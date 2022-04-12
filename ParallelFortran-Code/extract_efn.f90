PROGRAM arnoldi
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE transform
  USE init_jet_buffer
  USE pass1_linearized
  USE pass2
  USE solve
  USE diff_int
  IMPLICIT NONE

  INCLUDE 'hdf.f90'
  INCLUDE 'mpif.h'
  
  REAL(KIND=prec), DIMENSION(6) :: time
  INTEGER mpierr
  REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! Beginning of main program
  CALL initialize_mpi
  CALL read_parms
  CALL distribute_data
  CALL allocate_data_arrays
!  test_run = 1 !ESTIMATE OPTIMAL FFT, RATHER THAN SPENDING TIME MEASURING.
  CALL initialize_fft

  CALL define_coordinates
  CALL define_parameters

  linearized_run = 1
  CALL init_baseflow          ! READ IN BASEFLOW AND SET TIME STEP.
  dt = half*dt

  IF (xfringe+zfringe > 0) CALL init_buffer      ! INITIALIZE BUFFER REGION
  IF (linearized_run == 1) target_inflow = zero  ! PERTURBATIONS -> 0 IN BUFFER

  CALL extract_evectors(12)

  CALL deallocate_data_arrays
  CALL finalize_fft
  CALL finalize_mpi
  STOP 'Normal Termination'
! End of main program
CONTAINS
!=========================================================================
  SUBROUTINE extract_evectors(Nev)
    IMPLICIT NONE
    
    INTEGER Nev, i, j, jj, kk, lwork, info, n, nloc
    REAL(KIND=prec), DIMENSION(nx*local_nz_small*(ny+1)*3,Nev) :: v
    REAL(KIND=prec), DIMENSION(nx*local_nz_small*(ny+1)*3) :: wr, wi, rhs
    REAL(KIND=prec), DIMENSION(Nev,Nev)     :: m, h_arnoldi, vl, vr
    REAL(KIND=prec), DIMENSION(Nev)   :: lambdar, lambdai, res
    REAL(KIND=prec), DIMENSION(6*Nev) :: work
    REAL(KIND=prec) dlapy2, pdnorm2, tmp

    n    = nx*local_nz_small*(ny+1)*3
    nloc = nx*local_nz_small*(ny+1)*3

    ! READ IN ARNOLDI VECTORS.
    DO i = 1,Nev
       jj = 50 + i
       CALL read_evector(jj,lambdar(i),lambdai(i),v(1,i))
    END DO

    ! COMPUTE CORRELATION TENSOR AMONG ARNOLDI VECTORS.
    ! (ENSURE THAT ARNOLDI VECTORS ARE ORTHONORMAL.)
    m = zero
    DO i = 1,Nev
       DO j = 1,i
          tmp = SUM(v(1:nloc,i)*v(1:nloc,j))
          CALL MPI_ALLREDUCE(tmp,m(i,j),1,MPI_DOUBLE_PRECISION,MPI_SUM, &
               MPI_COMM_WORLD,mpierr)
       END DO
    END DO
    DO i = 1,Nev
       IF (me == 0) WRITE(*,930) (m(i,j),j=1,Nev)
       v(:,i) = v(:,i)/sqrt(m(i,i) + EPSILON(m(i,i)))
    END DO
    IF (me == 0) WRITE(*,*)
    930 FORMAT(100e14.6)

    ! COMPUTE PROJECTION OF RHS ONTO ARNOLDI VECTORS.
    H_arnoldi = zero
    DO i = 1,Nev
       CALL compute_rhs(v(1,i),rhs)
       DO j = 1,Nev
          tmp = SUM(rhs*v(1:nloc,j))
          CALL MPI_ALLREDUCE(tmp,h_arnoldi(j,i),1,MPI_DOUBLE_PRECISION, &
               MPI_SUM,MPI_COMM_WORLD,mpierr)
       END DO
    END DO
    DO i = 1,Nev
       IF (me == 0) WRITE(*,930) (H_arnoldi(i,j),j=1,Nev)
    END DO
    IF (me == 0) WRITE(*,*)

    ! COMPUTE SPECTRAL DECOMPOSITION OF H_ARNOLDI.
    lwork = 6*Nev
    CALL dgeev('N','V',Nev,H_arnoldi,Nev,lambdar,lambdai,vl,Nev,vr,Nev, &
         work,lwork,info)
    DO i = 1,Nev
       IF (me == 0) WRITE(*,*) lambdar(i), lambdai(i)
    END DO
    IF (me == 0) WRITE(*,*)
    DO i = 1,Nev
       IF (me == 0) WRITE(*,930) (vr(i,j),j=1,Nev)
    END DO
    IF (me == 0) WRITE(*,*)

    res = zero

    i = 1
    DO jj = 1,Nev
       IF (i > Nev) EXIT
       IF (lambdai(i) == zero) THEN
          ! REAL EIGENVECTOR
          wr = zero
          DO j = 1,Nev
             wr = wr + vr(j,i)*v(:,j)
          END DO
          CALL compute_rhs(wr,rhs)
          CALL daxpy(nloc,-lambdar(i),wr,1,rhs,1)
          res(i) = pdnorm2(MPI_COMM_WORLD,nloc,rhs,1)
          res(i) = res(i)/(ABS(lambdar(i)) + EPSILON(lambdar(i)))
          kk = 80 + i
          CALL write_evector(kk,lambdar(i),lambdai(i),wr)

          i = i + 1
       ELSE
          ! COMPLEX CONJUGATE EIGENVECTORS
          wr= zero
          wi= zero
          DO j = 1,Nev
             wr = wr + vr(j,i)  *v(:,j)
             wi = wi + vr(j,i+1)*v(:,j)
          END DO

          CALL compute_rhs(wr,rhs)
          CALL daxpy(nloc,-lambdar(i),wr,1,rhs,1)
          CALL daxpy(nloc, lambdai(i),wi,1,rhs,1)
          res(i) = pdnorm2(MPI_COMM_WORLD,nloc,rhs,1)

          CALL compute_rhs(wi,rhs)
          CALL daxpy(nloc,-lambdai(i),wr,1,rhs,1)
          CALL daxpy(nloc,-lambdar(i),wi,1,rhs,1)
          res(i) = dlapy2(res(i),pdnorm2(MPI_COMM_WORLD,nloc,rhs,1))
          res(i) = &
               res(i)/(dlapy2(lambdar(i),lambdai(i)) + EPSILON(lambdar(i)))
          kk = 80 + i
          CALL write_evector(kk,lambdar(i),lambdai(i),wr)
          kk = 81 + i
          CALL write_evector(kk,lambdar(i),-lambdai(i),wi)

          i = i + 2
       END IF
    END DO

    DO i = 1,Nev
       IF (me == 0) WRITE(*,990) i, lambdar(i), lambdai(i), res(i)
    END DO
990 FORMAT(i2,3e14.5)

  END SUBROUTINE extract_evectors
!=========================================================================
  SUBROUTINE compute_time_deriv(v,w)
    IMPLICIT NONE
    REAL(KIND=prec), DIMENSION(nx,local_nz_small,ny+1,3) :: v, w

    REAL(KIND=prec) t0, t1, scale, tmp1, tmp2
    REAL(KIND=prec), PARAMETER :: eps = 1.d-8
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
    INTEGER x,z,j,i

    ! INITIALIZE ARRAYS FOR VELOCITY FIELD AND RHS.
    data1 = zero
    data2 = zero

    ! SCALE INPUT VELOCITY FIELD AND PUT INTO data1 AND data2
    tmp1 = MAXVAL(ABS(v))
    CALL MPI_ALLREDUCE(tmp1,tmp2,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
         MPI_COMM_WORLD,mpierr)
    scale = eps/tmp2

    DO i = 1,3
       DO j = 1,ny+1
          CALL small_fft_real_to_complex(v(1,1,j,i),data1(1,1,j,i))
          SELECT CASE (i)
          CASE (2)
             IF (dy(j) /= zero) THEN
                data1(1:nz,1:local_nx,j,i) = &
                     (scale/sqrt(dy(j)))*data1(1:nz,1:local_nx,j,i)
             ELSE
                data1(1:nz,1:local_nx,j,i) = zero
             END IF
          CASE DEFAULT
             IF (dyh(j) /= zero) THEN
                data1(1:nz,1:local_nx,j,i) = &
                     (scale/sqrt(dyh(j)))*data1(1:nz,1:local_nx,j,i)
             ELSE
                data1(1:nz,1:local_nx,j,i) = zero
             END IF
          END SELECT
          data2(j,i,1:nz,1:local_nx) = data1(1:nz,1:local_nx,j,i)
       END DO
    END DO

!  COMPUTE INITIAL PRESSURE
    CALL nonlinear_linearized
    CALL solve_continuous_poisson
    ! RETURN INITIAL VELOCITY FIELD TO data1.
    DO X = 1,local_nx
       DO Z = 1,nz
          data1(z,x,1:ny+1,1:3+nscalar) = data2(1:ny+1,1:3+nscalar,z,x)
       END DO
    END DO
    IF (me == 0) WRITE(*,*) 'INITIALIZED PRESSURE'
    
    t0 = MPI_WTIME()
    t_total = zero
    DO step = 1,tend
       DO rk_step = 1,rk_last
!!$          CALL nonlinear
          CALL nonlinear_linearized
          CALL rhs_solve
       END DO
       T_TOTAL = T_TOTAL + DT
    END DO
    t1 = MPI_WTIME()
    IF (me==0) WRITE(*,*) 'FINAL TIME = ', t_total
    IF (me==0) write(*,999) t1-t0, output_step
999 format('WTIME = ',f12.6,'  OVER ', i6,'  STEPS')

    DO i = 1,3
       DO j = 1,ny+1
          SELECT CASE (i)
          CASE (2)
             data1(1:nz,1:local_nx,j,i) = &
                  (sqrt(dy(j))/scale)*data1(1:nz,1:local_nx,j,i)
          CASE DEFAULT
             data1(1:nz,1:local_nx,j,i) = &
                  (sqrt(dyh(j))/scale)*data1(1:nz,1:local_nx,j,i)
          END SELECT
          CALL small_fft_complex_to_real(data1(1,1,j,i),w(1,1,j,i))
       END DO
    END DO

  END SUBROUTINE compute_time_deriv
!=====================================================================
  SUBROUTINE solve_continuous_poisson
    IMPLICIT NONE

    COMPLEX(KIND=prec), DIMENSION(ny+1,3) :: cu, cru
    COMPLEX(KIND=prec), DIMENSION(ny+1)   :: du
    REAL(KIND=prec)    visc
    INTEGER            x, z

    DO X = 1,Klim
       DO Z = 1,Mlim

          !  READ IN THE VELOCITY FIELD AND NONLINEAR TERM.
          cu(1:ny+1,1:3)  = data2(1:ny+1,1:3,z,x)
          cru(1:ny+1,1:3) = - data1(z,x,1:ny+1,1:3)

          !  ADD THE VISCOUS TERM TO THE NONLINEAR TERM.         
          VISC = ONE/RE
          CALL d2_dy(ny+1,cu(1,1),du,d2h)
          cru(1:ny+1,1) = cru(1:ny+1,1) &
               + visc*(du(1:ny+1) - k2_me(z,x)*cu(1:ny+1,1))
          CALL d2_dy(ny,  cu(1,2),du,d2n)
          cru(1:ny,2) = cru(1:ny,2)  &
               + visc*(du(1:ny)   - k2_me(z,x)*cu(1:ny,2))
          CALL d2_dy(ny+1,cu(1,3),du,d2h)
          cru(1:ny+1,3) = cru(1:ny+1,3) &
               + visc*(du(1:ny+1) - k2_me(z,x)*cu(1:ny+1,3))

          IF (K2_ME(z,x) .EQ. ZERO) THEN
             !  SOLVE FOR PRESSURE IN THE MEAN MODE.
             data4(1:ny+1,z,x) = ZERO
          ELSE
             !  FOR EACH WAVENUMBER PAIR (X,Z), 
             !  TAKE THE DIVERGENCE OF THE RHS, 
             !  COMPUTE THE BOUNDARY CONDITIONS AND 
             !  SOLVE THE POISSON EQUATION FOR THE PRESSURE
             CALL strip_divergence(cru,data4(1:ny+1,z,x),one,x,z)
          END IF
       END DO
    END DO
  END SUBROUTINE solve_continuous_poisson
!=========================================================================
  SUBROUTINE compute_rhs(v,w)
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(nx,local_nz_small,ny+1,3) :: v, w

    COMPLEX(KIND=prec), DIMENSION(1:ny+1,3) :: cu
    COMPLEX(KIND=prec), DIMENSION(1:ny+1) :: d2u, cp
    COMPLEX(KIND=prec) :: ctmp
    REAL(KIND=prec) tmp_visc, tmp_diff, scale, tmp1, tmp2
    REAL(KIND=prec), PARAMETER :: eps = 1.d-8
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
    INTEGER x,z,j,i

    ! INITIALIZE ARRAYS FOR VELOCITY FIELD AND RHS.
    data1 = zero
    data2 = zero

    ! SCALE INPUT VELOCITY FIELD AND PUT INTO data1 AND data2
    tmp1 = MAXVAL(ABS(v))
    CALL MPI_ALLREDUCE(tmp1,tmp2,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
         MPI_COMM_WORLD,mpierr)
    scale = eps/tmp2
!    WRITE(*,*) me, scale
    DO i = 1,3
       DO j = 1,ny+1
          CALL small_fft_real_to_complex(v(1,1,j,i),data1(1,1,j,i))
          SELECT CASE (i)
          CASE (2)
             IF (dy(j) /= zero) THEN
                data1(1:nz,1:local_nx,j,i) = &
                     (scale/sqrt(dy(j)))*data1(1:nz,1:local_nx,j,i)
             ELSE
                data1(1:nz,1:local_nx,j,i) = zero
             END IF
          CASE DEFAULT
             IF (dyh(j) /= zero) THEN
                data1(1:nz,1:local_nx,j,i) = &
                     (scale/sqrt(dyh(j)))*data1(1:nz,1:local_nx,j,i)
             ELSE
                data1(1:nz,1:local_nx,j,i) = zero
             END IF
          END SELECT
          data2(j,i,1:nz,1:local_nx) = data1(1:nz,1:local_nx,j,i)
       END DO
    END DO

    ! COMPUTE LINEARIZED FORM OF NONLINEAR TERM.
    CALL nonlinear_linearized

    ! ADD VISCOUS PART OF RIGHT HAND SIDE AND REMOVE DIVERGENCE.
    tmp_visc = one/re
    DO x = 1,local_nx
       DO z = 1,nz
          ! COMPUTE VISCOUS TERM.
          CALL d2_dy(ny+1,data2(1,1,z,x),d2u,d2h)
          cu(1:ny+1,1) = - data1(z,x,1:ny+1,1) &
               + tmp_visc*(d2u(1:ny+1) - k2_me(z,x)*data2(1:ny+1,1,z,x))
          CALL d2_dy(ny,data2(1,2,z,x),d2u,d2n)
          cu(1:ny,2) = - data1(z,x,1:ny,2) &
               + tmp_visc*(d2u(1:ny) - k2_me(z,x)*data2(1:ny,2,z,x))
          CALL d2_dy(ny+1,data2(1,3,z,x),d2u,d2h)
          cu(1:ny+1,3) = - data1(z,x,1:ny+1,3) &
               + tmp_visc*(d2u(1:ny+1) - k2_me(z,x)*data2(1:ny+1,3,z,x))
          ! REMOVE DIVERGENCE
          IF (k2_me(z,x) == zero) THEN
             ctmp = DBLE(SUM(dy(1:ny)*cu(1:ny,2)))
             cu(1:ny,2) = ctmp
          ELSE
             CALL strip_divergence(cu,cp,one,x,z)
          END IF
          ! UPLOAD TO data1
          data1(z,x,1:ny+1,1:3) = cu
          data1(z,x,ny+1,2) = zero
          ! ADD DIFFUSIVE TERMS FOR SCALAR
          DO j = 1,nscalar
             tmp_diff = one/(re*pr(j))
             CALL d2_dy(ny+1,data2(1,3+j,z,x),d2u,d2h)
             data1(z,x,1:ny+1,3+j) = - data1(z,x,1:ny+1,3+j) &
                  + tmp_diff*(d2u(1:ny+1) - k2_me(z,x)*data2(1:ny+1,3+j,z,x))
          END DO
       END DO
    END DO

    DO i = 1,3
       DO j = 1,ny+1
          SELECT CASE (i)
          CASE (2)
             data1(1:nz,1:local_nx,j,i) = &
                  (sqrt(dy(j))/scale)*data1(1:nz,1:local_nx,j,i)
          CASE DEFAULT
             data1(1:nz,1:local_nx,j,i) = &
                  (sqrt(dyh(j))/scale)*data1(1:nz,1:local_nx,j,i)
          END SELECT
          CALL small_fft_complex_to_real(data1(1,1,j,i),w(1,1,j,i))
       END DO
    END DO

  END SUBROUTINE compute_rhs
!=========================================================================
  SUBROUTINE read_evector(n,cr,ci,evector)
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(nx,nz,ny+1)        :: data
    INTEGER, PARAMETER :: time_rank = 1, time_length = 6, vel_rank = 3
    INTEGER :: plane, i, status, dsgdata, vel_length(vel_rank)

    REAL(KIND=prec), DIMENSION(nx,local_nz_small,ny+1,3) :: evector
    REAL(KIND=prec) :: cr,ci
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
    INTEGER :: n
    CHARACTER(LEN=13) :: filename

    vel_length(1) = nx
    vel_length(2) = nz
    vel_length(3) = ny+1

    IF (me == 0) THEN
       ! GENERATE FILENAME ACCORDING TO EIGENVECTOR NUMBER.
       filename = 'evector' // char(48 + mod(int(n/10),10))// &
            char(48 + mod(n,10)) // '.hdf'
       WRITE(*,*) 'READING ', filename

       ! READ IN THE TIME
       ! TIME ARRAY HAS SIX ELEMENTS: time, real/imag part of growth rate,
       ! mean pressure gradient, reynolds number and prandtl number.
       status = dsgdata(filename,time_rank,time_length,time)
       cr = time(2)
       ci = time(3)
    END IF

    evector = zero
    DO i = 1,3+nscalar
       ! READ IN THE VELOCITY/SCALAR ARRAY
       IF (me == 0) THEN
          status = dsgdata(filename,vel_rank,vel_length,data)
          WRITE(*,*) 'READ FIELD ', i, ' WITH STATUS ', status
       END IF
       ! TRANSFORM INTO FOURIER SPACE AND STORE IN DATA1
       DO plane = 1,ny+1
          CALL MPI_SCATTERV( &
               data(1,1,plane),nx*nz_small_proc,nx*z_small_start, &
               MPI_DOUBLE_PRECISION, &
               evector(1,1,plane,i),nx*local_nz_small,MPI_DOUBLE_PRECISION, &
               0,MPI_COMM_WORLD,mpierr)
          SELECT CASE (i)
          CASE (2)
             IF (dy(plane) /= zero) THEN
                evector(1:nx,1:local_nz_small,plane,i) = &
                     sqrt(dy(plane))*evector(1:nx,1:local_nz_small,plane,i)
             END IF
          CASE DEFAULT
             IF (dyh(plane) /= zero) THEN
                evector(1:nx,1:local_nz_small,plane,i) = &
                     sqrt(dyh(plane))*evector(1:nx,1:local_nz_small,plane,i)
             END IF
          END SELECT 
      END DO
    END DO

  END SUBROUTINE read_evector
!=========================================================================
  SUBROUTINE write_evector(n,realpart,imagpart,evector)
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(nx,nz,ny+1) :: data
    REAL(KIND=prec), DIMENSION(6)          :: time

    REAL(KIND=prec), DIMENSION(nx,local_nz_small,ny+1,3) :: evector
    REAL(KIND=prec) :: vel_xcoord(nx), vel_zcoord(nz), realpart, imagpart
    REAL(KIND=prec) :: cr,ci    
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
    INTEGER, PARAMETER :: time_rank = 1, time_length = 6, vel_rank = 3
    INTEGER :: n, plane, i, x, z, status, dspdata, dssdims, dssnt, dsadata, &
         dssdisc, vel_length(vel_rank)
    CHARACTER(LEN=13) :: filename

    vel_length(1) = nx
    vel_length(2) = nz
    vel_length(3) = ny+1

    DO x = 1,nx
       vel_xcoord(x) = DBLE(x-1)*lx/nx
    END DO
    DO z = 1,nz
       vel_zcoord(z) = DBLE(z-1)*lz/nz
    END DO

    IF (me == 0) THEN
       ! GENERATE FILENAME ACCORDING TO EIGENVECTOR NUMBER.
       filename = 'evector' // char(48 + mod(int(n/10),10))// &
            char(48 + mod(n,10)) // '.hdf'
       WRITE(*,*) 'WRITING ', filename

       ! CONVERT COMPUTED GROWTH RATE = exp((cr + i*ci)*T_TOTAL) TO (cr,ci)
       ! WHICH ARE THE EIGENVALUES OF THE LINEARIZED NAVIER--STOKES EQNS.
       cr = realpart
       ci = imagpart
       WRITE(*,980) cr,ci
980    FORMAT('EIGENVALUE = ( ',F14.10,', ',F14.10,')')

       ! WRITE THE TIME
       ! TIME ARRAY HAS SIX ELEMENTS: efn #, real/imag part of growth rate,
       ! mean pressure gradient, reynolds number and prandtl number.
       time(1) = n
       time(2) = cr
       time(3) = ci

       ! SETUP THE HDF FILE WHICH WILL HOLD THE EIGENVECTOR.
       status = dssdims(time_rank,time_length)
       status = dssnt(DFNT_FLOAT64)
       status = dspdata(filename,time_rank,time_length,time)
    END IF

    DO i = 1,3+nscalar
       ! UNSCALE THE EIGENVECTOR (DIVIDE BY SQRT OF INTEGRATION WEIGHTS)
       ! AND GATHER ONTO MASTER NODE.
       DO plane = 1,ny+1
          SELECT CASE (i)
          CASE (2)
             IF (dy(plane) /= zero) THEN
                evector(1:nx,1:local_nz_small,plane,i) = &
                     evector(1:nx,1:local_nz_small,plane,i)/sqrt(dy(plane))
             ELSE
                evector(1:nx,1:local_nz_small,plane,i) = zero
             END IF
          CASE DEFAULT
             IF (dyh(plane) /= zero) THEN
                evector(1:nx,1:local_nz_small,plane,i) = &
                     evector(1:nx,1:local_nz_small,plane,i)/sqrt(dyh(plane))
             ELSE
                evector(1:nx,1:local_nz_small,plane,i) = zero
             END IF
          END SELECT
          CALL MPI_GATHERV( &
               evector(1,1,plane,i),nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,&
               data(1,1,plane),nx*nz_small_proc,nx*z_small_start, &
               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
       END DO
       ! WHEN FINISHED GATHERING ONE COMPONENT (i.e. U, V OR W) OF EVECTOR,
       ! SAVE THAT COMPONENT TO THE hdf FILE.
       IF (me == 0) THEN
          status = dssdims(vel_rank,vel_length)
          status = dssnt(DFNT_FLOAT64)
          status = dssdisc(1,vel_length(1),vel_xcoord) ! WRITE COORDINATES.
          status = dssdisc(2,vel_length(2),vel_zcoord)
          status = dssdisc(3,vel_length(3),yh)
          status = dsadata(filename,vel_rank,vel_length,data) ! WRITE EVECTOR.
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    END DO


  END SUBROUTINE write_evector
!==========================================================================
!
!     matrix vector subroutine
!
!     The matrix used is the 2 dimensional convection-diffusion 
!     operator discretized using central difference.
!
  subroutine av (nx, v, w)
    integer           nx, j, lo
    Double precision  v(nx*nx), w(nx*nx), one, h2
    parameter         (one = 1.0D+0)
    external          daxpy
    !
    !     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block 
    !     tridiagonal matrix
    !
    !                  | T -I          | 
    !                  |-I  T -I       |
    !             OP = |   -I  T       |
    !                  |        ...  -I|
    !                  |           -I T|
    !
    !     derived from the standard central difference discretization 
    !     of the 2 dimensional convection-diffusion operator 
    !     (Laplacian u) + rho*(du/dx) on a unit square with zero boundary 
    !     condition.
    !
    !     When rho*h/2 <= 1, the discrete convection-diffusion operator 
    !     has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX 
    !     eigenvalues.
    !
    !     The subroutine TV is called to compute y<---T*x.
    !
    !
    h2 = one / dble((nx+1)*(nx+1))
    !
    call tv(nx,v(1),w(1))
    call daxpy(nx, -one/h2, v(nx+1), 1, w(1), 1)
    !
    do j = 2, nx-1
       lo = (j-1)*nx
       call tv(nx, v(lo+1), w(lo+1))
       call daxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
       call daxpy(nx, -one/h2, v(lo+nx+1), 1, w(lo+1), 1)
    end do
    !
    lo = (nx-1)*nx
    call tv(nx, v(lo+1), w(lo+1))
    call daxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
    !
  end subroutine AV
  !=========================================================================
  subroutine tv (nx, x, y)
    !
    integer           nx, j 
    Double precision  x(nx), y(nx), h, h2, dd, dl, du
    !
    Double precision  one, zero, rho
    parameter         (one = 1.0D+0, zero = 0.0D+0, rho = 21.99D+0)
    !
    !     Compute the matrix vector multiplication y<---T*x
    !     where T is a nx by nx tridiagonal matrix with DD on the 
    !     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
    !
    !     When rho*h/2 <= 1, the discrete convection-diffusion operator 
    !     has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX 
    !     eigenvalues.
    !
    h   = one / dble(nx+1)
    h2  = one / dble((nx+1)*(nx+1))
    dd  = 4.0D+0 / h2
    dl  = -one / h2 - 0.5D+0*rho / h
    du  = -one / h2 + 0.5D+0*rho / h
    ! 
    y(1) =  dd*x(1) + du*x(2)
    do j = 2,nx-1
       y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
    end do
    y(nx) =  dl*x(nx-1) + dd*x(nx) 
  end subroutine tv
!=========================================================================
  SUBROUTINE create_initial_vector
    IMPLICIT NONE

    COMPLEX(KIND=prec), DIMENSION(1:ny+1,3) :: cu
    COMPLEX(KIND=prec), DIMENSION(1:ny+1)   :: cp
    REAL(KIND=prec) :: TEMP(NY+1,6), ushape(1:ny+1), vshape(1:ny)
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    INTEGER :: x, z

    data1 = zero
    data2 = zero
    data4 = zero

    ! INITIALIZE WHITE RANDOM VELOCITY FIELD.
    DO x = 1,local_nx
       DO z = 1,nz
!!$    DO x = 2,2
!!$       DO z = 1,1
          IF (k2_me(z,x) .ne. zero) THEN
             CALL RANDOM_NUMBER(TEMP)
             ushape(1:ny+1) = tanh(10.*(yh(1:ny+1)-y(1))) &
                  *tanh(10.*(yh(ny+1)-yh(1:ny+1)))
             vshape(1:ny)   = tanh(10.*(y(1:ny)-y(1)))    &
                  *tanh(10.*(y(ny)-y(1:ny)))
             data1(z,x,1:ny+1,1) = two*ushape(1:ny+1)          &
                  *cmplx(temp(1:ny+1,1)-half,temp(1:ny+1,2)-half)
             data1(z,x,1:ny,2)   = two*vshape(1:ny)            &
                  *cmplx(temp(1:ny,3)-half,temp(1:ny,4)-half)
             data1(z,x,1:ny+1,3) = two*ushape(1:ny+1)          &
                  *cmplx(temp(1:ny+1,5)-half,temp(1:ny+1,6)-half)
             ! REMOVE DIVERGENCE FROM INITIAL VECTOR AND UPLOAD TO data1
             cu = data1(z,x,1:ny+1,1:3)
             CALL strip_divergence(cu,cp,one,x,z)
             data1(z,x,1:ny+1,1:3) = cu
          ELSE
             ! ZERO OUT MEAN VELOCITY PROFILE SINCE THIS SIMULATION
             ! IS LINEARIZED ABOUT A BASEFLOW.
             data1(z,x,1:ny+1,1:3) = zero
          END IF
       END DO
    END DO
  END SUBROUTINE create_initial_vector
END PROGRAM arnoldi
