PROGRAM arnoldi
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE transform
  USE init_jet_buffer
!!$  USE pass1
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
!!$  tend = INT(one/dt)
  tend = MAX(INT(0.1d0/dt),10)
!!$  output_step = 1
  output_step = tend
  IF (me == 0) WRITE(*,*) 'DT = ', dt
  IF (me == 0) WRITE(*,*) tend, ' STEPS PER RHS CALL'

  IF (xfringe+zfringe > 0) CALL init_buffer      ! INITIALIZE BUFFER REGION
  IF (linearized_run == 1) target_inflow = zero  ! PERTURBATIONS -> 0 IN BUFFER

!!$  CALL read_snap
!!$  CALL compute_time_deriv
!!$  CALL write_ddt
!!$
  CALL arpack_call(nx*local_nz_small*(ny+1)*3,nx*local_nz_small*(ny+1)*3)

  CALL deallocate_data_arrays
  CALL finalize_fft
  CALL finalize_mpi
  STOP 'Normal Termination'
! End of main program
CONTAINS
!=========================================================================
  SUBROUTINE arpack_call(maxn,ldv)
    IMPLICIT NONE
    
    INCLUDE 'debug.h'
    INCLUDE 'stat.h'

    !     %---------------%
    !     | MPI INTERFACE |
    !     %---------------%
    integer           rc, nloc 
    !     %-----------------------------%
    !     | Define maximum dimensions   |
    !     | for all arrays.             |
    !     | MAXN:   Maximum dimension   |
    !     |         of the A allowed.   |
    !     | MAXNEV: Maximum NEV allowed |
    !     | MAXNCV: Maximum NCV allowed |
    !     %-----------------------------%
    !
    integer           maxn, maxnev, maxncv, ldv
    parameter         (maxnev=12, maxncv=25)
    !
    !     %--------------%
    !     | Local Arrays |
    !     %--------------%
    !
    integer           iparam(11), ipntr(14)
    logical           select(maxncv)
    Double precision  ax(maxn), resid(maxn), v(ldv,maxncv), &
                      d(maxncv,3), workd(3*maxn), workev(3*maxncv), &
                      workl(3*maxncv*maxncv+6*maxncv)
    !
    !     %---------------%
    !     | Local Scalars |
    !     %---------------%
    !
    character         bmat*1, which*2
    integer           ido, n, nev, ncv, lworkl, info, j, &
                      ierr, nconv, maxitr, ishfts, mode, jj, i, dsnum
    Double precision  tol, sigmar, sigmai
    logical           first, rvec
    !
    !     %----------------------------------------------%
    !     | Local Buffers needed for MPI communication |
    !     %----------------------------------------------%
    !
    Double precision  mv_buf(maxn)
    !
    !     %------------%
    !     | Parameters |
    !     %------------%
    !
    Double precision   zero, one
    parameter         (zero = 0.0D+0)
    parameter         (one  = 1.0D+0)
    !
    !     %-----------------------------%
    !     | BLAS & LAPACK routines used |
    !     %-----------------------------%
    !
    Double precision  dlapy2, pdnorm2, dnrm2
    external          dlapy2, daxpy, pdnorm2, dnrm2
    !
    !     %--------------------%
    !     | Intrinsic function |
    !     %--------------------%
    !
    intrinsic         abs
    double precision     tmpr, tmpi
    integer              status
    CHARACTER(LEN=13) :: filename
    !
    !     %-----------------------%
    !     | Executable Statements |
    !     %-----------------------%
    ! 
    ndigit = -3
    logfil = 6
    mnaupd = 2
    mnaup2 = 2
    !
    !     %--------------------------------------------------%
    !     | The number NX is the number of interior points   |
    !     | in the discretization of the 2-dimensional       |
    !     | convection-diffusion operator on the unit        |
    !     | square with zero Dirichlet boundary condition.   | 
    !     | The number N(=NX*NX) is the dimension of the     |
    !     | matrix.  A standard eigenvalue problem is        |
    !     | solved (BMAT = 'I').  NEV is the number of       |
    !     | eigenvalues to be approximated.  The user can    |
    !     | modify NX, NEV, NCV, WHICH to solve problems of  |
    !     | different sizes, and to get different parts of   |
    !     | the spectrum.  However, The following            |
    !     | conditions must be satisfied:                    |
    !     |                   N <= MAXN                      |
    !     |                 NEV <= MAXNEV                    |
    !     |           NEV + 2 <= NCV <= MAXNCV               | 
    !     %--------------------------------------------------% 
    !
    n     = nx*local_nz_small*(ny+1)*3 
    nloc  = nx*local_nz_small*(ny+1)*3 
    nev   = 12
    ncv   = 25
    if ( n .gt. maxn ) then
       print *, ' ERROR with _NDRV1: N is greater than MAXN '
       STOP
    else if ( nev .gt. maxnev ) then
       print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
       STOP 
    else if ( ncv .gt. maxncv ) then
       print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
       STOP 
    end if
    bmat  = 'I'
    which = 'LM'
    !
    !     %-----------------------------------------------------%
    !     | The work array WORKL is used in DNAUPD as           |  
    !     | workspace.  Its dimension LWORKL is set as          |
    !     | illustrated below.  The parameter TOL determines    |
    !     | the stopping criterion. If TOL<=0, machine          |
    !     | precision is used.  The variable IDO is used for    |
    !     | reverse communication, and is initially set to 0.   |
    !     | Setting INFO=0 indicates that a random vector is    |
    !     | generated in DNAUPD to start the Arnoldi iteration. | 
    !     %-----------------------------------------------------%
    !
    lworkl  = 3*ncv**2+6*ncv 
    tol    = 1.d-8
    ido    = 0
    info   = 0
    !
    !     %---------------------------------------------------%
    !     | This program uses exact shifts with respect to    |
    !     | the current Hessenberg matrix (IPARAM(1) = 1).    |
    !     | IPARAM(3) specifies the maximum number of Arnoldi |
    !     | iterations allowed.  Mode 1 of DNAUPD is used     |
    !     | (IPARAM(7) = 1). All these options can be changed |
    !     | by the user. For details see the documentation in |
    !     | DNAUPD.                                           |
    !     %---------------------------------------------------%
    !
    ishfts = 1
    maxitr = 2
    mode   = 1
    !
    iparam(1) = ishfts
    iparam(3) = maxitr 
    iparam(7) = mode
    !
    !      SET INITIAL VALUE FOR RESIDUAL.
    !      THIS ACTS AS A STARTING VALUE FOR THE ARNOLDI ALGORITHM
    !
    info   = 1
    ! INITIALIZE RESIDUAL VECTOR.
    resid = zero
    ! IF OLD ARNOLDI VECTORS EXIST, READ IN AND ADD TO INITIAL VECTOR.
    filename = 'evector51.hdf'
    status = dsnum(filename)
    IF (status >= 1) THEN
       IF (me == 0) WRITE(*,*) 'READING ', filename
       CALL read_evector(51,tmpr,tmpi,resid)
    ELSE
       ! CREATE DIVERGENCE-FREE INITIAL VECTOR FROM WHITE NOISE.
       CALL create_initial_vector  ! CREATE DIV-FREE INITIAL VECTOR.
       jj = 1
       ax = zero
       DO i = 1,3
          DO j = 1,ny+1
             SELECT CASE (i)
             CASE (2)
                data1(1:nz,1:local_nx,j,i) = &
                     sqrt(dy(j))*data1(1:nz,1:local_nx,j,i)
             CASE DEFAULT
                data1(1:nz,1:local_nx,j,i) = &
                     sqrt(dyh(j))*data1(1:nz,1:local_nx,j,i)
             END SELECT
             CALL small_fft_complex_to_real(data1(1,1,j,i),resid(jj))
             jj = jj + nx*local_nz_small
          END DO
       END DO
    END IF
    !     %-------------------------------------------%
    !     | M A I N   L O O P (Reverse communication) | 
    !     %-------------------------------------------%
    !
    DO jj = 1,1000000
       IF (mod(jj,100) == 0) WRITE(*,*) 'DONE MATRIX MULTIPLICATION NO. ', jj
       !
       !        %---------------------------------------------%
       !        | Repeatedly call the routine DNAUPD and take |
       !        | actions indicated by parameter IDO until    |
       !        | either convergence is indicated or maxitr   |
       !        | has been exceeded.                          |
       !        %---------------------------------------------%
       !
!!$       CALL dnaupd ( ido, bmat, n, which, nev, tol, resid,  &
!!$            ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,  &
!!$            info )
       call pdnaupd(MPI_COMM_WORLD, ido, bmat, nloc, which, nev, tol, resid, &
            ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
       !
!!$       ! IF AT maxitr OR EVERY 100 RHS COMPUTATIONS, SAVE STATE. 
!!$       IF ((ido == -2) .or. (mod(jj,100) == 0)) CALL write_saved_state
       !
       ! CHECK TO SEE IF WE ARE DONE COMPUTING THE RHS.
       IF ((ido .ne. -1) .and. (ido .ne. 1)) EXIT

       !  IF NOT, COMPUTE RHS AGAIN
       !         if (ido .eq. -1 .or. ido .eq. 1) then
       !
       !           %-------------------------------------------%
       !           | Perform matrix vector multiplication      |
       !           |                y <--- OP*x                |
       !           | The user should supply his/her own        |
       !           | matrix vector multiplication routine here |
       !           | that takes workd(ipntr(1)) as the input   |
       !           | vector, and return the matrix vector      |
       !           | product to workd(ipntr(2)).               | 
       !           %-------------------------------------------%
       !
       CALL compute_time_deriv(workd(ipntr(1)), workd(ipntr(2)))

       
       !
       !           %-----------------------------------------%
       !           | L O O P   B A C K to call DNAUPD again. |
       !           %-----------------------------------------%
       !
       !
    end do
    ! 
    !     %----------------------------------------%
    !     | Either we have convergence or there is |
    !     | an error.                              |
    !     %----------------------------------------%
    !
    if ( info .lt. 0 ) then
       !
       !        %--------------------------%
       !        | Error message, check the |
       !        | documentation in DNAUPD. |
       !        %--------------------------%
       !
       IF (me == 0) THEN
          print *, ' '
          print *, ' Error with _naupd, info = ', info
          print *, ' Check the documentation of _naupd'
          print *, ' '
       END IF
       !
    else 
       !
       !        %-------------------------------------------%
       !        | No fatal errors occurred.                 |
       !        | Post-Process using DNEUPD.                |
       !        |                                           |
       !        | Computed eigenvalues may be extracted.    |
       !        |                                           |
       !        | Eigenvectors may also be computed now if  |
       !        | desired.  (indicated by rvec = .true.)    |
       !        %-------------------------------------------%
       !
       rvec = .true.
       !
!!$       call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv,  &
!!$            sigmar, sigmai, workev, bmat, n, which, nev, tol,  &
!!$            resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
!!$            lworkl, ierr )
       call pdneupd (MPI_COMM_WORLD, rvec, 'A', select, d, d(1,2), v, ldv, &
            sigmar, sigmai, workev, bmat, nloc, which, nev, tol, &
            resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
            lworkl, ierr )
       !
       !        %-----------------------------------------------%
       !        | The real part of the eigenvalue is returned   |
       !        | in the first column of the two dimensional    |
       !        | array D, and the imaginary part is returned   |
       !        | in the second column of D.  The corresponding |
       !        | eigenvectors are returned in the first NEV    |
       !        | columns of the two dimensional array V if     |
       !        | requested.  Otherwise, an orthogonal basis    |
       !        | for the invariant subspace corresponding to   |
       !        | the eigenvalues in D is returned in V.        |
       !        %-----------------------------------------------%
       !
       if ( ierr .ne. 0) then
          !
          !           %------------------------------------%
          !           | Error condition:                   |
          !           | Check the documentation of DNEUPD. |
          !           %------------------------------------%
          !
          IF (me == 0) THEN
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
          END IF
          !
       else 
          !
          first  = .true.
          nconv  = iparam(5)
          do j=1, nconv
             !
             !               %---------------------------%
             !               | Compute the residual norm |
             !               |                           |
             !               |   ||  A*x - lambda*x ||   |
             !               |                           |
             !               | for the NCONV accurately  |
             !               | computed eigenvalues and  |
             !               | eigenvectors.  (iparam(5) |
             !               | indicates how many are    |
             !               | accurate to the requested |
             !               | tolerance)                |
             !               %---------------------------%
             !
             IF (d(j,2) .eq. zero)  THEN
                !
                !                  %--------------------%
                !                  | Ritz value is real |
                !                  %--------------------%
                !
                CALL compute_time_deriv(v(1,j),ax)
!!$                call av(nx, v(1,j), ax)
                CALL daxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
                d(j,3) = pdnorm2(MPI_COMM_WORLD, nloc, ax, 1)
!!$                d(j,3) = dnrm2(n, ax, 1)
                d(j,3) = d(j,3) / (abs(d(j,1)) + EPSILON(d(j,1)))
                !
             ELSEIF (first) THEN
                !
                !                  %------------------------%
                !                  | Ritz value is complex. |
                !                  | Residual of one Ritz   |
                !                  | value of the conjugate |
                !                  | pair is computed.      | 
                !                  %------------------------%
                !        
                CALL compute_time_deriv(v(1,j),ax)
                CALL daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                CALL daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                   d(j,3) = pdnorm2(MPI_COMM_WORLD, nloc, ax, 1)
!!$                d(j,3) = dnrm2(n, ax, 1)
                CALL compute_time_deriv(v(1,j),ax)
!!$                call av(nx, v(1,j+1), ax)
                CALL daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                CALL daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                d(j,3) = dlapy2(d(j,3), pdnorm2(MPI_COMM_WORLD,nloc,ax,1) )
!!$                d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                d(j,3) = d(j,3) / (dlapy2(d(j,1),d(j,2)) + EPSILON(d(j,1)))
                d(j+1,3) = d(j,3)
                first = .false.
             ELSE
                first = .true.
             END IF
             !
          END DO
          !
          !            %-----------------------------%
          !            | Display computed residuals. |
          !            %-----------------------------%
          !
             CALL pdmout(MPI_COMM_WORLD, 6, nconv, 3, d, maxncv, -6, &
                  'Ritz values (Real,Imag) and direct residuals')
!!$          call dmout(6, nconv, 3, d, maxncv, -6, &
!!$               'Ritz values (Real,Imag) and relative residuals')
       end if
       !
       !        %-------------------------------------------%
       !        | Print additional convergence information. |
       !        %-------------------------------------------%
       !
       IF (me == 0) THEN
          IF ( info .eq. 1) THEN
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
          ELSE IF ( info .eq. 3) THEN
             print *, ' ' 
             print *, ' No shifts could be applied during implicit', &
                  ' Arnoldi update, try increasing NCV.'
             print *, ' '
          END IF
          !
          print *, ' '
          print *, ' _NDRV1 '
          print *, ' ====== '
          print *, ' ' 
          print *, ' Size of the matrix is ', n
          print *, ' The number of Ritz values requested is ', nev
          print *, ' The number of Arnoldi vectors generated', &
               ' (NCV) is ', ncv
          print *, ' What portion of the spectrum: ', which
          print *, ' The number of converged Ritz values is ',  &
               nconv 
          print *, ' The number of Implicit Arnoldi update', &
               ' iterations taken is ', iparam(3)
          print *, ' The number of OP*x is ', iparam(9)
          print *, ' The convergence criterion is ', tol
          print *, ' '
          !
       END IF
    END IF
    !
    !
    ! WRITE EIGENVECTORS OUT TO A FILE.
    !
    DO i = 1,nconv
       CALL write_evector(i,d(i,1),d(i,2),v(1,i))
    END DO
    IF (nconv < nev) THEN
       DO i = 1,nev
          jj = 50 + i
          CALL write_evector(jj,one,one,v(1,i))
       END DO
    END IF

    !
    !     %---------------------------%
    !     | Done with program dndrv1. |
    !     %---------------------------%
    !
  END SUBROUTINE arpack_call
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
       cr = log(sqrt(realpart**2 + imagpart**2) + EPSILON(realpart))/t_total
       ci = atan2(imagpart,realpart)/t_total
!!$       cr = realpart
!!$       ci = imagpart
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
!=========================================================================
  SUBROUTINE create_initial_vector
    IMPLICIT NONE

    COMPLEX(KIND=prec), DIMENSION(1:ny+1,3) :: cu
    COMPLEX(KIND=prec), DIMENSION(1:ny+1)   :: cp
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
          CALL random_divfree_vector(0,x,z,1.d0,cu)
          IF (k2_me(z,x) .ne. zero) THEN
             ! REMOVE DIVERGENCE FROM INITIAL VECTOR AND UPLOAD TO data1
             CALL strip_divergence(cu,cp,one,x,z)
             data1(z,x,1:ny+1,1:3) = cu
          ELSE
             ! ZERO OUT MEAN VELOCITY PROFILE SINCE THIS SIMULATION
             ! IS LINEARIZED ABOUT A BASEFLOW.
             data1(z,x,1:ny+1,1) = DBLE(cu(1:ny+1,1))
             data1(z,x,1:ny+1,3) = DBLE(cu(1:ny+1,3))
          END IF
       END DO
    END DO

  END SUBROUTINE create_initial_vector
  !=====================================================================
  SUBROUTINE random_divfree_vector(filter,p,q,ampl,cu)
    IMPLICIT NONE
    COMPLEX(KIND=prec), DIMENSION(ny+1,3) :: cu
    INTEGER :: filter, p, q !! x- AND z-WAVENUMBERS OF THIS FOURIER MODE.
    REAL(KIND=prec) :: ampl !! APPROX AMPLITUDE OF VECTOR BEFORE FILTERING.

    REAL(KIND=prec), DIMENSION(ny+1,4) :: temp, temp2
    REAL(KIND=prec), DIMENSION(ny+1)   :: ushape
    REAL(KIND=prec), DIMENSION(ny)     :: vshape
    REAL(KIND=prec) :: tmp
    INTEGER         :: k

    cu = zero
    ushape = one
    vshape = one

    !! GENERATE RANDOM VECTORS
    temp = zero
    CALL RANDOM_NUMBER(TEMP)

    IF (filter .eq. 1) THEN
       !! SMOOTH NOISE IN WALL-NORMAL DIRECTION.
       DO k = 2,ny
          temp2(k,:) = 0.25d0*(temp(k-1,:)+temp(k+1,:)) + half*temp(k,:)
       END DO
       temp(2:ny,:) = temp2(2:ny,:)

       !! SCALE NOISE TO DECREASE AMPLITUDE NEXT TO WALLS.
       IF (bl_grid .eq. 0) THEN
          ushape = (tanh(10.d0*(yh-yh(1)))*tanh(10.d0*(yh(ny+1)-yh)))
          vshape = (tanh(10.d0*(y -y(1))) *tanh(10.d0*(y(ny)   -y)))**2
       ELSE
          ushape = tanh(10.d0*(yh-yh(1)))
          vshape = tanh(10.d0*(y -y(1)))**2
       END IF
    END IF

    !! INITIALIZE WITH RANDOM VERTICAL VELOCITY AND VERTICAL VORTICITY 
    !! PROFILES.  PLACE THESE IN cu(:,1) AND cu(:,2), RESPECTIVELY.
    !! COMPUTE STREAMWISE AND SPANWISE VELOCITY THAT WILL GIVE A 
    !! DIVERGENCE-FREE VELOCITY FIELD.  PLACE THESE IN cu(:,1) AND 
    !! cu(:,3), RESPECTIVELY.

    IF (k2_me(q,p) .ne. zero) THEN

       !! INITIALIZE VELOCITY/VORTICITY PROFILES WITH RANDOM AMPLITUDE/PHASE.
       temp = two*(temp-half)
       cu(1:ny+1,1) = ushape*CMPLX(temp(1:ny+1,1),temp(1:ny+1,2))
       cu(1:ny,2) = vshape*CMPLX(temp(1:ny,3),temp(1:ny,4))
       cu(ny+1,2) = zero

       !! COMPUTE SPANWISE VELOCITY PROFILE FIRST.
       cu(1,3) = zero
       cu(2:ny,3) = (ikx_me(p)*cu(2:ny,1) &
            + ikz(q)*dh(2:ny)*(cu(2:ny,2) - cu(1:ny-1,2)))/k2_me(q,p)
       cu(ny+1,3) = zero

       !! COMPUTE STREAMWISE VELOCITY PROFILE NEXT.
       cu(1,1) = zero
       cu(2:ny,1) = (- ikz(q)*cu(2:ny,1) &
            + ikx_me(p)*dh(2:ny)*(cu(2:ny,2) - cu(1:ny-1,2)))/k2_me(q,p)
       cu(ny+1,1) = zero

       !! COMPUTE 2-NORM OF RANDOM VECTOR, NORMALIZE AND SCALE BY ampl.
       tmp = SUM(dyh(2:ny)*(ABS(cu(2:ny,1)) + ABS(cu(2:ny,3)) &
            + half*ABS(cu(1:ny-1,2)) + half*ABS(cu(2:ny,2))))
       IF (tmp .gt. zero) cu = ampl*cu/tmp
       !! APPLY RAISED COSINE FILTER TO NOISE IN x/z-DIRECTIONS
       IF (filter .eq. 1) cu = cu*COS(half*pi*MAX(ABS(kx_me(p))/MAXVAL(kx), &
            ABS(kz(q))/MAXVAL(kz)))

    ELSE

       cu(1:ny+1,1) = ushape*CMPLX(temp(1:ny+1,1),zero)
       cu(1:ny+1,3) = ushape*CMPLX(temp(1:ny+1,2),zero)

    END IF

  END SUBROUTINE random_divfree_vector
END PROGRAM arnoldi
