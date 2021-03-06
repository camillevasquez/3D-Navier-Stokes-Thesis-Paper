MODULE pass1
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE transform
  USE pass_data
  USE solve
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nonlinear, advance_scalar

  INCLUDE 'mpif.h'

! miscellaneous loop indices.  
  INTEGER  i, j, k, x, z, nzl, mpierr, ierr

! ARRAYS FOR VELOCITY, SCALAR AND WALL-NORMAL FLUXES
  REAL(KIND=prec) :: deltax, deltaz
  REAL(KIND=prec), PARAMETER :: zero=0.0d0, half=0.5d0, one=1.0d0, two=2.0d0

CONTAINS
!=================================================================
!==============================NONLINEAR==========================
!=================================================================
  SUBROUTINE nonlinear
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz,1:3) :: u, up, yflux
    REAL(KIND=prec), DIMENSION(1:ny+1)               :: cfl

    REAL(KIND=prec) :: temp(1:nx2,1:local_nz), cfl_temp, dt_temp
    REAL(KIND=prec), DIMENSION(2) :: tmp_max, tmp2_max

    deltax = lx/DBLE(nx2)
    deltaz = lz/DBLE(nz2)

! INITIALIZE AN ABBREVIATION FOR local_nz
    nzl = local_nz

! INITIALIZE ARRAYS
    u   = zero
    up  = zero
    yflux = zero

! DOWNLOAD THE INITIAL PLANES OF THE VELOCITY AND SCALAR FIELDS.
    DO i = 1,3
       CALL xz_fft_complex_to_real(data1(1,1,1,i),u(1,1,i))
       CALL xz_fft_complex_to_real(data1(1,1,2,i),up(1,1,i))
    END DO
 
! COMPUTE WALL-NORMAL FLUXES AT BOTTOM BOUNDARY.
    yflux(1:nx2,1:nzl,1) = u(1:nx2,1:nzl,2)*u(1:nx2,1:nzl,1)
    temp = half*(u(1:nx2,1:nzl,2) + up(1:nx2,1:nzl,2))
    yflux(1:nx2,1:nzl,2) = temp*temp
    yflux(1:nx2,1:nzl,3) = u(1:nx2,1:nzl,2)*u(1:nx2,1:nzl,3)

! MOVE VARIABLES UP ONE PLANE.
    u  = up

! ZERO OUT THE RHS AT THE BOTTOM BOUNDARY FOR ALL VARIABLES.
    data1(1:nz,1:local_nx,1,1:3) = zero 

! COMPUTE THE NONLINEAR TERM AT EACH OF THE PLANES
    DO K = 2,NY-1

! LOAD THE NEW PLANE INTO UP.
       DO I = 1,3
          CALL xz_fft_complex_to_real(data1(1,1,k+1,i),up(1,1,i))
       END DO

       CALL compute_wall_parallel_flux(k)
       CALL compute_wall_normal_flux(k)

!  CALCULATE THE CFL NUMBER
       IF (((mod(step,1)==0).OR.(dt_flag==1)).and.(rk_step==1)) THEN
          CALL check_cfl(k)
       END IF

! SAVE VELOCITY AND SCALAR FIELDS BEFORE MOVING UP TO NEXT PLANE.
       u  = up

    END DO  ! END LOOP OVER PLANES

! COMPUTE THE FLUXES AND THE RESULTING NONLINEAR TERM AT y_{ny-1/2}
! (ONE HALF-PLANE BELOW THE BOUNDARY).

! DOWNLOAD VELOCITIES AT THE UPPER BOUNDARY.
    CALL xz_fft_complex_to_real(data1(1,1,ny+1,1),up(1,1,1))
    CALL xz_fft_complex_to_real(data1(1,1,ny+1,3),up(1,1,3))
    up(:,:,2) = u(:,:,2) ! u(:,:,2) ALREADY HOLDS v AT UPPER BOUNDARY

    CALL compute_wall_parallel_flux(ny)

! ZERO OUT RHS AT THE UPPER BOUNDARY.
    data1(1:nz,1:local_nx,ny,2) = zero
    DO j = 1,3
      data1(1:nz,1:local_nx,ny+1,j) = zero
    END DO
    
    !  CALCULATE THE CFL NUMBER
    IF (((mod(step,1)==0) .OR. (dt_flag == 1)) .and. (rk_step==1)) THEN
       CALL check_cfl(ny)
    END IF

    ! COMPUTE A NEW TIME STEP TO KEEP THE CFL NEAR ONE
    IF ((STEP .GE. 0) .AND. (RK_STEP .EQ. 1)) THEN
       CFL_TEMP = MAXVAL(CFL(2:ny))
       IF (nscalar .eq. 0) THEN
          !IF (me == 0) WRITE(*,*) 'MAX CFL = ', cfl_temp
       ELSE
          tmp_max(1) = MAXVAL(data_scalar)
          tmp_max(2) = -MINVAL(data_scalar)
          CALL MPI_ALLREDUCE(tmp_max,tmp2_max,2,MPI_DOUBLE_PRECISION, &
               MPI_MAX, MPI_COMM_WORLD, MPIERR)
          IF (me == 0) WRITE(*,989) cfl_temp, tmp2_max(1), -tmp2_max(2)
989       FORMAT('MAX CFL = ',E12.4,'  MAX SCALAR = ',F8.4, &
               '  MIN SCALAR = ',F8.4)
       END IF
       IF (DT_FLAG .EQ. 1) THEN
          DT_TEMP = DT*CFL_MAX/CFL_TEMP
          DT = DMIN1(DT_TEMP, DT + (DT_TEMP-DT)*DT/(ONE + DT))
       ENDIF
    ENDIF

  CONTAINS
    !=========================================================================
    SUBROUTINE compute_wall_parallel_flux(plane)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: plane
      REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz) :: fluxpu, fluxpw


      ! COMPUTE WALL-PARALLEL FLUXES AT yh(plane)

      ! WALL NORMAL FLUX OF U AND W AT y(plane).
      IF (plane < ny) THEN
         fluxpu(1:nx2,1:nzl) = u(1:nx2,1:nzl,2) &
              *half*(u(1:nx2,1:nzl,1) + up(1:nx2,1:nzl,1))
         fluxpw(1:nx2,1:nzl) = u(1:nx2,1:nzl,2) &
              *half*(u(1:nx2,1:nzl,3) + up(1:nx2,1:nzl,3))
      ELSE
         fluxpu(1:nx2,1:nzl) = u(1:nx2,1:nzl,2)*up(1:nx2,1:nzl,1)
         fluxpw(1:nx2,1:nzl) = u(1:nx2,1:nzl,2)*up(1:nx2,1:nzl,3)
      END IF

      ! TAKE WALL-NORMAL DERIVATE OF WALL-NORMAL FLUXES.  
      ! ADD BUFFER/FRINGE TERMS IF APPROPRIATE.
      yflux(1:nx2,1:nzl,1) =                                                 &
           dh(plane)*(fluxpu(1:nx2,1:nzl) - yflux(1:nx2,1:nzl,1))            &
           + buffer(1:nx2,1:nzl)*(u(1:nx2,1:nzl,1) - target_inflow(plane,1))
      yflux(1:nx2,1:nzl,3) =                                                 &
           dh(plane)*(fluxpw(1:nx2,1:nzl)-yflux(1:nx2,1:nzl,3))              &
           + buffer(1:nx2,1:nzl)*(u(1:nx2,1:nzl,3) - target_inflow(plane,3))

      ! CALL ROUTINE TO COMPUTE DIVERGENCE OF FLUXES AND PUT INTO DATA1.
      CALL compute_divergence_of_fluxes(      &
           u(1:nx2,1:nzl,1)*u(1:nx2,1:nzl,1), &  ! FLUX IN X-DIR.
           yflux(1:nx2,1:nzl,1),              &  ! D/DY of Y-FLUX
           u(1:nx2,1:nzl,3)*u(1:nx2,1:nzl,1), &  ! FLUX IN Z-DIR.
           data1(1,1,plane,1))                       ! ARRAY RETURNS RHS
      CALL compute_divergence_of_fluxes(      &
           u(1:nx2,1:nzl,1)*u(1:nx2,1:nzl,3), &  ! FLUX IN X-DIR.
           yflux(1:nx2,1:nzl,3),              &  ! D/DY of Y-FLUX
           u(1:nx2,1:nzl,3)*u(1:nx2,1:nzl,3), &  ! FLUX IN Z-DIR.
           data1(1,1,plane,3))                       ! ARRAY RETURNS RHS

      ! MOVE THE FLUX AT y(plane) INTO yflux.
      yflux(:,:,1) = fluxpu
      yflux(:,:,3) = fluxpw

    END SUBROUTINE compute_wall_parallel_flux
    !=========================================================================
    SUBROUTINE compute_wall_normal_flux(plane)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: plane
      REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz) :: fluxp, temp

      ! COMPUTE THE FLUXES FOR THE WALL-NORMAL VELOCITY AT y_k

      ! COMPUTE WALL-NORMAL FLUXES AT y(plane).
      temp = HALF*(u(1:nx2,1:nzl,2) + up(1:nx2,1:nzl,2))
      fluxp = temp*temp

      ! TAKE THE DIVERGENCE OF FLUXES AND PUT INTO DATA1.
      yflux(1:nx2,1:nzl,2)= dn(plane)*(fluxp(:,:) - yflux(:,:,2))          &
           + BUFFER(1:nx2,1:nzl)*(u(1:nx2,1:nzl,2) - TARGET_INFLOW(plane,2))
      CALL compute_divergence_of_fluxes(u(1:nx2,1:nzl,2)           &
           *(IU(1,plane)*u(1:nx2,1:nzl,1)+ IU(2,plane)*up(1:nx2,1:nzl,1)), &
           yflux(1:nx2,1:nzl,2),                                   &
           u(1:nx2,1:nzl,2)                                        &
           *(IU(1,plane)*u(1:nx2,1:nzl,3)+ IU(2,plane)*up(1:nx2,1:nzl,3)), &
           data1(1,1,plane,2))

      ! MOVE THE FLUX AT y(plane) INTO yflux.
      yflux(:,:,2) = fluxp

    END SUBROUTINE compute_wall_normal_flux
    !=======================================================================
    SUBROUTINE check_cfl(plane)
      IMPLICIT NONE

      INTEGER plane
      REAL(KIND=prec) TEMPA, TEMPB, TEMPC, big_cfl, poss_max
      REAL(KIND=prec), PARAMETER :: CFL_FOURIER =0.55d0,  CFL_FD =1.73d0, &
           DIFF_FOURIER=0.25d0,  DIFF_FD=2.51d0
      REAL(KIND=prec), PARAMETER :: zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0

      tempa = half*(MAXVAL(ABS(u(:,:,1)+up(:,:,1))))*DBLE(nx)/(cfl_fourier*lx)
      tempb = MAXVAL(ABS(u(:,:,2)))/(cfl_fd*min(hn(plane),hn(plane-1)))
      tempc = half*(MAXVAL(ABS(u(:,:,3)+up(:,:,3))))*DBLE(nz)/(cfl_fourier*lz)
      poss_max = tempa + tempb + tempc
      CALL MPI_ALLREDUCE(poss_max,big_cfl,1,MPI_DOUBLE_PRECISION, &
           MPI_MAX, MPI_COMM_WORLD, MPIERR)
      cfl(plane) = big_cfl*dt
      !         WRITE(*,*) STEP, RK_STEP, plane, CFL(plane)

!!$      IF (plane == 2) poss_max = zero
!!$      CALL MPI_ALLREDUCE(MAX(poss_max,MAXVAL(abs(u(:,:,1:3))),&
!!$           MAXVAL(abs(up(:,:,1:3)))),poss_max,&
!!$           1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,MPIERR)
!!$      IF ((plane==ny-1).and.(me==0)) WRITE(*,999) me,STEP,t_total,poss_max
!!$      999 FORMAT(i4,i6,2f12.4)

      IF (DT_FLAG .EQ. 1) THEN
         IF (CFL(plane) .GE. 3.0D+2) THEN             !DT IS TOO BIG  ... ABORT
            WRITE(*, 1000) STEP, RK_STEP, plane, CFL(plane)
            WRITE(*, 1010)
            CALL finalize_mpi
            STOP 'DT IS TOO BIG'
         ENDIF
      ELSE
         IF (CFL(plane) .GE. 1.15D0) THEN            !CFL IS TOO BIG ... ABORT
            WRITE(*,*) MAXVAL(ABS(u(:,:,1))), MAXVAL(ABS(u(:,:,2))), &
                 MAXVAL(ABS(u(:,:,3)))
            WRITE(*,*) MAXVAL(ABS(up(:,:,1))), MAXVAL(ABS(u(:,:,2))), &
                 MAXVAL(ABS(up(:,:,3)))
            WRITE(*, 1000) STEP, RK_STEP, plane, CFL(plane)
            WRITE(*, 1020)
            CALL finalize_mpi
            STOP 'CFL IS TOO BIG'
         ENDIF
      END IF
1000  FORMAT(' ',5X,'STEP = ',I6,' RK_STEP = ',I4,' plane = ',I4,&
           ' CFL = ',1PE15.4)
1010  FORMAT(' ','ABORTING ... TOO RAPID CHANGE OF DT REQUIRED')
1020  FORMAT(' ','ABORTING ... CFL IS TOO BIG')
    END SUBROUTINE check_cfl
! END OF nonlinear
  END SUBROUTINE nonlinear
  !=========================================================================
  SUBROUTINE advance_scalar(maxm)
    IMPLICIT NONE
    INTEGER maxm, mwork, maux, method(7), mthlim(1)

    REAL(KIND=prec), DIMENSION(-1:maxm+2,3,3) :: aux1, aux2, aux3
    REAL(KIND=prec), DIMENSION(-1:maxm+2,1)   :: qadd, fadd, q1d
    REAL(KIND=prec), DIMENSION(-1:maxm+2,1,2,-1:1) :: gadd, hadd
    REAL(KIND=prec), DIMENSION(-1:maxm+2)     :: dtdx1d, dtdy1d, dtdz1d
    REAL(KIND=prec), &
         DIMENSION(2*(nx+4)*(local_nz_small+4)*(ny+3) + 86*(maxm+4)) :: work

!!$    REAL(KIND=prec), DIMENSION(:,:,:), ALLOCATABLE:: aux1,aux2,aux3
!!$    REAL(KIND=prec), DIMENSION(:,:,:,:), ALLOCATABLE:: gadd,hadd
!!$    REAL(KIND=prec), DIMENSION(:,:), ALLOCATABLE  :: qadd, fadd, q1d
!!$    REAL(KIND=prec), DIMENSION(:), ALLOCATABLE :: dtdx1d,dtdy1d,dtdz1d,work

    REAL(KIND=prec), DIMENSION(-1:nx+2,-1:local_nz_small+2,-1:ny+1) :: qin,qout
    REAL(KIND=prec), DIMENSION(-1:nx+2,-1:local_nz_small+2,-1:ny+1,4):: fluxvel
    REAL(KIND=prec), DIMENSION(1:nx,1:local_nz_small) :: factor, rtmp
    REAL(KIND=prec), DIMENSION(0:nx+1,0:local_nz_small+1) :: rpad
    REAL(KIND=prec) :: tmpcfl, tmpcfl2(4)
    EXTERNAL rpn3, rpt3, rptt3

    deltax = lx/DBLE(nx)
    deltaz = lz/DBLE(nz)
    nzl = local_nz_small

    !! TRANSFORM CURRENT VELOCITY FIELD INTO fluxvel ARRAY.
    DO k = 1,ny-1
       CALL small_fft_complex_to_real_pad(data1(1,1,k+1,1),rpad)
       fluxvel(1:nx,1:nzl,k,1) = deltaz*hh(k+1) &
            *half*(rpad(0:nx-1,1:nzl)+rpad(1:nx,1:nzl))
    END DO
    CALL pass_slabs(nx,nzl,ny-1,2,fluxvel(-1,-1,1,1))

    DO k = 1,ny-1
       CALL small_fft_complex_to_real_pad(data1(1,1,k+1,3),rpad)
       fluxvel(1:nx,1:nzl,k,2) = deltax*hh(k+1) &
            *half*(rpad(1:nx,0:nzl-1)+rpad(1:nx,1:nzl))
    END DO
    CALL pass_slabs(nx,nzl,ny-1,2,fluxvel(-1,-1,1,2))

    DO k = 1,ny
       CALL small_fft_complex_to_real(data1(1,1,k,2),rtmp)
       fluxvel(1:nx,1:nzl,k,3) = rtmp(1:nx,1:nzl)*deltax*deltaz
    END DO
    CALL pass_slabs(nx,nzl,ny,2,fluxvel(-1,-1,1,3))

    ! SET UP CAPACITY MATRIX -- ACCOUNTS FOR VARIABLE GRID SPACING IN y-DIR.
    DO k = 2,ny
       fluxvel(-1:nx+2,-1:nzl+2,k-1,4) = deltax*hh(k)*deltaz
    END DO

    !! SET COEFFICIENT FOR DAMPING OF SCALAR IN BUFFER LAYER.
    factor = EXP(-half*dt*buffer_small)

!!$    maxm = MAX(nx,nzl,ny-1)
    mwork = 2*(nx+4)*(local_nz_small+4)*(ny+3) + 86*(maxm+4)
!!$    ALLOCATE(aux1(-1:maxm+2,3,3),aux2(-1:maxm+2,3,3),aux3(-1:maxm+2,3,3), &
!!$         qadd(-1:maxm+2,1), fadd(-1:maxm+2,1), q1d(-1:maxm+2,1), &
!!$         gadd(-1:maxm+2,1,2,-1:1), hadd(-1:maxm+2,1,2,-1:1), &
!!$         dtdx1d(-1:maxm+2), dtdy1d(-1:maxm+2), dtdz1d(-1:maxm+2), &
!!$         work(mwork), STAT=ierr)
!!$    IF (ierr .ne. 0) STOP 'CANNOT DECLARE ARRAYS IN CLAWPACK_CALL'

    DO i = 1,nscalar
       qin = zero
       IF (pr(i) .ne. zero) THEN
          DO k = 2,ny
             qin(1:nx,1:nzl,k-1) = data_scalar(1:nx,1:nzl,k,i)
          END DO
          CALL pass_slabs(nx,nzl,ny+3,2,qin(-1,-1,-1))
          CALL diffuse_scalar(dt,pr(i)*re,qin,qout)
       ELSE
          DO k = 2,ny
             qin(1:nx,1:nzl,k-1) = factor*data_scalar(1:nx,1:nzl,k,i) &
                  + (one-factor)*target_inflow(k,4)
          END DO
       END IF

       ! SET UP BOUNDARY CONDITIONS AT TOP WALLS.
       IF (bl_grid .gt. 0) THEN
          !! OUTFLOW BOUNDARY CONDITIONS AT TOP WALL.
          !! USE ZEROTH ORDER EXTRAPOLATION.
          qin(1:nx,1:nzl,ny)   = qin(1:nx,1:nzl,ny-1)
          qin(1:nx,1:nzl,ny+1) = qin(1:nx,1:nzl,ny-1)
          fluxvel(:,:,ny,1)   = fluxvel(:,:,ny-1,1) ! STREAMWISE VELOCITY
          fluxvel(:,:,ny+1,1) = fluxvel(:,:,ny-1,1)
          fluxvel(:,:,ny,2)   = fluxvel(:,:,ny-1,2) ! SPANWISE VELOCITY
          fluxvel(:,:,ny+1,2) = fluxvel(:,:,ny-1,2)
          fluxvel(:,:,ny+1,3) = fluxvel(:,:,ny,3)   ! WALL-NORMAL VELOCITY
          fluxvel(:,:,ny,4)   = fluxvel(:,:,ny-1,4) ! CAPACITY FUNCTION
          fluxvel(:,:,ny+1,4) = fluxvel(:,:,ny-1,4)
       ELSEIF (jetmag_top .eq. zero) THEN
          !! SOLID TOP WALL.
          qin(1:nx,1:nzl,ny)   = qin(1:nx,1:nzl,ny-1)
          qin(1:nx,1:nzl,ny+1) = qin(1:nx,1:nzl,ny-2)
          fluxvel(:,:,ny,1)   = fluxvel(:,:,ny-1,1) ! STREAMWISE VELOCITY
          fluxvel(:,:,ny+1,1) = fluxvel(:,:,ny-2,1)
          fluxvel(:,:,ny,2)   = fluxvel(:,:,ny-1,2) ! SPANWISE VELOCITY
          fluxvel(:,:,ny+1,2) = fluxvel(:,:,ny-2,2)
          fluxvel(:,:,ny+1,3) = - fluxvel(:,:,ny-1,3) ! WALL-NORMAL VEL.
          fluxvel(:,:,ny,4)   = fluxvel(:,:,ny-1,4) ! CAPACITY FUNCTION
          fluxvel(:,:,ny+1,4) = fluxvel(:,:,ny-2,4)
       ELSE
          !! TOP WALL INFLOW BC'S.
          qin(1:nx,1:nzl,ny)   = bc_top(1:nx,1:local_nz_small,4)
          qin(1:nx,1:nzl,ny+1) = bc_top(1:nx,1:local_nz_small,4)
          fluxvel(:,:,ny,1)   = fluxvel(:,:,ny-1,1) ! STREAMWISE VELOCITY
          fluxvel(:,:,ny+1,1) = fluxvel(:,:,ny-1,1)
          fluxvel(:,:,ny,2)   = fluxvel(:,:,ny-1,2) ! SPANWISE VELOCITY
          fluxvel(:,:,ny+1,2) = fluxvel(:,:,ny-1,2)
          fluxvel(:,:,ny+1,3) = fluxvel(:,:,ny,3) ! WALL-NORMAL VEL.
          ! INCLUDE TIME-VARIATION LATER
          fluxvel(:,:,ny,4)   = fluxvel(:,:,ny-1,4) ! CAPACITY FUNCTION
          fluxvel(:,:,ny+1,4) = fluxvel(:,:,ny-2,4)
       END IF

       !! SET UP BOTTOM WALL BOUNDARY CONDITIONS.
       IF (jetmag_bot .eq. zero) THEN
          !! SOLID BOTTOM WALL.
          qin(1:nx,1:nzl,0)  = qin(1:nx,1:nzl,1)
          qin(1:nx,1:nzl,-1) = qin(1:nx,1:nzl,2)
          fluxvel(:,:, 0,1) = fluxvel(:,:,1,1) ! STREAMWISE VELOCITY
          fluxvel(:,:,-1,1) = fluxvel(:,:,2,1)
          fluxvel(:,:, 0,2) = fluxvel(:,:,1,2) ! SPANWISE VELOCITY
          fluxvel(:,:,-1,2) = fluxvel(:,:,2,2)
          fluxvel(:,:, 0,3) = - fluxvel(:,:,1,4) ! WALL-NORMAL VEL.
          fluxvel(:,:,-1,3) = - fluxvel(:,:,2,3) 
          fluxvel(:,:, 0,4) = fluxvel(:,:,1,4) ! CAPACITY FUNCTION
          fluxvel(:,:,-1,4) = fluxvel(:,:,2,4)
       ELSE
          !! BOTTOM WALL INFLOW BC'S.
          qin(1:nx,1:nzl,0)  = bc_bot(1:nx,1:local_nz_small,4)
          qin(1:nx,1:nzl,-1) = bc_bot(1:nx,1:local_nz_small,4)
          fluxvel(:,:, 0,1) = fluxvel(:,:,1,1) ! STREAMWISE VELOCITY
          fluxvel(:,:,-1,1) = fluxvel(:,:,2,1)
          fluxvel(:,:, 0,2) = fluxvel(:,:,1,2) ! SPANWISE VELOCITY
          fluxvel(:,:,-1,2) = fluxvel(:,:,2,2)
          fluxvel(:,:, 0,3) = fluxvel(:,:,1,3) ! WALL-NORMAL VEL.
          fluxvel(:,:,-1,3) = fluxvel(:,:,1,3) ! INCLUDE TIME-VARIATION LATER
          fluxvel(:,:, 0,4) = fluxvel(:,:,1,4) ! CAPACITY FUNCTION
          fluxvel(:,:,-1,4) = fluxvel(:,:,2,4)
       END IF
       CALL pass_slabs(nx,nzl,ny+3,2,qin(-1,-1,-1)) ! PASS INFO TO NEIGHBORS

       qout = qin  !! SET qout EQUAL TO qin.

       method(1) = 0  !! FIXED TIME STEP
       method(2) = 2  !! USE SECOND ORDER CORRECTIONS.
       method(4) = 0  !! SUPPRESS PRINTING OF dt AND CFL NUMBER.
       method(5) = 0  !! NO SOURCE TERM.
       method(6) = 4  !! USE CAPACITY FUNCTION TO ACCOUNT FOR STRETCHED GRID
       method(7) = 3  !! THREE COMPONENTS OF VELOCITY DRIVE ADVECTION.

       mthlim(1) = 4  !! USE MC LIMITER (MONOTONIZED CENTERED-DIFFERENCE)

       maux = 3
       tmpcfl2 = zero
       tmpcfl2(1) = one

       IF (split_scalar .eq. 1) THEN
          method(3) = -1 !! GODUNOV SPLITTING
          CALL dimsp3(maxm,nx,nzl,ny-1,1,1,2,nx,nzl,ny-1,qin,qout, &
               fluxvel,one,one,one,dt,method,mthlim,tmpcfl,tmpcfl2, &
               qadd, fadd, gadd, hadd, q1d,   dtdx1d, dtdy1d, dtdz1d,      &
               aux1, aux2, aux3, maux, work, mwork, rpn3,   rpt3,   rptt3)
       ELSE
          method(3) = 22 !! 3D PROPAGATION OF INCREMENT AND CORRECTION WAVES.
          CALL step3(maxm,nx,nzl,ny-1,1,1,2,nx,nzl,ny-1,qin,qout, &
               fluxvel,one,one,one,dt,method,mthlim,tmpcfl, &
               qadd, fadd, gadd, hadd, q1d,   dtdx1d, dtdy1d, dtdz1d,      &
               aux1, aux2, aux3, maux, work, mwork, rpn3,   rpt3,   rptt3)
       END IF

       IF (pr(i) .ne. zero) THEN
          DO k = 1,ny+1
             data_scalar(1:nx,1:nzl,k,i) = qout(1:nx,1:nzl,k-1)
          END DO
       ELSE
          DO k = 1,ny+1
             data_scalar(1:nx,1:nzl,k,i) = factor*qout(1:nx,1:nzl,k-1) &
                  + (one-factor)*target_inflow(k,4)
          END DO
       END IF

    END DO

!!$    DEALLOCATE(aux1,aux2,aux3,qadd, fadd, q1d, gadd, hadd, &
!!$         dtdx1d, dtdy1d, dtdz1d, work, STAT=ierr)
!!$    IF (ierr .ne. 0) STOP 'CANNOT DE-ALLOCATE ARRAYS IN CLAWPACK_CALL'

  CONTAINS
    !=======================================================================
    SUBROUTINE diffuse_scalar(deltat,pe,qin,qout)
      IMPLICIT NONE
      REAL(KIND=prec) :: deltat, pe, diff, gamma, delta, factor, dx2, dz2
      REAL(KIND=prec), DIMENSION(-1:nx+2,-1:local_nz_small+2,-1:ny+1)::qin,qout
      REAL(KIND=prec), DIMENSION(0:nx+1,0:local_nz_small+1) :: &
           damping, topbc, botbc
      REAL(KIND=prec) :: factorx, factorz, factor2x, factor2z, factor2y, &
           rhs(0:ny,1), temp(1:nx,1:local_nz_small)

      nzl = local_nz_small

      gamma = one - one/SQRT(two)
      delta = - one/SQRT(two)
      damping(1:nx,1:nzl) = buffer_small
      CALL pass_neighbors(nx,nzl,1,damping)
      topbc(1:nx,1:nzl) = bc_top(1:nx,1:nzl,4)
      CALL pass_neighbors(nx,nzl,1,topbc)
      botbc(1:nx,1:nzl) = bc_bot(1:nx,1:nzl,4)
      CALL pass_neighbors(nx,nzl,1,botbc)

      !! FIRST SUBSTEP OF TWO-STEP RUNGE-KUTTA SCHEME.
      !! SEE ASCHER, RUUTH & SPITERI, APPLIED NUMERICAL MATHEMATICS
      !! VOLUME 25, pp. 151--167, 1997.  THE FOLLOWING SCHEME IS
      !! A SECOND-ORDER, TWO-STAGE EXPLICIT/IMPLICIT RUNGE-KUTTA SCHEME.
      !! THE VERTICAL DIFFUSIVITY AND BUFFER TERMS ARE TREATED IMPLICITLY
      !! AND THE HORIZONTAL DIFFUSIVITY IS TREATED EXPLICITLY.
      !! EXPLICIT PIECE OF FIRST STAGE.
      factorx = deltat*gamma/pe*(DBLE(nx)/lx)**2
      factorz = deltat*gamma/pe*(DBLE(nz)/lz)**2
      DO k = 1,ny-1
         DO z = 0,nzl+1
            DO x = 0,nx+1
               qout(x,z,k) = (one - two*factorx - two*factorz)*qin(x,z,k) &
                    + factorx*qin(x+1,z,k) + factorx*qin(x-1,z,k) &
                    + factorz*qin(x,z+1,k) + factorz*qin(x,z-1,k) &
                    + deltat*gamma*damping(x,z)*target_inflow(k+1,4)
            END DO
         END DO
      END DO

     !! IMPLICIT PIECE OF FIRST STAGE.
      factor = deltat*gamma/pe
      DO z = 0,nzl+1
         DO x = 0,nx+1
            rhs(0,1) = zero
            rhs(1:ny-1,1) = qout(x,z,1:ny-1)
            rhs(ny,1) = zero
            CALL solve_real_helmholtz(ny+1, 1, rhs(0,1), &
                 one + deltat*gamma*damping(x,z),factor,d2h, &
                 topbc(x,z),botbc(x,z),bl_grid,0)
            qout(x,z,0:ny) = rhs(:,1)
         END DO
      END DO

      !! EXPLICIT PIECE OF SECOND STAGE.
      !! TO THE OLD VALUE OF q (BEFORE FIRST STAGE), 
      !! ADD delta TIMES HORIZONTAL DIFFUSION TERM WITH q_n,
      factorx = delta*deltat*(DBLE(nx)/lx)**2/pe
      factorz = delta*deltat*(DBLE(nz)/lz)**2/pe
      !! ADD (one-delta) TIMES HORIZONTAL DIFFUSION TERM WITH INTERMEDIATE q,
      factor2x = (one-delta)*deltat*(DBLE(nx)/lx)**2/pe
      factor2z = (one-delta)*deltat*(DBLE(nz)/lz)**2/pe
      !! AND ADD (one-gamma) TIMES IMPLICIT TERMS WITH INTERMEDIATE q.
      factor2y  = (one-gamma)*deltat
      diff = one/pe
      DO k = 1,ny-1
         DO z = 1,nzl
            DO x = 1,nx
               temp(x,z) = (one - two*factorx - two*factorz)*qin(x,z,k) &
                    + factorx*qin(x+1,z,k) + factorx*qin(x-1,z,k) &
                    + factorz*qin(x,z+1,k) + factorz*qin(x,z-1,k) &
                    - (two*factor2x + two*factor2z)*qout(x,z,k) &
                    + factor2x*qout(x+1,z,k) + factor2x*qout(x-1,z,k) &
                    + factor2z*qout(x,z+1,k) + factor2z*qout(x,z-1,k) &
                    + factor2y*(diff*d2h(2,k+1) - damping(x,z))*qout(x,z,k) &
                    + factor2y*diff*d2h(1,k+1)*qout(x,z,k-1) &
                    + factor2y*diff*d2h(3,k+1)*qout(x,z,k+1) &
                    + deltat*damping(x,z)*target_inflow(k+1,4)
            END DO
         END DO
         qin(1:nx,1:nzl,k) = temp
      END DO

      !! IMPLICIT PIECE OF SECOND STAGE.
      factor = deltat*gamma/pe
      DO z = 1,nzl
         DO x = 1,nx
            rhs(0,1) = zero
            rhs(1:ny-1,1) = qin(x,z,1:ny-1)
            rhs(ny,1) = zero
            CALL solve_real_helmholtz(ny+1, 1, rhs(0,1), &
                 one + deltat*gamma*damping(x,z),factor,d2h, &
                 topbc(x,z),botbc(x,z),bl_grid,0)
            qin(x,z,0:ny) = rhs(:,1)
         END DO
      END DO

    END SUBROUTINE diffuse_scalar
  END SUBROUTINE advance_scalar
END MODULE pass1

!!$    !=======================================================================
!!$    SUBROUTINE compute_momentum_flux(plane)
!!$      IMPLICIT NONE
!!$
!!$      INTEGER, INTENT(IN) :: plane
!!$      INTEGER :: x,z
!!$      REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz) :: fluxpu, fluxpv, fluxpw
!!$      REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz,3) :: xflux, zflux
!!$      REAL(KIND=prec) :: temp
!!$
!!$      !! ASSUME plane < ny --- CALL compute_wall_parallel_flux for plane==ny
!!$
!!$      DO z = 1,nzl
!!$         DO x = 1,nx2
!!$            ! WALL NORMAL FLUX OF U AND W AT y(plane).
!!$            fluxpu(x,z) = u(x,z,2)*half*(u(x,z,1) + up(x,z,1))
!!$            fluxpw(x,z) = u(x,z,2)*half*(u(x,z,3) + up(x,z,3))
!!$            ! COMPUTE WALL-NORMAL FLUXES AT y(plane).
!!$            temp = HALF*(u(x,z,2) + up(x,z,2))
!!$            fluxpv(x,z) = temp*temp
!!$            
!!$            ! TAKE WALL-NORMAL DERIVATE OF WALL-NORMAL FLUXES.  
!!$            ! ADD BUFFER/FRINGE TERMS IF APPROPRIATE.
!!$            yflux(x,z,1) = dh(plane)*(fluxpu(x,z) - yflux(x,z,1))  &
!!$                 + buffer(x,z)*(u(x,z,1) - target_inflow(plane,1))
!!$            yflux(x,z,2) = dn(plane)*(fluxpv(x,z) - yflux(x,z,2))  &
!!$                 + buffer(x,z)*(u(x,z,2) - target_inflow(plane,2))
!!$            yflux(x,z,3) = dh(plane)*(fluxpw(x,z)-yflux(x,z,3))    &
!!$                 + buffer(x,z)*(u(x,z,3) - target_inflow(plane,3))
!!$
!!$            ! COMPUTE STREAMWISE FLUXES OF MOMENTUM
!!$            xflux(x,z,1) = u(x,z,1)*u(x,z,1)
!!$            xflux(x,z,2) = u(x,z,2)*(iu(1,plane)*u(x,z,1) &
!!$                                   + iu(2,plane)*up(x,z,1))
!!$            xflux(x,z,3) = u(x,z,3)*u(x,z,1)
!!$
!!$            ! COMPUTE SPANWISE FLUXES OF MOMENTUM
!!$            zflux(x,z,1) = u(x,z,1)*u(x,z,3)
!!$            zflux(x,z,2) = u(x,z,2)*(iu(1,plane)*u(x,z,3) &
!!$                                   + iu(2,plane)*up(x,z,3))
!!$            zflux(x,z,3) = u(x,z,3)*u(x,z,3)
!!$
!!$         END DO
!!$      END DO
!!$
!!$      ! CALL ROUTINE TO COMPUTE DIVERGENCE OF FLUXES AND PUT INTO DATA1.
!!$      CALL compute_divergence_of_fluxes( &
!!$           xflux(1,1,1),yflux(1,1,1),zflux(1,1,1),data1(1,1,plane,1))
!!$      CALL compute_divergence_of_fluxes( &
!!$           xflux(1,1,2),yflux(1,1,2),zflux(1,1,2),data1(1,1,plane,2))
!!$      CALL compute_divergence_of_fluxes(      &
!!$           xflux(1,1,3),yflux(1,1,3),zflux(1,1,3),data1(1,1,plane,3))
!!$
!!$      ! MOVE THE FLUX AT y(plane) INTO yflux.
!!$      yflux(:,:,1) = fluxpu
!!$      yflux(:,:,2) = fluxpv
!!$      yflux(:,:,3) = fluxpw
!!$
!!$    END SUBROUTINE compute_momentum_flux
