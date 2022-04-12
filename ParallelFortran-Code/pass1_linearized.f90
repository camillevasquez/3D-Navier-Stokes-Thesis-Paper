MODULE pass1_linearized
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE transform
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nonlinear_linearized, init_baseflow

  INCLUDE 'mpif.h'

! miscellaneous loop indices.  
  INTEGER  i, j, k, nzl, mpierr, ierr

CONTAINS
!=================================================================
!=========================INIT_BASEFLOW===========================
!=================================================================
  SUBROUTINE init_baseflow
    IMPLICIT NONE

!    INCLUDE 'hdf.f90'
    
    REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:) :: data
    REAL(KIND=prec) :: temp
    REAL(KIND=prec), DIMENSION(nx,local_nz_small) :: u1
    COMPLEX(KIND=prec), DIMENSION(nz,local_nx)    :: cu
    REAL(KIND=prec), PARAMETER :: zero=0.0, half=0.5, one=1.0, two=2.0
    INTEGER, DIMENSION(0:2) :: vel_length, vel_start, vel_edges, vel_stride
    INTEGER status, sds_id, sd_id, sfstart, sfselect, sfrdata, sfendacc, sfend
    CHARACTER(LEN=8) fname

    ALLOCATE(baseflow(nx2,local_nz,ny+1,3+nscalar),STAT=ierr)
    IF (ierr /= 0) STOP 'COULD NOT ALLOCATE ARRAY FOR BASEFLOW'

    fname = 'mean.hdf'
!    sd_id = sfstart(FNAME,DFACC_READ)
    IF (sd_id .eq. -1) THEN
       baseflow = zero
       DO k = 1,ny+1
          baseflow(1:nx2,1:local_nz,k,1) = (one - yh(k)**2)*cos(beta)
          baseflow(1:nx2,1:local_nz,k,3) = (one - yh(k)**2)*sin(beta)
       END DO
       ! SET TIME STEP
       dt = cfl_max*delta1/32.d0
    ELSE
       IF (me == 0) WRITE(*,*) 'READ IN BASEFLOW FROM mean.hdf'
       ! READ IN BASEFLOW VELOCITY FIELD.
       ALLOCATE(data(nx,nz,ny+1),STAT=ierr)
       DO i = 1,3
          IF (me == 0) THEN
             ! DESCRIBE VELOCITY FIELD
             vel_length(0) = nx
             vel_length(1) = nz
             vel_length(2) = ny+1
             vel_start(0:2) = 0
             vel_stride(0:2) = 1
             vel_edges(0:2) = vel_length(0:2)
             ! INITIALIZE ARRAY AND READ IN CURRENT COMPONENT OF VELOCITY.
             data = zero
 !            sds_id = sfselect(sd_id,i)
 !            status = sfrdata(sds_id,vel_start,vel_stride,vel_edges,data)
             IF (status .ne. 0) STOP 'CANNOT READ VELOCITY IN mean.hdf'
 !            status = sfendacc(sds_id)
             IF (i == 2) THEN
                ! SET DT ACCORDING TO CFL CONDITION AT THE BOTTOM WALL.
                temp = max(MAXVAL(ABS(data(:,:,1))), &
                     MAXVAL(ABS(data(:,:,ny))))
                WRITE(*,*) temp
                IF (temp > 0.1d0) THEN
                   dt = 0.75d0*1.73d0*(yh(2) - yh(1))/temp
                ELSE
                   dt = cfl_max*delta1/32.d0
                END IF
             END IF
          END IF
          DO k = 1,ny+1
             !! SCATTER VELOCITY FIELD, WHICH IS DEFINED ON SMALL NXxNZ GRID
             CALL MPI_SCATTERV( &
                  data(1,1,k),nz_small_proc*nx,z_small_start*nx, &
                  MPI_DOUBLE_PRECISION,u1,nz_small_proc(me)*nx, &
                  MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
             !! TRANSFORM INTO FOURIER SPACE
             CALL small_fft_real_to_complex(u1,cu)
             !! TRANSFORM BACK INTO PHYSICAL SPACE ON DE-ALIASED GRID
             CALL xz_fft_complex_to_real(cu,baseflow(1,1,k,i))
          END DO
       END DO
       DEALLOCATE(data)
       IF (nscalar .gt. 0) THEN
          ! READ IN BASEFLOW SCALAR FIELDS
          ALLOCATE(data(nx2,nz2,ny+1),STAT=ierr)
          DO i = 1,nscalar
             IF (me == 0) THEN
                ! DESCRIBE SCALAR FIELD
                vel_length(0) = nx2
                vel_length(1) = nz2
                vel_length(2) = ny+1
                vel_start(0:2) = 0
                vel_stride(0:2) = 1
                vel_edges(0:2) = vel_length(0:2)
                ! INITIALIZE ARRAY AND READ IN CURRENT SCALAR FIELD
                data = zero
 !               sds_id = sfselect(sd_id,3+i)
 !               status = sfrdata(sds_id,vel_start,vel_stride,vel_edges,data)
                IF (status .ne. 0) STOP 'CANNOT READ SCALAR IN mean.hdf'
 !               status = sfendacc(sds_id)
             END IF
             DO k = 1,ny+1
                !! SCATTER SCALAR FIELD, WHICH IS DEFINED ON LARGE NX2xNZ2 GRID
                CALL MPI_SCATTERV( &
                     data(1,1,k),nz_proc*nx2,z_start*nx2,MPI_DOUBLE_PRECISION,&
                     baseflow(1,1,k,3+i),nz_proc(me)*nx2,MPI_DOUBLE_PRECISION,&
                     0,MPI_COMM_WORLD,mpierr)
             END DO
          END DO
          DEALLOCATE(data)
       END IF
!       status = sfend(sd_id)
       CALL MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
    END IF

! ZERO OUT THE PRESSURE GRADIENT AND TARGET INFLOW PROFILES SINCE THIS
! SIMULATION IS LINEARIZED ABOUT AN EXISTING BASEFLOW.
    p_grad = zero
    
    IF (me == 0) WRITE(*,*) 'LINEARIZED SIMULATION'
    IF (me == 0) WRITE(*,*) 'CONSTANT DT = ', dt

  END SUBROUTINE init_baseflow
!=================================================================
!=======================NONLINEAR_LINEARIZED======================
!=================================================================
  SUBROUTINE nonlinear_linearized
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz,1:3+nscalar) :: u, up, yflux
    REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz,1:nscalar)   :: tm, tp2
    REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz) :: temp, temp_baseflow
    REAL(KIND=prec)  :: cfl(1:ny+1), CFL_TEMP, DT_TEMP
    REAL(KIND=prec), PARAMETER :: zero=0.0, half=0.5, one=1.0, two=2.0

! INITIALIZE AN ABBREVIATION FOR local_nz
    nzl = local_nz

! DOWNLOAD THE INITIAL PLANES OF THE VELOCITY AND SCALAR FIELDS.
    DO i = 1,3+nscalar
       CALL xz_fft_complex_to_real(data1(1,1,1,i),u(1,1,i))
       CALL xz_fft_complex_to_real(data1(1,1,2,i),up(1,1,i))
    END DO
 
    IF (nscalar > 0) THEN
       DO i = 1,nscalar
          CALL xz_fft_complex_to_real(data1(1,1,3,3+i),tp2(1,1,i))
       END DO
    END IF

! COMPUTE WALL-NORMAL FLUXES AT BOTTOM BOUNDARY.
    yflux(1:nx2,1:nzl,1) = &
           baseflow(1:nx2,1:nzl,1,2)*u(1:nx2,1:nzl,1) &
         + u(1:nx2,1:nzl,2)*baseflow(1:nx2,1:nzl,1,1)
    temp = half*(u(1:nx2,1:nzl,2) + up(1:nx2,1:nzl,2))
    temp_baseflow = half*(baseflow(1:nx2,1:nzl,1,2) &
                        + baseflow(1:nx2,1:nzl,2,2))
    yflux(1:nx2,1:nzl,2) = two*temp*temp_baseflow
    DO i = 3,3+nscalar
       yflux(1:nx2,1:nzl,i) = &
              baseflow(1:nx2,1:nzl,1,2)*u(1:nx2,1:nzl,i) &
            + u(1:nx2,1:nzl,2)*baseflow(1:nx2,1:nzl,1,i)
    END DO

! MOVE VARIABLES UP ONE PLANE.
    u  = up
    IF (nscalar > 0) tm(:,:,1:nscalar) = u(:,:,4:3+nscalar)

! ZERO OUT THE RHS AT THE BOTTOM BOUNDARY FOR ALL VARIABLES.
    DATA1(1:nz,1:local_nx,1,1:3+nscalar) = ZERO

! COMPUTE THE LINEARIZED NONLINEAR TERM AT EACH OF THE PLANES
    DO K = 2,NY-1

! LOAD THE NEW PLANES INTO UP AND TP2.
       DO I = 1,3
          CALL xz_fft_complex_to_real(data1(1,1,k+1,i),up(1,1,i))
       END DO

       IF (nscalar > 0) THEN
          up(1:nx2,1:nzl,4:3+nscalar) = tp2(1:nx2,1:nzl,1:nscalar)
          DO i = 1,nscalar
             CALL xz_fft_complex_to_real(data1(1,1,k+2,3+i),tp2(1,1,i))
          END DO
       END IF

       CALL compute_wall_parallel_flux(k)
       CALL compute_wall_normal_flux(k)
       IF (nscalar > 0) CALL compute_scalar_flux(k)
       
!  CALCULATE THE CFL NUMBER
       IF (((mod(step,1)==0) .OR. (dt_flag == 1))  &
            .and. (rk_step==1)) THEN
          CALL check_cfl(k)
       END IF

! SAVE VELOCITY AND SCALAR FIELDS BEFORE MOVING UP TO NEXT PLANE.
       u  = up
       IF (nscalar > 0) tm(:,:,1:nscalar) = u(:,:,4:3+nscalar)

    END DO  ! END LOOP OVER PLANES

! COMPUTE THE FLUXES AND THE RESULTING NONLINEAR TERM AT y_{ny-1/2}
! (ONE HALF-PLANE BELOW THE BOUNDARY).

! DOWNLOAD VELOCITIES AT THE UPPER BOUNDARY.
    CALL xz_fft_complex_to_real(data1(1,1,ny+1,1),up(1,1,1))
    CALL xz_fft_complex_to_real(data1(1,1,ny+1,3),up(1,1,3))
    up(:,:,2) = u(:,:,2) ! u(:,:,2) ALREADY HOLDS v AT UPPER BOUNDARY

    IF (nscalar > 0) THEN
       up(1:nx2,1:nzl,4:3+nscalar) = tp2(1:nx2,1:nzl,1:nscalar)
    END IF

    CALL compute_wall_parallel_flux(ny)
    IF (nscalar > 0) CALL compute_scalar_flux(ny)

! ZERO OUT RHS AT THE UPPER BOUNDARY.
    data1(1:nz,1:local_nx,ny,2) = zero
    DO j = 1,3+nscalar
      data1(1:nz,1:local_nx,ny+1,j) = zero
    END DO

!  CALCULATE THE CFL NUMBER
    IF (((mod(step,1)==0) .OR. (dt_flag == 1))  &
         .and. (rk_step==1)) THEN
       CALL check_cfl(ny)
    END IF

! COMPUTE A NEW TIME STEP TO KEEP THE CFL NEAR ONE
    IF ((STEP .GE. 0) .AND. (RK_STEP .EQ. 1)) THEN
       CFL_TEMP = MAXVAL(CFL(2:ny))
       IF ((me == 0) .and. ((step == 1) .or. (crossflow == 1) &
            .or. (mod(step,output_step) == 0))) THEN
          WRITE(*,*) 'MAX CFL = ', cfl_temp
       END IF
       IF (DT_FLAG .EQ. 1) THEN
          DT_TEMP = DT*CFL_MAX/CFL_TEMP
          DT = DMIN1(DT_TEMP, DT + (DT_TEMP-DT)*DT/(ONE + DT))
       ENDIF
    ENDIF

! END OF nonlinear_linearized
  CONTAINS
!=========================================================================
    SUBROUTINE compute_wall_parallel_flux(plane)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: plane
      REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz) :: fluxpu, fluxpw


! COMPUTE WALL-PARALLEL FLUXES AT yh(plane)

! WALL NORMAL FLUX OF U AND W AT y(plane).
      IF (plane < ny) THEN
         fluxpu(1:nx2,1:nzl) = &
              baseflow(1:nx2,1:nzl,plane,2)*half*(u(1:nx2,1:nzl,1)     &
              + up(1:nx2,1:nzl,1))    &
              + u(1:nx2,1:nzl,2)*half*(baseflow(1:nx2,1:nzl,plane,1)   &
              + baseflow(1:nx2,1:nzl,plane+1,1))
         fluxpw(1:nx2,1:nzl) = baseflow(1:nx2,1:nzl,plane,2)           &
              *half*(u(1:nx2,1:nzl,3) + up(1:nx2,1:nzl,3))           &
              + u(1:nx2,1:nzl,2)*half*(baseflow(1:nx2,1:nzl,plane,3)   &
              + baseflow(1:nx2,1:nzl,plane+1,3))
      ELSE
         fluxpu(1:nx2,1:nzl) = &
                baseflow(1:nx2,1:nzl,ny,2)*up(1:nx2,1:nzl,1)    &
              + u(1:nx2,1:nzl,2)*baseflow(1:nx2,1:nzl,ny+1,1)
         fluxpw(1:nx2,1:nzl) = &
                baseflow(1:nx2,1:nzl,ny,2)*up(1:nx2,1:nzl,3)    &
              + u(1:nx2,1:nzl,2)*baseflow(1:nx2,1:nzl,ny+1,3)
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
           two*baseflow(1:nx2,1:nzl,plane,1)*u(1:nx2,1:nzl,1), &! FLUX IN X-DIR
           yflux(1:nx2,1:nzl,1),              &  ! D/DY of Y-FLUX
           baseflow(1:nx2,1:nzl,plane,3)*u(1:nx2,1:nzl,1) &
           + u(1:nx2,1:nzl,3)*baseflow(1:nx2,1:nzl,plane,1), &  ! FLUX IN Z-DIR
           data1(1,1,plane,1))                       ! ARRAY RETURNS RHS
      CALL compute_divergence_of_fluxes(      &
           baseflow(1:nx2,1:nzl,plane,3)*u(1:nx2,1:nzl,1) &
           + u(1:nx2,1:nzl,3)*baseflow(1:nx2,1:nzl,plane,1), &  ! FLUX IN X-DIR
           yflux(1:nx2,1:nzl,3),              &  ! D/DY of Y-FLUX
           two*baseflow(1:nx2,1:nzl,plane,3)*u(1:nx2,1:nzl,3), &! FLUX IN Z-DIR
           data1(1,1,plane,3))                       ! ARRAY RETURNS RHS

! MOVE THE FLUX AT y(plane) INTO yflux.
      yflux(:,:,1) = fluxpu
      yflux(:,:,3) = fluxpw

    END SUBROUTINE compute_wall_parallel_flux
!=========================================================================
    SUBROUTINE compute_wall_normal_flux(plane)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: plane
      REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz) :: fluxp

! COMPUTE THE FLUXES FOR THE WALL-NORMAL VELOCITY AT y_k

! COMPUTE WALL-NORMAL FLUXES AT y(plane).
      temp = HALF*(u(1:nx2,1:nzl,2) + up(1:nx2,1:nzl,2))
      temp_baseflow = HALF*(baseflow(1:nx2,1:nzl,plane,2) &
                          + baseflow(1:nx2,1:nzl,plane+1,2))
      fluxp = two*temp*temp_baseflow

! TAKE THE DIVERGENCE OF FLUXES AND PUT INTO DATA1.
      yflux(1:nx2,1:nzl,2)= dn(plane)*(fluxp(:,:) - yflux(:,:,2))          &
           + BUFFER(1:nx2,1:nzl)*(u(1:nx2,1:nzl,2) - TARGET_INFLOW(plane,2))
      CALL compute_divergence_of_fluxes( &
             baseflow(1:nx2,1:nzl,plane,2)*(iu(1,plane)*u(1:nx2,1:nzl,1) &
                                          + iu(2,plane)*up(1:nx2,1:nzl,1)) &
           + u(1:nx2,1:nzl,2)*(iu(1,plane)*baseflow(1:nx2,1:nzl,plane,1) &
                             + iu(2,plane)*baseflow(1:nx2,1:nzl,plane+1,1)), &
           yflux(1:nx2,1:nzl,2),                                   &
             baseflow(1:nx2,1:nzl,plane,2)*(iu(1,plane)*u(1:nx2,1:nzl,3) &
                                          + iu(2,plane)*up(1:nx2,1:nzl,3)) &
           + u(1:nx2,1:nzl,2)*(iu(1,plane)*baseflow(1:nx2,1:nzl,plane,3) &
                             + iu(2,plane)*baseflow(1:nx2,1:nzl,plane+1,3)), &
           data1(1,1,plane,2))

! MOVE THE FLUX AT y(plane) INTO yflux.
      yflux(:,:,2) = fluxp

    END SUBROUTINE compute_wall_normal_flux
!=========================================================================
    SUBROUTINE compute_scalar_flux(plane)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: plane
      INTEGER  x, z, limit_scalar
      REAL(KIND=prec) T1, T2, T3, T4, T5, FLUXT, TMPD, TMPC, TMPU
      REAL(KIND=prec), DIMENSION(1:nx2,1:local_nz) :: fluxp

      REAL(KIND=prec), PARAMETER :: &
             QUICKU =-0.125d0,   CURVC = -1.0d0/12.0d0,  &
             QUICKC = 0.75d0,    CURVU = 1.0d0/24.0d0,   &
             QUICKD = 0.375d0,   SLOPE = 10.0d0

      limit_scalar = 1

! COMPUTE WALL-NORMAL FLUXES AT y(plane)
      DO i = 1,nscalar
         IF (plane == ny) THEN
            fluxp(1:nx2,1:nzl) = up(1:nx2,1:nzl,2)*up(1:nx2,1:nzl,3+i)
         ELSE
            IF (limit_scalar == 1) THEN
               DO z = 1,nzl
                  DO x = 1,nx2
                     IF (baseflow(X,Z,plane,2) .LE. ZERO) THEN
                        TMPD = U(X,Z,3+I)
                        TMPC = UP(X,Z,3+I)
                        TMPU = TP2(X,Z,I)
                        T5 = UPWIND(1,plane,2)*TMPD    &
                             + UPWIND(2,plane,2)*TMPC  &
                             + UPWIND(3,plane,2)*TMPU
                     ELSE
                        TMPD = UP(X,Z,3+I)
                        TMPC = U(X,Z,3+I)
                        TMPU = tm(X,Z,I)
                        T5 = UPWIND(1,plane,1)*TMPU    &
                             + UPWIND(2,plane,1)*TMPC  &
                             + UPWIND(3,plane,1)*TMPD
                     END IF
                     T1 = TMPC
                     T2 = TMPD
                     T3 = TMPU + SLOPE*TMPC - SLOPE*TMPU
                     T4 = DMAX1(DMIN1(T1,T2),DMIN1(T2,T3),DMIN1(T1,T3))
                     FLUXT = baseflow(X,Z,plane,2) &
                          *DMAX1(DMIN1(T1,T4),DMIN1(T4,T5),DMIN1(T1,T5))
                     IF (U(X,Z,2) .LE. ZERO) THEN
                        TMPD = baseflow(X,Z,plane,3+I)
                        TMPC = baseflow(X,Z,plane+1,3+I)
                        TMPU = baseflow(X,Z,plane+2,3+I)
                        T5 = UPWIND(1,plane,2)*TMPD    &
                             + UPWIND(2,plane,2)*TMPC  &
                             + UPWIND(3,plane,2)*TMPU
                     ELSE
                        TMPD = baseflow(X,Z,plane+1,3+I)
                        TMPC = baseflow(X,Z,plane,3+I)
                        TMPU = baseflow(X,Z,plane-1,I)
                        T5 = UPWIND(1,plane,1)*TMPU    &
                             + UPWIND(2,plane,1)*TMPC  &
                             + UPWIND(3,plane,1)*TMPD
                     END IF
                     T1 = TMPC
                     T2 = TMPD
                     T3 = TMPU + SLOPE*TMPC - SLOPE*TMPU
                     T4 = DMAX1(DMIN1(T1,T2),DMIN1(T2,T3),DMIN1(T1,T3))
                     FLUXT = fluxt &
                          + u(X,Z,2)*MAX(MIN(T1,T4),MIN(T4,T5),MIN(T1,T5))
                     fluxp(x,z) = FLUXT
                  END DO
               END DO
            ELSE
! APPLY QUADRATIC UPWIND INTERPOLATION OF SCALAR FLUX (QUICK SCHEME OF LEONARD)
               WHERE (U(1:nx2,1:nzl,2) <= zero)
                  fluxp(1:nx2,1:nzl) = baseflow(1:nx2,1:nzl,plane,2) &
                       *(upwind(1,plane,2)*u(1:nx2,1:nzl,3+i)   &
                       + upwind(2,plane,2)*up(1:nx2,1:nzl,3+i)  &
                       + upwind(3,plane,2)*tp2(1:nx2,1:nzl,i)) &
                       + u(1:nx2,1:nzl,2) &
                       *(upwind(1,plane,2)*baseflow(1:nx2,1:nzl,plane,3+i)   &
                       + upwind(2,plane,2)*baseflow(1:nx2,1:nzl,plane+1,3+i)  &
                       + upwind(3,plane,2)*baseflow(1:nx2,1:nzl,plane+2,3+i)) 
               ELSEWHERE
                  fluxp(1:nx2,1:nzl) = baseflow(1:nx2,1:nzl,plane,2) &
                       *(UPWIND(1,plane,1)*up(1:nx2,1:nzl,3+i)       &
                       + UPWIND(2,plane,1)*u(1:nx2,1:nzl,3+i)        &
                       + UPWIND(3,plane,1)*tm(1:nx2,1:nzl,i))        &
                       + u(1:nx2,1:nzl,2)                            &
                       *(UPWIND(1,plane,1)*baseflow(1:nx2,1:nzl,plane+1,3+i)  &
                       + UPWIND(2,plane,1)*baseflow(1:nx2,1:nzl,plane,3+i)   &
                       + UPWIND(3,plane,1)*baseflow(1:nx2,1:nzl,plane-1,3+i))
               END WHERE
            END IF
         END IF
         yflux(:,:,3+I) = DH(plane)*(fluxp - yflux(:,:,3+I)) &
              + buffer(1:nx2,1:nzl)                          &
              *(u(1:nx2,1:nzl,3+i) - target_inflow(plane,3+i))

! CALL ROUTINE TO COMPUTE DIVERGENCE OF FLUXES AND PUT INTO DATA1.
         CALL compute_divergence_of_fluxes(      &
              baseflow(1:nx2,1:nzl,plane,1)*u(1:nx2,1:nzl,3+i) &
              + u(1:nx2,1:nzl,1)*baseflow(1:nx2,1:nzl,plane,3+i), &  
              yflux(1,1,3+I),                    &
              baseflow(1:nx2,1:nzl,plane,3)*u(1:nx2,1:nzl,3+i) &  
              + u(1:nx2,1:nzl,3)*baseflow(1:nx2,1:nzl,plane,3+i), &  
              data1(1,1,plane,3+i))                ! RETURN RHS

! MOVE THE FLUX AT y(plane) INTO yflux.
         yflux(:,:,3+i) = fluxp
      END DO

    END SUBROUTINE compute_scalar_flux
!=========================================================================
    SUBROUTINE check_cfl(plane)
      IMPLICIT NONE
      
      INTEGER plane, x, z
      REAL(KIND=prec) TEMPA, TEMPB, TEMPC, big_cfl, poss_max, tmp
      REAL(KIND=prec), PARAMETER :: CFL_FOURIER =0.55d0,  CFL_FD =1.73d0, &
                                    DIFF_FOURIER=0.25d0,  DIFF_FD=2.51d0

      tempa = half*MAX(MAXVAL(ABS(u(:,:,1)+up(:,:,1))),  &
           MAXVAL(ABS(baseflow(:,:,plane,1)+baseflow(:,:,plane+1,1)))) &
           *DBLE(nx)/(cfl_fourier*lx)
      tempb = MAX(MAXVAL(ABS(u(:,:,2))),MAXVAL(ABS(baseflow(:,:,plane,2)))) &
           /(cfl_fd*min(hn(plane),hn(plane-1)))
      tempc = half*MAX(MAXVAL(ABS(u(:,:,3)+up(:,:,3))), &
           MAXVAL(ABS(baseflow(:,:,plane,3)+baseflow(:,:,plane+1,3)))) &
           *DBLE(nz)/(cfl_fourier*lz)
      poss_max = tempa + tempb + tempc
      CALL MPI_ALLREDUCE(poss_max,big_cfl,1,MPI_DOUBLE_PRECISION, &
           MPI_MAX, MPI_COMM_WORLD, MPIERR)
      cfl(plane) = big_cfl*dt

!!$      big_cfl = MAXVAL(tempa*half*abs(u(:,:,1)+up(:,:,1))   &
!!$           + tempb*abs(u(:,:,2)) + tempc*half*abs(u(:,:,3)+up(:,:,3)))
!         WRITE(*,*) STEP, RK_STEP, plane, CFL(plane)

!!$      IF (plane == 2) poss_max = zero
!!$      CALL MPI_ALLREDUCE(MAX(poss_max,MAXVAL(abs(u(:,:,1:3))),&
!!$           MAXVAL(abs(up(:,:,1:3)))),poss_max,&
!!$           1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,MPIERR)
!!$      IF ((plane==ny-1).and.(me==0)) WRITE(*,999) me,STEP,t_total,poss_max
!!$      999 FORMAT(i4,i6,2f12.4)

      IF (DT_FLAG .EQ. 1) THEN
         IF (CFL(plane) .GE. 3.0D+2) THEN    !DT IS TOO BIG  ... ABORT
            WRITE(*, 1000) STEP, RK_STEP, plane, CFL(plane)
            WRITE(*, 1010)
            CALL finalize_mpi
            STOP 'DT IS TOO BIG'
         ENDIF
      ELSE
         IF (CFL(plane) .GE. 1.15D0) THEN            !CFL IS TOO BIG ... ABORT
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
  END SUBROUTINE nonlinear_linearized
END MODULE pass1_linearized
