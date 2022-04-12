MODULE init_vel
  USE runparms
  USE gridparms
  USE data
  USE read_write
  USE diff_int
  USE init_jet_buffer
  USE solve
  USE pass1
!!$  USE pass1_dynamic2
  USE pass1_linearized
!!$  USE dynamic2
  USE pass2
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: define_initial_field
            
! miscellaneous loop indices.  
  INTEGER           i, j, k, x, z, plane

CONTAINS
!=================================================================
!=======================DEFINE INITIAL FIELD======================
!=================================================================
  SUBROUTINE DEFINE_INITIAL_FIELD
    IMPLICIT NONE

    COMPLEX(KIND=prec) :: cu(ny+1,3), TEMPU(NY+1), TEMPV(NY+1), &
                          topbc(1), botbc(1)
    REAL(KIND=prec)    :: tmp, alpha
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0, &
                                  EPS  = 1.0D-8,EPS2 = 1.0D-4

! INITIALIZE TIME TO ZERO.
    T_TOTAL = ZERO

!  INITIALIZE DATA1, DATA2 AND DATA4 WHICH HOLD THE VELOCITY, RHS,
!     AND PRESSURE FIELDS.
    data1 = cmplx(zero,zero)
    data2 = cmplx(zero,zero)
    data4 = cmplx(zero,zero)
   IF (nscalar .gt. 0) data_scalar = zero
!!$    data_rhs = zero
    
    IF (READ_VEL .EQ. 1) THEN

       CALL READ_VELOCITY
       IF (me == 0) WRITE(*,*) 'READ VELOCITY FIELD FROM savevel.out'

!         IF (STAT_WRITE .EQ. 2) THEN
!            OPEN (UNIT = 65, FILE = 'savestat.out', &
!                FORM = 'UNFORMATTED')
!            READ(65) STAT_TIME
!            READ(65) STAT_LAST
!            READ(65) STATS
!            READ(65) SPECTX
!            READ(65) SPECTZ
!            CLOSE(UNIT = 65)
!         ELSE
!            STAT_TIME = ZERO
!            STAT_LAST = T_TOTAL
!         END IF            

    ELSE

       IF (nscalar > 0) THEN
!!$          !! SPECIFY A GAUSSIAN SCALAR FIELD CENTERED IN MIDDLE OF DOMAIN.
!!$          DO k = 2,ny
!!$             DO z = 1,local_nz_small
!!$                DO x = 1,nx
!!$                   data_scalar(x,z,k,1) = EXP(-half*8.d0**2 &
!!$                        *((lx*DBLE(x-nx/2)/DBLE(nx2))**2 &
!!$                        + (lz*DBLE(z+local_z_small_start-nz/2)/DBLE(nz))**2&
!!$                        + (yh(k) - yh(ny/2+1))**2))
!!$                END DO
!!$             END DO
!!$          END DO

          !! SPECIFY A TOPHAT SCALAR FIELD IN MIDDLE OF DOMAIN.
          DO z = 1,local_nz_small
             IF ((z+local_z_small_start >= nz/4) .and. &
                  (z+local_z_small_start <= 3*nz/4)) THEN
                data_scalar(nx/4:3*nx/4,z,ny/4+1:3*ny/4+1,1) = one
             END IF
          END DO
       END IF

       IF (WHICH_TEST .EQ. 1) THEN

          IF (local_x_start == 0) THEN
             DATA1(1,1,1:ny+1,1) = (ONE - YH(1:ny+1)**2)*cos(beta)
             DATA1(1,1,1:ny+1,3) = (ONE - YH(1:ny+1)**2)*sin(beta)
          END IF

       ELSEIF (WHICH_TEST .EQ. 2) THEN

!!$          CALL INIT_BASEFLOW
!!$          linearized_run = 1
!!$
          IF (local_x_start == 0) THEN
             DATA1(1,1,1:ny+1,1) = ONE - YH(1:ny+1)**2
          END IF
          
          IF ((local_x_start+1 <= 2) .and. (local_x_start+local_nx >= 2)) THEN
             x = 2 - (local_x_start)
             z = 1
             CALL READ_EIGENFUNCTION(TEMPU, TEMPV)
             data1(z,x,1:ny+1,1) = disturb_ampl*tempu(1:ny+1)
             data1(z,x,1:ny,2)   = disturb_ampl*tempv(1:ny)

             tmp = SUM(dyh(2:ny)*(ABS(data1(z,x,2:ny,1))**2              &
                  + ABS(half*(data1(z,x,1:ny-1,2)+data1(z,x,2:ny,2)))**2))
             WRITE(80,990) T_TOTAL, tmp, P_GRAD
             990 FORMAT(F10.3,2E24.14)
          END IF
       
       ELSEIF (WHICH_TEST .EQ. 5) THEN

!!$          CALL INIT_BASEFLOW
!!$          linearized_run = 1
!!$
          IF (local_x_start == 0) THEN
             DATA1(1,1,1:ny+1,3) = ONE - YH(1:ny+1)**2
             CALL READ_EIGENFUNCTION(TEMPU, TEMPV)
             x = 1
             z = 2
             data1(2,1,1:ny+1,3) = disturb_ampl*tempu(1:ny+1)
             data1(2,1,1:ny,2)   = disturb_ampl*tempv(1:ny)

             data1(nz,1,1:ny+1,3) = disturb_ampl*CONJG(tempu(1:ny+1))
             data1(nz,1,1:ny,2)   = disturb_ampl*CONJG(tempv(1:ny))

             tmp = SUM(dyh(2:ny)*(ABS(data1(z,x,2:ny,1))**2              &
                  + ABS(half*(data1(z,x,1:ny-1,2)+data1(z,x,2:ny,2)))**2))
             WRITE(80,990) T_TOTAL, tmp, P_GRAD
          END IF
          
       ELSEIF (WHICH_TEST .EQ. 8) THEN

!!$          CALL INIT_BASEFLOW
!!$          linearized_run = 1
!!$
          IF (local_x_start == 0) THEN
             DATA1(1,1,1:ny+1,1) = (ONE - YH(1:ny+1)**2)*cos(beta)
             DATA1(1,1,1:ny+1,3) = (ONE - YH(1:ny+1)**2)*sin(beta)
          END IF
          
          IF ((local_x_start+1 <= 2) .and. (local_x_start+local_nx >= 2)) THEN
             x = 2 - (local_x_start)
             z = 2
             CALL READ_EIGENFUNCTION(TEMPU, TEMPV)
             data1(z,x,1:ny+1,1) = disturb_ampl*tempu(1:ny+1)*cos(beta)
             data1(z,x,1:ny+1,3) = disturb_ampl*tempu(1:ny+1)*sin(beta)
             data1(z,x,1:ny,2)   = disturb_ampl*tempv(1:ny)

             tmp = SUM(dyh(2:ny)*(ABS(data1(z,x,2:ny,1))**2              &
                  + ABS(half*(data1(z,x,1:ny-1,2)+data1(z,x,2:ny,2)))**2))
             WRITE(80,990) T_TOTAL, tmp, P_GRAD
          END IF
          
       ELSEIF (WHICH_TEST .EQ. 6) THEN

          p_grad = zero
          fix_mass = 0
          
          IF ((local_x_start+1 <= 2) .and. (local_x_start+local_nx >= 2)) THEN
             x = 2 - (local_x_start)
             z = 1
             data1(z,x,1:ny+1,1) = cmplx(-half*sin(pi*yh(1:ny+1)),zero)
             data1(z,x,1:ny,2)   = cmplx(zero,-half*cos(pi*y(1:ny)))
          END IF
          
       ELSEIF (WHICH_TEST .EQ. 3) THEN

          linearized_run = 1
          CALL INIT_BASEFLOW

! INITIALIZE WHITE RANDOM VELOCITY FIELD.
          DO x = 1,local_nx
             DO z = 1,nz
                CALL random_divfree_vector(0,x,z,disturb_ampl,cu)
                IF (k2_me(z,x) .ne. zero) THEN
                   data1(z,x,1:ny+1,1:3) = cu                   
                ELSE
                   ! ADD PERTURBATIONS TO MEAN FLOW IN STREAMWISE AND
                   ! SPANWISE DIRECTION.  THESE PERTURBATIONS ARE IN
                   ! ADDITION TO THE BASE FLOW INITIALIZED ABOVE.
                   data1(z,x,1:ny+1,1) = DBLE(cu(1:ny+1,1))
                   data1(z,x,1:ny+1,3) = DBLE(cu(1:ny+1,3))
                END IF
             END DO
          END DO

!          IF (local_x_start == 0) THEN
!             DATA1(1,1,1:ny+1,1) = ONE - YH(1:ny+1)**2
!          END IF
!            DO PLANE = 1,NY+1
!               DATA1(1,PLANE,3,1) = ZERO
!               DATA1(2,PLANE,3,3) = - EPS2*PI*DSIN(PI*YH(PLANE))/KZ(3)
!               DATA1(1,PLANE,NZ-1,1) = ZERO
!               DATA1(2,PLANE,NZ-1,3) = - EPS2*PI*DSIN(PI*YH(PLANE))/KZ(NZ-1)
!            END DO
!            DO PLANE = 1,NY
!               DATA1(1,PLANE,3,2) = EPS2*(ONE + DCOS(PI*Y(PLANE)))
!               DATA1(1,PLANE,NZ-1,2) = EPS2*(ONE + DCOS(PI*Y(PLANE)))
!            END DO            

       ELSEIF (WHICH_TEST .EQ. 4) THEN

          IF (local_x_start == 0) THEN
             DATA1(1,1,1:ny+1,1) = 0.75d0*(one - yh(1:ny+1)**8)
             IF ((crossflow == 1) .or. (bl_grid > 0)) THEN
                DATA1(1,1,1:ny+1,1) = ublasius(1:ny+1)
             END IF
          END IF

          DO x = 1,local_nx
             DO z = 1,nz
                IF (k2_me(z,x) .ne. zero) THEN
                   ! SPECIFY RANDOM INITIAL FIELD W/ k^-2 DECAY AT HIGH WAVENO.
                   CALL random_divfree_vector(0,x,z,disturb_ampl,cu)
                   data1(z,x,1:ny+1,1:3) = cu*min(one,k2_me(z,x)**(-1))
                END IF
             END DO
          END DO

       ELSEIF (WHICH_TEST .EQ. 9) THEN

          DO x = 1,local_nx
             DO z = 1,nz
                IF (k2_me(z,x) == zero) THEN
                   DATA1(z,x,1:ny+1,1) = 0.75d0*(one - yh(1:ny+1)**8)*COS(beta)
                   DATA1(z,x,1:ny+1,3) = 0.75d0*(one - yh(1:ny+1)**8)*SIN(beta)
                   IF ((crossflow == 1) .or. (bl_grid > 0)) THEN
                      DATA1(z,x,1:ny+1,1) = ublasius(1:ny+1)
                   END IF
                ELSE 
                   CALL random_divfree_vector(1,x,z,disturb_ampl,cu)
                   data1(z,x,1:ny+1,1:3) = cu
                END IF
             END DO
          END DO

       END IF

    END IF

    IF (me==0) WRITE(*,*) 'INITIALIZED VELOCITY'

! IF BL OR JET IN CROSSFLOW SIMULATION, USE BLASIUS MEAN FLOW
    IF (((bl_grid > 0) .or. (crossflow==1)) .and. (read_vel==0)) THEN
       IF ((me == 0) .and. (linearized_run == 0)) THEN
          DO plane = 1,ny+1
             data1(1,1,plane,1) = ublasius(plane)*COS(beta)
             data1(1,1,plane,3) = ublasius(plane)*SIN(beta)
          END DO
       END IF
    END IF

! IF JET IN CROSSFLOW SIMULATION OR FRINGE REGION IS REQUIRED, INITIALIZE THEM.
    IF (crossflow > 0) THEN
       !! CALL ROUTINE TO SET UP BOUNDARY CONDITIONS.
       CALL init_jet

       IF ((read_vel == 0) .and. (jetmag_bot /= zero)) THEN
          !! SET UP BOUNDARY CONDITIONS AT BOTTOM OF DOMAIN
          !! FILTER INITIAL CONDITION USING GAUSSIAN FILTER
          !! WITH FILTER WIDTH = MAX(deltax,deltay_wall,deltaz)
          alpha = (8.d0*MAX(lx/DBLE(nx),hn(1),lz/DBLE(nz)))**2/24.d0
          DO z = 1,nz
             DO x = 1,local_nx
                cu = zero
                topbc = zero
                botbc = swirl_bot*cbc_bot(z,x,1)
                CALL solve_helmholtz(ny+1,1,cu(1,1), &
                     one+alpha*k2_me(z,x),alpha,d2h,topbc,botbc,bl_grid,0)
                botbc = jetmag_bot*cbc_bot(z,x,2)
                CALL solve_helmholtz(ny,1,cu(1,2), &
                     one+alpha*k2_me(z,x),alpha,d2n,topbc,botbc,bl_grid,0)
                botbc = swirl_bot*cbc_bot(z,x,3)
                CALL solve_helmholtz(ny+1,1,cu(1,3), &
                     one+alpha*k2_me(z,x),alpha,d2h,topbc,botbc,bl_grid,0)
                data1(z,x,1:ny+1,1:3) = data1(z,x,1:ny+1,1:3) + cu
             END DO
          END DO
!!$          data1(1:nz,1:local_nx,1,1) = swirl_bot*cbc_bot(1:nz,1:local_nx,1)
!!$          data1(1:nz,1:local_nx,1,2) = jetmag_bot*cbc_bot(1:nz,1:local_nx,2)
!!$          data1(1:nz,1:local_nx,1,3) = swirl_bot*cbc_bot(1:nz,1:local_nx,3)
          IF (nscalar .gt. 0) THEN
             data_scalar(1:nx,1:local_nz_small,1,1) = &
                  bc_bot(1:nx,1:local_nz_small,4)
          END IF
       END IF

       IF ((read_vel == 0) .and. (jetmag_top /= zero)) THEN
          !! SET UP BOUNDARY CONDITIONS AT TOP OF DOMAIN
          !! FILTER INITIAL CONDITION USING GAUSSIAN FILTER
          !! WITH FILTER WIDTH = MAX(deltax,deltay_wall,deltaz)
          alpha = (8.d0*MAX(lx/DBLE(nx),hn(1),lz/DBLE(nz)))**2/24.d0
          DO z = 1,nz
             DO x = 1,local_nx
                cu = zero
                botbc = zero
                topbc = swirl_top*cbc_top(z,x,1)
                CALL solve_helmholtz(ny+1,1,cu(1,1), &
                     one+alpha*k2_me(z,x),alpha,d2h,topbc,botbc,bl_grid,0)
                topbc = jetmag_top*cbc_top(z,x,2)
                CALL solve_helmholtz(ny,1,cu(1,2), &
                     one+alpha*k2_me(z,x),alpha,d2n,topbc,botbc,bl_grid,0)
                topbc = swirl_top*cbc_top(z,x,3)
                CALL solve_helmholtz(ny+1,1,cu(1,3), &
                     one+alpha*k2_me(z,x),alpha,d2h,topbc,botbc,bl_grid,0)
                data1(z,x,1:ny+1,1:3) = data1(z,x,1:ny+1,1:3) + cu
             END DO
          END DO
!!$       data1(1:nz,1:local_nx,ny+1,1) = swirl_top*cbc_top(1:nz,1:local_nx,1)
!!$       data1(1:nz,1:local_nx,ny,2)   = jetmag_top*cbc_top(1:nz,1:local_nx,2)
!!$       data1(1:nz,1:local_nx,ny+1,3) = swirl_top*cbc_top(1:nz,1:local_nx,3)
          IF (nscalar .gt. 0) THEN
             data_scalar(1:nx,1:local_nz_small,1,1) = &
                  bc_top(1:nx,1:local_nz_small,4)
          END IF
       END IF
    END IF
    IF (xfringe+zfringe > 0) CALL init_buffer
    IF (linearized_run == 1) target_inflow = zero

!  REMOVE INITIAL DIVERGENCE FROM VELOCITY FIELD.
    CALL remove_initial_divergence
    IF (me==0) WRITE(*,*) 'REMOVED DIVERGENCE'

!  COMPUTE D_UM_DY, THE MEAN VELOCITY GRADIENT.
    IF (local_x_start == 0) THEN	
       D_UM_DY(1) = DBLE(DM(1,1)*DATA1(1,1,1,1) &
                       + DM(2,1)*DATA1(1,1,2,1) &
                       + DM(3,1)*DATA1(1,1,3,1))
       D_UM_DY(2:ny-1) = &
            DBLE(DN(2:ny-1)*(DATA1(1,1,3:ny,1) - DATA1(1,1,2:ny-1,1)))
       D_UM_DY(NY) = DBLE(DM(1,2)*DATA1(1,1,NY-1,1) &
                        + DM(2,2)*DATA1(1,1,NY,  1) &
                        + DM(3,2)*DATA1(1,1,NY+1,1))
    END IF

    DO X = 1,Klim
       DO Z = 1,Mlim
          data2(1:ny+1,1:3,z,x) = data1(z,x,1:ny+1,1:3)
       END DO
    END DO

! COMPUTE THE INITIAL PRESSURE FIELD FROM THE VELOCITY.
    CALL compute_initial_pressure

    IF (me==0) WRITE(*,*) 'INITIALIZED VELOCITY AND PRESSURE'

  CONTAINS
    !=====================================================================
    SUBROUTINE random_divfree_vector(filter,p,q,ampl,cu)
      IMPLICIT NONE
      COMPLEX(KIND=prec), DIMENSION(ny+1,3) :: cu
      INTEGER :: filter, p, q !! x- AND z-WAVENUMBERS OF THIS FOURIER MODE.
      REAL(KIND=prec) :: ampl !! APPROX AMPLITUDE OF VECTOR BEFORE FILTERING.
      
      REAL(KIND=prec) :: TEMP(NY+1,4), TEMP2(NY+1,4), &
                         ushape(1:ny+1), vshape(1:ny), TMP

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
  END SUBROUTINE define_initial_field
  !=====================================================================
  !===================COMPUTE INITIAL PRESSURE==========================
  !=====================================================================
  SUBROUTINE COMPUTE_INITIAL_PRESSURE
    IMPLICIT NONE

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    !  COMPUTE THE NONLINEAR TERM FOR THE INITIAL VELOCITY FIELD.

    step = -1
    RK_STEP = 3
    IF (me==0) WRITE(*,*) 'CALL NONLINEAR'
    IF ((READ_VEL .EQ. 0) .AND. (WHICH_TEST .EQ. 3)) THEN
!!!!!!!!!!!         CALL nonlinear_linearized
       CALL nonlinear            
    ELSE
       IF (subgrid==1) THEN
!!$          IF (lagr_avg == 1) THEN
!!$             IF (MAXVAL(data_dyn)<10.d0*EPSILON(1.d0)) THEN
!!$                CALL compute_dynamic_coef(filter_dim,0)
!!$             END IF
!!$          ELSEIF (lagr_avg==0) THEN
!!$             CALL compute_dynamic_coef(filter_dim,0)
!!$          END IF
!!$          CALL nonlinear_dynamic2
       ELSE
          CALL nonlinear
       END IF
    END IF
    IF (me==0) WRITE(*,*) 'NONLINEAR DONE'

    CALL solve_continuous_poisson

    ! RETURN INITIAL VELOCITY/SCALAR FIELD TO data1
    DO X = 1,Klim
       DO Z = 1,Mlim
          data1(z,x,1:ny+1,1:3) = data2(1:ny+1,1:3,z,x)
       END DO
    END DO

  CONTAINS
    !=====================================================================
    SUBROUTINE solve_continuous_poisson
      IMPLICIT NONE

      COMPLEX(KIND=prec), DIMENSION(ny+1,3) :: cu, cru
      COMPLEX(KIND=prec), DIMENSION(ny+1)   :: du
      REAL(KIND=prec)    visc

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
               !  FOR EACH WAVENUMBER PAIR (X,Z), TAKE THE DIVERGENCE OF THE RHS, COMPUTE 
               !    THE BOUNDARY CONDITIONS AND SOLVE THE POISSON EQUATION FOR THE PRESSURE
               CALL strip_divergence(cru,data4(1:ny+1,z,x),one,x,z)
            END IF
         END DO
      END DO
    END SUBROUTINE solve_continuous_poisson
  END SUBROUTINE compute_initial_pressure
  !=====================================================================
  !=============================REMOVE DIVERGENCE=======================
  !=====================================================================
  SUBROUTINE REMOVE_INITIAL_DIVERGENCE
    IMPLICIT NONE

    COMPLEX(KIND=prec)         :: CH(NY+1)
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    DO X = 1,Klim
       DO Z = 1,Mlim
          IF (K2_ME(z,x) .EQ. ZERO) THEN
! set the mean vertical velocity to that specified at the bottom wall.
             data1(z,x,2:ny,2) = data1(z,x,1,2)
          ELSE
! REMOVE THE DIVERGENCE FROM THE VELOCITY FIELD IN THIS WAVENUMBER PAIR.
             CALL strip_divergence(data1(z,x,1:ny+1,1:3),ch,one,x,z)
          END IF
       END DO
    END DO

  END SUBROUTINE remove_initial_divergence
  !=====================================================================
  !=============================READ EIGENFUNCTION====================
  !=====================================================================
  SUBROUTINE READ_EIGENFUNCTION(U,V)

    !   THIS SUBROUTINE READS IN THE FASTEST-GROWING EIGENFUNCTION 
    !     OF THE ORR-SOMMERFELD OPERATOR AT A REYNOLDS NUMBER OF 7500.
    !     THE EIGENFUNCTION WAS COMPUTED IN MATLAB USING A CHEBYSHEV
    !     COLLOCATION TECHNIQUE, SO WE INTERPOLATE THE EIGENFUNCTION
    !     (USING THE CHEBYSHEV POLYNOMIALS) ONTO THE STRETCHED FINITE
    !     DIFFERENCE GRID.

    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(OUT) :: U(NY+1), V(NY+1)

    INTEGER, PARAMETER :: NCHEB = 65
    REAL(KIND=prec)    :: N, UIN(NCHEB,4), UCHEB(NCHEB,4), T(NCHEB,NCHEB)
    REAL(KIND=prec)    :: SUM, UTMP(NY+1,4)
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    !  DEFINE N, A USEFUL CONSTANT.
    N = DBLE(NCHEB-1)

    !  INITIALIZE THE OUTPUT.
    u = zero
    v = zero

    !  READ IN THE EIGENFUNCTION OF THE ORR-SOMMERFELD EQUATION.
    OPEN(UNIT = 19, FILE = 'efn.m', STATUS = 'OLD', FORM = 'FORMATTED') 
    DO I = 1,4
       DO PLANE = 1,NCHEB
          READ(19,*) UIN(PLANE,I)
       END DO
    END DO
    CLOSE(UNIT = 19, STATUS = 'KEEP')

    !  WE NOW HAVE THE EIGENFUNCTION ON A CHEBYSHEV COLLOCATION GRID.
    !     WE NEED TO GET IT ON OUR (HYPERBOLIC TANGENT) STRETCHED GRID.
    !     TWO STEPS: 
    !
    !     1.  TRANSFORM FROM A PHYSICAL TO CHEBYSHEV COEFFICIENT SPACE.
    !         TO DO THIS, WE BUILD THE CHEBYSHEV TRANSFORM MATRIX, T.
    !     2.  INTERPOLATE THE EIGENFUNCTION ON OUR GRID USING THE 
    !         CHEBYSHEV POLYNOMIALS.

    !  BUILD THE TRANSFORM MATRIX, T.
    DO J = 1,NCHEB
       DO I = 1,NCHEB
          T(I,J) = (TWO/N)*DCOS(PI*DBLE((I-1)*(J-1))/N)
          IF ((I .EQ. 1) .OR. (I .EQ. NCHEB)) T(I,J) = HALF*T(I,J)
          IF ((J .EQ. 1) .OR. (J .EQ. NCHEB)) T(I,J) = HALF*T(I,J)
       END DO
    END DO

    !  TRANSFORM UIN INTO COEFFICIENT SPACE.
    ucheb = zero

    DO K = 1,4
       DO I = 1,NCHEB
          DO J = 1,NCHEB
             UCHEB(I,K) = UCHEB(I,K) + T(I,J)*UIN(J,K)
          END DO
       END DO
    END DO

    !  INTERPOLATE U ONTO THE HALF POINTS, YH(J).
    DO K = 1,2
       DO I = 2,NY
          UTMP(I,K) = ZERO
          DO J = 1,NCHEB
             UTMP(I,K) = UTMP(I,K) + COS(DBLE(J-1)*ACOS(YH(I)))*UCHEB(J,K)
          END DO
       END DO
    END DO

    !  INTERPOLATE V ONTO THE POINTS, Y(J), J = 1,...,NY.
    DO K = 3,4
       DO I = 1,NY
          UTMP(I,K) = ZERO
          DO J = 1,NCHEB
             UTMP(I,K) = UTMP(I,K) + COS(DBLE(J-1)*ACOS(Y(I)))*UCHEB(J,K)
          END DO
       END DO
    END DO

    !  PUT THIS INTO THE VELOCITY VECTOR.
    U(1) = ZERO
    DO PLANE = 2,NY
       U(PLANE) = CMPLX(UTMP(PLANE,1),UTMP(PLANE,2),KIND=prec)
    END DO
    U(NY+1) = ZERO

    DO PLANE = 1,NY
       V(PLANE) = CMPLX(UTMP(PLANE,3),UTMP(PLANE,4),KIND=prec)
    END DO

    !  CHECK THE DIVERGENCE OF THE EIGENFUNCTION
    U(1) = ZERO
    DO PLANE = 2,NY
       !         U(PLANE) = - (DH(PLANE)*(V(PLANE)-V(PLANE-1)))/IKX(2)
    END DO
    U(NY+1) = ZERO

    !  LAST THING: NORMALIZE SO THAT THE INITAIL DISTURBANCE
    !     HAS MAGNITUDE ONE.
    SUM = ZERO
    DO PLANE = 2,NY
       SUM = SUM + DYH(PLANE)*(ABS(U(PLANE))**2  &
            + ABS(HALF*V(PLANE)+HALF*V(PLANE-1))**2)
    END DO

    SUM = ONE/DSQRT(SUM)
    DO PLANE = 1,NY+1
       U(PLANE) = SUM*U(PLANE)
       V(PLANE) = SUM*V(PLANE)
    END DO

  END SUBROUTINE read_eigenfunction
END MODULE init_vel
