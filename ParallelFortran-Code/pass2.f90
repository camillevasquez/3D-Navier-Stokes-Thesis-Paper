MODULE pass2
  USE runparms
  USE gridparms
  USE data
  USE solve
  USE diff_int
  USE transform
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rhs_solve, jet_profile
  
  INCLUDE 'mpif.h'

  REAL(KIND=prec), PARAMETER, DIMENSION(3) :: &
         rka = (/  0.0,      -17.0/60.0, -5.0/12.0 /), &
         rkb = (/  8.0/15.0,   5.0/12.0,  3.0/ 4.0 /), &
         rkc = (/ 29.0/96.0,  -3.0/40.0,  1.0/ 6.0 /), &
         rkd = (/ 37.0/160.0,  5.0/24.0,  1.0/ 6.0 /)

  REAL(KIND=prec) tmp_time, top_jet_amplitude, bot_jet_amplitude, &
       top_swirl_amplitude, bot_swirl_amplitude

  INTEGER :: i, j, k, x, z, plane, mpierr

CONTAINS
!=====================================================================
!==============================RHS SOLVE==============================
!=====================================================================
  SUBROUTINE rhs_solve
    IMPLICIT NONE

    REAL(KIND=prec)  KE(local_nx,NZ), TMP_SUM, TEMP(NY+1), tmp_incr
    COMPLEX(KIND=prec), DIMENSION(ny+1,3) :: cu

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! SET THE TIME CORRESPONDING TO THE END OF THIS RK_STEP.
    tmp_time = t_total + dt*SUM(rkc(1:rk_step)+rkd(1:rk_step))
    tmp_incr = dt*(rkc(rk_step) + rkd(rk_step))

! SET THE (POSSIBLY TIME-VARYING) AMPLITUDE OF THE CROSSFLOW JET.
    IF (crossflow == 1) THEN
       CALL jet_profile(tmp_time,tmp_incr,cbc_top(1,1,2),cbc_bot(1,1,2))
       bot_swirl_amplitude = swirl_bot/(jetmag_bot + EPSILON(jetmag_bot))
       IF (bl_grid == 0) THEN
          top_swirl_amplitude = swirl_top/(jetmag_top + EPSILON(jetmag_top))
       ELSE
          top_jet_amplitude = zero
          top_swirl_amplitude = zero
       END IF
    ELSE
       top_jet_amplitude = one
       bot_jet_amplitude = one
       top_swirl_amplitude = one
       bot_swirl_amplitude = one
    END IF

!  LOOP THROUGH THE (X,Z) WAVENUMBER PAIRS INDIVIDUALLY
    DO X = 1,Klim
       DO Z = 1,Mlim

!  INITIALIZE THE VARIABLES.
          cu  = zero

!  DOWNLOAD THE OLD VELOCITY AND PRESSURE FIELDS FROM DATA1
! NEGATE THE NONLINEAR TERM (IN cu/data1) SINCE IT IS MOVED TO RIGHT HAND SIDE
          cu  = -data1(z,x,1:ny+1,1:3)

!  SET UP THE RHS AND COMPUTE THE NEW VELOCITY AND PRESSURE FIELDS.
          CALL PHASE1(cu,data2(1,1,z,x),data4(1,z,x))

!  SAVE THE VELOCITY AND PRESSURE FIELDS
          data1(z,x,1:ny+1,1:3) = cu
!!$          data2(1:ny+1,1:3,z,x) = cru
!!$          data4(1:ny+1,z,x) = cp

!!$          unit = 80 + me
!!$          WRITE(unit,880) z, x, (cu(j,4),j=1,3)
!!$          WRITE(unit,880) z, x, (data1(z,x,j,4),j=1,3)
880       FORMAT(2i4,6e14.6)

          
!  IF WE ARE IN THE MEAN MODE, UPDATE THE MEAN VELOCITY GRADIENT.
          IF (k2_me(z,x) == zero) THEN
             d_um_dy(1) = dble(DM(1,1)*cu(1,1)+DM(2,1)*cu(2,1)+DM(3,1)*cu(3,1))
             d_um_dy(2:ny-1) = dble(dn(2:ny-1)*(cu(3:ny,1) - cu(2:ny-1,1)))
             d_um_dy(ny) = &
                  dble(DM(1,2)*cu(NY-1,1)+DM(2,2)*cu(ny,1)+DM(3,2)*cu(ny+1,1))
          END IF

!  OUTPUT THE DATA.
          IF ((RK_STEP == RK_LAST) .AND. &
               ((step==tend).or.(MOD(STEP,OUTPUT_STEP) == 0))) THEN

! WRITE MEAN STREAMWISE VELOCITY PROFILE TO vel.out
!  AS WELL AS VOLUME FLOW, PRESSURE GRADIENT TO STANDARD OUTPUT.
!  AND TIME, TOP WALL SHEAR STRESS AND BOTTOM WALL SHEAR STRESS TO drag.out.
             IF (k2_me(z,x) == zero) THEN

                WRITE(12,910) (DBLE(CU(PLANE,1)),PLANE=1,NY+1)
                IF ((subgrid==1).and.(lagr_avg==0)) THEN
                   WRITE(10,910) (dyn_visc(PLANE,1),PLANE=1,NY)
                END IF
910             FORMAT(1000F16.12)

                temp(1:ny+1) = dble(cu(1:ny+1,1))
                CALL int_y(ny+1,temp,dyh,tmp_sum)
                WRITE(*,930) TMP_SUM, P_GRAD, STEP, T_TOTAL+dt
!		drag.out is writen with l2.out from read_write.f90
!               WRITE(11,"(F10.3,F12.6)") T_TOTAL+dt, D_UM_DY(1)- D_UM_DY(NY)
             END IF
          
! COMPUTE KINETIC ENERGY IN THIS WAVENUMBER AND STORE IN KE(X,Z).
             temp(1:ny+1) = abs(cu(1:ny+1,1))**2 + abs(cu(1:ny+1,3))**2
             temp(1)    = temp(1)    + abs(cu(1,2))**2
             temp(2:ny) = temp(2:ny) + abs(half*(cu(2:ny,2)+cu(1:ny-1,2)))**2
             temp(ny+1) = temp(ny+1) + abs(cu(ny,2))**2
             CALL int_y(ny+1,temp,dyh,ke(x,z))

!!$             IF ((x_start(me)+x == 2) .and. (z == 1) .and. &
!!$                  (READ_VEL == 0) .AND. (WHICH_TEST == 3)) THEN
!!$                WRITE(80,990) T_TOTAL+dt, KE(X,Z), &
!!$                     half*log(ke(x,z)+1.d-16)/(t_total+dt)
!!$             END IF
!!$              
             IF ((x_start(me)+x == 2) .and. (z == 1) .and. &
                  (READ_VEL == 0) .AND. (WHICH_TEST == 2)) THEN
                WRITE(80,990) T_TOTAL+dt, KE(X,Z), P_GRAD, &
                     half*log(ke(x,z)+1.d-16)/(t_total+dt)
!!$		WRITE(76,980) (DBLE(CU(J,1)),j=1,NY+1), &
!!$                     (aimag(CU(J,1)),j=1,NY+1),&
!!$                     (DBLE(CU(J,2)),j=1,NY), (aimag(CU(J,2)),j=1,NY)
             END IF

!  CHECK THE DIVERGENCE OF THE VELOCITY
             IF (DIVG_CHECK == 1) CALL CHECK_DIVERGENCE(CU,x,z)

          END IF
          
!!$          IF ((x_start(me)+x == 2) .and. (z == 1) .and. &
!!$               (STEP == tend) .and. (rk_step==rk_last) .and. &
!!$               (READ_VEL == 0) .AND. (WHICH_TEST == 6)) THEN
!!$                WRITE(80,990) T_TOTAL+dt, KE(X,Z), P_GRAD
!!$                do plane = 1,ny
!!$                   write(*,970) plane, &
!!$                        abs(cu(plane,1) + cmplx(half*sin(pi*yh(plane))&
!!$                        *exp(-two*pi**2*tmp_time/re),zero)), &
!!$                        abs(cu(plane,2) + cmplx(zero,half*cos(pi*y(plane)) &
!!$                        *exp(-two*pi**2*tmp_time/re))) ,&
!!$                        data4(plane,1,3)
!!$                end do
!!$                write(*,970) ny+1, &
!!$                     abs(dble(cu(ny+1,1))                                &
!!$                    +half*sin(pi*yh(ny+1))*exp(-two*pi**2*tmp_time/re)), &
!!$                        data4(ny+1,1,3)
!!$             END IF
              
! END LOOPS OVER WAVENUMBERS
       END DO
    END DO

    IF ((RK_STEP == RK_LAST) .AND. (MOD(STEP,OUTPUT_STEP) == 0)) THEN
       IF (local_x_start <= 5) THEN
          IF ((me == 0) .or. (read_vel == 0)) THEN
             DO Z = 1,MIN(NZ,10)
                WRITE(*,950) (KE(X,Z),X=1,MIN(local_nx,6))
             END DO
          END IF
       END IF
    END IF

    ! 920  FORMAT(F8.3,4F12.6)
    930  FORMAT('Q = ',F14.10,'  P GRAD = ',E14.6, &
              '  STEP = ',I6,'  TIME = ',F8.2)
    950  FORMAT(6E12.4)
    970  FORMAT(i4,4E18.10)
    980  FORMAT(1000E24.16)
    990  FORMAT(F10.3,2E24.14)

  END SUBROUTINE rhs_solve
!=====================================================================
!================================PHASE 1==============================
!=====================================================================
  SUBROUTINE PHASE1(CU, CRU, CP)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(INOUT), DIMENSION(NY+1,3) :: CU, CRU
    COMPLEX(KIND=prec), INTENT(INOUT), DIMENSION(NY+1)   :: CP

    COMPLEX(KIND=prec)  :: CV(NY+1,3), ch(ny+1), DU(ny+1), TOPBC(1), BOTBC(1)
    REAL(KIND=prec)     :: VIS_EXP, VIS_IMP, COEF, COEF2, RK_PRESS, ALPHA

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    995 format(12f12.6)
    996 format(12e12.4)

    VIS_IMP = RKD(RK_STEP)*(DT/RE)
    RK_PRESS = (RKC(RK_STEP)+RKD(RK_STEP))*DT

!  COMBINE THE EXPLICIT PART OF LINEAR TERMS WITH THE RHS
    COEF    = RKB(RK_STEP)*DT
    cv = cru + coef*cu

    IF (RK_STEP == 1) THEN
!  ADD EXPLICIT PART OF VISCOUS TERM.
       VIS_EXP = RKC(RK_STEP)*(DT/RE)
       CALL add_explicit_diffusive_terms(cv,cru)

!  ADD PRESSURE GRADIENT TO NON-ZERO WAVENUMBERS.
       IF (k2_me(z,x) /= zero) THEN
          CALL subtract_pressure_gradient(cv,cp,rk_press,x,z)
       END IF
    END IF

!  SOLVE THE HELMHOLTZ EQUATIONS FOR THE NEW VELOCITY FIELD.
    IF (K2_ME(Z,X) == ZERO) THEN
!  ADD THE MEAN PRESSURE GRADIENT TO THE RHS.
      IF (FIX_MASS == 1) CALL compute_p_grad(cu)
      cv(1:ny+1,1) = cv(1:ny+1,1) - cos(beta)*rk_press*p_grad
      cv(1:ny+1,3) = cv(1:ny+1,3) - sin(beta)*rk_press*p_grad

!  SOLVE EQUATIONS FOR MEAN MODE. 
!      TOPBC = top_swirl_amplitude*cbc_top(z,x,1)
!      IF (bl_grid > 0) TOPBC = zero
!      BOTBC = bot_swirl_amplitude*cbc_bot(z,x,1)
! AB
      TOPBC = cbc_top(z,x,1)
      BOTBC = cbc_bot(z,x,1)
!      write(*,*) "me in pass2=", me, cbc_top(2:3,2:3,1)
!      write(*,*) maxval(abs(cbc_top(:,:,1))), maxval(abs(cbc_top(:,:,2))), maxval(abs(cbc_top(:,:,3)))
 !     write(*,*) maxval(abs(cbc_bot(:,:,1))), maxval(abs(cbc_bot(:,:,2))), maxval(abs(cbc_bot(:,:,3)))

!stop
      CALL SOLVE_HELMHOLTZ(NY+1,1,cv(1,1),ONE,VIS_IMP,D2H,&
           TOPBC,BOTBC,bl_grid,0)

!      TOPBC = top_swirl_amplitude*cbc_top(z,x,3)
!      IF (bl_grid > 0) TOPBC = zero
!      BOTBC = bot_swirl_amplitude*cbc_bot(z,x,3)
      TOPBC = cbc_top(z,x,3)
      BOTBC = cbc_bot(z,x,3)
      CALL SOLVE_HELMHOLTZ(NY+1,1,cv(1,3),ONE,VIS_IMP,D2H,&
           TOPBC,BOTBC,bl_grid,0)

!  SOLVE FOR MEAN PRESSURE AND SET MEAN WALL-NORMAL VELOCITY.
!!$      TMP = ONE/RK_PRESS
!!$      CP(1) = ZERO
!!$      cp(2) = cp(1) + tmp*DBLE(cv(1,2))*hn(1) &
!!$           + tmp*DBLE(cv(2,2)-cv(1,2))*half*hn(1)**2/hh(2)
!!$      DO PLANE = 2,NY-1
!!$         CP(PLANE+1) = CP(PLANE) + TMP*DBLE(cv(PLANE,2))*HN(PLANE)
!!$      END DO
!!$      cp(ny+1) = cp(ny) + tmp*DBLE(cv(ny,2))*hn(ny) &
!!$           - tmp*DBLE(cv(ny,2)-cv(ny-1,2))*half*hn(ny)**2/hh(ny)
!!$      tmp = SUM(dyh(1:ny+1)*dble(cp(1:ny+1)))/SUM(dyh(1:ny+1))
!!$      cp = cp - cmplx(tmp,zero)
      cp(1:ny+1) = zero
!  SET BOTTOM BC FOR MEAN WALL-NORMAL VELOCITY.
!!$      cv(1,2) = bot_jet_amplitude*cbc_bot(z,x,2)
      cv(1,2) = cbc_bot(z,x,2)
      cv(2:ny,2) = cv(1,2)
      cv(ny+1,2) = zero ! ZERO OUT UNUSED VALUE IN ARRAY FOR V.

    ELSE

!  SOLVE EQUATIONS FOR NON-ZERO WAVENUMBERS.
!  ADVANCE MOMENTUM EQUATION FOR U AND W.
       ALPHA = ONE + VIS_IMP*K2_ME(Z,X)
!  AB
!       TOPBC = top_swirl_amplitude*cbc_top(z,x,1)
!       IF (bl_grid > 0) TOPBC = zero
!       BOTBC = bot_swirl_amplitude*cbc_bot(z,x,1)
!       IF ((x==2).and.(z==1).and.(read_vel == 0).and.(which_test == 6)) THEN
!          TOPBC = cmplx(-half,zero)*exp(-two*pi**2*tmp_time/re)
!          BOTBC = cmplx( half,zero)*exp(-two*pi**2*tmp_time/re)
!       END IF
      TOPBC = cbc_top(z,x,1)
      BOTBC = cbc_bot(z,x,1)
!      write(*,*) maxval(abs(cbc_top))
       CALL SOLVE_HELMHOLTZ(NY+1,1,cv(1,1),ALPHA,VIS_IMP,D2H, &
            TOPBC,BOTBC,bl_grid,0)
!       TOPBC = top_swirl_amplitude*cbc_top(z,x,3)
!       IF (bl_grid > 0) TOPBC = zero
!       BOTBC = bot_swirl_amplitude*cbc_bot(z,x,3)
      TOPBC = cbc_top(z,x,3)
      BOTBC = cbc_bot(z,x,3)
       CALL SOLVE_HELMHOLTZ(NY+1,1,cv(1,3),ALPHA,VIS_IMP,D2H, &
            TOPBC,BOTBC,bl_grid,0)
!  ADVANCE MOMENTUM EQUATION FOR V.
!!$       TOPBC = top_jet_amplitude*cbc_top(z,x,2)
!!$       IF (bl_grid > 0) TOPBC = zero
!!$       BOTBC = bot_jet_amplitude*cbc_bot(z,x,2)
       TOPBC = cbc_top(z,x,2)
       BOTBC = cbc_bot(z,x,2)
       ALPHA = ONE + VIS_IMP*K2_ME(Z,X)
       CALL SOLVE_HELMHOLTZ(NY,1,cv(1,2),ALPHA,VIS_IMP,D2N,&
            TOPBC,BOTBC,bl_grid,0)
       cv(ny+1,2) = zero ! ZERO OUT UNUSED VALUE IN ARRAY FOR V.

!  REMOVE DIVERGENCE FROM NEW VELOCITY FIELD AND UPDATE PRESSURE.
       CALL strip_divergence(cv,ch,rk_press,x,z)

!  NEW PRESSURE = OLD PRESSURE + (ONE - NU*DT/2*LAPLACIAN)*CORRECTION.
!     SHOULD MAKE PRESSURE FIELD SECOND ORDER ACCURATE IN TIME.
!     COMES FROM BROWN, CORTEZ & MINION, JCP 168, P. 483, 2001.
       ALPHA = ONE + vis_imp*K2_ME(Z,X)
       CALL d2_dy(ny+1,ch,du,d2h)
       cp = cp + alpha*ch - vis_imp*du
    END IF

    IF (RK_STEP .LT. RK_LAST) THEN

!  PUT THE NEW (SOON TO BE OLD) VELOCITY FIELD INTO CRU
!     ALONG WITH A PIECE OF THE OLD RHS 
       COEF2 = DT*RKA(RK_STEP+1)
       cru = cv + coef2*cu

!  ADD THE EXPLICIT PART OF THE VISCOUS TERM FOR THE NEXT SUBSTEP.
       VIS_EXP = (DT/RE)*RKC(RK_STEP+1)
       CALL add_explicit_diffusive_terms(cru,cv)

!  ADD PRESSURE GRADIENT FOR THE NEXT SUBSTEP TO NON-ZERO WAVENUMBERS.
       IF (k2_me(z,x) /= zero) THEN
          RK_PRESS = (RKC(RK_STEP+1)+RKD(RK_STEP+1))*DT
          CALL subtract_pressure_gradient(cru,cp,rk_press,x,z)
       END IF
          
    ELSE
!  IF IT'S THE LAST RK_STEP, PUT THE NEW VELOCITY FIELD
!     INTO BOTH CU AND CRU.
       CRU = CV
    END IF

! PUT THE NEW VELOCITY FIELD INTO CU FOR COMPUTING THE NONLINEAR TERM.
    CU = CV
  CONTAINS
!=====================================================================
    SUBROUTINE add_explicit_diffusive_terms(crhs,cvel)
      IMPLICIT NONE

      COMPLEX(KIND=prec), DIMENSION(1:ny+1,1:3) :: crhs, cvel
      COMPLEX(KIND=prec), DIMENSION(1:ny+1)     :: du

      CALL d2_dy(ny+1,cvel(1,1),du,d2h)
      crhs(:,1)    = crhs(:,1)    + vis_exp*(du(:)   -k2_me(z,x)*cvel(:,1))
      CALL d2_dy(ny,  cvel(1,2),du,d2n)
      crhs(1:ny,2) = crhs(1:ny,2) + vis_exp*(du(1:ny)-k2_me(z,x)*cvel(1:ny,2))
      CALL d2_dy(ny+1,cvel(1,3),du,d2h)
      crhs(:,3)    = crhs(:,3)    + vis_exp*(du(:)   -k2_me(z,x)*cvel(:,3))
    END SUBROUTINE add_explicit_diffusive_terms
  END SUBROUTINE phase1
!=====================================================================
!=========================COMPUTE P GRAD==============================
!=====================================================================
  SUBROUTINE COMPUTE_P_GRAD(CU)
    IMPLICIT NONE

!  THIS SUBROUTINE COMPUTES THE PRESSURE GRADIENT FOR THE ENTIRE
!  Y PLANE OF THE CURRENT POSITION.

!  INCLUDE ALL RELEVANT FILES

    COMPLEX(KIND=prec), INTENT(IN) :: CU(NY+1,3)
    REAL(KIND=prec)                   TMP
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

!  COMPUTE THE PRESSURE GRADIENT
    TMP = ZERO
    DO PLANE = 2,NY
      TMP = TMP + DYH(PLANE)*DBLE(CU(PLANE,1))
    END DO
    P_GRAD = - TMP + HALF*(D_UM_DY(NY) - D_UM_DY(1))/RE

  END SUBROUTINE compute_p_grad

!=====================================================================
!==========================CHECK DIVERGENCE===========================
!=====================================================================
  SUBROUTINE CHECK_DIVERGENCE(CU,p,q)
    IMPLICIT NONE

    INTEGER                            p, q, m, n
    COMPLEX(KIND=prec), INTENT(IN) ::  CU(NY+1,3)
    COMPLEX(KIND=prec)             ::  DV(NY+1), TEMPD, MAX_DIVG
    REAL(KIND=prec)                    POSS_MAX

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    CALL d_dy(ny,cu(1,2),ny+1,dv,dh)
    POSS_MAX = ZERO
    MAX_DIVG = ZERO
    DO j = 2,NY
       TEMPD = ikx_me(p)*cu(j,1) + dv(j) + ikz(q)*cu(j,3)
       IF (ABS(TEMPD) .GT. POSS_MAX) THEN
          POSS_MAX = ABS(TEMPD)
          MAX_DIVG = TEMPD
       END IF
    END DO
    
    DIVG(p,q) = MAX_DIVG

    IF ((p == local_nx) .AND. (q == NZ)  &
         .AND. (RK_STEP == RK_LAST) &
         .AND. (MOD(STEP,OUTPUT_STEP) == 0)) THEN
       WRITE(*,*) 'DIVERGENCE INFORMATION'
       WRITE(*,890) 'X', 'Z', 'MAX DIVERGENCE'
       DO n = 1,Mlim
          DO m = 1,Klim
             WRITE(*,900) m+local_x_start, n, DIVG(m,n)
          END DO
       END DO
       
    END IF

    890 FORMAT(2X,A1,3X,A1,10X,A15)
    900 FORMAT(2I4,' (',E13.6,',',E13.6,')')

  END SUBROUTINE check_divergence
!=================================================================
!=======================INITIALIZE JET BC'S=======================
!=================================================================
  SUBROUTINE jet_profile(time,time_incr,ctop,cbot)
    IMPLICIT NONE

    REAL(KIND=prec) time, time_incr
    COMPLEX(KIND=prec), DIMENSION(nz,local_nx) :: ctop,cbot

    REAL(KIND=prec)  xstart, xrise, xfall, xend, &
         frac_suct, suct_start, suct_rise, suct_fall

    REAL(KIND=prec), DIMENSION(nx,local_nz_small) :: vjet_top, vjet_bot
    COMPLEX(KIND=prec), DIMENSION(nz,local_nx) :: csuction,cvjet_top,cvjet_bot

    REAL(KIND=prec), DIMENSION(nx,local_nz_small) :: temp
    COMPLEX(KIND=prec), DIMENSION(nz,local_nx)    :: ctemp

    REAL(KIND=prec), DIMENSION(3,2,4) :: forcing_incr_sign
    REAL(KIND=prec), DIMENSION(nx) :: f
    REAL(KIND=prec)  x2, sum_suction, sum_top, sum_bot, tmp

    INTEGER  xi, xf, x, z, n

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! JET IS DEFINED IN THIS ROUTINE.  SMALL FORCING IS ADDED TO JET OUTFLOW
! VELOCITY TO PROMOTE BREAKUP OF JET.  THIS FORCING FOLLOWS THE DESIGN
! OF FREUND 2001, JFM 438, p. 281.
!     v = (1 - tanh(b*(r/r_0 - r_0/r)))
!     where b(theta,t) = 5 + \sum_m \sum_n A(m,n)*cos(St(m,n)*t + phi(m,n))
!                                                *cos(m*theta   + psi(m,n))
! HERE, m = 0..2 and n = 0..1.
! IF START OF RUN, INITIALIZE FORCING PARAMETERS.
    IF (me == 0) THEN
       tmp = one
       CALL RANDOM_NUMBER(forcing_incr_sign)
       jet_forcing_incr = jet_forcing_incr*SIGN(tmp,forcing_incr_sign-0.05d0)
       jet_forcing = jet_forcing + time_incr*jet_forcing_incr
       DO j = 1,2
          DO i = 1,3
             IF (jet_forcing(i,j,1) <= 0.01d0) THEN
                jet_forcing(i,j,1) = 0.01d0
                jet_forcing_incr(i,j,1) = ABS(jet_forcing_incr(i,j,1))
             ELSEIF (jet_forcing(i,j,1) >= 0.07d0) THEN
                jet_forcing(i,j,1) = 0.07d0
                jet_forcing_incr(i,j,1) = -ABS(jet_forcing_incr(i,j,1))
             END IF
             IF (jet_forcing(i,j,2) <= 0.1d0) THEN
                jet_forcing(i,j,2) = 0.1d0
                jet_forcing_incr(i,j,2) = ABS(jet_forcing_incr(i,j,2))
             ELSEIF (jet_forcing(i,j,2) >= 0.7d0) THEN
                jet_forcing(i,j,2) = 0.7d0
                jet_forcing_incr(i,j,2) = -ABS(jet_forcing_incr(i,j,2))
             END IF
          END DO
       END DO
    END IF
    CALL MPI_BCAST(jet_forcing,24,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, &
         mpierr)

! SET THE POSITION OF THE SUCTION WHICH BALANCES THE MASS FLUX OF THE JET
! IN A CLOSED CHANNEL.
    IF (bl_grid == 0) THEN
       CALL define_suction_profile
       CALL small_fft_real_to_complex(temp,csuction)
       DO n = 1,filter_jet_profile
          CALL filter_cosine(csuction)
       END DO
       IF (k2_me(1,1) .eq. zero) sum_suction = DBLE(csuction(1,1))
       CALL MPI_BCAST(sum_suction,1,MPI_DOUBLE_PRECISION, &
               0,MPI_COMM_WORLD,mpierr)
    END IF

    IF (jetmag_top /= zero) THEN
       !  define the jet profile for the top wall
       CALL define_jet_profile(xjet_top,zjet_top,jetmag_top,forcing_top, &
            for_freq_top,vjet_top)
       ! JET OUTFLOW
       CALL small_fft_real_to_complex(vjet_top,ctemp)
       cvjet_top = - ctemp
       DO n = 1,filter_jet_profile
          CALL filter_cosine(cvjet_top)
       END DO
       ! IF THE DOMAIN IS A CLOSED CHANNEL, ADD SUCTION TO JET PROFILE
       ! SO THAT THERE IS NO NET MASS FLOW ACROSS THE WALL.    
       IF (bl_grid == 0) THEN
          IF (k2_me(1,1) .eq. zero) sum_top = dble(cvjet_top(1,1))
          CALL MPI_BCAST(sum_top,1,MPI_DOUBLE_PRECISION, &
               0,MPI_COMM_WORLD,mpierr)
          cvjet_top = cvjet_top - (sum_top/sum_suction)*csuction
       END IF
       ctop = cvjet_top
    ELSE
       ctop = zero
    END IF

    IF (jetmag_bot /= zero) THEN
       !  DEFINE THE JET PROFILE FOR THE BOTTOM WALL
       CALL define_jet_profile(xjet_bot,zjet_bot,jetmag_bot,forcing_bot, &
            for_freq_bot,vjet_bot)
       ! JET OUTFLOW
       CALL small_fft_real_to_complex(vjet_bot,cvjet_bot)
       DO n = 1,filter_jet_profile
          CALL filter_cosine(cvjet_bot)
       END DO
       ! IF THE DOMAIN IS A CLOSED CHANNEL, ADD SUCTION TO JET PROFILE
       ! SO THAT THERE IS NO NET MASS FLOW ACROSS THE WALL.    
       IF (bl_grid == 0) THEN
          IF (k2_me(1,1) .eq. zero) sum_bot = DBLE(cvjet_bot(1,1))
          CALL MPI_BCAST(sum_bot,1,MPI_DOUBLE_PRECISION, &
               0,MPI_COMM_WORLD,mpierr)
          cvjet_bot = cvjet_bot - (sum_bot/sum_suction)*csuction
       END IF
       cbot= cvjet_bot
    ELSE
       cbot = zero
    END IF

!!$! SET UP TARGET VELOCITY PROFILE IN BUFFER REGION IF IN CHANNEL FLOW DOMAIN
!!$    XI = INT(NX2*(ONE-FRAC_xBUFF))
!!$    XF = INT(NX2*(ONE-HALF*FRAC_xBUFF))      ! ANOTHER OPTION: XF = NX2
!!$
!!$    IF (BL_GRID .EQ. 0) THEN
!!$       junk = zero
!!$       CALL xz_transform(cvjet_top,vtop_full,fftw_complex_to_real)
!!$       CALL xz_transform(cvjet_bot,vbot_full,fftw_complex_to_real)
!!$       CALL make_suction_target(xi,xf,vtop_full,junk,target_top)
!!$       CALL make_suction_target(xi,xf,junk,vbot_full,target_bot)
!!$    END IF

  CONTAINS
!========================================
    SUBROUTINE define_suction_profile
      IMPLICIT NONE

!!$ DEFINES A ONE-DIMENSIONAL PROFILE FOR THE SUCTION (VARYING IN THE
!!$ X-DIRECTION).  THIS PROFILE IS ZERO IN THE FLOW DOMAIN AND ONE IN
!!$ OR JUST UPSTREAM OF THE FRINGE REGION.  IT IS BASED ON THE PROFILE
!!$ SUGGESTED FOR DEFINING THE FRINGE REGION DAMPING FUNCTION IN:
!!$ NORDSTROM, NORDIN & HENNINGSON, SIAM J. SCI COMP, 1999.

      f(1:nx)  = zero
      IF (bl_grid == 0) THEN
         frac_suct  = 0.5d0*frac_xbuff
         suct_start = one - frac_xbuff
         suct_rise  = 0.25d0*frac_suct
         suct_fall  = 0.25d0*frac_suct

         xstart = lx*suct_start
         xrise  = xstart + lx*suct_rise
         xfall  = xstart + lx*frac_suct - lx*suct_fall
         xend   = xstart + lx*frac_suct

         IF (xend .gt. lx) THEN
            write(*,*) 'Suction must end inside box.'
            write(*,*) 'Reset start_suct and frac_suct ', &
                 'in subroutine INIT_JET'
            STOP 'Stopped in subroutine INIT_JET'
         END IF
!      write(*,*) xstart, xrise, xfall, xend
         DO x = 1,nx
            IF ((xcoord(x) .gt. xstart) .and. (xcoord(x) .lt. xend)) THEN
               IF (xcoord(x) .lt. xrise) THEN
                  x2 = (xcoord(x) - xstart)/(xrise-xstart)
                  f(x) = one/(one + dexp(one/x2 + one/(x2-one)))
               ELSEIF (xcoord(x) .gt. xfall) THEN
                  x2 = (xend - xcoord(x))/(xend - xfall)
                  f(x) = one/(one + dexp(one/x2 + one/(x2-one)))
               ELSE
                  f(x) = one
               END if
            ELSE
               f(x) = zero
            END if
!  make suction stronger near the middle of the buffer region
!     than at the entrance.
!         f(x) = f(x)*half*(one + (xcoord(x)-xstart)/(xend-xstart))
         END DO
      ELSE
         f = zero
      END IF

! return a two-dimensional matrix in physical space that defines the
! suction profile over the whole wall.
      DO z = 1,local_nz_small
         temp(1:nx,z) = f(1:nx)
      END DO

    END SUBROUTINE define_suction_profile
!========================================
    SUBROUTINE define_jet_profile(xjet,zjet,jetmag,forcing,for_freq,vjet)
      IMPLICIT NONE

      REAL(KIND=prec) xjet, zjet, jetmag, forcing, for_freq
      REAL(KIND=prec), DIMENSION(nx,local_nz_small) :: vjet
      REAL(KIND=prec)  tmpr, tmpr2, tmpr3, tmpr4, bjet, &
           xdist1, xdist2, zdist1, zdist2, tmpphi, v_ampl
      REAL(KIND=prec), DIMENSION(3,2), PARAMETER :: &
           forcing_waveno = RESHAPE( (/ 0.d0, 1.d0, 2.d0, 0.d0, 1.d0, 2.d0 /),&
                                     (/ 3, 2 /))

! COMPUTE AMPLITUDE OF JET AND SWIRL VELOCITIES
      v_ampl = jetmag*(one + forcing*SIN(two*pi*time*for_freq*jetmag/jet_diam))
      IF (me == 0) WRITE(*,999) step, rk_step, time, v_ampl
      999 format(2i4,2f12.8)

!  define the jet profile for the top wall
      DO z = 1,local_nz_small
         DO x = 1,nx
            xdist1 = xcoord(x) - xjet
            xdist2 = xcoord(x) - (lx+xjet)
            zdist1 = zcoord(z+local_z_small_start) - zjet
            zdist2 = zcoord(z+local_z_small_start) - (lz+zjet)
! COMPUTE DISTANCE OF CURRENT POINT FROM JET AND ITS PERIODIC EXTENSIONS
            tmpr  = ABS(CMPLX(xdist1,zdist1))*two/jet_diam
            tmpr2 = ABS(CMPLX(xdist2,zdist1))*two/jet_diam
            tmpr3 = ABS(CMPLX(xdist1,zdist2))*two/jet_diam
            tmpr4 = ABS(CMPLX(xdist2,zdist2))*two/jet_diam
! SET UP PROFILE OF JET OUTFLOW VELOCITY.  DEFINED POSITIVE HERE.
! BE SURE TO REVERSE SIGN FOR JET BLOWING INTO DOMAIN FROM TOP WALL.
            tmpr = MIN(tmpr,tmpr2,tmpr3,tmpr4)
            IF (tmpr == zero) THEN
               vjet(x,z) = v_ampl
            ELSE
               ! COMPUTE ANGLE WITH RESPECT TO JET CENTER AND X-AXIS
               tmpphi = atan2(zdist1,xdist1)
               bjet = 5.d0 + SUM(jet_forcing(:,:,1) &
                    *cos(jet_forcing(:,:,2)*jetmag*time/jet_diam &
                        + jet_forcing(:,:,3)) &
                    *cos(forcing_waveno*tmpphi + jet_forcing(:,:,4)))
               vjet(x,z) = v_ampl*half*(one - tanh(bjet*(tmpr - one/tmpr)))
            END IF
         END DO
      END DO

    END SUBROUTINE define_jet_profile
!========================================
    SUBROUTINE filter_cosine(cu_unfiltered)
      IMPLICIT NONE

      COMPLEX(KIND=prec), DIMENSION(nz,local_nx) :: cu_unfiltered
      REAL(KIND=prec) :: xfilter

! APPLIES RAISED COSINE FILTER TO VELOCITY/SCALAR PROFILE AT WALL.
! SEE CANUTO et al 1988, Spectral Methods in Fluid Dynamics, pp. 246--252.
      DO x = 1,local_nx
         cu_unfiltered(1:nz,x) = cu_unfiltered(1:nz,x) &
                 *half*(one + cos(pi*kx_me(x)/(DBLE(nx/2)*two*pi/lx))) &
                 *half*(one + cos(pi*kz(1:nz)/(DBLE(nz/2)*two*pi/lz)))
      END DO

    END SUBROUTINE filter_cosine
  END SUBROUTINE jet_profile

END MODULE pass2

