MODULE init_jet_buffer
  USE runparms
  USE gridparms
  USE data
  USE transform
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_jet, init_buffer

  INTEGER  i, k, n, x, z, xi, xf

CONTAINS
!=================================================================
!=======================INITIALIZE JET BC'S=======================
!=================================================================
  SUBROUTINE INIT_JET
    IMPLICIT NONE

    REAL(KIND=prec)  xstart, xrise, xfall, xend, &
         frac_suct, suct_start, suct_rise, suct_fall

    REAL(KIND=prec), DIMENSION(nx,nz) ::         &
         vjet_top, uswirl_top, wswirl_top, vjet_bot, uswirl_bot, wswirl_bot
    COMPLEX(KIND=prec), DIMENSION(nx/2+1,nz) :: csuction, &
         cujet_top, cvjet_top, cwjet_top, cujet_bot, cvjet_bot, cwjet_bot
    REAL(KIND=prec), DIMENSION(nx,nz) :: tjet_top, tjet_bot

    REAL(KIND=prec), DIMENSION(nx,nz)          :: temp
    COMPLEX(KIND=prec), DIMENSION(nx/2+1,nz)   :: ctemp

    REAL(KIND=prec), DIMENSION(nx2,nz2) :: junk, vtop_full, vbot_full, temp2
    REAL(KIND=prec), DIMENSION(nx) :: tmpx, f
    REAL(KIND=prec), DIMENSION(nz) :: tmpz
    REAL(KIND=prec)  x2
    REAL(KIND=prec)  sum_suction, sum_top, sum_bot
    REAL(KIND=prec), DIMENSION(nx2) :: tmpx2
    REAL(KIND=prec), DIMENSION(nz2) :: tmpz2

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! SET THE LOCATION OF THE JET IN THE X- AND Z-DIRECTIONS.
! THESE VALUES ARE SCALES IN SUBROUTINE read_parms
!!$    xjet_top = lx*xjet_top
!!$    xjet_bot = lx*xjet_bot
!!$    zjet_top = lz*zjet_top
!!$    zjet_bot = lz*zjet_bot

! INITIALIZE tmpx AND tmpz WHICH HOLD THE X- AND Z-COORDINATES.
    do x = 1,nx
       tmpx(x) = DBLE(x-1)*lx/DBLE(nx)
    end do
    do z = 1,nz
       tmpz(z) = DBLE(z-1)*lz/DBLE(nz)
    end do
    do x = 1,nx2
       tmpx2(x) = DBLE(x-1)*lx/DBLE(nx2)
    end do
    do z = 1,nz2
       tmpz2(z) = DBLE(z-1)*lz/DBLE(nz2)
    end do
         
! SET THE POSITION OF THE SUCTION WHICH BALANCES THE MASS FLUX OF THE JET
! IN A CLOSED CHANNEL.
    CALL define_suction_profile
    CALL output_fft_real_to_complex(temp,csuction)
    DO n = 1,filter_jet_profile
       CALL filter_cosine(csuction)
    END DO
    sum_suction = DBLE(csuction(1,1))

!  DEFINE JET BOUNDARY CONDITIONS FOR VELOCITY FIELD.
    IF (jetmag_top /= zero) THEN
       !  define the jet profile for the top wall
       CALL define_jet_velocity_profile(xjet_top,zjet_top,vjet_top, &
            uswirl_top,wswirl_top)
       ! JET OUTFLOW
       CALL output_fft_real_to_complex(vjet_top,ctemp)
       ! REVERSE SENSE OF VELOCITY SO THAT JET EXHAUSTS INTO DOMAIN.
       cvjet_top = - ctemp 
       DO n = 1,filter_jet_profile
          CALL filter_cosine(cvjet_top)
       END DO
       ! IF THE DOMAIN IS A CLOSED CHANNEL, ADD SUCTION TO JET PROFILE
       ! SO THAT THERE IS NO NET MASS FLOW ACROSS THE WALL.    
       sum_top = dble(cvjet_top(1,1))
       IF (bl_grid == 0) cvjet_top = cvjet_top - (sum_top/sum_suction)*csuction
       ! U AND W SWIRL
       CALL output_fft_real_to_complex(uswirl_top,cujet_top)
       DO n = 1,filter_jet_profile
          CALL filter_cosine(cujet_top)
       END DO
       CALL output_fft_real_to_complex(wswirl_top,cwjet_top)
       DO n = 1,filter_jet_profile
          CALL filter_cosine(cwjet_top)
       END DO

       ! PUT THE BOUNDARY CONDITIONS INTO COMPLEX ARRAYS THAT WILL BE USED 
       ! IN DEFINING BOUNDARY CONDITIONS FOR THE VELOCITY AND SCALAR FIELD.
       DO x = 1,local_nx
          cbc_top(1:nz,x,1) = cujet_top(local_x_start+x,1:nz)
          cbc_top(1:nz,x,2) = cvjet_top(local_x_start+x,1:nz)
          cbc_top(1:nz,x,3) = cwjet_top(local_x_start+x,1:nz)
       END DO

       ! TRANSFORM THE BOUNDARY CONDITIONS INTO PHYSICAL SPACE ON THE 
       ! DE-ALIASED GRID AND STORE IN ARRAYS THAT CAN BE USED FOR 
       ! BOUNDARY CONDITIONS WHEN SOLVING FOR VELOCITY/SCALAR PROFILE 
       ! WITH IMPLICIT WALL-NORMAL CONVECTION.
       CALL xz_transform(cujet_top,temp2,fftw_complex_to_real)
       bc_top(1:nx2,1:local_nz,1) = &
            temp2(1:nx2,local_z_start+1:local_z_start+local_nz)
       CALL xz_transform(cvjet_top,temp2,fftw_complex_to_real)
       bc_top(1:nx2,1:local_nz,2) = &
            temp2(1:nx2,local_z_start+1:local_z_start+local_nz)
       CALL xz_transform(cwjet_top,temp2,fftw_complex_to_real)
       bc_top(1:nx2,1:local_nz,3) = &
            temp2(1:nx2,local_z_start+1:local_z_start+local_nz)

       ! SCALAR BOUNDARY CONDITION
       IF (nscalar >= 1) THEN
          CALL define_jet_scalar_profile(xjet_top,zjet_top,tjet_top)
          DO i = 4,3+nscalar
             DO z = 1,local_nz_small
                bc_top(1:nx,z,i) = tjet_top(1:nx,local_z_small_start+z)
             END DO
          END DO
       END IF

    ELSE
       bc_top = zero
       cbc_top = zero
    END IF

    IF (jetmag_bot /= zero) THEN
       !  DEFINE THE JET PROFILE FOR THE BOTTOM WALL
       CALL define_jet_velocity_profile(xjet_bot,zjet_bot,vjet_bot, &
            uswirl_bot,wswirl_bot)
       ! JET OUTFLOW
       CALL output_fft_real_to_complex(vjet_bot,cvjet_bot)
       DO n = 1,filter_jet_profile
          CALL filter_cosine(cvjet_bot)
       END DO
       ! IF THE DOMAIN IS A CLOSED CHANNEL, ADD SUCTION TO JET PROFILE
       ! SO THAT THERE IS NO NET MASS FLOW ACROSS THE WALL.    
       sum_bot = dble(cvjet_bot(1,1))
       IF (bl_grid == 0) cvjet_bot = cvjet_bot - (sum_bot/sum_suction)*csuction
       ! U AND W SWIRL
       CALL output_fft_real_to_complex(uswirl_bot,cujet_bot)
       DO n = 1,filter_jet_profile
          CALL filter_cosine(cujet_bot)
       END DO
       CALL output_fft_real_to_complex(wswirl_bot,cwjet_bot)
       DO n = 1,filter_jet_profile
          CALL filter_cosine(cwjet_bot)
       END DO

       ! PUT THE BOUNDARY CONDITIONS INTO COMPLEX ARRAYS THAT WILL BE USED 
       ! IN DEFINING BOUNDARY CONDITIONS FOR THE VELOCITY AND SCALAR FIELD.
       DO x = 1,local_nx
          cbc_bot(1:nz,x,1) = cujet_bot(local_x_start+x,1:nz)
          cbc_bot(1:nz,x,2) = cvjet_bot(local_x_start+x,1:nz)
          cbc_bot(1:nz,x,3) = cwjet_bot(local_x_start+x,1:nz)
       END DO

       ! TRANSFORM THE BOUNDARY CONDITIONS INTO PHYSICAL SPACE ON THE DE-ALIASED
       ! GRID AND STORE IN ARRAYS THAT CAN BE USED FOR BOUNDARY CONDITIONS
       ! WHEN SOLVING FOR VELOCITY/SCALAR PROFILE WITH IMPLICIT WALL-NORMAL 
       ! CONVECTION.
       CALL xz_transform(cujet_bot,temp2,fftw_complex_to_real)
       bc_bot(1:nx2,1:local_nz,1) = &
            temp2(1:nx2,local_z_start+1:local_z_start+local_nz)
       CALL xz_transform(cvjet_bot,temp2,fftw_complex_to_real)
       bc_bot(1:nx2,1:local_nz,2) = &
            temp2(1:nx2,local_z_start+1:local_z_start+local_nz)
       CALL xz_transform(cwjet_bot,temp2,fftw_complex_to_real)
       bc_bot(1:nx2,1:local_nz,3) = &
            temp2(1:nx2,local_z_start+1:local_z_start+local_nz)

       ! SCALAR BOUNDARY CONDITION
       IF (nscalar >= 1) THEN
          CALL define_jet_scalar_profile(xjet_bot,zjet_bot,tjet_bot)
          DO i = 4,3+nscalar
             DO z = 1,local_nz_small
                bc_bot(1:nx,z,i) = tjet_bot(1:nx,local_z_small_start+z)
             END DO
          END DO
       END IF
    ELSE
       bc_bot = zero
       cbc_bot = zero
    END IF

! SET UP TARGET VELOCITY PROFILE IN BUFFER REGION IF IN CHANNEL FLOW DOMAIN
    XI = INT(NX2*(ONE-FRAC_xBUFF))
    XF = INT(NX2*(ONE-HALF*FRAC_xBUFF))      ! ANOTHER OPTION: XF = NX2

    IF (BL_GRID .EQ. 0) THEN
       junk = zero
       CALL xz_transform(cvjet_top,vtop_full,fftw_complex_to_real)
       CALL xz_transform(cvjet_bot,vbot_full,fftw_complex_to_real)
       CALL make_suction_target(xi,xf,vtop_full,junk,target_top)
       CALL make_suction_target(xi,xf,junk,vbot_full,target_bot)
    END IF

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
            IF ((tmpx(x) .gt. xstart) .and. (tmpx(x) .lt. xend)) THEN
               IF (tmpx(x) .lt. xrise) THEN
                  x2 = (tmpx(x) - xstart)/(xrise-xstart)
                  f(x) = one/(one + dexp(one/x2 + one/(x2-one)))
               ELSEIF (tmpx(x) .gt. xfall) THEN
                  x2 = (xend - tmpx(x))/(xend - xfall)
                  f(x) = one/(one + dexp(one/x2 + one/(x2-one)))
               ELSE
                  f(x) = one
               END if
            ELSE
               f(x) = zero
            END if
!  make suction stronger near the middle of the buffer region
!     than at the entrance.
!         f(x) = f(x)*half*(one + (tmpx(x)-xstart)/(xend-xstart))
         END DO
      ELSE
         f = zero
      END IF

! return a two-dimensional matrix in physical space that defines the
! suction profile over the whole wall.
      DO z = 1,nz
         temp(1:nx,z) = f(1:nx)
      END DO

    END SUBROUTINE define_suction_profile
!========================================
    SUBROUTINE define_jet_velocity_profile(xjet,zjet,vjet,uswirl,wswirl)
      IMPLICIT NONE

      REAL(KIND=prec) xjet, zjet, xdist1, xdist2, zdist1, zdist2, tmpphi
      REAL(KIND=prec), DIMENSION(nx,nz) :: vjet, uswirl, wswirl
      REAL(KIND=prec)  tmpr, tmpr2, tmpr3, tmpr4

!  define the jet profile for the top wall
      DO z = 1,nz
         DO x = 1,nx
            xdist1 = tmpx(x) - xjet
            xdist2 = tmpx(x) - (lx+xjet)
            zdist1 = tmpz(z) - zjet
            zdist2 = tmpz(z) - (lz+zjet)
! COMPUTE DISTANCE OF CURRENT POINT FROM JET AND ITS PERIODIC EXTENSIONS
            tmpr  = SQRT(xdist1**2 + zdist1**2)*two/jet_diam
            tmpr2 = SQRT(xdist2**2 + zdist1**2)*two/jet_diam
            tmpr3 = SQRT(xdist1**2 + zdist2**2)*two/jet_diam
            tmpr4 = SQRT(xdist2**2 + zdist2**2)*two/jet_diam
! SET UP PROFILE OF JET OUTFLOW VELOCITY.  DEFINED POSITIVE HERE.
! BE SURE TO REVERSE SIGN FOR JET BLOWING INTO DOMAIN FROM TOP WALL.
!!$! SET UP PROFILE OF PASSIVE SCALAR ON TOP WALL.  
            tmpr = min(tmpr,tmpr2,tmpr3,tmpr4)
            IF (tmpr == zero) THEN
               vjet(x,z) = one
               uswirl(x,z) = zero
               wswirl(x,z) = zero
            ELSE
               vjet(x,z) = half*(one - tanh(5.d0*(tmpr - one/tmpr)))
               ! COMPUTE ANGLE OF CURRENT POINT WITH RESPECT TO CENTER 
               ! OF JET AND AND X-AXIS.
               tmpphi = atan2(tmpz(z)-zjet,tmpx(x)-xjet)
               uswirl(x,z) = -tmpr*sin(tmpphi)*vjet(x,z)
               wswirl(x,z) =  tmpr*cos(tmpphi)*vjet(x,z)
            END IF
         END DO
      END DO

    END SUBROUTINE define_jet_velocity_profile
!========================================
    SUBROUTINE define_jet_scalar_profile(xjet,zjet,tjet)
      IMPLICIT NONE

      REAL(KIND=prec) xjet, zjet, xdist1, xdist2, zdist1, zdist2
      REAL(KIND=prec), DIMENSION(nx,nz) :: tjet
      REAL(KIND=prec)  tmpr, tmpr2, tmpr3, tmpr4

!  define the jet profile for the top wall
      DO z = 1,nz
         DO x = 1,nx
            xdist1 = tmpx(x) - xjet
            xdist2 = tmpx(x) - (lx+xjet)
            zdist1 = tmpz(z) - zjet
            zdist2 = tmpz(z) - (lz+zjet)
! COMPUTE DISTANCE OF CURRENT POINT FROM JET AND ITS PERIODIC EXTENSIONS
            tmpr  = SQRT(xdist1**2 + zdist1**2)*two/jet_diam
            tmpr2 = SQRT(xdist2**2 + zdist1**2)*two/jet_diam
            tmpr3 = SQRT(xdist1**2 + zdist2**2)*two/jet_diam
            tmpr4 = SQRT(xdist2**2 + zdist2**2)*two/jet_diam
! SET UP PROFILE OF JET OUTFLOW VELOCITY.  DEFINED POSITIVE HERE.
! BE SURE TO REVERSE SIGN FOR JET BLOWING INTO DOMAIN FROM TOP WALL.
!!$! SET UP PROFILE OF PASSIVE SCALAR ON TOP WALL.  
            tmpr = min(tmpr,tmpr2,tmpr3,tmpr4)
            IF (tmpr == zero) THEN
               tjet(x,z) = one
            ELSE
               tjet(x,z) = half*(one - tanh(5.d0*(tmpr - one/tmpr)))
            END IF
         END DO
      END DO

    END SUBROUTINE define_jet_scalar_profile
!========================================
    SUBROUTINE filter_cosine(cu_unfiltered)
      IMPLICIT NONE

      COMPLEX(KIND=prec), DIMENSION(nx/2+1,nz) :: cu_unfiltered

! APPLIES RAISED COSINE FILTER TO VELOCITY/SCALAR PROFILE AT WALL.
! SEE CANUTO et al 1988, Spectral Methods in Fluid Dynamics, pp. 246--252.
      DO z = 1,nz
         DO x = 1,nx/2+1
            cu_unfiltered(x,z) = cu_unfiltered(x,z) &
                 *half*(one + cos(pi*kx(x)/(DBLE(nx/2)*two*pi/lx)))  &
                 *half*(one + cos(pi*kz(z)/(DBLE(nz/2)*two*pi/lz)))
         END DO
      END DO

    END SUBROUTINE filter_cosine
  END SUBROUTINE init_jet
!=================================================================
!===============INITIALIZE BUFFER/FRINGE REGION===================
!=================================================================
  SUBROUTINE INIT_BUFFER
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(nx2,nz2) :: buffer_full
    REAL(KIND=prec), DIMENSION(nx)      :: tmpx, fx, fx2
    REAL(KIND=prec), DIMENSION(nz)      :: tmpz, fz, fza, fzb, fz2, fza2, fzb2
    REAL(KIND=prec) xstart, xrise, xfall, xend, zstart, zrise, zfall, zend, tmp
    
    REAL(KIND=prec)    :: temp(nx,nz)
    COMPLEX(KIND=prec) :: ctemp(nx/2+1,nz)

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    do x = 1,nx
       tmpx(x) = DBLE(x-1)*lx/DBLE(nx)
    end do
         
    do z = 1,nz
       tmpz(z) = DBLE(z-1)*lz/DBLE(nz)
    end do
         
    IF (xfringe == 1) THEN
       xstart = lx*(one-frac_xbuff)
       xrise  = xstart + lx*buff_xrise
       xfall  = lx  - lx*buff_xfall
       xend   = lx
       CALL define_fringe_region(nx,fx,fx2,tmpx,xstart,xrise,xfall,xend)
    ELSE
       fx = zero
       fx2 = zero
    END IF
   
    IF (zfringe == 1) THEN
       zstart = -one
       zrise  = -one
       zfall  = zero
       zend   = frac_zbuff*lz
       CALL define_fringe_region(nz,fza,fza2,tmpz,zstart,zrise,zfall,zend)

       zstart = lz - frac_zbuff*lz
       zrise  = lz 
       zfall  = lz + one
       zend   = lz + one
       CALL define_fringe_region(nz,fzb,fzb2,tmpz,zstart,zrise,zfall,zend)

! THIS COMPOSITION OF THE TWO FRINGE REGIONS IN Z RELIES ON THEM BEING DISJOINT
       fz = fza + fzb
       fz2 = fza2 + fzb2
    ELSE
       fz = zero
       fz2 = zero
    END IF

! COMBINE THE FRINGE REGIONS IN x AND z INTO A SMOOTH AND CONTINUOUS PROFILE
    DO z = 1,nz
       DO x = 1,nx
          IF ((fx(x)<one).and.(fz(z)<one).and.(fx(x)*fz(z) > zero)) THEN
             tmp = sqrt(fx2(x)**2 + fz2(z)**2)
             IF (tmp < one) THEN
                temp(x,z) = one/(one + exp(one/tmp + one/(tmp-one)))
             ELSE
                temp(x,z) = one
             END IF
          ELSE
             temp(x,z) = max(fx(x),fz(z))
          END IF
       END DO
    END DO
    DO z = 1,local_nz_small
       buffer_small(1:nx,z) = fringe_gain*temp(1:nx,local_z_small_start+z)
    END DO

    CALL output_fft_real_to_complex(temp,ctemp)
    CALL xz_transform(ctemp,buffer_full,fftw_complex_to_real)
    DO z = 1,local_nz
       buffer(1:nx2,z) = fringe_gain*buffer_full(1:nx2,local_z_start+z)
    END DO
    target_inflow = zero    
    IF (bl_grid .gt. 0) THEN
       target_inflow(1:ny+1,1) = ublasius(1:ny+1)
    ELSE
       target_inflow(1:ny+1,1) = 0.75d0*(one - yh(1:ny+1)**8)*COS(beta)
       target_inflow(1:ny+1,3) = 0.75d0*(one - yh(1:ny+1)**8)*SIN(beta)
    END IF

!!$    IF (me == 0) THEN
!!$       DO z = 1,ny+1
!!$          WRITE(*,999) (target_inflow(z,x),x=1,3+nscalar)
!!$       END DO
!!$       WRITE(*,*)
!!$    END IF
999 FORMAT(4e12.4)

  CONTAINS
!====================================================================
    SUBROUTINE define_fringe_region(n,f,f2,t,tstart,trise,tfall,tend)
      IMPLICIT NONE
      
      INTEGER :: n
      REAL(KIND=prec), DIMENSION(n) :: f, f2, t
      REAL(KIND=prec)               :: tstart, trise, tfall, tend, t2

      f = zero
      f2 = zero
      DO i = 1,n
         IF (t(i) > tstart) THEN
            IF (t(i) < trise) THEN
                t2 = (t(i) - tstart)/(trise-tstart)
                f(i) = one/(one + dexp(one/t2 + one/(t2-one)))
                f2(i) = t2
             ELSEIF ((t(i) > tfall) .and. (t(i) < tend)) THEN
                t2 = (tend - t(i))/(tend - tfall)
                f(i) = one/(one + dexp(one/t2 + one/(t2-one)))
                f2(i) = t2
             ELSEIF (t(i) >= tend) THEN
                f(i) = zero
             ELSE
                f(i) = one
             END IF
          ELSE
             f(i) = zero
          END IF
       END DO

     END SUBROUTINE define_fringe_region
  END SUBROUTINE INIT_BUFFER
!=================================================================
!===============INITIALIZE BUFFER/FRINGE REGION===================
!=================================================================
  SUBROUTINE MAKE_SUCTION_TARGET(XI, XF, VTOPBC_IN, VBOTBC_IN, VEL_TARGET)
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(nx2,nz2), INTENT(IN)     :: VTOPBC_IN, VBOTBC_IN
    REAL(KIND=prec), DIMENSION(nx2,ny+1,2), INTENT(OUT) :: VEL_TARGET

    REAL(KIND=prec)    xi_top, xi_bot
    INTEGER :: xi, xf

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! INITIALIZE TARGET VELOCITY IN FRINGE REGION TO ZERO.
    VEL_TARGET = zero

    DO k = 1,ny
       DO x = xf-1,xi,-1
          xi_top = (y(k)-y(ny))/(y(1)-y(ny))
          xi_bot = (y(k)-y(1)) /(y(ny)-y(1))
          vel_target(x,k,2) = vtopbc_in(x,1)*half*(one + cos(pi*xi_top)) &
               + vbotbc_in(x,1)*half*(one + cos(pi*xi_bot))
       END DO
    END DO

    vel_target(xi:xf,1,1) = zero
    DO k = 2,ny
       DO x = xf-1,xi,-1
          vel_target(x,k,1) = vel_target(x+1,k,1) &
               + half*delta1*dh(k)*(vel_target(x+1,k,2)-vel_target(x+1,k-1,2))&
               + half*delta1*dh(k)*(vel_target(x,k,2)-vel_target(x,k-1,2))
       END DO
    END DO
    vel_target(xi:xf,ny+1,1) = zero

  END SUBROUTINE MAKE_SUCTION_TARGET
END MODULE init_jet_buffer

!!$            vjet(x,z) =   0.99d0*(EXP(-tmpr**8)   &
!!$                 + EXP(-tmpr2**8)  &
!!$                 + EXP(-tmpr3**8)  &
!!$                 + EXP(-tmpr4**8)) &
!!$                 + 0.01D0*(EXP(-half*tmpr**2)     &
!!$                 + EXP(-half*tmpr2**2)    &
!!$                 + EXP(-half*tmpr3**2)    &
!!$                 + EXP(-half*tmpr4**2)) 
!!$            tjet(x,z) =  0.99d0*(EXP(-(0.8d0*tmpr)**8)   &
!!$                 + EXP(-(0.8d0*tmpr2)**8)  &
!!$                 + EXP(-(0.8d0*tmpr3)**8)  &
!!$                 + EXP(-(0.8d0*tmpr4)**8)) &
!!$                 + 0.01D0*(EXP(-half*(0.8d0*tmpr)**2)  &
!!$                 + EXP(-half*(0.8d0*tmpr2)**2) &
!!$                 + EXP(-half*(0.8d0*tmpr3)**2) &
!!$                 + EXP(-half*(0.8d0*tmpr4)**2))

