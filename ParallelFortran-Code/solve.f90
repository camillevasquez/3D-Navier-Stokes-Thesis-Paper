MODULE solve
  USE runparms
  USE gridparms
  USE diff_int
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: solve_helmholtz, solve_real_helmholtz, solve_p, &
       strip_divergence, subtract_pressure_gradient

CONTAINS
!=====================================================================
!==========================STRIP DIVERGENCE===========================
!=====================================================================
  SUBROUTINE strip_divergence(cu,cp,coef,p,q)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(INOUT) :: cu(1:ny+1,3)
    COMPLEX(KIND=prec), INTENT(OUT)   :: cp(1:ny+1)
    REAL(KIND=prec), INTENT(IN)       :: coef
    INTEGER, INTENT(IN)               :: p, q  ! p, q are wavenumbers in x & z

    REAL(KIND=prec)            :: tmp, alpha
    COMPLEX(KIND=prec)         :: dv(1:ny+1), topbc, botbc
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! compute divergence of cu.  scale by a factor of 1/coef.
    tmp = one/coef
    CALL d_dy(ny,cu(1,2),ny+1,dv,dh)
    cp(1)    = zero
    cp(2:ny) = tmp*(ikx_me(p)*cu(2:ny,1) + dv(2:ny) + ikz(q)*cu(2:ny,3))
    cp(ny+1) = zero

! solve for the pressure correction (cp) that will make cu divergence-free.
    topbc = zero
    botbc = zero
    alpha = -k2_me(q,p)
    CALL solve_p(cp,alpha,topbc,botbc)

! subtract the gradient of cp (scaled by coef) from cu to make it div-free.
    CALL subtract_pressure_gradient(cu,cp,coef,p,q)

!!$    IF ((p==2).and.(q==1)) THEN
!!$       CALL d_dy(ny,cu(1,2),ny+1,dv,dh)
!!$       DO i = 2,ny
!!$          write(*,995) cu(i,1), dv(i), ikx_me(p)*cu(i,1) + dv(i) 
!!$       END DO
!!$       write(*,*)
!!$       995 FORMAT(6e12.4)
!!$    END IF
  END SUBROUTINE strip_divergence
!=====================================================================
!==================SUBTRACT PRESSURE GRADIENT=========================
!=====================================================================
  SUBROUTINE subtract_pressure_gradient(crhs,cpress,coef,p,q)
    IMPLICIT NONE

    COMPLEX(KIND=prec), DIMENSION(1:ny+1,1:3), INTENT(INOUT) :: crhs
    COMPLEX(KIND=prec), DIMENSION(1:ny+1), INTENT(IN)        :: cpress
    REAL(KIND=prec), INTENT(IN)                              :: coef
    INTEGER, INTENT(IN) :: p, q   ! p AND q INDEX THE X- AND Z-WAVENUMBERS

    COMPLEX(KIND=prec), DIMENSION(1:ny) :: dp

!!$    CALL d_dy(ny+1,cpress,ny,dp,dn)
!!$    crhs(1:ny+1,1) = crhs(1:ny+1,1) - coef*ikx_me(p)*cpress(1:ny+1)
!!$!    crhs(1:ny,2)   = crhs(1:ny,2)   - coef*dp(1:ny)
!!$    crhs(1:ny+1,3) = crhs(1:ny+1,3) - coef*ikz(q)*cpress(1:ny+1)
!!$    crhs(1,2) = crhs(1,2) - coef*SUM(dm(1:3,1)*cpress(1:3))
!!$    crhs(2:ny-1,2) = crhs(2:ny-1,2) - coef*dp(2:ny-1)
!!$    crhs(ny,2) = crhs(ny,2) - coef*SUM(dm(1:3,2)*cpress(ny-1:ny+1))

    CALL d_dy(ny+1,cpress,ny,dp,dn)
    crhs(2:ny,1)   = crhs(2:ny,1)   - coef*ikx_me(p)*cpress(2:ny)
    crhs(2:ny-1,2) = crhs(2:ny-1,2) - coef*dp(2:ny-1)
    crhs(2:ny,3)   = crhs(2:ny,3)   - coef*ikz(q)*cpress(2:ny)

  END SUBROUTINE subtract_pressure_gradient
!=====================================================================
!=======================SOLVE REAL HELMHOLTZ==========================
!=====================================================================
  SUBROUTINE SOLVE_REAL_HELMHOLTZ(N, NRHS, F, ALPHA, DIFF, D2, &
          TOPBC, BOTBC, TOP_NEUMANN, BOT_NEUMANN)
    IMPLICIT NONE

    INTEGER, INTENT(IN)            :: N, NRHS
    REAL(KIND=prec), INTENT(INOUT) :: F(N,NRHS)

    INTEGER, INTENT(IN)            :: TOP_NEUMANN, BOT_NEUMANN
    REAL(KIND=prec), INTENT(IN)    :: ALPHA, DIFF, D2(3,N)
    REAL(KIND=prec), INTENT(IN)    :: TOPBC(NRHS), BOTBC(NRHS)

    INTEGER            I
    REAL(KIND=prec)    A(N), B(N), C(N)
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0


!  SET UP THE LEFT HAND SIDE MATRIX 
!     (DIFF*D2 - ALPHA*I)

! SET BOTTOM BOUNDARY CONDITION.
    IF (BOT_NEUMANN .EQ. 0) THEN
       A(1) = ZERO
       B(1) = ONE
       C(1) = ZERO
    ELSE
       A(1) = ZERO
       B(1) = - ONE
       C(1) =   ONE
    END IF

    a(2:n-1) = - diff*d2(1,2:n-1) 
    b(2:n-1) = - diff*d2(2,2:n-1) + alpha
    c(2:n-1) = - diff*d2(3,2:n-1) 

!  SET UP TOP BOUNDARY CONDITION.
    IF (TOP_NEUMANN .EQ. 0) THEN
       A(N) = ZERO
       B(N) = ONE
       C(N) = ZERO
    ELSE
       A(N) = - ONE
       B(N) =   ONE
       C(N) =  ZERO
    END IF

!  APPLY BOUNDARY CONDITIONS ON RHS
    F(1,1:nrhs) = BOTBC(1:nrhs)
    F(N,1:nrhs) = TOPBC(1:nrhs)

!  SOLVE TRIDIAGONAL SYSTEM USING GAUSSIAN ELIMINATION.
!  FORWARD ELIMINATION
    DO I = 2,N
       B(I-1) = B(I-1)**(-1)
       C(I-1) = C(I-1)*B(I-1)
       F(I-1,1:nrhs) = F(I-1,1:nrhs)*B(I-1)

       B(I) = B(I) - A(I)*C(I-1)
       F(I,1:nrhs) = F(I,1:nrhs) - A(I)*F(I-1,1:nrhs)
    END DO

! SOLVE FOR THE LAST ELEMENT IN THE VECTOR
    B(N) = B(N)**(-1)
    F(N,1:nrhs) = F(N,1:nrhs)*B(N)

! BACK SUBSTITUTION.
    DO I = N-1,1,-1
       F(I,1:nrhs) = F(I,1:nrhs) - C(I)*F(I+1,1:nrhs)
    END DO

  END SUBROUTINE solve_real_helmholtz
!=====================================================================
!==========================SOLVE HELMHOLTZ============================
!=====================================================================
  SUBROUTINE SOLVE_HELMHOLTZ(N, NRHS, F, ALPHA, DIFF, D2, &
          TOPBC, BOTBC, TOP_NEUMANN, BOT_NEUMANN)
    IMPLICIT NONE

    INTEGER, INTENT(IN)               :: N, NRHS
    COMPLEX(KIND=prec), INTENT(INOUT) :: F(N,NRHS)

    INTEGER, INTENT(IN)               :: TOP_NEUMANN, BOT_NEUMANN
    REAL(KIND=prec), INTENT(IN)       :: ALPHA, DIFF, D2(3,N)
    COMPLEX(KIND=prec), INTENT(IN)    :: TOPBC(NRHS), BOTBC(NRHS)

    INTEGER            I
    REAL(KIND=prec)    A(N), B(N), C(N)
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0


!  SET UP THE LEFT HAND SIDE MATRIX 
!     (DIFF*D2 - ALPHA*I)

! SET BOTTOM BOUNDARY CONDITION.
    IF (BOT_NEUMANN .EQ. 0) THEN
       A(1) = ZERO
       B(1) = ONE
       C(1) = ZERO
       f(1,:) = botbc(:)
    ELSE
       IF (n == ny+1) THEN
! MODIFY SECOND EQUATION TO INCORPORATE BOUNDARY CONDITION IN VISCOUS FLUX
          a(2) =       zero
          b(2) = one + diff*dh(2)*dn(2)
          c(2) =     - diff*dh(2)*dn(2)
          f(2,:) = f(2,:) - diff*dh(2)*botbc(:)

! LET FIRST ROW OF MATRIX HOLD NEUMANN BOUNDARY CONDITION.
! THIS ALLOWS THE APPLICATION OF A SECOND-ORDER ACCURATE NEUMANN BC.
! USE THE SECOND EQUATION TO ELIMINATE THE THIRD ENTRY IN THE FIRST ROW.
          a(1) = dm(3,1)/c(2)
          b(1) = dm(1,1) - a(1)*a(2)
          c(1) = dm(2,1) - a(1)*b(2)
          f(1,:) = botbc(:)   - a(1)*f(2,:)
       ELSEIF (n == ny) THEN
! USE SECOND-ORDER ACCURATE ONE-SIDED DIFFERENCE FORMULA FOR NEUMANN BC.
! ELIMINATE THIRD ENTRY IN TOP ROW USING SECOND EQUATION.
          A(1) = dm(3,3)/c(2)
          B(1) = dm(1,3) - a(1)*a(2)
          C(1) = dm(2,3) - a(1)*b(2)
          F(1,:) = BOTBC(:)   - a(1)*f(2,:)
       END IF
    END IF

    a(2:n-1) = - diff*d2(1,2:n-1) 
    b(2:n-1) = - diff*d2(2,2:n-1) + alpha
    c(2:n-1) = - diff*d2(3,2:n-1) 

!  SET UP TOP BOUNDARY CONDITION.
    IF (TOP_NEUMANN .EQ. 0) THEN
       A(N) = ZERO
       B(N) = ONE
       C(N) = ZERO
       F(N,1:nrhs) = TOPBC(1:nrhs)
    ELSE
       IF (n == ny+1) THEN
! MODIFY SECOND LAST EQUATION TO INCORPORATE BOUNDARY CONDITION IN VISCOUS FLUX
          a(n-1) =     - diff*dh(n-1)*dn(n-2)
          b(n-1) = one + diff*dh(n-1)*dn(n-2)
          c(n-1) = zero
          f(n-1,:) = f(n-1,:) + diff*dh(n-1)*topbc(:)

! LET LAST ROW OF MATRIX HOLD NEUMANN BOUNDARY CONDITION.       
! THIS ALLOWS THE APPLICATION OF A SECOND-ORDER ACCURATE NEUMANN BC.
! USE THE SECOND LAST EQUATION TO ELIMINATE THE THIRD LAST ENTRY IN LAST ROW.
          c(n) = dm(1,2)/a(n-1)
          a(n) = dm(2,2) - c(n)*b(n-1)
          b(n) = dm(3,2) - c(n)*c(n-1)
          f(n,:) = topbc(:) - c(n)*f(n-1,:)
       ELSEIF (n == ny) THEN
! USE SECOND-ORDER ACCURATE ONE-SIDED DIFFERENCE FORMULA FOR NEUMANN BC.
! ELIMINATE THIRD LAST ENTRY IN BOTTOM ROW USING SECOND LAST EQUATION.
          C(N) = dm(1,4)/a(n-1)
          A(N) = dm(2,4) - c(n)*b(n-1)
          B(N) = dm(3,4) - c(n)*c(n-1)
          F(N,:) = TOPBC(:) - c(n)*f(n-1,:)
       END IF

    END IF

!  APPLY BOUNDARY CONDITIONS ON RHS
!!$    F(1,1:nrhs) = BOTBC(1:nrhs)
!!$    F(N,1:nrhs) = TOPBC(1:nrhs)

!  SOLVE TRIDIAGONAL SYSTEM USING GAUSSIAN ELIMINATION.
!  FORWARD ELIMINATION
    DO I = 2,N
       B(I-1) = B(I-1)**(-1)
       C(I-1) = C(I-1)*B(I-1)
       F(I-1,1:nrhs) = F(I-1,1:nrhs)*B(I-1)

       B(I) = B(I) - A(I)*C(I-1)
       F(I,1:nrhs) = F(I,1:nrhs) - A(I)*F(I-1,1:nrhs)
    END DO

! SOLVE FOR THE LAST ELEMENT IN THE VECTOR
    B(N) = B(N)**(-1)
    F(N,1:nrhs) = F(N,1:nrhs)*B(N)

! BACK SUBSTITUTION.
    DO I = N-1,1,-1
       F(I,1:nrhs) = F(I,1:nrhs) - C(I)*F(I+1,1:nrhs)
    END DO

  END SUBROUTINE solve_helmholtz

!=====================================================================
!==========================SOLVE PRESSURE=============================
!=====================================================================
  SUBROUTINE SOLVE_P(P, ALPHA, TOPBC, BOTBC)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(INOUT) :: P(NY+1)
    REAL(KIND=prec), INTENT(IN)       :: ALPHA
    COMPLEX(KIND=prec), INTENT(IN)    :: TOPBC, BOTBC

    INTEGER                             j
    REAL(KIND=prec), DIMENSION(NY+1) :: A, B, C
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    A(1) =   ZERO
    B(1) = - DN(1)
    C(1) =   DN(1)

!!$    A(2) =   ZERO
!!$    B(2) = - DH(2)*DN(2) + ALPHA
!!$    C(2) =   DH(2)*DN(2)

    A(2) =   DH(2)*DN(1)
    B(2) = - DH(2)*DN(1) - DH(2)*DN(2) + ALPHA
    C(2) =                 DH(2)*DN(2)

    DO J = 3,ny-1
       A(J) = D2H(1,J)
       B(J) = D2H(2,J) + ALPHA
       C(J) = D2H(3,J)
    END DO

    A(NY) =   DH(NY)*DN(NY-1)
    B(NY) = - DH(NY)*DN(NY-1) - DH(NY)*DN(NY) + ALPHA
    C(NY) =                     DH(NY)*DN(NY)

!!$    A(NY) =   DH(NY)*DN(NY-1)
!!$    B(NY) = - DH(NY)*DN(NY-1) + ALPHA
!!$    C(NY) =   ZERO

    A(NY+1) = - DN(NY)
    B(NY+1) =   DN(NY)
    C(NY+1) =   ZERO

!  APPLY BOUNDARY CONDITIONS ON RHS

    P(1)    = BOTBC
    P(2)    = P(2)  ! + DH(2)*BOTBC
    P(NY)   = P(NY) ! - DH(NY)*TOPBC
    P(NY+1) = TOPBC

!  INVERT TRIDIAGONAL SYSTEM ON RHS.
!     ELIMINATE LOWER DIAGONAL

    DO J = 2,NY+1
       B(J-1) = B(J-1)**(-1)
       C(J-1) = C(J-1)*B(J-1)
       P(J-1) = P(J-1)*B(J-1)

       B(J) = B(J) - A(J)*C(J-1)
       P(J) = P(J) - A(J)*P(J-1)
    END DO

    P(NY+1) = P(NY+1)*B(NY+1)**(-1)

    DO J = NY,1,-1
       P(J) = P(J) - C(J)*P(J+1)
    END DO

    P(1) = (BOTBC - DM(2,1)*P(2) - DM(3,1)*P(3))/DM(1,1)
    P(NY+1) = (TOPBC - DM(1,2)*P(NY-1) - DM(2,2)*P(NY))/DM(3,2)

  END SUBROUTINE solve_p

END MODULE solve 
  
