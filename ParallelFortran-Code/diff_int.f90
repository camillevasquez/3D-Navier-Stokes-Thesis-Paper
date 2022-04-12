MODULE diff_int
  USE runparms
  USE gridparms
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: d_dy, d2_dy, int_y

CONTAINS
!=====================================================================
!===================================D DY==============================
!=====================================================================
  SUBROUTINE d_dy(n1,f,n2,df,d)
    IMPLICIT NONE
    INTEGER             n1, n2
    REAL(KIND=prec)    d(1:n2)
    COMPLEX(KIND=prec) f(1:n1), df(1:n2)

    IF (n1 == ny) THEN
      df(1)      = d(1)     *(f(2)    - f(1))
      df(2:n2-1) = d(2:n2-1)*(f(2:n1) - f(1:n1-1))
      df(n2)     = d(n2)    *(f(n1)   - f(n1-1))
    ELSEIF (n1 == ny+1) THEN
      df(1:n2) = d(1:n2)*(f(2:n1) - f(1:n1-1))
    END IF     

  END SUBROUTINE d_dy

!=====================================================================
!==================================D2 DY==============================
!=====================================================================
  SUBROUTINE d2_dy(n,f,df,d)
    IMPLICIT NONE
    INTEGER             n
    COMPLEX(KIND=prec) f(1:n), df(1:n)
    REAL(KIND=prec)    d(1:3,1:n)

    df(1)     = d(1,1)    *f(1)     + d(2,1)    *f(2)     + d(3,1)    *f(3)
    df(2:n-1) = d(1,2:n-1)*f(1:n-2) + d(2,2:n-1)*f(2:n-1) + d(3,2:n-1)*f(3:n)
    df(n)     = d(1,n)    *f(n-2)   + d(2,n)    *f(n-1)   + d(3,n)    *f(n)

!!$    IF (n == ny+1) THEN
!!$       df(1) = SUM(d2m(1:4,1)*f(1:4))
!!$       df(ny+1) = SUM(d2m(1:4,2)*f(ny-2:ny+1))
!!$    END IF

  END SUBROUTINE d2_dy
!=====================================================================
!===================================D DY==============================
!=====================================================================
  SUBROUTINE int_y(n,f,dy,fint)
    IMPLICIT NONE
    INTEGER             n
    REAL(KIND=prec)    f(1:n), dy(1:n), fint

    fint = sum(f(1:n)*dy(1:n))

  END SUBROUTINE int_y

END MODULE diff_int
