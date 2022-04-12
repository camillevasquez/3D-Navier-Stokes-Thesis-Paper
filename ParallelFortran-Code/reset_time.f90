PROGRAM reset_time
  USE runparms
  USE gridparms
  USE setup
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  INTEGER               i, j, jj, k, x, z, ierr, mpierr, nx2_test, nz2_test
  REAL(KIND=prec)       t1, t0

  INTEGER           NX_TEST, NY_TEST, NZ_TEST, status(MPI_STATUS_SIZE)
  REAL(KIND=prec)  LX_TEST, LZ_TEST, RE_TEST, DT_TEST
  LOGICAL           ABORT
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:) :: slice, slice2, temp
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:)     :: temp_in, temp_out
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:)        :: filterx, filterz
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:)   :: &
       tempx1, tempx2, tempy1, tempy2, tempz1, tempz2
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:) :: &
       scalar_slice, scalar_slice2, temp_xz_scalar
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:) :: &
       scalar, scalar_xz_interp

  REAL(KIND=prec), PARAMETER :: zero=0.0, half=0.5, one=1.0, two=2.0

! Beginning of main program
  CALL initialize_mpi
  CALL read_parms
  CALL define_coordinates
  CALL define_parameters

  ALLOCATE(filterx(nx/2+1),filterz(nz),STAT=ierr)
  DO z = 1,nz/2+1
     filterz(z) = cos(half*pi*DBLE(z-1)/DBLE(nz/2))
!!$     WRITE(*,*) filterz(z)
  END DO
  DO z = nz/2+2,nz
     filterz(z) = cos(half*pi*DBLE(nz+1-z)/DBLE(nz/2))
!!$     WRITE(*,*) filterz(z)
  END DO
  WRITE(*,*)
  DO x = 1,nx/2+1
     filterx(x) = cos(half*pi*DBLE(x-1)/DBLE(nx/2))
!!$     WRITE(*,*) filterx(x)
  END DO
  WRITE(*,*)

  OPEN (UNIT=14,FILE='savevel.old',FORM='UNFORMATTED') 
  READ(14) NX_TEST,NY_TEST,NZ_TEST,LX_TEST,LZ_TEST,RE_TEST,DT_TEST,T_TOTAL

  ALLOCATE(slice(1:nz_test,1:ny_test+1,1:3), &
       temp(1:nz,1:ny_test+1,1:3), &
       slice2(1:nz,1:ny+1,1:3), &
       temp_in(ny_test+1),temp_out(ny+1),STAT=ierr)

  WRITE(*,*) NX_TEST,NY_TEST,NZ_TEST 
  WRITE(*,*) LX_TEST,LZ_TEST,RE_TEST
  WRITE(*,*) DT_TEST,T_TOTAL
  t_total = 0.d0
  
  OPEN (UNIT=15,FILE='savevel.new',FORM='UNFORMATTED') 
  WRITE(15) NX,NY,NZ,LX,LZ,RE,DT_TEST,T_TOTAL
  
!$$ READ IN, INTERPOLATE AND WRITE OUT VELOCITY FIELD.
  DO x = 1,nx/2+1
     WRITE (*,*) 'X = ', x, '  OUT OF ', nx/2
! READ OFF y-z SLICE FROM savevel.old
     IF (x <= nx_test/2+1) THEN
        read(14) slice
     ELSE
        slice = zero
     END IF
! SET UP ARRAY WITH CORRECT NUMBER OF z WAVENUMBERS
     temp = zero
     IF (nz_test >= nz) THEN
        DO z = 1,nz/2+1
           temp(z,:,:)  = slice(z,:,:)*filterx(x)*filterz(z)
        END DO
        DO z = nz/2+2,nz
           temp(z,:,:) = slice(z+nz_test-nz,:,:)*filterx(x)*filterz(z)
        END DO
     ELSE
        DO z = 1,nz_test/2+1
           temp(z,:,:) = slice(z,:,:)*filterx(x)*filterz(z)
        END DO
        DO z = nz-nz_test/2+2,nz
           temp(z,:,:) = slice(z+nz_test-nz,:,:)*filterx(x)*filterz(z)
        END DO
     END IF
     DO i = 1,3
        DO z = 1,nz
           temp_in = temp(z,1:ny_test+1,i)
           SELECT CASE (i)
           CASE (2)
              CALL interpolate(ny_test,ny,temp_in,temp_out)
              temp_out(ny+1) = zero
           CASE DEFAULT
              CALL interpolate(ny_test+1,ny+1,temp_in,temp_out)
           END SELECT
           slice2(z,1:ny+1,i) = temp_out
        END DO
     END DO
     WRITE(15) slice2
  END DO
  DEALLOCATE(slice,temp,slice2,temp_in,temp_out,STAT=ierr)

  WRITE(*,*) 'BEGIN SCALAR CONVERSION'

!$$ READ IN, INTERPOLATE AND WRITE OUT SCALAR FIELD.
  ALLOCATE(scalar_slice(nx_test,nz_test), scalar_slice2(nx,nz), STAT=ierr)

  DO i = 1,nscalar
!$$ READ OLD SCALAR FIELD PLANE BY PLANE, INTERPOLATE PLANES ONTO
!$$ NEW GRID AND STORE.
     DO k = 1,ny_test+1
        READ(14) scalar_slice
        scalar_slice2 = scalar_slice
        WRITE(15) scalar_slice2
     END DO
  END DO

!$$ CLOSE INPUT AND OUTPUT FILES.
  CLOSE (UNIT=14)
  CLOSE (UNIT=15)

!$$ SHUT DOWN MPI BEFORE STOPPING.
  CALL finalize_mpi
  STOP 'Normal Termination'
! End of main program
CONTAINS
!================================================
  SUBROUTINE interpolate(n1,n2,in,out)
    IMPLICIT NONE
    INTEGER n1,n2, i, nlast
    COMPLEX(KIND=prec), DIMENSION(n1) :: in, temp
    COMPLEX(KIND=prec), DIMENSION(n2) :: out
    REAL(KIND=prec), DIMENSION(n1) :: x1
    REAL(KIND=prec), DIMENSION(n2) :: x2
    
    DO i = 1,n1
       x1(i) = DBLE(i-1)/DBLE(n1-1)
    END DO

    DO i = 1,n2
       x2(i) = DBLE(i-1)/DBLE(n2-1)
    END DO

    IF (n1 >= 3*n2/2) THEN
! FILTER INPUT PROFILE
       DO i = 2,n1-1
          temp(i) = 0.25d0*in(i-1) + 0.5d0*in(i) + 0.25d0*in(i+1)
       END DO
       in(2:n1-1) = temp(2:n1-1)
    END IF

    IF (n1 >= 3*n2) THEN
! FILTER INPUT PROFILE
       DO i = 2,n1-1
          temp(i) = 0.25d0*in(i-1) + 0.5d0*in(i) + 0.25d0*in(i+1)
       END DO
       in(2:n1-1) = temp(2:n1-1)
    END IF

    out(1) = in(1)
    nlast = 1
    DO i = 2,n1
       DO WHILE (x2(nlast+1) < x1(i)) 
          nlast = nlast + 1
          out(nlast) = in(i-1) + (x2(nlast) - x1(i-1))*(in(i)-in(i-1)) &
               /(x1(i) - x1(i-1))
       END DO
    END DO
    out(n2) = in(n1)

  END SUBROUTINE interpolate
!================================================
  SUBROUTINE interpolate_real(n1,n2,in,out)
    IMPLICIT NONE
    INTEGER n1,n2, i, nlast
    REAL(KIND=prec), DIMENSION(n1) :: in, temp
    REAL(KIND=prec), DIMENSION(n2) :: out
    REAL(KIND=prec), DIMENSION(n1) :: x1
    REAL(KIND=prec), DIMENSION(n2) :: x2
    
    DO i = 1,n1
       x1(i) = DBLE(i-1)/DBLE(n1-1)
    END DO

    DO i = 1,n2
       x2(i) = DBLE(i-1)/DBLE(n2-1)
    END DO

    IF (n1 >= 3*n2/2) THEN
! FILTER INPUT PROFILE
       DO i = 2,n1-1
          temp(i) = 0.25d0*in(i-1) + 0.5d0*in(i) + 0.25d0*in(i+1)
       END DO
       in(2:n1-1) = temp(2:n1-1)
    END IF

    IF (n1 >= 3*n2) THEN
! FILTER INPUT PROFILE
       DO i = 2,n1-1
          temp(i) = 0.25d0*in(i-1) + 0.5d0*in(i) + 0.25d0*in(i+1)
       END DO
       in(2:n1-1) = temp(2:n1-1)
    END IF

    out(1) = in(1)
    nlast = 1
    DO i = 2,n1
       DO WHILE (x2(nlast+1) < x1(i)) 
          nlast = nlast + 1
          out(nlast) = in(i-1) + (x2(nlast) - x1(i-1))*(in(i)-in(i-1)) &
               /(x1(i) - x1(i-1))
       END DO
    END DO
    out(n2) = in(n1)

  END SUBROUTINE interpolate_real
END PROGRAM reset_time
