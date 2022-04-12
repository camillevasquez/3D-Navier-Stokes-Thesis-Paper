PROGRAM test
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE init_vel
  USE pass1
  USE pass2
  USE transform
  USE read_write
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(14,32)
  
  INTEGER               i, j, x, z, ierr, mpierr
  REAL(KIND=prec)       t1, t0

  INTEGER           NX_TEST, NY_TEST, NZ_TEST, status(MPI_STATUS_SIZE)
  REAL(KIND=prec)  LX_TEST, LZ_TEST, RE_TEST, DT_TEST
  LOGICAL           ABORT
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:) :: slice, slice2

! Beginning of main program
  CALL init
  CALL define_parameters
  CALL define_coordinates

  ALLOCATE(slice(1:nz,1:ny+1,1:3+nscalar),slice2(1:nz,1:ny+1,1:3+nscalar), &
       STAT=ierr)

  OPEN (UNIT=14,FILE='newvel.1',STATUS='OLD',FORM='UNFORMATTED') 
  READ(14) NX_TEST,NY_TEST,NZ_TEST,LX_TEST,LZ_TEST,RE_TEST,DT_TEST,T_TOTAL
  
  OPEN (UNIT=15,FILE='newvel.5',STATUS='OLD',FORM='UNFORMATTED') 
  READ(15) NX_TEST,NY_TEST,NZ_TEST,LX_TEST,LZ_TEST,RE_TEST,DT_TEST,T_TOTAL
  
  DO x = 1,nx/2
     read(14) slice
     read(15) slice2

     write(*,*) x, MAXVAL(ABS(slice-slice2))
  END DO

  CALL finalize
! End of main program
END PROGRAM test
