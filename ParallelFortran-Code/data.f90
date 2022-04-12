MODULE data
  USE runparms
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: data1, data2, data4, data_dyn, dyn_visc, &
       data_scalar, data_rhs, baseflow

! data arrays
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:,:) :: data1, data2
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:)   :: data4
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:,:)    :: data_scalar, data_rhs

! array to hold dynamic model coefficients
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:,:) :: data_dyn
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:)     :: dyn_visc

! array to hold baseflow for linearized simulations
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:,:) :: baseflow

END MODULE data
