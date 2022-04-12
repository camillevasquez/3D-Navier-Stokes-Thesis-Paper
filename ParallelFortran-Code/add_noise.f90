PROGRAM add_noise
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE pass_data
  USE read_write
  USE transform
  USE solve
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(14,32)
  
  INTEGER           i, j, k, x, z, plane, ierr
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:) :: cu
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:)   :: cdummy
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:)      :: ushape, vshape
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:)    :: temp, temp2
  REAL(KIND=prec) :: tmp

  REAL(KIND=prec), PARAMETER :: zero=0.0d0, half=0.5d0, one=1.0d0, two=2.0d0

! Beginning of main program
  CALL initialize_mpi
  CALL read_parms
  CALL distribute_data
  CALL allocate_data_arrays

  CALL define_coordinates
  CALL define_parameters

  CALL read_velocity

  ALLOCATE(cu(ny+1,3),cdummy(ny+1),ushape(ny+1),vshape(ny),temp(ny+1,6), &
       temp2(ny+1,6), STAT=ierr)
  DO i = 1,3
     DO z = 1,nz
        DO x = 1,local_nx
           IF (k2_me(z,x) .ne. zero) THEN
              ! ADD NOISE IN NON-ZERO WAVENUMBERS.
              cu = zero
              CALL RANDOM_NUMBER(TEMP)
              temp2 = 1.d-4*two*(temp-half) &
                   *half*(one + cos(pi*kx_me(x)/MAXVAL(kx))) &
                   *half*(one + cos(pi*kz(z)   /MAXVAL(kz)))
              ! FILTER NOISE IN WALL-NORMAL DIRECTION.
              temp(1,1:6) = temp2(1,1:6)
              DO k = 2,ny
                 temp(k,:) = 0.25d0*(temp2(k-1,:)+two*temp2(k,:)+temp2(k+1,:))
              END DO
              temp(ny+1,:) = temp2(ny+1,:)
              ! SHAPE NOISE TO SUPPRESS NOISE NEAR BOUNDARIES
              ushape = (tanh(yh-yh(1))*tanh(yh(ny+1)-yh))**2
              vshape = (tanh(y  -y(1))*tanh(y(ny)   -y))**2
              cu(1:ny+1,1) = ushape*cmplx(temp(1:ny+1,1),temp(1:ny+1,2))
              cu(1:ny,2)   = vshape*cmplx(temp(1:ny,3),  temp(1:ny,4))
              cu(1:ny+1,3) = ushape*cmplx(temp(1:ny+1,5),temp(1:ny+1,6))
              ! REMOVE DIVERGENCE FROM NOISE.
              CALL strip_divergence(cu,cdummy,one,x,z)
              ! ADD TO VELOCITY FIELD.
              data1(z,x,1:ny+1,1:3) = data1(z,x,1:ny+1,1:3) + cu
           END IF
        END DO
     END DO
  END DO
  DEALLOCATE(cu,cdummy,ushape,vshape,temp,STAT=ierr)

  CALL save_velocity

!$$ SHUT DOWN MPI BEFORE STOPPING.
  CALL deallocate_data_arrays
  CALL finalize_mpi
  STOP 'Normal Termination'
! End of main program
END PROGRAM add_noise
