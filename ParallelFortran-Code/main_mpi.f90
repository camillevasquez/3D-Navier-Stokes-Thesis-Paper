PROGRAM mpiflow
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE diff_int
  USE init_vel
  USE init_jet_buffer
  USE pass1
!!$  USE pass1_dynamic
  USE pass1_linearized
!!$  USE dynamic2
!!$  USE pass1_dynamic2
  USE pass2
  USE transform
  USE read_write
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  INTEGER              j, x, z, ierr
  REAL(KIND=prec)      t1, t0
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:) :: u, v
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:) :: tmp

  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:,:) :: junk

! Beginning of main program
  CALL initialize_mpi
  CALL read_parms
  CALL distribute_data
  CALL allocate_data_arrays
  CALL initialize_fft
  IF (me .eq. 0) WRITE(*,*) 'INITIALIZED FFTS'

  CALL define_coordinates
  IF (me .eq. 0) WRITE(*,*) 'DEFINED COORDINATES'
  CALL define_parameters
  IF (me .eq. 0) WRITE(*,*) 'DEFINED PARAMETERS'
  CALL define_initial_field

  IF (vel_out == 1) CALL output_hdf_dfsd

  CALL timestep

  IF (save_vel==1) CALL save_velocity
  CALL deallocate_data_arrays
  CALL finalize_fft
  CALL finalize_mpi
  write(*,*) 'Normal Termination'
  STOP
! End of main program
CONTAINS
!====================================================
  SUBROUTINE timestep
    IMPLICIT NONE
    REAL(KIND=prec) ::temp, temp2
    INTEGER mpierr, maxm

    IF (me==0) THEN
       WRITE(*,*) 'INITIAL DT = ', dt
       WRITE(*,*) 'INITIAL TIME = ', t_total
    END IF
    
    maxm = MAX(nx,local_nz_small,ny-1)

    call l2_norm

    t0 = MPI_WTIME()
    DO step = 1,tend

!!en       CALL advance_scalar(maxm)
       call boundary
!!$       IF (subgrid == 1) CALL compute_dynamic_coef(filter_dim,lagr_avg)
       DO rk_step = 1,rk_last

          IF (subgrid == 1) THEN
!!$             CALL nonlinear_dynamic2
          ELSEIF (linearized_run == 1) THEN
             CALL nonlinear_linearized
          ELSE
             CALL nonlinear
          END IF
          CALL rhs_solve

       END DO
       T_TOTAL = T_TOTAL + DT

       IF (mod(step,output_step) == 0) THEN
          t1 = MPI_WTIME()
          IF (me==0) write(*,999) t1-t0, output_step
          999 format('WTIME = ',f12.6,'  OVER ', i6,'  STEPS')
          t0 = MPI_WTIME()
       END IF

       IF ((mod(step,output_step)==0).and.(vel_out==1))  CALL output_hdf_dfsd
       IF ((mod(step,output_step)==0).and.(snap_out/=0)) CALL output_snapshot
       IF ((mod(step,save_step)==0) .and. (save_vel==1)) CALL save_velocity
       IF (mod(step,l2_step) == 0)     call l2_norm


       IF (mod(step,output_step) == 0) THEN
          t1 = MPI_WTIME()
          IF (me==0) write(*,998) t1-t0
          998 format('OUTPUT TIME = ',f12.6)
          t0 = MPI_WTIME()
       END IF
       
       IF ((crossflow > 0) .and. (nscalar > 0)) THEN
          temp = MAXVAL(data_scalar)
          CALL MPI_REDUCE(temp,temp2,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
               0,MPI_COMM_WORLD,mpierr)
          IF (me == 0) WRITE(*,*) 'MAX SCALAR VALUE = ', temp2
       END IF

    END DO
    990 format(6e12.4)
    995 format(12e12.4)


    IF (which_test == 6) THEN
       ALLOCATE(u(1:nx2,1:nz2),v(1:nx2,1:nz2),tmp(1:ny+1),STAT=ierr)
       IF (ierr==0) WRITE(*,*) 'BEGIN ERROR CHECK'
       do j = 1,ny+1
          call xz_fft_complex_to_real(data1(1:nz,1:local_nx,j,1),u)
          do z = 1,nz2
             do x = 1,nx2
                v(x,z) = u(x,z) + cos(2.d0*pi*dble(x-1)/dble(nx2)) &
                     *sin(pi*yh(j))*exp(-2.d0*pi**2*t_total/re)
             end do
          end do
          tmp(j) = MAXVAL(ABS(v(1:nx2,1:nz2)))
       end do
       write(*,*) 'MAX ERROR IN U = ', MAXVAL(tmp)

       do j = 1,ny
          call xz_fft_complex_to_real(data1(1:nz,1:local_nx,j,2),u)
          do z = 1,nz2
             do x = 1,nx2
                v(x,z) = u(x,z) - sin(2.d0*pi*dble(x-1)/dble(nx2)) &
                     *cos(pi*y(j))*exp(-2.d0*pi**2*t_total/re)
             end do
          end do
          v = ABS(v)
          tmp(j) = MAXVAL(v)
       end do
       write(*,*) 'MAX ERROR IN V = ', MAXVAL(tmp(1:ny))
       
       do j = 1,ny+1
          call xz_fft_complex_to_real(data4(j,1:nz,1:local_nx),u)
          do z = 1,nz2
             do x = 1,nx2
                v(x,z) = u(x,z) +0.25d0*(cos(4.d0*pi*dble(x-1)/dble(nx2)) &
                     + cos(2.d0*pi*yh(j)))*exp(-4.d0*pi**2*t_total/re)
             end do
          end do
          tmp(j) = MAXVAL(ABS(v(1:nx2,1:nz2)))
       end do
       write(*,*) 'MAX ERROR IN P = ', MAXVAL(tmp)
       DEALLOCATE(u,v,tmp)
    END IF

    IF (me==0) WRITE(*,*) 'FINAL TIME = ', t_total


END SUBROUTINE timestep


END PROGRAM mpiflow

