PROGRAM test_jet_forcing
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE transform
  USE pass2
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  INTEGER, PARAMETER :: prec=SELECTED_REAL_KIND(14,32)
  
  INTEGER              j, x, z, ierr, mpierr
  REAL(KIND=prec)      t1, t0
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:) :: u, v
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:) :: tmp

  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:,:) :: junk

! Beginning of main program
  CALL initialize_mpi
  CALL read_parms
  CALL distribute_data
  CALL initialize_fft

  CALL define_coordinates
  CALL define_parameters

  t_total = 0.d0
  dt = 0.001d0
  CALL jet_profile(t_total,dt,cbc_top(1,1,2),cbc_bot(1,1,2))
  CALL write_jet

  DO j = 1,30
     t_total = t_total + dt
     WRITE(*,*) j
     t0 = MPI_WTIME()
     CALL jet_profile(t_total,dt,cbc_top(1,1,2),cbc_bot(1,1,2))
     t1 = MPI_WTIME()
     IF (me==0) write(*,999) t1-t0
999  format('WTIME = ',f12.6,'  OVER ONE CALL')
     CALL write_jet
  END DO
  STOP 'Normal Termination'

CONTAINS
  SUBROUTINE write_jet
    IMPLICIT NONE

    INCLUDE 'hdf.f90'
    
    INTEGER status, dssdims, dssnt, dssdisc, dsadata
    INTEGER jet_rank, jet_length(2)
    REAL, DIMENSION(nx+1)      :: jet_xcoord
    REAL, DIMENSION(nz+1)      :: jet_zcoord
    REAL, DIMENSION(nx+1,nz+1) :: data
    REAL(KIND=prec), DIMENSION(nx,nz) :: u2
    REAL(KIND=prec), DIMENSION(nx,local_nz_small) :: u1
    
    jet_rank = 2
    jet_length(1) = nx+1
    jet_length(2) = nz+1
    jet_xcoord = xcoord
    jet_zcoord = zcoord
    CALL small_fft_complex_to_real(cbc_bot(1,1,2),u1)
    CALL MPI_GATHERV(u1,nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,   &
         u2,nx*nz_small_proc(0:nproc-1),nx*z_small_start(0:nproc-1), &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
    IF (me == 0) THEN
       data(1:nx,1:nz)   = u2
       data(1:nx,nz+1)   = data(1:nx,1)
       data(nx+1,1:nz+1) = data(1,1:nz+1)

       status = dssdims(jet_rank,jet_length)   ! CREATE AN ARRAY
       status = dssnt(DFNT_FLOAT32)            ! SPECIFY DATA TYPE
       status = dssdisc(1,jet_length(1),jet_xcoord) ! DIMENSIONS
       status = dssdisc(2,jet_length(2),jet_zcoord)
       status = dsadata('jet.hdf',jet_rank,jet_length,data) ! WRITE ARRAY
    END IF

  END SUBROUTINE write_jet
END PROGRAM test_jet_forcing
