PROGRAM convert_float
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE transform
  USE init_jet_buffer
  USE pass1_linearized
  USE pass2
  USE solve
  USE diff_int
  IMPLICIT NONE

  INCLUDE 'hdf.f90'
  INCLUDE 'mpif.h'
  
  INTEGER, PARAMETER :: prec = SELECTED_REAL_KIND(14,32)
  REAL(KIND=prec), DIMENSION(6) :: time
  REAL(KIND=prec) :: cr, ci
  INTEGER mpierr, ierr, x, z, jj, plane, i, j
  REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! Beginning of main program
  CALL initialize_mpi
  CALL read_parms
  CALL distribute_data
  CALL allocate_data_arrays
!  test_run = 1 !ESTIMATE OPTIMAL FFT, RATHER THAN SPENDING TIME MEASURING.
  CALL initialize_fft

  CALL define_coordinates
  CALL define_parameters

  linearized_run = 1
  CALL init_baseflow          ! READ IN BASEFLOW AND SET TIME STEP.
  dt = half*dt

  IF (xfringe+zfringe > 0) CALL init_buffer      ! INITIALIZE BUFFER REGION
  IF (linearized_run == 1) target_inflow = zero  ! PERTURBATIONS -> 0 IN BUFFER

  CALL read_evector
  DO i = 1,3
     DO j = 1,ny+1
        data2(j,i,1:nz,1:local_nx) = data1(1:nz,1:local_nx,j,i)
     END DO
  END DO
  CALL nonlinear_linearized
  CALL solve_continuous_poisson
  CALL output_hdf_dfsd

  CALL deallocate_data_arrays
  CALL finalize_fft
  CALL finalize_mpi
  STOP 'Normal Termination'
! End of main program
CONTAINS
!=====================================================================
  SUBROUTINE solve_continuous_poisson
    IMPLICIT NONE

    COMPLEX(KIND=prec), DIMENSION(ny+1,3) :: cu, cru
    COMPLEX(KIND=prec), DIMENSION(ny+1)   :: du
    REAL(KIND=prec)    visc
    INTEGER            x, z

    DO X = 1,Klim
       DO Z = 1,Mlim

          !  READ IN THE VELOCITY FIELD AND NONLINEAR TERM.
          cu(1:ny+1,1:3)  = data2(1:ny+1,1:3,z,x)
          cru(1:ny+1,1:3) = - data1(z,x,1:ny+1,1:3)

          !  ADD THE VISCOUS TERM TO THE NONLINEAR TERM.         
          VISC = ONE/RE
          CALL d2_dy(ny+1,cu(1,1),du,d2h)
          cru(1:ny+1,1) = cru(1:ny+1,1) &
               + visc*(du(1:ny+1) - k2_me(z,x)*cu(1:ny+1,1))
          CALL d2_dy(ny,  cu(1,2),du,d2n)
          cru(1:ny,2) = cru(1:ny,2)  &
               + visc*(du(1:ny)   - k2_me(z,x)*cu(1:ny,2))
          CALL d2_dy(ny+1,cu(1,3),du,d2h)
          cru(1:ny+1,3) = cru(1:ny+1,3) &
               + visc*(du(1:ny+1) - k2_me(z,x)*cu(1:ny+1,3))

          IF (K2_ME(z,x) .EQ. ZERO) THEN
             !  SOLVE FOR PRESSURE IN THE MEAN MODE.
             data4(1:ny+1,z,x) = ZERO
          ELSE
             !  FOR EACH WAVENUMBER PAIR (X,Z), 
             !  TAKE THE DIVERGENCE OF THE RHS, 
             !  COMPUTE THE BOUNDARY CONDITIONS AND 
             !  SOLVE THE POISSON EQUATION FOR THE PRESSURE
             CALL strip_divergence(cru,data4(1:ny+1,z,x),one,x,z)
          END IF
       END DO
    END DO
  END SUBROUTINE solve_continuous_poisson
!=========================================================================
  SUBROUTINE read_evector
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(nx2,nz2,ny+1) :: data
    REAL(KIND=prec), DIMENSION(nx2,local_nz) :: u1
    INTEGER, PARAMETER :: time_rank = 1, time_length = 6, vel_rank = 3
    INTEGER :: plane, i, status, dsrref, dsgdata, vel_length(vel_rank)

    REAL(KIND=prec), DIMENSION(nx2,local_nz,ny+1,3) :: evector
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
    INTEGER :: n, x, z, jj
    CHARACTER(LEN=13) :: filename

    vel_length(1) = nx2
    vel_length(2) = nz2
    vel_length(3) = ny+1

    IF (me == 0) THEN
       ! GENERATE FILENAME ACCORDING TO EIGENVECTOR NUMBER.
       filename = 'in.hdf'
       WRITE(*,*) 'READING ', filename

       ! READ IN THE TIME
       ! TIME ARRAY HAS SIX ELEMENTS: time, real/imag part of growth rate,
       ! mean pressure gradient, reynolds number and prandtl number.
       status = dsrref(filename,2)
       status = dsgdata(filename,time_rank,time_length,time)
       cr = time(2)
       ci = time(3)
    END IF

    DO i = 1,3+nscalar
       ! READ IN THE VELOCITY/SCALAR ARRAY
       IF (me == 0) THEN
          status = dsrref(filename,2+i)
          status = dsgdata(filename,vel_rank,vel_length,data)
          WRITE(*,*) 'READ FIELD ', i, ' WITH STATUS ', status
       END IF
       ! TRANSFORM INTO FOURIER SPACE AND STORE IN DATA1
       DO plane = 1,ny+1
          CALL MPI_SCATTERV( &
               data(1,1,plane),nx2*nz_proc,nx2*z_start,MPI_DOUBLE_PRECISION, &
               u1,nx2*local_nz,MPI_DOUBLE_PRECISION, &
               0,MPI_COMM_WORLD,mpierr)
          CALL xz_fft_real_to_complex(u1,data1(1,1,plane,i))
      END DO
    END DO
! ZERO OUT THE EXTRA PLANE IN THE VERTICAL VELOCITY.
    data1(1:nz,1:local_nx,ny+1,2) = zero

  END SUBROUTINE read_evector
!=================================================================
!=========================OUTPUT HDF DFSD=========================
!=================================================================
  SUBROUTINE output_hdf_dfsd
    IMPLICIT NONE
    
    INCLUDE 'hdf.f90'

    COMPLEX(KIND=prec), DIMENSION(nz,local_nx)    :: cu
    REAL(KIND=prec), DIMENSION(nx,local_nz_small) :: u1
    REAL(KIND=prec), DIMENSION(nx,nz)             :: u2
    REAL(KIND=prec) :: temp(4+nscalar), poss_max, temp2(ny,1+ndynscalar), &
         temp3(ny,1+ndynscalar)
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    REAL, DIMENSION(:,:,:), ALLOCATABLE :: data
    REAL, DIMENSION(:), ALLOCATABLE :: vel_xcoord, vel_ycoord, vel_zcoord, &
         tmp_dyncoef
    REAL :: TIME_OUT(4)
         
    INTEGER, DIMENSION(:), ALLOCATABLE  :: xindex, yindex, zindex
    INTEGER xlength, ylength,zlength

! DECLARE HDF FUNCTIONS
    INTEGER  dsadata, dssdims, dssdisc, dspdata, dssdast, dsslens
    INTEGER  dssnt, dsnum, dswref
! DECLARE HDF VARIABLES
    INTEGER, PARAMETER :: vel_rank = 3
    INTEGER, DIMENSION(0:vel_rank-1)::vel_length
    INTEGER  time_rank, time_length
    INTEGER  num_datasets, fname_rank, fname_length, status
    INTEGER  subgrid_flags(2)
    CHARACTER(LEN=4) time_name
    CHARACTER(LEN=1) vel_name(1:5), dim_name(0:2)
    CHARACTER(LEN=12) FNAME

! SET UP LENGTH OF OUTPUT ARRAYS AND INDEXING INTO FULL DATA ARRAYS.
    xlength = nx/xskip + 1
    ylength = ny/yskip + 1
    zlength = nz/zskip + 1
    vel_length(0) = xlength
    vel_length(1) = ylength
    vel_length(2) = zlength
    ALLOCATE(vel_xcoord(xlength), xindex(xlength),  &
         vel_ycoord(ylength), yindex(ylength),  &
         vel_zcoord(zlength), zindex(zlength), STAT=ierr)
    CALL generate_coord_index(nx,xskip,xlength,xcoord,vel_xcoord,xindex)
    CALL generate_coord_index(ny,yskip,ylength,yh,    vel_ycoord,yindex)
    CALL generate_coord_index(nz,zskip,zlength,zcoord,vel_zcoord,zindex)

    IF (me==0) THEN
       ALLOCATE(data(xlength,ylength,zlength), STAT=ierr)
! DESCRIBE VELOCITY FIELD
       vel_name(1) = 'u'
       vel_name(2) = 'v'
       vel_name(3) = 'w'
       vel_name(4) = 't'
       vel_name(5) = 's'
! DESCRIBE COORDINATES
       dim_name(0) = 'x'
       dim_name(1) = 'y'
       dim_name(2) = 'z'
! DESCRIBE VELOCITY FIELD
       time_rank = 1
       time_length = 4
       time_name = 'time'
       time_out(1) = t_total
       time_out(2) = REAL(cr)
       time_out(3) = REAL(ci)
       time_out(4) = p_grad

       FNAME = 'out.hdf'
! WRITE OUT TIME STAMP AS FIRST ENTRY IN FILE.
       status = dswref(fname,2)
       status = dssdims(time_rank,time_length)
       status = dssnt(DFNT_FLOAT32)
       status = dsslens(4,4,5,4)
       status = dssdast('time','hL/U','f10.5','none')
       status = dspdata(fname,time_rank,time_length,time_out)

    END IF

    DO i = 1,3+nscalar
       poss_max = zero
       DO jj = 1,ylength

          plane = yindex(jj)
          IF (i == 2) THEN
             IF (plane == 1) THEN
                cu = data1(1:nz,1:local_nx,1,i)
             ELSEIF (plane == ny+1) THEN
                cu = data1(1:nz,1:local_nx,ny,i)
             ELSE
                cu = half*(data1(1:nz,1:local_nx,plane-1,i) &
                     + data1(1:nz,1:local_nx,plane,i))
             END IF
          ELSE
             cu = data1(1:nz,1:local_nx,plane,i)
          END IF

          CALL small_fft_complex_to_real(cu,u1)
          CALL MPI_GATHERV(u1,nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,   &
               u2,nx*nz_small_proc(0:nproc-1),nx*z_small_start(0:nproc-1), &
               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
          poss_max = MAX(MAXVAL(ABS(u1)),poss_max)
          IF (me == 0) THEN
             DO z = 1,zlength-1
                DO x = 1,xlength-1
                   data(x,jj,z) = u2(xindex(x),zindex(z))
                END DO
             END DO
             data(xlength,jj,1:zlength-1) = data(1,jj,1:zlength-1)
             data(1:xlength,jj,zlength)   = data(1:xlength,jj,1)
          END IF
       END DO
       CALL MPI_Reduce(poss_max,temp(i),1,MPI_DOUBLE_PRECISION, &
            MPI_MAX,0,MPI_COMM_WORLD,mpierr)

! WRITE THE CURRENT COMPONENT OF VELOCITY/SCALAR INTO THE HDF FILE
       IF (me == 0) THEN
          status = dswref(fname,2+i)
          status = dssdims(vel_rank,vel_length)   ! CREATE AN ARRAY
          status = dssnt(DFNT_FLOAT32)            ! SPECIFY DATA TYPE
          status = dsslens(1,4,5,9)               ! SPECIFY LABELS
          status = dssdast(vel_name(i),'U','f10.5','cartesian')
          status = dssdisc(1,vel_length(0),vel_xcoord) ! DIMENSIONS
          status = dssdisc(2,vel_length(1),vel_ycoord)
          status = dssdisc(3,vel_length(2),vel_zcoord)
          status = dsadata(fname,vel_rank,vel_length,data) ! WRITE ARRAY
       END IF

    END DO

! SAVE THE PRESSURE TO THE DATA FILE AS WELL.
    DO jj = 1,ylength
       poss_max = zero
       plane = yindex(jj)
       cu = data4(plane,1:nz,1:local_nx)
       CALL small_fft_complex_to_real(cu,u1)
       CALL MPI_GATHERV(u1,nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,   &
            u2,nx*nz_small_proc(0:nproc-1),nx*z_small_start(0:nproc-1), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
          poss_max = MAX(MAXVAL(ABS(u1)),poss_max)
       IF (me == 0) THEN
          DO z = 1,zlength-1
             DO x = 1,xlength-1
                data(x,jj,z) = u2(xindex(x),zindex(z))
             END DO
          END DO
          data(xlength,jj,1:zlength-1) = data(1,jj,1:zlength-1)
          data(1:xlength,jj,zlength)   = data(1:xlength,jj,1)
       END IF
    END DO
    CALL MPI_Reduce(poss_max,temp(4+nscalar),1,MPI_DOUBLE_PRECISION, &
            MPI_MAX,0,MPI_COMM_WORLD,mpierr)

! WRITE THE PRESSURE TO THE DATA FILE.
    IF (me == 0) THEN
! CREATE AN ARRAY FOR THE CURRENT COMPONENT OF VELOCITY/SCALAR.
       status = dswref(fname,6+nscalar)
       status = dssdims(vel_rank,vel_length)
       status = dssnt(DFNT_FLOAT32)
       status = dsslens(1,4,5,9)
       status = dssdast('p','hU','f10.5','cartesian')
       status = dssdisc(1,vel_length(0),vel_xcoord)
       status = dssdisc(2,vel_length(1),vel_ycoord)
       status = dssdisc(3,vel_length(2),vel_zcoord)
       status = dsadata(fname,vel_rank,vel_length,data)
    END IF

! WRITE MAX VALUES OF VELOCITY/PRESSURE TO STANDARD OUTPUT & DEALLOCATE VARS
    IF (me == 0) THEN
       WRITE(*,997) (temp(i),i=1,4+nscalar)
       WRITE(*,998) dt
       WRITE(*,*)
997    FORMAT('MAX VALUES =',6E12.4)
998    FORMAT('DT         =',e12.4)

       DEALLOCATE(data,vel_xcoord,vel_ycoord,vel_zcoord,xindex,yindex,zindex,&
            STAT=ierr)
    END IF
       
  END SUBROUTINE output_hdf_dfsd
!=============================================================================
  SUBROUTINE generate_coord_index(n,skip,length,coord,vel_coord,vel_index)
    IMPLICIT NONE

    INTEGER n, skip, length
    REAL(KIND=prec), DIMENSION(1:n+1) :: coord
    REAL, DIMENSION(1:length) :: vel_coord
    INTEGER, DIMENSION(1:length) :: vel_index

    INTEGER i, j

    j = 1
    DO i = 1,n/skip
       vel_index(i) = j
       vel_coord(i) = coord(j)
       j = j + skip
    END DO
    vel_index(length) = n+1
    vel_coord(length) = coord(n+1)

  END SUBROUTINE generate_coord_index
END PROGRAM convert_float
