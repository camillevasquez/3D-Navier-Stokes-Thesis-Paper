MODULE read_write
  USE runparms
  USE gridparms
  USE data
  USE setup
  USE transform
  USE pass_data
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_velocity, save_velocity, output_snapshot, output_hdf_dfsd,&
       l2_norm, boundary
  
  INCLUDE 'mpif.h'

  REAL(KIND=prec) :: ts,tf
  INTEGER  i, j, k, n, x, z, jj, plane, ierr, mpierr, ioerror
CONTAINS
!=================================================================
!=========================READ VELOCITY===========================
!=================================================================
  SUBROUTINE read_velocity
    IMPLICIT NONE

    INTEGER    NX_TEST, NY_TEST, NZ_TEST, status(MPI_STATUS_SIZE), &
         LAGR_TEST, ndynscalar_test, bundle(3)
    REAL(KIND=prec)  LX_TEST, LZ_TEST, RE_TEST, DT_TEST
    LOGICAL           ABORT
    COMPLEX(KIND=prec), DIMENSION(1:nz,1:ny+1,1:3) :: slice
    REAL(KIND=prec), DIMENSION(1:nx,1:nz)          :: scalar_slice
    REAL(KIND=prec), DIMENSION(:,:,:), ALLOCATABLE :: dyn_slice

!  OPEN THE INPUT DATA FILE AND READ IN THE DATA.
    IF (me==0) THEN
       OPEN (UNIT=14,FILE='savevel.out',STATUS='OLD',FORM='UNFORMATTED') 
       READ(14) NX_TEST,NY_TEST,NZ_TEST,LX_TEST,LZ_TEST,RE_TEST,DT_TEST,T_TOTAL

!       write(*,*) NX_TEST,NY_TEST,NZ_TEST,LX_TEST,LZ_TEST,RE_TEST,DT_TEST,T_TOTAL
!       write(*,*) NX,NY,NZ,LX,LZ,RE,DT
       
!  TEST THE VARIABLES TO AVOID DIMENSION CONFLICTS.
       abort = .false.
       IF (nx .NE. nx_test) CALL inconsistent_parameter('NX',nx_test,nx)
       IF (ny .NE. ny_test) CALL inconsistent_parameter('NY',ny_test,ny)
       IF (nz .NE. nz_test) CALL inconsistent_parameter('NZ',nz_test,nz)
       IF (ABS(lx-lx_test) > lx*EPSILON(lx)) THEN
          CALL inconsistent_real_parameter('LX',lx_test,lx)
       END IF
       IF (ABS(lz-lz_test) > lz*EPSILON(lz)) THEN
          CALL inconsistent_real_parameter('LZ',lz_test,lz) 
       END IF
       IF (ABS(re-re_test) > re*EPSILON(re)) THEN
!          CALL inconsistent_real_parameter('RE',re_test,re)
          IF (me==0) THEN
             WRITE(*,990)  re, re_test
990          FORMAT('REYNOLDS NUMBER CHANGED TO',f12.4,'  FROM ',f12.4)
          END IF
       END IF
    END IF

    CALL MPI_BCAST(ABORT,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    
    IF (abort) THEN
       CALL finalize_fft
       CALL finalize_mpi
       CALL deallocate_data_arrays
       STOP 'INCONSISTENT PARAMETERS'
    END IF

! INITIALIZE AND PASS TIME STEP AND INITIAL TIME
    IF (DT_FLAG == 1) THEN
       IF (me == 0) DT = DT_TEST
       CALL MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
    END IF
    CALL MPI_BCAST(t_total,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
       
    IF (me == 0) WRITE(*,*) 'READ OLD RUN PARAMETERS'

! READ DATA FIELD, SLICE-BY-SLICE IN THE X-DIRECTION.
    IF (me == 0) THEN
       DO x = 1,nx_proc(me)        ! READ SLICE IN ON SERVER NODE
          READ(14) slice           ! PUT SLICE INTO DATA1 
          data1(1:nz,x,1:ny+1,1:3) = slice(1:nz,1:ny+1,1:3)
       END DO
       WRITE(*,*) 'DONE READING OLD VELOCITY FOR PROCESSOR ', me
       DO n = 1,nproc-1
          ts = MPI_WTIME()
          DO x = 1,nx_proc(n)
             READ(14) slice      ! READ SLICE IN ON SERVER NODE
             CALL MPI_SEND(slice,nz*(ny+1)*3,MPI_DOUBLE_COMPLEX, &
                  n,x_start(n)+x,MPI_COMM_WORLD,mpierr)
          END DO
          tf = MPI_WTIME()
          WRITE(*,*) 'PASSED OLD VELOCITY TO PROC ', n, ' IN ', REAL(tf-ts)
       END DO
    ELSE
       DO x = 1,nx_proc(me)
          CALL MPI_RECV(slice,nz*(ny+1)*3,MPI_DOUBLE_COMPLEX, &
               0,x_start(me)+x,MPI_COMM_WORLD,status,mpierr)
          data1(1:nz,x,1:ny+1,1:3) = slice(1:nz,1:ny+1,1:3)
       END DO
    END IF

! READ SCALAR DATA FIELD, SLICE-BY-SLICE IN THE X-DIRECTION.
    DO i = 1,nscalar
       DO k = 1,ny+1
          IF (me == 0) READ(14) scalar_slice
          CALL MPI_SCATTERV(scalar_slice,nx*nz_small_proc,nx*z_small_start,&
               MPI_DOUBLE_PRECISION, &
               data_scalar(1,1,k,i),nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,&
               0, MPI_COMM_WORLD, mpierr)
       END DO
    END DO

    IF (me .eq. 0) WRITE(*,*) 'READ SCALAR FIELD'

    IF (subgrid == 1) THEN
! CHECK TO SEE WHETHER THE LAGRANGIAN DYNAMIC COEFFICIENTS ARE INCLUDED IN 
! THE FILE AND SHARE INFO WITH OTHER PROCESSORS.
       IF (me == 0) THEN
          READ(UNIT=14,IOSTAT=ioerror) lagr_test, ndynscalar_test
          bundle(1) = ioerror
          bundle(2) = lagr_test
          bundle(3) = ndynscalar_test
       END IF
       CALL MPI_BCAST(bundle,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
       IF (me /= 0) THEN
          ioerror = bundle(1)
          lagr_test = bundle(2)
          ndynscalar_test = bundle(3)
       END IF
! IF LAGRANGIAN COEFFICIENTS ARE INCLUDED< READ THEM IN.
       IF (ioerror == 0) THEN
          CALL MPI_BCAST(lagr_test,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
          CALL MPI_BCAST(ndynscalar_test,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
          IF (me == 0) WRITE(*,*) lagr_test, ndynscalar_test
          IF (lagr_test == 1) THEN
! READ DATA FIELD, SLICE-BY-SLICE IN THE z-DIRECTION.
             ALLOCATE(dyn_slice(1:nx,1:ny,1:2+2*ndynscalar_test),STAT=ierr)
             IF (me == 0) THEN
                DO z = 1,nz_small_proc(me)! READ SLICE IN ON SERVER NODE
                   READ(14) dyn_slice     ! PUT SLICE INTO DATA1 ON SERVER NODE
                   data_dyn(1:nx,z,1:ny,1:2+2*ndynscalar_test) &
                        = dyn_slice(1:nx,1:ny,1:2+2*ndynscalar_test)
                END DO

                DO n = 1,nproc-1
                   DO z = 1,nz_small_proc(n)
                      READ(14) dyn_slice  ! READ SLICE IN ON SERVER NODE
                      CALL MPI_SEND(dyn_slice,nx*ny*(2+2*ndynscalar_test), &
                           MPI_DOUBLE_PRECISION,n,z_small_start(n)+z, &
                           MPI_COMM_WORLD,mpierr)
                   END DO
                END DO
             ELSE
                DO z = 1,nz_small_proc(me)
                   CALL MPI_RECV(dyn_slice,nx*ny*(2+2*ndynscalar_test), &
                        MPI_DOUBLE_PRECISION,0,z_small_start(me)+z,&
                        MPI_COMM_WORLD,status,mpierr)
                   data_dyn(1:nx,z,1:ny,1:2+2*ndynscalar_test) &
                        = dyn_slice(1:nx,1:ny,1:2+2*ndynscalar_test)
                END DO
             END IF
             DEALLOCATE(dyn_slice,STAT=ierr)
! WHERE DATA IS SPLIT ACROSS PROCESSORS, PASS INFORMATION ABOUT 
! DYNAMIC COEFFICIENTS AT NEIGHBORING POINTS.  THIS IS USED WHEN
! THE DYNAMIC COEFFICIENTS ARE "ADVECTED" BY THE LOCAL VELOCITY.
             DO i = 1,2+2*ndynscalar_test
                DO plane = 1,ny
                   CALL pass_neighbors( &
                        nx,local_nz_small,1,data_dyn(0,0,plane,i))
                END DO
             END DO

          END IF
       END IF
    END IF

    IF (me == 0) CLOSE (UNIT = 14, STATUS = 'KEEP')

  CONTAINS
!------------------------------------------------------------------
    SUBROUTINE inconsistent_parameter(which,what_test,what)
      IMPLICIT NONE
      CHARACTER(LEN=2) which
      INTEGER          what_test, what

      write(*,1000)
      write(*,1010) which, what_test, what

      abort = .true.

1000  FORMAT(' ',' ')
1010  FORMAT(' ','ABORTING ... ',A2,' HAS BEEN CHANGED FROM ',I4,' TO ',I4)
    END SUBROUTINE inconsistent_parameter
!------------------------------------------------------------------
    SUBROUTINE inconsistent_real_parameter(which,what_test,what)
      IMPLICIT NONE
      CHARACTER(LEN=2) which
      REAL(kind=prec)  what_test, what

      write(*,1000)
      write(*,1020) which, what_test, what
      abort = .true.

1000  FORMAT(' ',' ')
1020  FORMAT(' ','ABORTING ... ',A2,' HAS BEEN CHANGED FROM ',&
           F10.5,' TO ',F10.5)
    END SUBROUTINE inconsistent_real_parameter
  END SUBROUTINE read_velocity
!=================================================================
!=========================SAVE VELOCITY===========================
!=================================================================
  SUBROUTINE SAVE_VELOCITY
    IMPLICIT NONE
    INTEGER  status(MPI_STATUS_SIZE), LAGR_TEST, ndynscalar_test
    COMPLEX(KIND=prec), DIMENSION(nz,ny+1,3) :: slice
    REAL(KIND=prec), DIMENSION(1:nx,1:nz)  :: scalar_slice
    REAL(KIND=prec), DIMENSION(:,:,:), ALLOCATABLE   :: dyn_slice

!  OPEN THE NEW DATA FILE AND WRITE IN THE DATA.
    IF (me == 0) THEN
       OPEN (UNIT=13,FILE='newvel.out', FORM='UNFORMATTED')
       REWIND(13)
       WRITE(13) NX, NY, NZ, LX, LZ, RE, DT, T_TOTAL
    END IF

! WRITE DATA FIELD, SLICE-BY-SLICE IN THE X-DIRECTION.
    DO n = 0,nproc-1
       DO x = 1,nx_proc(n)
          IF ((n /= 0) .and. (me == n)) THEN  ! SEND SLICE TO SERVER NODE
             slice(1:nz,1:ny+1,1:3) = data1(1:nz,x,1:ny+1,1:3)
             CALL MPI_SEND(slice,nz*(ny+1)*3,MPI_DOUBLE_COMPLEX, &
                  0,n*x,MPI_COMM_WORLD,mpierr)
          ELSEIF (me == 0) THEN
             IF (n == 0) THEN        ! COPY DATA1 INTO SLICE ON SERVER NODE
                slice(1:nz,1:ny+1,1:3) = data1(1:nz,x,1:ny+1,1:3)
             ELSE                    ! RECV SLICE FROM NODE n
                CALL MPI_RECV(slice,nz*(ny+1)*3,MPI_DOUBLE_COMPLEX, &
                     n,n*x,MPI_COMM_WORLD,status,mpierr)
             END IF
             write(13) slice         ! WRITE OUT SLICE ON SERVER NODE
          END IF
       END DO
    END DO

! WRITE SCALAR DATA FIELD, SLICE-BY-SLICE IN THE X-DIRECTION.
    DO i = 1,nscalar
       DO k = 1,ny+1
          CALL MPI_GATHERV( &
               data_scalar(1,1,k,i),nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,&
               scalar_slice,nx*nz_small_proc,nx*z_small_start, &
               MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          IF (me == 0) WRITE(13) scalar_slice
       END DO
    END DO

! WRITE OUT DYNAMIC COEFFICIENTS IF LAGRANGIAN AVERAGING IS USED.
    IF ((subgrid == 0) .or. (lagr_avg == 0)) THEN
       lagr_test = 0
       ndynscalar_test = 0
       write(13) lagr_test,ndynscalar_test
    ELSEIF ((subgrid == 1) .and. (lagr_avg == 1)) THEN
       write(13) lagr_avg,ndynscalar       

! WRITE DATA FIELD, SLICE-BY-SLICE IN THE X-DIRECTION.
       ALLOCATE(dyn_slice(1:nx,1:ny,1:2+2*ndynscalar),STAT=ierr)
       DO n = 0,nproc-1
          DO z = 1,nz_small_proc(n)
             IF ((n /= 0) .and. (me == n)) THEN  ! SEND SLICE TO SERVER NODE
                dyn_slice(1:nx,1:ny,1:2+2*ndynscalar)  &
                     = data_dyn(1:nx,z,1:ny,1:2+2*ndynscalar)
                CALL MPI_SEND(dyn_slice,nx*ny*(2+2*ndynscalar), &
                     MPI_DOUBLE_PRECISION,0,n*z,MPI_COMM_WORLD,mpierr)
             ELSEIF (me == 0) THEN
                IF (n == 0) THEN  ! COPY DATA_DYN INTO DYN_SLICE ON SERVER NODE
                   dyn_slice(1:nx,1:ny,1:2+2*ndynscalar) = &
                        data_dyn(1:nx,z,1:ny,1:2+2*ndynscalar)
                ELSE                    ! RECV SLICE FROM NODE n
                   CALL MPI_RECV(dyn_slice,nx*ny*(2+2*ndynscalar), &
                        MPI_DOUBLE_PRECISION,n,n*z,MPI_COMM_WORLD,&
                        status,mpierr)
                END IF
                write(13) dyn_slice ! WRITE OUT dyn_slice ON SERVER NODE
             END IF
          END DO
       END DO
       DEALLOCATE(dyn_slice,STAT=ierr)
    END IF

    IF (me == 0) CLOSE (UNIT = 13, STATUS = 'KEEP')

  END SUBROUTINE SAVE_VELOCITY
!=================================================================
!=========================OUTPUT HDF DFSD=========================
!=================================================================
  SUBROUTINE output_hdf_dfsd
    IMPLICIT NONE
    

    COMPLEX(KIND=prec), DIMENSION(nz,local_nx)   :: cu
    REAL(KIND=prec), DIMENSION(:,:), ALLOCATABLE :: u1,u2
    REAL(KIND=prec) :: temp(4+nscalar), tmp2(4+nscalar)
    REAL(KIND=prec) :: poss_max, poss_min, xfine(nx2+1),zfine(nz2+1)
    REAL(KIND=prec) :: temp2(ny,1+ndynscalar), temp3(ny,1+ndynscalar)
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
    
    REAL :: data_vel_out(nx,ny+1,nz,3)

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
    CHARACTER(LEN=4) time_name
    CHARACTER(LEN=1) vel_name(1:5), dim_name(0:2)
    CHARACTER(LEN=14) FNAME
    CHARACTER(LEN=30) vtk_string
    character :: lf*1

    lf = achar(10) ! line feed character

    ALLOCATE(u1(nx,local_nz_small),u2(nx,nz),STAT=ierr)
    IF (ierr /= 0) WRITE(*,*) 'CANNOT ALLOCATE u1 AND u2'
! SET UP LENGTH OF OUTPUT ARRAYS AND INDEXING INTO FULL DATA ARRAYS (nx x nz)
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
       time_out(2) = d_um_dy(1)
       time_out(3) = d_um_dy(ny)
       time_out(4) = p_grad

! CHECK NUMBER OF ENTRIES IN NAMEHOLDER HDF FILE.
!       fname="t_xxx.hdf"
       fname="t_xxx.vtk"
       write(fname(3:5),"(i3.3)") step/output_step

!       fname_rank = 1
!       fname_length = 9
!       status = dssdims(fname_rank,fname_length)
!       status = dssnt(DFNT_CHAR8)
!       status = dsadata('datafiles.hdf',fname_rank,fname_length,fname)

! WRITE OUT TIME STAMP AS FIRST ENTRY IN FILE.
 !      status = dssdims(time_rank,time_length)
 !      status = dssnt(DFNT_FLOAT32)
 !      status = dsslens(4,4,5,4)
 !      status = dssdast('time','hL/U','f10.5','none')
 !      status = dspdata(fname,time_rank,time_length,time_out)

! the vtk file output is something like this
       open(unit=12,file=fname,access='stream',status='replace',convert='BIG_ENDIAN')
       write(12) "# vtk DataFile Version 2.0"//lf
       write(12) "(This is header line) vtk rectilinear grid with binary data"//lf
       write(12) 'BINARY'//lf
       write(12) 'DATASET RECTILINEAR_GRID'//lf
       vtk_string='DIMENSIONS'
       write(vtk_string(12:26),"(I5,I5,I5)") nx, ny+1, nz
       write(12) vtk_string//lf
       vtk_string='X_COORDINATES     float'
       write(vtk_string(15:17),"(I3)") nx
       write(12) vtk_string//lf
       write(12)  (real(lx*dble(x-1)/dble(nx)),x = 1,nx), lf
       vtk_string='Y_COORDINATES     float'
       write(vtk_string(15:17),"(I3)") ny+1
       write(12) vtk_string//lf
       write(12)  (real(yh(plane)),plane = 1,ny+1), lf
       vtk_string='Z_COORDINATES     float'
       write(vtk_string(15:17),"(I3)") nz
       write(12) vtk_string//lf
       write(12)  (real(lz*dble(z-1)/dble(nz)),z=1,nz), lf
!       write(*,*)  "Z_COORDINATES", (lz*dble(z-1)/dble(nz),z=1,nz)

    END IF

    DO i = 1,3
       IF (me == 0) poss_max = zero
       IF (me == 0) poss_min = zero
       u1 = zero
       u2 = zero
       DO jj = 1,ylength

          plane = yindex(jj)
          IF (i == 2) THEN
             IF (plane == 1) THEN
                cu = data1(1:nz,1:local_nx,1,i)
             ELSEIF (plane == ny+1) THEN
                cu = data1(1:nz,1:local_nx,ny,i)
             ELSE
                IF (yskip == 1) THEN
                   cu = half*(data1(1:nz,1:local_nx,plane-1,i) &
                        + data1(1:nz,1:local_nx,plane,i))
                ELSE
                   cu = 0.125d0*data1(1:nz,1:local_nx,plane-2,i) &
                      + 0.375d0*data1(1:nz,1:local_nx,plane-1,i) &
                      + 0.375d0*data1(1:nz,1:local_nx,plane,i) &
                      + 0.125d0*data1(1:nz,1:local_nx,plane+1,i)
                END IF
             END IF
          ELSE
             IF ((plane == 1) .or. (plane == ny+1)) THEN
                cu = data1(1:nz,1:local_nx,plane,i)
             ELSE
                IF (yskip == 1) THEN
                   cu = data1(1:nz,1:local_nx,plane,i)
                ELSE
                   cu = 0.25d0*data1(1:nz,1:local_nx,plane-1,i) &
                      + 0.50d0*data1(1:nz,1:local_nx,plane,i) &
                      + 0.25d0*data1(1:nz,1:local_nx,plane+1,i)
                END IF
             END IF
          END IF
          IF (xskip .gt. 1) CALL filterx(cu)
          IF (zskip .gt. 1) CALL filterz(cu)
          CALL small_fft_complex_to_real(cu,u1)
          CALL MPI_GATHERV(u1,nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,   &
               u2,nx*nz_small_proc(0:nproc-1),nx*z_small_start(0:nproc-1), &
               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
          IF (me == 0) THEN
             poss_max = MAX(MAXVAL(u2),poss_max)
             poss_min = MIN(MINVAL(u2),poss_min)
             CALL upload_data(jj)
          END IF
       END DO
       IF (me == 0) temp(i) = poss_max
       IF (me == 0) tmp2(i) = poss_min
!!$       CALL MPI_Reduce(poss_max,temp(i),1,MPI_DOUBLE_PRECISION, &
!!$            MPI_MAX,0,MPI_COMM_WORLD,mpierr)

! WRITE THE CURRENT COMPONENT OF VELOCITY/SCALAR INTO THE HDF FILE
       IF (me == 0) THEN
!           status = dssdims(vel_rank,vel_length)   ! CREATE AN ARRAY
!           status = dssnt(DFNT_FLOAT32)            ! SPECIFY DATA TYPE
!           status = dsslens(1,4,5,9)               ! SPECIFY LABELS
!           status = dssdast(vel_name(i),'U','f10.5','cartesian')
!           status = dssdisc(1,vel_length(0),vel_xcoord) ! DIMENSIONS
!           status = dssdisc(2,vel_length(1),vel_ycoord)
!           status = dssdisc(3,vel_length(2),vel_zcoord)
!           status = dsadata(fname,vel_rank,vel_length,data) ! WRITE ARRAY
! this is where vtk output should happen

          data_vel_out(:,:,:,i)=real(data)

       END IF

    END DO

    IF (me == 0) then
       vtk_string='POINT_DATA '
       write(vtk_string(12:20),"(I9)") nx*(ny+1)*nz
       write(12) vtk_string//lf
       write(12) 'VECTORS flowvec float'//lf
       write(12) ((((data_vel_out(x,plane,z,i),i=1,3),x=1,nx),plane=1,ny+1),z=1,nz),lf
    
! close the vtk file
       close(12)
    end IF


! SAVE THE PRESSURE TO THE DATA FILE AS WELL.
    IF (me == 0) poss_max = zero
    IF (me == 0) poss_min = zero
    u1 = zero
    u2 = zero
    DO jj = 1,ylength
       plane = yindex(jj)
       cu = data4(plane,1:nz,1:local_nx)
       CALL small_fft_complex_to_real(cu,u1)
       CALL MPI_GATHERV(u1,nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,   &
            u2,nx*nz_small_proc(0:nproc-1),nx*z_small_start(0:nproc-1), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
       IF (me == 0) THEN
          poss_max = MAX(MAXVAL(u2),poss_max)
          poss_min = MIN(MINVAL(u2),poss_min)
          CALL upload_data(jj)
       END IF
    END DO
    IF (me == 0) temp(4+nscalar) = poss_max
    IF (me == 0) tmp2(4+nscalar) = poss_min
!!$    CALL MPI_Reduce(poss_max,temp(4+nscalar),1,MPI_DOUBLE_PRECISION, &
!!$            MPI_MAX,0,MPI_COMM_WORLD,mpierr)

! WRITE THE PRESSURE TO THE DATA FILE.
    IF (me == 0) THEN
! CREATE AN ARRAY FOR THE PRESSURE
!        status = dssdims(vel_rank,vel_length)
!        status = dssnt(DFNT_FLOAT32)
!        status = dsslens(1,4,5,9)
!        status = dssdast('p','hU','f10.5','cartesian')
!        status = dssdisc(1,vel_length(0),vel_xcoord)
!        status = dssdisc(2,vel_length(1),vel_ycoord)
!        status = dssdisc(3,vel_length(2),vel_zcoord)
!        status = dsadata(fname,vel_rank,vel_length,data)
! this is where vtk output should happen
    END IF

! OUTPUT SCALAR FIELD TO HDF FILE.
    DO i = 1,nscalar
       IF (me == 0) poss_max = zero
       IF (me == 0) poss_min = zero
       u1 = zero
       u2 = zero
       DO jj = 1,ylength
          plane = yindex(jj)
          CALL MPI_GATHERV( &
               data_scalar(1,1,plane,i),nx*nz_small_proc(me), &
               MPI_DOUBLE_PRECISION,&
               u2,nx*nz_small_proc,nx*z_small_start,MPI_DOUBLE_PRECISION, &
               0, MPI_COMM_WORLD, mpierr)
          IF (me == 0) THEN
             poss_max = MAX(MAXVAL(u2),poss_max)
             poss_min = MIN(MINVAL(u2),poss_min)
             CALL upload_data(jj)
          END IF
       END DO
       IF (me == 0) temp(3+i) = poss_max
       IF (me == 0) tmp2(3+i) = poss_min

! WRITE THE SCALAR INTO THE HDF FILE
       IF (me == 0) THEN
!           status = dssdims(vel_rank,vel_length)   ! CREATE AN ARRAY
!           status = dssnt(DFNT_FLOAT32)            ! SPECIFY DATA TYPE
!           status = dsslens(1,4,5,9)               ! SPECIFY LABELS
!           status = dssdast(vel_name(i),'U','f10.5','cartesian')
!           status = dssdisc(1,vel_length(0),vel_xcoord) ! DIMENSIONS
!           status = dssdisc(2,vel_length(1),vel_ycoord)
!           status = dssdisc(3,vel_length(2),vel_zcoord)
!           status = dsadata(fname,vel_rank,vel_length,data) ! WRITE ARRAY
! this is where vtk output should happen
       END IF
    END DO

    IF (subgrid == 1) THEN
       IF (lagr_avg == 0) THEN
          
! WRITE THE DYNAMIC COEFFICIENT TO THE DATA FILE.
          IF (me == 0) THEN
             ALLOCATE(tmp_dyncoef(ylength),STAT=ierr)
             DO i = 1,1+ndynscalar
                DO jj = 1,ylength
                   tmp_dyncoef(jj) = dyn_visc(yindex(jj),i)
                END DO
! CREATE AN ARRAY FOR THE CURRENT COMPONENT OF VELOCITY/SCALAR.
!                 status = dssdims(1,vel_length(1))
!                 status = dssnt(DFNT_FLOAT32)
!                 status = dssdisc(1,vel_length(1),vel_ycoord)
!                 status = dsadata(fname,1,vel_length(1),tmp_dyncoef)
! this is where vtk output should happen
                temp2(1:ny,i) = dyn_visc(1:ny,i)
             END DO
             DEALLOCATE(tmp_dyncoef)
          END IF


       ELSEIF (lagr_avg == 1) THEN
          DO i = 1,2+2*ndynscalar
             DO jj = 1,ylength
                plane = yindex(jj)
                IF (plane == 1) THEN
                   u1 = data_dyn(1:nx,1:local_nz_small,1,i)
                ELSEIF (plane == ny+1) THEN
                   u1 = data_dyn(1:nx,1:local_nz_small,ny,i)
                ELSE
                   u1 = half*(data_dyn(1:nx,1:local_nz_small,plane-1,i)&
                        + data_dyn(1:nx,1:local_nz_small,plane,i))
                END IF

                CALL MPI_GATHERV(                 &
                     u1, nx*nz_small_proc(me), MPI_DOUBLE_PRECISION, &
                     u2, nx*nz_small_proc,     nx*z_small_start,     &
                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
                IF (me == 0) CALL upload_data(jj)
             END DO

! WRITE THE DYNAMIC COEFFICIENT TO THE DATA FILE.
             IF (me == 0) THEN
! CREATE AN ARRAY FOR THE CURRENT COMPONENT OF VELOCITY/SCALAR.
!                 status = dssdims(vel_rank,vel_length)
!                 status = dssnt(DFNT_FLOAT32)
!                 status = dsslens(2,4,5,9)
!                 SELECT CASE (i)
!                 CASE (1)
!                    status = dssdast('lm','hU','f10.5','cartesian')
!                 CASE (2)
!                    status = dssdast('mm','hU','f10.5','cartesian')
!                 CASE (3)
!                    status = dssdast('lt','hU','f10.5','cartesian')
!                 CASE (4)
!                    status = dssdast('mt','hU','f10.5','cartesian')
!                 CASE (5)
!                    status = dssdast('ls','hU','f10.5','cartesian')
!                 CASE (6)
!                    status = dssdast('ms','hU','f10.5','cartesian')
!                 END SELECT
!                 status = dssdisc(1,vel_length(0),vel_xcoord)
!                 status = dssdisc(2,vel_length(1),vel_ycoord)
!                 status = dssdisc(3,vel_length(2),vel_zcoord)
!                 status = dsadata(fname,vel_rank,vel_length,data)
! this is where vtk output should happen
             END IF
          END DO
          
          DO i = 1,1+ndynscalar
             DO jj = 1,ny
                temp3(jj,i) = MAXVAL(data_dyn(1:nx,1:local_nz_small,jj,2*i-1) &
                     /data_dyn(1:nx,1:local_nz_small,jj,2*i))
             END DO
          END DO
          CALL MPI_Reduce(temp3,temp2,ny*(1+ndynscalar),MPI_DOUBLE_PRECISION, &
               MPI_MAX, 0, MPI_COMM_WORLD,mpierr)
       END IF
    END IF


! WRITE MAX VALUES OF VELOCITY/PRESSURE TO STANDARD OUTPUT & DEALLOCATE VARS
    IF (me == 0) THEN
       WRITE(*,996) (temp(i),i=1,4+nscalar)
       WRITE(*,997) (tmp2(i),i=1,4+nscalar)
       WRITE(*,998) dt
       WRITE(*,*)
996    FORMAT('MAX VALUES =',6E12.4)
997    FORMAT('MIN VALUES =',6E12.4)
998    FORMAT('DT         =',e12.4)

       IF (subgrid == 1) THEN
          DO PLANE = 1,ny
             WRITE(*,999) plane, (temp2(plane,i),i=1,1+ndynscalar)
          END DO
          WRITE(*,*)
       END IF
999    FORMAT('PLANE = ',i4,'  MAX DYN COEF = ',4F12.8)
       DEALLOCATE(data,STAT=ierr)
    END IF
    DEALLOCATE(vel_xcoord,vel_ycoord,vel_zcoord,xindex,yindex,zindex,&
         STAT=ierr)
    DEALLOCATE(u1,u2)

  CONTAINS
!=============================================================================
    SUBROUTINE filterz(cvel)
      IMPLICIT NONE
      COMPLEX(KIND=prec), DIMENSION(nz,local_nx) :: cvel
      REAL(KIND=prec) :: kmax
      INTEGER x,z

      kmax = DBLE(nz/2)*two*pi/lz
      DO x = 1,local_nx
         DO z = 1,nz
            cvel(z,x) = cvel(z,x)*half*(one + COS(pi*kz(z)/kmax))
         END DO
      END DO

    END SUBROUTINE filterz
!=============================================================================
    SUBROUTINE filterx(cvel)
      IMPLICIT NONE
      COMPLEX(KIND=prec), DIMENSION(nz,local_nx) :: cvel
      REAL(KIND=prec) :: kmax
      INTEGER x,z

      kmax = DBLE(nx/2)*two*pi/lx
      DO x = 1,local_nx
         DO z = 1,nz
            cvel(z,x) = cvel(z,x)*half*(one + COS(pi*kx_me(x)/kmax))
         END DO
      END DO

    END SUBROUTINE filterx
!=============================================================================
    SUBROUTINE upload_data(plane)
      IMPLICIT NONE
      INTEGER plane

      DO z = 1,zlength-1
         DO x = 1,xlength-1
            data(x,plane,z) = u2(xindex(x),zindex(z))
         END DO
      END DO
      data(xlength,plane,1:zlength-1) = data(1,plane,1:zlength-1)
      data(1:xlength,plane,zlength)   = data(1:xlength,plane,1)

!!$      data(1:nx,jj,1:nz)   = u2(1:nx,1:nz)
!!$      data(nx+1,jj,1:nz)   = data(1,jj,1:nz)
!!$      data(1:nx+1,jj,nz+1) = data(1:nx+1,jj,1)

    END SUBROUTINE upload_data
!=============================================================================
    SUBROUTINE generate_coord_index(n,skip,length,coord,vel_coord,vel_index)
      IMPLICIT NONE

      INTEGER n, skip, length
      REAL(KIND=prec), DIMENSION(1:n+1) :: coord
      REAL, DIMENSION(1:length) :: vel_coord
      INTEGER, DIMENSION(1:length) :: vel_index

      j = 1
      DO i = 1,n/skip
         vel_index(i) = j
         vel_coord(i) = coord(j)
         j = j + skip
       END DO
       vel_index(length) = n+1
       vel_coord(length) = coord(n+1)

     END SUBROUTINE generate_coord_index
  END SUBROUTINE output_hdf_dfsd
!=================================================================
!=======================OUTPUT SNAPSHOT===========================
!=================================================================
  SUBROUTINE output_snapshot
    IMPLICIT NONE
!    INCLUDE 'hdf.f90'
    REAL(KIND=prec) :: vel_xcoord(nx2), vel_ycoord(ny+1), vel_zcoord(nz2)
    REAL(KIND=prec), DIMENSION(6) :: time
    ! DECLARE HDF VARIABLES
    INTEGER, PARAMETER :: vel_rank = 3
    INTEGER, DIMENSION(0:vel_rank-1)::vel_length
    INTEGER  time_rank, time_length
    CHARACTER(LEN=4) time_name
    CHARACTER(LEN=1) vel_name(1:5), dim_name(0:2)

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    CALL setup_dims_time

    IF (snap_out == 1) THEN
       CALL output_snap('snap') ! SAVE VELOCITY/SCALAR FIELD TO A FILE.
       IF (subgrid == 1) CALL output_snap('visc') ! SAVE DYNAMIC MODEL COEF.
    ELSEIF (snap_out == 2) THEN
       CALL compute_mean ! COMPUTE RUNNING MEAN OF VELOCITY/SCALAR/DYN. COEF.
       CALL compute_uu   ! COMPUTE RUNNING MEAN OF VELOCITY/SCALAR/DYN. COEF.
    END IF

  CONTAINS
!=========================================================================
    SUBROUTINE output_snap(name)
      IMPLICIT NONE
      
      REAL(KIND=prec), DIMENSION(:,:,:), ALLOCATABLE :: data
      COMPLEX(KIND=prec), DIMENSION(nz,local_nx) :: cu
      REAL(KIND=prec), DIMENSION(nx2,local_nz)   :: u1
      REAL(KIND=prec), DIMENSION(nx2,nz2)         :: u2

      REAL(KIND=prec) :: temp(4+nscalar), poss_max

      ! DECLARE HDF FUNCTIONS
      INTEGER  dsadata, dssdims, dssdisc, dspdata, dssdast, dsslens
      INTEGER  dssnt, dsnum
      ! DECLARE HDF VARIABLES
      INTEGER  num_datasets, fname_rank, fname_length, status
      CHARACTER(LEN=16) FNAME, fnameholder
      CHARACTER(LEN=4) name

    IF (me == 0) THEN
       ALLOCATE(data(nx2,nz2,ny+1), STAT=ierr)
       IF (ierr /= 0) THEN
          STOP 'CANNOT ALLOCATE ARRAY FOR DOUBLE PRECISION SNAPSHOT'
       END IF
    END IF

      ! SET UP ARRAY INFO.
      IF (me==0) THEN
         ! CHECK NUMBER OF ENTRIES IN NAMEHOLDER HDF FILE.
         fname_rank = 1
         fname_length = 16
         fnameholder = name // 'files.hdf'
!         num_datasets = dsnum(fnameholder)
         write(*,*) 'NUMBER OF OLD DATASETS = ', num_datasets
         IF (num_datasets <= 0) THEN
            num_datasets = 1
         ELSE
            num_datasets = num_datasets+1
         END IF
         FNAME = name // &
              char(48 + int(t_total/1000))// &
              char(48 + int(mod(int(t_total/100),10)))// & 
              char(48 + mod(int(t_total/10),10))// &
              char(48 + mod(int(t_total),10))//&
              '_'// &
              char(48 + mod(int(t_total*10),10))//&
              char(48 + mod(int(t_total*100),10))//&
              char(48 + mod(int(t_total*1000),10))//&
              '.hdf'
!          status = dssdims(fname_rank,fname_length)
!          status = dssnt(DFNT_CHAR8)
!          status = dsadata(fnameholder,fname_rank,fname_length,fname)

!          ! WRITE OUT TIME STAMP AS FIRST ENTRY IN FILE.
!          status = dssdims(time_rank,time_length)
!          status = dssnt(DFNT_FLOAT64)
!          status = dsslens(4,4,5,4)
!          status = dssdast('time','hL/U','f10.5','none')
!          status = dspdata(fname,time_rank,time_length,time)
! this is where vtk output should happen

      END IF

      IF ((name == 'snap') .or. (name == 'd_dt')) THEN

         DO i = 1,3
            IF (me == 0) poss_max = zero
            DO plane = 1,ny+1
               cu = data1(1:nz,1:local_nx,plane,i)
               CALL xz_fft_complex_to_real(cu,u1)
               CALL MPI_GATHERV( &
                    u1,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
                    u2,nx2*nz_proc(0:nproc-1),nx2*z_start(0:nproc-1), &
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
               IF (me == 0) THEN
                  poss_max = MAX(MAXVAL(u2),MAXVAL(-u2),poss_max)
                  data(:,:,plane) = u2
               END IF
            END DO
            IF (me == 0) temp(i) = poss_max
!!$            CALL MPI_Reduce(poss_max,temp(i),1,MPI_DOUBLE_PRECISION,&
!!$                 MPI_MAX,0,MPI_COMM_WORLD,mpierr)

            ! WRITE THE CURRENT COMPONENT OF VELOCITY/SCALAR INTO THE HDF FILE
            IF (me == 0) THEN
!                status = dssdims(vel_rank,vel_length)   ! CREATE AN ARRAY
!                status = dssnt(DFNT_FLOAT64)            ! SPECIFY DATA TYPE
!                status = dsslens(1,4,5,9)               ! SPECIFY LABELS
!                status = dssdast(vel_name(i),'U','f10.5','cartesian')
!                status = dssdisc(1,vel_length(0),vel_xcoord) ! DIMENSIONS
!                status = dssdisc(2,vel_length(1),vel_zcoord)
!                status = dssdisc(3,vel_length(2),vel_ycoord)
!                status = dsadata(fname,vel_rank,vel_length,data) ! WRITE ARRAY
! this is where vtk output should happen
            END IF
         END DO

         DO i = 1,nscalar
            IF (me == 0) poss_max = zero
            DO plane = 1,ny+1
               CALL MPI_GATHERV( &
                    data_scalar(1,1,plane,i),nx*nz_small_proc(me), &
                    MPI_DOUBLE_PRECISION,   &
                    u2,nx*nz_small_proc,nx*z_small_start,&
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
               IF (me == 0) THEN
                  poss_max = MAX(MAXVAL(u2),MAXVAL(-u2),poss_max)
                  data(:,:,plane) = u2
               END IF
            END DO
            IF (me == 0) temp(3+i) = poss_max
!!$            CALL MPI_Reduce(poss_max,temp(i),1,MPI_DOUBLE_PRECISION,&
!!$                 MPI_MAX,0,MPI_COMM_WORLD,mpierr)

            ! WRITE THE CURRENT COMPONENT OF VELOCITY/SCALAR INTO THE HDF FILE
            IF (me == 0) THEN
!                status = dssdims(vel_rank,vel_length)   ! CREATE AN ARRAY
!                status = dssnt(DFNT_FLOAT64)            ! SPECIFY DATA TYPE
!                status = dsslens(1,4,5,9)               ! SPECIFY LABELS
!                status = dssdast(vel_name(3+i),'U','f10.5','cartesian')
!                status = dssdisc(1,vel_length(0),vel_xcoord) ! DIMENSIONS
!                status = dssdisc(2,vel_length(1),vel_zcoord)
!                status = dssdisc(3,vel_length(2),vel_ycoord)
!                status = dsadata(fname,vel_rank,vel_length,data) ! WRITE ARRAY
            END IF
         END DO

         ! WRITE MAX VALUES OF VELOCITY/PRESSURE TO STANDARD OUTPUT
         IF (me == 0) THEN
            WRITE(*,998) (temp(i),i=1,3+nscalar)
            WRITE(*,999) dt
998         FORMAT('MAX VALUES =',6E12.4)
999         FORMAT('DT         =',e12.4)
         END IF

      ELSEIF (name == 'visc') THEN
         IF (lagr_avg == 0) THEN

            ! WRITE THE DYNAMIC COEFFICIENT TO THE DATA FILE.
            IF (me == 0) THEN
               DO i = 1,1+ndynscalar
                  ! CREATE AN ARRAY FOR THE CURRENT DYNAMIC MODEL COEF.
!                   status = dssdims(1,vel_length(2))
!                   status = dssnt(DFNT_FLOAT64)
!                   status = dssdisc(1,vel_length(2),vel_ycoord)
!                   status = dsadata(fname,1,vel_length(2),dyn_visc(1,i))
! this is where vtk output should happen
               END DO
            END IF


         ELSEIF (lagr_avg == 1) THEN
            ! IF LAGRANGIAN AVERAGING IS USED, nx2 == nx and nz2 == nz.
            DO i = 1,2+2*ndynscalar
               DO plane = 1,ny
                  u1 = data_dyn(1:nx2,1:local_nz,plane,i)
                  CALL MPI_GATHERV(                 &
                       u1, nx2*nz_proc(me), MPI_DOUBLE_PRECISION, &
                       u2, nx2*nz_proc,     nx2*z_start,     &
                       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
                  IF (me == 0) data(:,:,plane) = u2
               END DO
               data(:,:,ny+1) = zero

               ! WRITE THE DYNAMIC COEFFICIENT TO THE DATA FILE.
               IF (me == 0) THEN
                  ! CREATE AN ARRAY FOR THE CURRENT VELOCITY/SCALAR COMPONENT.
!                   status = dssdims(vel_rank,vel_length)
!                   status = dssnt(DFNT_FLOAT64)
!                   status = dsslens(2,4,5,9)
!                   SELECT CASE (i)
!                   CASE (1)
!                      status = dssdast('lm','hU','f10.5','cartesian')
!                   CASE (2)
!                      status = dssdast('mm','hU','f10.5','cartesian')
!                   CASE (3)
!                      status = dssdast('lt','hU','f10.5','cartesian')
!                   CASE (4)
!                      status = dssdast('mt','hU','f10.5','cartesian')
!                   CASE (5)
!                      status = dssdast('ls','hU','f10.5','cartesian')
!                   CASE (6)
!                      status = dssdast('ms','hU','f10.5','cartesian')
!                   END SELECT
!                   status = dssdisc(1,vel_length(0),vel_xcoord)
!                   status = dssdisc(2,vel_length(1),vel_zcoord)
!                   status = dssdisc(3,vel_length(2),vel_ycoord)
!                   status = dsadata(fname,vel_rank,vel_length,data)
! this is where vtk output should happen
               END IF

            END DO
         END IF
      END IF

      IF (me == 0) DEALLOCATE(data)

    END SUBROUTINE output_snap
!=========================================================================
    SUBROUTINE compute_mean
      IMPLICIT NONE
      
!      INCLUDE 'hdf.f90'
    
      REAL(KIND=prec), DIMENSION(:,:,:), ALLOCATABLE :: data
      COMPLEX(KIND=prec), DIMENSION(nz,local_nx) :: cu
      REAL(KIND=prec), DIMENSION(:,:), ALLOCATABLE :: u1, u2

      REAL(KIND=prec) :: temp(4+nscalar), poss_max, tempa, tempb

      ! DECLARE HDF FUNCTIONS
      INTEGER   sfstart, sfcreate, sfwdata, sfendacc, sfend, sfselect, sfrdata
      INTEGER   sfdimid, sfsdscale, sfsdmname
      ! DECLARE HDF VARIABLES
      INTEGER, PARAMETER :: vel_rank = 3
      INTEGER, DIMENSION(0:vel_rank-1)::vel_length, vel_start, &
                                        vel_stride, vel_edges
      INTEGER  time_rank, time_length, time_start, time_stride, time_edges
      INTEGER  sd_id, sds_id, status, dim_id
      REAL(KIND=prec) :: vel_xcoord(nx), vel_ycoord(ny+1), vel_zcoord(nz), &
           vel_xcoord2(nx2), vel_zcoord2(nz2)
      INTEGER  num_datasets, time_ref
      CHARACTER(LEN=8) FNAME
      CHARACTER(LEN=4) time_name
      CHARACTER(LEN=1) u_name, dim_name(0:2)

      ALLOCATE(u1(nx,local_nz_small),u2(nx,nz),STAT=ierr)
      IF (me == 0) THEN
         ALLOCATE(data(nx,nz,ny+1),STAT=ierr)
         IF (ierr /= 0) STOP 'CANNOT ALLOCATE ARRAY FOR MEAN FIELDS'

         ! DESCRIBE VELOCITY FIELD
         vel_length(0) = nx
         vel_length(1) = nz
         vel_length(2) = ny+1
         vel_start(0:2) = 0
         vel_stride(0:2) = 1
         vel_edges(0:2) = vel_length(0:2)

         ! DESCRIBE VELOCITY FIELD
         time_rank = 1
         time_length = 6
         time_start = 0
         time_stride = 1
         time_edges = time_length
         time_name = 'time'
!!$         time_out  = real(t_total)

         ! DESCRIBE COORDINATES
         dim_name(0) = 'x'
         dim_name(1) = 'y'
         dim_name(2) = 'z'
         DO i = 1,nx
            vel_xcoord(i) = real(i-1)*lx/real(nx)
         END DO
         vel_ycoord = real(yh)
         DO i = 1,nz
            vel_zcoord(i) = real(i-1)*lz/real(nz)
         END DO
      END IF

      IF (me == 0) THEN
         fname = 'mean.hdf'
         WRITE(*,*) 'COMPUTING ENSEMBLE-AVERAGED MEAN'
!         sd_id = sfstart(FNAME,DFACC_WRITE)
         IF (sd_id .eq. -1) THEN
!            sd_id = sfstart(FNAME,DFACC_CREATE)
!            sds_id = sfcreate(sd_id,time_name,DFNT_FLOAT64, &
!                 time_rank,time_length)
            time(1) = zero
         ELSE
!            sds_id = sfselect(sd_id,0)
!            status = sfrdata(sds_id,time_start,time_stride,time_edges,time)
         END IF
         num_datasets = INT(time(1))
         WRITE(*,*) 'NUMBER OF OLD DATASETS IN mean.hdf = ', num_datasets
         tempa =     one/(time(1) + one) ! WEIGHT OF NEW VELOCITY FIELD
         tempb = time(1)/(time(1) + one) ! WEIGHT OF OLD RUNNING AVERAGE
         time(1) = time(1) + one ! time(1) HOLDS NO. OF FIELDS INCLUDED IN AVG.
!!$         write(*,991) tempa, tempb, time(1)
!!$991      FORMAT(3f8.3)

         ! WRITE OUT TIME STAMP AS FIRST ENTRY IN FILE.
!          status = sfwdata(sds_id,time_start,time_stride,time_edges,time)
!          status = sfendacc(sds_id)
! this is where vtk output should happen
      END IF

      DO i = 1,3
         IF (me == 0) THEN
            ! CREATE AN ARRAY FOR THE CURRENT COMPONENT OF VELOCITY/SCALAR.
            SELECT CASE (i)
            CASE (1)
               u_name = 'u'
            CASE (2) 
               u_name = 'v'
            CASE (3) 
               u_name = 'w'
            CASE (4) 
               u_name = 't'
            CASE (5) 
               u_name = 's'
            END SELECT
            IF (num_datasets .eq. 0) THEN
               !! IF mean.hdf DOESN'T EXIST, CREATE DATASET FOR CURRENT FIELD
                data = zero
!                sds_id = sfcreate(sd_id,u_name,DFNT_FLOAT64,vel_rank,vel_length)
            ELSE
               !! IF mean.hdf EXISTS, READ DATASET FOR CURRENT FIELD.
               data = zero
!               sds_id = sfselect(sd_id,i)
!               status = sfrdata(sds_id,vel_start,vel_stride,vel_edges,data)
               IF (status .ne. 0) WRITE(*,*) 'CANNOT READ ', u_name, &
                    ' FIELD IN mean.hdf'
            END IF
         END IF

         IF (me == 0) poss_max = zero
         DO plane = 1,ny+1
            cu = data1(1:nz,1:local_nx,plane,i)
            CALL small_fft_complex_to_real(cu,u1)
            CALL MPI_GATHERV( &
                 u1,nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,   &
                 u2,nx*nz_small_proc,nx*z_small_start,MPI_DOUBLE_PRECISION, &
                 0,MPI_COMM_WORLD,mpierr)
            IF (me == 0) THEN
               data(:,:,plane) = tempa*u2 + tempb*data(:,:,plane)
               poss_max = MAX(MAXVAL(ABS(data(:,:,plane))),poss_max)
            END IF
         END DO
         IF (me == 0) temp(i) = poss_max

         ! WRITE THE CURRENT COMPONENT OF VELOCITY/SCALAR INTO THE HDF FILE
         IF (me == 0) THEN
!            status = sfwdata(sds_id,vel_start,vel_stride,vel_edges,data)
!            status = sfendacc(sds_id)
         END IF
      END DO

      DO i = 1,nscalar
         IF (me == 0) THEN
            ! CREATE AN ARRAY FOR THE CURRENT COMPONENT OF VELOCITY/SCALAR.
            SELECT CASE (i)
            CASE (1)
               u_name = 't'
            CASE (2) 
               u_name = 's'
            CASE (3) 
               u_name = 'r'
            END SELECT
            IF (num_datasets <= 0) THEN
               data = zero
 !              sds_id = sfcreate(sd_id,u_name,DFNT_FLOAT64,vel_rank,vel_length)
            ELSE
               data = zero
 !              sds_id = sfselect(sd_id,3+i)
 !              status = sfrdata(sds_id,vel_start,vel_stride,vel_edges,data)
               IF (status .ne. 0) WRITE(*,*) 'CANNOT READ ',u_name, &
                    ' IN mean.hdf'
            END IF
         END IF

         IF (me == 0) poss_max = zero
         DO plane = 1,ny+1
            CALL MPI_GATHERV( &
                 data_scalar(1,1,plane,i),nx*nz_small_proc(me), &
                 MPI_DOUBLE_PRECISION,   &
                 u2,nx*nz_small_proc,nx*z_small_start,MPI_DOUBLE_PRECISION, &
                 0,MPI_COMM_WORLD,mpierr)
            IF (me == 0) THEN
               data(:,:,plane) = tempa*u2 + tempb*data(:,:,plane)
               poss_max = MAX(MAXVAL(ABS(data(:,:,plane))),poss_max)
            END IF
         END DO
         IF (me == 0) temp(3+i) = poss_max

         ! WRITE THE CURRENT COMPONENT OF SCALAR INTO THE HDF FILE
         IF (me == 0) THEN
!            status = sfwdata(sds_id,vel_start,vel_stride,vel_edges,data)
!            status = sfendacc(sds_id)
         END IF
      END DO

      ! DEALLOCATE THE DATA ARRAY.
      DEALLOCATE(u1,u2)
      IF (me == 0) DEALLOCATE(data)

      IF (num_datasets .le. 0) THEN
         !! ADD y-COORDINATE AS LAST ENTRY IN FILE.
         !! IT HAS SAME RANK AS TIME BUT LENGTH OF ny+1
!         sds_id = sfcreate(sd_id,'ycoord',DFNT_FLOAT64,time_rank,vel_length(2))
!         status = sfwdata(sds_id,time_start,time_stride,vel_length(2), &
!              vel_ycoord)
!         status = sfendacc(sds_id)
      END IF

      ! WRITE MAX VALUES OF VELOCITY/PRESSURE TO STANDARD OUTPUT
      IF (me == 0) THEN
         WRITE(*,998) (temp(i),i=1,3+nscalar)
!!$         WRITE(*,999) dt
998      FORMAT('MAX VALUES =',6E12.4)
999      FORMAT('DT         =',e12.4)
      END IF

      ! CLOSE THE HDF FILE
!      IF (me == 0) status = sfend(sd_id)

    END SUBROUTINE compute_mean
!=========================================================================
    SUBROUTINE compute_uu
      IMPLICIT NONE
      
!      INCLUDE 'hdf.f90'
    
      REAL(KIND=prec), DIMENSION(:,:,:), ALLOCATABLE :: data
      COMPLEX(KIND=prec), DIMENSION(nz,local_nx) :: cu
      REAL(KIND=prec), DIMENSION(:,:), ALLOCATABLE :: u1, u2

      REAL(KIND=prec) :: temp(4+nscalar), poss_max, tempa, tempb

      ! DECLARE HDF FUNCTIONS
      INTEGER   sfstart, sfcreate, sfwdata, sfendacc, sfend, sfselect, sfrdata
      INTEGER   sfdimid, sfsdscale, sfsdmname
      ! DECLARE HDF VARIABLES
      INTEGER, PARAMETER :: vel_rank = 3
      INTEGER, DIMENSION(0:vel_rank-1)::vel_length, vel_start, &
                                        vel_stride, vel_edges
      INTEGER  time_rank, time_length, time_start, time_stride, time_edges
      INTEGER  sd_id, sds_id, status, dim_id
      REAL(KIND=prec) :: vel_xcoord(nx), vel_ycoord(ny+1), vel_zcoord(nz), &
           vel_xcoord2(nx2), vel_zcoord2(nz2)
      INTEGER  num_datasets, time_ref
      CHARACTER(LEN=6) FNAME
      CHARACTER(LEN=4) time_name
      CHARACTER(LEN=1) u_name, dim_name(0:2)

      ALLOCATE(u1(nx,local_nz_small),u2(nx,nz),STAT=ierr)
      IF (me == 0) THEN
         ALLOCATE(data(nx,nz,ny+1),STAT=ierr)
         IF (ierr /= 0) STOP 'CANNOT ALLOCATE ARRAY FOR UU FIELDS'

         ! DESCRIBE VELOCITY FIELD
         vel_length(0) = nx
         vel_length(1) = nz
         vel_length(2) = ny+1
         vel_start(0:2) = 0
         vel_stride(0:2) = 1
         vel_edges(0:2) = vel_length(0:2)

         ! DESCRIBE VELOCITY FIELD
         time_rank = 1
         time_length = 6
         time_start = 0
         time_stride = 1
         time_edges = time_length
         time_name = 'time'
!!$         time_out  = real(t_total)

         ! DESCRIBE COORDINATES
         dim_name(0) = 'x'
         dim_name(1) = 'y'
         dim_name(2) = 'z'
         DO i = 1,nx
            vel_xcoord(i) = real(i-1)*lx/real(nx)
         END DO
         vel_ycoord = real(yh)
         DO i = 1,nz
            vel_zcoord(i) = real(i-1)*lz/real(nz)
         END DO
      END IF

      IF (me == 0) THEN
         fname = 'uu.hdf'
         WRITE(*,*) 'COMPUTING ENSEMBLE-AVERAGED VARIANCE'
!         sd_id = sfstart(FNAME,DFACC_WRITE)
         IF (sd_id .eq. -1) THEN
!            sd_id = sfstart(FNAME,DFACC_CREATE)
!            sds_id = sfcreate(sd_id,time_name,DFNT_FLOAT64, &
!                 time_rank,time_length)
            time(1) = zero
         ELSE
!            sds_id = sfselect(sd_id,0)
!            status = sfrdata(sds_id,time_start,time_stride,time_edges,time)
         END IF
         num_datasets = INT(time(1))
         WRITE(*,*) 'NUMBER OF OLD DATASETS IN uu.hdf = ', num_datasets
         tempa =     one/(time(1) + one) ! WEIGHT OF NEW VELOCITY FIELD
         tempb = time(1)/(time(1) + one) ! WEIGHT OF OLD RUNNING AVERAGE
         time(1) = time(1) + one ! time(1) HOLDS NO. OF FIELDS INCLUDED IN AVG.
!!$         write(*,991) tempa, tempb, time(1)
!!$991      FORMAT(3f8.3)

         ! WRITE OUT TIME STAMP AS FIRST ENTRY IN FILE.
!         status = sfwdata(sds_id,time_start,time_stride,time_edges,time)
!         status = sfendacc(sds_id)
      END IF

      DO i = 1,3
         IF (me == 0) THEN
            ! CREATE AN ARRAY FOR THE CURRENT COMPONENT OF VELOCITY
            SELECT CASE (i)
            CASE (1)
               u_name = 'u'
            CASE (2) 
               u_name = 'v'
            CASE (3) 
               u_name = 'w'
            END SELECT
            IF (num_datasets .eq. 0) THEN
               !! IF uu.hdf DOESN'T EXIST, CREATE DATASET FOR CURRENT FIELD
               data = zero
!               sds_id = sfcreate(sd_id,u_name,DFNT_FLOAT64,vel_rank,vel_length)
            ELSE
               !! IF uu.hdf EXISTS, READ DATASET FOR CURRENT FIELD.
               data = zero
!               sds_id = sfselect(sd_id,i)
!               status = sfrdata(sds_id,vel_start,vel_stride,vel_edges,data)
               IF (status .ne. 0) WRITE(*,*) 'CANNOT READ ', u_name, &
                    ' FIELD IN uu.hdf'
            END IF
         END IF

         IF (me == 0) poss_max = zero
         DO plane = 1,ny+1
            cu = data1(1:nz,1:local_nx,plane,i)
            CALL small_fft_complex_to_real(cu,u1)
            CALL MPI_GATHERV( &
                 u1,nx*nz_small_proc(me),MPI_DOUBLE_PRECISION,   &
                 u2,nx*nz_small_proc(0:nproc-1),nx*z_small_start(0:nproc-1), &
                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
            IF (me == 0) THEN
               data(:,:,plane) = tempa*u2*u2 + tempb*data(:,:,plane)
               poss_max = MAX(MAXVAL(ABS(data(:,:,plane))),poss_max)
            END IF
         END DO
         IF (me == 0) temp(i) = poss_max

         ! WRITE THE CURRENT COMPONENT OF VELOCITY/SCALAR INTO THE HDF FILE
         IF (me == 0) THEN
!            status = sfwdata(sds_id,vel_start,vel_stride,vel_edges,data)
!            status = sfendacc(sds_id)
         END IF
      END DO

      DO i = 1,nscalar
         IF (me == 0) THEN
            ! CREATE AN ARRAY FOR THE CURRENT COMPONENT OF VELOCITY/SCALAR.
            SELECT CASE (i)
            CASE (1)
               u_name = 't'
            CASE (2) 
               u_name = 's'
            CASE (3) 
               u_name = 'r'
            END SELECT
            IF (num_datasets <= 0) THEN
               data = zero
!               sds_id = sfcreate(sd_id,u_name,DFNT_FLOAT64,vel_rank,vel_length)
            ELSE
               data = zero
!               sds_id = sfselect(sd_id,3+i)
!               status = sfrdata(sds_id,vel_start,vel_stride,vel_edges,data)
               IF (status .ne. 0) WRITE(*,*) 'CANNOT READ ',u_name, &
                    ' IN uu.hdf'
            END IF
         END IF

         IF (me == 0) poss_max = zero
         DO plane = 1,ny+1
            CALL MPI_GATHERV( &
                 data_scalar(1,1,plane,i),nx*nz_small_proc(me), &
                 MPI_DOUBLE_PRECISION,   &
                 u2,nx*nz_small_proc,nx*z_small_start, &
                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
            IF (me == 0) THEN
               data(:,:,plane) = tempa*u2*u2 + tempb*data(:,:,plane)
               poss_max = MAX(MAXVAL(ABS(data(:,:,plane))),poss_max)
            END IF
         END DO
         IF (me == 0) temp(3+i) = poss_max

         ! WRITE THE CURRENT COMPONENT OF VELOCITY/SCALAR INTO THE HDF FILE
         IF (me == 0) THEN
!            status = sfwdata(sds_id,vel_start,vel_stride,vel_edges,data)
!            status = sfendacc(sds_id)
         END IF
      END DO

      ! DEALLOCATE THE DATA ARRAY.
      DEALLOCATE(u1,u2)
      IF (me == 0) DEALLOCATE(data)

      IF (num_datasets .le. 0) THEN
         !! ADD y-COORDINATE AS LAST ENTRY IN FILE.
         !! IT HAS SAME RANK AS TIME BUT LENGTH OF ny+1
!         sds_id = sfcreate(sd_id,'ycoord',DFNT_FLOAT64,time_rank,vel_length(2))
!         status = sfwdata(sds_id,time_start,time_stride,vel_length(2), &
!              vel_ycoord)
!         status = sfendacc(sds_id)
      END IF

      ! WRITE MAX VALUES OF VELOCITY/PRESSURE TO STANDARD OUTPUT
      IF (me == 0) THEN
         WRITE(*,998) (temp(i),i=1,3+nscalar)
!!$         WRITE(*,999) dt
998      FORMAT('MAX VALUES =',6E12.4)
999      FORMAT('DT         =',e12.4)
      END IF

      ! CLOSE THE HDF FILE
!      IF (me == 0) status = sfend(sd_id)

    END SUBROUTINE compute_uu
!===================================================================
    SUBROUTINE setup_dims_time
      IMPLICIT NONE

      ! SET UP LENGTH OF OUTPUT ARRAYS
      vel_length(0) = nx2
      vel_length(1) = nz2
      vel_length(2) = ny+1
      ! SET UP COORDINATES.
      DO x = 1,nx2
         vel_xcoord(x) = DBLE(x-1)*lx/DBLE(nx2)
      END DO
      vel_ycoord = yh
      DO z = 1,nz2
         vel_zcoord(z) = DBLE(z-1)*lz/DBLE(nz2)
      END DO

      IF (me == 0) THEN
         ! DESCRIBE VELOCITY FIELD
         vel_name(1) = 'u'
         vel_name(2) = 'v'
         vel_name(3) = 'w'
         vel_name(4) = 't'
         vel_name(5) = 's'
         ! DESCRIBE COORDINATES
         dim_name(0) = 'x'
         dim_name(1) = 'z'
         dim_name(2) = 'y'
         ! DESCRIBE VELOCITY FIELD
         time_rank = 1
         time_length = 6
         time_name = 'time'
         time = zero ! INITIALIZE THE VECTOR HOLDING TIME, SHEAR, PGRAD, RE, PR
         time(1) = t_total
         time(2) = lx
         time(3) = lz
         time(4) = DBLE(bl_grid)
         time(5) = re
         IF (nscalar > 0) time(6) = pr(1)
      END IF

    END SUBROUTINE setup_dims_time
  END SUBROUTINE OUTPUT_SNAPSHOT
!=================================================================


!=================================================================
    subroutine l2_norm
      implicit none
      REAL(KIND=prec) l2, u1(nx2,local_nz), u2(nx2,nz2)
      integer plane,i
!      
      l2=0.D0
!    
      do plane=1,ny+1
         do i=1,3
            CALL xz_fft_complex_to_real(data1(1,1,plane,i),u1)
            call MPI_GATHERV(u1,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
                 u2,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
                 0,MPI_COMM_WORLD,mpierr)
            
            if (me==0) then
               if(i==1)  u2=u2-(1.0D0-y(plane)**2)  ! U steady is not zero
               l2=l2+sum(u2**2)*dy(plane)
            end if
         end do
     end do
      
      if (me==0) then
!     write l2.out and drag.out
         l2=sqrt(l2/dble(nx2)/dble(nz2)/2.0)
         write(16,*) t_total, " , " , l2
         call flush(16)
	 WRITE(11,*) T_TOTAL, " , " ,  D_UM_DY(1)- D_UM_DY(NY)
         call flush(11)
      end if
      
    end subroutine l2_norm


    subroutine boundary
    implicit none
    real(KIND=prec) ul2(nx2,local_nz), ul3(nx2,local_nz), &
         u2(nx2,local_nz),u3(nx2,local_nz),ubot(nx2,local_nz),utop(nx2,local_nz) 
!      
    ubot=0.0D0; utop=0.0D0;
    cbc_top=0.0D0; cbc_bot=0.0D0;

    if (ctrlu.ne.0.0) then
!! U at bottom    
       call xz_fft_complex_to_real(data1(1:nz,1:local_nx,2,1),ul2)
!       call MPI_GATHERV(ul2,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
!            u2,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
!            0,MPI_COMM_WORLD,mpierr)
       call xz_fft_complex_to_real(data1(1:nz,1:local_nx,3,1),ul3)
!       call MPI_GATHERV(ul3,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
!            u3,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
!            0,MPI_COMM_WORLD,mpierr)
       ubot=(dm(2,1)*ul2+dm(3,1)*ul3-2.0D0)/(ctrlu-dm(1,1))
       call xz_fft_real_to_complex(ubot,cbc_bot(1:nz,1:local_nx,1))
!       write(*,*) shape(cbc_bot),nz,local_nx, nx_proc(me)
!       CALL MPI_BCAST(cbc_bot(1:nz,1:nx_proc(me),1),nz*nx_proc(me),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
!stop

! U at top    
       call xz_fft_complex_to_real(data1(1:nz,1:local_nx,ny,1),ul2)
!       call MPI_GATHERV(ul2,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
!            u2,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
!            0,MPI_COMM_WORLD,mpierr)
       call xz_fft_complex_to_real(data1(1:nz,1:local_nx,ny-1,1),ul3)
!       call MPI_GATHERV(ul3,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
!            u3,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
!            0,MPI_COMM_WORLD,mpierr)
       utop=-(dm(2,2)*ul2+dm(1,2)*ul3+2.0D0)/(ctrlu+dm(3,2))
       call xz_fft_real_to_complex(utop,cbc_top(1:nz,1:local_nx,1))
!       write(*,*) "me in boundary=", me, maxval(utop)
!       CALL MPI_BCAST(cbc_top(1:nz,1:local_nx,1),nz*local_nx,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
    end if

    if (ctrlw.ne.0.0) then
!! w at bottom    
       call xz_fft_complex_to_real(data1(1:nz,1:local_nx,2,3),ul2)
!       call MPI_GATHERV(ul2,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
!            u2,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
!            0,MPI_COMM_WORLD,mpierr)
       call xz_fft_complex_to_real(data1(1:nz,1:local_nx,3,3),ul3)
!       call MPI_GATHERV(ul3,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
!            u3,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
!            0,MPI_COMM_WORLD,mpierr)
       ubot=(dm(2,1)*ul2+dm(3,1)*ul3)/(ctrlw-dm(1,1))
       call xz_fft_real_to_complex(ubot,cbc_bot(1:nz,1:local_nx,3))
!       CALL MPI_BCAST(cbc_bot(1:nz,1:local_nx,3),nz*local_nx,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)

! w at top
       call xz_fft_complex_to_real(data1(1:nz,1:local_nx,ny,3),ul2)
!       call MPI_GATHERV(ul2,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
!            u2,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
!            0,MPI_COMM_WORLD,mpierr)
       call xz_fft_complex_to_real(data1(1:nz,1:local_nx,ny-1,3),ul3)
!       call MPI_GATHERV(ul3,nx2*nz_proc(me),MPI_DOUBLE_PRECISION,   &
!            u3,nx2*nz_proc,nx2*z_start, MPI_DOUBLE_PRECISION, &
!            0,MPI_COMM_WORLD,mpierr)
       utop=-(dm(2,2)*ul2+dm(1,2)*ul3)/(ctrlw+dm(3,2))
       call xz_fft_real_to_complex(utop,cbc_top(1:nz,1:local_nx,3))
!       CALL MPI_BCAST(cbc_top(1:nz,1:local_nx,3),nz*local_nx,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpierr)
    end if

  end subroutine boundary

END MODULE read_write

