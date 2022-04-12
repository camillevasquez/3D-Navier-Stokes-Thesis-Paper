MODULE setup
  USE runparms
  USE gridparms
  USE data
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_mpi, read_parms, &
       distribute_data, allocate_data_arrays, &
       finalize_mpi, deallocate_data_arrays, &
       define_parameters, define_coordinates
            
  INCLUDE  'mpif.h'

! miscellaneous loop indices.  
  INTEGER           i, x, z, plane, mpierr, ierr


CONTAINS
!=================================================================
!=======================INITIALIZE MPI============================
!=================================================================
  SUBROUTINE initialize_mpi
    IMPLICIT NONE

    CALL MPI_INIT(mpierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,mpierr)
    CALL mpi_COMM_RANK(MPI_COMM_WORLD,me,mpierr)
    IF (me == 0) THEN
      write(*,*) 'I am processor ', me+1, 'out of ', nproc
    END IF

  END SUBROUTINE initialize_mpi
!=================================================================
!===========================FINALIZE MPI==========================
!=================================================================
  SUBROUTINE finalize_mpi
    IMPLICIT NONE

    CALL mpi_finalize(mpierr)

  END SUBROUTINE finalize_mpi
!=================================================================
!===========================READ_PARMS============================
!=================================================================
  SUBROUTINE read_parms
    IMPLICIT NONE

    INTEGER            scale_by_pi
    REAL(KIND=prec) :: aspect_ratio

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    pi = ATAN2(zero,-one)

! Read information from parameters file, flow.dat    
    OPEN(unit=1, file='flow.dat', form='FORMATTED')
    READ(1,*) re
    DO i = 1,3
       READ(1,*)
    END DO
    read(1,*) ctrlu
    read(1,*) ctrlw
    READ(1,*) lx
    READ(1,*) aspect_ratio
    READ(1,*) scale_by_pi
    IF (scale_by_pi == 1) THEN
       lx = lx*pi
    END IF
    lz = lx/aspect_ratio

    READ(1,*) nx
    READ(1,*) ny
    READ(1,*) nz
    READ(1,*) nx2
    READ(1,*) nz2
    READ(1,*) nscalar
    ALLOCATE(pr(nscalar),pr_t(nscalar),STAT=ierr)
    READ(1,*) (pr(i),i=1,nscalar)
    READ(1,*) (pr_t(i),i=1,nscalar)
    READ(1,*) split_scalar
    READ(1,*) subgrid
    READ(1,*) filter_dim
    READ(1,*) lagr_avg
    READ(1,*) dynamic_scalar
    READ(1,*) tend
    READ(1,*) output_step
    READ(1,*) save_step
    READ(1,*) l2_step
    READ(1,*) read_vel
    READ(1,*) save_vel
    READ(1,*) vel_out
    READ(1,*)
    READ(1,*) snap_out
    READ(1,*)
    READ(1,*) which_test
    READ(1,*) disturb_ampl    
    DO i = 1,10
       READ(1,*)
    END DO
    READ(1,*) test_run
    READ(1,*) fix_mass
    READ(1,*) dt_flag
    READ(1,*) divg_check
    READ(1,*) beta
    beta = beta*asin(1.d0)
    DO i = 1,2
       READ(1,*)
    END DO

! READ IN OUTPUT PARAMETERS
    READ(1,*)
    READ(1,*) xskip
    READ(1,*) yskip
    READ(1,*) zskip
    READ(1,*) filter_output
    
! READ IN JET IN CROSSFLOW PARAMETERS
    READ(1,*)
    READ(1,*) bl_grid
    READ(1,*) ly
    READ(1,*) ly_half
    READ(1,*) delta_bl
    READ(1,*) nbl
    READ(1,*)
    READ(1,*) crossflow
    READ(1,*) filter_jet_profile
    READ(1,*) 
    READ(1,*) jet_diam
    READ(1,*) forcing_phase
    READ(1,*)
    READ(1,*)
    READ(1,*) xjet_top
    READ(1,*) zjet_top
    READ(1,*) jetmag_top
    READ(1,*) forcing_top
    READ(1,*) for_freq_top
    READ(1,*) swirl_top
    READ(1,*)
    READ(1,*) xjet_bot
    READ(1,*) zjet_bot
    READ(1,*) jetmag_bot
    READ(1,*) forcing_bot
    READ(1,*) for_freq_bot
    READ(1,*) swirl_bot

! READ IN BUFFER REGION PARAMETERS
    READ(1,*)
    READ(1,*) fringe_gain
    READ(1,*) 
    READ(1,*) xfringe
    READ(1,*) frac_xbuff
    READ(1,*) buff_xrise
    READ(1,*) buff_xfall
    READ(1,*)
    READ(1,*) zfringe
    READ(1,*) frac_zbuff
    READ(1,*) buff_zrise
    READ(1,*) buff_zfall

    CLOSE(unit=1)

!  IF THIS IS A TEST RUN, SET PARAMETERS TO APPROPRIATE VALUES
    IF (read_vel == 0) THEN
       IF (which_test==2) THEN
          re = 7500.d0
          lx = 2.d0*pi
          beta = 0.d0
          fix_mass = 0
          crossflow = 0
          bl_grid = 0
          xfringe = 0
          zfringe = 0
       ELSEIF (which_test==5) THEN
          re = 7500.d0
          lz = 2.d0*pi
          beta = pi/2.d0
          fix_mass = 0
          crossflow = 0
          bl_grid = 0
          xfringe = 0
          zfringe = 0
       ELSEIF (which_test==8) THEN
          re = 7500.d0
          lx = 2.d0*sqrt(2.d0)*pi
          lz = lx
          beta = pi/4.d0
          fix_mass = 0
          crossflow = 0
          bl_grid = 0
          xfringe = 0
          zfringe = 0
       ELSEIF (which_test == 6) THEN
          re = 10.d0
          lx = 2.d0
          fix_mass = 0
          crossflow = 0
          bl_grid = 0
          xfringe = 0
          zfringe = 0
       END IF
    END IF

! NO DE-ALIASING IF USING THE LAGRANGIAN-AVERAGED DYNAMIC MODEL.
    IF ((subgrid == 1) .and. (lagr_avg == 1)) THEN
       nx2 = nx
       nz2 = nz
    END IF

    IF (me == 0) THEN
      WRITE(*,900) re
      WRITE(*,*)
      DO i = 1,nscalar
         WRITE(*,907) i, split_scalar
         WRITE(*,908) pr(i), pr_t(i)
         WRITE(*,*)
      END DO
      WRITE(*,*) 'BOX DIMENSIONS AND GRID SIZE'
      WRITE(*,901) nx, nx2, lx
      WRITE(*,902) ny
      WRITE(*,903) nz, nz2, lz
      WRITE(*,*)
      WRITE(*,*) 'RUN PARAMETERS'
      WRITE(*,904) tend, read_vel, nscalar, fix_mass
      WRITE(*,905) output_step, save_vel, which_test, dt_flag
      WRITE(*,906) save_step, vel_out, test_run
      IF (vel_out == 1) THEN
         WRITE(*,909) xskip, yskip, zskip, filter_output
      END IF
      WRITE(*,*)

      IF (crossflow == 1) THEN
         WRITE(*,910)          
         WRITE(*,911) bl_grid, jet_diam, filter_jet_profile
         WRITE(*,*)
         IF (jetmag_top /= 0.d0) THEN
            WRITE(*,912) 
            WRITE(*,914) xjet_top, zjet_top
            WRITE(*,915) jetmag_top, swirl_top
            WRITE(*,916) forcing_top, for_freq_top
            WRITE(*,*)
         END IF
         IF (jetmag_bot /= 0.d0) THEN
            WRITE(*,913) 
            WRITE(*,914) xjet_bot, zjet_bot
            WRITE(*,915) jetmag_bot, swirl_bot
            WRITE(*,916) forcing_bot, for_freq_bot
            WRITE(*,*)
         END IF
         IF (jetmag_top*jetmag_bot /= 0.d0) THEN
            WRITE(*,917) forcing_phase
            WRITE(*,*)
         END IF
      END IF
      IF (xfringe == 1) THEN
         WRITE(*,920) 
         WRITE(*,922) fringe_gain, frac_xbuff
         WRITE(*,923) buff_xrise, buff_xfall
         WRITE(*,*)
      END IF
      IF (zfringe == 1) THEN
         WRITE(*,921) 
         WRITE(*,922) fringe_gain, frac_zbuff
         WRITE(*,923) buff_zrise, buff_zfall
         WRITE(*,*)
      END IF
   END IF

900 format('Reynolds number = ',f10.2)
901 format('Nx = ',i4,'   Nx2 = ',i4,'   Lx = ',f6.2)
902 format('Ny = ',i4)
903 format('Nz = ',i4,'   Nz2 = ',i4,'   Lz = ',f6.2)
904 format('tend        = ',i8,'   read vel = ',i2,&
         '   nscalar    = ',i2,'   fix mass = ',i2)
905 format('output step = ',i8,'   save vel = ',i2,&
         '   which test = ',i2,'   dt flag  = ',i2)
906 format('save step   = ',i8,'   vel out  = ',i2,'   test run   = ',i2)
907 format('Scalar  no = ',i2,'          dimensional splitting    = ',i2)
908 format('Prandtl no = ',f10.4,     '  Turbulent Prandtl no = ',f10.4)
909 format('OUTPUT PARAMETERS  skip = (',3i4,') filter_output = ',i2)
910 format('Jet in crossflow simulation')
911 format('Boundary layer domain = ',i2,'  jet diameter = ',f6.2,&
         '  Filter jet profile = ',i2)
912 format('Top wall crossflow jet')
913 format('Bottom wall crossflow jet')
914 format('x location        = ',f6.2,'  z location      = ',f6.2)
915 format('velocity ratio    = ',f6.2,'  swirl number    = ',f6.2)
916 format('Forcing amplitude = ',f6.2,'  Strouhal number = ',f6.2)
917 format('relative phase of jet forcing = ',f6.2)
920 format('Fringe region parameters in x-direction')
921 format('Fringe region parameters in z-direction')
922 format('damping coef = ',f6.2,'  fraction of box length = ',f6.2)
923 format('rise length  = ',f6.2,'  fall length            = ',f6.2)

  END SUBROUTINE read_parms
!=================================================================
!==============DISTRIBUTE DATA AMONG PROCESSORS===================
!=================================================================
  SUBROUTINE distribute_data
    IMPLICIT NONE

    INTEGER  xleft, zleft, z_small_left, n

! ALLOCATE GRID AMONG THE PROCESSORS.
! GRID IS SPLIT ALONG Z-DIRECTION IN PHYSICAL SPACE,
! ALONG X-DIRECTION IN FOURIER SPACE.
! See how many vectors will be left over after being divided evenly.
    ALLOCATE(nx_proc(0:nproc-1),nz_proc(0:nproc-1),nz_small_proc(0:nproc-1), &
             x_start(0:nproc-1),z_start(0:nproc-1),z_small_start(0:nproc-1), &
             scount(0:nproc-1),       scount_pad(0:nproc-1), &
             scount_small(0:nproc-1), scount_small_pad(0:nproc-1), &
             rcount(0:nproc-1),       rcount_pad(0:nproc-1), &
             rcount_small(0:nproc-1), rcount_small_pad(0:nproc-1), &
             sdispl(0:nproc-1),       sdispl_pad(0:nproc-1), &
             sdispl_small(0:nproc-1), sdispl_small_pad(0:nproc-1), &
             rdispl(0:nproc-1),       rdispl_pad(0:nproc-1), &
             rdispl_small(0:nproc-1), rdispl_small_pad(0:nproc-1), &
             STAT=ierr)
    xleft = mod(nx/2+1,nproc)
    zleft = mod(nz2,nproc)
    z_small_left = mod(nz,nproc)
    IF (me==0) WRITE(*,*) xleft, zleft, z_small_left
    DO n = 0,nproc-1
      IF (n < xleft) THEN
        nx_proc(n) = (nx/2+1-xleft)/nproc + 1
      ELSE
        nx_proc(n) = (nx/2+1-xleft)/nproc
      END IF
      IF (n < zleft) THEN
        nz_proc(n) = (nz2-zleft)/nproc + 1
      ELSE
        nz_proc(n) = (nz2-zleft)/nproc
      END IF
      IF (n < z_small_left) THEN
        nz_small_proc(n) = (nz-z_small_left)/nproc + 1
      ELSE
        nz_small_proc(n) = (nz-z_small_left)/nproc
      END IF
    END DO

    x_start(0) = 0
    z_start(0) = 0
    z_small_start(0) = 0
    DO n = 1,nproc-1
      x_start(n) = x_start(n-1) + nx_proc(n-1)
      z_start(n) = z_start(n-1) + nz_proc(n-1)
      z_small_start(n) = z_small_start(n-1) + nz_small_proc(n-1)
    END DO

! Compute maximum block size.  Use this blocksize for all messages 
! to permit use of the (more likely optimized) MPI_Alltoall command.
    scount = nz_proc(me)*nx_proc      !SCOUNT FOR P2F TRANSFORM, RCOUNT FOR F2P
    rcount = nz_proc    *nx_proc(me)  !RCOUNT FOR P2F TRANSFORM, SCOUNT FOR F2P

    scount_pad = (nz_proc(me)+2)*nx_proc
    rcount_pad = (nz_proc    +2)*nx_proc(me)

    scount_small = nz_small_proc(me)*nx_proc
    rcount_small = nz_small_proc    *nx_proc(me)

    scount_small_pad = (nz_small_proc(me)+2)*nx_proc
    rcount_small_pad = (nz_small_proc    +2)*nx_proc(me)

    sdispl = 0      !SDISPL FOR P2F TRANSFORM, RDISPL FOR F2P
    rdispl = 0      !RDISPL FOR P2F TRANSFORM, SDISPL FOR F2P

    sdispl_pad = 0
    rdispl_pad = 0

    sdispl_small = 0
    rdispl_small = 0

    sdispl_small_pad = 0
    rdispl_small_pad = 0

    DO n = 1,nproc-1
       sdispl(n) = sdispl(n-1) + scount(n-1)
       rdispl(n) = rdispl(n-1) + rcount(n-1)

       sdispl_pad(n) = sdispl_pad(n-1) + scount_pad(n-1)
       rdispl_pad(n) = rdispl_pad(n-1) + rcount_pad(n-1)

       sdispl_small(n) = sdispl_small(n-1) + scount_small(n-1)
       rdispl_small(n) = rdispl_small(n-1) + rcount_small(n-1)

       sdispl_small_pad(n) = sdispl_small_pad(n-1) + scount_small_pad(n-1)
       rdispl_small_pad(n) = rdispl_small_pad(n-1) + rcount_small_pad(n-1)
    END DO

! Assign local values so that data arrays may be declared.
    local_nx      = nx_proc(me)
    local_x_start = x_start(me)
    local_nz      = nz_proc(me)
    local_z_start = z_start(me)
    local_nz_small      = nz_small_proc(me)
    local_z_small_start = z_small_start(me)
    total_local_size = MAX(SUM(scount),SUM(rcount))
    padded_local_size = max(SUM(scount_pad),SUM(rcount_pad), &
         (2+maxval(nz_proc))*(nx2/2+1))
    small_local_size  = max(SUM(scount_small_pad),SUM(rcount_small_pad), &
         (2+maxval(nz_small_proc))*(nx/2+1))

! Write out partitioning of data to standard output.
    IF (me == 0) THEN
       WRITE(*,921)
       DO n = 0,nproc-1
          WRITE(*,922) n, nx_proc(n), x_start(n), nz_proc(n), z_start(n), &
                          nz_small_proc(n), z_small_start(n)
       END DO
       WRITE(*,*)    
       WRITE(*,924)
       DO n = 0,nproc-1
          WRITE(*,926) n, scount(n), sdispl(n), scount_pad(n), sdispl_pad(n)
       END DO
       WRITE(*,*)    
       WRITE(*,925)
       DO n = 0,nproc-1
          WRITE(*,926) n, scount_small(n), sdispl_small(n), &
               scount_small_pad(n), sdispl_small_pad(n)
       END DO
       WRITE(*,*)    
    END IF
921 FORMAT('Proc   nx_proc   x_start   nz_proc    z_start   nz_small  z_small')
922 FORMAT(i4,i7,3x,i7,3x,i7,3x,i7,3x,i8,3x,i7)
923 FORMAT(4i8)
924 FORMAT('Proc   scountc   sdispl  scount_pad  sdispl_pad')
925 FORMAT('Proc   scountc   sdispl  scount_pad  sdispl_pad')
926 FORMAT(i4,i7,3x,i7,3x,i7,3x,i7)

! SET ARRAY TO INDICATE WHICH PROCESSOR HOLDS A PARTICULAR x-z PLANE.
    left = me - 1
    right = me + 1
    IF (left < 0) left = nproc-1
    IF (right > nproc-1) right = 0
!!$    write(*,*) me, left, right
!!$    WRITE(*,*) 
    
  END SUBROUTINE distribute_data
!=================================================================
!==================ALLOCATE DATA ARRAYS===========================
!=================================================================
  SUBROUTINE allocate_data_arrays
    IMPLICIT NONE
    REAL(KIND=prec), PARAMETER :: zero = 0.d0
    INTEGER recl_vel

! Allocate data arrays accordingly.
    ALLOCATE(data1(nz,local_nx,ny+1,3), &
             data2(ny+1,3,nz,local_nx),  &
             data4(ny+1,nz,local_nx),            &
             divg(local_nx,nz), STAT=ierr)
    IF (ierr /= 0) THEN
       WRITE(*,*) 'NOT ENOUGH SPACE TO ALLOCATE DATA ARRAYS ON n',me+1
       STOP
    ELSE
      data1 = zero
      data2 = zero
      data4 = zero
      divg = zero
    END IF

    IF (nscalar > 0) THEN
       ALLOCATE(data_scalar(nx,local_nz_small,ny+1,1:nscalar), STAT=ierr)
       IF (ierr /= 0) THEN
          WRITE(*,*) 'NOT ENOUGH SPACE TO ALLOCATE SCALAR ARRAYS ON n',me+1
          STOP
       ELSE
          data_scalar = zero
       END IF
    END IF

    IF (me ==0) THEN
       recl_vel = (ny+1)*17
       OPEN (UNIT=11,FILE='drag.out',FORM='FORMATTED')
       OPEN (UNIT=16,FILE='l2.out',FORM='FORMATTED',position='rewind')
       OPEN (UNIT=12,FILE='vel.out', FORM='FORMATTED',RECL=recl_vel)
    END IF

! Allocate arrays to hold dynamic coefficients.
    IF (subgrid == 1) THEN

! ALWAYS DECLARE ARRAY FOR PLANE-AVERAGED EDDY VISCOSITY.
! IT CAN BE USED TO INITIALIZE LAGRANGIAN EDDY VISCOSITY.
       IF (dynamic_scalar == 1) THEN
          ALLOCATE(dyn_visc(1:ny,1:1+nscalar),STAT=ierr)
          ndynscalar = nscalar
       ELSE
          ALLOCATE(dyn_visc(1:ny,1:1),STAT=ierr)
          ndynscalar = 0
       END IF
       IF (ierr/=0) WRITE(*,*) 'FAILED TO DECLARE ARRAY FOR DYN VISCOSITY'
       dyn_visc = zero

! IF LAGRANGIAN AVERAGING OF EDDY VISCOSITY WILL BE USED, DECLARE ARRAYS
! FOR AVERAGED DYNAMIC MODEL COEFFICIENTS.
       IF (lagr_avg==1) THEN
          IF (dynamic_scalar == 1) THEN
             ALLOCATE(data_dyn(0:nx+1,0:local_nz_small+1,1:ny,1:2+2*nscalar), &
                  STAT=ierr)
          ELSE
             ALLOCATE(data_dyn(0:nx+1,0:local_nz_small+1,1:ny,1:2),STAT=ierr)
          END IF
          IF (ierr/=0) WRITE(*,*) 'FAILED TO DECLARE ARRAY FOR DYN VISCOSITY'
          data_dyn = zero
       ELSEIF (lagr_avg==0) THEN
! OPEN FILE TO STORE PLANE-AVERAGED EDDY VISCOSITY TIME HISTORY.         
          IF (me==0) THEN
             OPEN(unit=10,file='dyn_visc.out',RECL=recl_vel,FORM='FORMATTED')
          END IF
       END IF
    END IF

    IF (me == 0) WRITE(*,*) 'DATA ARRAYS ALLOCATED'

  END SUBROUTINE allocate_data_arrays
!=================================================================
!==========================DEALLOCATE DATA ARRAYS=================
!=================================================================
  SUBROUTINE deallocate_data_arrays
    IMPLICIT NONE

    DEALLOCATE(data1,data2,data4,divg)
    IF (nscalar > 0) DEALLOCATE(data_scalar)
    DEALLOCATE(nx_proc,nz_proc,x_start,z_start)
    DEALLOCATE(kx_me,kx,kz,ikx_me,ikx,ikz,k2_me)
    DEALLOCATE(y,yh,dn,dh,dm,hn,hh,d2h,d2n,d_um_dy,sum_y,dyh,dy,iu,upwind)
! De-allocate arrays for dynamic coefficients.
    IF ((subgrid == 1).and.(lagr_avg==1)) DEALLOCATE(data_dyn)
    IF (linearized_run == 1) DEALLOCATE(baseflow)

    IF (me==0) THEN
       CLOSE (UNIT=11)
       CLOSE (UNIT=12)
       CLOSE (UNIT=16)
       IF ((subgrid==1).and.(lagr_avg==0)) CLOSE (UNIT=10)
    END IF

  END SUBROUTINE deallocate_data_arrays
!=================================================================
!=======================DEFINE COORDINATES========================
!=================================================================
  SUBROUTINE DEFINE_COORDINATES
    IMPLICIT NONE

    REAL(KIND=prec)  STRETCH, STRETCH2, TEMP1, TEMP2, TEMP3, TEMP4, XI
    PARAMETER (STRETCH = 2.75D0)
    PARAMETER (STRETCH2 = 1.4D0)

    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0, &
                                  three=3.0D0

    ALLOCATE(y(1:ny),      hn(1:ny),          sum_y(1:ny),          &
             yh(1:ny+1),   hh(1:ny+1),        dyh(1:ny+1),          &
             xcoord(nx+1), zcoord(nz+1),      dy(1:ny+1),           &
             dn(1:ny),     d2h(1:3,1:ny+1),   iu(1:2,1:ny),         &
             dh(1:ny+1),   d2n(1:3,1:ny),     upwind(1:3,1:ny,1:2), &
             dm(1:3,1:4),  d2m(1:4,1:2),      d_um_dy(1:ny), metric(1:3,1:ny),&
             jacobian(ny), filter_delta(ny),  filter_coef(3,2),     &
             covariant_gij(1:3,1:3,1:ny), STAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(*,*) 'NOT ENOUGH SPACE TO ALLOCATE GRID ARRAYS ON PROCESSOR ',&
                  me+1
      STOP
    END IF

    ALLOCATE(cbc_bot(nz,local_nx,3+nscalar),bc_bot(nx2,local_nz,3+nscalar), &
         cbc_top(nz,local_nx,3+nscalar),bc_top(nx2,local_nz,3+nscalar), &
         buffer(nx2,local_nz), buffer_small(nx,local_nz_small), &
         target_inflow(ny+1,3+nscalar), &
         target_top(nx2,ny+1,2), target_bot(nx2,ny+1,2),ublasius(ny+1), &
         STAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(*,*) 'NOT ENOUGH SPACE TO ALLOCATE BC ARRAYS ON PROCESSOR ',&
                  me+1
      STOP
    END IF

!  INITIALIZE BC ARRAYS.
    cbc_bot = zero
    cbc_top = zero
    bc_bot = zero
    bc_top = zero
    buffer = zero
    buffer_small = zero
    target_inflow = zero
    target_top = zero
    target_bot = zero

!  SET UP THE GRID SPACING IN THE X AND Z DIRECTIONS
    DELTA1 = LX/DBLE(NX2)
    DELTA3 = LZ/DBLE(NZ2)
    DO x = 1,nx+1
       xcoord(x) = dble(x-1)*lx/dble(nx)
    END DO
    DO z = 1,nz+1
       zcoord(z) = dble(z-1)*lz/dble(nz)
    END DO

!  COMPUTE THE Y-COORDINATE
    IF ((crossflow == 0) .and. (bl_grid == 0)) THEN

       DO PLANE = 1,NY
          XI = TWO*DBLE(PLANE-1)/DBLE(NY-1) - ONE
          Y(PLANE) = TANH(STRETCH*XI)/TANH(STRETCH)
          !     Y(PLANE) = TANH(STRETCH2*XI + (STRETCH2**3/THREE)*XI**3)
          !     .         /TANH(STRETCH2 + STRETCH2**3/THREE)
       END DO

       IF ((read_vel == 0) .and. (which_test == 6)) THEN
          y(1) = -half
          DO plane = 2,ny-1
             y(plane) = y(1) + dble(plane-1)/dble(ny-1)
          END DO
          y(ny) = half
       END IF

       ! MAKE CHANNEL FLOW GRID PERFECTLY SYMMETRIC ABOUT CENTERLINE.
       !    do j = ny,ny/2+1,-1
       !       y(j) = - y(ny+1-j)
       !    end do

       YH(1) = Y(1)
       yh(2:ny) = half*(y(1:ny-1) + y(2:ny))
       YH(NY+1) = Y(NY)

    ELSEIF (bl_grid == 1) THEN
       CALL blasius_grid(y,yh,ublasius)
    ELSEIF (bl_grid == 2) THEN
       CALL blasius_grid2(y,yh,ublasius)
    ELSEIF ((bl_grid == 3) .or. (crossflow == 1)) THEN
       CALL blasius_grid3(y,yh,ublasius)
    END IF

    IF (me == 0) CALL output_grids

!  DEFINE SPACING AND FIRST DIVERATIVE OPERATORS, 
!  DH TO COMPUTE DERIVATIVE AT HALF-POINT LOCATIONS (OF V)
!  DN TO COMPUTE DERIVATIVE AT Y_N (OF U, W AND P)
    HH(1)    = Y(2)    - Y(1)
    HH(2:NY) = Y(2:NY) - Y(1:NY-1)
    HH(NY+1) = Y(NY)   - Y(NY-1)    
    DH(1:NY+1) = ONE/HH(1:NY+1)

    HN(1:NY) = YH(2:NY+1) - YH(1:NY)
    DN(1:NY) = ONE/HN(1:NY)

!  DEFINE INTERPOLATION OPERATOR FOR U AND W, TO BE USED IN
!  NONLINEAR TERMS FOR V.
    IU(1,1:NY) = HALF*HH(1:NY)/HN(1:NY)
    IU(2,1:NY) = HALF*HH(2:NY+1)/HN(1:NY)

!  DEFINE A DIFFERENCE OPERATOR FOR EVALUATING DU/DY AT THE
!     WALL.  THESE PERMIT GOOD APPROXIMATION OF THE SECOND
!     DERIVATIVE OF U AT Y_1/2 AND Y_N-1/2.
!!$      DM(1,1) = - TWO*TWO/(HN(1)+HN(2))
!!$      DM(2,1) =   THREE*DN(2)
!!$      DM(3,1) =   DN(2)*(HN(2) - THREE*HN(1))/(HN(1)+HN(2))
!!$
!!$      DM(1,2) = - DN(NY-1)*(HN(NY-1) - THREE*HN(NY))/(HN(NY)+HN(NY-1))
!!$      DM(2,2) = - THREE*DN(NY-1)
!!$      DM(3,2) =   TWO*TWO/(HN(NY)+HN(NY-1))

    DM(1,1) = - DN(1)*(TWO*HN(1) + HN(2))/(HN(1) + HN(2))
    DM(2,1) =   DN(1)*DN(2)*(HN(1)+HN(2))
    DM(3,1) = - DN(2)*HN(1)/(HN(1)+HN(2))

    DM(1,2) =   DN(NY-1)*HN(NY)/(HN(NY-1)+HN(NY))
    DM(2,2) = - DN(NY)*DN(NY-1)*(HN(NY)+HN(NY-1))
    DM(3,2) =   DN(NY)*(TWO*HN(NY)+HN(NY-1))/(HN(NY)+HN(NY-1))

! DEFINE A SECOND-ORDER ONE-SIDED DIFFERENCE OPERATOR
! FOR USE IN APPLYING NEUMANN BOUNDARY CONDITIONS TO THE VERTICAL VELOCITY
! THIS WILL GUARANTEE A SECOND-ORDER METHOD.
    DM(1,3) = (y(1)-y(2) + y(1)-y(3))/((y(1)-y(2))*(y(1)-y(3)))
    DM(2,3) = (y(1)-y(3))/((y(2)-y(1))*(y(2)-y(3)))
    DM(3,3) = (y(1)-y(2))/((y(3)-y(1))*(y(3)-y(2)))

    DM(1,4) = (y(ny)-y(ny-1))/((y(ny-2)-y(ny))*(y(ny-2)-y(ny-1)))
    DM(2,4) = (y(ny)-y(ny-2))/((y(ny-1)-y(ny))*(y(ny-1)-y(ny-2)))
    DM(3,4) = (y(ny)-y(ny-1)+y(ny)-y(ny-2))/((y(ny)-y(ny-1))*(y(ny)-y(ny-2)))

!  DEFINE A SECOND DIFFERENCE OPERATOR FOR EVALUATING DU/DY AT THE
!     WALL.  THESE PERMIT GOOD APPROXIMATION OF THE SECOND
!     DERIVATIVE OF U AT Y_1/2 AND Y_N-1/2.
    temp1 = yh(1) - yh(2)
    temp2 = yh(1) - yh(3)
    temp3 = yh(1) - yh(4)
    D2M(1,1) = two*(temp1+temp2+temp3)/(temp1*temp2*temp3)
    D2M(2,1) = two*(temp2+temp3)/((yh(2)-yh(1))*(yh(2)-yh(3))*(yh(2)-yh(4)))
    D2M(3,1) = two*(temp1+temp3)/((yh(3)-yh(1))*(yh(3)-yh(2))*(yh(3)-yh(4)))
    D2M(4,1) = two*(temp1+temp2)/((yh(4)-yh(1))*(yh(4)-yh(2))*(yh(4)-yh(3)))

    temp1 = yh(ny+1)-yh(ny)
    temp2 = yh(ny+1)-yh(ny-1)
    temp3 = yh(ny+1)-yh(ny-2)
    D2M(1,2) = two*(temp1+temp2) &
         /((yh(ny-2)-yh(ny+1))*(yh(ny-2)-yh(ny))*(yh(ny-2)-yh(ny-1)))
    D2M(2,2) = two*(temp1+temp3) &
         /((yh(ny-1)-yh(ny+1))*(yh(ny-1)-yh(ny))*(yh(ny-1)-yh(ny-2)))
    D2M(3,2) = two*(temp2+temp3) &
         /((yh(ny)-yh(ny+1))*(yh(ny)-yh(ny-1))*(yh(ny)-yh(ny-2)))
    D2M(4,2) = two*(temp1+temp2+temp3)/(temp1*temp2*temp3)

!  SAME FOR SECOND DERIVATIVE OPERATORS.
    DO PLANE = 2,NY-1
       D2N(1,PLANE) =   DN(PLANE)*DH(PLANE)
       D2N(2,PLANE) = - DN(PLANE)*DH(PLANE) - DN(PLANE)*DH(PLANE+1)
       D2N(3,PLANE) =                         DN(PLANE)*DH(PLANE+1)

       D2H(1,PLANE) =   DH(PLANE)*DN(PLANE-1)
       D2H(2,PLANE) = - DH(PLANE)*DN(PLANE-1) - DH(PLANE)*DN(PLANE)
       D2H(3,PLANE) =                           DH(PLANE)*DN(PLANE)
    END DO
    D2H(1,2) =               - DH(2)*DM(1,1)
    D2H(2,2) = - DH(2)*DN(2) - DH(2)*DM(2,1) 
    D2H(3,2) =   DH(2)*DN(2) - DH(2)*DM(3,1)

    D2H(1,NY) =   DH(NY)*DN(NY-1) + DH(NY)*DM(1,2)
    D2H(2,NY) = - DH(NY)*DN(NY-1) + DH(NY)*DM(2,2)
    D2H(3,NY) =                     DH(NY)*DM(3,2)
    DO I = 1,3
       D2H(I,1) = D2H(I,2)
       D2N(I,1) = D2N(I,2)

       D2N(I,NY)   = D2N(I,NY-1)
       D2H(I,NY+1) = D2H(I,NY)
    END DO
    
!  DEFINE THE COEFFICIENTS FOR INTEGRATION IN Y, USING THE
!     HALF-POINTS IN THE INTERIOR OF THE DOMAIN.
!     THIS IS USED IN INTEGRATING THE RHS OF THE MEAN
!     VELOCITY EQUATION TO PRESERVE CONSTANT MASS FLUX.
!!$    TEMPH = TWO/DBLE(NY-1)
!!$    DO J = 1,NY+1
!!$       DYH(J) = TEMPH*(STRETCH/DTANH(STRETCH))* &
!!$           (ONE - (DTANH(STRETCH)*YH(J))**2)
!!$    END DO
!!$    DYH(1) = ZERO
!!$    DYH(2) = DYH(2)*17.D0/24.D0
!!$    DYH(3) = DYH(3)*33.D0/24.D0
!!$    DYH(4) = DYH(4)*22.D0/24.D0
!!$    DYH(NY-2) = DYH(NY-2)*22.D0/24.D0
!!$    DYH(NY-1) = DYH(NY-1)*33.D0/24.D0
!!$    DYH(NY)   = DYH(NY)  *17.D0/24.D0
!!$    DYH(NY+1) = ZERO

! USE MIDPOINT METHOD FOR INTEGRATION.
    dyh(1)    = zero
    dyh(2:ny) = hh(2:ny)
    dyh(ny+1) = zero

    dy(1)      = half*hh(2)
    dy(2:ny-1) = half*(hh(2:ny-1)+hh(3:ny))
    dy(ny)     = half*hh(ny)
    dy(ny+1)   = zero

!  DEFINE INTERPOLATION OPERATORS FOR UPWINDING OF THE CONSERVED SCALAR.  
!  UPWIND(I,PLANE,1) WILL BE USED IF THE V(PLANE) IS POSITIVE, 
!     UPWIND(I,PLANE,2) OTHERRWISE.

! INITIALIZE ARRAY.
    upwind = zero

!  DEFINE THE INTERPOLATION OPERATOR FOR THE CELL FACE AT Y(PLANE)
!     IF THE VELOCITY THERE, V(PLANE), IS POSITIVE.
    UPWIND(1,1,1) = ONE

!  DEFINE THE INTERPOLATION OPERATOR FOR THE CELL FACE AT Y(PLANE)
!     IF THE VELOCITY THERE, V(PLANE), IS NEGATIVE.
    TEMP1 = Y(1) - YH(1)
    TEMP2 = Y(1) - YH(2)
    TEMP3 = Y(1) - YH(3)
    TEMP4 = Y(1) - YH(4)

    UPWIND(1,1,2) = TEMP2*TEMP3/((TEMP2-TEMP1)*(TEMP3-TEMP1))
    UPWIND(2,1,2) = TEMP1*TEMP3/((TEMP1-TEMP2)*(TEMP3-TEMP2))
    UPWIND(3,1,2) = TEMP1*TEMP2/((TEMP1-TEMP3)*(TEMP2-TEMP3))

!  REPEAT FOR THE REST OF THE PLANES IN THE BOX.
    DO PLANE = 2,NY-1
       TEMP1 = Y(PLANE) - YH(PLANE-1)
       TEMP2 = Y(PLANE) - YH(PLANE)
       TEMP3 = Y(PLANE) - YH(PLANE+1)
       TEMP4 = Y(PLANE) - YH(PLANE+2)
       UPWIND(1,PLANE,1) = TEMP2*TEMP3/((TEMP2-TEMP1)*(TEMP3-TEMP1))
       UPWIND(2,PLANE,1) = TEMP1*TEMP3/((TEMP1-TEMP2)*(TEMP3-TEMP2))
       UPWIND(3,PLANE,1) = TEMP1*TEMP2/((TEMP1-TEMP3)*(TEMP2-TEMP3))
       UPWIND(1,PLANE,2) = TEMP3*TEMP4/((TEMP3-TEMP2)*(TEMP4-TEMP2))
       UPWIND(2,PLANE,2) = TEMP2*TEMP4/((TEMP2-TEMP3)*(TEMP4-TEMP3))
       UPWIND(3,PLANE,2) = TEMP2*TEMP3/((TEMP2-TEMP4)*(TEMP3-TEMP4))
    END DO
! JUST BELOW TOP WALL.
    TEMP1 = Y(NY) - YH(NY-2)
    TEMP2 = Y(NY) - YH(NY-1)
    TEMP3 = Y(NY) - YH(NY)
    TEMP4 = Y(NY) - YH(NY+1)
    UPWIND(1,NY,1) = TEMP3*TEMP4/((TEMP3-TEMP2)*(TEMP4-TEMP2))
    UPWIND(2,NY,1) = TEMP2*TEMP4/((TEMP2-TEMP3)*(TEMP4-TEMP3))
    UPWIND(3,NY,1) = TEMP2*TEMP3/((TEMP2-TEMP4)*(TEMP3-TEMP4))
    UPWIND(3,NY,2) = ONE

!  DEFINE THE GRID STRETCHING METRICS -- USEFUL IN LES WHERE FILTERING
!  CAN BE PERFORMED IN COMPUTATIONAL SPACE, RATHER THAN PHYSICAL SPACE.

!  Following the suggestion of jordan, a large-eddy simulation method
!     in generalized curvilinear coordinates, jcp 148:322--340, 1998,
!     filtering is performed on a uniform grid in transformed coordinates.
!     The grid spacing is uniform in the transformed coordinates
!     with spacing one in all directions.  Filtering is then performed on
!     the contravariant contravariant quantities, e.g. in constructing L12
!     the following quantity is needed:
!
!           filter( jacobian*metric(2,.)*u(.,1)*u(.,2))

    metric(1,1:ny) = dble(nx)/lx
    metric(3,1:ny) = dble(nz)/lz

    metric(2,1)      = one/(y(2)-y(1))
    metric(2,2:ny-1) = two/(Y(3:ny)-Y(1:ny-2))
    metric(2,ny)     = one/(y(ny)-y(ny-1))

    jacobian(1:ny)   = (metric(1,1:ny)*metric(2,1:ny)*metric(3,1:ny))**(-1)
    filter_delta(1:ny) = jacobian(1:ny)**(1.d0/3.d0)

    covariant_gij(1:3,1:3,1:ny) = zero

    covariant_gij(1,1,1:ny) = metric(1,1:ny)**(-2)
    covariant_gij(2,2,1:ny) = metric(2,1:ny)**(-2)
    covariant_gij(3,3,1:ny) = metric(3,1:ny)**(-2)

!!$    DO plane = 1,ny
!!$       WRITE(*,999) ((covariant_gij(x,z,plane),x=1,3),z=1,3)
!!$999    FORMAT(9e10.2)
!!$    END DO
!!$    STOP 'IN SETUP'

    CALL COMPUTE_FILTER_COEF(FILTER_COEF(1,1),FILTER_COEF(1,2), &
         filter_alpha,2.0D0,ONE,ONE)

  CONTAINS
!==============================================================
    FUNCTION DELTA_ANISO(DX,DY,DZ)
      IMPLICIT NONE

! USE IF FILTERING QUANTITIES IN PHYSICAL SPACE.  ADAPTS VALUE OF
! DELTA FOR ANISOTROPIC GRIDS FOLLOWING SCOTTI, MENEVEAU & LILLY, 
! PHYS FLU (199?)
      REAL(KIND=prec)  DELTA_ANISO, DX, DY, DZ, A1, A2
      
      A1 = DMIN1(DX,DY,DZ)/DMAX1(DX,DY,DZ)
      A2 = DMAX1(DMIN1(DX,DY),DMIN1(DX,DZ),DMIN1(DY,DZ))/DMAX1(DX,DY,DZ)
      DELTA_ANISO = (DX*DY*DZ)**(1.D0/3.D0)*DCOSH(DSQRT((4.D0/27.D0) &
           *(DLOG(A1)**2 - DLOG(A1)*DLOG(A2) + DLOG(A2)**2)))

    END FUNCTION DELTA_ANISO
!==============================================================
    SUBROUTINE COMPUTE_FILTER_COEF(COEFA,COEFB,ALPHA,DELTA,H0,H1)
      IMPLICIT NONE
          
      REAL(KIND=prec) :: COEFA(3), COEFB(3), ALPHA, DELTA, H0, H1
      REAL(KIND=prec), PARAMETER :: zero=0.0, half=0.5, one=1.0, two=2.0
      
! THE COMPUTATION OF THE FILTER WIDTH alpha IS DERIVED FROM A SUGGESTION
! OF TOM LUND IN THE CTR ANNUAL REPORT:
!@techreport{Lund97,
!author="Lund, T.~S.",
!title="On the use of discrete filters for large eddy simulation",
!institution="Center for Turbulence Research",
!series="Annual Research Briefs", pages="83--95", year=1997}
! AVAILABLE ONLINE AT: http://ctr-sgi1.stanford.edu/CTR/ResBriefs97/lund.ps.gz

! FILTER COEFFICIENTS COMPUTED SO THAT THIS BOX FILTER HAS A 
! ZEROTH MOMENT OF ONE, A ZERO FIRST MOMENT AND A SECOND MOMENT 
! OF DELTA**2/8.  THIS YIELDS TRAPEZOIDAL RULE WHEN H0 = H1 = DELTA/2.
      COEFA(1) = DELTA**2/(8.D0*H0*(H0+H1))
      COEFA(3) = DELTA**2/(8.D0*H1*(H0+H1))
      COEFA(2) = ONE - COEFA(1) - COEFA(3)

!  THIS SET OF COEFFICIENTS WILL GIVE SIMPSON'S RULE.
!!$        COEFA(1) = DELTA**2/(12.D0*H0*(H0+H1)) 
!!$        COEFA(3) = DELTA**2/(12.D0*H1*(H0+H1))
!!$        COEFA(2) = ONE - COEFA(1) - COEFA(3)

      ALPHA = DSQRT(12.D0*(COEFA(1) + COEFA(3)))

      COEFB(1) = (ALPHA**2 - ONE)*H1/(12.d0*(H0+H1))
      COEFB(3) = (ALPHA**2 - ONE)*H0/(12.d0*(H0+H1))
      COEFB(2) = ONE - COEFB(1) - COEFB(3)

    END SUBROUTINE COMPUTE_FILTER_COEF
!============================================
    SUBROUTINE blasius_grid(yblas,yhblas,ublas)
      IMPLICIT NONE

      REAL(KIND=prec) ::  yblas(ny), yhblas(ny+1), ublas(ny+1)
      INTEGER           I, J, PLANE
      REAL(KIND=prec)  :: F(3), DF(3), COEFA, COEFB, COEFC, UTMP(NY+1), H

! SET INITIAL CONDITION FOR SECOND DERIVATIVE IN BLASIUS PROFILE INTEGRATION
! AND DEFINITION OF DISPLACEMENT THICKNESS IN TERMS OF ETA, THE BLASIUS
! SIMILARITY VARIABLE.
      REAL(KIND=prec), PARAMETER :: &
           D2F = 0.33205733621291D0, &
           DISPL_THICKNESS = 1.72D0

!  SET UP THE PARAMETERS FOR THE RUNGE-KUTTA TIME INTEGRATION.
!     WE ARE USING A FIVE-STEP, FOURTH ORDER RUNGE-KUTTA SCHEME
!     FOR THE EXPLICIT TERMS (NONLINEAR AND ANY EDDY VISOSITY TERMS).
!     THE (CONSTANT COEFFICIENT) VISCOUS TERMS ARE COMPUTED
!     IMPLICITLY WITH A CRANK-NICOLSON SCHEME WITHIN EACH RK SUBSTEP.
!     FOR DETAILS ON THE RK SCHEME, LOOK FOR A PAPER:
!          MH CARPENTER AND CA KENNEDY, A FOURTH-ORDER 2N-STORAGE
!                RUNGE-KUTTA SCHEME.  (UNPUBLISHED)
!   Coefficients taken from Wilson, Demuren & Carpenter,
!     High-order compact schemes for numerical simulation of
!     incompressible flows, ICASE Report 98-13, 1998.  

      REAL(KIND=prec), DIMENSION(5), PARAMETER :: &
           rka = (/ 0.0D0, -0.41789047D0, -1.19215169D0, &
                    -1.69778469D0, -1.51418344D0 /),   &
           rkb = (/ 0.14965902D0, 0.37921031D0, 0.82295502D0, &
                    0.69945045D0, 0.15305724D0 /)
 
      REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

      COEFA = (LY - TWO*DELTA_BL)/DBLE(NY-NBL)
      COEFB = DELTA_BL/(DISPL_THICKNESS*DBLE(NBL))
      COEFC = DISPL_THICKNESS/DELTA_BL      

      IF (BL_GRID .EQ. 0) THEN

!  SET UP INITIAL CONDITION FOR INTEGRATION OF BLASIUS PROFILE.
         f = zero
         df = zero
         f(3) = d2f

!  SET INITIAL VALUE OF Y AND INTEGRATE TO FIND OTHER Y.  THEY'RE
!     ROUGHLY EVENLY SPACE IN THE BOUNDARY LAYER ACCORDING TO
!     VELOCITY INCREMENTS, I.E. UBLAS(J) - UBLAS(J-1) ~ CONSTANT.
!     OUTSIDE THE BL, THEY'RE EVENLY SPACED WITH A SMOOTH TRANSITION
!     IN BETWEEN USING A HYBERBOLIC TANGENT.

         yblas(1) = zero
         DO plane = 1,ny/2

!  COMPUTE THE Y VALUE OF THE NEXT POINT, USING A SIMPLE INTEGRATION RULE.

            IF (coefc*yblas(plane) .lt. 20.d0) THEN
               yblas(plane+1) = yblas(plane) + coefa*TANH(coefb/(coefa*f(3)))
            ELSE
               YBLAS(PLANE+1) = YBLAS(PLANE) + COEFA            
            END IF

!  INTEGRATE THE BLASIUS PROFILE UP TO YBLAS(PLANE+1)
!  USE 100 STEPS OF RK4 SCHEME BETWEEN EACH PAIR OF
!     GRID POINTS IN THE WALL-NORMAL DIRECTION.

            H = COEFC*(YBLAS(PLANE+1) - YBLAS(PLANE))/100.d0
            DO J = 1,100
               DO RK_STEP = 1,5
                  DF(1) = RKA(RK_STEP)*DF(1) + F(2)
                  DF(2) = RKA(RK_STEP)*DF(2) + F(3)
                  DF(3) = RKA(RK_STEP)*DF(3) - HALF*F(1)*F(3)

                  F(1) = F(1) + H*RKB(RK_STEP)*DF(1)
                  F(2) = F(2) + H*RKB(RK_STEP)*DF(2)
                  F(3) = F(3) + H*RKB(RK_STEP)*DF(3)
               END DO
            END DO
         END DO
         yblas(ny/2+1) = two*yblas(ny/2) - yblas(ny/2-1)
         DO plane = ny/2+2,ny
            yblas(plane) = yblas(plane-1) &
                 + (yblas(ny-plane+2) - yblas(ny-plane+1))
         END DO

      ELSE

!  SET UP INITIAL CONDITION FOR INTEGRATION OF BLASIUS PROFILE.
         DO I = 1,3
            F(I) = ZERO
            DF(I) = ZERO
         END DO
         F(3) = D2F

!  SET INITIAL VALUE OF Y AND INTEGRATE TO FIND OTHER Y.  THEY'RE
!     ROUGHLY EVENLY SPACE IN THE BOUNDARY LAYER ACCORDING TO
!     VELOCITY INCREMENTS, I.E. UBLAS(J) - UBLAS(J-1) ~ CONSTANT.
!     OUTSIDE THE BL, THEY'RE EVENLY SPACED WITH A SMOOTH TRANSITION
!     IN BETWEEN USING A HYBERBOLIC TANGENT.

         YBLAS(1) = ZERO
         DO PLANE = 1,NY-1

!  COMPUTE THE Y VALUE OF THE NEXT POINT, USING A SIMPLE INTEGRATION RULE.

            IF (COEFC*YBLAS(PLANE) .LT. 20.D0) THEN
               YBLAS(PLANE+1) = YBLAS(PLANE) + COEFA*DTANH(COEFB/(COEFA*F(3)))
            ELSE
               YBLAS(PLANE+1) = YBLAS(PLANE) + COEFA            
            END IF

!  INTEGRATE THE BLASIUS PROFILE UP TO YBLAS(PLANE+1)
!  USE 100 STEPS OF RK4 SCHEME BETWEEN EACH PAIR OF
!     GRID POINTS IN THE WALL-NORMAL DIRECTION.

            H = COEFC*(YBLAS(PLANE+1) - YBLAS(PLANE))/100.d0
            DO J = 1,100
               DO RK_STEP = 1,5
                  DF(1) = RKA(RK_STEP)*DF(1) + F(2)
                  DF(2) = RKA(RK_STEP)*DF(2) + F(3)
                  DF(3) = RKA(RK_STEP)*DF(3) - HALF*F(1)*F(3)

                  F(1) = F(1) + H*RKB(RK_STEP)*DF(1)
                  F(2) = F(2) + H*RKB(RK_STEP)*DF(2)
                  F(3) = F(3) + H*RKB(RK_STEP)*DF(3)
               END DO
            END DO

         END DO

      END IF

! NORMALIZE THE Y-COORDINATE SO THAT yblas(NY) == ly
      yblas = ly*yblas/MAXVAL(yblas)

!  NOW, WITH THE VALUES OF YBLAS(PLANE) -- THE V POINTS -- SPECIFIED,
!     COMPUTE THE VALUES OF YHBLAS(PLANE) -- THE U,W,P POINTS.
      YHBLAS(1) = YBLAS(1)
      DO PLANE = 2,NY
         YHBLAS(PLANE) = HALF*(YBLAS(PLANE) + YBLAS(PLANE-1))
      END DO
      YHBLAS(NY+1) = YBLAS(NY)

!  SET UP INITIAL CONDITION FOR INTEGRATION OF BLASIUS PROFILE.
      DO I = 1,3
         F(I) = ZERO
         DF(I) = ZERO
      END DO
      F(3) = D2F

!  COMPUTE THE BLASIUS PROFILE ON THE STAGGERED POINTS.

      UBLAS(1) = F(2)
      DO PLANE = 1,NY
      
!  INTEGRATE THE BLASIUS PROFILE UP TO YBLAS(PLANE+1)
!  USE 100 STEPS OF RK4 SCHEME BETWEEN EACH PAIR OF
!     GRID POINTS IN THE WALL-NORMAL DIRECTION.

         H = COEFC*(YHBLAS(PLANE+1) - YHBLAS(PLANE))/100.d0
         DO J = 1,100
            DO RK_STEP = 1,5
               DF(1) = RKA(RK_STEP)*DF(1) + F(2)
               DF(2) = RKA(RK_STEP)*DF(2) + F(3)
               DF(3) = RKA(RK_STEP)*DF(3) - HALF*F(1)*F(3)

               F(1) = F(1) + H*RKB(RK_STEP)*DF(1)
               F(2) = F(2) + H*RKB(RK_STEP)*DF(2)
               F(3) = F(3) + H*RKB(RK_STEP)*DF(3)
            END DO
         END DO

         UBLAS(PLANE+1) = F(2)

      END DO

      IF (BL_GRID .EQ. 0) THEN
         DO PLANE = 1,NY+1
            UTMP(PLANE) = UBLAS(PLANE)
         END DO
         DO PLANE = 1,NY+1
            UBLAS(PLANE) = UTMP(PLANE)*UTMP(NY+2-PLANE)
         END DO
      END IF

    END SUBROUTINE BLASIUS_GRID
!============================================
    SUBROUTINE blasius_grid2(yblas,yhblas,ublas)
      IMPLICIT NONE

      REAL(KIND=prec) :: yblas(ny), yhblas(ny+1), ublas(ny+1), tmp_dy(ny-1)
      REAL(KIND=prec) :: F(3), DF(3), UTMP(NY+1)
      REAL(KIND=prec) :: COEF1, COEF2, coef3, coefc, tmp_h
      INTEGER           I, J, PLANE, n0

! SET INITIAL CONDITION FOR SECOND DERIVATIVE IN BLASIUS PROFILE INTEGRATION
! AND DEFINITION OF DISPLACEMENT THICKNESS IN TERMS OF ETA, THE BLASIUS
! SIMILARITY VARIABLE.
      REAL(KIND=prec), PARAMETER :: &
           D2F = 0.33205733621291D0, &
           DISPL_THICKNESS = 1.72D0

!  SET UP THE PARAMETERS FOR THE RUNGE-KUTTA TIME INTEGRATION.
!     WE ARE USING A FIVE-STEP, FOURTH ORDER RUNGE-KUTTA SCHEME
!     FOR THE EXPLICIT TERMS (NONLINEAR AND ANY EDDY VISOSITY TERMS).
!     THE (CONSTANT COEFFICIENT) VISCOUS TERMS ARE COMPUTED
!     IMPLICITLY WITH A CRANK-NICOLSON SCHEME WITHIN EACH RK SUBSTEP.
!     FOR DETAILS ON THE RK SCHEME, LOOK FOR A PAPER:
!          MH CARPENTER AND CA KENNEDY, A FOURTH-ORDER 2N-STORAGE
!                RUNGE-KUTTA SCHEME.  (UNPUBLISHED)
!   Coefficients taken from Wilson, Demuren & Carpenter,
!     High-order compact schemes for numerical simulation of
!     incompressible flows, ICASE Report 98-13, 1998.  

      REAL(KIND=prec), DIMENSION(5), PARAMETER :: &
           rka = (/ 0.0D0, -0.41789047D0, -1.19215169D0, &
                    -1.69778469D0, -1.51418344D0 /),   &
           rkb = (/ 0.14965902D0, 0.37921031D0, 0.82295502D0, &
                    0.69945045D0, 0.15305724D0 /)
 
      REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

      IF (bl_grid == 0) THEN

! GRID IS DENSE AT UPPER AND BOTTOM BOUNDARIES.
! SET COEFFICIENTS FOR GRID STRETCHING.
! COEF1 IS THE (UNIFORM) GRID SPACING NEAR THE WALL
! COEF2 IS THE APPROXIMATE ASYMPTOTIC SPACING AT INFINITY
! COEF3 IS THE INVERSE WIDTH OF THE TRANSITION REGION BETWEEN THE TWO SPACINGS
         COEF1 = delta_bl/DBLE(nbl)
         COEF2 = (half*LY - 6.D0*DELTA_BL)/DBLE(NY/2-NBL)
         coef3 = 0.64d0
         n0 = 3*nbl
         DO i = 1,ny/2
            tmp_dy(i) = coef1 + half*coef2*(one &
                 + TANH(coef3*(DBLE(i)/DBLE(n0) - DBLE(n0)/DBLE(i))))
            tmp_dy(ny-i) = tmp_dy(i)
         END DO
         tmp_dy = tmp_dy*ly/SUM(tmp_dy)

! DEFINE THE COORDINATE LOCATIONS.
         yblas(1)  = zero
         yhblas(1) = zero
         DO i = 1,ny-1
            yblas(i+1)  = yblas(i) + tmp_dy(i)
            yhblas(i+1) = yblas(i) + half*tmp_dy(i)
         END DO
         yh(ny+1) = y(ny)

      ELSE

! SET COEFFICIENTS FOR GRID STRETCHING.
! COEF1 IS THE (UNIFORM) GRID SPACING NEAR THE WALL
! COEF2 IS THE APPROXIMATE ASYMPTOTIC SPACING AT INFINITY
! COEF3 IS THE INVERSE WIDTH OF THE TRANSITION REGION BETWEEN THE TWO SPACINGS
         COEF1 = delta_bl/DBLE(nbl)
         COEF2 = (LY - 6.D0*DELTA_BL)/DBLE(NY-NBL)
         coef3 = 0.62d0
         n0 = 3*nbl
         DO i = 1,ny-1
            tmp_dy(i) = coef1 + half*coef2*(one &
                 + TANH(coef3*(DBLE(i)/DBLE(n0) - DBLE(n0)/DBLE(i))))
         END DO
         tmp_dy = tmp_dy*ly/SUM(tmp_dy)

! DEFINE THE COORDINATE LOCATIONS.
         yblas(1)  = zero
         yhblas(1) = zero
         DO i = 1,ny-1
            yblas(i+1)  = yblas(i) + tmp_dy(i)
            yhblas(i+1) = yblas(i) + half*tmp_dy(i)
         END DO
         yh(ny+1) = y(ny)

      END IF

      COEFC = DISPL_THICKNESS/DELTA_BL      

!  SET UP INITIAL CONDITION FOR INTEGRATION OF BLASIUS PROFILE.
      f  = zero
      df = zero
      f(3) = d2f

!  COMPUTE THE BLASIUS PROFILE ON THE STAGGERED POINTS.
      ublas(1) = f(2)
      DO plane = 1,ny
!  INTEGRATE THE BLASIUS PROFILE UP TO YBLAS(PLANE+1)
!  USE 100 STEPS OF RK4 SCHEME BETWEEN EACH PAIR OF
!     GRID POINTS IN THE WALL-NORMAL DIRECTION.
         tmp_h = 0.01d0*coefc*(yhblas(plane+1) - yhblas(plane))
         DO j = 1,100
            DO rk_step = 1,5
               df(1) = rka(rk_step)*df(1) + f(2)
               df(2) = rka(rk_step)*df(2) + f(3)
               df(3) = rka(rk_step)*df(3) - half*f(1)*f(3)
               f = f + tmp_h*rkb(rk_step)*df
            END DO
         END DO
         ublas(plane+1) = f(2)
      END DO

      IF (BL_GRID .EQ. 0) THEN
         DO PLANE = 1,NY+1
            UTMP(PLANE) = UBLAS(PLANE)
         END DO
         DO PLANE = 1,NY+1
            UBLAS(PLANE) = UTMP(PLANE)*UTMP(NY+2-PLANE)
         END DO
      END IF

    END SUBROUTINE blasius_grid2
!============================================
    SUBROUTINE blasius_grid3(yblas,yhblas,ublas)
      IMPLICIT NONE

      REAL(KIND=prec) :: yblas(ny), yhblas(ny+1), ublas(ny+1), tmp_dy(ny-1)
      REAL(KIND=prec) :: F(3), DF(3), UTMP(NY+1)
      REAL(KIND=prec) :: COEF1, COEF2, coef3, coefc, tmp_h, h, a, b
      INTEGER           I, J, PLANE, n0

! SET INITIAL CONDITION FOR SECOND DERIVATIVE IN BLASIUS PROFILE INTEGRATION
! AND DEFINITION OF DISPLACEMENT THICKNESS IN TERMS OF ETA, THE BLASIUS
! SIMILARITY VARIABLE.
      REAL(KIND=prec), PARAMETER :: &
           D2F = 0.33205733621291D0, &
           DISPL_THICKNESS = 1.72D0

!  SET UP THE PARAMETERS FOR THE RUNGE-KUTTA TIME INTEGRATION.
!     WE ARE USING A FIVE-STEP, FOURTH ORDER RUNGE-KUTTA SCHEME
!     FOR THE EXPLICIT TERMS (NONLINEAR AND ANY EDDY VISOSITY TERMS).
!     THE (CONSTANT COEFFICIENT) VISCOUS TERMS ARE COMPUTED
!     IMPLICITLY WITH A CRANK-NICOLSON SCHEME WITHIN EACH RK SUBSTEP.
!     FOR DETAILS ON THE RK SCHEME, LOOK FOR A PAPER:
!          MH CARPENTER AND CA KENNEDY, A FOURTH-ORDER 2N-STORAGE
!                RUNGE-KUTTA SCHEME.  (UNPUBLISHED)
!   Coefficients taken from Wilson, Demuren & Carpenter,
!     High-order compact schemes for numerical simulation of
!     incompressible flows, ICASE Report 98-13, 1998.  

      REAL(KIND=prec), DIMENSION(5), PARAMETER :: &
           rka = (/ 0.0D0, -0.41789047D0, -1.19215169D0, &
                    -1.69778469D0, -1.51418344D0 /),   &
           rkb = (/ 0.14965902D0, 0.37921031D0, 0.82295502D0, &
                    0.69945045D0, 0.15305724D0 /)
 
      REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

      IF (bl_grid == 0) THEN

! GRID IS DENSE AT UPPER AND BOTTOM BOUNDARIES.
! SET COEFFICIENTS FOR GRID STRETCHING.
! COEF1 IS THE (UNIFORM) GRID SPACING NEAR THE WALL
! COEF2 IS THE APPROXIMATE ASYMPTOTIC SPACING AT INFINITY
! COEF3 IS THE INVERSE WIDTH OF THE TRANSITION REGION BETWEEN THE TWO SPACINGS
         COEF1 = delta_bl/DBLE(nbl)
         COEF2 = (half*LY - 6.D0*DELTA_BL)/DBLE(NY/2-NBL)
         coef3 = 0.64d0
         n0 = 3*nbl
         DO i = 1,ny/2
            tmp_dy(i) = coef1 + half*coef2*(one &
                 + TANH(coef3*(DBLE(i)/DBLE(n0) - DBLE(n0)/DBLE(i))))
            tmp_dy(ny-i) = tmp_dy(i)
         END DO
         tmp_dy = tmp_dy*ly/SUM(tmp_dy)

! DEFINE THE COORDINATE LOCATIONS.
         yblas(1)  = zero
         yhblas(1) = zero
         DO i = 1,ny-1
            yblas(i+1)  = yblas(i) + tmp_dy(i)
            yhblas(i+1) = yblas(i) + half*tmp_dy(i)
         END DO
         yh(ny+1) = y(ny)

      ELSE

         ! DEFINE GRID.
         ! GRID IS STRETCHED USING MAPPING: y = A*(1+eta)/(B-eta)
         ! WHERE eta \in [-1,1] AND y \in [0,L].  HERE
         ! A = L*h/(L-2*h) AND B = L/(L-2*h).  THE POINT eta = 0
         ! IS MAPPED TO y = h, SO THAT HALF OF THE GRID POINTS ARE
         ! WITHIN THE INTERVAL [0,h].
         IF (ly_half .ne. 0.5d0*ly) THEN
            a = ly*ly_half/(ly-two*ly_half)
            b = ly/(ly-two*ly_half)
            DO i = 1,ny
               yblas(i) = -one + two*DBLE(i-1)/DBLE(ny-1)
            END DO
            yblas = a*(1 + yblas)/(b - yblas)
            yblas(1)  = zero
         ELSE
            DO i = 1,ny
               yblas(i) = zero + ly*DBLE(i-1)/DBLE(ny-1)
            END DO
         END IF
         yhblas(1)    = yblas(1)
         yhblas(2:ny) = half*(yblas(2:ny)+yblas(1:ny-1))
         yhblas(ny+1) = yblas(ny)

      END IF

      COEFC = DISPL_THICKNESS/DELTA_BL      

!  SET UP INITIAL CONDITION FOR INTEGRATION OF BLASIUS PROFILE.
      f  = zero
      df = zero
      f(3) = d2f

!  COMPUTE THE BLASIUS PROFILE ON THE STAGGERED POINTS.
      ublas(1) = f(2)
      DO plane = 1,ny
!  INTEGRATE THE BLASIUS PROFILE UP TO YBLAS(PLANE+1)
!  USE 100 STEPS OF RK4 SCHEME BETWEEN EACH PAIR OF
!     GRID POINTS IN THE WALL-NORMAL DIRECTION.
         tmp_h = 0.01d0*coefc*(yhblas(plane+1) - yhblas(plane))
         DO j = 1,100
            DO rk_step = 1,5
               df(1) = rka(rk_step)*df(1) + f(2)
               df(2) = rka(rk_step)*df(2) + f(3)
               df(3) = rka(rk_step)*df(3) - half*f(1)*f(3)
               f = f + tmp_h*rkb(rk_step)*df
            END DO
         END DO
         ublas(plane+1) = f(2)
      END DO

      IF (BL_GRID .EQ. 0) THEN
         DO PLANE = 1,NY+1
            UTMP(PLANE) = UBLAS(PLANE)
         END DO
         DO PLANE = 1,NY+1
            UBLAS(PLANE) = UTMP(PLANE)*UTMP(NY+2-PLANE)
         END DO
      END IF

    END SUBROUTINE blasius_grid3
  END SUBROUTINE define_coordinates
!=====================================================================
!=========================OUTPUT GRIDS================================
!=====================================================================
  SUBROUTINE output_grids
    IMPLICIT NONE

!  OUTPUT GRID FILE FOR READING INTO MATLAB.
    open (unit=19, file='makegrid.m', form='FORMATTED')
    write(19,*) 'y = ['
    DO plane = 1,ny
       write(19,990) y(plane)
    END DO
    write(19,*) '];'
    write(19,*) 'yh = ['
    do plane = 1,ny+1
       write(19,990) yh(plane)
    end do
    write(19,*) '];'
    write(19,*) 'x = ['
    do x = 1,nx+1
       write(19,990) lx*dble(x-1)/dble(nx)
    end do
    write(19,*) '];'
    write(19,*) 'z = ['
    do z = 1,nz+1
       write(19,990) lz*dble(z-1)/dble(nz)
    end do
    write(19,*) '];'
    close (unit=19)
    990 format(1f16.8)

  END SUBROUTINE output_grids
!=================================================================
!=======================DEFINE PARAMETERS=========================
!=================================================================
  SUBROUTINE define_parameters
    IMPLICIT none
    REAL(KIND=prec) :: tmp
    REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    CALL random_seed ! INITIALIZE RANDOM NUMBER GENERATOR.
    DO i = 1,100*me
       CALL RANDOM_NUMBER(tmp)
    END DO

    pi = ATAN2(zero,-one)

    linearized_run = 0

    dt      =  0.25d0*lx/DBLE(nx)
    cfl_max =  0.95d0
!!$    IF ((read_vel == 0) .and. (which_test == 9)) dt = lx/DBLE(tend)
    IF (dt_flag == 1) THEN
       IF (subgrid == 1) cfl_max = 0.5d0
       IF ((subgrid == 1) .and. (lagr_avg == 1)) cfl_max = 0.3d0
    END IF
    p_grad  = -2.0d0/re
    IF (xfringe == 1) p_grad = zero
    IF ((read_vel == 0) .and. (mod(which_test,3) == 2) .and. (tend >= 40)) THEN
         DT = 2.0d0*(2.0d0*PI/0.24989154D0)/DBLE(TEND)
    END IF
    IF ((read_vel == 0) .and. (which_test == 6)) THEN
       dt = 0.5d0/dble(ny)
       tend = int(1.d0/dt)
!       dt = 1.0d0/dble(tend)
    END IF
    IF (crossflow==1) THEN
       IF (read_vel == 0) THEN
          dt = 0.40d0*1.73d0*(yh(2) - yh(1)) &
               /MAX(jetmag_bot*(one+forcing_bot),jetmag_top*(one+forcing_top))
       ELSE
          dt = 0.6d0*1.73d0*(yh(2) - yh(1)) &
               /MAX(jetmag_bot*(one+forcing_bot),jetmag_top*(one+forcing_top))
       END IF
    END IF
    IF (me == 0) WRITE(*,*) 'INITIAL DT = ', dt

! Allocate arrays for wavenumbers and compute.
    IF (nx == 1) THEN
      ALLOCATE(kx_me(1:1), kz(1:nz/2), ikx_me(1:1), ikz(1:nz/2), &
             kx(1:1), ikx(1:1), k2_me(1:nz/2,1:1), STAT=ierr)
      kx_me(1) = 0.0
      ikx_me(1) = 0.0
      kx(1) = 0.0
      ikx(1) = 0.0
      DO z = 1,nz/2
        kz(z) = DBLE(z-1)*2.0*pi/lz
        ikz(z) = CMPLX(0.,kz(z),KIND=prec)
        k2_me(z,x) = kz(z)*kz(z)
      END DO
      Klim = 1
      Mlim = nz/2

    ELSEIF (nz == 1) THEN
      ALLOCATE(kx_me(1:local_nx),  kz(1:1),  kx(1:nx/2+1), &
              ikx_me(1:local_nx), ikz(1:1), ikx(1:nx/2+1), &
               k2_me(1:1,1:local_nx), STAT=ierr)
      DO x = 1,local_nx
        kx_me(x)  = dble(local_x_start+x-1)*2.d0*pi/lx
        ikx_me(x) = cmplx(0.d0,kx_me(x),KIND=prec)
        k2_me(z,x) = kx_me(x)*kx_me(x)
      END DO
      DO x = 1,nx/2+1
        kx(x)  = dble(x-1)*2.d0*pi/lx
        ikx(x) = cmplx(0.d0,kx(x),KIND=prec)
      END DO
      kz(1) = 0.d0
      ikz(1) = 0.d0
      Klim = local_nx
      Mlim = 1

    ELSE
      ALLOCATE(kx_me(1:local_nx),  kx(1:nx/2+1),  kz(1:nz), &
              ikx_me(1:local_nx), ikx(1:nx/2+1), ikz(1:nz), &
               k2_me(1:nz,1:local_nx), STAT=ierr)
      DO x = 1,local_nx
        kx_me(x)  = dble(local_x_start+x-1)*2.d0*pi/lx
        ikx_me(x) = cmplx(0.d0,kx_me(x),KIND=prec)
      END DO
      DO x = 1,nx/2+1
        kx(x)  = dble(x-1)*2.d0*pi/lx
        ikx(x) = cmplx(0.d0,kx(x),KIND=prec)
      END DO
      DO z = 1,nz/2+1
        kz(z)  = DBLE(z-1)*2.d0*pi/lz
        ikz(z) = CMPLX(0.d0,kz(z),KIND=prec)
      END DO
      DO z = nz,nz/2+2,-1
        kz(z)  = -kz(nz+2-z)
        ikz(z) = CMPLX(0.d0,kz(z),KIND=prec)
      END DO
      DO x = 1,local_nx
         DO z = 1,nz
          k2_me(z,x) = kx_me(x)*kx_me(x) + kz(z)*kz(z)
        END DO
      END DO
      Klim = local_nx
      Mlim = nz
    END IF

! IF CROSSFLOW JET RUN, MULTIPLY JET POSITIONS BY BOX LENGTH
! AND DECLARE VARIABLES FOR JET FORCING AT INPUT (AFTER FREUND 2001, JFM 438).
    IF (crossflow == 1) THEN
       ! SET THE LOCATION OF THE JET IN THE X- AND Z-DIRECTIONS.
       xjet_top = lx*xjet_top
       xjet_bot = lx*xjet_bot
       zjet_top = lz*zjet_top
       zjet_bot = lz*zjet_bot
       
       ! ALLOCATE VARIABLE TO HOLD JET FORCING PARAMETERS AND INCREMENT.
       ALLOCATE(jet_forcing(3,2,4),jet_forcing_incr(3,2,4),STAT=ierr)
       jet_forcing      = zero
       jet_forcing_incr = zero
       IF (me == 0) THEN
          !! INITIALIZE JET FORCING
          CALL RANDOM_NUMBER(jet_forcing)
          jet_forcing(:,:,1) = 0.01d0 + 0.06d0*jet_forcing(:,:,1) ! AMPLITUDE
          jet_forcing(:,:,2) = 0.1d0  + 0.6d0 *jet_forcing(:,:,2) ! STROUHAL #
          jet_forcing(:,:,3) = two*pi*jet_forcing(:,:,3)          ! PHI
          jet_forcing(:,:,4) = two*pi*jet_forcing(:,:,4)          ! PSI

          !! INITIALIZE JET FORCING INCREMENT -- PROVIDES SIGN FOR CHANGE 
          !! OF JET FORCING PARAMETERS AT EACH STEP.
          CALL RANDOM_NUMBER(jet_forcing_incr)
          jet_forcing_incr = SIGN(0.1d0,jet_forcing_incr-half)
       END IF
       CALL MPI_BCAST(jet_forcing,24,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, &
            mpierr)
    END IF

  END SUBROUTINE define_parameters
END MODULE setup    



