MODULE gridparms
  USE runparms
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: lx, lz,   delta1,   delta3, divg, ctrlu, ctrlw,             &
            ikx, ikx_me, ikz, kx, kx_me, kz, k2_me, xcoord, zcoord,     &
            y,  yh,   hn, hh,   dn, dh, d2n,    d2h,  iu, dm, d2m,      &
            sum_y,    d_um_dy,  upwind, xzfour, xzphys, left, right,    &
            filter_alpha, filter_delta, filter_coef, metric, jacobian,  &
            covariant_gij, &
            dyh, dy, &
            nx, ny,   nz, nx2,  nz2,    Klim,   Mlim,                   &
            local_nx, local_x_start, nx_proc, x_start,        &
            local_nz, local_z_start, nz_proc, z_start,          &
            local_nz_small, local_z_small_start, nz_small_proc, z_small_start,&
            total_local_size, padded_local_size, small_local_size, &
            xskip, yskip, zskip, filter_output, &
            cbc_bot, cbc_top, bc_top, bc_bot, &
            target_inflow, buffer, target_top, target_bot, ublasius,    &
            jet_forcing, jet_forcing_incr, buffer_small, &
            scount, scount_pad, scount_small, scount_small_pad,      &
            rcount, rcount_pad, rcount_small, rcount_small_pad,      &
            sdispl, sdispl_pad, sdispl_small, sdispl_small_pad,      &
            rdispl, rdispl_pad, rdispl_small, rdispl_small_pad
  
  REAL(KIND=prec)  lx, lz, delta1, delta3, filter_alpha, ctrlu, ctrlw

! grid parameters and indices
  INTEGER    nx, ny, nz, nx2, nz2, Klim, Mlim, &
             padded_local_size, left, right
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: xzfour, xzphys

  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:)    :: kx,  kz, y, yh, hn, hh,  &
                                            kx_me, dn, dh, dyh, sum_y, dy,  &
                                            d_um_dy, filter_delta, jacobian,&
                                            xcoord, zcoord
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:)  :: k2_me, d2n, d2h, iu, dm, &
                     d2m, metric, filter_coef
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:):: upwind, covariant_gij, &
       jet_forcing, jet_forcing_incr

  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:)   :: ikx, ikz, ikx_me
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:) :: divg

! These variables describe the partition of data arrays across the processors.
  INTEGER  local_nx,local_x_start, local_nz, local_z_start, total_local_size, &
       local_nz_small, local_z_small_start, small_local_size
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: &
       nx_proc, nz_proc, x_start, z_start, nz_small_proc, z_small_start, &
       scount, scount_pad, scount_small, scount_small_pad, &
       rcount, rcount_pad, rcount_small, rcount_small_pad, &
       sdispl, sdispl_pad, sdispl_small, sdispl_small_pad, &
       rdispl, rdispl_pad, rdispl_small, rdispl_small_pad

! THESE ARRAYS ALLOW THE SPECIFICATION OF INHOMOGENEOUS BOUNDARY CONDITIONS
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:)    :: bc_top, bc_bot
  COMPLEX(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:) :: cbc_top, cbc_bot

! FRINGE REGION ARRAYS
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:)   :: target_inflow, buffer, &
       buffer_small
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:,:,:) :: target_top, target_bot

! BLASIUS VELOCITY PROFILE
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:)     :: ublasius

! OUTPUT PARAMETERS
  INTEGER  xskip, yskip, zskip, filter_output

END MODULE gridparms





