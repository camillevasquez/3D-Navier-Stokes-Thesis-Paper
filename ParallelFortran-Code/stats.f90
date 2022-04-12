PROGRAM statistics
  USE runparms
  USE gridparms
  USE setup
  IMPLICIT NONE
  
  INCLUDE 'hdf.f90'

  INTEGER, PARAMETER :: prec = SELECTED_REAL_KIND(14,32)

  INTEGER i, j, k, n, nn, x, z, ierr
  REAL(KIND=prec) wgt_old, wgt_new, scale, dwall(3,2)

  REAL(KIND=prec), DIMENSION(:,:), ALLOCATABLE :: moments, moments_vor, &
       budget_vel, budget_scalar, budget_vor, fluxes, taylor_scale
  REAL(KIND=prec), DIMENSION(:), ALLOCATABLE   :: yplus, temp
  REAL(KIND=prec), DIMENSION(:,:), ALLOCATABLE :: diff, dmean, d2mean, &
       omega, domega

  REAL  time(4), u_tau, re_tau, drag(2), press_grad, time_in
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: u

  REAL(KIND=prec), DIMENSION(:,:), ALLOCATABLE       :: temp_u
  REAL(KIND=prec), DIMENSION(:,:,:), ALLOCATABLE     :: vor
  REAL(KIND=prec), DIMENSION(:,:,:,:), ALLOCATABLE   :: sij, grad_u, grad_vor
  REAL(KIND=prec), DIMENSION(:,:,:,:,:), ALLOCATABLE :: d2_u

  INTEGER dsgdata, dsnum, dsfirst, dsgnt, dsgdims, status, dsrref

  INTEGER, PARAMETER :: vel_rank = 3
  INTEGER, DIMENSION(vel_rank):: vel_length
  INTEGER  time_rank, time_length, num_datasets, fname_rank, fname_length
  INTEGER  nfields
  CHARACTER(LEN=15) fname
  CHARACTER(LEN=15) output_file

  REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! READ IN RUN PARAMETERS AND SET UP COORDINATES.
  CALL read_parms
  test_run = 1
  CALL initialize_fft
  CALL define_parameters
  CALL define_coordinates

! ALLOCATE DATA ARRAY AND STATISTICS ARRAY
  ny = ny + 1
  ALLOCATE(u(1:nx+1,1:ny,1:nz+1,4+nscalar),fluxes(ny,3+nscalar), &
       moments(ny,4*(3+nscalar)),moments_vor(ny,12), &
       budget_vel(ny,5), budget_vor(ny,7),budget_scalar(ny,4*nscalar), &
       taylor_scale(ny,1+nscalar), yplus(1:ny/2), diff(3,ny), STAT=ierr)
  moments = zero
  fluxes = zero
  moments_vor = zero
  budget_vel = zero
  budget_vor = zero
  budget_scalar = zero
  taylor_scale = zero

  ALLOCATE(dmean(ny,3+nscalar),d2mean(ny,3+nscalar),omega(ny,3),domega(ny,3),&
       sij(1:nx,1:nz,1:3,1:3),grad_u(1:nx,1:nz,1:3+nscalar,1:3), &
       grad_vor(nx,nz,3,3),d2_u(nx,nz,3,3,3),vor(nx,nz,3), temp(ny), &
       temp_u(1:nx,1:nz), STAT=ierr)

! SET UP DIFFERENTIATION COEFFICIENT IN WALL-NORMAL DIRECTION.
  CALL setup_diff(3,yh(1),yh(1),diff(1,1))
  DO j = 2,ny-1
     CALL setup_diff(3,yh(j),yh(j-1),diff(1,j))
  END DO
  CALL setup_diff(3,yh(ny),yh(ny-2),diff(1,ny))

!!$  diff(1,1) = one/(y(1) - y(2)) + one/(y(1) - y(3))
!!$  diff(2,1) = (y(1) - y(3))/((y(2)-y(1))*(y(2)-y(3)))
!!$  diff(3,1) = (y(1) - y(2))/((y(3)-y(1))*(y(3)-y(2)))
!!$
!!$  diff(1,2:ny-1) = - half*dh(2:ny-1)
!!$  diff(2,2:ny-1) =   half*dh(2:ny-1) - half*dh(3:ny)
!!$  diff(3,2:ny-1) =                     half*dh(3:ny)
!!$
!!$  diff(1,ny) = (y(ny) - y(ny-1))/((y(ny-2)-y(ny))*(y(ny-2)-y(ny-1)))
!!$  diff(2,ny) = (y(ny) - y(ny-2))/((y(ny-1)-y(ny))*(y(ny-1)-y(ny-2)))
!!$  diff(3,ny) = one/(y(ny) - y(ny-2)) + one/(y(ny) - y(ny-1))

! SET PARAMETERS FOR READING DATAFILE NAMES.
  fname_rank = 1
  fname_length = 12

! DESCRIBE TIME
  time_rank = 1
  time_length = 4

! DESCRIBE VELOCITY
  vel_length(1) = nx+1
  vel_length(2) = ny
  vel_length(3) = nz+1

! FIND THE NUMBER OF DATASETS IN THE ENSEMBLE.
  num_datasets = dsnum('datafiles.hdf')

  DO n = 1,num_datasets

     wgt_old = DBLE(n-1)/DBLE(n)
     wgt_new =       one/DBLE(n)

! READ IN THE FILENAME OF THE nTH DATA SET.
     status = dsrref('datafiles.hdf',n+1)
     status = dsgdata('datafiles.hdf',fname_rank,fname_length,fname)

! READ IN THE TIME
     status = dsgdata(fname,time_rank,time_length,time)
     WRITE(*,*) n, fname, '  TIME = ', time(1)
     drag(1:2) = wgt_old*drag(1:2) + wgt_new*time(2:3)
     press_grad = wgt_old*press_grad + wgt_new*time(4)

! READ IN THE VELOCITY/SCALAR ARRAY
     DO i = 1,3+nscalar
        status = dsgdata(fname,vel_rank,vel_length,u(:,:,:,i))
     END DO

! COMPUTE MEAN VELOCITY AND SCALAR PROFILES.
     scale = one/DBLE(nx*nz)
     nn = 3+nscalar
     DO i = 1,3+nscalar
        DO j = 1,ny
           moments(j,i) = wgt_old*moments(j,i) &
                + wgt_new*scale*SUM(u(1:nx,j,1:nz,i))
        END DO
     END DO

  END DO

! COMPUTE GRADIENT OF MEAN PROFILE ALONG WITH MEAN VORTICITY
! AND THE GRADIENT OF MEAN VORTICITY.
  DO i = 1,nn
     dmean(1,i) = SUM(diff(1:3,1)*moments(1:3,i))
     d2mean(1,i) = SUM(d2n(1:3,1)*moments(1:3,i))
     DO j = 2,ny-1
        dmean(j,i) = SUM(diff(1:3,j)*moments(j-1:j+1,i))
        d2mean(j,i) = SUM(d2n(1:3,j)*moments(j-1:j+1,i))
     END DO
     dmean(ny,i) = SUM(diff(1:3,ny)*moments(ny-2:ny,i))
     d2mean(ny,i) = SUM(d2n(1:3,ny)*moments(ny-2:ny,i))
  END DO

  omega(:,1) = dmean(:,3)
  omega(:,2) = zero
  omega(:,3) = - dmean(:,1)

  domega(:,1) = d2mean(:,3)
  domega(:,2) = zero
  domega(:,3) = - d2mean(:,1)


  moments_vor(:,1:3) = omega

! COMPUTE MEAN VELOCITY AND VARIANCE OF FLUCTUATIONS.
  re_tau = sqrt(re*SUM(diff(1:3,1)*moments(1:3,1)))
  u_tau  = sqrt(SUM(diff(1:3,1)*moments(1:3,1))/re)
  write(*,*) re_tau, u_tau, drag(1), press_grad
  re_tau = sqrt(re*drag(1))
  u_tau  = sqrt(drag(1)/re)
  write(*,*) re_tau, u_tau, drag(1), press_grad

  DO n = 1,num_datasets

     wgt_old = DBLE(n-1)/DBLE(n)
     wgt_new =       one/DBLE(n)

! READ IN THE FILENAME OF THE nTH DATA SET.
     status = dsrref('datafiles.hdf',n+1)
     status = dsgdata('datafiles.hdf',fname_rank,fname_length,fname)
     
! READ IN THE TIME
     status = dsgdata(fname,time_rank,time_length,time)
     WRITE(*,*) n, fname, '  TIME = ', time(1)
     drag(1:2) = wgt_old*drag(1:2) + wgt_new*time(2:3)
     press_grad = wgt_old*press_grad + wgt_new*time(4)

! READ IN THE VELOCITY/SCALAR/PRESSURE ARRAY
     nn = 3+nscalar
     DO i = 1,nn+1
        status = dsgdata(fname,vel_rank,vel_length,u(:,:,:,i))
     END DO

!!$     DO z = 1,nz
!!$        WRITE(20,*) (u(nx/2,k,z,1),k=1,ny/2)
!!$     END DO

! COMPUTE VARIANCE, SKEWNESS AND FLATNESS OF VELOCITY FIELD
     scale = one/DBLE(nx*nz)

! COMPUTE DERIVATIVES OF FLUCTUATING VELOCITY/SCALAR FIELDS PLANE-BY-PLANE
!  AND COMPUTE STATISTICS/BUDGETS ON THESE PLANES.

     DO k = 1,ny

! FIRST COMPUTE DERIVATIVES AT THE BOTTOM WALL:
        DO i = 1,3+nscalar
           temp_u = u(1:nx,k,1:nz,i) - moments(k,i)
           CALL compute_xz_derivatives(temp_u, &
                grad_u(1:nx,1:nz,i,1),grad_u(1:nx,1:nz,i,3))
           IF (k==1) THEN
              grad_u(1:nx,1:nz,i,2) = &
                     diff(1,1)*(u(1:nx,1,1:nz,i)-moments(1,i)) &
                   + diff(2,1)*(u(1:nx,2,1:nz,i)-moments(2,i)) &
                   + diff(3,1)*(u(1:nx,3,1:nz,i)-moments(3,i)) 
           ELSEIF (k == ny) THEN
              grad_u(1:nx,1:nz,i,2) = &
                     diff(1,ny)*(u(1:nx,ny-2,1:nz,i)-moments(ny-2,i)) &
                   + diff(2,ny)*(u(1:nx,ny-1,1:nz,i)-moments(ny-1,i)) &
                   + diff(3,ny)*(u(1:nx,ny,  1:nz,i)-moments(ny,  i)) 
           ELSE              
              grad_u(1:nx,1:nz,i,2) = &
                     diff(1,k)*(u(1:nx,k-1,1:nz,i)-moments(k-1,i)) &
                   + diff(2,k)*(u(1:nx,k,  1:nz,i)-moments(k,  i)) &
                   + diff(3,k)*(u(1:nx,k+1,1:nz,i)-moments(k+1,i)) 
           END IF
           IF (i <= 3) THEN ! COMPUTE SECOND PARTIALS OF VELOCITY
              DO j = 1,3    ! USED IN COMPUTING VORTICITY GRADIENTS
                 CALL compute_xz_derivatives(grad_u(1:nx,1:nz,i,j),   &
                      d2_u(1:nx,1:nz,i,j,1),d2_u(1:nx,1:nz,i,j,3))
              END DO
              IF (k==1) THEN
                 d2_u(1:nx,1:nz,i,2,2) = &
                        d2n(1,1)*(u(1:nx,1,1:nz,i)-moments(1,i)) &
                      + d2n(2,1)*(u(1:nx,2,1:nz,i)-moments(2,i)) &
                      + d2n(3,1)*(u(1:nx,3,1:nz,i)-moments(3,i)) 
              ELSEIF (k == ny) THEN
                 d2_u(1:nx,1:nz,i,2,2) = &
                        d2n(1,ny)*(u(1:nx,ny-2,1:nz,i)-moments(ny-2,i)) &
                      + d2n(2,ny)*(u(1:nx,ny-1,1:nz,i)-moments(ny-1,i)) &
                      + d2n(3,ny)*(u(1:nx,ny,  1:nz,i)-moments(ny,  i)) 
              ELSE              
                 d2_u(1:nx,1:nz,i,2,2) = &
                        d2n(1,k)*(u(1:nx,k-1,1:nz,i)-moments(k-1,i)) &
                      + d2n(2,k)*(u(1:nx,k,  1:nz,i)-moments(k,  i)) &
                      + d2n(3,k)*(u(1:nx,k+1,1:nz,i)-moments(k+1,i)) 
              END IF
              d2_u(1:nx,1:nz,i,1,2) = d2_u(1:nx,1:nz,i,2,1)
              d2_u(1:nx,1:nz,i,3,2) = d2_u(1:nx,1:nz,i,2,3)
           END IF
        END DO

        CALL compute_velocity_moments(k)
        CALL compute_stress_tensor
        CALL compute_vorticity
        CALL compute_vorticity_moments(k)
        CALL compute_velocity_budget(k)
        CALL compute_vorticity_budget(k)
        CALL compute_scalar_budget(k)
        CALL compute_taylor_scale(k)
     
     END DO

  END DO

! NORMALIZE SKEWNESS AND FLATNESS OF VELOCITY/SCALAR AND VORTICITY FIELDS.
  nn = 3+nscalar
  DO i = 1,3+nscalar
     moments(1:ny,2*nn+i) = moments(1:ny,2*nn+i) &
          /(moments(1:ny,nn+i)**(1.5d0) + EPSILON(moments(1:ny,nn+i)))
     moments(1:ny,3*nn+i) = moments(1:ny,3*nn+i) &
          /(moments(1:ny,nn+i)**2 + EPSILON(moments(1:ny,nn+i)))
  END DO
  DO i = 1,3
     moments_vor(1:ny,2*nn+i) = moments_vor(1:ny,2*nn+i) &
          /(moments_vor(1:ny,nn+i)**(1.5d0) + EPSILON(moments_vor(1:ny,nn+i)))
     moments_vor(1:ny,3*nn+i) = moments_vor(1:ny,3*nn+i) &
          /(moments_vor(1:ny,nn+i)**2 + EPSILON(moments_vor(1:ny,nn+i)))
  END DO

! COMPUTE DERIVATIVES OF TRANSPORT TERMS IN BUDGETS
  DO i = 3,5
     CALL compute_y_derivative(budget_vel(1:ny,i),temp(1:ny))
     budget_vel(1:ny,i) = temp(1:ny)     
  END DO
  DO i = 1,nscalar
     DO j = 3,4
        CALL compute_y_derivative(budget_scalar(1:ny,4*(i-1)+j),temp(1:ny))
        budget_scalar(1:ny,4*(i-1)+j) = temp(1:ny)     
     END DO
  END DO
  DO i = 6,7
     CALL compute_y_derivative(budget_vor(1:ny,i),temp(1:ny))
     budget_vor(1:ny,i) = temp(1:ny)     
  END DO

! COMPUTE REYNOLDS NUMBER BASED ON TAYLOR MICROSCALE
  taylor_scale(1:ny,1) = re*moments(1:ny,nn+1)/sqrt(taylor_scale(1:ny,1) &
       +EPSILON(taylor_scale(1:ny,1)))
  DO i = 1,nscalar
     taylor_scale(1:ny,1+i) = &
          re*pr(i)*moments(1:ny,nn+3+i)/sqrt(taylor_scale(1:ny,1+i) &
          +EPSILON(taylor_scale(1:ny,1+i)))
  END DO

  yplus = re_tau*(one+yh(1:ny/2))

! WRITE OUT MEAN PROFILES
  output_file = 'meanvel.out'
  CALL save_formatted_output(output_file,yh,    moments(1:ny,1:nn))
  output_file = 'meanplus.out'
  CALL save_formatted_output(output_file,yplus,moments(1:ny/2,1:nn)/u_tau)
  output_file = 'fluct.out'
  CALL save_formatted_output(output_file,yplus,&
       moments(1:ny/2,nn+1:2*nn)/u_tau**2)
  output_file = 'skewness.out'
  CALL save_formatted_output(output_file,yplus,moments(1:ny/2,2*nn+1:3*nn))
  output_file = 'flatness.out'
  CALL save_formatted_output(output_file,yplus,moments(1:ny/2,3*nn+1:4*nn))
  output_file = 'fluxes.out'
  CALL save_formatted_output(output_file,yplus,fluxes(1:ny/2,1:nn)/u_tau**2)
  output_file = 'vor_moments.out'
  CALL save_formatted_output(output_file,yplus,moments_vor(1:ny/2,1:12))
  output_file = 'budget_vel.out'
  CALL save_formatted_output(output_file,yh,budget_vel)
  output_file = 'budget_vor.out'
  CALL save_formatted_output(output_file,yh,budget_vor)
  IF (nscalar >= 1) THEN
     output_file = 'budget_t.out'
     CALL save_formatted_output(output_file,yh,budget_scalar)
  END IF
  output_file = 'r_lambda.out'
  CALL save_formatted_output(output_file,yh,taylor_scale)

CONTAINS
!=================================================
  SUBROUTINE save_formatted_output(filename,coord,stats)
    IMPLICIT NONE
    CHARACTER(LEN=15)        :: filename    
    REAL(KIND=prec), DIMENSION(:)   :: coord
    REAL(KIND=prec), DIMENSION(:,:) :: stats

    open (unit = 12, file = filename, form = 'FORMATTED')
    DO j = 1,SIZE(coord)
       WRITE(12,990) coord(j), (stats(j,i),i=1,SIZE(stats(1,:)))
990    FORMAT(18f20.10)
    END DO
    close (unit = 12)
  END SUBROUTINE save_formatted_output
!=================================================
  SUBROUTINE setup_diff(n,z0,z,d)
    IMPLICIT NONE

    INTEGER            n
    REAL(KIND=prec) :: z0, z(n), d(n)

    INTEGER         i, j, k
    REAL(KIND=prec) temp

    d = zero
    DO i = 1,n
       DO j = 1,n
          IF (j /= i) THEN
             temp = one
             DO k = 1,n
                IF ((k /= j) .and. (k /= i)) THEN
                   temp = temp*(z0-z(k))/(z(i)-z(k))
                END IF
             END DO
             d(i) = d(i) + temp/(z(i) - z(j))
          END IF
       END DO
    END DO

  END SUBROUTINE setup_diff
!=================================================
  SUBROUTINE compute_xz_derivatives(vel,dvdx,dvdz)
    IMPLICIT NONE

    REAL(KIND=prec), DIMENSION(nx,nz) :: vel, dvdx, dvdz
    COMPLEX(KIND=prec), DIMENSION(nx/2+1,nz) :: cvel, cdvdx, cdvdz

! TRANSFORM VELOCITY INTO FOURIER SPACE.
    CALL rfftwnd_f77_one_real_to_complex(output_p2f_plan,vel,cvel)
    cvel = cvel/DBLE(nx*nz)

! COMPUTE X- AND Z-DERIVATIVES IN FOURIER SPACE
    DO z = 1,nz
       cdvdx(1:nx/2,z) = ikx(1:nx/2)*cvel(1:nx/2,z)
       cdvdz(1:nx/2,z) = ikz(z)     *cvel(1:nx/2,z)
    END DO

! TRANSFORM X-DERIVATIVE BACK INTO PHYSICAL SPACE
    cdvdx(nx/2+1,1:nz)   = zero
    cdvdx(1:nx/2,nz/2+1) = zero
    CALL rfftwnd_f77_one_complex_to_real(output_f2p_plan,cdvdx,dvdx)

! TRANSFORM Z-DERIVATIVE BACK INTO PHYSICAL SPACE
    cdvdz(nx/2+1,1:nz)   = zero
    cdvdz(1:nx/2,nz/2+1) = zero
    CALL rfftwnd_f77_one_complex_to_real(output_f2p_plan,cdvdz,dvdz)

  END SUBROUTINE compute_xz_derivatives
!==============================================
  SUBROUTINE compute_stress_tensor
    IMPLICIT NONE

     DO j = 1,3
        DO i = 1,3
           sij(1:nx,1:nz,i,j) = &
                half*(grad_u(1:nx,1:nz,i,j)+grad_u(1:nx,1:nz,j,i))
        END DO
     END DO

  END SUBROUTINE compute_stress_tensor
!==============================================
  SUBROUTINE compute_vorticity
    IMPLICIT NONE

    vor(1:nx,1:nz,1) = grad_u(1:nx,1:nz,3,2) - grad_u(1:nx,1:nz,2,3)
    vor(1:nx,1:nz,2) = grad_u(1:nx,1:nz,1,3) - grad_u(1:nx,1:nz,3,1)
    vor(1:nx,1:nz,3) = grad_u(1:nx,1:nz,2,1) - grad_u(1:nx,1:nz,1,2)

    DO i = 1,3
       grad_vor(1:nx,1:nz,1,i) = d2_u(1:nx,1:nz,3,2,i) - d2_u(1:nx,1:nz,2,3,i)
       grad_vor(1:nx,1:nz,2,i) = d2_u(1:nx,1:nz,1,3,i) - d2_u(1:nx,1:nz,3,1,i)
       grad_vor(1:nx,1:nz,3,i) = d2_u(1:nx,1:nz,2,1,i) - d2_u(1:nx,1:nz,1,2,i)
    END DO

  END SUBROUTINE compute_vorticity
!==============================================
  SUBROUTINE compute_velocity_moments(plane)
    IMPLICIT NONE
    INTEGER plane

    scale = one/DBLE(nx*nz)
    DO i = 1,3+nscalar
       moments(plane,nn+i) = wgt_old*moments(plane,nn+i) &
            + wgt_new*scale*SUM((u(1:nx,plane,1:nz,i)-moments(plane,i))**2)
       moments(plane,2*nn+i) = wgt_old*moments(plane,2*nn+i) &
            + wgt_new*scale*SUM((u(1:nx,plane,1:nz,i)-moments(plane,i))**3)
       moments(plane,3*nn+i) = wgt_old*moments(plane,3*nn+i) &
            + wgt_new*scale*SUM((u(1:nx,plane,1:nz,i)-moments(plane,i))**4)
       fluxes(plane,i) = wgt_old*fluxes(plane,i) &
            + wgt_new*scale*SUM((u(1:nx,plane,1:nz,i)-moments(plane,i)) &
            *(u(1:nx,plane,1:nz,2)-moments(plane,2)))
    END DO

  END SUBROUTINE compute_velocity_moments
!==============================================
  SUBROUTINE compute_vorticity_moments(plane)
    IMPLICIT NONE
    INTEGER plane

    scale = one/DBLE(nx*nz)
    DO i = 1,3
       moments_vor(plane,3+i) = wgt_old*moments_vor(plane,3+i) &
            + wgt_new*scale*SUM(vor(1:nx,1:nz,i)**2)
       moments_vor(plane,2*3+i) = wgt_old*moments_vor(plane,2*3+i) &
            + wgt_new*scale*SUM(vor(1:nx,1:nz,i)**3)
       moments_vor(plane,3*3+i) = wgt_old*moments_vor(plane,3*3+i) &
            + wgt_new*scale*SUM(vor(1:nx,1:nz,i)**4)
    END DO

  END SUBROUTINE compute_vorticity_moments
!==============================================
  SUBROUTINE compute_velocity_budget(plane)
    IMPLICIT NONE
    INTEGER plane
    REAL(KIND=prec) tmp_sum

    scale = one/DBLE(nx*nz)

! NOTE THAT WE WILL HAVE TO COMPUTE THE y-DERIVATIVE OF THE TRANSPORT/
! PRESSURE QUANTITIES AT THE END OF THE PROGRAM TO COMPUTE THE TRUE
! TRANSPORT.  THE PRODUCTION WILL ALSO BE COMPUTED THEN FROM THE FLUXES
! AND MEAN VELOCITY PROFILE.

! PRODUCTION
    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum &
               + u(x,plane,z,2)*SUM(u(x,plane,z,1:3)*dmean(plane,1:3))
       END DO
    END DO
    budget_vel(plane,1) = wgt_old*budget_vel(plane,1) - wgt_new*scale*tmp_sum

! DISSIPATION
    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum + SUM(sij(x,z,1:3,1:3)**2)
       END DO
    END DO
    budget_vel(plane,2) = wgt_old*budget_vel(plane,2) &
         + wgt_new*scale*(two/re)*tmp_sum

! TRANSPORT
    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum + u(x,plane,z,2)*SUM(u(x,plane,z,1:3)**2)
       END DO
    END DO
    budget_vel(plane,3) = wgt_old*budget_vel(plane,3) &
         - wgt_new*scale*half*tmp_sum
    
! VISCOUS TRANSPORT
    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum + SUM(u(x,plane,z,1:3)*sij(x,z,1:3,2))
       END DO
    END DO
    budget_vel(plane,4) = wgt_old*budget_vel(plane,4) &
         + wgt_new*scale*(two/re)*tmp_sum
    
! PRESSURE WORK
    tmp_sum = SUM(u(1:nx,plane,1:nz,2)*u(1:nx,plane,1:nz,nn+1))
    budget_vel(plane,5) = wgt_old*budget_vel(plane,5) &
         - wgt_new*scale*tmp_sum
    
  END SUBROUTINE compute_velocity_budget
!==============================================
  SUBROUTINE compute_scalar_budget(plane)
    IMPLICIT NONE
    INTEGER plane
    REAL(KIND=prec) tmp_sum

    scale = one/DBLE(nx*nz)

! NOTE THAT WE WILL HAVE TO COMPUTE THE y-DERIVATIVE OF THE TRANSPORT/
! PRESSURE QUANTITIES AT THE END OF THE PROGRAM TO COMPUTE THE TRUE
! TRANSPORT.  THE PRODUCTION WILL ALSO BE COMPUTED THEN FROM THE FLUXES
! AND MEAN VELOCITY PROFILE.

    DO i = 1,nscalar

! PRODUCTION
       tmp_sum = SUM(u(1:nx,plane,1:nz,3+i)*u(1:nx,plane,1:nz,2))&
            *dmean(plane,3+i)
       budget_scalar(plane,4*i-3) = wgt_old*budget_scalar(plane,4*i-3) &
            - wgt_new*scale*tmp_sum

! DISSIPATION
       tmp_sum = SUM(grad_u(1:nx,1:nz,3+i,1:3)**2)
       budget_scalar(plane,4*i-2) = wgt_old*budget_scalar(plane,4*i-2)&
            + wgt_new*scale*tmp_sum/(re*pr(i))

! TRANSPORT
       tmp_sum = SUM(u(1:nx,plane,1:nz,3+i)**2*u(1:nx,plane,1:nz,2))
       budget_scalar(plane,4*i-1) = wgt_old*budget_scalar(plane,4*i-1)&
            - wgt_new*scale*half*tmp_sum
    
! VISCOUS TRANSPORT
       tmp_sum = SUM(u(1:nx,plane,1:nz,3+i)*grad_u(1:nx,1:nz,3+i,2))
       budget_scalar(plane,4*i) = wgt_old*budget_scalar(plane,4*i) &
            + wgt_new*scale*tmp_sum/(re*pr(i))
    
    END DO

  END SUBROUTINE compute_scalar_budget
!==============================================
  SUBROUTINE compute_vorticity_budget(plane)
    IMPLICIT NONE
    INTEGER plane
    REAL(KIND=prec) tmp_sum

    scale = one/DBLE(nx*nz)

! NOTE THAT WE WILL HAVE TO COMPUTE THE y-DERIVATIVE OF THE TRANSPORT/
! PRESSURE QUANTITIES AT THE END OF THE PROGRAM TO COMPUTE THE TRUE
! TRANSPORT.  THE PRODUCTION WILL ALSO BE COMPUTED THEN FROM THE FLUXES
! AND MEAN VELOCITY PROFILE.

! GRADIENT PRODUCTION
    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum &
               + u(x,plane,z,2)*SUM(vor(x,z,1:3)*domega(plane,1:3))
       END DO
    END DO
    budget_vor(plane,1) = wgt_old*budget_vor(plane,1) - wgt_new*scale*tmp_sum

    tmp_sum = zero
    DO i = 1,3    
       DO z = 1,nz
          DO x = 1,nx
             tmp_sum = tmp_sum + SUM(vor(x,z,i)*vor(x,z,1:3)*sij(x,z,i,1:3))
          END DO
       END DO
    END DO
    budget_vor(plane,2) = wgt_old*budget_vor(plane,2) &
         - wgt_new*scale*tmp_sum

    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum + vor(x,z,2)*SUM(vor(x,z,1:3)*dmean(plane,1:3))
       END DO
    END DO
    budget_vor(plane,3) = wgt_old*budget_vor(plane,3) &
         - wgt_new*scale*tmp_sum

    tmp_sum = zero
    DO i = 1,3    
       DO z = 1,nz
          DO x = 1,nx
             tmp_sum = tmp_sum+vor(x,z,i)*SUM(omega(plane,1:3)*sij(x,z,i,1:3))
          END DO
       END DO
    END DO
    budget_vor(plane,4) = wgt_old*budget_vor(plane,4) &
         - wgt_new*scale*tmp_sum

! DISSIPATION
    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum + SUM(grad_vor(x,z,1:3,1:3)**2)
       END DO
    END DO
    budget_vor(plane,5) = wgt_old*budget_vor(plane,5) &
         + wgt_new*scale*tmp_sum/re

! TRANSPORT
    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum + u(x,plane,z,2)*SUM(vor(x,z,1:3)**2)
       END DO
    END DO
    budget_vor(plane,6) = wgt_old*budget_vor(plane,6) &
         - wgt_new*scale*half*tmp_sum
    
! VISCOUS TRANSPORT
    tmp_sum = zero
    DO z = 1,nz
       DO x = 1,nx
          tmp_sum = tmp_sum + SUM(vor(x,z,1:3)*grad_vor(x,z,1:3,2))
       END DO
    END DO
    budget_vor(plane,7) = wgt_old*budget_vor(plane,7) &
         + wgt_new*scale*(one/re)*tmp_sum
    
  END SUBROUTINE compute_vorticity_budget
!==============================================
  SUBROUTINE compute_taylor_scale(plane)
    IMPLICIT NONE
    INTEGER plane

    scale = one/DBLE(nx*nz)
    taylor_scale(plane,1) = wgt_old*taylor_scale(plane,1) &
         + wgt_new*scale*SUM(grad_u(1:nx,1:nz,1,1)**2)
    DO i = 1,nscalar
       taylor_scale(plane,1+i) = wgt_old*taylor_scale(plane,1+i) &
         + wgt_new*scale*SUM(grad_u(1:nx,1:nz,1,3+i)**2)
    END DO

  END SUBROUTINE compute_taylor_scale
!=========================================
  SUBROUTINE compute_y_derivative(in, out)
    IMPLICIT NONE

    INTEGER kk
    REAL(KIND=prec), DIMENSION(ny) :: in, out

    out(1) = SUM(diff(1:3,1)*in(1:3))
    DO kk = 2,ny-1
       out(kk) = SUM(diff(1:3,kk)*in(kk-1:kk+1))       
    END DO
    out(ny) = SUM(diff(1:3,ny)*in(ny-2:ny))

  END SUBROUTINE compute_y_derivative
!==============================================
END PROGRAM statistics

