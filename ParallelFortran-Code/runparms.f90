MODULE runparms
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: re,  pr,  pr_t, pi, nscalar,  split_scalar,      prec,     &
            t_total,     dt,          cfl_max,   p_grad,     beta,      &
            tend,        step,        rk_step,   rk_last,    dt_flag,   &
            output_step, save_step,   read_vel,  save_vel,   vel_out,   &
            which_test,  test_run,    fix_mass,  divg_check, nproc, me, &
            subgrid,     filter_dim,  lagr_avg,  dynamic_scalar, ndynscalar, &
            crossflow,   bl_grid,     filter_jet_profile,    forcing_phase, &
            jetmag_top,  jetmag_bot,  swirl_top,    swirl_bot,  jet_diam,  &
            forcing_top, forcing_bot, for_freq_top, for_freq_bot, &
            xjet_top,    xjet_bot,    zjet_top,  zjet_bot, &
            fringe_gain, xfringe,     zfringe,   frac_xbuff, frac_zbuff, &
            buff_xrise, buff_zrise, buff_xfall, buff_zfall, &
            ly, ly_half, delta_bl, nbl, linearized_run, snap_out, & 
            disturb_ampl, l2_step
            
! precision of floating point numbers.
  INTEGER, PARAMETER :: prec = SELECTED_REAL_KIND(14,32)

  REAL(KIND=prec)  re, lx, lz, beta
  REAL(KIND=prec)  t_total, dt, pi, p_grad, cfl_max, disturb_ampl
  REAL(KIND=prec), ALLOCATABLE, DIMENSION(:) :: pr, pr_t

  INTEGER, PARAMETER :: rk_last = 3
  INTEGER               tend,         step,      rk_step,   &
                        output_step,  save_step, l2_step, &
                        read_vel,     save_vel,  vel_out,   &
                        which_test,   test_run,  nscalar,   &
                        fix_mass,     dt_flag,   divg_check,&
                        split_scalar, subgrid,   filter_dim,&
                        lagr_avg,     dynamic_scalar, ndynscalar, &
                        linearized_run, snap_out

! mpi variables
  INTEGER    nproc, me

! JET IN CROSSFLOW/BOUNDARY LAYER DOMAIN PARAMETERS
  INTEGER    bl_grid, nbl
  REAL(KIND=prec)  ly, ly_half, delta_bl

! JET IN CROSSFLOW PARAMETERS
  INTEGER    crossflow, filter_jet_profile
  REAL(KIND=prec)  jet_diam, jetmag_top, jetmag_bot, swirl_top, swirl_bot, &
       forcing_top, forcing_bot, forcing_phase, for_freq_top, for_freq_bot, &
       xjet_top, xjet_bot, zjet_top, zjet_bot

! FRINGE/BUFFER REGION PARAMETERS
  INTEGER         xfringe, zfringe
  REAL(KIND=prec) fringe_gain, frac_xbuff, frac_zbuff, &
            buff_xrise, buff_zrise, buff_xfall, buff_zfall
END MODULE runparms



