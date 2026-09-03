module nonlinear_solution
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE, DHALF, DTWO
  use types_module, only: sol_set
  use utility_module, only: st_mpi, slope_sign_num, nan_recv_num, maxstep_num, satlim_num
  use initial_module, only: st_ctrl
  use read_input, only: len_scal
  use check_condition, only: st_out_fnum
  use set_cell, only: ncalc
  use prep_calculation, only: st_time
  use allocate_solution, only: nreg_num, array_var
  use calc_function, only: qext_sum, calc_func
  use calc_simulation, only: calc_l2norm2
#ifdef MPI_MSG
  use mpi_utility, only: mpisum_val
#endif

  implicit none
  private
  public :: allocate_nonlin, calc_numsol

  ! -- local
  integer(I4), parameter :: MAXSTEP_RUN_MAX = 5
  real(DP), parameter :: VARMAX = 1.00E+03_DP, XMAX = 1.00E+04_DP
  real(DP), parameter :: XMAX_INV = 1.00E-04_DP
  real(DP), allocatable :: new_func(:), jacvec(:)

  contains

  subroutine allocate_nonlin()
  !*********************************************************************************************
  ! allocate_nonlin -- Allocate for nonlinear solution work arrays
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    allocate(new_func(ncalc), jacvec(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      new_func(i) = DZERO ; jacvec(i) = DZERO
    end do
    !$omp end parallel do

  end subroutine allocate_nonlin

  subroutine calc_numsol(st_kryl, st_amgt, st_coef, st_sol)
  !*********************************************************************************************
  ! calc_numsol -- Calculate numerical solution
  !*********************************************************************************************
    ! -- modules
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use constval_module, only: VARLEN
    use types_module, only: kryl_set, amgt_set, coef_set
    use utility_module, only: log_fnum, write_err_stop
    use initial_module, only: st_sim, st_out_step
    use make_linearsystem, only: make_matvec
    use check_simulation, only: check_abserrmax, check_residual
    use linear_solution, only: solve_linalg, in_iter
#ifdef MPI_MSG
    use mpi_solve, only: check_mpimaxerr, bcast_convinfo
#endif
    ! -- inout
    type(kryl_set), intent(inout) :: st_kryl
    type(amgt_set), intent(inout) :: st_amgt
    type(coef_set), intent(inout) :: st_coef
    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
    integer(I4) :: out_iter
    integer(I4) :: max_num, conv_fnum
    integer(I4) :: back_iter, beta_iter, maxstep_run
    character(VARLEN) :: cxyzn
    real(DP) :: max_var, max_unk, check_val
    real(DP) :: conv_dmat, conv_rhs, conv_head, conv_var
    real(DP) :: l2norm_new, l2norm_pre, l2norm_jac, lambda, eater, gradient, max_step
    real(DP) :: res_l1, qext_l1, mass_error
    logical :: back_flag, res_flag, maxs_flag
#ifdef MPI_MSG
    real(DP) :: sum_l2
    real(DP) :: var_max, unk_max, var_abs_max
#endif
    ! -- format
    10 format(//1x,"CURRENT TIME : ",es12.5,1x,"(",a,")",20x,"TIME STEP : ",&
              es12.5,1x,"(SEC)",/,1x,95("-"),/,1x,&
              " OUTER INNER BACK BETA    MAXIMUM           MAXIMUM   DIAGONAL RIGHT HAND    &
              &UNKNOWN       MASS",/,1x,&
              "                           CHANGE              CELL     MATRIX     VECTOR      &
              &VALUE      ERROR",/,1x,95("-"))
    11 format(1X,2(i6),2(i5),es11.3,a18,4(es11.3))
    12 format(1X,"Didn't converge due to maximum value or change")
    13 format(1X,"Stop due to maximum value or change in backtracking")
    14 format(1X,"Stop due to maximum number of nonlinear iteration")
    15 format(1X,"Didn't converge in steady state calculation")
    16 format(1X,"RESIDUAL  SUM OF |F| = ",es11.3)
    !-------------------------------------------------------------------------------------------
    conv_fnum = st_out_fnum%conv ; eater = DHALF ; maxstep_run = 0
    res_l1 = DZERO ; qext_l1 = DZERO ; mass_error = DZERO
    ! -- Set for backtracking (backtr)
      call set_backtr(st_sol, max_step)

    outer_loop : do out_iter = 1, st_ctrl%maxout_iter
      st_time%out_iter = out_iter

      if (st_time%out_iter == 1) then
        if (st_mpi%rank == 0) then
          write(conv_fnum,10) st_time%now_time, trim(st_sim%cal_unit), st_time%delt
        end if
        !$omp parallel do private(i)
        do i = 1, ncalc
          new_func(i) = DZERO
        end do
        !$omp end parallel do
      else
        l2norm_pre = l2norm_new
      end if

      if (st_ctrl%picard_iter < 0) then
        st_time%form_switch = 0
      else if (st_time%out_iter > st_ctrl%picard_iter) then
        st_time%form_switch = 1
      else
        st_time%form_switch = 0
      end if

      back_iter = 0 ; beta_iter = 0
      back_flag = .false.

      ! -- Reset coefficients matrix and constant vector (matvec)
        call reset_matvec()

      if (st_sim%sim_type == -1) then
        ! -- Calculate surface water level (surfw)
          call calc_surfw(st_sol)
      end if

      ! -- Make coefficients matrix and constant vector (matvec)
        call make_matvec(st_coef, st_sol)

      if (st_time%out_iter == 1) then
        ! -- Calculate l2 norm square (resl2norm2)
          call calc_l2norm2(1, array_var(1)%rhs, l2norm_new)
#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
          ! -- Sum value for MPI (val)
            call mpisum_val(l2norm_new, "initial function l2-norm", sum_l2)
          l2norm_new = sum_l2
        end if
#endif
        l2norm_pre = l2norm_new
      end if

      !$omp parallel do private(i)
      do i = 1, nreg_num
        st_sol%head_pre(i) = st_sol%head_new(i)
        st_sol%head_change(i) = DZERO
      end do
      !$omp end parallel do

      st_time%conv_flag = .false.
      if (st_sim%sim_type /= -1) then
        st_ctrl%errtol = eater
      else
        st_ctrl%errtol = XMAX_INV**3
      end if
      ! -- Solve linear algebra (linalg)
        call solve_linalg(l2norm_pre, st_sol%head_change, st_kryl, st_amgt, l2norm_jac)

      !$omp parallel do private(i)
      do i = 1, nreg_num
        st_sol%head_new(i) = st_sol%head_pre(i) + st_sol%head_change(i)
      end do
      !$omp end parallel do

      ! -- Check absolute error max norm
        call check_abserrmax(st_sol%head_new, st_sol%head_pre, max_var, max_unk, max_num)

#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Check mpi max error (mpimaxerr)
          call check_mpimaxerr(max_var, max_unk, var_abs_max, unk_max, var_max)
        check_val = var_abs_max*len_scal ; max_unk = unk_max*len_scal
      else
        check_val = abs(max_var)*len_scal ; max_unk = abs(max_unk)*len_scal
        var_max = max_var
      end if
#else
      check_val = abs(max_var)*len_scal ; max_unk = abs(max_unk)*len_scal
#endif
      if (st_time%conv_flag .and. .not. ieee_is_nan(max_unk)) then
          st_time%conv_flag = .false.
        if (check_val <= st_ctrl%criteria .and. max_unk < XMAX) then
          st_time%conv_flag = .true.
        end if
      else if (ieee_is_nan(max_unk) .and. st_sim%sim_type /= -1) then
        if (st_mpi%rank == 0) then
          write(log_fnum,'(a)') "Nan detected."
        end if
        exit outer_loop
      else if (ieee_is_nan(max_unk) .and. st_sim%sim_type == -1) then
        if (st_mpi%rank == 0) then
          write(conv_fnum,15)
        end if
        call write_err_stop("Nan detected in the steady state calculation.")
      end if

      if (.not. st_time%conv_flag .and. st_time%form_switch == 1) then
        ! -- Run backtracking (backtr)
          call run_backtr(back_iter, back_flag, beta_iter, maxs_flag, l2norm_new, l2norm_pre,&
                          l2norm_jac, lambda, gradient, max_step, st_coef%stod, new_func,&
                          st_sol)
        if (maxs_flag) then
          maxstep_run = maxstep_run + 1
        else
          maxstep_run = 0
        end if
        if (maxstep_run == MAXSTEP_RUN_MAX) then
          back_flag = .true.
        end if
        ! -- Check absolute error max norm
          call check_abserrmax(st_sol%head_new, st_sol%head_pre, max_var, max_unk, max_num)
#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
          ! -- Check mpi max error (mpimaxerr)
            call check_mpimaxerr(max_var, max_unk, var_abs_max, unk_max, var_max)
          check_val = var_abs_max*len_scal ; max_unk = unk_max*len_scal
        else
          check_val = abs(max_var)*len_scal ; max_unk = abs(max_unk)*len_scal
          var_max = max_var
        end if
#else
        check_val = abs(max_var)*len_scal ; max_unk = abs(max_unk)*len_scal
#endif
        if ((check_val >= VARMAX .or. max_unk >= XMAX) .and. st_sim%sim_type /= -1) then
          back_flag = .true.
        else if (check_val <= st_ctrl%criteria .and. max_unk < XMAX .and.&
                 .not. (back_flag .and. st_sim%sim_type /= -1)) then
          st_time%conv_flag = .true.
        else if (st_sim%sim_type == -1) then
          back_flag = .false.
        end if
        if (.not. st_time%conv_flag .and. .not. back_flag) then
          ! -- Set Eisenstat-Walker forcing term (eise_walk)
            call set_eise_walk(l2norm_new, l2norm_pre, l2norm_jac, gradient, eater)
        end if
      end if

      ! -- Calculate function value (func)
        call calc_func(st_sol%stor_old, st_sol%stor_new, st_sol%surf_head, st_sol%head_new,&
                       st_sol%srat_new, st_sol%rel_perm, st_sol%surf_rati, new_func)
      ! -- Calculate l2 norm square (resl2norm2)
        call calc_l2norm2(1, new_func, l2norm_new)
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(l2norm_new, "new function l2-norm", sum_l2)
        l2norm_new = sum_l2
      end if
      if (max_var == var_max) then
        cxyzn = get_cnum(max_num)
      else
        cxyzn = ""
      end if
      conv_dmat = array_var(1)%dmat(max_num)*len_scal**2
      conv_rhs = array_var(1)%rhs(max_num)*len_scal**3
      conv_head = st_sol%head_new(max_num)*len_scal
      ! -- Bcast converge information (convinfo)
        call bcast_convinfo(cxyzn, conv_dmat, conv_rhs, conv_head, max_var)
#else
      cxyzn = get_cnum(max_num)
      conv_dmat = array_var(1)%dmat(max_num)*len_scal**2
      conv_rhs = array_var(1)%rhs(max_num)*len_scal**3
      conv_head = st_sol%head_new(max_num)*len_scal
#endif
      conv_var = max_var*len_scal
      res_l1 = DZERO
      !$omp parallel do private(i) reduction(+:res_l1)
      do i = 1, ncalc
        res_l1 = res_l1 + abs(new_func(i))
      end do
      !$omp end parallel do
      qext_l1 = qext_sum
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(res_l1, "residual l1-norm", sum_l2)
        res_l1 = sum_l2
        ! -- Sum value for MPI (val)
          call mpisum_val(qext_l1, "external flux l1-norm", sum_l2)
        qext_l1 = sum_l2
      end if
#endif
      res_l1 = res_l1*len_scal**3 ; qext_l1 = qext_l1*len_scal**3
      if (qext_l1 > DZERO) then
        mass_error = res_l1/qext_l1
      else
        mass_error = DZERO
      end if
      if (st_mpi%rank == 0) then
        write(conv_fnum,11) st_time%out_iter, in_iter, back_iter, beta_iter, conv_var,&
                            trim(adjustl(cxyzn)), conv_dmat, conv_rhs, conv_head, mass_error
      end if

      if (st_time%conv_flag .and. st_ctrl%conv_type == 1) then
        ! -- Check residual convergence (residual)
          call check_residual(new_func, st_sol%stor_new, res_flag)
        st_time%conv_flag = res_flag
      end if

      ! check outer_loop
      if (st_time%conv_flag) then
        ! -- Calculate surface water level (surfw)
          call calc_surfw(st_sol)
        if (st_out_step%rest == DZERO) then
          ! -- Write restart file (rest)
            call write_rest(st_sol%head_new)
        else if (mod(st_time%current_t,st_out_step%rest) == 0) then
          ! -- Write restart file (rest)
            call write_rest(st_sol%head_new)
        end if
        exit outer_loop
      else if (back_flag .and. st_sim%sim_type /= -1) then
        if (st_mpi%rank == 0) then
          write(conv_fnum,13)
        end if
        exit outer_loop
      else if (st_time%out_iter == st_ctrl%maxout_iter .and. st_sim%sim_type == -1) then
        if (st_mpi%rank == 0) then
          write(conv_fnum,15)
        end if
        call write_err_stop("Steady state calculation didn't converge.")
      else if (st_time%out_iter == st_ctrl%maxout_iter) then
        if (st_mpi%rank == 0) then
          write(conv_fnum,14)
        end if
        exit outer_loop
      else if ((abs(conv_var) >= VARMAX .or. max_unk >= XMAX) .and. st_sim%sim_type /= -1) then
        if (st_mpi%rank == 0) then
          write(conv_fnum,12)
        end if
        exit outer_loop
      end if

    end do outer_loop

    if (st_mpi%rank == 0) then
      write(conv_fnum,16) res_l1
    end if


  end subroutine calc_numsol

  subroutine reset_matvec
  !*********************************************************************************************
  ! reset_matvec -- Reset coefficients matrix and constant vector
  !*********************************************************************************************
    ! -- modules
    use set_cell, only: amg_setflag
    use allocate_solution, only: crs_index, pro_var, res_var
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    array_var(1)%lumat(:) = DZERO
    array_var(1)%dmat(:) = DZERO
    array_var(1)%rhs(:) = DZERO

    if (st_ctrl%precon_type == 1 .and. amg_setflag == 1) then
      do i = 2, st_ctrl%nlevel
        deallocate(array_var(i)%dmat, array_var(i)%rhs, array_var(i)%lumat, array_var(i)%x)
        deallocate(crs_index(i)%offrow, crs_index(i)%offind)
        deallocate(pro_var(i)%pindex, pro_var(i)%poffrow, pro_var(i)%pval)
        deallocate(res_var(i)%rindex, res_var(i)%roffrow, res_var(i)%rval)
      end do
    end if

  end subroutine reset_matvec

  subroutine calc_surfw(st_sol)
  !*********************************************************************************************
  ! calc_surfw -- Calculate surface water level
  !*********************************************************************************************
    ! -- modules
    use set_cell, only: ncals
!    use make_cell, only: surf_elev
!    use prep_calculation, only: surf_bott
    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, ncals
!      if (st_sol%head_new(i) <= surf_elev(i)) then
!        st_sol%surf_head(i) = st_sol%head_new(i)
!      else
!        st_sol%surf_head(i) = surf_elev(i)
!      end if
      ! all surface head = surf_elev
!      st_sol%surf_head(i) = surf_elev(i)
      ! all surface head = surf_bottom
!      st_sol%surf_head(i) = surf_bott(i)
      ! all surface head = st_sol%head_new
      st_sol%surf_head(i) = st_sol%head_new(i)
    end do
    !$omp end parallel do

  end subroutine calc_surfw

  subroutine set_backtr(st_sol, maxstep)
  !*********************************************************************************************
  ! set_backtr -- Set for backtracking
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    type(sol_set), intent(in) :: st_sol
    real(DP), intent(out) :: maxstep
    ! -- local
    real(DP) :: l2_xnew
#ifdef MPI_MSG
    real(DP) :: sum_l2
#endif
    !-------------------------------------------------------------------------------------------
    ! -- Calculate l2 norm square (resl2norm2)
      call calc_l2norm2(1, st_sol%head_new, l2_xnew)
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_xnew, "previous function l2-norm", sum_l2)
      l2_xnew = sum_l2
    end if
#endif

    maxstep = sqrt(l2_xnew)*VARMAX
    if (maxstep < DONE) then
      maxstep = DONE
    end if

  end subroutine set_backtr

  subroutine set_eise_walk(l2_new, l2_pre, l2_jac, grad, eta)
  !*********************************************************************************************
  ! set_eise_walk -- Set Eisenstat-Walker forcing term
  !*********************************************************************************************
    ! -- modules
    use kind_module, only: SP
    use utility_module, only: log_fnum
    ! -- inout
    real(DP), intent(in) :: l2_new, l2_pre, l2_jac, grad
    real(DP), intent(inout) :: eta
    ! -- local
    real(DP), parameter :: ETA_MAX = 0.9_DP
    real(DP), parameter :: ETA_MIN = 1.0E-4_DP
    real(DP), parameter :: ETA_ALPHA = (1.0_DP+sqrt(5.0_DP))*DHALF
    real(DP), parameter :: DESCENT_TOL = epsilon(1.00_SP)
    real(DP) :: eta_safe, l2_line, lin_l2norm
    !-------------------------------------------------------------------------------------------
    eta_safe = eta**ETA_ALPHA
    l2_line = l2_pre + DTWO*grad + l2_jac
    if (l2_line < -DESCENT_TOL*l2_pre) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Warning!! Negative linear model norm in the forcing term."
      end if
    end if
    lin_l2norm = sqrt(max(DZERO, l2_line))
    eta = abs(sqrt(l2_new) - lin_l2norm)/sqrt(l2_pre)

    if (eta_safe < 0.1_DP) then
      eta_safe = DZERO
    end if

    eta = max(eta, eta_safe)
    eta = max(eta, ETA_MIN)
    eta = min(eta, ETA_MAX)

  end subroutine set_eise_walk

  subroutine run_backtr(backi, backf, betai, maxsf, l2_new, l2_pre, l2_jac, lam, grad, maxstep,&
                        stod, new_f, st_sol)
  !*********************************************************************************************
  ! run_backtr -- Run backtracking
  !*********************************************************************************************
    ! -- modules
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use make_cell, only: st_geom
    use calc_function, only: calc_vecjacf
#ifdef MPI_MSG
    use mpi_utility, only: mpimax_val, mpimin_val
#endif
    ! -- inout
    integer(I4), intent(inout) :: backi, betai
    logical, intent(inout) :: backf
    logical, intent(out) :: maxsf
    real(DP), intent(inout) :: l2_new, l2_jac, lam
    real(DP), intent(in) :: l2_pre, maxstep, stod(:)
    real(DP), intent(out) :: grad
    real(DP), intent(inout) :: new_f(:)
    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
    real(DP) :: l2_new2, l2_pnorm, slope, av, bv, rhs1, rhs2, root, step_len
    real(DP) :: f1_pre, f1_new
    real(DP) :: lam2, temp_lam, lam_inv, lam2_inv, del_lam, lam_max, lam_min
    real(DP) :: lam_length, lam_base, lam_diff, lam_incr, sql2_pnorm, maxpnorm
    real(DP) :: lam_maxi, dsat_scale, dsat_cell
    real(DP) :: alpha_cond, beta_cond
    real(DP), parameter :: BACK_ALPHA = 1.00E-4_DP, BACK_BETA = 0.9_DP, MAXSTEP_RATIO = 0.99_DP
    real(DP), parameter :: DSAT_STEP_FRAC = 0.9_DP
    real(DP), parameter :: DIVERGE_LIMIT = huge(1.00_DP)*0.1_DP
    real(DP), parameter :: STEP_TOL = 1.00E-07_DP
#ifdef MPI_MSG
    real(DP) :: sum_l2, max_val, min_val
#endif
    !-------------------------------------------------------------------------------------------
    l2_new2 = l2_new ; lam = DONE ; lam2 = DZERO ; maxpnorm = DONE ; lam_length = DZERO
    maxsf = .false.
    !$omp parallel do private(i)
    do i = 1, ncalc
      jacvec(i) = DZERO
    end do
    !$omp end parallel do
    ! -- Calculate function value (func)
      call calc_func(st_sol%stor_old, st_sol%stor_new, st_sol%surf_head, st_sol%head_new,&
                     st_sol%srat_new, st_sol%rel_perm, st_sol%surf_rati, new_f)
    ! -- Calculate l2 norm square (resl2norm2)
      call calc_l2norm2(1, st_sol%head_change, l2_pnorm)
    ! -- Calculate vector by jacobi-free (vecjocf)
      call calc_vecjacf(1, st_sol%head_change, st_sol%stor_old, st_sol%stor_new,&
                        st_sol%surf_head, st_sol%head_pre, st_sol%srat_new, st_sol%rel_perm,&
                        st_sol%surf_rati, jacvec)

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_pnorm, "change function l2-norm", sum_l2)
      l2_pnorm = sum_l2
    end if
#endif

    sql2_pnorm = sqrt(l2_pnorm)

    lam_maxi = DONE
    if (st_ctrl%expd_type == 1 .and. sql2_pnorm > DZERO) then
      lam_maxi = maxstep/sql2_pnorm
    end if

    if (sql2_pnorm > maxstep) then
      maxpnorm = maxstep/sql2_pnorm
    end if

    if (st_ctrl%dsat_max > DZERO) then
      dsat_scale = DONE
      !$omp parallel do private(i, dsat_cell) reduction(min:dsat_scale)
      do i = 1, ncalc
        dsat_cell = abs(stod(i))*st_time%delt/st_geom%cell_vol(i)*abs(st_sol%head_change(i))
        if (dsat_cell > st_ctrl%dsat_max) then
          dsat_scale = min(dsat_scale, DSAT_STEP_FRAC*st_ctrl%dsat_max/dsat_cell)
        end if
      end do
      !$omp end parallel do
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- MIN value for MPI (val)
          call mpimin_val(dsat_scale, "saturation limit ratio", min_val)
        dsat_scale = min_val
      end if
#endif
      if (dsat_scale < DONE) then
        satlim_num = satlim_num + 1
        maxpnorm = min(maxpnorm, dsat_scale)
      end if
    end if

    if (maxpnorm < DONE) then
      !$omp parallel do private(i)
      do i = 1, ncalc
        st_sol%head_change(i) = st_sol%head_change(i)*maxpnorm
      end do
      !$omp end parallel do
      sql2_pnorm = sql2_pnorm*maxpnorm
      lam_maxi = DONE
    end if
    step_len = sql2_pnorm

    l2_new = DZERO ; slope = DZERO ; l2_jac = DZERO
    !$omp parallel
    !$omp do private(i) reduction(+:l2_new, slope, l2_jac)
    do i = 1, ncalc
      l2_new = l2_new + new_f(i)*new_f(i)
      slope = slope + array_var(1)%rhs(i)*jacvec(i)*maxpnorm
      l2_jac = l2_jac + jacvec(i)*jacvec(i)
    end do
    !$omp end do

    !$omp do private(i, temp_lam) reduction(max:lam_length)
    do i = 1, ncalc
      temp_lam = abs(st_sol%head_change(i))/(DONE + abs(st_sol%head_pre(i)))
      if (temp_lam > lam_length) then
        lam_length = temp_lam
      end if
    end do
    !$omp end do
    !$omp end parallel

    l2_jac = l2_jac*maxpnorm*maxpnorm


#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_new, "current function l2-norm", sum_l2)
      l2_new = sum_l2
      ! -- Sum value for MPI (val)
        call mpisum_val(slope, "slope function l2-norm", sum_l2)
      slope = sum_l2
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_jac, "jacobi vector l2-norm", sum_l2)
      l2_jac = sum_l2
      ! -- MAX value for MPI (val)
        call mpimax_val(lam_length, "lambda length", max_val)
      lam_length = max_val
    end if
#endif
    if (lam_length > DZERO) then
      lam_min = STEP_TOL/lam_length
    else
      lam_min = DZERO
    end if
    if (slope >= DZERO) then
      slope_sign_num = slope_sign_num + 1
      if (slope > DZERO) then
        slope = -slope
      end if
    end if
    f1_pre = DHALF*l2_pre
    grad = slope
    back_aloop: do
      ! -- Calculate function and l2norm2 (func2norm)
        call calc_funcl2norm(lam, l2_new, new_f, st_sol)
      backi = backi + 1
      if (.not. ieee_is_finite(l2_new) .or. l2_new > DIVERGE_LIMIT) then
        nan_recv_num = nan_recv_num + 1
        if (lam < lam_min) then
          ! -- Calculate function and l2norm2 (func2norm)
            call calc_funcl2norm(DZERO, l2_new, new_f, st_sol)
          backf = .true.
          return
        end if
        lam = DHALF*lam
        cycle back_aloop
      end if
      f1_new = DHALF*l2_new
      alpha_cond = f1_pre + BACK_ALPHA*lam*slope
      if (f1_new <= alpha_cond) then
        exit back_aloop
      end if
      if (lam < lam_min) then
        ! -- Calculate function and l2norm2 (func2norm)
          call calc_funcl2norm(DZERO, l2_new, new_f, st_sol)
        backf = .true.
        return
      end if

      if (lam == DONE) then
        temp_lam = -slope/(DTWO*(f1_new-f1_pre-slope))
      else
        rhs1 = f1_new - f1_pre - lam*slope
        rhs2 = l2_new2 - f1_pre - lam2*slope
        lam_inv = DONE/(lam**2) ; lam2_inv = DONE/(lam2**2)
        del_lam = DONE/(lam - lam2)
        av = (rhs1*lam_inv - rhs2*lam2_inv)*del_lam
        bv = (-lam2*rhs1*lam_inv + lam*rhs2*lam2_inv)*del_lam
        if (av == 0) then
          temp_lam = -slope/(DTWO*bv)
        else
          root = bv*bv - 3.0_DP*av*slope
          if (root < DZERO) then
            temp_lam = DHALF*lam
          else if (bv <= DZERO) then
            temp_lam = (-bv + sqrt(root))/(3.0_DP*av)
          else
            temp_lam = -slope/(bv + sqrt(root))
          end if
          if (temp_lam > DHALF*lam) then
            temp_lam = DHALF*lam
          end if
        end if
      end if
      lam2 = lam
      l2_new2 = f1_new
      lam = max(temp_lam, 0.1_DP*lam)
    end do back_aloop

    alpha_cond = DHALF*l2_pre + BACK_ALPHA*lam*slope
    beta_cond = DHALF*l2_pre + BACK_BETA*lam*slope
    if (DHALF*l2_new < beta_cond) then
      if (lam == DONE .and. lam_maxi > DONE) then
        lam_max = lam_maxi
        b1_loop: do
          if (DHALF*l2_new > alpha_cond .or. DHALF*l2_new >= beta_cond .or. &
              lam >= lam_max) then
            exit b1_loop
          end if
          lam2 = lam ; l2_new2 = DHALF*l2_new ; lam = min(DTWO*lam, lam_max)
          ! -- Calculate function and l2norm2 (func2norm)
            call calc_funcl2norm(lam, l2_new, new_f, st_sol)
          alpha_cond = DHALF*l2_pre + BACK_ALPHA*lam*slope
          beta_cond = DHALF*l2_pre + BACK_BETA*lam*slope
        end do b1_loop
      end if

      if (lam < DONE .or. (lam > DONE .and. DHALF*l2_new > alpha_cond)) then
        lam_base = min(lam, lam2) ; lam_diff = abs(lam2 - lam)
        b2_loop: do
          if (DHALF*l2_new <= alpha_cond .and. (DHALF*l2_new >= beta_cond .or. &
              lam_diff <= lam_min)) then
            exit b2_loop
          end if
          lam_incr = DHALF*lam_diff ; lam = lam_base + lam_incr
          ! -- Calculate function and l2norm2 (func2norm)
            call calc_funcl2norm(lam, l2_new, new_f, st_sol)
          alpha_cond = DHALF*l2_pre + BACK_ALPHA*lam*slope
          beta_cond = DHALF*l2_pre + BACK_BETA*lam*slope

          if (DHALF*l2_new > alpha_cond) then
            lam_diff = lam_incr
          else if (DHALF*l2_new < beta_cond) then
            lam_base = lam ; lam_diff = lam_diff - lam_incr
          end if
          if (lam_diff == DZERO) then
            exit b2_loop
          end if
        end do b2_loop

        if (DHALF*l2_new < beta_cond .or.&
            (lam_diff < lam_min .and. DHALF*l2_new > alpha_cond)) then
          ! -- Calculate function and l2norm2 (func2norm)
            call calc_funcl2norm(lam_base, l2_new, new_f, st_sol)
          betai = betai + 1
        end if
        if (betai == 10) then
          backf = .true.
          return
        end if
      end if
    end if

    if (lam*sql2_pnorm > MAXSTEP_RATIO*maxstep) then
      maxsf = .true.
      maxstep_num = maxstep_num + 1
    end if

    grad = slope*lam
    l2_jac = l2_jac*lam*lam

  end subroutine run_backtr

  function get_cnum(calc_num) result(char_cell)
  !*********************************************************************************************
  ! get_cnum -- Get cell number
  !*********************************************************************************************
    ! -- modules
    use utility_module, only: get_ilen, conv_i2s
    use set_cell, only: get_calc_grid
    ! -- inout
    integer(I4), intent(in) :: calc_num
    ! -- local
    integer(I4) :: i_num, j_num, k_num
    character(:), allocatable :: cx_num, cy_num, cz_num
    character(:), allocatable :: char_cell
    !-------------------------------------------------------------------------------------------
    ! -- Get calculation number from grid number (calc_grid)
      call get_calc_grid(calc_num, i_num, j_num, k_num)

    allocate(character(get_ilen(i_num)) :: cx_num)
    allocate(character(get_ilen(j_num)) :: cy_num)
    allocate(character(get_ilen(k_num)) :: cz_num)

    call conv_i2s(i_num, cx_num) ; call conv_i2s(j_num, cy_num) ; call conv_i2s(k_num, cz_num)

    char_cell = "("//cx_num//","//cy_num//","//cz_num//")"

  end function get_cnum

  subroutine calc_funcl2norm(in_lam, l2_new, new_f, st_sol)
  !*********************************************************************************************
  ! calc_funcl2norm -- Calculate function and l2norm2
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: in_lam
    real(DP), intent(out) :: l2_new
    real(DP), intent(inout) :: new_f(:)
    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
#ifdef MPI_MSG
    real(DP) :: sum_l2
#endif
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, nreg_num
      st_sol%head_new(i) = st_sol%head_pre(i) + in_lam*st_sol%head_change(i)
    end do
    !$omp end parallel do

    ! -- Calculate function value (func)
      call calc_func(st_sol%stor_old, st_sol%stor_new, st_sol%surf_head, st_sol%head_new,&
                     st_sol%srat_new, st_sol%rel_perm, st_sol%surf_rati, new_f)
    ! -- Calculate l2 norm square (resl2norm2)
      call calc_l2norm2(1, new_f, l2_new)

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_new, "backtracking new function l2-norm", sum_l2)
      l2_new = sum_l2
    end if
#endif

  end subroutine calc_funcl2norm

  subroutine write_rest(rest_head)
  !*********************************************************************************************
  ! write_rest -- Write restart file
  !*********************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_write, only: write_mpi_rest
#endif
    ! -- inout
    real(DP), intent(in) :: rest_head(:)
    ! -- local
    integer(I4) :: i, rest_fnum
    !-------------------------------------------------------------------------------------------
    rest_fnum = st_out_fnum%rest
#ifdef MPI_MSG
    i = 0
    ! -- Write mpi restart value (mpi_rest)
      call write_mpi_rest(rest_fnum, st_time%now_time, len_scal, rest_head)
#else
    rewind(rest_fnum)
    write(rest_fnum) real(st_time%now_time, kind=DP)
    write(rest_fnum) (rest_head(i)*len_scal, i = 1, ncalc)
#endif

  end subroutine write_rest

end module nonlinear_solution
