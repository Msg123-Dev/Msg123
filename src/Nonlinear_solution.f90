module nonlinear_solution
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DHALF, DONE, DTWO, VARMAX
  use initial_module, only: criteria
  use read_input, only: len_scal
  use check_condition, only: st_out_fnum
  use set_cell, only: ncalc
  use prep_calculation, only: out_iter
  use allocate_solution, only: nreg_num, array_var, head_new, head_pre, head_change
  use time_module, only: now_time
  use calc_simulation, only: calc_l2norm2
  use calc_function, only: calc_func
#ifdef MPI_MSG
  use initial_module, only: pro_totn
  use mpi_utility, only: mpisum_val
#endif

  implicit none
  private
  public :: calc_numsol

  contains

  subroutine calc_numsol()
  !*********************************************************************************************
  ! calc_numsol -- Calculate numerical solution
  !*********************************************************************************************
    ! -- modules
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use constval_module, only: XMAX, XMAX_INV, VARLEN
    use utility_module, only: log_fnum
    use initial_module, only: st_sim, maxout_iter, errtol, st_out_step, my_rank
    use prep_calculation, only: current_t, delt, conv_flag, form_switch
    use make_linearsystem, only: make_matvec
    use check_simulation, only: check_abserrmax
    use linear_solution, only: solve_linalg, in_iter
#ifdef MPI_MSG
    use mpi_solve, only: check_mpimaxerr, bcast_convinfo
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4) :: max_num, conv_fnum
    integer(I4) :: back_iter, beta_iter
    character(VARLEN) :: cxyzn
    real(DP) :: max_var, max_unk, check_val
    real(DP) :: conv_dmat, conv_rhs, conv_head, conv_var
    real(DP) :: l2norm_new, l2norm_pre, l2norm_jac, lambda, eater, gradient, max_step
    real(DP), allocatable :: new_func(:)
    logical :: back_flag
#ifdef MPI_MSG
    real(DP) :: sum_l2
    real(DP) :: var_max, unk_max, var_abs_max
#endif
    ! -- format
    10 format(//1x,"CURRENT TIME : ",es12.5,1x,"(",a,")",20x,"TIME STEP : ",&
              es12.5,1x,"(SEC)",/,1x,84("-"),/,1x,&
              " OUTER INNER  BACK     MAXIMUM           MAXIMUM    DIAGONAL  RIGHT HAND     &
              &UNKNOWN",/,1x,&
              "                        CHANGE              CELL      MATRIX      VECTOR&
              &       VALUE",/,1x,84("-"))
    11 format(1X,3(i6),es12.3,a18,4(es12.3))
    12 format(1X,"Didn't converge due to maximum value or change")
    13 format(1X,"Stop due to maximum value or change in backtracking")
    14 format(1X,"Stop due to maximum number of nonlinear iteration")
    15 format(1X,"Didn't converge in steady state calculation")
    !-------------------------------------------------------------------------------------------
    conv_fnum = st_out_fnum%conv
    eater = DHALF
    ! -- Set for backtracking (backtr)
      call set_backtr(max_step)

    outer_loop : do out_iter = 1, maxout_iter

      if (out_iter == 1) then
        if (my_rank == 0) then
          write(conv_fnum,10) now_time, trim(st_sim%cal_unit), delt
        end if
        form_switch = 1
        allocate(new_func(ncalc))
        !$omp parallel do private(i)
        do i = 1, ncalc
          new_func(i) = DZERO
        end do
        !$omp end parallel do
      else
        l2norm_pre = l2norm_new
      end if

      back_iter = 0 ; beta_iter = 0
      back_flag = .false.

      ! -- Reset coefficients matrix and constant vector (matvec)
        call reset_matvec()

      if (st_sim%sim_type == -1) then
        ! -- Calculate surface water level (surfw)
          call calc_surfw()
      end if

      ! -- Make coefficients matrix and constant vector (matvec)
        call make_matvec()

      if (out_iter == 1) then
        ! -- Calculate l2 norm square (resl2norm2)
          call calc_l2norm2(1, array_var(1)%rhs, l2norm_new)
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Sum value for MPI (val)
            call mpisum_val(l2norm_new, "initial function l2-norm", sum_l2)
          l2norm_new = sum_l2
        end if
#endif
        l2norm_pre = l2norm_new
      end if

      !$omp parallel do private(i)
      do i = 1, nreg_num
        head_pre(i) = head_new(i)
        head_change(i) = DZERO
      end do
      !$omp end parallel do

      conv_flag = .false.
      if (st_sim%sim_type /= -1) then
        errtol = eater
      else
        errtol = XMAX_INV**3
      end if
      ! -- Solve linear algebra (linalg)
        call solve_linalg(l2norm_pre, head_change, l2norm_jac)

      !$omp parallel do private(i)
      do i = 1, nreg_num
        head_new(i) = head_pre(i) + head_change(i)
      end do
      !$omp end parallel do

      ! -- Check absolute error max norm
        call check_abserrmax(head_new, head_pre, max_var, max_unk, max_num)

#ifdef MPI_MSG
      if (pro_totn /= 1) then
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
      if (conv_flag .and. .not. ieee_is_nan(max_unk)) then
          conv_flag = .false.
        if (check_val <= criteria .and. max_unk < XMAX) then
          conv_flag = .true.
        end if
      else if (ieee_is_nan(max_unk) .and. st_sim%sim_type /= -1) then
        if (my_rank == 0) then
          write(log_fnum,'(a)') "Nan detected."
        end if
        exit outer_loop
      else if (ieee_is_nan(max_unk) .and. st_sim%sim_type == -1) then
        st_sim%sim_type = 0
        current_t = DZERO
        if (my_rank == 0) then
          write(log_fnum,'(a)') "Steady state calculation changes to timestep calculation."
          write(conv_fnum,15)
        end if
        exit outer_loop
      end if

      if (.not. conv_flag .and. form_switch == 1) then
        ! -- Run backtracking (backtr)
          call run_backtr(back_iter, back_flag, beta_iter, l2norm_new, l2norm_pre, l2norm_jac,&
                          lambda, gradient, max_step, new_func)
        ! -- Check absolute error max norm
          call check_abserrmax(head_new, head_pre, max_var, max_unk, max_num)
#ifdef MPI_MSG
        if (pro_totn /= 1) then
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
        else if (check_val <= criteria .and. max_unk < XMAX) then
          conv_flag = .true.
        else if (st_sim%sim_type == -1) then
          back_flag = .false.
        end if
        if (.not. conv_flag .and. .not. back_flag) then
          ! -- Set Eisenstat-Walker forcing term (eise_walk)
            call set_eise_walk(lambda, l2norm_new, l2norm_pre, l2norm_jac, gradient, eater)
        end if
      end if

      ! -- Calculate function value (func)
        call calc_func(head_new, new_func)
      ! -- Calculate l2 norm square (resl2norm2)
        call calc_l2norm2(1, new_func, l2norm_new)
#ifdef MPI_MSG
      if (pro_totn /= 1) then
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
      conv_head = head_new(max_num)*len_scal
      ! -- Bcast converge information (convinfo)
        call bcast_convinfo(cxyzn, conv_dmat, conv_rhs, conv_head, max_var)
#else
      cxyzn = get_cnum(max_num)
      conv_dmat = array_var(1)%dmat(max_num)*len_scal**2
      conv_rhs = array_var(1)%rhs(max_num)*len_scal**3
      conv_head = head_new(max_num)*len_scal
#endif
      conv_var = max_var*len_scal
      if (my_rank == 0) then
        write(conv_fnum,11) out_iter, in_iter, back_iter, conv_var, trim(adjustl(cxyzn)),&
                            conv_dmat, conv_rhs, conv_head
      end if
      ! check outer_loop
      if (conv_flag) then
        ! -- Calculate surface water level (surfw)
          call calc_surfw()
        if (st_out_step%rest == DZERO) then
          ! -- Write restart file (rest)
            call write_rest(head_new)
        else if (mod(current_t,st_out_step%rest) == 0) then
          ! -- Write restart file (rest)
            call write_rest(head_new)
        end if
        exit outer_loop
      else if (back_flag .and. st_sim%sim_type /= -1) then
        if (my_rank == 0) then
          write(conv_fnum,13)
        end if
        exit outer_loop
      else if (out_iter == maxout_iter .and. st_sim%sim_type == -1) then
        st_sim%sim_type = 0
        current_t = DZERO
        if (my_rank == 0) then
          write(log_fnum,'(a)') "Steady state calculation changes to timestep calculation."
          write(conv_fnum,15)
        end if
        exit outer_loop
      else if (out_iter == maxout_iter) then
        if (my_rank == 0) then
          write(conv_fnum,14)
        end if
        exit outer_loop
      else if ((abs(conv_var) >= VARMAX .or. max_unk >= XMAX) .and. st_sim%sim_type /= -1) then
        if (my_rank == 0) then
          write(conv_fnum,12)
        end if
        exit outer_loop
      end if

    end do outer_loop

    deallocate(new_func)

  end subroutine calc_numsol

  subroutine reset_matvec
  !*********************************************************************************************
  ! reset_matvec -- Reset coefficients matrix and constant vector
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: precon_type, nlevel
    use set_cell, only: amg_setflag
    use allocate_solution, only: crs_index, pro_var, res_var
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    array_var(1)%lumat(:) = DZERO
    array_var(1)%dmat(:) = DZERO
    array_var(1)%rhs(:) = DZERO

    if (precon_type == 1 .and. amg_setflag == 1) then
      do i = 2, nlevel
        deallocate(array_var(i)%dmat, array_var(i)%rhs, array_var(i)%lumat, array_var(i)%x)
        deallocate(crs_index(i)%offrow, crs_index(i)%offind)
        deallocate(pro_var(i)%pindex, pro_var(i)%poffrow, pro_var(i)%pval)
        deallocate(res_var(i)%rindex, res_var(i)%roffrow, res_var(i)%rval)
      end do
    end if

  end subroutine reset_matvec

  subroutine calc_surfw()
  !*********************************************************************************************
  ! calc_surfw -- Calculate surface water level
  !*********************************************************************************************
    ! -- modules
    use set_cell, only: ncals
!    use make_cell, only: surf_elev
!    use prep_calculation, only: surf_bott
    use allocate_solution, only: surf_head
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, ncals
!      if (head_new(i) <= surf_elev(i)) then
!        surf_head(i) = head_new(i)
!      else
!        surf_head(i) = surf_elev(i)
!      end if
      ! all surface head = surf_elev
!      surf_head(i) = surf_elev(i)
      ! all surface head = surf_bottom
!      surf_head(i) = surf_bott(i)
      ! all surface head = head_new
      surf_head(i) = head_new(i)
    end do
    !$omp end parallel do

  end subroutine calc_surfw

  subroutine set_backtr(maxstep)
  !*********************************************************************************************
  ! set_backtr -- Set for backtracking
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(out) :: maxstep
    ! -- local
    real(DP) :: l2_xnew
#ifdef MPI_MSG
    real(DP) :: sum_l2
#endif
    !-------------------------------------------------------------------------------------------
    ! -- Calculate l2 norm square (resl2norm2)
      call calc_l2norm2(1, head_new, l2_xnew)
#ifdef MPI_MSG
    if (pro_totn /= 1) then
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

  subroutine set_eise_walk(lam, l2_new, l2_pre, l2_jac, grad, eta)
  !*********************************************************************************************
  ! set_eise_walk -- Set Eisenstat-Walker forcing term
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: lam, l2_new, l2_pre, l2_jac, grad
    real(DP), intent(inout) :: eta
    ! -- local
    real(DP), parameter :: ETA_MAX = 0.9_DP
    real(DP), parameter :: ETA_MIN = 1.0E-4_DP
    real(DP), parameter :: ETA_ALPHA = (1.0_DP+sqrt(5.0_DP))*DHALF
    real(DP) :: eta_safe, minus1, l2_line, lin_l2norm
    !-------------------------------------------------------------------------------------------
    eta_safe = eta**ETA_ALPHA
    minus1 = DONE - lam
    l2_line = minus1*minus1*l2_pre + DTWO*lam*minus1*grad + lam*lam*l2_jac
    lin_l2norm = sqrt(max(DZERO, l2_line))
    eta = abs(sqrt(l2_new) - lin_l2norm)/sqrt(l2_pre)

    if (eta_safe < 0.1_DP) then
      eta_safe = DZERO
    end if

    eta = max(eta, eta_safe)
    eta = max(eta, ETA_MIN)
    eta = min(eta, ETA_MAX)

  end subroutine set_eise_walk

  subroutine run_backtr(backi, backf, betai, l2_new, l2_pre, l2_jac, lam, grad, maxstep, new_f)
  !*********************************************************************************************
  ! run_backtr -- Run backtracking
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: MACHI_EPS
    use calc_function, only: calc_vecjacf
#ifdef MPI_MSG
    use mpi_utility, only: mpimax_val
#endif
    ! -- inout
    integer(I4), intent(inout) :: backi, betai
    logical, intent(inout) :: backf
    real(DP), intent(inout) :: l2_new, l2_jac, lam
    real(DP), intent(in) :: l2_pre, maxstep
    real(DP), intent(out) :: grad
    real(DP), intent(inout) :: new_f(:)
    ! -- local
    integer(I4) :: i
    real(DP) :: l2_new2, l2_pnorm, slope, av, bv, rhs1, rhs2, root, step_len, step_tol
    real(DP) :: lam2, temp_lam, lam_inv, lam2_inv, del_lam, lam_max, lam_min
    real(DP) :: lam_length, lam_base, lam_diff, lam_incr, sql2_pnorm, maxpnorm
    real(DP) :: alpha_cond, beta_cond
    real(DP), parameter :: BACK_ALPHA = 1.00E-4_DP
    real(DP), parameter :: BACK_BETA = 0.9_DP
    real(DP), allocatable :: jacvec(:)
#ifdef MPI_MSG
    real(DP) :: sum_l2, max_val
#endif
    !-------------------------------------------------------------------------------------------
    l2_new2 = l2_new ; lam = DONE ; lam2 = DONE ; maxpnorm = DONE ; lam_length = DONE
    allocate(jacvec(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      jacvec(i) = DZERO
    end do
    !$omp end parallel do
    ! -- Calculate function value (func)
      call calc_func(head_new, new_f)
    ! -- Calculate l2 norm square (resl2norm2)
      call calc_l2norm2(1, head_change, l2_pnorm)
    ! -- Calculate vector by jacobi-free (vecjocf)
      call calc_vecjacf(1, head_pre, head_change, jacvec)

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_pnorm, "change function l2-norm", sum_l2)
      l2_pnorm = sum_l2
    end if
#endif

    sql2_pnorm = sqrt(l2_pnorm)
    if (sql2_pnorm > maxstep) then
      maxpnorm = maxstep/sql2_pnorm
      !$omp parallel do private(i)
      do i = 1, ncalc
        head_change(i) = head_change(i)*maxpnorm
      end do
      !$omp end parallel do
      sql2_pnorm = maxstep
    end if
    step_len = sql2_pnorm

    l2_new = DZERO ; slope = DZERO
    !$omp parallel
    !$omp do private(i) reduction(+:l2_new, slope)
    do i = 1, ncalc
      l2_new = l2_new + new_f(i)*new_f(i)
      slope = slope + array_var(1)%rhs(i)*jacvec(i)*maxpnorm
    end do
    !$omp end do

    !$omp do private(i, temp_lam) reduction(max:lam_length)
    do i = 1, ncalc
      if (head_new(i) /= DZERO) then
        temp_lam = abs(head_change(i))/abs(head_new(i))
        if (temp_lam > lam_length) then
          lam_length = temp_lam
        end if
      end if
    end do
    !$omp end do
    !$omp end parallel

    l2_jac = l2_jac*maxpnorm*maxpnorm

    deallocate(jacvec)

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_new, "current function l2-norm", sum_l2)
      l2_new = sum_l2
      ! -- Sum value for MPI (val)
        call mpisum_val(slope, "slope function l2-norm", sum_l2)
      slope = sum_l2
      ! -- MAX value for MPI (val)
        call mpimax_val(lam_length, "lambda length", max_val)
      lam_length = max_val
    end if
#endif
    step_tol = MACHI_EPS**(2.0_DP/3.0_DP)
    lam_min = step_tol/lam_length ; alpha_cond = DHALF*l2_pre + BACK_ALPHA*lam*slope
    grad = slope
    back_aloop: do
      if (DHALF*l2_new <= alpha_cond) then
        exit back_aloop
      end if
      backi = backi + 1
      ! -- Calculate function and l2norm2 (func2norm)
        call calc_funcl2norm(lam, l2_new, new_f)

      if (lam == DONE) then
        temp_lam = -slope/(DTWO*(DHALF*l2_new-DHALF*l2_pre-slope))
      else
        rhs1 = DHALF*l2_new - DHALF*l2_pre - lam*slope
        rhs2 = l2_new2 - DHALF*l2_pre - lam2*slope
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
      l2_new2 = DHALF*l2_new
      lam = max(temp_lam, 0.1_DP*lam)
      if (lam < lam_min) then
        backf = .true.
        return
      end if
      alpha_cond = DHALF*l2_pre + BACK_ALPHA*lam*slope
    end do back_aloop

    alpha_cond = DHALF*l2_pre + BACK_ALPHA*lam*slope
    beta_cond = DHALF*l2_pre + BACK_BETA*lam*slope
    if (DHALF*l2_new < beta_cond) then
      if (lam == DONE .and. sql2_pnorm < step_len) then
        lam_max = step_len/sql2_pnorm
        b1_loop: do
          if (DHALF*l2_new > alpha_cond .or. DHALF*l2_new >= beta_cond .or. &
              lam >= lam_max) then
            exit b1_loop
          end if
          lam2 = lam ; l2_new2 = DHALF*l2_new ; lam = min(DTWO*lam, lam_max)
          ! -- Calculate function and l2norm2 (func2norm)
            call calc_funcl2norm(lam, l2_new, new_f)
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
            call calc_funcl2norm(lam, l2_new, new_f)
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

        if (DHALF*l2_new < beta_cond) then
          ! -- Calculate function and l2norm2 (func2norm)
            call calc_funcl2norm(lam_base, l2_new, new_f)
          betai = betai + 1
        end if
        if (betai == 10) then
          backf = .true.
          return
        end if
      end if
    end if

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

  subroutine calc_funcl2norm(in_lam, l2_new, new_f)
  !*********************************************************************************************
  ! calc_funcl2norm -- Calculate function and l2norm2
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: in_lam
    real(DP), intent(out) :: l2_new
    real(DP), intent(inout) :: new_f(:)
    ! -- local
    integer(I4) :: i
#ifdef MPI_MSG
    real(DP) :: sum_l2
#endif
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, nreg_num
      head_new(i) = head_pre(i) + in_lam*head_change(i)
    end do
    !$omp end parallel do

    ! -- Calculate function value (func)
      call calc_func(head_new, new_f)
    ! -- Calculate l2 norm square (resl2norm2)
      call calc_l2norm2(1, new_f, l2_new)

#ifdef MPI_MSG
    if (pro_totn /= 1) then
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
      call write_mpi_rest(rest_fnum, now_time, len_scal, rest_head)
#else
    rewind(rest_fnum)
    write(rest_fnum) real(now_time, kind=DP)
    write(rest_fnum) (rest_head(i)*len_scal, i = 1, ncalc)
#endif

  end subroutine write_rest

end module nonlinear_solution
