module nonlinear_solution
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DHALF, DONE, DTWO, VARMAX
  use initial_module, only: criteria
  use read_input, only: len_scal
  use check_condition, only: st_out_fnum
  use set_cell, only: ncalc
  use prep_calculation, only: out_iter
  use allocate_solution, only: array_var, head_new, head_pre, head_change
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

  ! -- local
  integer(I4) :: back_iter, back_flag, beta_iter
  real(DP) :: maxvar
  real(DP) :: l2_new, l2_pre, l2_jnorm, lam, eta, gradient, maxstep
  real(DP), allocatable :: new_f(:)
#ifdef MPI_MSG
  real(DP) :: sum_l2
#endif
  contains

  subroutine calc_numsol()
  !***************************************************************************************
  ! calc_numsol -- Calculate numerical solution
  !***************************************************************************************
    ! -- modules
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
    integer(I4) :: maxnun, conv_fnum
    character(VARLEN) :: cxyzn
    real(DP) :: maxunk, check_val
    real(DP) :: conv_dmat, conv_rhs, conv_head, conv_var
#ifdef MPI_MSG
  real(DP) :: var_max, unk_max, var_abs_max
#endif
    ! -- format
    10 format(//1x,"CURRENT TIME : ",es12.5,1x,"(",a,")",20x,"TIME STEP : ",&
              es12.5,1x,"(SEC)",/,1x,84("-"),/,1x,&
              " OUTER INNER  BACK     MAXIMUM           MAXIMUM    DIAGONAL  RIGHT HAND     UNKNOWN",/,1x,&
              "                        CHANGE              CELL      MATRIX      VECTOR       VALUE",/,1x,84("-"))
    11 format(1X,3(i6),es12.3,a18,4(es12.3))
    12 format(1X,"Didn't converge due to maximum value or change")
    13 format(1X,"Stop due to maximum value or change in backtracking")
    14 format(1X,"Stop due to maximum number of nonlinear iteration")
    15 format(1X,"Didn't converge in steady state calculation")
    !-------------------------------------------------------------------------------------
    conv_fnum = st_out_fnum%conv
    ! -- Set for backtracking (backtr)
      call set_backtr()

    outer_loop : do out_iter = 1, maxout_iter

      if (out_iter == 1) then
        if (my_rank == 0) then
          write(conv_fnum,10) now_time, trim(st_sim%cal_unit), delt
        end if
        form_switch = 1
        allocate(new_f(ncalc))
        !$omp parallel workshare
        new_f(:) = DZERO
        !$omp end parallel workshare
      else
        l2_pre = l2_new
      end if

      back_flag = 0 ; back_iter = 0 ; beta_iter = 0

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
          call calc_l2norm2(1, array_var(1)%rhs, l2_new)
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Sum value for MPI (val)
            call mpisum_val(l2_new, "initial function l2-norm", sum_l2)
          l2_new = sum_l2
        end if
#endif
        l2_pre = l2_new ; eta = DHALF
      end if

      !$omp parallel workshare
      head_pre(:) = head_new(:) ; head_change(:) = DZERO
      !$omp end parallel workshare

      conv_flag = 0
      if (st_sim%sim_type /= -1) then
        errtol = eta
      else
        errtol = XMAX_INV**3
      end if
      ! -- Solve linear algebra (linalg)
        call solve_linalg(l2_pre, head_change, l2_jnorm)

      !$omp parallel workshare
      head_new(:) = head_pre(:) + head_change(:)
      !$omp end parallel workshare

      ! -- Check absolute error max norm
        call check_abserrmax(head_new, head_pre, maxvar, maxunk, maxnun)

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Check mpi max error (mpimaxerr)
          call check_mpimaxerr(maxvar, maxunk, var_abs_max, unk_max, var_max)
        check_val = var_abs_max*len_scal ; maxunk = unk_max*len_scal
      else
        check_val = abs(maxvar)*len_scal ; maxunk = abs(maxunk)*len_scal
        var_max = maxvar
      end if
#else
      check_val = abs(maxvar)*len_scal ; maxunk = abs(maxunk)*len_scal
#endif
      if (conv_flag /= 0 .and. maxunk == maxunk) then
        conv_flag = 0
        if (check_val <= criteria .and. maxunk < XMAX) then
          conv_flag = 1
        end if
      else if (maxunk /= maxunk .and. st_sim%sim_type /= -1) then
        if (my_rank == 0) then
          write(log_fnum,'(a)') "Nan detected."
        end if
        exit outer_loop
      else if (maxunk /= maxunk .and. st_sim%sim_type == -1) then
        st_sim%sim_type = 0
        current_t = DZERO
        if (my_rank == 0) then
          write(log_fnum,'(a)') "Steady state calculation changes to timestep calculation."
          write(conv_fnum,15)
        end if
        exit outer_loop
      end if

      if (conv_flag /= 1 .and. form_switch == 1) then
        ! -- Run backtracking (backtr)
          call run_backtr()
        ! -- Check absolute error max norm
          call check_abserrmax(head_new, head_pre, maxvar, maxunk, maxnun)
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Check mpi max error (mpimaxerr)
            call check_mpimaxerr(maxvar, maxunk, var_abs_max, unk_max, var_max)
          check_val = var_abs_max*len_scal ; maxunk = unk_max*len_scal
        else
          check_val = abs(maxvar)*len_scal ; maxunk = abs(maxunk)*len_scal
          var_max = maxvar
        end if
#else
        check_val = abs(maxvar)*len_scal ; maxunk = abs(maxunk)*len_scal
#endif
        if ((check_val >= VARMAX .or. maxunk >= XMAX) .and. st_sim%sim_type /= -1) then
          back_flag = 1
        else if (check_val <= criteria .and. maxunk < XMAX) then
          conv_flag = 1
        else if (st_sim%sim_type == -1) then
          back_flag = 0
        end if
        if (conv_flag /= 1 .and. back_flag /= 1) then
          ! -- Set eater (eater)
            call set_eater()
        end if
      end if

!      if (conv_flag /= 1 .and. form_switch == 1 .and. back_flag /= 1) then
!        ! -- Run under relax using cooley underrelaxation (relaxcooly)
!          call run_relax_cooly(head_new, head_pre)
!        ! -- Run under relax using delta-bar-delta underrelaxation (relaxdelta)
!          call run_relax_delta(head_new, head_pre)

!        ! -- Check absolute error max norm
!          call check_abserrmax(head_new, head_pre, maxvar, maxunk, maxnun)
!#ifdef MPI_MSG
!        if (pro_totn /= 1) then
!          ! -- Check mpi max error (mpimaxerr)
!            call check_mpimaxerr(maxvar, maxunk, var_abs_max, unk_max, var_max)
!          check_val = var_abs_max*len_scal ; maxunk = unk_max*len_scal
!        else
!          check_val = abs(maxvar)*len_scal ; maxunk = abs(maxunk)*len_scal
!        end if
!#else
!        check_val = abs(maxvar)*len_scal ; maxunk = abs(maxunk)*len_scal
!#endif
!        if (check_val <= criteria .and. maxunk < XMAX) then
!          conv_flag = 1
!        end if
!      end if

      ! -- Calculate function value (func)
        call calc_func(head_new, new_f)
      ! -- Calculate l2 norm square (resl2norm2)
        call calc_l2norm2(1, new_f, l2_new)
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(l2_new, "new function l2-norm", sum_l2)
        l2_new = sum_l2
      end if
      if (maxvar == var_max) then
        ! -- Get cell number (cnum)
          call get_cnum(maxnun, cxyzn)
      else
        cxyzn = ""
      end if
      conv_dmat = array_var(1)%dmat(maxnun)*len_scal**2
      conv_rhs = array_var(1)%rhs(maxnun)*len_scal**3
      conv_head = head_new(maxnun)*len_scal
      ! -- Bcast converge information (convinfo)
        call bcast_convinfo(cxyzn, conv_dmat, conv_rhs, conv_head, maxvar)
#else
      ! -- Get cell number (cnum)
        call get_cnum(maxnun, cxyzn)
      conv_dmat = array_var(1)%dmat(maxnun)*len_scal**2
      conv_rhs = array_var(1)%rhs(maxnun)*len_scal**3
      conv_head = head_new(maxnun)*len_scal
#endif
      conv_var = maxvar*len_scal
      if (my_rank == 0) then
        write(conv_fnum,11) out_iter, in_iter, back_iter, conv_var, trim(adjustl(cxyzn)),&
                            conv_dmat, conv_rhs, conv_head
      end if
      ! check outer_loop
      if (conv_flag == 1) then
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
      else if (back_flag == 1 .and. st_sim%sim_type /= -1) then
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
      else if ((abs(conv_var) >= VARMAX .or. maxunk >= XMAX) .and. st_sim%sim_type /= -1) then
        if (my_rank == 0) then
          write(conv_fnum,12)
        end if
        exit outer_loop
      end if

    end do outer_loop

    deallocate(new_f)

  end subroutine calc_numsol

  subroutine reset_matvec
  !***************************************************************************************
  ! reset_matvec -- Reset coefficients matrix and constant vector
  !***************************************************************************************
    ! -- modules
    use initial_module, only: precon_type, nlevel
    use set_cell, only: amg_setflag
    use allocate_solution, only: crs_index, pro_var, res_var
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel workshare
    array_var(1)%lumat(:) = DZERO
    array_var(1)%dmat(:) = DZERO
    array_var(1)%rhs(:) = DZERO
    !$omp end parallel workshare

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
  !***************************************************************************************
  ! calc_surfw -- Calculate surface water level
  !***************************************************************************************
    ! -- modules
    use set_cell, only: ncals
!    use make_cell, only: surf_elev
!    use set_condition, only: surf_area
!    use prep_calculation, only: surf_bott
    use allocate_solution, only: surf_head
!    use calc_function, only: calc_mass
    ! -- inout

    ! -- local
    integer(I4) :: i
!    real(DP), allocatable :: sst_ms(:), sco_ms(:), sse_ms(:), swe_ms(:)
!    real(DP), allocatable :: sre_ms(:), ssu_ms(:), sri_ms(:), sla_ms(:)
!    real(DP) :: surf_ms
    !-------------------------------------------------------------------------------------
!    allocate(sst_ms(ncalc), sco_ms(ncalc), sse_ms(ncalc), swe_ms(ncalc))
!    allocate(sre_ms(ncals), ssu_ms(ncals), sri_ms(ncals), sla_ms(ncals))
!    !$omp parallel workshare
!    sst_ms(:) = DZERO ; sco_ms(:) = DZERO ; sse_ms(:) = DZERO ; swe_ms(:) = DZERO
!    sre_ms(:) = DZERO ; ssu_ms(:) = DZERO ; sri_ms(:) = DZERO ; sla_ms(:) = DZERO
!    !$omp end parallel workshare
!    ! -- Calculate massbalance (mass)
!      call calc_mass(1, head_new, head_pre, sst_ms, sco_ms, sse_ms, swe_ms, sre_ms, ssu_ms, sri_ms, sla_ms)

    !$omp parallel do private(i)
!    !$omp parallel do private(i, surf_ms)
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
      ! all surface head = calculated from massbalance
!      if (seal_cflag(i) /= 1) then
!        surf_ms = sre_ms(i) + swe_ms(i) + sst_ms(i) + sco_ms(i) + ssu_ms(i) + sri_ms(i) + sla_ms(i)
!        if (surf_area(i) /= DZERO) then
!          surf_head(i) = head_new(i) - surf_ms/surf_area(i)
!        else
!          surf_head(i) = head_new(i)
!        end if
!      else
!        surf_head(i) = head_new(i)
!      end if
    end do
    !$omp end parallel do

!    deallocate(sst_ms, sco_ms, sre_ms, swe_ms, sre_ms, ssu_ms, sri_ms, sla_ms)

  end subroutine calc_surfw

  subroutine set_backtr
  !***************************************************************************************
  ! set_backtr -- Set for backtracking
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    real(DP) :: l2_xnew
    !-------------------------------------------------------------------------------------
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

  subroutine set_eater
  !***************************************************************************************
  ! set_eater -- Set eater
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    real(DP) :: eta_max = 0.9_DP, eta_min = 1.0E-4_DP, eta_safe, minus1
    real(DP) :: eta_alpha = (DONE+sqrt(5.0_DP))/DTWO, lin_l2norm
    !-------------------------------------------------------------------------------------
    eta_safe = eta**eta_alpha
    minus1 = DONE - lam
    lin_l2norm = sqrt(minus1*minus1*l2_pre + DTWO*lam*minus1*gradient + lam*lam*l2_jnorm)
    eta = abs(sqrt(l2_new) - lin_l2norm)/sqrt(l2_pre)

    if (eta_safe < 0.1_DP) then
      eta_safe = DZERO
    end if

    eta = max(eta, eta_safe)
    eta = max(eta, eta_min)
    eta = min(eta, eta_max)

  end subroutine set_eater

  subroutine run_backtr
  !***************************************************************************************
  ! run_backtr -- Run backtracking
  !***************************************************************************************
    ! -- modules
    use constval_module, only: MACHI_EPS
    use calc_function, only: calc_vecjacf
#ifdef MPI_MSG
    use mpi_utility, only: mpimax_val
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i
    real(DP) :: l2_new2, l2_pnorm, slope, av, bv, rhs1, rhs2, root, step_len, step_tol
    real(DP) :: lam2, temp_lam, lam_inv, lam2_inv, del_lam, lam_max, lam_min
    real(DP) :: lam_length, lam_base, lam_diff, lam_incr, sql2_pnorm, maxpnorm
    real(DP) :: alpha_cond, beta_cond
    real(DP) :: back_alpha = 1.00E-4_DP, back_beta = 0.9_DP
    real(DP), allocatable :: jacvec(:)
#ifdef MPI_MSG
    real(DP) :: max_val
#endif
    !-------------------------------------------------------------------------------------
    l2_new2 = l2_new ; lam = DONE ; lam2 = DONE ; maxpnorm = DONE ; lam_length = DONE
    allocate(jacvec(ncalc))
    !$omp parallel workshare
    jacvec(:) = DZERO
    !$omp end parallel workshare
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

    !$omp parallel
    !$omp workshare
    l2_new = dot_product(new_f(1:ncalc), new_f(1:ncalc))
    slope = dot_product(array_var(1)%rhs(1:ncalc), jacvec(1:ncalc))*maxpnorm
    !$omp end workshare

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

    l2_jnorm = l2_jnorm*maxpnorm*maxpnorm

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
    lam_min = step_tol/lam_length ; alpha_cond = DHALF*l2_pre + back_alpha*lam*slope
    gradient = slope
    back_aloop: do
      if (DHALF*l2_new <= alpha_cond) then
        exit back_aloop
      end if
      back_iter = back_iter + 1
      ! -- Calculate function and l2norm2 (func2norm)
        call calc_funcl2norm(lam)

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
        back_flag = 1
        return
      end if
      alpha_cond = DHALF*l2_pre + back_alpha*lam*slope
    end do back_aloop

    alpha_cond = DHALF*l2_pre + back_alpha*lam*slope
    beta_cond = DHALF*l2_pre + back_beta*lam*slope
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
            call calc_funcl2norm(lam)
          alpha_cond = DHALF*l2_pre + back_alpha*lam*slope
          beta_cond = DHALF*l2_pre + back_beta*lam*slope
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
            call calc_funcl2norm(lam)
          alpha_cond = DHALF*l2_pre + back_alpha*lam*slope
          beta_cond = DHALF*l2_pre + back_beta*lam*slope

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
            call calc_funcl2norm(lam_base)
          beta_iter = beta_iter + 1
        end if
        if (beta_iter == 10) then
          back_flag = 1
          return
        end if
      end if
    end if

  end subroutine run_backtr

!  subroutine run_relax_cooly(x_new, x_pre)
!  !***************************************************************************************
!  ! run_relax_cooly -- Run under relax using cooley underrelaxation
!  !***************************************************************************************
!    ! -- modules
!
!    ! -- inout
!    real(DP), intent(inout) :: x_new(:), x_pre(:)
!    ! -- local
!    integer(I4) :: i
!    real(DP), save :: old_relax, old_maxvar
!    real(DP) :: relax, esk, esk_abs, del_x
!    !-------------------------------------------------------------------------------------
!    if (out_iter == 1) then
!      relax = DONE
!      old_relax = DONE
!      old_maxvar = maxvar
!    else
!      esk = maxvar/(old_maxvar*old_relax)
!      esk_abs = abs(esk)
!      if (esk < -DONE) then
!        relax = DHALF/esk_abs
!      else
!        relax = (DONE*3 + esk)/(DONE*3 + esk_abs)
!      end if
!    end if
!
!    old_relax = relax
!    old_maxvar = maxvar
!
!    if (relax < DONE) then
!      !$omp parallel do private(i, del_x)
!      do i = 1, ncalc
!        del_x = x_new(i) - x_pre(i)
!        x_new(i) = x_pre(i) + relax*del_x
!      end do
!      !$omp end parallel do
!    end if
!
!  end subroutine run_relax_cooly

!  subroutine run_relax_delta(x_new, x_pre)
!  !***************************************************************************************
!  ! run_relax_delta -- Run under relax using delta-bar-delta under-relaxation
!  !***************************************************************************************
!    ! -- modules
!
!    ! -- inout
!    real(DP), intent(inout) :: x_new(:), x_pre(:)
!    ! -- local
!    integer(I4) :: i
!    real(DP), allocatable, save :: old_w(:), weight_x(:), old_del_x(:)
!    real(DP) :: w, mom, theta = 0.97_DP, kappa = 1.0E-5_DP
!    real(DP) :: gamma = DZERO, momentum = 0.1_DP, del_x
!    !-------------------------------------------------------------------------------------
!    if (.not. allocated(old_w)) then
!      allocate(old_w(ncalc), weight_x(ncalc), old_del_x(ncalc))
!      !$omp parallel workshare
!      old_w(:) = DZERO ; weight_x(:) = DZERO ; old_del_x(:) = DZERO
!      !$omp end parallel workshare
!    end if
!
!    !$omp parallel do private(i, w, mom, del_x)
!    do i = 1, ncalc
!      del_x = x_new(i) - x_pre(i)
!      if (out_iter == 1) then
!        old_w(i) = DONE
!        old_del_x(i) = DZERO
!        weight_x(i) = DZERO
!      end if
!
!      w = old_w(i)
!
!      if (old_del_x(i)*del_x < DZERO) then
!        w = theta*old_w(i)
!      else
!        w = old_w(i) + kappa
!      end if
!
!      if (w > DONE) then
!        w = DONE
!      end if
!
!      old_w(i) = w
!
!      if (out_iter == 1) then
!        weight_x(i) = del_x
!      else
!        weight_x(i) = (DONE - gamma)*del_x + gamma*weight_x(i)
!      end if
!
!      old_del_x(i) = del_x
!
!      mom = DZERO
!      if (out_iter > 4) then
!        mom = momentum
!      end if
!      del_x = del_x*w + mom*weight_x(i)
!      x_new(i) = x_pre(i) + del_x
!    end do
!    !$omp end parallel do
!
!  end subroutine run_relax_delta

  subroutine get_cnum(calc_num, char_cell)
  !***************************************************************************************
  ! get_cnum -- Get cell number
  !***************************************************************************************
    ! -- modules
    use utility_module, only: conv_i2s
    use set_cell, only: get_calc_grid
    ! -- inout
    integer(I4), intent(in) :: calc_num
    character(*), intent(out) :: char_cell
    ! -- local
    integer(I4) :: i_num, j_num, k_num
    character(:), allocatable :: cx_num, cy_num, cz_num
    !-------------------------------------------------------------------------------------
    ! -- Get calculation number from grid number (calc_grid)
      call get_calc_grid(calc_num, i_num, j_num, k_num)

    cx_num = conv_i2s(i_num) ; cy_num = conv_i2s(j_num) ; cz_num = conv_i2s(k_num)

    char_cell = "("//cx_num//","//cy_num//","//cz_num//")"

  end subroutine get_cnum

  subroutine calc_funcl2norm(in_lam)
  !***************************************************************************************
  ! calc_funcl2norm -- Calculate function and l2norm2
  !***************************************************************************************
    ! -- modules
    use allocate_solution, only: nreg_num
    ! -- inout
    real(DP), intent(in) :: in_lam
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
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
  !***************************************************************************************
  ! write_rest -- Write restart file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_write, only: write_mpi_rest
#endif
    ! -- inout
    real(DP), intent(in) :: rest_head(:)
    ! -- local
    integer(I4) :: i, rest_fnum
    !-------------------------------------------------------------------------------------
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
