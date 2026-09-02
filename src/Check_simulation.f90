module check_simulation
  ! -- modules
  use kind_module, only: I4, DP
  use initial_module, only: st_sim
  use prep_calculation, only: st_time

  implicit none
  private
  public :: check_insol, check_abserrmax, check_residual, check_stepnorm
  public :: check_outtiming, check_lastts
  integer(I4), public :: write_flag, lasttime_flag

  ! -- local

  contains

  subroutine check_insol(bsum, rsum)
  !*********************************************************************************************
  ! check_insol -- Check inner solution
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_ctrl
    ! -- inout
    real(DP), intent(in) :: bsum, rsum
    ! -- local
    real(DP) :: rerr
    !-------------------------------------------------------------------------------------------
    rerr = sqrt(rsum/bsum)

    if (rerr <= st_ctrl%errtol) then
      st_time%conv_flag = .true.
    else
      st_time%conv_flag = .false.
    end if

  end subroutine check_insol

  subroutine check_abserrmax(x_new, x_pre, chmax, xmax, nunmax)
  !*********************************************************************************************
  ! check_abserrmax -- Check absolute error max norm
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DZERO
    use set_cell, only: ncalc
    ! -- inout
    real(DP), intent(in) :: x_new(:), x_pre(:)
    real(DP), intent(out) :: chmax, xmax
    integer(I4), intent(out) :: nunmax
    ! -- local
    integer(I4) :: i, nun
    integer(I4) :: nan_flag
    real(DP) :: ach, aval, bch, bval, xdif, xdif_abs, x_abs
    real(DP) :: change, absolute
    !-------------------------------------------------------------------------------------------
    nan_flag = 0 ; nun = 0
    aval = DZERO ; bval = DZERO ; bch = DZERO ; ach = DZERO
    change = DZERO ; absolute = DZERO
    do i = 1, ncalc
      if (x_new(i) /= x_new(i)) then
        nan_flag = 1
      end if
      x_abs = abs(x_new(i))
      xdif = x_new(i) - x_pre(i)
      xdif_abs = abs(xdif)
      if (xdif_abs >= ach) then
        bch = xdif
        ach = xdif_abs
        nun = i
      end if
      if (x_abs >= aval) then
        aval = x_abs
        bval = x_new(i)
      end if
    end do

    if (nan_flag == 0) then
      if (ach >= change) then
        chmax = bch
        nunmax = nun
        change = ach
      end if
      if (aval >= absolute) then
        xmax = bval
        absolute = aval
      end if
    else
      chmax = DZERO
      xmax = transfer(-1_8,DZERO)
      nunmax = 0
    end if

  end subroutine check_abserrmax

  subroutine check_residual(funcv, stnew, res_flag)
  !*********************************************************************************************
  ! check_residual -- Check residual convergence criteria
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DZERO
#ifdef MPI_MSG
    use utility_module, only: st_mpi
    use mpi_utility, only: mpisum_val
#endif
    use initial_module, only: st_ctrl
    use read_input, only: len_scal
    use set_cell, only: ncalc
    use make_cell, only: st_geom
    ! -- inout
    real(DP), intent(in) :: funcv(:), stnew(:)
    logical, intent(out) :: res_flag
    ! -- local
    integer(I4) :: i, viol_num
    real(DP), parameter :: ACC_FLOOR = 1.00E-15_DP
    real(DP) :: vol_scal, res_val, acc_val
    logical :: abs_flag, rel_flag, pass_flag
#ifdef MPI_MSG
    integer(I4) :: sum_viol
#endif
    !-------------------------------------------------------------------------------------------
    res_flag = .true.
    abs_flag = st_ctrl%res_abs_tol > DZERO ; rel_flag = st_ctrl%res_rel_tol > DZERO
    if (.not. abs_flag .and. .not. rel_flag) then
      return
    end if

    viol_num = 0 ; vol_scal = len_scal**3
    !$omp parallel do private(i, res_val, acc_val, pass_flag) reduction(+:viol_num)
    do i = 1, ncalc
      res_val = abs(funcv(i))*vol_scal
      acc_val = abs(stnew(i))*st_time%delt_inv*st_geom%cell_vol(i)*vol_scal
      pass_flag = .false.
      if (abs_flag .and. res_val <= st_ctrl%res_abs_tol) then
        pass_flag = .true.
      else if (rel_flag .and. acc_val <= ACC_FLOOR) then
        pass_flag = .true.
      else if (rel_flag .and. res_val <= st_ctrl%res_rel_tol*acc_val) then
        pass_flag = .true.
      end if
      if (.not. pass_flag) then
        viol_num = viol_num + 1
      end if
    end do
    !$omp end parallel do

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(viol_num, "residual criteria violation", sum_viol)
      viol_num = sum_viol
    end if
#endif
    if (viol_num /= 0) then
      res_flag = .false.
    end if

  end subroutine check_residual

  subroutine check_stepnorm(x_new, x_pre, stepmax)
  !*********************************************************************************************
  ! check_stepnorm -- Check scaled step max norm
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DZERO, DONE
#ifdef MPI_MSG
    use utility_module, only: st_mpi
    use mpi_utility, only: mpimax_val
#endif
    use read_input, only: len_scal
    use set_cell, only: ncalc
    ! -- inout
    real(DP), intent(in) :: x_new(:), x_pre(:)
    real(DP), intent(out) :: stepmax
    ! -- local
    integer(I4) :: i
    real(DP) :: new_val, chg_val, stp_val
#ifdef MPI_MSG
    real(DP) :: max_step
#endif
    !-------------------------------------------------------------------------------------------
    stepmax = DZERO
    !$omp parallel do private(i, new_val, chg_val, stp_val) reduction(max:stepmax)
    do i = 1, ncalc
      new_val = abs(x_new(i))*len_scal
      chg_val = abs(x_new(i) - x_pre(i))*len_scal
      stp_val = chg_val/(DONE + new_val)
      stepmax = max(stepmax, stp_val)
    end do
    !$omp end parallel do

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Max value for MPI (val)
        call mpimax_val(stepmax, "scaled step max norm", max_step)
      stepmax = max_step
    end if
#endif

  end subroutine check_stepnorm

  subroutine check_outtiming()
  !*********************************************************************************************
  ! check_outtiming -- check output timing
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_step_flag
    ! -- inout

    ! -- local
    integer(I4) :: step_flag
    !-------------------------------------------------------------------------------------------
    step_flag = st_step_flag%rech + st_step_flag%well + st_step_flag%seal +&
                st_step_flag%prec + st_step_flag%evap + st_step_flag%riwl +&
                st_step_flag%riwd + st_step_flag%ribl + st_step_flag%ride +&
                st_step_flag%riwi + st_step_flag%lawl + st_step_flag%lawd +&
                st_step_flag%labl + st_step_flag%laar

    if (step_flag > 0) then
      write_flag = 1
    else if (st_time%current_t >= st_sim%end_time) then
      write_flag = 1
    else if (st_time%conv_flag .and. st_sim%sim_type /= 1) then
      write_flag = 1
    else
      write_flag = 0
    end if

  end subroutine check_outtiming

  subroutine check_lastts()
  !*********************************************************************************************
  ! check_lastts -- Check last time step
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------------
    if (st_time%current_t >= st_sim%end_time) then
      lasttime_flag = 1
    else if (st_time%conv_flag .and. st_sim%sim_type == -1) then
      lasttime_flag = 1
    else
      lasttime_flag = 0
    end if

  end subroutine check_lastts

end module check_simulation
