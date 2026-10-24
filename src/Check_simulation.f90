module check_simulation
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, UROUND
  use utility_module, only: st_mpi
  use initial_module, only: st_sim, st_ctrl
  use read_input, only: len_scal
  use set_cell, only: ncalc
  use make_cell, only: st_geom
  use prep_calculation, only: st_time
  use allocate_solution, only: crs_index, array_var

  implicit none
  private
  public :: check_insol, check_abserrmax, check_residual, calc_res_weight
  public :: check_outtiming, check_lastts
  integer(I4), public :: write_flag, lasttime_flag

  ! -- local
  real(DP), parameter :: ROUND_FACT = 4.00E+00_DP

  contains

  subroutine check_insol(bsum, rsum)
  !*********************************************************************************************
  ! check_insol -- Check inner solution
  !*********************************************************************************************
    ! -- modules
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

  subroutine check_residual(funcv, fscal, res_flag, res_fact, stnew, stold, hnew)
  !*********************************************************************************************
  ! check_residual -- Check residual convergence criteria
  !*********************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_utility, only: mpisum_val
#endif
    ! -- inout
    real(DP), intent(in) :: funcv(:), fscal(:), res_fact
    logical, intent(out) :: res_flag
    real(DP), intent(in), optional :: stnew(:), stold(:), hnew(:)
    ! -- local
    integer(I4) :: i, k, viol_num
    real(DP) :: vol_scal, res_val, scal_val, abs_tol, rel_tol, round_val, acc_val
    logical :: abs_flag, rel_flag, pass_flag, round_flag, stor_flag
#ifdef MPI_MSG
    integer(I4) :: sum_viol
#endif
    !-------------------------------------------------------------------------------------------
    res_flag = .true.
    abs_flag = st_ctrl%res_abs_tol > DZERO ; rel_flag = st_ctrl%res_rel_tol > DZERO
    if (.not. abs_flag .and. .not. rel_flag) then
      return
    end if
    round_flag = present(hnew) ; stor_flag = present(stnew) .and. present(stold)

    viol_num = 0 ; vol_scal = len_scal**3 ; acc_val = calc_acc_floor(fscal)
    abs_tol = st_ctrl%res_abs_tol*res_fact ; rel_tol = st_ctrl%res_rel_tol*res_fact
    !$omp parallel do private(i, k, res_val, scal_val, round_val, pass_flag) &
    !$omp             reduction(+:viol_num)
    do i = 1, ncalc
      res_val = abs(funcv(i))*vol_scal
      scal_val = fscal(i)*vol_scal
      pass_flag = .false.
      if (abs_flag .and. res_val <= abs_tol) then
        pass_flag = .true.
      else if (rel_flag .and. scal_val <= acc_val) then
        pass_flag = .true.
      else if (rel_flag .and. res_val <= rel_tol*scal_val) then
        pass_flag = .true.
      else if (round_flag) then
        round_val = abs(array_var(1)%dmat(i))*abs(hnew(i))
        if (stor_flag) then
          round_val = max(abs(stnew(i)), abs(stold(i)))*st_geom%cell_vol(i)*st_time%delt_inv&
                      + round_val
        end if
        do k = crs_index(1)%offind(i-1)+1, crs_index(1)%offind(i)
          round_val = round_val + abs(array_var(1)%lumat(k))*abs(hnew(crs_index(1)%offrow(k)))
        end do
        pass_flag = res_val <= ROUND_FACT*UROUND*round_val*vol_scal
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

  subroutine calc_res_weight(fscal, stnew, stold, hnew, res_wgt)
  !*********************************************************************************************
  ! calc_res_weight -- Calculate residual weight from the residual convergence criteria
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: fscal(:), stnew(:), stold(:), hnew(:)
    real(DP), intent(out) :: res_wgt(:)
    ! -- local
    integer(I4) :: i, k
    real(DP) :: vol_scal, scal_val, abs_tol, rel_tol, round_val, allow_val, acc_val
    logical :: abs_flag, rel_flag
    !-------------------------------------------------------------------------------------------
    abs_flag = st_ctrl%res_abs_tol > DZERO ; rel_flag = st_ctrl%res_rel_tol > DZERO
    vol_scal = len_scal**3 ; acc_val = calc_acc_floor(fscal)
    abs_tol = st_ctrl%res_abs_tol ; rel_tol = st_ctrl%res_rel_tol
    !$omp parallel do private(i, k, scal_val, round_val, allow_val)
    do i = 1, ncalc
      scal_val = fscal(i)*vol_scal
      round_val = max(abs(stnew(i)), abs(stold(i)))*st_geom%cell_vol(i)*st_time%delt_inv&
                  + abs(array_var(1)%dmat(i))*abs(hnew(i))
      do k = crs_index(1)%offind(i-1)+1, crs_index(1)%offind(i)
        round_val = round_val + abs(array_var(1)%lumat(k))*abs(hnew(crs_index(1)%offrow(k)))
      end do
      allow_val = ROUND_FACT*UROUND*round_val*vol_scal
      if (abs_flag) then
        allow_val = allow_val + abs_tol
      end if
      if (rel_flag) then
        allow_val = allow_val + rel_tol*scal_val
      end if
      if (rel_flag .and. scal_val <= acc_val) then
        res_wgt(i) = DZERO
      else if (allow_val > DZERO) then
        res_wgt(i) = vol_scal/allow_val
      else
        res_wgt(i) = DZERO
      end if
    end do
    !$omp end parallel do

  end subroutine calc_res_weight

  function calc_acc_floor(fscal) result(acc_val)
  !*********************************************************************************************
  ! calc_acc_floor -- Calculate the floor of the residual scale relative to the largest cell
  !*********************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_utility, only: mpimax_val
#endif
    ! -- inout
    real(DP), intent(in) :: fscal(:)
    ! -- local
    integer(I4) :: i
    real(DP), parameter :: ACC_REL = 1.00E-09_DP
    real(DP) :: smax
    real(DP) :: acc_val
#ifdef MPI_MSG
    real(DP) :: all_max
#endif
    !-------------------------------------------------------------------------------------------
    smax = DZERO
    !$omp parallel do private(i) reduction(max:smax)
    do i = 1, ncalc
      smax = max(smax, fscal(i))
    end do
    !$omp end parallel do
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Max value for MPI (val)
        call mpimax_val(smax, "residual scale maximum", all_max)
      smax = all_max
    end if
#endif
    acc_val = ACC_REL*smax*len_scal**3

  end function calc_acc_floor

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
    step_flag = st_step_flag%rech + st_step_flag%well + st_step_flag%seal + st_step_flag%prec +&
                st_step_flag%evap + st_step_flag%riwl + st_step_flag%riwd + st_step_flag%ribl +&
                st_step_flag%ride + st_step_flag%riwi + st_step_flag%rile + st_step_flag%ribk +&
                st_step_flag%ribt + st_step_flag%lawl + st_step_flag%lawd + st_step_flag%labl +&
                st_step_flag%laar + st_step_flag%labk + st_step_flag%labt + st_step_flag%rfra

    if (step_flag > 0) then
      write_flag = 1
    else if (st_time%current_t >= st_sim%end_time) then
      write_flag = 1
    else if (st_time%conv_flag .and. (st_sim%sim_type /= 1 .or. st_ctrl%outstep_type == 1)) then
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
