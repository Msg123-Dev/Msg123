module check_simulation
  ! -- modules
  use kind_module, only: I4, DP
  use initial_module, only: st_sim
  use prep_calculation, only: current_t, conv_flag

  implicit none
  private
  public :: check_insol, check_abserrmax, check_outtiming, check_lastts
  integer(I4), public :: write_flag, lasttime_flag

  ! -- local

  contains

  subroutine check_insol(bsum, rsum)
  !***************************************************************************************
  ! check_insol -- Check inner solution
  !***************************************************************************************
    ! -- modules
    use initial_module, only: errtol
    ! -- inout
    real(DP), intent(in) :: bsum, rsum
    ! -- local
    real(DP) :: rerr
    !-------------------------------------------------------------------------------------
    rerr = sqrt(rsum/bsum)

    if (rerr <= errtol) then
      conv_flag = 1
    else
      conv_flag = 0
    end if

  end subroutine check_insol

  subroutine check_abserrmax(x_new, x_pre, chmax, xmax, nunmax)
  !***************************************************************************************
  ! check_abserrmax -- Check absolute error max norm
  !***************************************************************************************
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
    !-------------------------------------------------------------------------------------
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

  subroutine check_outtiming()
  !***************************************************************************************
  ! check_outtiming -- check output timing
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_step_flag
    ! -- inout

    ! -- local
    integer(I4) :: step_flag
    !-------------------------------------------------------------------------------------
    step_flag = st_step_flag%rech + st_step_flag%well + st_step_flag%seal +&
                st_step_flag%prec + st_step_flag%evap + st_step_flag%riwl +&
                st_step_flag%riwd + st_step_flag%ribl + st_step_flag%ride +&
                st_step_flag%riwi + st_step_flag%lawl + st_step_flag%lawd +&
                st_step_flag%labl + st_step_flag%laar

    if (step_flag > 0) then
      write_flag = 1
    else if (current_t >= st_sim%end_time) then
      write_flag = 1
    else if (conv_flag == 1 .and. st_sim%sim_type /= 1) then
      write_flag = 1
    else
      write_flag = 0
    end if

  end subroutine check_outtiming

  subroutine check_lastts()
  !***************************************************************************************
  ! check_lastts -- Check last time step
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    if (current_t >= st_sim%end_time) then
      lasttime_flag = 1
    else if (conv_flag == 1 .and. st_sim%sim_type == -1) then
      lasttime_flag = 1
    else
      lasttime_flag = 0
    end if

  end subroutine check_lastts

end module check_simulation
