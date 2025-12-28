module calc_boundary
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SZERO, DZERO
  use set_cell, only: ncals
  use make_cell, only: surf_elev
  use set_condition, only: set_bound2calc
  use assign_boundary, only: read_rech, rech_cflag

  implicit none
  private
  public :: calc_reprev, conv_rech2calc, count_rivecalc, count_lakecalc
  public :: calc_wlbd, calc_blld, calc_wlsl, calc_blsl, calc_lsurf, calc_rivea
  integer(I4), allocatable, public :: rech2cals(:), rive2cals(:), lake2cals(:)
  real(DP), allocatable, public :: calc_rech(:)
  real(DP), allocatable, public :: rive_head(:), rive_bott(:), rive_area(:)
  real(DP), allocatable, public :: lake_head(:), lake_bott(:), lake_area(:)

  ! -- local

  contains

  subroutine calc_reprev(rec_num)
  !***************************************************************************************
  ! calc_reprev -- Calculate recharge from precipitation and evapotranspiration
  !***************************************************************************************
    ! -- modules
    use constval_module, only: SNOVAL
    use initial_module, only: st_prec, st_evap
    use read_input, only: len_scal_inv
    use assign_boundary, only: read_prec, read_evap
    ! -- inout
    integer(I4), intent(inout) :: rec_num
    ! -- local
    integer(I4) :: i
    real(SP) :: no_val, prec_noval, evap_noval
    !-------------------------------------------------------------------------------------
    rec_num = 0 ; no_val = SNOVAL*len_scal_inv
    prec_noval = no_val*st_prec%uni_conv
    evap_noval = no_val*st_evap%uni_conv
    allocate(rech_cflag(ncals), read_rech(ncals))
    !$omp parallel
    !$omp workshare
    rech_cflag(:) = 0 ; read_rech(:) = SZERO
    !$omp end workshare

    !$omp do private(i) reduction(+:rec_num)
    do i = 1, ncals
      if (read_prec(i) > prec_noval .and. read_evap(i) > evap_noval) then
        read_rech(i) = read_prec(i) - read_evap(i)
        rech_cflag(i) = 1 ; rec_num = rec_num + 1
      end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine calc_reprev

  subroutine conv_rech2calc(rec_cnum)
  !***************************************************************************************
  ! conv_rech2calc -- Convert recharge to calculation
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: rec_cnum
    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(rech2cals(rec_cnum), calc_rech(rec_cnum))
    !$omp parallel workshare
    rech2cals(:) = 0 ; calc_rech(:) = DZERO
    !$omp end parallel workshare
    call set_bound2calc(ncals, rech_cflag, read_rech, rech2cals, calc_rech)
    deallocate(read_rech, rech_cflag)
    allocate(read_rech(rec_cnum))
    !$omp parallel workshare
    read_rech(:) = real(calc_rech(:), kind=SP)
    !$omp end parallel workshare

  end subroutine conv_rech2calc

  subroutine calc_wlbd(bl_flag, bl_calc, de_flag, de_calc, wl_flag, wl_calc, wl_num)
  !***************************************************************************************
  ! calc_wlbd -- Calculate water level from bottom level and water depth
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: bl_flag(:), de_flag(:)
    real(SP), intent(in) :: bl_calc(:), de_calc(:)
    integer(I4), intent(out) :: wl_flag(:)
    real(SP), intent(out) :: wl_calc(:)
    integer(I4), intent(out) :: wl_num
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    wl_num = 0
    !$omp parallel
    !$omp workshare
    wl_flag(:) = 0 ; wl_calc(:) = SZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      if (bl_flag(i) == 1 .and. de_flag(i) == 1) then
        wl_calc(i) = bl_calc(i) + de_calc(i)
        wl_flag(i) = 1
      end if
    end do
    !$omp end do

    !$omp workshare
    wl_num = sum(wl_flag)
    !$omp end workshare
    !$omp end parallel

  end subroutine calc_wlbd

  subroutine calc_blld(wl_flag, wl_calc, de_flag, de_calc, bl_flag, bl_calc, wb_num)
  !***************************************************************************************
  ! calc_blld -- Calculate bottom level from water level and water depth
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wl_flag(:), de_flag(:)
    real(SP), intent(in) :: wl_calc(:), de_calc(:)
    integer(I4), intent(out) :: bl_flag(:)
    real(SP), intent(out) :: bl_calc(:)
    integer(I4), intent(out) :: wb_num
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    wb_num = 0
    !$omp parallel
    !$omp workshare
    bl_flag(:) = 0 ; bl_calc(:) = SZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      if (wl_flag(i) == 1 .and. de_flag(i) == 1) then
        bl_calc(i) = wl_calc(i) - de_calc(i)
        bl_flag(i) = 1
      end if
    end do
    !$omp end do

    !$omp workshare
    wb_num = sum(bl_flag)
    !$omp end workshare
    !$omp end parallel

  end subroutine calc_blld

  subroutine calc_wlsl(wd_flag, wd_calc, wl_flag, wl_calc)
  !***************************************************************************************
  ! calc_wlsl -- Calculate water level from surface level
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wd_flag(:)
    real(SP), intent(in) :: wd_calc(:)
    integer(I4), intent(out) :: wl_flag(:)
    real(SP), intent(out) :: wl_calc(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel
    !$omp workshare
    wl_flag(:) = 0 ; wl_calc(:) = SZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      if (wd_flag(i) == 1) then
        wl_calc(i) = real(surf_elev(i) + wd_calc(i), kind=SP)
        wl_flag(i) = 1
      end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine calc_wlsl

  subroutine calc_blsl(de_flag, de_calc, bl_flag, bl_calc, bl_num)
  !***************************************************************************************
  ! calc_blsl -- Calculate bottom level from surface level
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: de_flag(:)
    real(SP), intent(in) :: de_calc(:)
    integer(I4), intent(out) :: bl_flag(:)
    real(SP), intent(out) :: bl_calc(:)
    integer(I4), intent(out) :: bl_num
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    bl_num = 0
    !$omp parallel
    !$omp workshare
    bl_flag(:) = 0 ; bl_calc(:) = SZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      if (de_flag(i) == 1) then
        bl_calc(i) = real(surf_elev(i) - de_calc(i), kind=SP)
        bl_flag(i) = 1
      end if
    end do
    !$omp end do

    !$omp workshare
    bl_num = sum(bl_flag)
    !$omp end workshare
    !$omp end parallel

  end subroutine calc_blsl

  subroutine calc_lsurf(in_flag, out_flag, surf_out, out_num)
  !***************************************************************************************
  ! calc_lsurf -- Calculate level from surface
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: in_flag(:)
    integer(I4), intent(out) :: out_flag(:)
    real(SP), intent(out) :: surf_out(:)
    integer(I4), intent(out) :: out_num
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    out_num = 0
    !$omp parallel
    !$omp workshare
    out_flag(:) = 0 ; surf_out(:) = SZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      if (in_flag(i) == 1) then
        surf_out(i) = real(surf_elev(i), kind=SP)
        out_flag(i) = 1
      end if
    end do
    !$omp end do

    !$omp workshare
    out_num = sum(out_flag)
    !$omp end workshare
    !$omp end parallel

  end subroutine calc_lsurf

  subroutine calc_rivea(wi_flag, le_flag, riv_wi, riv_le, ar_flag, riv_ar, riar_num)
  !***************************************************************************************
  ! calc_rivea -- Calculate river area
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wi_flag(:), le_flag(:)
    real(SP), intent(in) :: riv_wi(:), riv_le(:)
    integer(I4), intent(out) :: ar_flag(:)
    real(SP), intent(out) :: riv_ar(:)
    integer(I4), intent(out) :: riar_num
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    riar_num = 0
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      if (wi_flag(i) == 1 .and. le_flag(i) == 1) then
        riv_ar(i) = riv_wi(i)*riv_le(i)
        ar_flag(i) = 1
      end if
    end do
    !$omp end do

    !$omp workshare
    riar_num = sum(ar_flag)
    !$omp end workshare
    !$omp end parallel

  end subroutine calc_rivea

  subroutine count_rivecalc(wl_flag, bl_flag, ar_flag, riv_wi, riv_bl, riv_ar, riv_cnum)
  !***************************************************************************************
  ! count_rivecalc -- Count river calculation
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wl_flag(:), bl_flag(:), ar_flag(:)
    real(SP), intent(in) :: riv_wi(:), riv_bl(:), riv_ar(:)
    integer(I4), intent(inout) :: riv_cnum
    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: rive_cflag(:)
    !-------------------------------------------------------------------------------------
    riv_cnum = 0
    allocate(rive_cflag(ncals))
    !$omp parallel
    !$omp workshare
    rive_cflag(:) = 0
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      if (wl_flag(i) == 1 .and. bl_flag(i) == 1 .and. ar_flag(i) == 1) then
        rive_cflag(i) = 1
      end if
    end do
    !$omp end do

    !$omp workshare
    riv_cnum = sum(rive_cflag)
    !$omp end workshare
    !$omp end parallel

    allocate(rive2cals(riv_cnum))
    !$omp parallel workshare
    rive2cals(:) = 0
    !$omp end parallel workshare

    if (riv_cnum > 0) then
      allocate(rive_head(riv_cnum), rive_bott(riv_cnum), rive_area(riv_cnum))
      !$omp parallel workshare
      rive_head(:) = DZERO ; rive_bott(:) = DZERO ; rive_area(:) = DZERO
      !$omp end parallel workshare
      call set_bound2calc(ncals, rive_cflag, riv_wi, rive2cals, rive_head)
      call set_bound2calc(ncals, rive_cflag, riv_bl, rive2cals, rive_bott)
      call set_bound2calc(ncals, rive_cflag, riv_ar, rive2cals, rive_area)
    end if

    deallocate(rive_cflag)

  end subroutine count_rivecalc

  subroutine count_lakecalc(wl_flag, bl_flag, ar_flag, lak_wi, lak_bl, lak_ar, lak_cnum)
  !***************************************************************************************
  ! count_lakecalc -- Count lake calculation cell
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wl_flag(:), bl_flag(:), ar_flag(:)
    real(SP), intent(in) :: lak_wi(:), lak_bl(:), lak_ar(:)
    integer(I4), intent(inout) :: lak_cnum
    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: lake_cflag(:)
    !-------------------------------------------------------------------------------------
    lak_cnum = 0
    allocate(lake_cflag(ncals))
    !$omp parallel
    !$omp workshare
    lake_cflag(:) = 0
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      if (wl_flag(i) == 1 .and. bl_flag(i) == 1 .and. ar_flag(i) == 1) then
        lake_cflag(i) = 1
      end if
    end do
    !$omp end do

    !$omp workshare
    lak_cnum = sum(lake_cflag)
    !$omp end workshare
    !$omp end parallel

    allocate(lake2cals(lak_cnum))
    !$omp parallel workshare
    lake2cals(:) = 0
    !$omp end parallel workshare

    if (lak_cnum > 0) then
      allocate(lake_head(lak_cnum), lake_bott(lak_cnum), lake_area(lak_cnum))
      !$omp parallel workshare
      lake_head(:) = DZERO ; lake_bott(:) = DZERO ; lake_area(:) = DZERO
      !$omp end parallel workshare
      call set_bound2calc(ncals, lake_cflag, lak_wi, lake2cals, lake_head)
      call set_bound2calc(ncals, lake_cflag, lak_bl, lake2cals, lake_bott)
      call set_bound2calc(ncals, lake_cflag, lak_ar, lake2cals, lake_area)
    end if

    deallocate(lake_cflag)

  end subroutine count_lakecalc

end module calc_boundary
