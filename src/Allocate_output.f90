module allocate_output
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO
  use set_cell, only: ncalc, ncals
  use assign_calc, only: msout_tnum

  implicit none
  private
  public :: allocate_outvar
  public :: allocate_wtab, allocate_mass, allocate_velc, allocate_rivr, allocate_lakr
  public :: allocate_sufr, allocate_dunr, allocate_sear, allocate_recr, allocate_welr

  type, public :: st_msout
    real(DP), allocatable :: sto(:), con(:), sea(:), wel(:), rec(:), sur(:), riv(:), lak(:), tot(:)
  end type st_msout
  type(st_msout), public :: st_msloc, st_msglo

  real(DP), allocatable, public :: wtable(:)
  real(DP), allocatable, public :: pointv(:,:), facev(:,:)
  real(DP), allocatable, public :: roff_rive(:), roff_lake(:), roff_surf(:), roff_dunn(:)
  real(DP), allocatable, public :: res_seal(:), res_rech(:), res_well(:)
  integer(I4), allocatable, public :: res_snum(:), res_rnum(:), res_wnum(:)
  real(DP), public :: rive_sumtime, lake_sumtime, surf_sumtime, dunn_sumtime
  character(:), allocatable, public :: ms_head

  contains

  subroutine allocate_outvar()
  !***************************************************************************************
  ! allocate_outvar -- Allocate output variable
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_in_type, st_out_type, in_type, out_type
    use assign_calc, only: massout_name
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    ! water table file
    if (st_out_type%wtab == out_type(2)) then
      ! -- Allocate water table (wtab)
        call allocate_wtab()
    end if

    ! massbalance file
    if (st_out_type%mass == out_type(1)) then
      if (st_in_type%mass /= in_type(7) .and. st_in_type%mass /= in_type(1)) then
        allocate(massout_name(msout_tnum))
        !$omp parallel
        !$omp workshare
        massout_name(:) = ""
        !$omp end workshare

        !$omp do private(i)
        do i = 1, msout_tnum
          write(massout_name(i),*) i
        end do
        !$omp end do
        !$omp end parallel
      else if (st_in_type%mass == in_type(7)) then
        allocate(massout_name(1))
        massout_name(1) = "whole area"
      end if

      do i = 1, msout_tnum
        if (i == 1) then
          ms_head = trim(adjustl(massout_name(i)))
        else
          ms_head = trim(adjustl(ms_head))//",,,,,,,,,"//trim(adjustl(massout_name(i)))
        end if
      end do
      ! -- Allocate massbalance (mass)
        call allocate_mass()
    end if

    ! velocity file
    if (st_out_type%velc == out_type(3)) then
      ! -- Allocate velocity (velc)
        call allocate_velc()
    end if

    ! river runoff file
    if (st_out_type%rivr == out_type(2)) then
      ! -- Allocate river runoff (rivr)
        call allocate_rivr()
    end if

    ! lake runoff file
    if (st_out_type%lakr == out_type(2)) then
      ! -- Allocate lake runoff (lakr)
        call allocate_lakr()
    end if

    ! surface runoff file
    if (st_out_type%sufr == out_type(2)) then
      ! -- Allocate surface runoff (sufr)
        call allocate_sufr()
    end if

    ! dunne runoff file
    if (st_out_type%dunr == out_type(2)) then
      ! -- Allocate dunne runoff (dunr)
        call allocate_dunr()
    end if

    ! sea results file
    if (st_out_type%seal == out_type(3)) then
      ! -- Allocate seal results (sear)
        call allocate_sear()
    end if

    ! rech results file
    if (st_out_type%rech == out_type(2)) then
      ! -- Allocate recharge results (recr)
        call allocate_recr()
    end if

    ! well pumping results file
    if (st_out_type%well == out_type(3)) then
      ! -- Allocate well results (welr)
        call allocate_welr()
    end if

  end subroutine allocate_outvar

  subroutine allocate_wtab()
  !***************************************************************************************
  ! allocate_wtab -- Allocate water table
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(wtable(ncals))
    !$omp parallel workshare
    wtable(:) = DZERO
    !$omp end parallel workshare

  end subroutine allocate_wtab

  subroutine allocate_mass()
  !***************************************************************************************
  ! allocate_mass -- Allocate massbalance
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(st_msloc%sto(ncalc), st_msloc%con(ncalc), st_msloc%sea(ncalc))
    allocate(st_msloc%wel(ncalc), st_msloc%rec(ncals), st_msloc%sur(ncals))
    allocate(st_msloc%riv(ncals), st_msloc%lak(ncals))
    allocate(st_msglo%sto(msout_tnum), st_msglo%con(msout_tnum), st_msglo%sea(msout_tnum))
    allocate(st_msglo%wel(msout_tnum), st_msglo%rec(msout_tnum), st_msglo%sur(msout_tnum))
    allocate(st_msglo%riv(msout_tnum), st_msglo%lak(msout_tnum), st_msglo%tot(msout_tnum))
    !$omp parallel workshare
    st_msloc%sto(:) = DZERO ; st_msloc%con(:) = DZERO ; st_msloc%sea(:) = DZERO
    st_msloc%wel(:) = DZERO ; st_msloc%rec(:) = DZERO ; st_msloc%sur(:) = DZERO
    st_msloc%riv(:) = DZERO ; st_msloc%lak(:) = DZERO
    st_msglo%sto(:) = DZERO ; st_msglo%con(:) = DZERO ; st_msglo%sea(:) = DZERO
    st_msglo%wel(:) = DZERO ; st_msglo%rec(:) = DZERO ; st_msglo%sur(:) = DZERO
    st_msglo%riv(:) = DZERO ; st_msglo%lak(:) = DZERO ; st_msglo%tot(:) = DZERO
    !$omp end parallel workshare

  end subroutine allocate_mass

  subroutine allocate_velc()
  !***************************************************************************************
  ! allocate_velc -- Allocate velocity
  !***************************************************************************************
    ! -- modules
    use constval_module, only: FACE
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(pointv(ncalc,3), facev(ncalc,FACE))
    !$omp parallel workshare
    pointv(:,:) = DZERO ; facev(:,:) = DZERO
    !$omp end parallel workshare

  end subroutine allocate_velc

  subroutine allocate_rivr()
  !***************************************************************************************
  ! allocate_rivr -- Allocate river runoff
  !***************************************************************************************
    ! -- modules
    use set_boundary, only: rive_num
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(roff_rive(rive_num))
    !$omp parallel workshare
    roff_rive(:) = DZERO
    !$omp end parallel workshare

    rive_sumtime = DZERO

  end subroutine allocate_rivr

  subroutine allocate_lakr()
  !***************************************************************************************
  ! allocate_lakr -- Allocate lake runoff
  !***************************************************************************************
    ! -- modules
    use set_boundary, only: lake_num
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(roff_lake(lake_num))
    !$omp parallel workshare
    roff_lake(:) = DZERO
    !$omp end parallel workshare

    lake_sumtime = DZERO

  end subroutine allocate_lakr

  subroutine allocate_sufr()
  !***************************************************************************************
  ! allocate_sufr -- Allocate surface runoff
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(roff_surf(ncals))
    !$omp parallel workshare
    roff_surf(:) = DZERO
    !$omp end parallel workshare

    surf_sumtime = DZERO

  end subroutine allocate_sufr

  subroutine allocate_dunr()
  !***************************************************************************************
  ! allocate_dunr -- Allocate dunne runoff
  !***************************************************************************************
    ! -- modules
    use set_boundary, only: rech_num
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(roff_dunn(rech_num))
    !$omp parallel workshare
    roff_dunn(:) = DZERO
    !$omp end parallel workshare

    dunn_sumtime = DZERO

  end subroutine allocate_dunr

  subroutine allocate_sear()
  !***************************************************************************************
  ! allocate_sear -- Allocate sea results
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(res_seal(ncalc), res_snum(ncalc))
    !$omp parallel workshare
    res_seal(:) = DZERO ; res_snum(:) = 0
    !$omp end parallel workshare

  end subroutine allocate_sear

  subroutine allocate_recr()
  !***************************************************************************************
  ! allocate_recr -- Allocate recharge results
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(res_rech(ncals), res_rnum(ncals))
    !$omp parallel workshare
    res_rech(:) = DZERO ; res_rnum(:) = 0
    !$omp end parallel workshare

  end subroutine allocate_recr

  subroutine allocate_welr()
  !***************************************************************************************
  ! allocate_welr -- Allocate well results
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(res_well(ncalc), res_wnum(ncalc))
    !$omp parallel workshare
    res_well(:) = DZERO ; res_wnum(:) = 0
    !$omp end parallel workshare

  end subroutine allocate_welr

end module allocate_output
