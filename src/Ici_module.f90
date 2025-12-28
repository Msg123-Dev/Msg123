module ici_module
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: DZERO
  use initial_module, only: my_rank, st_grid, st_out_type, out_type
  use read_input, only: len_scal
  use set_cell, only: ncals, ncalc, get_calc_grid
  use set_condition, only: rech_area
  use prep_calculation, only: now_date
  use allocate_output, only: wtable
  use calc_output, only: calc_wtable
  use palmtime, only: palm_TimeStart, palm_TimeEnd
  use ici_api, only: ici_put_data
  use mpi

  implicit none
  private
  public :: init_ici, set_mapt, put_initv, get_var, alloc_outvar, put_var, fin_ici

  ! -- local
  integer(I4) :: int_mat, int_cama, sum_cals
  integer(I4), allocatable :: xy2cals(:)
  real(DP) :: mat_noval, cama_noval
  real(DP), allocatable, save :: water_in(:), roff_a(:), roff_b(:), evap_s(:), evap_v(:)
  character(:), allocatable :: my_comp, my_grid, ici_file
  logical :: coupled_matsiro = .false., coupled_cama = .false.
  logical :: mat_get = .false., mat_put = .false., cama_get = .false., cama_put = .false.

  contains

  subroutine init_ici()
  !***************************************************************************************
  ! init_ici -- Initialize ici
  !***************************************************************************************
    ! -- module
    use mpi_initfin, only: my_comm
    use initial_module, only: pro_totn
    use mpi_set, only: bcast_ici_set
    use ici_api, only: ici_init, ici_get_comm_local, ici_get_numpe_local,&
                       ici_get_irank_local, ici_is_coupled
    use palmtime, only: palm_TimeInit
    ! -- inout

    ! -- local
    integer(I4) :: ierr, errcode
    logical :: mpi_init_check
    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      ! -- Read ici main file (ici_main)
        call read_ici_main()
    end if

    ! -- Bcast ici set (ici_set)
      call bcast_ici_set(my_comp, my_grid, ici_file, mat_get, mat_put, cama_get, cama_put)

    call ici_init(my_comp, ici_file)
    call palm_TimeInit("MSG", comm=ici_get_comm_local())

    call palm_TimeStart('Initialize')

    pro_totn = ici_get_numpe_local()
    my_rank = ici_get_irank_local()
    my_comm = ici_get_comm_local()

    coupled_matsiro  = ici_is_coupled('MATSIRO')
    coupled_cama  = ici_is_coupled('CAMA')

    call palm_TimeEnd('Initialize')

  end subroutine init_ici

  subroutine set_mapt()
  !***************************************************************************************
  ! set_mapt -- Set mapping table
  !***************************************************************************************
    ! -- module
    use ici_api, only: ici_init_time
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    call palm_TimeStart('Setup')

    ! -- Define grid (grid)
      call define_grid()

    call ici_init_time(now_date)

    ! -- Get exchange interval (exint)
      call get_exint()

    call palm_TimeEnd('Setup')

  end subroutine set_mapt

  subroutine get_var()
  !***************************************************************************************
  ! get_var -- Get variables
  !***************************************************************************************
    ! -- module
    use constval_module, only: SZERO
    use initial_module, only: st_rech, st_riwd, st_step_flag
    use prep_calculation, only: current_t, inter_time
    use assign_boundary, only: read_rech, rech_cflag
    use set_boundary, only: cflag_riv, criv, rech_num, rivnum
    use ici_api, only: ici_set_time, ici_get_data, ici_get_get_fill_value
    ! -- inout

    ! -- local
    real(DP), allocatable :: rive_wd(:)
    logical :: wain_get, runa_get, runb_get, evas_get, evav_get, rivd_get
    !-------------------------------------------------------------------------------------
    if (current_t == DZERO) then
      call ici_set_time(now_date, min(int_mat, int_cama))
    else
      call ici_set_time(now_date, nint(inter_time))
    end if

    if (coupled_matsiro .and. mat_get) then
      if (.not. allocated(read_rech)) then
        rech_num = ncals
        allocate(read_rech(ncals))
        !$omp parallel workshare
        read_rech(:) = SZERO
        !$omp end parallel workshare
      end if
      if (.not. allocated(rech_cflag)) then
        allocate(rech_cflag(ncals))
        !$omp parallel workshare
        rech_cflag(:) = 1
        !$omp end parallel workshare
      end if
      call ici_get_data("water_input", water_in, IS_GET_OK=wain_get)
      call ici_get_data("runoff_all", roff_a, IS_GET_OK=runa_get)
      call ici_get_data("runoff_base", roff_b, IS_GET_OK=runb_get)
      call ici_get_data("evap_soil", evap_s, IS_GET_OK=evas_get)
      call ici_get_data("evap_vegt", evap_v, IS_GET_OK=evav_get)
      if (current_t == DZERO) then
        mat_noval = ici_get_get_fill_value('water_input')
        if (wain_get .and. runa_get .and. runb_get .and. evas_get .and. evav_get) then
          call calc_infil(water_in, roff_a, roff_b, evap_s, evap_v, read_rech)
          st_rech%etime = real(int_mat, kind=SP)
        else if (my_rank == 0) then
          call write_err_stop("Not get initial variables from MATSIRO.")
        end if
      else if (st_step_flag%rech == 1) then
        if (wain_get .or. runa_get .or. runb_get .or. evas_get .or. evav_get) then
          call calc_infil(water_in, roff_a, roff_b, evap_s, evap_v, read_rech)
          st_rech%etime = st_rech%etime + real(int_mat, kind=SP)
        end if
      end if
    end if

    if (coupled_cama .and. cama_get) then
      if (.not. allocated(criv%wd)) then
        allocate(criv%wd(ncals))
        !$omp parallel workshare
        criv%wd(:) = SZERO
        !$omp end parallel workshare
      end if
      if (.not. allocated(cflag_riv%wd)) then
        allocate(cflag_riv%wd(ncals))
        !$omp parallel workshare
        cflag_riv%wd(:) = 1
        !$omp end parallel workshare
      end if
      allocate(rive_wd(ncals))
      !$omp parallel workshare
      rive_wd(:) = DZERO
      !$omp end parallel workshare
      call ici_get_data("rive_wdep", rive_wd, IS_GET_OK=rivd_get)
      if (current_t == DZERO) then
        cama_noval = ici_get_get_fill_value('rive_wdep')
        if (rivd_get) then
          call conv_rive(rive_wd, cflag_riv%wd, criv%wd, rivnum%wd)
          st_riwd%etime = real(int_cama, kind=SP)
        else if (my_rank == 0) then
          call write_err_stop("Not get initial variables from CaMa-Flood.")
        end if
      else if (st_step_flag%riwd == 1) then
        if (rivd_get) then
          call conv_rive(rive_wd, cflag_riv%wd, criv%wd, rivnum%wd)
          st_riwd%etime = st_riwd%etime + real(int_cama, kind=SP)
        end if
      end if
      deallocate(rive_wd)
    end if

  end subroutine get_var

  subroutine alloc_outvar()
  !***************************************************************************************
  ! alloc_outvar -- Allocate output variables
  !***************************************************************************************
    ! -- module
    use allocate_output, only: allocate_wtab, allocate_rivr, allocate_lakr, allocate_sufr
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    if (coupled_matsiro .and. mat_put) then
      if (st_out_type%wtab /= out_type(2)) then
        ! -- Allocate water table (wtab)
          call allocate_wtab()
      end if
    end if

    if (coupled_cama .and. cama_put) then
      if (st_out_type%rivr /= out_type(2)) then
        ! -- Allocate river runoff (rivr)
          call allocate_rivr()
      end if
      if (st_out_type%lakr /= out_type(2)) then
        ! -- Allocate lake runoff (lakr)
          call allocate_lakr()
      end if
      if (st_out_type%sufr /= out_type(2)) then
        ! -- Allocate surface runoff (sufr)
          call allocate_sufr()
      end if
    end if

  end subroutine alloc_outvar

  subroutine put_var()
  !***************************************************************************************
  ! put_var -- Put variables
  !***************************************************************************************
    ! -- module
    use initial_module, only: st_in_type
    use allocate_solution, only: head_new, srat_new
    use allocate_output, only: roff_rive, roff_lake, roff_surf
    use calc_output, only: calc_rivr_off, calc_lakr_off, calc_sufr_off
    ! -- inout

    ! -- local
    integer(I4) :: i, s, nz
    real(DP), allocatable :: head_ici(:), wtab_dep(:)
    real(DP), allocatable :: rive_flux(:), lake_flux(:), surf_flux(:), tot_flux(:)
    !-------------------------------------------------------------------------------------
    allocate(head_ici(ncalc))
    !$omp parallel workshare
    head_ici(:) = DZERO
    !$omp end parallel workshare

    call change_ici_put(head_new, head_ici)
    call ici_put_data("hyd_head", reshape([head_ici(:)*len_scal], [ncals,st_grid%nz]))

    if (coupled_matsiro .and. mat_put) then
      if (st_out_type%wtab /= out_type(2)) then
        ! -- Calculate water table (wtable)
          call calc_wtable(head_new, srat_new)
      end if
      allocate(wtab_dep(ncals))
      ! -- Calculate water table depth (wtab_depth)
        call calc_wtab_depth(wtable, wtab_dep)

      call ici_put_data("wtab_depth", wtab_dep(:)*len_scal)

      deallocate(wtab_dep)
    end if

    if (coupled_cama .and. cama_put) then
      allocate(rive_flux(ncals), lake_flux(ncals), surf_flux(ncals), tot_flux(ncals))
      !$omp parallel workshare
      rive_flux(:) = DZERO ; lake_flux(:) = DZERO ; surf_flux(:) = DZERO
      tot_flux(:) = DZERO
      !$omp end parallel workshare
      if (st_out_type%rivr /= out_type(2) .and. st_in_type%rive >= 0) then
        ! -- Calculate river runoff (rive_roff)
          call calc_rivr_off()
        ! -- Calculate river flux (rive_flux)
          call calc_rive_flux(roff_rive, rive_flux)
      end if
      if (st_out_type%lakr /= out_type(2) .and. st_in_type%lake >= 0) then
        ! -- Calculate lake runoff (lake_roff)
          call calc_lakr_off()
        ! -- Calculate lake flux (lake_flux)
          call calc_lake_flux(roff_lake, lake_flux)
      end if
      if (st_out_type%sufr /= out_type(2)) then
        ! -- Calculate surface runoff (surf_roff)
          call calc_sufr_off()
        ! -- Calculate surface flux (surf_flux)
          call calc_surf_flux(roff_surf, surf_flux)
      end if
      !$omp parallel workshare
      tot_flux(:) = rive_flux(:) + lake_flux(:) + surf_flux(:)
      !$omp end parallel workshare

      call ici_put_data("runoff_total", tot_flux(:)*len_scal)

      !$omp parallel workshare
      roff_rive(:) = DZERO ; roff_lake(:) = DZERO ; roff_surf(:) = DZERO
      !$omp end parallel workshare

      deallocate(tot_flux, rive_flux, lake_flux, surf_flux)
    end if

  end subroutine put_var

  subroutine fin_ici()
  !***************************************************************************************
  ! fin_ici -- Finalize ici
  !***************************************************************************************
    ! -- module
    use ici_api, only: ici_finalize
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    call ici_finalize(.true., .false.)

  end subroutine fin_ici

  subroutine read_ici_main()
  !***************************************************************************************
  ! read_ici_main -- Set ici mail file
  !***************************************************************************************
    ! -- module
    use utility_module, only: open_new_rtxt
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: main_fnum
    character(:), allocatable :: main_file
    namelist/set_ici/my_comp, my_grid, ici_file, mat_get, mat_put, cama_get, cama_put
    !-------------------------------------------------------------------------------------
    main_file = 'msg.main'
    ! -- Open new read text file (new_rtxt)
      call open_new_rtxt(1, 1, main_file, "msg123 main", main_fnum)

    ierr = 0
    read(unit=main_fnum,nml=set_ici,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading Integrated Land Simulator(ILS) section in main file.")
    end if

    call close_file(main_fnum)

  end subroutine read_ici_main

  subroutine define_grid()
  !***************************************************************************************
  ! define_grid -- Define grid
  !***************************************************************************************
    ! -- module
    use ici_api, only: ici_def_grid, ici_end_grid_def
    use set_cell, only: loc2unk_ij
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: msg_idex(:)
    !-------------------------------------------------------------------------------------
    allocate(msg_idex(ncals))
    !$omp parallel
    !$omp workshare
    msg_idex(:) = 0
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      msg_idex(i) = loc2unk_ij(i)
    end do
    !$omp end do
    !$omp end parallel

    call ici_def_grid(my_grid, ncals, 1, 1, msg_idex)
    call ici_end_grid_def()

    deallocate(msg_idex, loc2unk_ij)

  end subroutine define_grid

  subroutine put_initv()
  !***************************************************************************************
  ! put_initv -- Put initial variables
  !***************************************************************************************
    ! -- module
    use assign_calc, only: read_init
    ! -- inout

    ! -- local
    integer(I4) :: i, xn, yn, zn, xy
    integer(I4) :: coun_cals
    real(DP), allocatable :: init_val(:), dummy_var(:)
    !-------------------------------------------------------------------------------------
    allocate(init_val(ncalc))
    !$omp parallel workshare
    init_val(:) = DZERO
    !$omp end parallel workshare

    call make_ici_put()
    call change_ici_put(read_init, init_val)
    call ici_put_data("hyd_head", reshape([init_val(:)*len_scal], [ncals,st_grid%nz]))

    deallocate(init_val)

    if (coupled_matsiro .and. mat_put) then
      allocate(dummy_var(ncals))
      !$omp parallel workshare
      dummy_var(:) = DZERO
      !$omp end parallel workshare
      call ici_put_data("wtab_depth", dummy_var)

      deallocate(dummy_var)
    end if

    if (coupled_cama .and. cama_put) then
      allocate(dummy_var(ncals))
      !$omp parallel workshare
      dummy_var(:) = DZERO
      !$omp end parallel workshare
      call ici_put_data("runoff_total", dummy_var)

      deallocate(dummy_var)
    end if

  end subroutine put_initv

  subroutine get_exint()
  !***************************************************************************************
  ! get_exint -- Get exchange interval
  !***************************************************************************************
    ! -- module
    use constval_module, only: VARLEN
    use ici_api, only: ici_get_exchange_interval_get, ici_get_num_of_configuration,&
                       ici_get_num_of_exchange, ici_get_get_data_name
    ! -- inout

    ! -- local
    integer(I4) :: i, j, count
    integer(I4) :: int_winp, int_runa, int_runb, int_evas, int_evav, int_riwd
    character(:), allocatable :: data_name
    character(9), allocatable :: exch_name(:)
    character(VARLEN), allocatable :: exdata_name(:)
    !-------------------------------------------------------------------------------------
    allocate(exch_name(6), exdata_name(6))
    exch_name = ["water_inp", "runoff_al", "runoff_ba", "evap_soil", "evap_vegt", "rive_wdep"]
    count = 0
    do i = 1, ici_get_num_of_configuration()
      do j = 1, ici_get_num_of_exchange(i)
        data_name = ici_get_get_data_name(i,j)
        if (any(data_name(1:9) == exch_name)) then
          count = count + 1
          exdata_name(count) = trim(data_name)
        end if
      end do
    end do

    if (coupled_matsiro .and. mat_get) then
      call ici_get_exchange_interval_get("MSG", trim(exdata_name(1)), int_winp)
      call ici_get_exchange_interval_get("MSG", trim(exdata_name(2)), int_runa)
      call ici_get_exchange_interval_get("MSG", trim(exdata_name(3)), int_runb)
      call ici_get_exchange_interval_get("MSG", trim(exdata_name(4)), int_evas)
      call ici_get_exchange_interval_get("MSG", trim(exdata_name(5)), int_evav)
      int_mat = min(int_winp, int_runa, int_runb, int_evas, int_evav)

      allocate(water_in(ncals), roff_a(ncals), roff_b(ncals))
      allocate(evap_s(ncals), evap_v(ncals))
      !$omp parallel workshare
      water_in(:) = DZERO ; roff_a(:) = DZERO ; roff_b(:) = DZERO
      evap_s(:) = DZERO ; evap_v(:) = DZERO
      !$omp end parallel workshare
    end if

    if (coupled_cama .and. cama_get) then
      call ici_get_exchange_interval_get("MSG", trim(exdata_name(6)), int_riwd)
      int_cama = int_riwd
    end if

    deallocate(exch_name, exdata_name)

  end subroutine get_exint

  subroutine make_ici_put()
  !***************************************************************************************
  ! make_ici_put -- Make for ici put
  !***************************************************************************************
    ! -- module

    ! -- inout

    ! -- local
    integer(I4) :: i, xn, yn, zn, xy
    integer(I4) :: coun_cals
    !-------------------------------------------------------------------------------------
    allocate(xy2cals(st_grid%nx*st_grid%ny))
    !$omp parallel workshare
    xy2cals(:) = 0
    !$omp end parallel workshare

    coun_cals = 0
    do i = 1, ncalc
      call get_calc_grid(i, xn, yn, zn)
      if (zn == 1) then
        coun_cals = coun_cals + 1 ; xy = st_grid%nx*(yn-1) + xn
        xy2cals(xy) = coun_cals
      end if
    end do

    sum_cals = coun_cals

  end subroutine make_ici_put

  subroutine change_ici_put(in_head, out_head)
  !***************************************************************************************
  ! change_ici_put -- Change for ici put
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: in_head(:)
    real(DP), intent(out) :: out_head(:)
    ! -- local
    integer(I4) :: i, xn, yn, zn, xy
    integer(I4) :: coun_calc
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i,  xn, yn, zn, xy, coun_calc)
    do i = 1, ncalc
      call get_calc_grid(i, xn, yn, zn)
      xy = st_grid%nx*(yn-1) + xn ; coun_calc = xy2cals(xy) + sum_cals*(zn-1)
      out_head(coun_calc) = in_head(i)
    end do
    !$omp end parallel do

  end subroutine change_ici_put

  subroutine calc_infil(winpt, runofa, runofb, evapso, evapve, infilt)
  !***************************************************************************************
  ! calc_infil -- Calculation infiltration
  !***************************************************************************************
    ! -- module
    use read_input, only: len_scal_inv
    ! -- inout
    real(DP), intent(in) :: winpt(:), runofa(:), runofb(:), evapso(:), evapve(:)
    real(SP), intent(out) :: infilt(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, ncals
      if (winpt(i)/=mat_noval .and. runofa(i)/=mat_noval .and. runofb(i)/=mat_noval .and.&
          evapso(i)/=mat_noval .and. evapve(i)/=mat_noval) then
        infilt(i) = real((winpt(i)+runofa(i)-runofb(i)-evapso(i)-evapve(i))*len_scal_inv, kind=SP)
      end if
    end do
    !$omp end parallel do

  end subroutine calc_infil

  subroutine conv_rive(inwd, wdflag, outwd, wdnum)
  !***************************************************************************************
  ! conv_rive -- Convert river
  !***************************************************************************************
    ! -- module
    use read_input, only: len_scal_inv
    ! -- inout
    real(DP), intent(in) :: inwd(:)
    integer(I4), intent(out) :: wdflag(:)
    real(SP), intent(out) :: outwd(:)
    integer(I4), intent(out) :: wdnum
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, ncals
      if (inwd(i) /= cama_noval) then
        outwd(i) = real(inwd(i)*len_scal_inv, kind=SP)
        wdflag(i) = 1
      end if
    end do
    !$omp end parallel do
    wdnum = sum(wdflag)

  end subroutine conv_rive

  subroutine calc_wtab_depth(wtab, wtab_depth)
  !***************************************************************************************
  ! calc_wtab_depth -- Calculate water table depth
  !***************************************************************************************
    ! -- module
    use make_cell, only: surf_elev
    ! -- inout
    real(DP), intent(in) :: wtab(:)
    real(DP), intent(out) :: wtab_depth(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel
    !$omp workshare
    wtab_depth(:) = DZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      wtab_depth(i) = surf_elev(i) - wtab(i)
      if (wtab_depth(i) < DZERO) then
        wtab_depth(i) = DZERO
      end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine calc_wtab_depth

  subroutine calc_rive_flux(rivoff, rivflux)
  !***************************************************************************************
  ! calc_rive_flux -- Calculate river flux
  !***************************************************************************************
    ! -- module
    use calc_boundary, only: rive2cals
    use set_boundary, only: rive_num
    use allocate_output, only: rive_sumtime
    ! -- inout
    real(DP), intent(in) :: rivoff(:)
    real(DP), intent(out) :: rivflux(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, s)
    do i = 1, rive_num
      s = rive2cals(i)
      rivflux(i) = rivoff(i)/rech_area(s)/rive_sumtime
    end do
    !$omp end parallel do

    rive_sumtime = DZERO

  end subroutine calc_rive_flux

  subroutine calc_lake_flux(lakoff, lakflux)
  !***************************************************************************************
  ! calc_lake_flux -- Calculate lake flux
  !***************************************************************************************
    ! -- module
    use calc_boundary, only: lake2cals
    use set_boundary, only: lake_num
    use allocate_output, only: lake_sumtime
    ! -- inout
    real(DP), intent(in) :: lakoff(:)
    real(DP), intent(out) :: lakflux(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, s)
    do i = 1, lake_num
      s = lake2cals(i)
      lakflux(i) = lakoff(i)/rech_area(s)/lake_sumtime
    end do
    !$omp end parallel do

    lake_sumtime = DZERO

  end subroutine calc_lake_flux

  subroutine calc_surf_flux(sufoff, sufflux)
  !***************************************************************************************
  ! calc_surf_flux -- Calculate surface flux
  !***************************************************************************************
    ! -- module
    use allocate_output, only: surf_sumtime
    ! -- inout
    real(DP), intent(in) :: sufoff(:)
    real(DP), intent(out) :: sufflux(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, ncals
      sufflux(i) = sufoff(i)/rech_area(i)/surf_sumtime
    end do
    !$omp end parallel do

    surf_sumtime = DZERO

  end subroutine calc_surf_flux

end module ici_module
