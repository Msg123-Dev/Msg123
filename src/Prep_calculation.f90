module prep_calculation
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: DZERO
  use initial_module, only: st_in_type, st_in_path, in_type
  use set_cell, only: ncalc
#ifdef MPI_MSG
  use initial_module, only: pro_totn
  use mpi_set, only: cals_r4view, calc_r4view
#endif

  implicit none
  private
  public :: prepare_calc
  integer(I4), public :: now_date(6), out_iter, conv_flag, form_switch
  real(SP), public :: current_t, inter_time
  real(DP), public :: delt, delt_inv
  real(DP), allocatable, public :: surf_top(:)

  ! -- local

  contains

  subroutine prepare_calc()
  !***************************************************************************************
  ! prepare_calc -- Prepare calculation
  !***************************************************************************************
    ! -- modules
    use constval_module, only: SZERO, DONE, MACHI_EPS
    use initial_module, only: newper, newper_inv, st_sim
    use read_input, only: len_scal_inv
    use check_condition, only: check_outf_cond
#ifdef MPI_MSG
    use mpi_set, only: bcast_calc_ftype, bcast_sim_val, bcast_out_type
#endif
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    ! -- Set grid information (grid_info)
      call set_grid_info()

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast calculation file type (calc_ftype)
        call bcast_calc_ftype()
      ! -- Bcast simulation value (sim_val)
        call bcast_sim_val()
    end if
#endif

    ! -- Set retention information (retn_info)
      call set_retn_info()

    ! -- Set parameter information (parm_info)
      call set_parm_info()

    ! -- Set geography information (geog_info)
      call set_geog_info()

    ! -- Set initial information (init_info)
      call set_init_info()

    ! -- Set massbalance information (mass_info)
      call set_mass_info()

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast output file type (out_type)
        call bcast_out_type()
    end if
#endif

    ! -- Check output file condition (outf_cond)
      call check_outf_cond()

    current_t = SZERO ; delt = DZERO ; delt_inv = DZERO
    newper = MACHI_EPS*len_scal_inv ; newper_inv = DONE/newper
    now_date(:) = st_sim%sta_date(:) ; out_iter = 0 ; conv_flag = 0 ; form_switch = 0

  end subroutine prepare_calc

  subroutine set_grid_info()
  !***************************************************************************************
  ! set_grid_info -- Set grid information
  !***************************************************************************************
    ! -- modules
    use initial_module, only: my_rank
    use initial_module, only: st_grid
    use read_input, only: read_grid_file
    use make_cell, only: make_cell_info
    ! -- inout

    ! -- local
    integer(I4) :: nx, ny, nz
    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      nx = st_grid%nx ; ny = st_grid%ny ; nz = st_grid%nz
      ! -- Read grid file (grid_file)
        call read_grid_file(st_grid%fnum, nx, ny, nz, st_in_type%grid)
    end if

    ! -- Make cell information (cell_info)
      call make_cell_info()

  end subroutine set_grid_info

  subroutine set_retn_info()
  !***************************************************************************************
  ! set_retn_info -- Set retention information
  !***************************************************************************************
    ! -- modules
    use open_file, only: open_in_retnf
    use check_condition, only: check_calc_retn
    use assign_calc, only: assign_retnv, read_reta, read_retn, read_resi
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
#ifdef MPI_MSG
    ! -- Open input retention file (in_retnf)
      call open_in_retnf(st_in_type%retn, st_in_path%retn, calc_r4view)
#else
    ! -- Open input retention file (in_retnf)
      call open_in_retnf(st_in_type%retn, st_in_path%retn)
#endif
    ! -- Assign retention value (retnv)
      call assign_retnv(st_in_type%retn)

    ! -- Check calculation retention value (calc_retn)
      call check_calc_retn(ncalc, read_reta, read_retn, read_resi)

  end subroutine set_retn_info

  subroutine set_parm_info()
  !***************************************************************************************
  ! set_parm_info -- Set parameter information
  !***************************************************************************************
    ! -- modules
    use open_file, only: open_in_parmf
    use check_condition, only: check_calc_parm
    use assign_calc, only: assign_parmv, read_ksx, read_ksy, read_ksz, read_ss, read_poro
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
#ifdef MPI_MSG
    ! -- Open input parameter file (in_parmf)
      call open_in_parmf(st_in_type%parm, st_in_path%parm, calc_r4view)
#else
    ! -- Open input parameter file (in_parmf)
      call open_in_parmf(st_in_type%parm, st_in_path%parm)
#endif
    ! -- Assign parameter value (parav)
      call assign_parmv(st_in_type%parm)

    ! -- Check calculation parameter value (calc_para)
      call check_calc_parm(ncalc, read_ksx, read_ksy, read_ksz, read_ss, read_poro)

  end subroutine set_parm_info

  subroutine set_geog_info()
  !***************************************************************************************
  ! set_geog_info -- Set geography information
  !***************************************************************************************
    ! -- modules
    use open_file, only: open_in_geogf
    use make_cell, only: surf_elev
    use set_cell, only: ncals
    use assign_calc, only: assign_geogv, surf_bott, surf_reli, geog_num
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    allocate(surf_bott(ncals), surf_top(ncals), surf_reli(ncals))
    !$omp parallel workshare
    surf_bott(:) = surf_elev(:) ; surf_top(:) = surf_elev(:) ; surf_reli(:) = DZERO
    !$omp end parallel workshare

    if (st_in_type%geog == in_type(0)) then
#ifdef MPI_MSG
      ! -- Open input geography file (in_geogf)
        call open_in_geogf(st_in_path%geog, cals_r4view)
#else
      ! -- Open input geography file (in_geogf)
        call open_in_geogf(st_in_path%geog)
#endif
      ! -- Assign geography value (geogv)
        call assign_geogv()
    end if

    if (geog_num /= 0) then
      !$omp parallel do private(i)
      do i = 1, ncals
        surf_top(i) = surf_bott(i) + surf_reli(i)
      end do
      !$omp end parallel do
    end if

  end subroutine set_geog_info

  subroutine set_init_info()
  !***************************************************************************************
  ! set_init_info -- Set initial information
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_in_unit
    use open_file, only: open_in_initf
    use check_condition, only: check_calc_init
    use assign_calc, only: assign_initv, read_init
#ifdef MPI_MSG
!    use mpi_utility, only: bcast_path_unit
    use mpi_utility, only: bcast_file
    use mpi_set, only: cals_r4hview, calc_r4hview, rest_view
#endif
    ! -- inout

    ! -- local
#ifdef MPI_MSG
    integer(I4) :: init_view
#endif
    !-------------------------------------------------------------------------------------
    if (st_in_type%init /= in_type(7)) then
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast file (file)
          call bcast_file(st_in_path%init, st_in_unit%init, "initial")
      end if

      if (st_in_type%init == in_type(4)) then
        if (len_trim(adjustl(st_in_unit%init)) == 0) then
          init_view = cals_r4view
        else
          init_view = cals_r4hview
        end if
      else if (st_in_type%init == in_type(6)) then
        if (len_trim(adjustl(st_in_unit%init)) == 0) then
          init_view = calc_r4view
        else
          init_view = calc_r4hview
        end if
      else if (st_in_type%init == in_type(0)) then
        init_view = rest_view
      else
        init_view = 0
      end if

      ! -- Open input initial file (in_initf)
        call open_in_initf(st_in_type%init, st_in_path%init, st_in_unit%init, init_view)
#else
      ! -- Open input initial file (in_initf)
        call open_in_initf(st_in_type%init, st_in_path%init, st_in_unit%init)
#endif

    end if

    ! -- Assign initial value (initv)
      call assign_initv(st_in_type%init, st_in_unit%init)

    ! -- Check calculation initial value (calc_init)
      call check_calc_init(ncalc, read_init)

  end subroutine set_init_info

  subroutine set_mass_info()
  !***************************************************************************************
  ! set_mass_info -- Set massbalance information
  !***************************************************************************************
    ! -- modules
    use open_file, only: open_in_massf
    use assign_calc, only: assign_massv
#ifdef MPI_MSG
!    use mpi_utility, only: bcast_file_path
    use mpi_utility, only: bcast_file
    use mpi_set, only: cals_i4view, calc_i4view
#endif
    ! -- inout

    ! -- local
    integer(I4), allocatable :: all_mass_type(:)
    logical, allocatable :: mass_mask(:)
    !-------------------------------------------------------------------------------------
    allocate(all_mass_type(6), mass_mask(6))
    all_mass_type(:) = [in_type(1), in_type(3:7)]
    mass_mask(:) = (st_in_type%mass == all_mass_type(:))

    if (any(mass_mask)) then
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast file (file)
          call bcast_file(st_in_path%mass, "input massbalance")
      end if

      if (st_in_type%mass == in_type(4)) then
        ! -- Open input massbalance file for output (in_massf)
          call open_in_massf(st_in_type%mass, st_in_path%mass, cals_i4view)
      else if (st_in_type%mass == in_type(6)) then
        ! -- Open input massbalance file for output (in_massf)
          call open_in_massf(st_in_type%mass, st_in_path%mass, calc_i4view)
      else
        ! -- Open input massbalance file for output (in_massf)
          call open_in_massf(st_in_type%mass, st_in_path%mass)
      end if
#else
      ! -- Open input massbalance file for output (in_massf)
        call open_in_massf(st_in_type%mass, st_in_path%mass)
#endif
      ! -- Assign massbalance value (massv)
        call assign_massv(st_in_type%mass)
    end if

  end subroutine set_mass_info

end module prep_calculation
