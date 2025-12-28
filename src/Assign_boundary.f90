module assign_boundary
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: DZERO, SNOVAL
  use initial_module, only: pro_totn, my_rank, in_type
  use utility_module, only: close_file
  use read_module, only: read_clasf
  use read_input, only: len_scal_inv
  use set_condition, only: set_clas2calc, set_2dfile2cals
#ifdef MPI_MSG
  use mpi_read, only: close_mpi_file
  use mpi_set, only: bcast_clas_val
#endif

  implicit none
  private
  public :: assign_sealv, assign_surfbv, assign_wellv, assign_rilav

  real(SP), allocatable, public :: read_seal(:), read_rech(:), read_well(:)
  real(SP), allocatable, public :: read_prec(:), read_evap(:)
  real(DP), allocatable, public :: well_top(:), well_bott(:), calc_well(:)
  integer(I4), allocatable, public :: rech_cflag(:), prec_cflag(:), evap_cflag(:)

  ! -- local

  contains

  subroutine assign_sealv(seal_ftype)
  !***************************************************************************************
  ! assign_sealv -- Assign sea level value
  !***************************************************************************************
    ! -- modules
    use utility_module, only: conv_i2s, write_err_stop
    use initial_module, only: st_seal, st_in_type
    use read_module, only: read_3dpointf
    use open_file, only: st_intse
    use set_cell, only: ncell, seal_cnum, no_ncals, no_ncalc
    use set_condition, only: set_point2seal, set_2dfile2seal, set_3dfile2seal,&
                             set_bound2calc
#ifdef MPI_MSG
    use mpi_utility, only: bcast_char
    use mpi_set, only: bcast_3dpoint
#endif
    ! -- inout
    integer(I4), intent(in) :: seal_ftype
    ! -- local
    integer(I4) :: sealn, intse_type, temp_ftype
    integer(I4), allocatable :: seal_cflag(:), seal2cell(:)
    integer(I4), allocatable :: seal_all_type(:)
    real(DP), allocatable :: cell_seal(:)
    character(:), allocatable :: err_mes, rank_str
    logical, allocatable :: seal_all_mask(:)
    !-------------------------------------------------------------------------------------
    if (st_seal%totn > 0) then
      allocate(read_seal(ncell), seal_cflag(ncell))
      allocate(seal_all_type(7), seal_all_mask(7))
      !$omp parallel workshare
      read_seal(:) = SNOVAL ; seal_cflag(:) = 0
      seal_all_type(:) = [in_type(1:7)]
      seal_all_mask(:) = (st_in_type%seal == seal_all_type(:))
      !$omp end parallel workshare

      if (seal_ftype == in_type(1)) then
        allocate(st_seal%value(st_seal%totn), st_seal%name(st_seal%totn))
        !$omp parallel workshare
        st_seal%value(:) = SNOVAL ; st_seal%name(:) = ""
        !$omp end parallel workshare
        if (my_rank == 0) then
          call read_clasf(st_seal%fnum, st_seal%totn, st_seal%name, st_seal%value)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          call bcast_clas_val(st_seal%totn, st_seal%name, st_seal%value)
        end if
#endif
        call set_clas2calc(st_seal%totn, st_seal%name, st_seal%value, read_seal,&
                           seal_cflag, sealn)
        deallocate(st_seal%value, st_seal%name)
      else if (seal_ftype == in_type(2)) then
        allocate(st_seal%value(st_seal%totn))
        allocate(st_seal%i(st_seal%totn), st_seal%j(st_seal%totn), st_seal%k(st_seal%totn))
        !$omp parallel workshare
        st_seal%value(:) = SNOVAL
        st_seal%i(:) = 0 ; st_seal%j(:) = 0 ; st_seal%k(:) = 0
        !$omp end parallel workshare
        if (my_rank == 0) then
          call read_3dpointf(st_seal%fnum, st_seal%totn, st_seal%i, st_seal%j, st_seal%k,&
                             st_seal%value)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          call bcast_3dpoint(st_seal%totn, st_seal%i, st_seal%j, st_seal%k, st_seal%value)
        end if
#endif
        call set_point2seal(st_seal%totn, st_seal%i, st_seal%j, st_seal%k, st_seal%value,&
                            read_seal, seal_cflag, sealn)
        deallocate(st_seal%i, st_seal%j, st_seal%k, st_seal%value)
      else
        if (seal_ftype == in_type(7)) then
          intse_type = st_intse%type
          temp_ftype = intse_type
        else
          temp_ftype = seal_ftype
        end if

        if (temp_ftype == in_type(3) .or. temp_ftype == in_type(4)) then
          call set_2dfile2seal(st_seal%fnum, seal_ftype, intse_type, no_ncals, SNOVAL,&
                               read_seal, seal_cflag, sealn)
        else if (temp_ftype == in_type(5) .or. temp_ftype == in_type(6)) then
          call set_3dfile2seal(st_seal%fnum, seal_ftype, intse_type, no_ncalc, SNOVAL,&
                               read_seal, seal_cflag, sealn)
        end if
        if (seal_ftype == in_type(7)) then
          if (intse_type == in_type(3) .or. intse_type == in_type(5)) then
            if (my_rank == 0) then
              call close_file(st_seal%fnum)
            end if
          else if (intse_type == in_type(4) .or. intse_type == in_type(6)) then
#ifdef MPI_MSG
            call close_mpi_file(st_seal%fnum)
#else
            call close_file(st_seal%fnum)
#endif
          end if
        end if
      end if

      if (seal_cnum /= sealn .or. my_rank == 0) then
        if (seal_cnum /= sealn) then
          rank_str = conv_i2s(my_rank)
          err_mes = "The number of sea grids is different in rank "//rank_str//"."
        end if
#ifdef MPI_MSG
        call bcast_char(err_mes)
#endif
        if (my_rank == 0 .and. allocated(err_mes)) then
          call write_err_stop(err_mes)
        end if
        if (allocated(err_mes)) then
          deallocate(err_mes)
        end if
        if (allocated(rank_str)) then
          deallocate(rank_str)
        end if
      end if

      if (any(seal_all_mask)) then
        allocate(seal2cell(sealn), cell_seal(sealn))
        !$omp parallel workshare
        seal2cell(:) = 0 ; cell_seal(:) = DZERO
        !$omp end parallel workshare
        call set_bound2calc(ncell, seal_cflag, read_seal, seal2cell, cell_seal)
        deallocate(read_seal)
        allocate(read_seal(sealn))
        !$omp parallel workshare
        read_seal(:) = real(cell_seal(:), kind=SP)*len_scal_inv
        !$omp end parallel workshare
        deallocate(seal_cflag, seal2cell, cell_seal)
      end if
      deallocate(seal_all_type, seal_all_mask)
    end if

  end subroutine assign_sealv

  subroutine assign_surfbv(sb_ftype, int_ftype, sb_st, sb_num, sb_cflag, read_sb)
  !***************************************************************************************
  ! assign_surfbv -- Assign surface boundary value
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_surfb
    ! -- inout
    integer(I4), intent(in) :: sb_ftype, int_ftype
    type(st_surfb), intent(inout) :: sb_st
    integer(I4), intent(out) :: sb_num
    integer(I4), intent(out) :: sb_cflag(:)
    real(SP), intent(out) :: read_sb(:)
    ! -- local

    !-------------------------------------------------------------------------------------
    sb_num = 0

    if (sb_st%totn > 0) then
      if (sb_ftype == in_type(1)) then
        allocate(sb_st%value(sb_st%totn), sb_st%name(sb_st%totn))
        !$omp parallel workshare
        sb_st%value(:) = SNOVAL ; sb_st%name(:) = ""
        !$omp end parallel workshare
        if (my_rank == 0) then
          call read_clasf(sb_st%fnum, sb_st%totn, sb_st%name, sb_st%value)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          call bcast_clas_val(sb_st%totn, sb_st%name, sb_st%value)
        end if
#endif
        call set_clas2calc(sb_st%totn, sb_st%name, sb_st%value, read_sb, sb_cflag, sb_num)
        deallocate(sb_st%value, sb_st%name)
      else
        call set_2dfile2cals(sb_st%fnum, sb_ftype, int_ftype, SNOVAL, read_sb, sb_cflag, sb_num)
        if (sb_ftype == in_type(7)) then
          if (int_ftype == in_type(3)) then
            if (my_rank == 0) then
              call close_file(sb_st%fnum)
            end if
          else if (int_ftype == in_type(4)) then
#ifdef MPI_MSG
            call close_mpi_file(sb_st%fnum)
#else
            call close_file(sb_st%fnum)
#endif
          end if
        end if
      end if

      !$omp parallel workshare
      read_sb(:) = read_sb(:)*len_scal_inv*sb_st%uni_conv
      !$omp end parallel workshare

    end if

  end subroutine assign_surfbv

  subroutine assign_wellv(well_ftype, weks_ftype, weke_ftype, num_well)
  !***************************************************************************************
  ! assign_wellv -- Assign well value
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_well
    use open_file, only: st_intwe
    use set_cell, only: ncalc
    use set_condition, only: set_point2well, set_2dwell, set_well2index, set_3dfile2well,&
                             set_well3d2index, set_wellprop
    ! -- inout
    integer(I4), intent(in) :: well_ftype, weks_ftype, weke_ftype
    integer(I4), intent(inout) :: num_well
    ! -- local
    integer(I4) :: intwe_type
    integer(I4), allocatable :: well2calc(:)
    integer(I4), allocatable :: type_2d(:), type_3d(:)
    !-------------------------------------------------------------------------------------
    if (st_well%totn > 0) then
      intwe_type = st_intwe%type
      allocate(type_2d(2), type_3d(2))
      type_2d(:) = [in_type(3), in_type(4)] ; type_3d(:) = [in_type(5), in_type(6)]
      if (any(well_ftype == type_3d) .or. any(intwe_type == type_3d)) then
        allocate(well2calc(ncalc))
        !$omp parallel workshare
        well2calc(:) = 0
        !$omp end parallel workshare
      end if

      if (well_ftype == in_type(2)) then
        ! -- Set well from point file (point2well)
          call set_point2well(st_well%fnum, num_well)
      else if (any(well_ftype == type_2d)) then
        ! -- Set well from 2d file (2dwell)
          call set_2dwell(st_well%fnum, well_ftype, 0, weks_ftype, weke_ftype, num_well)
      else if (any(well_ftype == type_3d)) then
        ! -- Set well value from 3d text file (3dfile2well)
          call set_3dfile2well(st_well%fnum, well_ftype, 0, num_well, well2calc)
      else if (well_ftype == in_type(7)) then
        if (any(intwe_type == type_2d)) then
          ! -- Set well from 2d file (2dwell)
            call set_2dwell(st_well%fnum, well_ftype, intwe_type, weks_ftype, weke_ftype, num_well)
        else if (any(intwe_type == type_3d)) then
          ! -- Set well value from 3d text file (3dfile2well)
            call set_3dfile2well(st_well%fnum, well_ftype, intwe_type, num_well, well2calc)
        end if
#ifdef MPI_MSG
        call close_mpi_file(st_well%fnum)
#else
        call close_file(st_well%fnum)
#endif
      end if

      if (well_ftype == in_type(2) .or. any(well_ftype == type_2d)) then
        ! -- Set well index (well2index)
          call set_well2index(num_well)
      else if (any(well_ftype == type_3d)) then
        ! -- Set well index from 3d well (well3d2index)
          call set_well3d2index(num_well, well2calc)
          deallocate(well2calc)
      else if (well_ftype == in_type(7)) then
        if (any(intwe_type == type_2d)) then
          ! -- Set well index (well2index)
            call set_well2index(num_well)
        else if (any(intwe_type == type_3d)) then
          ! -- Set well index from 3d well (well3d2index)
            call set_well3d2index(num_well, well2calc)
            deallocate(well2calc)
        end if
      end if

      if (well_ftype == in_type(2) .or. any(well_ftype == type_2d) .or. &
          any(well_ftype == type_3d) .or. well_ftype == in_type(7)) then
        allocate(well_top(num_well), well_bott(num_well), read_well(num_well))
        !$omp parallel workshare
        well_top(:) = DZERO ; well_bott(:) = DZERO ; read_well(:) = st_well%value(:)
        !$omp end parallel workshare
        call set_wellprop(num_well, well_top, well_bott)
        allocate(calc_well(ncalc))
        deallocate(st_well%value)
        !$omp parallel workshare
        calc_well(:) = DZERO
        read_well(:) = read_well(:)*(len_scal_inv**3)*st_well%uni_conv
        !$omp end parallel workshare
      end if
      deallocate(type_2d, type_3d)
    end if

  end subroutine assign_wellv

  subroutine assign_rilav(rl_ftype, lake_aflag, rl_st, rl_num, rl_cflag, calc_rl)
  !***************************************************************************************
  ! assign_rilav -- Assign river and lake value
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_surfw
    use read_module, only: read_2dpointf
    use set_condition, only: set_point2surf
#ifdef MPI_MSG
    use mpi_set, only: bcast_2dpoint
#endif
    ! -- inout
    integer(I4), intent(in) :: rl_ftype, lake_aflag
    type(st_surfw), intent(inout) :: rl_st
    integer(I4), intent(out) :: rl_num
    integer(I4), intent(out) :: rl_cflag(:)
    real(SP), intent(out) :: calc_rl(:)
    ! -- local
    real(SP) :: nodim_unit
    !-------------------------------------------------------------------------------------
    rl_num = 0

    if (rl_st%totn > 0) then
      if (rl_ftype == in_type(1)) then
        allocate(rl_st%value(rl_st%totn), rl_st%name(rl_st%totn))
        !$omp parallel workshare
        rl_st%value(:) = SNOVAL ; rl_st%name(:) = ""
        !$omp end parallel workshare
        if (my_rank == 0) then
          call read_clasf(rl_st%fnum, rl_st%totn, rl_st%name, rl_st%value)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          call bcast_clas_val(rl_st%totn, rl_st%name, rl_st%value)
        end if
#endif
        call set_clas2calc(rl_st%totn, rl_st%name, rl_st%value, calc_rl, rl_cflag, rl_num)
        deallocate(rl_st%value, rl_st%name)
      else if (rl_ftype == in_type(2)) then
        allocate(rl_st%value(rl_st%totn))
        allocate(rl_st%i(rl_st%totn), rl_st%j(rl_st%totn))
        !$omp parallel workshare
        rl_st%value(:) = SNOVAL
        rl_st%i(:) = 0 ; rl_st%j(:) = 0
        !$omp end parallel workshare
        if (my_rank == 0) then
          call read_2dpointf(rl_st%fnum, rl_st%totn, rl_st%i, rl_st%j, rl_st%value)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          call bcast_2dpoint(rl_st%totn, rl_st%i, rl_st%j, rl_st%value)
        end if
#endif
        call set_point2surf(rl_st%totn, rl_st%i, rl_st%j, rl_st%value, calc_rl, rl_cflag,&
                            rl_num)
        deallocate(rl_st%i, rl_st%j, rl_st%value)
      else
        call set_2dfile2cals(rl_st%fnum, rl_ftype, rl_st%inttype, SNOVAL, calc_rl,&
                             rl_cflag, rl_num)
        if (rl_ftype == in_type(7)) then
          if (rl_st%inttype == in_type(3)) then
            if (my_rank == 0) then
              call close_file(rl_st%fnum)
            end if
          else if (rl_st%inttype == in_type(4)) then
#ifdef MPI_MSG
            call close_mpi_file(rl_st%fnum)
#else
            call close_file(rl_st%fnum)
#endif
          end if
        end if
      end if

      if (lake_aflag /= 1) then
        nodim_unit = len_scal_inv
      else
        nodim_unit = len_scal_inv*len_scal_inv
      end if
      !$omp parallel workshare
      calc_rl(:) = calc_rl(:)*nodim_unit
      !$omp end parallel workshare

    end if

  end subroutine assign_rilav

end module assign_boundary
