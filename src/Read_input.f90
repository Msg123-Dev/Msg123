module read_input
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SZERO, SONE, DNOVAL, CHALEN, TIMELEN
  use utility_module, only: open_new_rtxt, close_file, write_logf, write_success,&
                            write_err_stop, get_days
  use initial_module, only: st_sim, st_in_type, st_in_path, st_in_unit, in_type, unit_list

  implicit none
  private
  public :: read_main_file, read_grid_file
  integer(I4), public :: temp_maxinn_iter
  real(SP), public :: len_scal, len_scal_inv
  real(DP), allocatable, public :: glob_x(:,:), glob_y(:,:), glob_z(:,:,:)

  ! -- local
  integer(I4) :: main_fnum
  integer(I4) :: sdate(6), edate(6)
  integer(I4) :: grid_check, grid_xnum, grid_ynum, grid_znum

  contains

  subroutine read_main_file()
  !***************************************************************************************
  ! read_main_file -- Read main file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    character(:), allocatable :: main_file
    !-------------------------------------------------------------------------------------
    main_file = 'msg123.main'

    ! -- Open new read text file (new_rtxt)
      call open_new_rtxt(1, 1, main_file, "msg123 main", main_fnum)

    ! -- Read simulation name list (sim_list)
      call read_sim_list()

    ! -- Read calculation time name list (calc_time_list)
      call read_calc_time_list()

    ! -- Read calculation region name list (calc_reg_list)
      call read_calc_reg_list()

    ! -- Read solution name list (sol_list)
      call read_sol_list()

    ! -- Read grid name list (grid_list)
      call read_grid_list()

    ! -- Read input setting (inp_set)
      call read_inp_set()

    ! -- Read output setting (out_set)
      call read_out_set()

  end subroutine read_main_file

  subroutine read_inp_set()
  !***************************************************************************************
  ! read_inp_set -- Read input setting
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    character(CHALEN) :: time_input_file
    namelist/set_time_input/time_input_file
    !-------------------------------------------------------------------------------------
    ! -- Read retention & parameter name list (retn_parm_list)
      call read_retn_parm_list()

    ! -- Read initial name list (init_list)
      call read_init_list()

    ierr = 0 ; time_input_file = ""
    read(unit=main_fnum,nml=set_time_input,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading timeseries input section in main file.")
    end if

    ! -- Read timeseries input name list (tinp_list)
      call read_tinp_list(trim(adjustl(time_input_file)))

    ! -- Read geography name list (geog_list)
      call read_geog_list()

    ! -- Read water table name list (wtab_list)
      call read_wtab_list()

    ! -- Read massbalance name list (mass_list)
      call read_mass_list()

  end subroutine read_inp_set

  subroutine read_out_set()
  !***************************************************************************************
  ! read_out_set -- Read output setting
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: len_char, tot_vari
    character(CHALEN) :: out_list
    character(1), allocatable :: temp_char(:)
    character(4), allocatable :: out_vari(:)
    namelist/set_out_vari/out_list
    !-------------------------------------------------------------------------------------
    ierr = 0 ; out_list = ""
    read(unit=main_fnum,nml=set_out_vari,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading output variable section in main file.")
    end if

    if (len_trim(adjustl(out_list)) /= 0) then
      len_char = len_trim(adjustl(out_list))

      allocate(temp_char(len_char))
      temp_char(:) = transfer(out_list, ' ', size = len_char)
      tot_vari = count(temp_char(:) == ",") + 1

      allocate(out_vari(tot_vari+3))
      ! -- Get output variables (out_vari)
        call get_out_vari(tot_vari, trim(adjustl(out_list)), out_vari)
      deallocate(temp_char)
    else
      tot_vari = 3
      allocate(out_vari(3))
      call write_logf("Output default variables, hydraulic head & restart.")
      out_vari(:) = [character(4) :: 'conv', 'head', 'rest']
    end if

    ! -- Set output (output)
      call set_output(tot_vari, out_vari)

    deallocate(out_vari)

  end subroutine read_out_set

  subroutine read_sim_list()
  !***************************************************************************************
  ! read_sim_list -- Read simulation name list
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: sim_type
    character(CHALEN) :: sim_name
    namelist/set_simulation/sim_type, sim_name
    !-------------------------------------------------------------------------------------
    ierr = 0 ; sim_type = -2 ; sim_name = ""
    read(unit=main_fnum,nml=set_simulation,iostat=ierr)

    if (ierr /= 0) then
      call write_err_stop("While reading simulation section in main file.")
    end if

    if (sim_type == -1) then
      call write_logf("Simulation type is steady-state.")
      sdate(:) = 0
    else if (sim_type == 0) then
      call write_logf("Simulation type is long-term-steady.")
      sdate(:) = 0
    else if (sim_type == 1) then
      call write_logf("Simulation type is transient-simulation.")
    else
      call write_err_stop("Specified wrong simulation type.")
    end if

    st_sim%sim_type = sim_type ; st_sim%sim_name = trim(adjustl(sim_name))

  end subroutine read_sim_list

  subroutine read_calc_time_list()
  !***************************************************************************************
  ! read_calc_time_list -- Read calculation time name list
  !***************************************************************************************
    ! -- modules
    use utility_module, only: conv_unit
    use initial_module, only: my_rank
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: stime_type
    real(SP) :: end_time, calc_multi
    character(TIMELEN) :: calc_unit
    namelist/set_calc_time/stime_type, sdate, edate, end_time, calc_unit
    !-------------------------------------------------------------------------------------
    ierr = 0 ; stime_type = -1 ; end_time = SZERO ; calc_multi = SZERO ; calc_unit = ""
    read(unit=main_fnum,nml=set_calc_time,iostat=ierr)

    if (ierr /= 0) then
      call write_err_stop("While reading calculation time section in main file.")
    end if

    if (len_trim(calc_unit) == 0) then
      call write_err_stop("Specify calculation time unit in main file.")
    else
      ! -- Convert unit (unit)
        call conv_unit(my_rank, calc_unit, "main file", sdate, calc_multi)
      if (st_sim%sim_type == 0) then
        end_time = end_time*calc_multi
      else if (st_sim%sim_type == 1) then
        ! -- Check calculation date (date)
          call check_date(sdate, edate)
        ! -- Calculate end time (etime)
          call calc_etime(sdate, edate, end_time)
      end if
    end if

    st_sim%res_type = stime_type ; st_sim%sta_date = sdate ; st_sim%end_date = edate
    st_sim%end_time = end_time ; st_sim%cal_fact = calc_multi
    st_sim%cal_unit = calc_unit

  end subroutine read_calc_time_list

  subroutine read_calc_reg_list()
  !***************************************************************************************
  ! read_calc_reg_list -- Read calculation region name list
  !***************************************************************************************
    ! -- modules
    use initial_module, only: noclas_flag
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: calc_type, calcreg_neib, clas_fnum
    character(CHALEN) :: calcreg_name, inact_name, clas_file
    namelist/set_calc_reg/calc_type, calcreg_neib, calcreg_name, inact_name, clas_file
    !-------------------------------------------------------------------------------------
    ierr = 0 ; calc_type = -1 ; calcreg_neib = -1
    calcreg_name = "" ; inact_name = "" ; clas_file = ""
    read(unit=main_fnum,nml=set_calc_reg,iostat=ierr)

    if (ierr /= 0) then
      call write_err_stop("While reading calculation region section in main file.")
    end if

    if (len_trim(adjustl(clas_file)) /= 0) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, trim(adjustl(clas_file)), "classification", clas_fnum)
      ! -- Read classification file (clas_file)
        call read_clas_file(clas_fnum)
    else
      call write_logf("Set not to use classification file.")
      noclas_flag = 1
      if (calc_type == in_type(1)) then
        call write_err_stop("Set to use classification file in calc_type.")
      end if
    end if

    st_sim%reg_type = calc_type ; st_sim%reg_neib = calcreg_neib
    st_sim%reg_name = trim(adjustl(calcreg_name))
    st_sim%inact_name = trim(adjustl(inact_name))

  end subroutine read_calc_reg_list

  subroutine read_clas_file(file_num)
  !***************************************************************************************
  ! read_clas_file -- Read classification file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_clas
    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    integer(I4) :: i, j, ierr
    integer(I4) :: max_clas
    !-------------------------------------------------------------------------------------
    ierr = 0
    read(unit=file_num,fmt=*,iostat=ierr) st_clas%totn

    if (ierr /= 0) then
      call write_err_stop("While reading header in classification file.")
    else if (st_clas%totn <= 0) then
      call write_err_stop("Class number is below 0. Please check the number.")
    end if

    allocate(st_clas%name(st_clas%totn), st_clas%num(st_clas%totn))
    !$omp parallel workshare
    st_clas%name(:) = "" ; st_clas%num(:) = 0
    !$omp end parallel workshare
    max_clas = 0

    do i = 1, st_clas%totn
      read(unit=file_num,fmt=*,iostat=ierr) st_clas%name(i), st_clas%num(i)
      if (ierr /= 0) then
        call write_err_stop("While reading header in classification file.")
      end if
      do j = 1, st_clas%num(i)
        read(unit=file_num,fmt=*,iostat=ierr)
        if (ierr /= 0) then
          call write_err_stop("While reading in classification file.")
        end if
      end do
      if (st_clas%num(i) > max_clas) then
        max_clas = st_clas%num(i)
      end if
    end do

    allocate(st_clas%i(max_clas,st_clas%totn), st_clas%j(max_clas,st_clas%totn))
    allocate(st_clas%k(max_clas,st_clas%totn))
    !$omp parallel workshare
    st_clas%i(:,:) = 0 ; st_clas%j(:,:) = 0 ; st_clas%k(:,:) = 0
    !$omp end parallel workshare

    rewind(unit=file_num,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reopening in classification file.")
    end if

    read(unit=file_num,fmt=*) st_clas%totn
    do j = 1, st_clas%totn
      read(unit=file_num,fmt=*) st_clas%name(j), st_clas%num(j)
      do i = 1, st_clas%num(j)
        read(unit=file_num,fmt=*) st_clas%i(i,j), st_clas%j(i,j), st_clas%k(i,j)
      end do
    end do

    call write_success("Read classification file", file_num)

    call close_file(file_num)

  end subroutine read_clas_file

  subroutine read_sol_list()
  !***************************************************************************************
  ! read_sol_list -- Read solution name list
  !***************************************************************************************
    ! -- modules
    use initial_module, only: maxout_iter, maxinn_iter, criteria, precon_type, nlevel
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    real(SP) :: init_step, incr_multi, max_step
    namelist/set_solution/init_step, incr_multi, max_step, maxout_iter, criteria,&
                          maxinn_iter, precon_type
    !-------------------------------------------------------------------------------------
    ierr = 0 ; init_step = SZERO ; incr_multi = SZERO ; max_step = SZERO
    read(unit=main_fnum,nml=set_solution,iostat=ierr)

    if (ierr /= 0) then
      call write_err_stop("While reading solution section in main file.")
    else if (maxout_iter <= 0) then
      call write_err_stop("Input a positive value for maximum number of outer iteration.")
    else if (maxinn_iter <= 0) then
      call write_err_stop("Input a positive value for maximum number of inner iteration.")
    end if
    temp_maxinn_iter = maxinn_iter

    nlevel = 1
    if (precon_type /= 0 .and. precon_type /= 1) then
      call write_err_stop("Input a valid value for preconditoner type.")
    end if

    st_sim%ini_step = init_step*st_sim%cal_fact
    st_sim%max_step = max_step*st_sim%cal_fact
    st_sim%inc_fact = incr_multi

    if (precon_type == 1) then
      call read_amg_parm(main_fnum)
    end if

  end subroutine read_sol_list

  subroutine read_amg_parm(file_num)
  !***************************************************************************************
  ! read_amg_parm -- Read amg parameter
  !***************************************************************************************
    ! -- modules
    use initial_module, only: amg_nlevel, maxvcy_iter, max_sweep, jac_omega, amg_theta
    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    integer(I4) :: ierr
    namelist/set_amg/amg_nlevel, maxvcy_iter, max_sweep, jac_omega, amg_theta
    !-------------------------------------------------------------------------------------
    ierr = 0
    read(unit=file_num,nml=set_amg,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading amg section in main file.")
    end if

  end subroutine read_amg_parm

  subroutine read_grid_list()
  !***************************************************************************************
  ! read_grid_list -- Read grid name list
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_grid
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: gridx, gridy, gridz, gridxyz, grid_type
    character(CHALEN) :: grid_file
    namelist/set_grid/gridx, gridy, gridz, grid_type, grid_file
    !-------------------------------------------------------------------------------------
    ierr = 0
    gridx = 0 ; gridy = 0 ; gridz = 0 ; gridxyz = 0 ; grid_type = -1 ; grid_file = ""
    read(unit=main_fnum,nml=set_grid,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading grid section in main file.")
    else if (gridx <= 0 .or. gridy <= 0 .or. gridz <= 0) then
      call write_err_stop("Input a positive value for grid number.")
    end if

    if (len_trim(adjustl(grid_file)) /= 0) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, trim(adjustl(grid_file)), "grid", st_grid%fnum)
    else
      call write_err_stop("Specify grid file path in main file.")
    end if

    st_grid%nx = gridx ; st_grid%ny = gridy ; st_grid%nz = gridz
    gridxyz = gridx*gridy*gridz ; st_grid%nxyz = gridxyz
    st_in_type%grid = grid_type ; st_in_path%grid = trim(adjustl(grid_file))

  end subroutine read_grid_list

  subroutine read_grid_file(file_num, nx, ny, nz, gtype)
  !***************************************************************************************
  ! read_grid_file -- Read grid file
  !***************************************************************************************
    ! -- modules
    use constval_module, only: DZERO
    use utility_module, only: open_new_rbin
    use read_module, only: read_2dtxt, read_2dbin, read_3dtxt, read_3dbin
    ! -- inout
    integer(I4), intent(in) :: file_num, nx, ny, nz, gtype
    ! -- local
    integer(I4) :: i, j, k, ierr
    integer(I4) :: grid_xtype, grid_ytype, grid_ztype
    integer(I4) :: xyz_ftype(3), xyz_fnum(3)
    real(SP) :: xcorner, ycorner, zcorner
    real(SP) :: xwidth, ywidth, zwidth
    real(SP), allocatable :: read_xy(:,:), read_z(:,:,:)
    character(1) :: xyz_mess(3)
    character(CHALEN) :: in_xpath, in_ypath, in_zpath
    character(:), allocatable :: xyz_path, err_mes
    namelist/ingrid_type/grid_xtype, grid_ytype, grid_ztype
    namelist/ingrid_path/in_xpath, in_ypath, in_zpath
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(glob_x(nx+1,ny+1), glob_y(nx+1,ny+1))
    allocate(glob_z(nx+1,ny+1,nz+1))
    !$omp parallel workshare
    glob_x(:,:) = DZERO ; glob_y(:,:) = DZERO ; glob_z(:,:,:) = DZERO
    !$omp end parallel workshare

    if (gtype == in_type(0)) then
      read(unit=file_num,nml=ingrid_type,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading grid type section in grid compile file.")
      end if

      in_xpath = "" ; in_ypath = "" ; in_zpath = ""
      read(unit=file_num,nml=ingrid_path,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading grid file path section in grid compile file.")
      end if

      xyz_ftype(1) = grid_xtype ; xyz_ftype(2) = grid_ytype ; xyz_ftype(3) = grid_ztype
      xyz_mess(1) = "x" ; xyz_mess(2) = "y" ; xyz_mess(3) = "z"

      allocate(read_xy(nx+1,ny+1))
      allocate(read_z(nx+1,ny+1,nz+1))
      !$omp parallel workshare
      read_xy(:,:) = SZERO ; read_z(:,:,:) = SZERO
      !$omp end parallel workshare

      do i = 1, 3
        select case (i)
        case(1)
          xyz_path = trim(adjustl(in_xpath))
        case(2)
          xyz_path = trim(adjustl(in_ypath))
        case(3)
          xyz_path = trim(adjustl(in_zpath))
        end select

        if (xyz_ftype(i) == in_type(3)) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, xyz_path, "grid "//xyz_mess(i), xyz_fnum(i))
          if (i == 1 .or. i == 2) then
            ! -- Read 2D text file (2dtxt)
              call read_2dtxt(xyz_fnum(i), nx+1, ny+1, read_xy)
          else
            call write_err_stop("Specified wrong file number in gridz type.")
          end if
        else if (xyz_ftype(i) == in_type(4)) then
          ! -- Open new read binary file (new_btxt)
            call open_new_rbin(1, 1, xyz_path, "grid "//xyz_mess(i), xyz_fnum(i))
          if (i == 1 .or. i == 2) then
            ! -- Read 2D binary file (2dbin)
              call read_2dbin(xyz_fnum(i), nx+1, ny+1, SZERO, read_xy)
          else
            call write_err_stop("Specified wrong file number in gridz type.")
          end if
        else if (xyz_ftype(i) == in_type(5)) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, xyz_path, "grid "//xyz_mess(i), xyz_fnum(i))
          if (i == 3) then
            ! -- Read 3D text file (3dtxt)
              call read_3dtxt(xyz_fnum(i), nx+1, ny+1, nz+1, read_z)
          else
            call write_err_stop("Specified wrong file number in gridx or gridy type.")
          end if
        else if (xyz_ftype(i) == in_type(6)) then
          ! -- Open new read binary file (new_btxt)
            call open_new_rbin(1, 1, xyz_path, "grid "//xyz_mess(i), xyz_fnum(i))
          if (i == 3) then
            ! -- Read 3D binary file (3dbin)
              call read_3dbin(xyz_fnum(i), nx+1, ny+1, nz+1, SZERO, read_z)
          else
            call write_err_stop("Specified wrong file number in gridx or gridy type.")
          end if
        end if
        call close_file(xyz_fnum(i))
        deallocate(xyz_path)

        select case (i)
        case(1)
          glob_x = real(read_xy, kind=DP)
        case(2)
          glob_y = real(read_xy, kind=DP)
        case(3)
          glob_z = real(read_z, kind=DP)
        end select
      end do

      deallocate(read_xy, read_z)

    else if (gtype == in_type(7)) then
      read(unit=file_num,fmt=*,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading header in grid file.")
      end if
      read(unit=file_num,fmt=*,iostat=ierr) xcorner, xwidth
      if (ierr /= 0) then
        call write_err_stop("While reading x value in grid file.")
      end if
      read(unit=file_num,fmt=*,iostat=ierr) ycorner, ywidth
      if (ierr /= 0) then
        call write_err_stop("While reading y value in grid file.")
      end if
      read(unit=file_num,fmt=*,iostat=ierr) zcorner, zwidth
      if (ierr /= 0) then
        call write_err_stop("While reading z value in grid file.")
      end if

      !$omp parallel
      !$omp do private(i)
      do i = 1, nx+1
        glob_x(i,:) = real(xcorner + xwidth*(i-1), kind=DP)
      end do
      !$omp end do

      !$omp do private(j)
      do j = 1, ny+1
        glob_y(:,j) = real(ycorner - ywidth*(j-1), kind=DP)
      end do
      !$omp end do

      !$omp do private(k)
      do k = 1, nz+1
        glob_z(:,:,k) = real(zcorner + zwidth*(k-1), kind=DP)
      end do
      !$omp end do
      !$omp end parallel
    else
      call write_err_stop("Specify correct number for grid_type.")
    end if

    ! -- Check grid location (gridloc)
      call check_grid_loc(nx, ny, nz)

    err_mes = ""

    if (grid_check == 0) then
      call write_success("Read grid file", file_num)
    else
      if (grid_xnum /= 0) then
        write(err_mes,'(a,i0,a)') "X direction is reversed at ", grid_xnum, " point."
        call write_logf(err_mes)
      end if
      if (grid_ynum /= 0) then
        write(err_mes,'(a,i0,a)') "Y direction is reversed at ", grid_ynum, " point."
        call write_logf(err_mes)
      end if
      if (grid_znum /= 0) then
        write(err_mes,'(a,i0,a)') "Z direction is reversed at ", grid_znum, " point."
        call write_logf(err_mes)
      end if
      call write_err_stop("Check the direction in grid file.")
    end if

    deallocate(err_mes)

    call close_file(file_num)

  end subroutine read_grid_file

  subroutine read_retn_parm_list()
  !***************************************************************************************
  ! read_retn_parm_list -- Read retention & parameter name list
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: retn_type, parm_type
    integer(I4), allocatable :: all_retn_parm_type(:)
    character(CHALEN) :: retn_file, parm_file
    logical, allocatable :: retn_mask(:), parm_mask(:)
    namelist/set_retn_parm/retn_type, parm_type, retn_file, parm_file
    !-------------------------------------------------------------------------------------
    ierr = 0 ; retn_type = -1 ; parm_type = -1 ; retn_file = "" ; parm_file = ""
    read(unit=main_fnum,nml=set_retn_parm,iostat=ierr)

    if (ierr /= 0) then
      call write_err_stop("While reading retention & parameter section in main file.")
    end if

    allocate(all_retn_parm_type(2), retn_mask(2), parm_mask(2))
    all_retn_parm_type(:) = in_type(0:1)
    retn_mask(:) = (retn_type /= all_retn_parm_type(:))
    parm_mask(:) = (parm_type /= all_retn_parm_type(:))

    if (all(retn_mask)) then
      call write_err_stop("Specify correct number for retention file type.")
    else if (all(parm_mask)) then
      call write_err_stop("Specify correct number for parameter file type.")
    else if (len_trim(adjustl(retn_file)) == 0) then
      call write_err_stop("Specify retention file path in main file.")
    else if (len_trim(adjustl(parm_file)) == 0) then
      call write_err_stop("Specify parameter file path in main file.")
    end if

    st_in_type%retn = retn_type ; st_in_type%parm = parm_type
    st_in_path%retn = trim(adjustl(retn_file))
    st_in_path%parm = trim(adjustl(parm_file))

    deallocate(all_retn_parm_type)
    deallocate(retn_mask, parm_mask)

  end subroutine read_retn_parm_list

  subroutine read_init_list()
  !***************************************************************************************
  ! read_init_list -- Read initial name list
  !***************************************************************************************
    ! -- modules
    use constval_module, only: SNOVAL
    use initial_module, only: st_init
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: init_type
    integer(I4), allocatable :: all_init_type(:)
    real(SP) :: init_dept
    character(CHALEN) :: init_file
    character(TIMELEN) :: init_unit
    logical, allocatable :: init_mask(:)
    namelist/set_init/init_type, init_file, init_dept, init_unit
    !-------------------------------------------------------------------------------------
    ierr = 0 ; init_type = -1 ; init_dept = DNOVAL ; init_file = "" ; init_unit = ""
    read(unit=main_fnum,nml=set_init,iostat=ierr)

    if (ierr /= 0) then
      call write_err_stop("While reading initial condition section in main file.")
    end if

    allocate(all_init_type(6), init_mask(6))
    all_init_type(:) = [in_type(0), in_type(3:7)]
    init_mask(:) = (init_type /= all_init_type(:))

    if (all(init_mask)) then
      call write_err_stop("Specify correct number for initial file type.")
    else if (len_trim(adjustl(init_file)) == 0 .and. init_type /= in_type(7)) then
      call write_err_stop("Specify initial file path in main file.")
    else if (init_type == in_type(7) .and. init_dept == SNOVAL) then
      call write_err_stop("Specify depth from surface elevation in main file.")
    else if (init_type == in_type(0) .and. st_sim%res_type == 1) then
      if (all(unit_list(:) /= init_unit)) then
        call write_err_stop("Specify initial unit in main file.")
      end if
    end if

    st_in_type%init = init_type ; st_in_path%init = trim(adjustl(init_file))
    st_init%depth = init_dept ; st_in_unit%init = trim(adjustl(init_unit))

    deallocate(all_init_type)
    deallocate(init_mask)

  end subroutine read_init_list

  subroutine read_tinp_list(file_path)
  !***************************************************************************************
  ! read_tinp_list -- Read timeseries input name list
  !***************************************************************************************
    ! -- modules

    ! -- inout
    character(*), intent(in) :: file_path
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: tinp_fnum
    integer(I4) :: seal_type, rech_type, well_type, weks_type, weke_type
    integer(I4) :: rive_type, lake_type, prec_type, evap_type
    integer(I4), allocatable :: all_tinp_type(:)
    character(CHALEN) :: seal_path, rech_path, well_path, weks_path, weke_path
    character(CHALEN) :: rive_path, lake_path, prec_path, evap_path
    character(TIMELEN) :: seal_unit, rech_unit, well_unit, prec_unit, evap_unit
    namelist/tinp_type/seal_type, rech_type, well_type, weks_type, weke_type,&
                       rive_type, lake_type, prec_type, evap_type
    namelist/tinp_path/seal_path, rech_path, well_path, weks_path, weke_path,&
                       rive_path, lake_path, prec_path, evap_path
    namelist/tinp_unit/seal_unit, rech_unit, well_unit, prec_unit, evap_unit
    !-------------------------------------------------------------------------------------
    ! -- Open new read text file (new_rtxt)
      call open_new_rtxt(1, 1, file_path, "timeseries input", tinp_fnum)

    ierr = 0
    seal_type = -1 ; rech_type = -1 ; well_type = -1 ; weks_type = -1 ; weke_type = -1
    rive_type = -1 ; lake_type = -1 ; prec_type = -1 ; evap_type = -1
    read(unit=tinp_fnum,nml=tinp_type,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading file type section in timeseries input file.")
    end if

    seal_path = "" ; rech_path = "" ; well_path = "" ; weks_path = "" ; weke_path = ""
    rive_path = "" ; lake_path = "" ; prec_path = "" ; evap_path = ""
    read(unit=tinp_fnum,nml=tinp_path,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading file path section in timeseries input file.")
    end if

    seal_unit = "" ; rech_unit = "" ; well_unit = "" ; prec_unit = "" ; evap_unit = ""
    read(unit=tinp_fnum,nml=tinp_unit,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading file unit section in timeseries input file.")
    end if

    allocate(all_tinp_type(7))
    all_tinp_type(:) = [in_type(1:7)]

    if (seal_type > 0) then
      ! -- Check timeseries input file (tinp)
        call check_tinp(seal_type, all_tinp_type, trim(adjustl(seal_path)), "sea level", seal_unit)
    else
      call write_logf("Set not to use sea level boundary.")
    end if

    st_in_type%seal = seal_type ; st_in_path%seal = trim(adjustl(seal_path))
    st_in_unit%seal = seal_unit
    deallocate(all_tinp_type)

    allocate(all_tinp_type(4))
    all_tinp_type(:) = [in_type(1), in_type(3:4), in_type(7)]

    if (rech_type > 0) then
      ! -- Check timeseries input file (tinp)
        call check_tinp(rech_type, all_tinp_type, trim(adjustl(rech_path)), "recharge", rech_unit)
    else
      call write_logf("Set not to use recharge boundary.")
    end if

    st_in_type%rech = rech_type ; st_in_path%rech = trim(adjustl(rech_path))
    st_in_unit%rech = rech_unit
    deallocate(all_tinp_type)

    allocate(all_tinp_type(6))
    all_tinp_type(:) = [in_type(2:7)]

    if (well_type > 0) then
    ! -- Check timeseries input file (tinp)
      call check_tinp(well_type, all_tinp_type, trim(adjustl(well_path)), "well", well_unit)
    else
      call write_logf("Set not to use well boundary.")
    end if

    st_in_type%well = well_type ; st_in_path%well = trim(adjustl(well_path))
    st_in_unit%well = well_unit
    deallocate(all_tinp_type)

    if (well_type == in_type(3) .or. well_type == in_type(4)) then
      if (weks_type /= in_type(3) .and. weks_type /= in_type(4)) then
        call write_err_stop("Specify correct number for well start in timeseries input file.")
      else if (weke_type /= in_type(3) .and. weke_type /= in_type(4)) then
        call write_err_stop("Specify correct number for well end in timeseries input file.")
      else
        allocate(all_tinp_type(2))
        all_tinp_type(:) = [in_type(3:4)]
        ! -- Check timeseries input file (tinp)
          call check_tinp(weks_type, all_tinp_type, trim(adjustl(weks_path)), "well start")

        st_in_type%weks = weks_type ; st_in_path%weks = trim(adjustl(weks_path))

        ! -- Check timeseries input file (tinp)
          call check_tinp(weke_type, all_tinp_type, trim(adjustl(weke_path)), "well end")

        st_in_type%weke = weke_type ; st_in_path%weke = trim(adjustl(weke_path))
        deallocate(all_tinp_type)
      end if
    else
      if (weks_type > 0) then
        call write_logf("Ignored well start in timeseries input file.")
      else if (weke_type > 0) then
        call write_logf("Ignored well end in timeseries input file.")
      end if
    end if

    allocate(all_tinp_type(1))
    all_tinp_type(:) = [in_type(0)]
    if (rive_type >= 0) then
      ! -- Check timeseries input file (tinp)
        call check_tinp(rive_type, all_tinp_type, trim(adjustl(rive_path)), "river")
    else
      call write_logf("Set not to use river water boundary.")
    end if

    st_in_type%rive = rive_type ; st_in_path%rive = trim(adjustl(rive_path))
    deallocate(all_tinp_type)

    allocate(all_tinp_type(1))
    all_tinp_type(:) = [in_type(0)]
    if (lake_type >= 0) then
      ! -- Check timeseries input file (tinp)
        call check_tinp(lake_type, all_tinp_type, trim(adjustl(lake_path)), "lake")
    else
      call write_logf("Set not to use lake water boundary.")
    end if

    st_in_type%lake = lake_type ; st_in_path%lake = trim(adjustl(lake_path))
    deallocate(all_tinp_type)

    allocate(all_tinp_type(4))
    all_tinp_type(:) = [in_type(1), in_type(3:4), in_type(7)]
    if (prec_type > 0) then
      ! -- Check timeseries input file (tinp)
        call check_tinp(prec_type, all_tinp_type, trim(adjustl(prec_path)), "precipitation", prec_unit)
    else
      call write_logf("Set not to use precipitation boundary.")
    end if

    st_in_type%prec = prec_type ; st_in_path%prec = trim(adjustl(prec_path))
    st_in_unit%prec = prec_unit
    deallocate(all_tinp_type)

    allocate(all_tinp_type(4))
    all_tinp_type(:) = [in_type(1), in_type(3:4), in_type(7)]
    if (evap_type > 0) then
      ! -- Check timeseries input file (tinp)
        call check_tinp(evap_type, all_tinp_type, trim(adjustl(evap_path)), "evapotranspiration", evap_unit)
    else
      call write_logf("Set not to use evapotranspiration boundary.")
    end if

    st_in_type%evap = evap_type ; st_in_path%evap = trim(adjustl(evap_path))
    st_in_unit%evap = evap_unit

    deallocate(all_tinp_type)

  end subroutine read_tinp_list

  subroutine read_geog_list()
  !***************************************************************************************
  ! read_geog_list -- Read geography name list
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: geog_type
    character(CHALEN) :: geog_file
    namelist/set_geog/geog_type, geog_file
    !-------------------------------------------------------------------------------------
    ierr = 0 ; geog_type = -1 ; geog_file = ""
    read(unit=main_fnum,nml=set_geog,iostat=ierr)

    if (geog_type /= in_type(0)) then
      call write_logf("Set not to use geography function.")
    else
      if (len_trim(adjustl(geog_file)) == 0) then
        call write_err_stop("Specify geography file path in main file.")
      end if

      st_in_type%geog = geog_type ; st_in_path%geog = trim(adjustl(geog_file))
    end if

  end subroutine read_geog_list

  subroutine read_wtab_list()
  !***************************************************************************************
  ! read_wtab_list -- Read water table name list
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: wtab_type
    namelist/set_wtab/wtab_type
    !-------------------------------------------------------------------------------------
    ierr = 0 ; wtab_type = -1
    read(unit=main_fnum,nml=set_wtab,iostat=ierr)

    if (wtab_type /= in_type(7)) then
      call write_logf("Set not to use water table depth output function.")
    end if

    st_in_type%wtab = wtab_type

  end subroutine read_wtab_list

  subroutine read_mass_list()
  !***************************************************************************************
  ! read_mass_list -- Read massbalance name list
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: mass_type
    integer(I4), allocatable :: all_mass_type(:)
    character(CHALEN) :: mass_file
    logical, allocatable :: mass_mask(:)
    namelist/set_mass/mass_type, mass_file
    !-------------------------------------------------------------------------------------
    ierr = 0 ; mass_type = -1 ; mass_file = ""
    allocate(all_mass_type(4), mass_mask(4))
    all_mass_type(:) = [in_type(3:6)]
    mass_mask(:) = (mass_type == all_mass_type(:))

    read(unit=main_fnum,nml=set_mass,iostat=ierr)

    if (any(mass_mask) .and. len_trim(adjustl(mass_file)) == 0) then
      call write_err_stop("Specify massbalance file path in main file.")
    end if

    st_in_type%mass = mass_type ; st_in_path%mass = trim(adjustl(mass_file))

    deallocate(all_mass_type)
    deallocate(mass_mask)

  end subroutine read_mass_list

  subroutine check_date(stad, endd)
  !***************************************************************************************
  ! check_date -- Check calculation date
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: stad(:), endd(:)
    ! -- local
    integer(I4) :: staday, endday
    logical :: date_flag = .false.
    !-------------------------------------------------------------------------------------
    staday = get_days(stad(1), stad(2)) ; endday = get_days(endd(1), endd(2))
    if (stad(1) <= 0 .or. endd(1) <= 0) then
      date_flag = .true.
    else if (stad(1) > endd(1)) then
      date_flag = .true.
    else if (stad(2) < 1 .or. stad(2) > 12 .or. endd(2) < 1 .or. endd(2) > 12) then
      date_flag = .true.
    else if (stad(3) < 1 .or. stad(3) > staday .or. endd(2) < 1 .or. endd(2) > endday) then
      date_flag = .true.
    else if (stad(4) < 0 .or. stad(4) >= 24 .or. endd(4) < 0 .or. endd(4) >= 24) then
      date_flag = .true.
    else if (stad(5) < 0 .or. stad(5) >= 60 .or. endd(5) < 0 .or. endd(5) >= 60) then
      date_flag = .true.
    else if (stad(6) < 0 .or. stad(6) >= 60 .or. endd(6) < 0 .or. endd(6) >= 60) then
      date_flag = .true.
    else if (stad(1) == endd(1)) then
      if (stad(2) > endd(2)) then
        date_flag = .true.
      else if (stad(2) == endd(2)) then
        if (stad(3) > endd(3)) then
          date_flag = .true.
        else if (stad(3) == endd(3)) then
          if (stad(4) > endd(4)) then
            date_flag = .true.
          else if (stad(4) == endd(4)) then
            if (stad(5) > endd(5)) then
              date_flag = .true.
            else if (stad(5) == endd(5)) then
              if (stad(6) >= endd(6)) then
                date_flag = .true.
              end if
            end if
          end if
        end if
      end if
    end if

    if (date_flag) then
      call write_err_stop("Specified wrong date. Please check the time unit.")
    end if

  end subroutine check_date

  subroutine check_grid_loc(nx, ny, nz)
  !***************************************************************************************
  ! check_grid_loc -- Check grid location
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: nx, ny, nz
    ! -- local
    integer(I4) :: i, j, k
    integer(I4), allocatable :: gxnum(:,:), gynum(:,:), gznum(:,:,:)
    real(DP) :: minz, maxz, real_scal
    !-------------------------------------------------------------------------------------
    grid_check = 0 ; grid_xnum = 0 ; grid_ynum = 0 ; grid_znum = 0
    allocate(gxnum(nx+1,ny+1), gynum(nx+1,ny+1))
    allocate(gznum(nx+1,ny+1,nz+1))
    !$omp parallel
    !$omp workshare
    gxnum(:,:) = 0 ; gynum(:,:) = 0 ; gznum(:,:,:) = 0
    !$omp end workshare

    !$omp do private(i, j)
    do j = 1, ny+1
      do i = 2, nx+1
        if (glob_x(i,j) < glob_x(i-1,j)) then
          gxnum(i,j) = 1
        end if
      end do
    end do
    !$omp end do

    !$omp do private(i, j)
    do j = 2, ny+1
      do i = 1, nx+1
        if (glob_y(i,j) > glob_y(i,j-1)) then
          gynum(i,j) = 1
        end if
      end do
    end do
    !$omp end do

    !$omp do private(i, j, k)
    do k = 2, nz+1
      do j = 1, ny+1
        do i = 1, nx+1
          if (glob_z(i,j,k) > glob_z(i,j,k-1)) then
            gznum(i,j,k) = 1
          end if
        end do
      end do
    end do
    !$omp end do

    minz = -DNOVAL ; maxz = DNOVAL
    !$omp do private(i, j, k) reduction(min:minz) reduction(max:maxz)
    do k = 1, nz+1
      do j = 1, ny+1
        do i = 1, nx+1
          if (glob_z(i,j,k) < minz) then
            minz = glob_z(i,j,k)
          end if
          if (glob_z(i,j,k) > maxz) then
            maxz = glob_z(i,j,k)
          end if
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    grid_xnum = sum(gxnum) ; grid_ynum = sum(gynum) ; grid_znum = sum(gznum)
    deallocate(gxnum, gynum, gznum)

    grid_check = grid_xnum + grid_ynum + grid_znum

    real_scal = log10(max(abs(maxz),abs(minz)))
    if ((real_scal - int(real_scal)) > SZERO) then
      len_scal = 10**(int(real_scal)+1)
    else
      len_scal = 10**(int(real_scal))
    end if

    len_scal_inv = SONE/len_scal

    !$omp parallel workshare
    glob_x(:,:) = glob_x(:,:)*len_scal_inv ; glob_y(:,:) = glob_y(:,:)*len_scal_inv
    glob_z(:,:,:) = glob_z(:,:,:)*len_scal_inv
    !$omp end parallel workshare

  end subroutine check_grid_loc

  subroutine check_tinp(tinp_type, all_type, tinp_path, tinp_name, tinp_unit)
  !***************************************************************************************
  ! check_tinp -- Check timeseries input file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: tinp_type
    integer(I4), intent(in) :: all_type(:)
    character(*), intent(in) :: tinp_path, tinp_name
    character(*), intent(in), optional :: tinp_unit
    ! -- local
    character(:), allocatable :: err_mes
    logical, allocatable :: tinp_mask(:)
    !-------------------------------------------------------------------------------------
    allocate(tinp_mask(size(all_type)))
    !$omp parallel workshare
    tinp_mask(:) = (tinp_type /= all_type(:))
    !$omp end parallel workshare

    if (all(tinp_mask)) then
      err_mes = "Specify correct number for "//tinp_name//" in timeseries input file."
      call write_err_stop(err_mes)
    else if (len_trim(tinp_path) == 0) then
      err_mes = "Specify "//tinp_name//" file path in timeseries input file."
      call write_err_stop(err_mes)
    else if (present(tinp_unit)) then
      if (all(unit_list(:) /= tinp_unit) .and. len_trim(tinp_unit) /= 0) then
        err_mes = "Specify correct "//tinp_name//" unit in timeseries input file."
        call write_err_stop(err_mes)
      end if
    end if

    if (allocated(err_mes)) then
      deallocate(err_mes)
    end if

    deallocate(tinp_mask)

  end subroutine check_tinp

  subroutine calc_etime(stad, endd, etime)
  !***************************************************************************************
  ! calc_etime -- Calculate end time
  !***************************************************************************************
    ! -- modules
    use constval_module, only: MINSEC, HOURSEC, DAYSEC
    ! -- inout
    integer(I4), intent(inout) :: stad(:), endd(:)
    real(SP), intent(out) :: etime
    ! -- local
    integer(I4) :: stayear, stamonth, endmonth, temp_day
    !-------------------------------------------------------------------------------------
    etime = SZERO
    ! plus second
    etime = etime + (endd(6)-stad(6))*SONE
    if (endd(6) < stad(6)) then
      etime = etime + MINSEC ; endd(5) = endd(5) - 1
      if (endd(5) < 0) then
        endd(5) = 59 ; endd(4) = endd(4) - 1
      end if
    end if
    !plus minute
    etime = etime + (endd(5)-stad(5))*MINSEC
    if (endd(5) < stad(5)) then
      etime = etime + HOURSEC ; endd(4) = endd(4) - 1
      if (endd(4) < 0) then
        endd(4) = 23 ; endd(3) = endd(3) - 1
      end if
    end if
    !plus hour
    etime = etime + (endd(4)-stad(4))*HOURSEC
    if (endd(4) < stad(4)) then
      etime = etime + DAYSEC ; endd(3) = endd(3) - 1
    end if
    !plus day
    if (endd(3) == 0 .or. endd(3) < stad(3)) then
      endd(2) = endd(2) - 1
      if (endd(2) == 0) then
        endd(2) = 12 ; endd(1) = endd(1) - 1
      end if
      temp_day = get_days(endd(1), endd(2)) + endd(3)
    else
      temp_day = endd(3)
    end if
    etime = etime + (temp_day-stad(3))*DAYSEC

    !plus month&year
    stayear = stad(1) ; stamonth = stad(2) ; endmonth = 12
    year_loop: do while (stayear <= endd(1))
      if (stayear == endd(1)) then
        endmonth = endd(2) - 1
      end if
      month_loop: do while (stamonth <= endmonth)
        temp_day = get_days(stayear, stamonth)
        etime = etime + temp_day*DAYSEC
        stamonth = stamonth + 1
      end do month_loop
      stamonth = 1 ; stayear = stayear + 1
    end do year_loop

  end subroutine calc_etime

  subroutine get_out_vari(outv_num, out_vari, vari_name)
  !***************************************************************************************
  ! get_out_vari -- Get output variables
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(inout) :: outv_num
    character(*), intent(in) :: out_vari
    character(*), intent(inout) :: vari_name(:)
    ! -- local
    integer(I4) :: i
    integer(I4) :: first_pos, sta_pos, end_pos
    !-------------------------------------------------------------------------------------
    vari_name(1:3) = ['conv', 'head', 'rest']
    first_pos = index(out_vari, ",")
    sta_pos = 1 ; end_pos = first_pos

    do i = 1, outv_num-1
      vari_name(i+3) = out_vari(sta_pos:end_pos-1) ; sta_pos = end_pos + 1
      end_pos = index(out_vari(sta_pos:), ",") + sta_pos - 1
    end do
    vari_name(outv_num+3) = out_vari(sta_pos:)

    outv_num = outv_num + 3

  end subroutine get_out_vari

  subroutine set_output(outv_num, vari_name)
  !***************************************************************************************
  ! set_output -- Set output
  !***************************************************************************************
    ! -- modules
    use initial_module, only: out_type, st_out_type, st_out_path, st_out_unit, st_out_time
    ! -- inout
    integer(I4), intent(in) :: outv_num
    character(*), intent(in) :: vari_name(:)
    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: head_time, rest_time, srat_time, wtab_time, mass_time
    integer(I4) :: velc_time, rivr_time, lakr_time, sufr_time, dunr_time
    integer(I4) :: seal_time, well_time, rech_time
    character(:), allocatable :: str_sim_type, str_sim_name
    character(TIMELEN) :: head_unit, rest_unit, srat_unit, wtab_unit, mass_unit
    character(TIMELEN) :: velc_unit, rivr_unit, lakr_unit, sufr_unit, dunr_unit
    character(TIMELEN) :: seal_unit, well_unit, rech_unit
    namelist/set_out_unit/head_unit, rest_unit, srat_unit, wtab_unit, mass_unit,&
                          velc_unit, rivr_unit, lakr_unit, sufr_unit, dunr_unit,&
                          seal_unit, well_unit, rech_unit
    namelist/set_out_time/head_time, rest_time, srat_time, wtab_time, mass_time,&
                          velc_time, rivr_time, lakr_time, sufr_time, dunr_time,&
                          seal_time, well_time, rech_time
    !-------------------------------------------------------------------------------------
    ierr = 0
    head_unit = "" ; rest_unit = "" ; srat_unit = "" ; wtab_unit = "" ; mass_unit = ""
    velc_unit = "" ; rivr_unit = "" ; lakr_unit = "" ; sufr_unit = "" ; dunr_unit = ""
    seal_unit = "" ; well_unit = "" ; rech_unit = ""
    read(unit=main_fnum,nml=set_out_unit,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading output unit section in main file.")
    end if

    read(unit=main_fnum,nml=set_out_time,iostat=ierr)
    if (ierr /= 0) then
      call write_err_stop("While reading output interval time section in main file.")
    end if

    if (st_sim%sim_type == -1) then
      str_sim_type = "stat"
    else if (st_sim%sim_type == 0) then
      str_sim_type = "long"
    else if (st_sim%sim_type == 1) then
      str_sim_type = "tran"
    end if

    str_sim_name = str_sim_type//"_"//st_sim%sim_name

    deallocate(str_sim_type)

    do i = 1, outv_num
      select case (vari_name(i))
      case('conv')
        st_out_path%conv = "conv_"//str_sim_name//".txt"
      case('head')
        st_out_path%head = "head_"//str_sim_name//".bin"
        st_out_unit%head = head_unit ; st_out_time%head = head_time
      case('rest')
        st_out_path%rest = "rest_"//str_sim_name//".bin"
        st_out_unit%rest = rest_unit ; st_out_time%rest = rest_time
      case('srat')
        st_out_type%srat = out_type(3)
        st_out_path%srat = "srat_"//str_sim_name//".bin"
        st_out_unit%srat = srat_unit ; st_out_time%srat = srat_time
      case('wtab')
        st_out_type%wtab = out_type(2)
        st_out_path%wtab = "wtab_"//str_sim_name//".bin"
        st_out_unit%wtab = wtab_unit ; st_out_time%wtab = wtab_time
      case('mass')
        st_out_type%mass = out_type(1)
        st_out_path%mass = "mass_"//str_sim_name//".csv"
        st_out_unit%mass = mass_unit ; st_out_time%mass = mass_time
      case('velc')
        st_out_type%velc = out_type(3)
        st_out_path%velx = "velx_"//str_sim_name//".bin"
        st_out_path%vely = "vely_"//str_sim_name//".bin"
        st_out_path%velz = "velz_"//str_sim_name//".bin"
        st_out_unit%velc = velc_unit ; st_out_time%velc = velc_time
      case('rivr')
        st_out_type%rivr = out_type(2)
        st_out_path%rivr = "rivr_"//str_sim_name//".bin"
        st_out_unit%rivr = rivr_unit ; st_out_time%rivr = rivr_time
      case('lakr')
        st_out_type%lakr = out_type(2)
        st_out_path%lakr = "lakr_"//str_sim_name//".bin"
        st_out_unit%lakr = lakr_unit ; st_out_time%lakr = lakr_time
      case('sufr')
        st_out_type%sufr = out_type(2)
        st_out_path%sufr = "sufr_"//str_sim_name//".bin"
        st_out_unit%sufr = sufr_unit ; st_out_time%sufr = sufr_time
      case('dunr')
        st_out_type%dunr = out_type(2)
        st_out_path%dunr = "dunr_"//str_sim_name//".bin"
        st_out_unit%dunr = dunr_unit ; st_out_time%dunr = dunr_time
      case('seal')
        st_out_type%seal = out_type(3)
        st_out_path%seal = "seal_"//str_sim_name//".bin"
        st_out_unit%seal = seal_unit ; st_out_time%seal = seal_time
      case('well')
        st_out_type%well = out_type(3)
        st_out_path%well = "well_"//str_sim_name//".bin"
        st_out_unit%well = well_unit ; st_out_time%well = well_time
      case('rech')
        st_out_type%rech = out_type(2)
        st_out_path%rech = "rech_"//str_sim_name//".bin"
        st_out_unit%rech = rech_unit ; st_out_time%rech = rech_time
      case('calg')
        st_out_type%calg = out_type(1)
        st_out_path%calg = "calg_"//str_sim_name//".csv"
      end select
    end do

    deallocate(str_sim_name)

  end subroutine set_output

end module read_input
