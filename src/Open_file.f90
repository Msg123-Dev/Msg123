module open_file
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SZERO, SINFI, HOURSEC, CHALEN, TIMELEN
  use utility_module, only: conv_unit, close_file, open_new_rtxt, open_new_rbin,&
                            open_new_wtxt, write_err_stop, write_success, write_logf
  use initial_module, only: st_sim, st_grid, st_step_flag, pro_totn, my_rank, in_type
  use read_module, only: skip_file, skip_file_int
#ifdef MPI_MSG
  use mpi_utility, only: bcast_val, bcast_file, bcast_extr_set
  use utility_module, only: log_fnum
  use mpi_read, only: open_mpi_read_file, set_real4_fview, set_int4_fview,&
                      skip_mpi_file, skip_mpi_file_int
#endif

  implicit none
  private
  public :: open_in_retnf, open_in_parmf, open_in_geogf, open_in_initf
  public :: open_in_sealf, open_in_rechf, open_in_wellf, open_in_wlayf
  public :: open_in_precf, open_in_evapf, open_in_rivef, open_in_lakef, open_in_massf
  public :: open_out_convf, open_out_binf, open_out_massf
  integer(I4), public :: wells_fnum, welle_fnum, inmas_fnum

  type, public :: st_intf
    integer(I4) :: type, fnum
    real(SP) :: step
  end type st_intf
  type(st_intf), public :: st_intse, st_intre, st_intwe, st_intpr, st_intev

  ! -- local

  contains

  subroutine open_in_retnf(retn_type, retn_path, view_calc)
  !***************************************************************************************
  ! open_in_retnf -- Open input retention file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_retf_type, st_retn_path, st_retn_fnum, st_retn
    ! -- inout
    integer(I4), intent(in) :: retn_type
    character(*), intent(in) :: retn_path
    integer(I4), intent(in), optional :: view_calc
    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: retn_num, temp_view
    integer(I4) :: retn_ftype(3)
    integer(I4) :: vana_type, vann_type, resi_type
    character(CHALEN) :: in_vapath, in_vnpath, in_repath
    character(:), allocatable :: mess_vana, mess_vann, mess_resi
    character(:), allocatable :: path_retn, mess_retn
    namelist/inretn_type/vana_type, vann_type, resi_type
    namelist/inretn_path/in_vapath, in_vnpath, in_repath
    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, retn_path, "input retention", st_retn%fnum)
      if (retn_type == in_type(0)) then
        in_vapath = "" ; in_vnpath = "" ; in_repath = ""
        read(unit=st_retn%fnum,nml=inretn_type,iostat=ierr)
        if (ierr /= 0) then
          call write_err_stop("While reading file type section in retention file.")
        end if
        read(unit=st_retn%fnum,nml=inretn_path,iostat=ierr)
        if (ierr /= 0) then
          call write_err_stop("While reading file path section in retention file.")
        end if
        mess_vana = "input van genuchten parameter alpha"
        mess_vann = "input van genuchten parameter n"
        mess_resi = "input residual water content"
      else if (retn_type  == in_type(1)) then
        read(unit=st_retn%fnum,fmt=*,iostat=ierr) st_retn%totn
        if (ierr /= 0) then
          call write_err_stop("While reading retention number in classification file.")
        else if (st_retn%totn <= 0) then
          call write_err_stop("Retention number is below 0.")
        end if
        read(unit=st_retn%fnum,fmt=*,iostat=ierr)
        if (ierr /= 0) then
          call write_err_stop("While reading retention number in classification file.")
        end if
      end if
    end if

    if (present(view_calc)) then
      temp_view = view_calc
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast scalar value (val)
        call bcast_val(st_retn%totn, "retention number")
      if (retn_type == in_type(0)) then
        ! -- Bcast file (file)
          call bcast_file(vana_type, in_vapath, mess_vana)
          call bcast_file(vann_type, in_vnpath, mess_vann)
          call bcast_file(resi_type, in_repath, mess_resi)
      end if
    end if
#endif

    if (retn_type == in_type(0)) then
      retn_ftype(1) = vana_type ; retn_ftype(2) = vann_type ; retn_ftype(3) = resi_type
      st_retf_type%vana = vana_type ; st_retf_type%vann = vann_type
      st_retf_type%resi = resi_type
      st_retn_path%vana = trim(adjustl(in_vapath))
      st_retn_path%vann = trim(adjustl(in_vnpath))
      st_retn_path%resi = trim(adjustl(in_repath))

      do i = 1, 3
        retn_num = 0
        select case (i)
        case(1)
          path_retn = trim(adjustl(in_vapath)) ; mess_retn = mess_vana
          deallocate(mess_vana)
        case(2)
          path_retn = trim(adjustl(in_vnpath)) ; mess_retn = mess_vann
          deallocate(mess_vann)
        case(3)
          path_retn = trim(adjustl(in_repath)) ; mess_retn = mess_resi
          deallocate(mess_resi)
        end select

        if (retn_ftype(i) == in_type(3) .or. retn_ftype(i) == in_type(5)) then
          if (my_rank == 0) then
            ! -- Open new read text file (new_rtxt)
              call open_new_rtxt(1, 1, path_retn, mess_retn, retn_num)
          end if

#ifdef MPI_MSG
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(retn_num, "retention file number")
          end if
#endif

        else if (retn_ftype(i) == in_type(4) .or. retn_ftype(i) == in_type(6)) then
#ifdef MPI_MSG
          ! -- Open mpi read file (mpi_read_file)
            call open_mpi_read_file(1, 1, path_retn, mess_retn, retn_num)
          ! -- Set real4 file view (real4_fview)
            call set_real4_fview(retn_num, temp_view, mess_retn)
#else
          ! -- Open new read binary file (new_rbin)
            call open_new_rbin(1, 1, path_retn, mess_retn, retn_num)
#endif
        end if

        select case (i)
        case(1)
          st_retn_fnum%vana = retn_num
        case(2)
          st_retn_fnum%vann = retn_num
        case(3)
          st_retn_fnum%resi = retn_num
        end select

      end do
      deallocate(path_retn, mess_retn)
    end if

  end subroutine open_in_retnf

  subroutine open_in_parmf(parm_type, parm_path, view_calc)
  !***************************************************************************************
  ! open_in_parmf -- Open input parameter file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_parf_type, st_parm_path, st_parm_fnum, st_parm
    ! -- inout
    integer(I4), intent(in) :: parm_type
    character(*), intent(in) :: parm_path
    integer(I4), intent(in), optional :: view_calc
    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: parm_num, temp_view
    integer(I4) :: parm_ftype(5)
    integer(I4) :: pakx_type, paky_type, pakz_type, pass_type, pats_type
    character(CHALEN) :: in_kxpath, in_kypath, in_kzpath, in_sspath, in_tspath
    character(:), allocatable :: mess_kx, mess_ky, mess_kz, mess_ss, mess_ts
    character(:), allocatable :: path_parm, mess_parm
    namelist/inparm_type/pakx_type, paky_type, pakz_type, pass_type, pats_type
    namelist/inparm_path/in_kxpath, in_kypath, in_kzpath, in_sspath, in_tspath
    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, parm_path, "input parameter", st_parm%fnum)
      if (parm_type == in_type(0)) then
        in_kxpath = "" ; in_kypath = "" ; in_kzpath = "" ; in_sspath = "" ; in_tspath = ""
        read(unit=st_parm%fnum,nml=inparm_type,iostat=ierr)
        if (ierr /= 0) then
          call write_err_stop("While reading file type section in parameter file.")
        end if
        read(unit=st_parm%fnum,nml=inparm_path,iostat=ierr)
        if (ierr /= 0) then
          call write_err_stop("While reading file path section in parameter file.")
        end if
        mess_kx = "input saturated hydraulic conductivity in x direction"
        mess_ky = "input saturated hydraulic conductivity in y direction"
        mess_kz = "input saturated hydraulic conductivity in z direction"
        mess_ss = "input specific storage"
        mess_ts = "input porosity"
      else if (parm_type == in_type(1)) then
        read(unit=st_parm%fnum,fmt=*,iostat=ierr) st_parm%totn
        if (ierr /= 0) then
          call write_err_stop("While reading parameter number in classification file.")
        else if (st_parm%totn <= 0) then
          call write_err_stop("Parameter number is below 0.")
        end if
        read(unit=st_parm%fnum,fmt=*,iostat=ierr)
        if (ierr /= 0) then
          call write_err_stop("While reading parameter number in classification file.")
        end if
      end if
    end if

    if (present(view_calc)) then
      temp_view = view_calc
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast scalar value (val)
        call bcast_val(st_parm%totn, "parameter number")
      if (parm_type == in_type(0)) then
        ! -- Bcast file (file)
          call bcast_file(pakx_type, in_kxpath, mess_kx)
          call bcast_file(paky_type, in_kypath, mess_ky)
          call bcast_file(pakz_type, in_kzpath, mess_kz)
          call bcast_file(pass_type, in_sspath, mess_ss)
          call bcast_file(pats_type, in_tspath, mess_ts)
      end if
    end if
#endif

    if (parm_type == in_type(0)) then
      parm_ftype(1) = pakx_type ; parm_ftype(2) = paky_type ; parm_ftype(3) = pakz_type
      parm_ftype(4) = pass_type ; parm_ftype(5) = pats_type
      st_parf_type%pakx = pakx_type ; st_parf_type%paky = paky_type
      st_parf_type%pakz = pakz_type ; st_parf_type%pass = pass_type
      st_parf_type%pats = pats_type
      st_parm_path%pakx = trim(adjustl(in_kxpath))
      st_parm_path%paky = trim(adjustl(in_kypath))
      st_parm_path%pakz = trim(adjustl(in_kzpath))
      st_parm_path%pass = trim(adjustl(in_sspath))
      st_parm_path%pats = trim(adjustl(in_tspath))

      do i = 1, 5
        parm_num = 0
        select case (i)
        case(1)
          path_parm = trim(adjustl(in_kxpath)) ; mess_parm = mess_kx
          deallocate(mess_kx)
        case(2)
          path_parm = trim(adjustl(in_kypath)) ; mess_parm = mess_ky
          deallocate(mess_ky)
        case(3)
          path_parm = trim(adjustl(in_kzpath)) ; mess_parm = mess_kz
          deallocate(mess_kz)
        case(4)
          path_parm = trim(adjustl(in_sspath)) ; mess_parm = mess_ss
          deallocate(mess_ss)
        case(5)
          path_parm = trim(adjustl(in_tspath)) ; mess_parm = mess_ts
          deallocate(mess_ts)
        end select

        if (parm_ftype(i) == in_type(3) .or. parm_ftype(i) == in_type(5)) then
          if (my_rank == 0) then
            ! -- Open new read text file (new_rtxt)
              call open_new_rtxt(1, 1, path_parm, mess_parm, parm_num)
          end if

#ifdef MPI_MSG
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(parm_num, "parameter file number")
          end if
#endif

        else if (parm_ftype(i) == in_type(4) .or. parm_ftype(i) == in_type(6)) then
#ifdef MPI_MSG
          ! -- Open mpi read file (mpi_read_file)
            call open_mpi_read_file(1, 1, path_parm, mess_parm, parm_num)
          ! -- Set real4 file view (real4_fview)
            call set_real4_fview(parm_num, temp_view, mess_parm)
#else
          ! -- Open new read binary file (new_rbin)
            call open_new_rbin(1, 1, path_parm, mess_parm, parm_num)
#endif
        end if

        select case (i)
        case(1)
          st_parm_fnum%pakx = parm_num
        case(2)
          st_parm_fnum%paky = parm_num
        case(3)
          st_parm_fnum%pakz = parm_num
        case(4)
          st_parm_fnum%pass = parm_num
        case(5)
          st_parm_fnum%pats = parm_num
        end select

      end do
      deallocate(path_parm, mess_parm)
    end if

  end subroutine open_in_parmf

  subroutine open_in_geogf(geog_path, geog_view)
  !***************************************************************************************
  ! open_in_geogf -- Open geography file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_geof_type, st_geog_path, st_geog_fnum
    ! -- inout
    character(*), intent(in) :: geog_path
    integer(I4), intent(in), optional :: geog_view
    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: geog_num, temp_view
    integer(I4) :: geog_ftype(3)
    integer(I4) :: geoz_type, geor_type, geoa_type
    character(CHALEN) :: in_gzpath, in_grpath, in_gapath
    character(:), allocatable :: mess_geoz, mess_geor, mess_geoa
    character(:), allocatable :: path_geog, mess_geog
    namelist/ingeog_type/geoz_type, geor_type, geoa_type
    namelist/ingeog_path/in_gzpath, in_grpath, in_gapath
    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, geog_path, "input geography", geog_num)

      read(unit=geog_num,nml=ingeog_type,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading file type section in geography file.")
      end if
      in_gzpath = "" ; in_grpath = "" ;  in_gapath = ""
      read(unit=geog_num,nml=ingeog_path,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading file path section in geography file.")
      end if

      call close_file(geog_num)
    end if

    if (present(geog_view)) then
      temp_view = geog_view
    end if

    mess_geoz = "input minium z elevation"
    mess_geor = "input relief"
    mess_geoa = "input geography parameter"

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast file (file)
        call bcast_file(geoz_type, in_gzpath, mess_geoz)
        call bcast_file(geor_type, in_grpath, mess_geor)
        call bcast_file(geoa_type, in_gapath, mess_geoa)
    end if
#endif
    geog_ftype(1) = geoz_type ; geog_ftype(2) = geor_type ; geog_ftype(3) = geoa_type
    st_geof_type%geoz = geoz_type ; st_geof_type%geor = geor_type
    st_geof_type%geoa = geoa_type
    st_geog_path%geoz = trim(adjustl(in_gzpath))
    st_geog_path%geor = trim(adjustl(in_grpath))
    st_geog_path%geoa = trim(adjustl(in_gapath))

    do i = 1, 3
      geog_num = 0
      select case (i)
      case(1)
        path_geog = trim(adjustl(in_gzpath)) ; mess_geog = mess_geoz
        deallocate(mess_geoz)
      case(2)
        path_geog = trim(adjustl(in_grpath)) ; mess_geog = mess_geor
        deallocate(mess_geor)
      case(3)
        path_geog = trim(adjustl(in_gapath)) ; mess_geog = mess_geoa
        deallocate(mess_geoa)
      end select

      if (geog_ftype(i) == in_type(3)) then
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, path_geog, mess_geog, geog_num)
        end if

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(geog_num, "geography file number")
        end if
#endif

      else if (geog_ftype(i) == in_type(4)) then
#ifdef MPI_MSG
        ! -- Open mpi read file (mpi_read_file)
          call open_mpi_read_file(1, 1, path_geog, mess_geog, geog_num)
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(geog_num, temp_view, mess_geog)
#else
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, path_geog, mess_geog, geog_num)
#endif
      end if

      select case (i)
      case(1)
        st_geog_fnum%geoz = geog_num
      case(2)
        st_geog_fnum%geor = geog_num
      case(3)
        st_geog_fnum%geoa = geog_num
      end select

    end do

    deallocate(path_geog, mess_geog)

  end subroutine open_in_geogf

  subroutine open_in_initf(init_type, init_path, init_unit, view_rest)
  !***************************************************************************************
  ! open_in_initf -- Open input initial file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_init
#ifdef MPI_MSG
    use mpi_read, only: set_real8_fview
#endif
    ! -- inout
    integer(I4), intent(in) :: init_type
    character(*), intent(inout) :: init_path, init_unit
    integer(I4), intent(in), optional :: view_rest
    ! -- local
    integer(I4) :: temp_view
    integer(I4), allocatable :: txt_init_type(:)
    character(:), allocatable :: err_mes
    logical, allocatable :: init_mask(:)
    !-------------------------------------------------------------------------------------
    allocate(txt_init_type(2), init_mask(2))
    txt_init_type(:) = [in_type(3), in_type(5)]
    init_mask(:) = (init_type == txt_init_type(:))

    err_mes = "input initial"
    if (any(init_mask)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, init_path, err_mes, st_init%fnum)
      end if
    else if (init_type /= in_type(7)) then
      if (present(view_rest)) then
        temp_view = view_rest
      end if
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, init_path, err_mes, st_init%fnum)
      ! -- Set real8 file view (real8_fview)
        call set_real8_fview(st_init%fnum, temp_view, err_mes)
#else
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, init_path, err_mes, st_init%fnum)
#endif
    end if

    if (len_trim(adjustl(init_unit)) == 0) then
      st_init%multi = st_sim%cal_fact
    else
      ! -- Convert unit (unit)
        call conv_unit(my_rank, init_unit, err_mes, st_sim%sta_date, st_init%multi)
    end if

    deallocate(txt_init_type)
    deallocate(init_mask)

  end subroutine open_in_initf

  subroutine open_in_sealf(seal_type, seal_path, seal_unit, seal_view)
  !***************************************************************************************
  ! open_in_sealf -- Open input sea level file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_seal
    ! -- inout
    integer(I4), intent(in) :: seal_type
    character(*), intent(inout) :: seal_path, seal_unit
    integer(I4), intent(in), optional :: seal_view
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: temp_view
    integer(I4) :: intse_fnum, intse_type
    integer(I4), allocatable :: txt_seal_type(:), bin_seal_type(:)
    real(SP) :: intse_step, intse_end
    character(CHALEN) :: intsep
    character(:), allocatable :: vname, err_mes
    logical, allocatable :: txt_seal_mask(:), bin_seal_mask(:)
    !-------------------------------------------------------------------------------------
    allocate(txt_seal_type(2), txt_seal_mask(2))
    allocate(bin_seal_type(2), bin_seal_mask(2))
    txt_seal_type(:) = [in_type(3), in_type(5)]
    bin_seal_type(:) = [in_type(4), in_type(6)]
    txt_seal_mask(:) = (seal_type == txt_seal_type(:))
    bin_seal_mask(:) = (seal_type == bin_seal_type(:))

    if (present(seal_view)) then
      temp_view = seal_view
    end if

    vname = "sea level"
    st_intse%fnum = 0 ; st_intse%type = 0 ; st_intse%step = SZERO

    if (any(txt_seal_mask) .or. seal_type == in_type(1) .or. seal_type == in_type(2)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, seal_path, "input "//vname, st_seal%fnum)
      end if
    else if (any(bin_seal_mask)) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, seal_path, "input "//vname, st_seal%fnum)
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(st_seal%fnum, temp_view, "input "//vname)
#else
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, seal_path, "input "//vname, st_seal%fnum)
#endif
    else if (seal_type == in_type(7)) then
      err_mes = vname//" time interval"
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, seal_path, "input "//vname, intse_fnum)
        read(unit=intse_fnum,fmt=*,iostat=ierr) intse_type, intse_step, intse_end
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        read(unit=intse_fnum,fmt='(a)',iostat=ierr) intsep
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        st_intse%fnum = intse_fnum
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast extra setting (extr_set)
          call bcast_extr_set(intse_type, intse_step, intse_end, intsep, err_mes)
      end if
#endif
      txt_seal_mask(:) = (intse_type == txt_seal_type(:))
      bin_seal_mask(:) = (intse_type == bin_seal_type(:))
      if (any(txt_seal_mask)) then
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, trim(adjustl(intsep)), "input "//err_mes, st_seal%fnum)
        end if
      else if (any(bin_seal_mask)) then
#ifdef MPI_MSG
        ! -- Open mpi read file (mpi_read_file)
          call open_mpi_read_file(1, 1, trim(adjustl(intsep)), "input "//err_mes, st_seal%fnum)
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_seal%fnum, temp_view, "input "//err_mes)
#else
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, trim(adjustl(intsep)), "input "//err_mes, st_seal%fnum)
#endif
      else
        if (my_rank == 0) then
          call write_err_stop("Specified wrong number in "//err_mes//" file.")
        end if
      end if
      st_intse%type = intse_type ; st_intse%step = intse_step
      deallocate(err_mes)
    end if

    if (len_trim(adjustl(seal_unit)) == 0) then
      if (my_rank == 0) then
        call write_logf("Not specified "//vname//" time unit. Simulation time unit is used")
      end if
      st_seal%multi = SINFI ; st_seal%etime = SINFI
    else
      ! -- Convert unit (unit)
        call conv_unit(my_rank, seal_unit, "input "//vname, st_sim%sta_date, st_seal%multi)
    end if

    st_seal%totn = st_grid%nxyz

    if (st_seal%multi /= SINFI) then
      if (any(txt_seal_mask) .or. seal_type == in_type(1) .or. seal_type == in_type(2)) then
        if (my_rank == 0) then
          call skip_file(seal_type, st_seal%fnum, vname, st_seal%multi, st_seal%totn, st_seal%etime)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(st_seal%totn, "sea level number")
        end if
#endif
      else if (any(bin_seal_mask)) then
#ifdef MPI_MSG
        call skip_mpi_file(seal_type, st_seal%fnum, temp_view, vname, st_seal%multi, st_seal%etime)
#else
        call skip_file(seal_type, st_seal%fnum, vname, st_seal%multi, st_seal%totn, st_seal%etime)
#endif
      else if (seal_type == in_type(7)) then
#ifdef MPI_MSG
        if (any(txt_seal_mask)) then
          if (my_rank == 0) then
            call skip_file_int(intse_type, intse_fnum, st_seal%fnum, vname, st_seal%multi, intse_step, intse_end, st_seal%etime)
          end if
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(st_seal%fnum, vname//" file number")
          end if
        else if (any(bin_seal_mask)) then
          call skip_mpi_file_int(intse_fnum, vname, st_seal%multi, intse_step, intse_end, st_seal%fnum, st_seal%etime)
        end if
#else
        call skip_file_int(intse_type, intse_fnum, st_seal%fnum, vname, st_seal%multi, intse_step, intse_end, st_seal%etime)
#endif
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast scalar value (val)
          call bcast_val(st_seal%etime, vname//" end time value")
      end if
#endif
    else
      if (seal_type == in_type(1) .or. seal_type == in_type(2)) then
        if (my_rank == 0) then
          read(unit=st_seal%fnum,fmt=*,iostat=ierr) st_seal%totn
          if (ierr /= 0) then
            call write_err_stop("While reading the total number in "//vname//" file.")
          end if
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(st_seal%totn, vname//" total number")
        end if
#endif
      end if
    end if

    if (st_seal%etime > st_sim%end_time) then
      st_seal%etime = st_sim%end_time
    end if

    st_step_flag%seal = 0

    deallocate(vname)
    deallocate(txt_seal_type, bin_seal_type)
    deallocate(txt_seal_mask, bin_seal_mask)

  end subroutine open_in_sealf

  subroutine open_in_rechf(rech_type, rech_path, rech_unit, rech_view)
  !***************************************************************************************
  ! open_in_rechf -- Open input recharge file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_rech
    ! -- inout
    integer(I4), intent(in) :: rech_type
    character(*), intent(inout) :: rech_path, rech_unit
    integer(I4), intent(in), optional :: rech_view
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: temp_view
    integer(I4) :: intre_fnum, intre_type
    real(SP) :: intre_step, intre_end
    character(CHALEN) :: intrep
    character(:), allocatable :: vname, err_mes
    !-------------------------------------------------------------------------------------
    if (present(rech_view)) then
      temp_view = rech_view
    end if

    vname = "recharge"
    st_intre%fnum = 0 ; st_intre%type = 0 ; st_intre%step = SZERO

    if (rech_type == in_type(1) .or. rech_type == in_type(3)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, rech_path, "input "//vname, st_rech%fnum)
      end if
    else if (rech_type == in_type(4)) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, rech_path, "input "//vname, st_rech%fnum)
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(st_rech%fnum, temp_view, "input "//vname)
#else
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, rech_path, "input "//vname, st_rech%fnum)
#endif
    else if (rech_type == in_type(7)) then
      intre_type = 0 ; err_mes = vname//" time interval"
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, rech_path, "input "//vname, intre_fnum)
        read(unit=intre_fnum,fmt=*,iostat=ierr) intre_type, intre_step, intre_end
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        read(unit=intre_fnum,fmt='(a)',iostat=ierr) intrep
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        st_intre%fnum = intre_fnum
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast extra setting (extr_set)
          call bcast_extr_set(intre_type, intre_step, intre_end, intrep, err_mes)
      end if
#endif
      if (intre_type == in_type(3)) then
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, trim(adjustl(intrep)), "input "//err_mes, st_rech%fnum)
        end if
      else if (intre_type == in_type(4)) then
#ifdef MPI_MSG
        ! -- Open mpi read file (mpi_read_file)
          call open_mpi_read_file(1, 1, trim(adjustl(intrep)), "input "//err_mes, st_rech%fnum)
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_rech%fnum, temp_view, "input "//err_mes)
#else
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, trim(adjustl(intrep)), "input "//err_mes, st_rech%fnum)
#endif
      else
        if (my_rank == 0) then
          call write_err_stop("Specified wrong number in "//err_mes//" file.")
        end if
      end if
      st_intre%type = intre_type ; st_intre%step = intre_step
      deallocate(err_mes)
    end if

    st_rech%uni_conv = 1.0E-3_SP/HOURSEC

    if (len_trim(adjustl(rech_unit)) == 0) then
      if (my_rank == 0) then
        call write_logf("Not specified "//vname//" time unit. Simulation time unit is used")
      end if
      st_rech%multi = SINFI ; st_rech%etime = SINFI
    else
      ! -- Convert unit (unit)
        call conv_unit(my_rank, rech_unit, "input "//vname, st_sim%sta_date, st_rech%multi)
    end if

    st_rech%totn = st_grid%nx*st_grid%ny

    if (st_rech%multi /= SINFI) then
      if (rech_type == in_type(1) .or. rech_type == in_type(3)) then
        if (my_rank == 0) then
          call skip_file(rech_type, st_rech%fnum, vname, st_rech%multi, st_rech%totn, st_rech%etime)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(st_rech%totn, vname//" total number")
        end if
#endif
      else if (rech_type == in_type(4)) then
#ifdef MPI_MSG
        call skip_mpi_file(rech_type, st_rech%fnum, temp_view, vname, st_rech%multi, st_rech%etime)
#else
        call skip_file(rech_type, st_rech%fnum, vname, st_rech%multi, st_rech%totn, st_rech%etime)
#endif

      else if (rech_type == in_type(7)) then
#ifdef MPI_MSG
        if (intre_type == in_type(3)) then
          if (my_rank == 0) then
            call skip_file_int(intre_type, intre_fnum, st_rech%fnum, vname, st_rech%multi, intre_step, intre_end, st_rech%etime)
          end if
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(st_rech%fnum, vname//" file number")
          end if
        else if (intre_type == in_type(4)) then
          call skip_mpi_file_int(intre_fnum, vname, st_rech%multi, intre_step, intre_end, st_rech%fnum, st_rech%etime)
        end if
#else
        call skip_file_int(intre_type, intre_fnum, st_rech%fnum, vname, st_rech%multi, intre_step, intre_end, st_rech%etime)
#endif
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast scalar value (val)
          call bcast_val(st_rech%etime, vname//" end time value")
      end if
#endif
    else
      if (rech_type == in_type(1)) then
        if (my_rank == 0) then
          read(unit=st_rech%fnum,fmt=*,iostat=ierr) st_rech%totn
          if (ierr /= 0) then
            call write_err_stop("While reading the total number in "//vname//" file.")
          end if
        end if
#ifdef MPI_MSG
        ! -- Bcast scalar value (val)
          call bcast_val(st_rech%totn, vname//" total number")
#endif
      end if
    end if

    if (st_rech%etime > st_sim%end_time) then
      st_rech%etime = st_sim%end_time
    end if

    st_step_flag%rech = 0

    deallocate(vname)

  end subroutine open_in_rechf

  subroutine open_in_wellf(well_type, well_path, well_unit, well_view)
  !***************************************************************************************
  ! open_in_wellf -- Open input well file
  !***************************************************************************************
    ! -- modules
    use constval_module, only: SONE, DAYSEC
    use initial_module, only: st_well
    ! -- inout
    integer(I4), intent(in) :: well_type
    character(*), intent(inout) :: well_path, well_unit
    integer(I4), intent(in), optional :: well_view
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: temp_view
    integer(I4) :: intwe_fnum, intwe_type
    integer(I4), allocatable :: txt_well_type(:), bin_well_type(:)
    real(SP) :: intwe_step, intwe_end
    character(CHALEN) :: intwep
    character(:), allocatable :: vname, err_mes
    logical, allocatable :: txt_well_mask(:), bin_well_mask(:)
    !-------------------------------------------------------------------------------------
    allocate(txt_well_type(2), txt_well_mask(2))
    allocate(bin_well_type(2), bin_well_mask(2))
    txt_well_type(:) = [in_type(3), in_type(5)]
    bin_well_type(:) = [in_type(4), in_type(6)]
    txt_well_mask(:) = (well_type == txt_well_type(:))
    bin_well_mask(:) = (well_type == bin_well_type(:))

    if (present(well_view)) then
      temp_view = well_view
    end if

    vname = "well"
    st_intwe%fnum = 0 ; st_intwe%type = 0 ; st_intwe%step = SZERO

    if (any(txt_well_mask) .or. well_type == in_type(2)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, well_path, "input "//vname, st_well%fnum)
      end if
    else if (any(bin_well_mask)) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, well_path, "input "//vname, st_well%fnum)
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(st_well%fnum, temp_view, "input "//vname)
#else
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, well_path, "input "//vname, st_well%fnum)
#endif
    else if (well_type == in_type(7)) then
      err_mes = vname//" time interval"
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, well_path, "input "//vname, intwe_fnum)
        read(unit=intwe_fnum,fmt=*,iostat=ierr) intwe_type, intwe_step, intwe_end
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        read(unit=intwe_fnum,fmt='(a)',iostat=ierr) intwep
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        st_intwe%fnum = intwe_fnum
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast extra setting (extr_set)
          call bcast_extr_set(intwe_type, intwe_step, intwe_end, intwep, err_mes)
      end if
#endif
      txt_well_mask(:) = (txt_well_type(:) /= intwe_type)
      bin_well_mask(:) = (bin_well_type(:) /= intwe_type)
      if (any(txt_well_mask)) then
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, trim(adjustl(intwep)), "input "//err_mes, st_well%fnum)
        end if
      else if (any(bin_well_mask)) then
#ifdef MPI_MSG
        ! -- Open mpi read file (mpi_read_file)
          call open_mpi_read_file(1, 1, trim(adjustl(intwep)), "input "//err_mes, st_well%fnum)
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_well%fnum, temp_view, "input "//err_mes)
#else
        ! -- Open new read binary file (new_btxt)
          call open_new_rbin(1, 1, trim(adjustl(intwep)), "input "//err_mes, st_well%fnum)
#endif
      else
        if (my_rank == 0) then
          call write_err_stop("Specified wrong number in "//err_mes//" file.")
        end if
      end if
      st_intwe%type = intwe_type ; st_intwe%step = intwe_step
      deallocate(err_mes)
    end if

    st_well%uni_conv = SONE/DAYSEC

    if (len_trim(adjustl(well_unit)) == 0) then
      if (my_rank == 0) then
        call write_logf("Not specified "//vname//" time unit. Simulation time unit is used")
      end if
      st_well%multi = SINFI ; st_well%etime = SINFI
    else
      ! -- Convert unit (unit)
        call conv_unit(my_rank, well_unit, "input "//vname, st_sim%sta_date, st_well%multi)
    end if

    st_well%totn = st_grid%nxyz

    if (st_well%multi /= SINFI) then
      if (any(txt_well_mask) .or. well_type == in_type(2)) then
        if (my_rank == 0) then
          call skip_file(well_type, st_well%fnum, vname, st_well%multi, st_well%totn, st_well%etime)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(st_well%totn, vname//" total number")
        end if
#endif
      else if (any(bin_well_mask)) then
#ifdef MPI_MSG
        call skip_mpi_file(well_type, st_well%fnum, temp_view, vname, st_well%multi, st_well%etime)
#else
        call skip_file(well_type, st_well%fnum, vname, st_well%multi, st_well%totn, st_well%etime)
#endif
      else if (well_type == in_type(7)) then
#ifdef MPI_MSG
        if (any(txt_well_mask)) then
          if (my_rank == 0) then
            call skip_file_int(intwe_type, intwe_fnum, st_well%fnum, vname, st_well%multi, intwe_step, intwe_end, st_well%etime)
          end if
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(st_well%fnum, vname//" file number")
          end if
        else if (any(bin_well_mask)) then
          call skip_mpi_file_int(intwe_fnum, vname, st_well%multi, intwe_step, intwe_end, st_well%fnum, st_well%etime)
        end if
#else
        call skip_file_int(intwe_type, intwe_fnum, st_well%fnum, vname, st_well%multi, intwe_step, intwe_end, st_well%etime)
#endif
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast scalar value (val)
          call bcast_val(st_well%etime, vname//" end time value")
      end if
#endif
    else
      if (well_type == in_type(2)) then
        if (my_rank == 0) then
          read(unit=st_well%fnum,fmt=*,iostat=ierr) st_well%totn
          if (ierr /= 0) then
            call write_err_stop("While reading the total number in "//vname//" file.")
          end if
        end if
#ifdef MPI_MSG
        ! -- Bcast scalar value (val)
          call bcast_val(st_well%totn, vname//" total number")
#endif
      end if
    end if

    if (st_well%etime > st_sim%end_time) then
      st_well%etime = st_sim%end_time
    end if

    st_step_flag%well = 0

    deallocate(vname)
    deallocate(txt_well_type, bin_well_type)
    deallocate(txt_well_mask, bin_well_mask)

  end subroutine open_in_wellf

  subroutine open_in_wlayf(weks_type, weke_type, wells_path, welle_path, wlay_view)
  !***************************************************************************************
  ! open_in_wlayf -- Open input well layer file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
!    use mpi_utility, only: bcast_file_path
    use mpi_utility, only: bcast_file
#endif
    ! -- inout
    integer(I4), intent(in) :: weks_type, weke_type
    character(*), intent(inout) :: wells_path, welle_path
    integer(I4), intent(in), optional :: wlay_view
    ! -- local
    integer(I4) :: temp_view
    character(:), allocatable :: vname
    !-------------------------------------------------------------------------------------
    if (present(wlay_view)) then
      temp_view = wlay_view
    end if

    vname = "well 2d"

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast file (file)
        call bcast_file(wells_path, vname//" start layer")
        call bcast_file(welle_path, vname//" end layer")
    end if
#endif

    if (weks_type == in_type(3)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, wells_path, "input "//vname//" start layer", wells_fnum)
      end if
    end if

    if (weke_type == in_type(3)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, welle_path, "input "//vname//" end layer", welle_fnum)
      end if
    end if

    if (weks_type == in_type(4)) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, wells_path, "input "//vname//" start layer", wells_fnum)
      ! -- Set int4 file view (int4_fview)
        call set_int4_fview(wells_fnum, temp_view, "input "//vname//" start layer")
#else
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, wells_path, "input "//vname//" start layer", wells_fnum)
#endif
    end if

    if (weke_type == in_type(4)) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, welle_path, "input "//vname//" end layer", welle_fnum)
      ! -- Set int4 file view (int4_fview)
        call set_int4_fview(welle_fnum, temp_view, "input "//vname//" end layer")
#else
      if (present(wlay_view)) then
        temp_view = wlay_view
      end if
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, welle_path, "input "//vname//" end layer", welle_fnum)
#endif
    end if

    deallocate(vname)

  end subroutine open_in_wlayf

  subroutine open_in_precf(prec_type, prec_path, prec_unit, prec_view)
  !***************************************************************************************
  ! open_in_precf -- Open input precipitation file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_prec
    ! -- inout
    integer(I4), intent(in) :: prec_type
    character(*), intent(inout) :: prec_path, prec_unit
    integer(I4), intent(in), optional :: prec_view
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: temp_view
    integer(I4) :: intpr_fnum, intpr_type
    real(SP) :: intpr_step, intpr_end
    character(CHALEN) :: intprp
    character(:), allocatable :: vname, err_mes
    !-------------------------------------------------------------------------------------
    if (present(prec_view)) then
      temp_view = prec_view
    end if

    vname = "precipitation"
    st_intpr%fnum = 0 ; st_intpr%type = 0 ; st_intpr%step = SZERO

    if (prec_type == in_type(1) .or. prec_type == in_type(3)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, prec_path, "input "//vname, st_prec%fnum)
      end if
    else if (prec_type == in_type(4)) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, prec_path, "input "//vname, st_prec%fnum)
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(st_prec%fnum, temp_view, "input "//vname)
#else
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, prec_path, "input "//vname, st_prec%fnum)
#endif
    else if (prec_type == in_type(7)) then
      intpr_type = 0 ; err_mes = vname//" time interval"
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, prec_path, "input "//vname, intpr_fnum)
        read(unit=intpr_fnum,fmt=*,iostat=ierr) intpr_type, intpr_step, intpr_end
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        read(unit=intpr_fnum,fmt='(a)',iostat=ierr) intprp
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        st_intpr%fnum = intpr_fnum
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast extra setting (extr_set)
          call bcast_extr_set(intpr_type, intpr_step, intpr_end, intprp, err_mes)
      end if
#endif
      if (intpr_type == in_type(3)) then
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, trim(adjustl(intprp)), "input "//vname, st_prec%fnum)
        end if
      else if (intpr_type == in_type(4)) then
#ifdef MPI_MSG
        ! -- Open mpi read file (mpi_read_file)
          call open_mpi_read_file(1, 1, trim(adjustl(intprp)), "input "//err_mes, st_prec%fnum)
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_prec%fnum, temp_view, "input "//err_mes)
#else
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, trim(adjustl(intprp)), "input "//err_mes, st_prec%fnum)
#endif
      else
        if (my_rank == 0) then
          call write_err_stop("Specified wrong number in "//err_mes//" file.")
        end if
      end if
      st_intpr%type = intpr_type ; st_intpr%step = intpr_step
      deallocate(err_mes)
    end if

    st_prec%uni_conv = 1.0E-3_SP/HOURSEC

    if (len_trim(adjustl(prec_unit)) == 0) then
      if (my_rank == 0) then
        call write_logf("Not specified "//vname//" time unit. Simulation time unit is used")
      end if
      st_prec%multi = SINFI ; st_prec%etime = SINFI
    else
      ! -- Convert unit (unit)
        call conv_unit(my_rank, prec_unit, "input "//vname, st_sim%sta_date, st_prec%multi)
    end if

    st_prec%totn = st_grid%nx*st_grid%ny

    if (st_prec%multi /= SINFI) then
      if (prec_type == in_type(1) .or. prec_type == in_type(3)) then
        if (my_rank == 0) then
          call skip_file(prec_type, st_prec%fnum, vname, st_prec%multi, st_prec%totn, st_prec%etime)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(st_prec%totn, vname//" total number")
        end if
#endif
      else if (prec_type == in_type(4)) then
#ifdef MPI_MSG
        call skip_mpi_file(prec_type, st_prec%fnum, temp_view, vname, st_prec%multi, st_prec%etime)
#else
        call skip_file(prec_type, st_prec%fnum, vname, st_prec%multi, st_prec%totn, st_prec%etime)
#endif
      else if (prec_type == in_type(7)) then
#ifdef MPI_MSG
        if (intpr_type == in_type(3)) then
          if (my_rank == 0) then
            call skip_file_int(intpr_type, intpr_fnum, st_prec%fnum, vname, st_prec%multi, intpr_step, intpr_end, st_prec%etime)
          end if
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(st_prec%fnum, vname//" file number")
          end if
        else if (intpr_type == in_type(4)) then
          call skip_mpi_file_int(intpr_fnum, vname, st_prec%multi, intpr_step, intpr_end, st_prec%fnum, st_prec%etime)
        end if
#else
        call skip_file_int(intpr_type, intpr_fnum, st_prec%fnum, vname, st_prec%multi, intpr_step, intpr_end, st_prec%etime)
#endif
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast scalar value (val)
          call bcast_val(st_prec%etime, vname//" end time value")
      end if
#endif
    else
      if (prec_type == in_type(1)) then
        if (my_rank == 0) then
          read(unit=st_prec%fnum,fmt=*,iostat=ierr) st_prec%totn
          if (ierr /= 0) then
            call write_err_stop("While reading the total number in "//vname//" file.")
          end if
        end if
#ifdef MPI_MSG
        ! -- Bcast scalar value (val)
          call bcast_val(st_prec%totn, vname//" total number")
#endif
      end if
    end if

    if (st_prec%etime > st_sim%end_time) then
      st_prec%etime = st_sim%end_time
    end if

    st_step_flag%prec = 0

    deallocate(vname)

  end subroutine open_in_precf

  subroutine open_in_evapf(evap_type, evap_path, evap_unit, evap_view)
  !***************************************************************************************
  ! open_in_evapf -- Open input evapotranspiration file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_evap
    ! -- inout
    integer(I4), intent(in) :: evap_type
    character(*), intent(inout) :: evap_path, evap_unit
    integer(I4), intent(in), optional :: evap_view
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: temp_view
    integer(I4) :: intev_fnum, intev_type
    real(SP) :: intev_step, intev_end
    character(CHALEN) :: intevp
    character(:), allocatable :: vname, err_mes
    !-------------------------------------------------------------------------------------
    if (present(evap_view)) then
      temp_view = evap_view
    end if

    vname = "evapotranspiration"
    st_intev%fnum = 0 ; st_intev%type = 0 ; st_intev%step = SZERO

    if (evap_type == in_type(1) .or. evap_type == in_type(3)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, evap_path, "input "//vname, st_evap%fnum)
      end if
    else if (evap_type == in_type(4)) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, evap_path, "input "//vname, st_evap%fnum)
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(st_evap%fnum, temp_view, "input "//vname)
#else
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, evap_path, "input "//vname, st_evap%fnum)
#endif
    else if (evap_type == in_type(7)) then
      intev_type = 0 ; err_mes = vname//" time interval"
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, evap_path, "input "//vname, intev_fnum)
        read(unit=intev_fnum,fmt=*,iostat=ierr) intev_type, intev_step, intev_end
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        read(unit=intev_fnum,fmt='(a)',iostat=ierr) intevp
        if (ierr /= 0) then
          call write_err_stop("While reading list file in "//vname//" file.")
        end if
        st_intev%fnum = intev_fnum
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast extra setting (extr_set)
          call bcast_extr_set(intev_type, intev_step, intev_end, intevp, err_mes)
      end if
#endif
      if (intev_type == in_type(3)) then
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, trim(adjustl(intevp)), "input "//err_mes, st_evap%fnum)
        end if
      else if (intev_type == in_type(4)) then
#ifdef MPI_MSG
        ! -- Open mpi read file (mpi_read_file)
          call open_mpi_read_file(1, 1, trim(adjustl(intevp)), "input "//err_mes, st_evap%fnum)
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_evap%fnum, temp_view, "input "//err_mes)
#else
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, trim(adjustl(intevp)), "input "//err_mes, st_evap%fnum)
#endif
      else
        if (my_rank == 0) then
          call write_err_stop("Specified wrong number in "//err_mes//" file.")
        end if
      end if
      st_intev%type = intev_type ; st_intev%step = intev_step
      deallocate(err_mes)
    end if

    st_evap%uni_conv = 1.0E-3_SP/HOURSEC

    if (len_trim(adjustl(evap_unit)) == 0) then
      if (my_rank == 0) then
        call write_logf("Not specified "//vname//" time unit. Simulation time unit is used")
      end if
      st_evap%multi = SINFI ; st_evap%etime = SINFI
    else
      ! -- Convert unit (unit)
        call conv_unit(my_rank, evap_unit, "input "//vname, st_sim%sta_date, st_evap%multi)
    end if

    st_evap%totn = st_grid%nx*st_grid%ny

    if (st_evap%multi /= SINFI) then
      if (evap_type == in_type(1) .or. evap_type == in_type(3)) then
        if (my_rank == 0) then
          call skip_file(evap_type, st_evap%fnum, vname, st_evap%multi, st_evap%totn, st_evap%etime)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(st_evap%totn, vname//" total number")
        end if
#endif
      else if (evap_type == in_type(4)) then
#ifdef MPI_MSG
        call skip_mpi_file(evap_type, st_evap%fnum, temp_view, vname, st_evap%multi, st_evap%etime)
#else
        call skip_file(evap_type, st_evap%fnum, vname, st_evap%multi, st_evap%totn, st_evap%etime)
#endif
      else if (evap_type == in_type(7)) then
#ifdef MPI_MSG
        if (intev_type == in_type(3)) then
          if (my_rank == 0) then
            call skip_file_int(intev_type, intev_fnum, st_evap%fnum, vname, st_evap%multi, intev_step, intev_end, st_evap%etime)
          end if
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(st_evap%fnum, vname//" file number")
          end if
        else if (intev_type == in_type(4)) then
          call skip_mpi_file_int(intev_fnum, vname, st_evap%multi, intev_step, intev_end, st_evap%fnum, st_evap%etime)
        end if
#else
        call skip_file_int(intev_type, intev_fnum, st_evap%fnum, vname, st_evap%multi, intev_step, intev_end, st_evap%etime)
#endif
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast scalar value (val)
          call bcast_val(st_evap%etime, vname//" end time value")
      end if
#endif
    else
      if (evap_type == in_type(1)) then
        if (my_rank == 0) then
          read(unit=st_evap%fnum,fmt=*,iostat=ierr) st_evap%totn
          if (ierr /= 0) then
            call write_err_stop("While reading the total number in "//vname//" file.")
          end if
        end if
#ifdef MPI_MSG
        ! -- Bcast scalar value (val)
          call bcast_val(st_evap%totn, vname//" total number")
#endif
      end if
    end if

    if (st_evap%etime > st_sim%end_time) then
      st_evap%etime = st_sim%end_time
    end if

    st_step_flag%evap = 0

    deallocate(vname)

  end subroutine open_in_evapf

  subroutine open_in_rivef(rive_path, nohv, hv, riwlv, riwdv, riblv, ridev, riwiv, rilev)
  !***************************************************************************************
  ! open_in_rivef -- Read input river file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_rivf_type, st_riwl, st_riwd, st_ribl, st_ride, st_riwi,&
                              st_rile
    ! -- inout
    character(*), intent(in) :: rive_path
    integer(I4), intent(in) :: nohv, hv
    integer(I4), intent(out) :: riwlv, riwdv, riblv, ridev, riwiv, rilev
    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: rive_fnum, rive_totn, intri_type, intri_num, rive_stepflag
    integer(I4) :: riwl_type, riwd_type, ribl_type, ride_type, riwi_type, rile_type
    integer(I4) :: riven(6), rivet(6)
    integer(I4) :: rive_view
    integer(I4), allocatable :: riv_txt_type(:)
    real(SP) :: intri_end, rive_multi, rive_etime, intri_step
    character(CHALEN) :: riwl_path, riwd_path, ribl_path, ride_path
    character(CHALEN) :: riwi_path, rile_path, intri_path
    character(TIMELEN) :: riwl_unit, riwd_unit, ribl_unit, ride_unit
    character(TIMELEN) :: riwi_unit, rile_unit
    character(:), allocatable :: mess_wl, mess_wd, mess_bl, mess_de, mess_wi, mess_le
    character(:), allocatable :: path_rive, mess_rive, err_mes, unit_rive
    logical, allocatable :: rive_mask(:)
    namelist/inrive_type/riwl_type, riwd_type, ribl_type, ride_type, riwi_type, rile_type
    namelist/inrive_path/riwl_path, riwd_path, ribl_path, ride_path, riwi_path, rile_path
    namelist/inrive_unit/riwl_unit, riwd_unit, ribl_unit, ride_unit, riwi_unit, rile_unit
    !-------------------------------------------------------------------------------------
    riwl_type = 0 ; riwd_type = 0 ; ribl_type = 0 ; ride_type = 0 ; riwi_type = 0
    rile_type = 0
    riwl_path = "" ; riwd_path = "" ; ribl_path = "" ; ride_path = "" ; riwi_path = ""
    rile_path = ""
    riwl_unit = "" ; riwd_unit = "" ; ribl_unit = "" ; ride_unit = "" ; riwi_unit = ""
    rile_unit = ""
    riwlv = 0 ; riwdv = 0 ; riblv = 0 ; ridev = 0 ; riwiv = 0 ; rilev = 0

    if (my_rank == 0) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, rive_path, "input river", rive_fnum)
      read(unit=rive_fnum,nml=inrive_type,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading file type section in river file.")
      end if
      read(unit=rive_fnum,nml=inrive_path,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading file path section in river file.")
      end if
      read(unit=rive_fnum,nml=inrive_unit,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading file unit section in river file.")
      end if
      call close_file(rive_fnum)
    end if

    mess_wl = "input river water level"
    mess_wd = "input river water depth"
    mess_bl = "input river bottom level"
    mess_de = "input river depth"
    mess_wi = "input river width"
    mess_le = "input river length"

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast file (file)
        call bcast_file(riwl_type, riwl_path, riwl_unit, mess_wl)
        call bcast_file(riwd_type, riwd_path, riwd_unit, mess_wd)
        call bcast_file(ribl_type, ribl_path, ribl_unit, mess_bl)
        call bcast_file(ride_type, ride_path, ride_unit, mess_de)
        call bcast_file(riwi_type, riwi_path, riwi_unit, mess_wi)
        call bcast_file(rile_type, rile_path, rile_unit, mess_le)
    end if
#endif

    rivet(1) = riwl_type ; rivet(2) = riwd_type ; rivet(3) = ribl_type
    rivet(4) = ride_type ; rivet(5) = riwi_type ; rivet(6) = rile_type
    riven(1) = 0 ; riven(2) = 0 ; riven(3) = 0
    riven(4) = 0 ; riven(5) = 0 ; riven(6) = 0

    allocate(riv_txt_type(3), rive_mask(3))
    riv_txt_type(:) = [in_type(1:3)]

    do i = 1, 6
      rive_totn = 0
      rive_mask(:) = (rivet(i) /= riv_txt_type(:))
      if (all(rive_mask) .and. rivet(i) /= in_type(4) .and. rivet(i) /= in_type(7)) then
        cycle
      end if

      rive_view = 0 ; intri_type = 0 ; err_mes = "" ; path_rive = "" ; mess_rive = ""

      select case (i)
      case(1)
        path_rive = trim(adjustl(riwl_path)) ; mess_rive = mess_wl ; unit_rive = riwl_unit
        deallocate(mess_wl)
      case(2)
        path_rive = trim(adjustl(riwd_path)) ; mess_rive = mess_wd ; unit_rive = riwd_unit
        deallocate(mess_wd)
      case(3)
        path_rive = trim(adjustl(ribl_path)) ; mess_rive = mess_bl ; unit_rive = ribl_unit
        deallocate(mess_bl)
      case(4)
        path_rive = trim(adjustl(ride_path)) ; mess_rive = mess_de ; unit_rive = ride_unit
        deallocate(mess_de)
      case(5)
        path_rive = trim(adjustl(riwi_path)) ; mess_rive = mess_wi ; unit_rive = riwi_unit
        deallocate(mess_wi)
      case(6)
        path_rive = trim(adjustl(rile_path)) ; mess_rive = mess_le ; unit_rive = rile_unit
        deallocate(mess_le)
      case default

      end select

      if (any(rivet(i) == riv_txt_type(:))) then
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, path_rive, mess_rive, riven(i))
        end if
      else if (rivet(i) == in_type(4)) then
#ifdef MPI_MSG
        ! -- Open mpi read file (mpi_read_file)
          call open_mpi_read_file(1, 1, path_rive, mess_rive, riven(i))
        if (len_trim(adjustl(unit_rive)) == 0) then
          rive_view = nohv ; intri_type = 0
        else
          rive_view = hv ; intri_type = rivet(i)
        end if
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(riven(i), rive_view, mess_rive)
#else
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, path_rive, mess_rive, riven(i))
        rive_view = nohv ; rive_view = hv
#endif
      else if (rivet(i) == in_type(7)) then
        err_mes = mess_rive//" time interval"
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, path_rive, mess_rive, riven(i))
          read(unit=riven(i),fmt=*,iostat=ierr) intri_type, intri_step, intri_end
          if (ierr /= 0) then
            call write_err_stop("While reading list file in "//mess_rive//" file.")
          end if
          read(unit=riven(i),fmt='(a)',iostat=ierr) intri_path
          if (ierr /= 0) then
            call write_err_stop("While reading list file in "//mess_rive//" file.")
          end if
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast extra setting (extr_set)
            call bcast_extr_set(intri_type, intri_step, intri_end, intri_path, err_mes)
        end if
#endif
        if (intri_type == in_type(3)) then
          if (my_rank == 0) then
            ! -- Open new read text file (new_rtxt)
              call open_new_rtxt(1, 1, trim(adjustl(intri_path)), err_mes, intri_num)
          end if
        else if (intri_type == in_type(4)) then
#ifdef MPI_MSG
          ! -- Open mpi read file (mpi_read_file)
            call open_mpi_read_file(1, 1, trim(adjustl(intri_path)), err_mes, intri_num)
          rive_view = nohv
          ! -- Set real4 file view (real4_fview)
            call set_real4_fview(intri_num, rive_view, err_mes)
#else
          ! -- Open new read binary file (new_rbin)
            call open_new_rbin(1, 1, trim(adjustl(intri_path)), err_mes, intri_num)
          rive_view = nohv
#endif
        else
          if (my_rank == 0) then
            call write_err_stop("Specified wrong number in "//err_mes//" file.")
          end if
        end if
      end if

      if (len_trim(adjustl(unit_rive)) == 0) then
        if (my_rank == 0) then
          err_mes = "Not specified time unit in "//mess_rive//". Same value is used in simulation"
          call write_logf(err_mes)
        end if
        rive_multi = SINFI ; rive_etime = SINFI
      else
        ! -- Convert unit (unit)
          call conv_unit(my_rank, unit_rive, mess_rive, st_sim%sta_date, rive_multi)
      end if

      rive_totn = st_grid%nx*st_grid%ny

      if (rive_multi /= SINFI) then
        if (any(rivet(i) == riv_txt_type(:))) then
          if (my_rank == 0) then
            call skip_file(rivet(i), riven(i), mess_rive, rive_multi, rive_totn, rive_etime)
          end if
#ifdef MPI_MSG
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(rive_totn, mess_rive//" total number")
          end if
#endif
        else if (rivet(i) == in_type(4)) then
#ifdef MPI_MSG
          call skip_mpi_file(rivet(i), riven(i), rive_view, mess_rive, rive_multi, rive_etime)
#else
          call skip_file(rivet(i), riven(i), mess_rive, rive_multi, rive_totn, rive_etime)
#endif
        else if (rivet(i) == in_type(7)) then
#ifdef MPI_MSG
          if (intri_type == in_type(3)) then
            if (my_rank == 0) then
              call skip_file_int(intri_type, riven(i), intri_num, mess_rive, rive_multi, intri_step, intri_end, rive_etime)
            end if
            if (pro_totn /= 1) then
              ! -- Bcast scalar value (val)
                call bcast_val(intri_num, mess_rive//" file number")
            end if
          else if (intri_type == in_type(4)) then
            call skip_mpi_file_int(riven(i), mess_rive, rive_multi, intri_step, intri_end, intri_num, rive_etime)
          end if
#else
          call skip_file_int(intri_type, riven(i), intri_num, mess_rive, rive_multi, intri_step, intri_end, rive_etime)
#endif
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(rive_etime, mess_rive//" end time value")
        end if
#endif
      else
        if (rivet(i) == in_type(1) .or. rivet(i) == in_type(2)) then
          if (my_rank == 0) then
            read(unit=riven(i),fmt=*,iostat=ierr) rive_totn
            if (ierr /= 0) then
              err_mes = "While reading the total number in "//mess_rive//" file."
              call write_err_stop(err_mes)
            end if
          end if
#ifdef MPI_MSG
          ! -- Bcast scalar value (val)
            call bcast_val(rive_totn, mess_rive//" total number")
#endif
        end if
      end if

      if (rive_etime > st_sim%end_time) then
        rive_etime = st_sim%end_time
      end if

      rive_stepflag = 0

      select case (i)
      case(1)
        st_riwl%multi = rive_multi ; st_riwl%totn = rive_totn
        st_riwl%etime = rive_etime ; st_step_flag%riwl = rive_stepflag
        st_riwl%inttype = intri_type
#ifdef MPI_MSG
        riwlv = rive_view
#endif
        if (rivet(i) == in_type(7)) then
          st_riwl%fnum = intri_num ; st_riwl%intpath = trim(adjustl(intri_path))
          st_riwl%intstep = intri_step ; st_riwl%intfnum = riven(i)
        else
          st_riwl%fnum = riven(i)
        end if
      case(2)
        st_riwd%multi = rive_multi ; st_riwd%totn = rive_totn
        st_riwd%etime = rive_etime ; st_step_flag%riwd = rive_stepflag
        st_riwd%inttype = intri_type
#ifdef MPI_MSG
        riwdv = rive_view
#endif
        if (rivet(i) == in_type(7)) then
          st_riwd%fnum = intri_num ; st_riwd%intpath = trim(adjustl(intri_path))
          st_riwd%intstep = intri_step ; st_riwd%intfnum = riven(i)
        else
          st_riwd%fnum = riven(i)
        end if
      case(3)
        st_ribl%multi = rive_multi ; st_ribl%totn = rive_totn
        st_ribl%etime = rive_etime ; st_step_flag%ribl = rive_stepflag
        st_ribl%inttype = intri_type
#ifdef MPI_MSG
        riblv = rive_view
#endif
        if (rivet(i) == in_type(7)) then
          st_ribl%fnum = intri_num ; st_ribl%intpath = trim(adjustl(intri_path))
          st_ribl%intstep = intri_step ; st_ribl%intfnum = riven(i)
        else
          st_ribl%fnum = riven(i)
        end if
      case(4)
        st_ride%multi = rive_multi ; st_ride%totn = rive_totn
        st_ride%etime = rive_etime ; st_step_flag%ride = rive_stepflag
        st_ride%inttype = intri_type
#ifdef MPI_MSG
        ridev = rive_view
#endif
        if (rivet(i) == in_type(7)) then
          st_ride%fnum = intri_num ; st_ride%intpath = trim(adjustl(intri_path))
          st_ride%intstep = intri_step ; st_ride%intfnum = riven(i)
        else
          st_ride%fnum = riven(i)
        end if
      case(5)
        st_riwi%multi = rive_multi ; st_riwi%totn = rive_totn
        st_riwi%etime = rive_etime ; st_step_flag%riwi = rive_stepflag
        st_riwi%inttype = intri_type
#ifdef MPI_MSG
        riwiv = rive_view
#endif
        if (rivet(i) == in_type(7)) then
          st_riwi%fnum = intri_num ; st_riwi%intpath = trim(adjustl(intri_path))
          st_riwi%intstep = intri_step ; st_riwi%intfnum = riven(i)
        else
          st_riwi%fnum = riven(i)
        end if
      case(6)
        st_rile%multi = rive_multi ; st_rile%totn = rive_totn
        st_rile%etime = rive_etime ; st_step_flag%rile = rive_stepflag
        st_rile%inttype = intri_type
#ifdef MPI_MSG
        rilev = rive_view
#endif
        if (rivet(i) == in_type(7)) then
          st_rile%fnum = intri_num ; st_rile%intpath = trim(adjustl(intri_path))
          st_rile%intstep = intri_step ; st_rile%intfnum = riven(i)
        else
          st_rile%fnum = riven(i)
        end if
      end select

    end do

    st_rivf_type%wlev = riwl_type ; st_rivf_type%wdep = riwd_type
    st_rivf_type%blev = ribl_type ; st_rivf_type%dept = ride_type
    st_rivf_type%widt = riwi_type ; st_rivf_type%leng = rile_type

    deallocate(riv_txt_type)
    deallocate(rive_mask)


  end subroutine open_in_rivef

  subroutine open_in_lakef(lake_path, nohv, hv, lawlv, lawdv, lablv, laarv)
  !***************************************************************************************
  ! open_in_lakef -- Open input lake file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_lakf_type, st_lawl, st_lawd, st_labl, st_laar
    ! -- inout
    character(*), intent(in) :: lake_path
    integer(I4), intent(in) :: nohv, hv
    integer(I4), intent(inout) :: lawlv, lawdv, lablv, laarv
    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: lake_fnum, lake_totn, intla_type, intla_num, lake_stepflag
    integer(I4) :: laken(4), laket(4)
    integer(I4) :: lake_view
    integer(I4), allocatable :: lak_txt_type(:)
    real(SP) :: intla_end, lake_multi, lake_etime, intla_step
    character(CHALEN) :: lawl_path, lawd_path, labl_path, laar_path, intla_path
    character(TIMELEN) :: lawl_unit, lawd_unit, labl_unit, laar_unit
    character(:), allocatable :: mess_wl, mess_wd, mess_bl, mess_ar
    character(:), allocatable :: path_lake, mess_lake, err_mes, unit_lake
    integer(I4) :: lawl_type, lawd_type, labl_type, laar_type
    logical, allocatable :: lake_mask(:)
    namelist/inlake_type/lawl_type, lawd_type, labl_type, laar_type
    namelist/inlake_path/lawl_path, lawd_path, labl_path, laar_path
    namelist/inlake_unit/lawl_unit, lawd_unit, labl_unit, laar_unit
    !-------------------------------------------------------------------------------------
    lawl_type = 0 ; lawd_type = 0 ; labl_type = 0 ; laar_type = 0
    lawl_path = "" ; lawd_path = "" ; labl_path = "" ; laar_path = ""
    lawl_unit = "" ; lawd_unit = "" ; labl_unit = "" ; laar_unit = ""
    lawlv = 0 ; lawdv = 0 ; lablv = 0 ; laarv = 0

    if (my_rank == 0) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, lake_path, "input lake", lake_fnum)
      read(unit=lake_fnum,nml=inlake_type,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading file type section in lake file.")
      end if
      read(unit=lake_fnum,nml=inlake_path,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading file path section in lake file.")
      end if
      read(unit=lake_fnum,nml=inlake_unit,iostat=ierr)
      if (ierr /= 0) then
        call write_err_stop("While reading file unit section in lake file.")
      end if
      call close_file(lake_fnum)
    end if

    mess_wl = "input lake water level"
    mess_wd = "input lake water depth"
    mess_bl = "input lake bottom level"
    mess_ar = "input lake area"

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast file (file)
        call bcast_file(lawl_type, lawl_path, lawl_unit, mess_wl)
        call bcast_file(lawd_type, lawd_path, lawd_unit, mess_wd)
        call bcast_file(labl_type, labl_path, labl_unit, mess_bl)
        call bcast_file(laar_type, laar_path, laar_unit, mess_ar)
    end if
#endif

    laket(1) = lawl_type ; laket(2) = lawd_type ; laket(3) = labl_type
    laket(4) = laar_type
    laken(1) = 0 ; laken(2) = 0 ; laken(3) = 0 ; laken(4) = 0

    allocate(lak_txt_type(3), lake_mask(3))
    lak_txt_type(:) = [in_type(1:3)]

    do i = 1, 4
      lake_totn = 0
      lake_mask(:) = (laket(i) /= lak_txt_type(:))
      if (all(lake_mask) .and. laket(i) /= in_type(4) .and. laket(i) /= in_type(7)) then
        cycle
      end if

      lake_view = 0 ; intla_type = 0 ; err_mes = "" ; path_lake = "" ; mess_lake = ""

      select case (i)
      case(1)
        path_lake = trim(adjustl(lawl_path)) ; mess_lake = mess_wl ; unit_lake = lawl_unit
        deallocate(mess_wl)
      case(2)
        path_lake = trim(adjustl(lawd_path)) ; mess_lake = mess_wd ; unit_lake = lawd_unit
        deallocate(mess_wd)
      case(3)
        path_lake = trim(adjustl(labl_path)) ; mess_lake = mess_bl ; unit_lake = labl_unit
        deallocate(mess_bl)
      case(4)
        path_lake = trim(adjustl(laar_path)) ; mess_lake = mess_ar ; unit_lake = laar_unit
        deallocate(mess_ar)
      case default
        path_lake = ""
      end select

      if (any(laket(i) == lak_txt_type(:))) then
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, path_lake, mess_lake, laken(i))
        end if
      else if (laket(i) == in_type(4)) then
#ifdef MPI_MSG
        ! -- Open mpi read file (mpi_read_file)
          call open_mpi_read_file(1, 1, path_lake, mess_lake, laken(i))
        if (len_trim(adjustl(unit_lake)) == 0) then
          lake_view = nohv ; intla_type = 0
        else
          lake_view = hv ; intla_type = laket(i)
        end if
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(laken(i), lake_view, mess_lake)
#else
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, path_lake, mess_lake, laken(i))
        lake_view = nohv ; lake_view = hv
#endif
      else if (laket(i) == in_type(7)) then
        err_mes = mess_lake//" time interval"
        if (my_rank == 0) then
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 1, path_lake, mess_lake, laken(i))
          read(unit=laken(i),fmt=*,iostat=ierr) intla_type, intla_step, intla_end
          if (ierr /= 0) then
            call write_err_stop("While reading list file in "//mess_lake//" file.")
          end if
          read(unit=laken(i),fmt='(a)',iostat=ierr) intla_path
          if (ierr /= 0) then
            call write_err_stop("While reading list file in "//mess_lake//" file.")
          end if
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast extra setting (extr_set)
            call bcast_extr_set(intla_type, intla_step, intla_end, intla_path, err_mes)
        end if
#endif
        if (intla_type == in_type(3)) then
          if (my_rank == 0) then
            ! -- Open new read text file (new_rtxt)
              call open_new_rtxt(1, 1, trim(adjustl(intla_path)), err_mes, intla_num)
          end if
        else if (intla_type == in_type(4)) then
#ifdef MPI_MSG
          ! -- Open mpi read file (mpi_read_file)
            call open_mpi_read_file(1, 1, trim(adjustl(intla_path)), err_mes, intla_num)
          lake_view = nohv
          ! -- Set real4 file view (real4_fview)
            call set_real4_fview(intla_num, lake_view, err_mes)
#else
          ! -- Open new read binary file (new_rbin)
            call open_new_rbin(1, 1, trim(adjustl(intla_path)), err_mes, intla_num)
          lake_view = nohv
#endif
        else
          if (my_rank == 0) then
            call write_err_stop("Specified wrong number in "//err_mes//" file.")
          end if
        end if
      end if

      if (len_trim(adjustl(unit_lake)) == 0) then
        if (my_rank == 0) then
          err_mes = "Not specified time unit in "//mess_lake//". Same value is used in simulation"
          call write_logf(err_mes)
        end if
        lake_multi = SINFI ; lake_etime = SINFI
      else
        ! -- Convert unit (unit)
          call conv_unit(my_rank, unit_lake, mess_lake, st_sim%sta_date, lake_multi)
      end if

      lake_totn = st_grid%nx*st_grid%ny

      if (lake_multi /= SINFI) then
        if (any(laken(i) == lak_txt_type(:))) then
          if (my_rank == 0) then
            call skip_file(laket(i), laken(i), mess_lake, lake_multi, lake_totn, lake_etime)
          end if
#ifdef MPI_MSG
          if (pro_totn /= 1) then
            ! -- Bcast scalar value (val)
              call bcast_val(lake_totn, mess_lake//" total number")
          end if
#endif
        else if (laket(i) == in_type(4)) then
#ifdef MPI_MSG
          call skip_mpi_file(laket(i), laken(i), lake_view, mess_lake, lake_multi, lake_etime)
#else
          call skip_file(laket(i), laken(i), mess_lake, lake_multi, lake_totn, lake_etime)
#endif
        else if (laket(i) == in_type(7)) then
#ifdef MPI_MSG
          if (intla_type == in_type(3)) then
            if (my_rank == 0) then
              call skip_file_int(intla_type, laken(i), intla_num, mess_lake, lake_multi, intla_step, intla_end, lake_etime)
            end if
            if (pro_totn /= 1) then
              ! -- Bcast scalar value (val)
                call bcast_val(intla_num, mess_lake//" file number")
            end if
          else if (intla_type == in_type(4)) then
            call skip_mpi_file_int(laken(i), mess_lake, lake_multi, intla_step, intla_end, intla_num, lake_etime)
          end if
#else
          call skip_file_int(intla_type, laken(i), intla_num, mess_lake, lake_multi, intla_step, intla_end, lake_etime)
#endif
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(lake_etime, mess_lake//" end time value")
        end if
#endif
      else
        if (laket(i) == in_type(1) .or. laket(i) == in_type(2)) then
          if (my_rank == 0) then
            read(unit=laken(i),fmt=*,iostat=ierr) lake_totn
            if (ierr /= 0) then
              err_mes = "While reading the total number in "//mess_lake//" file."
              call write_err_stop(err_mes)
            end if
          end if
#ifdef MPI_MSG
          ! -- Bcast scalar value (val)
            call bcast_val(lake_totn, mess_lake//" total number")
#endif
        end if
      end if

      if (lake_etime > st_sim%end_time) then
        lake_etime = st_sim%end_time
      end if

      lake_stepflag = 0

      select case (i)
      case(1)
        st_lawl%multi = lake_multi ; st_lawl%totn = lake_totn
        st_lawl%etime = lake_etime ; st_step_flag%lawl = lake_stepflag
        st_lawl%inttype = intla_type
#ifdef MPI_MSG
        lawlv = lake_view
#endif
        if (laket(i) == in_type(7)) then
          st_lawl%fnum = intla_num ; st_lawl%intpath = trim(adjustl(intla_path))
          st_lawl%intstep = intla_step ; st_lawl%intfnum = laken(i)
        else
          st_lawl%fnum = laken(i)
        end if
      case(2)
        st_lawd%multi = lake_multi ; st_lawd%totn = lake_totn
        st_lawd%etime = lake_etime ; st_step_flag%lawd = lake_stepflag
        st_lawd%inttype = intla_type
#ifdef MPI_MSG
        lawdv = lake_view
#endif
        if (laket(i) == in_type(7)) then
          st_lawd%fnum = intla_num ; st_lawd%intpath = trim(adjustl(intla_path))
          st_lawd%intstep = intla_step ; st_lawd%intfnum = laken(i)
        else
          st_lawd%fnum = laken(i)
        end if
      case(3)
        st_labl%multi = lake_multi ; st_labl%totn = lake_totn
        st_labl%etime = lake_etime ; st_step_flag%labl = lake_stepflag
        st_labl%inttype = intla_type
#ifdef MPI_MSG
        lablv = lake_view
#endif
        if (laket(i) == in_type(7)) then
          st_labl%fnum = intla_num ; st_labl%intpath = trim(adjustl(intla_path))
          st_labl%intstep = intla_step ; st_labl%intfnum = laken(i)
        else
          st_labl%fnum = laken(i)
        end if
      case(4)
        st_laar%multi = lake_multi ; st_laar%totn = lake_totn
        st_laar%etime = lake_etime ; st_step_flag%laar = lake_stepflag
        st_laar%inttype = intla_type
#ifdef MPI_MSG
        laarv = lake_view
#endif
        if (laket(i) == in_type(7)) then
          st_laar%fnum = intla_num ; st_laar%intpath = trim(adjustl(intla_path))
          st_laar%intstep = intla_step ; st_laar%intfnum = laken(i)
        else
          st_laar%fnum = laken(i)
        end if
      end select

    end do

    st_lakf_type%wlev = lawl_type ; st_lakf_type%wdep = lawd_type
    st_lakf_type%blev = labl_type ; st_lakf_type%area = laar_type

    deallocate(lak_txt_type)
    deallocate(lake_mask)

  end subroutine open_in_lakef

  subroutine open_in_massf(in_mass_type, in_mass_path, view_mass)
  !***************************************************************************************
  ! open_in_massf -- Open input massbalance file for output
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4) , intent(in) :: in_mass_type
    character(*), intent(in) :: in_mass_path
    integer(I4), intent(in), optional :: view_mass
    ! -- local
    integer(I4) :: temp_view
    integer(I4), allocatable :: txt_mass_type(:)
    logical, allocatable :: txt_mass_mask(:)
    !-------------------------------------------------------------------------------------
    allocate(txt_mass_type(3), txt_mass_mask(3))
    txt_mass_type(:) = [in_type(1), in_type(3), in_type(5)]

    if (present(view_mass)) then
      temp_view = view_mass
    end if

    if (any(in_mass_type == txt_mass_type)) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, in_mass_path, "input massbalance", inmas_fnum)
      end if
    else if (in_mass_type == in_type(4) .or. in_mass_type == in_type(6)) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(1, 1, in_mass_path, "input massbalance", inmas_fnum)
      ! -- Set int4 file view (int4_fview)
        call set_int4_fview(inmas_fnum, temp_view, "input massbalance")
#else
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, in_mass_path, "input massbalance", inmas_fnum)
#endif
    end if

    deallocate(txt_mass_type)
    deallocate(txt_mass_mask)

  end subroutine open_in_massf

  subroutine open_out_convf(out_conv_path, out_conv_fnum)
  !***************************************************************************************
  ! open_out_convf -- Open output convergence file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    character(*), intent(in) :: out_conv_path
    integer(I4), intent(inout) :: out_conv_fnum
    ! -- local

    !-------------------------------------------------------------------------------------
    ! -- Open new read text file (new_rtxt)
      call open_new_wtxt(out_conv_path, "output convergence", out_conv_fnum)

  end subroutine open_out_convf

  subroutine open_out_binf(out_file, out_tint, out_path, out_unit, out_chra, out_trel)
  !***************************************************************************************
  ! open_out_binf -- Open output binary file
  !***************************************************************************************
    ! -- modules
    use utility_module, only: open_new_wbin
#ifdef MPI_MSG
    use mpi_read, only: open_mpi_write_file
#endif
    ! -- inout
    integer(I4), intent(inout) :: out_file, out_tint
    character(*), intent(inout) :: out_path, out_unit, out_chra
    real(SP), intent(out) :: out_trel
    ! -- local
#ifdef MPI_MSG
    integer(I4) :: out_fh
#endif
    !-------------------------------------------------------------------------------------
#ifdef MPI_MSG
    if (pro_totn /= 1) then
    ! -- Bcast file (file)
      call bcast_file(out_path, out_unit, out_chra)
    ! -- Bcast scalar value (val)
      call bcast_val(out_tint, out_chra//" time number")
    end if
    ! -- Open mpi write file (mpi_write file)
      call open_mpi_write_file(out_path, out_chra, out_fh)
    out_file = out_fh
#else
    ! -- Open new write binary file (new_wbin)
      call open_new_wbin(out_path, out_chra, out_file)
#endif

    if (out_tint /= 0) then
      ! -- Convert unit (unit)
        call conv_unit(my_rank, out_unit, out_chra, st_sim%sta_date, out_trel)
    end if

    out_trel = out_trel*out_tint

  end subroutine open_out_binf

  subroutine open_out_massf(out_mass_time, out_mass_path, out_mass_unit, out_mass_fnum)
  !***************************************************************************************
  ! open_out_massf -- Open output massbalance file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_out_step
    ! -- inout
    integer(I4), intent(inout) :: out_mass_time
    character(*), intent(inout) :: out_mass_path, out_mass_unit
    integer(I4), intent(out) :: out_mass_fnum
    ! -- local
    character(:), allocatable :: out_chra
    !-------------------------------------------------------------------------------------
    out_chra = "output massbalance"
    if (my_rank == 0) then
      ! -- Open new write text file (new_wtxt)
        call open_new_wtxt(out_mass_path, out_chra, out_mass_fnum)
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast file (file)
        call bcast_file(out_mass_path, out_mass_unit, out_chra)
      ! -- Bcast scalar value (val)
        call bcast_val(out_mass_fnum, out_chra//" file number")
      ! -- Bcast scalar value (val)
        call bcast_val(out_mass_time, out_chra//" time number")
    end if
#endif

    if (out_mass_time /= 0) then
      ! -- Convert unit (unit)
        call conv_unit(my_rank, out_mass_unit, out_chra, st_sim%sta_date, st_out_step%mass)
    end if

    st_out_step%mass = st_out_step%mass*out_mass_time

  end subroutine open_out_massf

end module open_file
