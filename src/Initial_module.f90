module initial_module
  ! -- modules
  use kind_module, only: I4
  use constval_module, only: TIMELEN
  use types_module, only: sim_set, ctrl_set, schm_set, grid_set, ftype_in, retn_in, parm_in
  use types_module, only: geog_in
  use types_module, only: rive_in, lake_in, path_in, retn_path, parm_path, geog_path, retn_fnum
  use types_module, only: parm_fnum, geog_fnum, unit_in, clas_set, retn_set, parm_set, init_set
  use types_module, only: seal_set, surfb_set, well_set, surfw_set, step_flag, ftype_out
  use types_module, only: path_out, unit_out, out_time, out_step
  use utility_module, only: log_fnum

  implicit none
  private
  public :: init_msg, init_log

  !file type
  integer(I4), allocatable, public :: in_type(:), out_type(:)
  !time unit
  character(TIMELEN), allocatable, public :: unit_list(:)

  ! main file
  type(sim_set), public :: st_sim

  !input solution file
  type(ctrl_set), public :: st_ctrl

  !input scheme file
  type(schm_set), public :: st_schm

  type(grid_set), public :: st_grid

  ! input ftype, path, unit
  type(ftype_in), public :: st_in_type
  type(retn_in), public :: st_retf_type
  type(parm_in), public :: st_parf_type
  type(geog_in), public :: st_geof_type
  type(rive_in), public :: st_rivf_type
  type(lake_in), public :: st_lakf_type

  type(path_in), public :: st_in_path
  type(retn_path), public :: st_retn_path
  type(parm_path), public :: st_parm_path
  type(geog_path), public :: st_geog_path

  type(retn_fnum), public :: st_retn_fnum
  type(parm_fnum), public :: st_parm_fnum
  type(geog_fnum), public :: st_geog_fnum

  type(unit_in), public :: st_in_unit

  ! input structure
  type(clas_set), public :: st_clas
  type(retn_set), public :: st_retn
  type(parm_set), public :: st_parm
  type(init_set), public :: st_init
  type(seal_set), public :: st_seal
  type(surfb_set), public :: st_rech, st_prec, st_evap
  type(well_set), public :: st_well
  type(surfw_set), public :: st_riwl, st_riwd, st_ribl, st_ride, st_riwi, st_rile
  type(surfw_set), public :: st_lawl, st_lawd, st_labl, st_laar

  type(step_flag), public :: st_step_flag

  ! output file number, path, unit
  type(ftype_out), public :: st_out_type
  type(path_out), public :: st_out_path
  type(unit_out), public :: st_out_unit
  type(out_time), public :: st_out_time
  type(out_step), public :: st_out_step

  ! -- local

  contains

  subroutine init_msg(in_stime)
  !*********************************************************************************************
  ! init_msg -- Initialize msg
  !*********************************************************************************************
    ! -- module
    use utility_module, only: st_mpi
#ifdef MPI_MSG
    use mpi_initfin, only: init_mpi
#endif
    ! -- inout
    integer(I4), intent(in) :: in_stime(:)
    ! -- local

    !-------------------------------------------------------------------------------------------
    st_mpi%totn = 1 ; st_mpi%rank = 0 ; st_mpi%comm = 0
    if (st_mpi%rank == 0) then
      ! -- Initialize log file (log)
        call init_log(in_stime)
    end if

#ifdef MPI_MSG
    ! -- Initialize mpi (mpi)
      call init_mpi(log_fnum, st_mpi)
#endif

#ifdef ICI
    ! -- Initialize ici (ici)
      call init_ici()
#endif

    ! -- Initialize variables (var)
      call init_var()

    ! -- Initialize OpenMP (omp)
    !$ call init_omp()

  end subroutine init_msg

  subroutine init_log(in_stime)
  !*********************************************************************************************
  ! init_log -- Initialize log file
  !*********************************************************************************************
    ! -- module
    use utility_module, only: open_new_wtxt
    ! -- inout
    integer(I4), intent(in) :: in_stime(:)
    ! -- local
    integer(I4) :: i
    character(:), allocatable :: log_file
    10 format(/"Run start date and time(yyyy/mm/dd hh:mm:ss) : ",i4,"/",i2.2,"/",i2.2,1x,&
              i2,":",i2.2,":",i2.2,/)
    !-------------------------------------------------------------------------------------------
    allocate(character(0) :: log_file)
    log_file = 'msg123_log.txt'
    ! -- Open new read text file (new_rtxt)
      call open_new_wtxt(log_file, "msg123 log", log_fnum)

    write(log_fnum,10) (in_stime(i), i = 1, 3), (in_stime(i), i = 5, 7)

    deallocate(log_file)

  end subroutine init_log

  subroutine init_var()
  !*********************************************************************************************
  ! init_var -- Initialize variables
  !*********************************************************************************************
    ! -- module
    use kind_module, only: SP, DP
    use constval_module, only: CHALEN, INF_SPEC, INF_CLAS, INF_POIN, INF_2DTX, INF_2DBI
    use constval_module, only: INF_3DTX, INF_3DBI, INF_EXTR, OUTF_TABL, OUTF_2DBI, OUTF_3DBI
    use constval_module, only: SZERO, SINFI, DZERO, DONE, MACHI_EPS
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------------
    ! input file type
    st_in_type%grid = -1 ; st_in_type%retn = -1 ; st_in_type%parm = -1
    st_in_type%geog = -1 ; st_in_type%init = -1 ; st_in_type%seal = -1
    st_in_type%rech = -1 ; st_in_type%well = -1 ; st_in_type%weks = -1
    st_in_type%weke = -1 ; st_in_type%rive = -1 ; st_in_type%lake = -1
    st_in_type%prec = -1 ; st_in_type%evap = -1 ; st_in_type%wtab = -1
    st_in_type%mass = -1

    st_retf_type%vana = -1 ; st_retf_type%vann = -1 ; st_retf_type%resi = -1
    st_parf_type%pakx = -1 ; st_parf_type%paky = -1 ; st_parf_type%pakz = -1
    st_parf_type%pass = -1 ; st_parf_type%pats = -1
    st_geof_type%geoz = -1 ; st_geof_type%geor = -1 ; st_geof_type%geoa = -1
    st_rivf_type%wlev = -1 ; st_rivf_type%wdep = -1 ; st_rivf_type%blev = -1
    st_rivf_type%dept = -1 ; st_rivf_type%widt = -1 ; st_rivf_type%leng = -1
    st_lakf_type%wlev = -1 ; st_lakf_type%wdep = -1 ; st_lakf_type%blev = -1
    st_lakf_type%area = -1

    ! total number of each variable
    st_clas%totn = 0 ; st_retn%totn = 0 ; st_parm%totn = 0 ; st_seal%totn = 0
    st_rech%totn = 0 ; st_well%totn = 0 ; st_prec%totn = 0 ; st_evap%totn = 0
    st_riwl%totn = 0 ; st_riwd%totn = 0 ; st_ribl%totn = 0 ; st_ride%totn = 0
    st_riwi%totn = 0 ; st_rile%totn = 0 ; st_lawl%totn = 0 ; st_lawd%totn = 0
    st_labl%totn = 0 ; st_laar%totn = 0

    ! end time of each variable
    st_seal%etime = SINFI ; st_rech%etime = SINFI ; st_well%etime = SINFI
    st_prec%etime = SINFI ; st_evap%etime = SINFI ; st_riwl%etime = SINFI
    st_riwd%etime = SINFI ; st_ribl%etime = SINFI ; st_ride%etime = SINFI
    st_riwi%etime = SINFI ; st_rile%etime = SINFI ; st_lawl%etime = SINFI
    st_lawd%etime = SINFI ; st_labl%etime = SINFI ; st_laar%etime = SINFI

    ! step flag of each variable
    st_step_flag%seal = 0 ; st_step_flag%rech = 0 ; st_step_flag%well = 0
    st_step_flag%prec = 0 ; st_step_flag%evap = 0 ; st_step_flag%riwl = 0
    st_step_flag%riwd = 0 ; st_step_flag%ribl = 0 ; st_step_flag%ride = 0
    st_step_flag%riwi = 0 ; st_step_flag%rile = 0 ; st_step_flag%lawl = 0
    st_step_flag%lawd = 0 ; st_step_flag%labl = 0 ; st_step_flag%laar = 0

    ! output time of each variable
    st_out_step%head = SZERO ; st_out_step%rest = SZERO ; st_out_step%srat = SZERO
    st_out_step%wtab = SZERO ; st_out_step%mass = SZERO ; st_out_step%velc = SZERO
    st_out_step%rivr = SZERO ; st_out_step%lakr = SZERO ; st_out_step%sufr = SZERO
    st_out_step%dunr = SZERO ; st_out_step%well = SZERO ; st_out_step%rech = SZERO

    ! output type of each variable
    st_out_type%srat = 0 ; st_out_type%wtab = 0 ; st_out_type%mass = 0
    st_out_type%velc = 0 ; st_out_type%rivr = 0 ; st_out_type%lakr = 0
    st_out_type%sufr = 0 ; st_out_type%dunr = 0 ; st_out_type%seal = 0
    st_out_type%rech = 0 ; st_out_type%well = 0 ; st_out_type%calg = 0

    ! input solution file
    st_ctrl%tstep_type = 0 ; st_ctrl%maxout_iter = 20 ; st_ctrl%picard_iter = 0
    st_ctrl%maxinn_iter = 10 ; st_ctrl%precon_type = 0 ; st_ctrl%nlevel = 0
    st_ctrl%maxvcy_iter = 1 ; st_ctrl%amg_nlevel = 5 ; st_ctrl%max_sweep = 1
    st_ctrl%criteria = 1.00E-03_DP ; st_ctrl%errtol = DZERO
    st_ctrl%res_abs_tol = DZERO ; st_ctrl%res_rel_tol = DZERO
    st_ctrl%dilu_shift = DZERO ; st_ctrl%dsat_max = DZERO ; st_ctrl%expd_type = 0
    st_ctrl%conv_type = 0
    st_ctrl%newper = MACHI_EPS ; st_ctrl%newper_inv = DONE/st_ctrl%newper
    st_ctrl%jac_omega = 0.67_SP ; st_ctrl%amg_theta = 0.05_SP

    ! input scheme file
    st_schm%krpos_type = 0 ; st_schm%stor_type = 0
    st_schm%abyd_type = 0 ; st_schm%abyd_ratio = DZERO

    ! time unit
    unit_list = ["SEC", "MIN", "HOU", "DAY", "YEA"]

    !file type
    allocate(in_type(0:7), out_type(3))
    in_type(0:3) = [INF_SPEC, INF_CLAS, INF_POIN, INF_2DTX]
    in_type(4:7) = [INF_2DBI, INF_3DTX, INF_3DBI, INF_EXTR]
    out_type(:) = [OUTF_TABL, OUTF_2DBI, OUTF_3DBI]

    st_sim%sim_name = repeat(' ', CHALEN)
    st_sim%cal_unit = repeat(' ', TIMELEN)
    st_sim%reg_name = repeat(' ', CHALEN)
    st_sim%inact_name = repeat(' ', CHALEN)

    st_in_path%grid = repeat(' ', CHALEN)
    st_in_path%retn = repeat(' ', CHALEN)
    st_in_path%parm = repeat(' ', CHALEN)
    st_in_path%geog = repeat(' ', CHALEN)
    st_in_path%init = repeat(' ', CHALEN)
    st_in_path%seal = repeat(' ', CHALEN)
    st_in_path%rech = repeat(' ', CHALEN)
    st_in_path%well = repeat(' ', CHALEN)
    st_in_path%weks = repeat(' ', CHALEN)
    st_in_path%weke = repeat(' ', CHALEN)
    st_in_path%rive = repeat(' ', CHALEN)
    st_in_path%lake = repeat(' ', CHALEN)
    st_in_path%prec = repeat(' ', CHALEN)
    st_in_path%evap = repeat(' ', CHALEN)
    st_in_path%mass = repeat(' ', CHALEN)

    st_in_unit%init = repeat(' ', TIMELEN)
    st_in_unit%seal = repeat(' ', TIMELEN)
    st_in_unit%rech = repeat(' ', TIMELEN)
    st_in_unit%well = repeat(' ', TIMELEN)
    st_in_unit%prec = repeat(' ', TIMELEN)
    st_in_unit%evap = repeat(' ', TIMELEN)

    st_out_path%conv = repeat(' ', CHALEN)
    st_out_path%head = repeat(' ', CHALEN)
    st_out_path%rest = repeat(' ', CHALEN)
    st_out_path%srat = repeat(' ', CHALEN)
    st_out_path%wtab = repeat(' ', CHALEN)
    st_out_path%mass = repeat(' ', CHALEN)
    st_out_path%velx = repeat(' ', CHALEN)
    st_out_path%vely = repeat(' ', CHALEN)
    st_out_path%velz = repeat(' ', CHALEN)
    st_out_path%rivr = repeat(' ', CHALEN)
    st_out_path%lakr = repeat(' ', CHALEN)
    st_out_path%sufr = repeat(' ', CHALEN)
    st_out_path%dunr = repeat(' ', CHALEN)
    st_out_path%seal = repeat(' ', CHALEN)
    st_out_path%well = repeat(' ', CHALEN)
    st_out_path%rech = repeat(' ', CHALEN)
    st_out_path%calg = repeat(' ', CHALEN)

    st_out_unit%head = repeat(' ', TIMELEN)
    st_out_unit%rest = repeat(' ', TIMELEN)
    st_out_unit%srat = repeat(' ', TIMELEN)
    st_out_unit%wtab = repeat(' ', TIMELEN)
    st_out_unit%mass = repeat(' ', TIMELEN)
    st_out_unit%velc = repeat(' ', TIMELEN)
    st_out_unit%rivr = repeat(' ', TIMELEN)
    st_out_unit%lakr = repeat(' ', TIMELEN)
    st_out_unit%sufr = repeat(' ', TIMELEN)
    st_out_unit%dunr = repeat(' ', TIMELEN)
    st_out_unit%seal = repeat(' ', TIMELEN)
    st_out_unit%well = repeat(' ', TIMELEN)
    st_out_unit%rech = repeat(' ', TIMELEN)

  end subroutine init_var

  !$ subroutine init_omp()
  !*********************************************************************************************
  ! init_omp -- Initialize OpenMP
  !*********************************************************************************************
    ! -- module
    !$ use omp_lib
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------------
    ! set internal control variables
    !$ call OMP_SET_SCHEDULE(1, 0)
    !$ call OMP_SET_DYNAMIC(.false.)
    !$ call OMP_SET_MAX_ACTIVE_LEVELS(1)

  !$ end subroutine init_omp

end module initial_module
