module initial_module
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: VARLEN, TIMELEN
  use utility_module, only: log_fnum

  implicit none
  private
  public :: init_msg, init_log
  integer(I4), public :: pro_totn = 1, my_rank = 0
  integer(I4), public :: noclas_flag = 0

  ! main file
  type :: sim_set
    integer(I4) :: sim_type, res_type, reg_type, reg_neib
    integer(I4) :: sta_date(6), end_date(6)
    character(:), allocatable :: sim_name, cal_unit, reg_name, inact_name
    real(SP) :: end_time, ini_step, max_step, inc_fact, cal_fact
  end type sim_set
  type(sim_set), public :: st_sim

  type :: grid_set
    integer(I4) :: nx, ny, nz, nxyz
    integer(I4) :: fnum
  end type grid_set
  type(grid_set), public :: st_grid

  ! input ftype, path, unit
  type :: ftype_in
    integer(I4) :: grid, retn, parm, geog, init
    integer(I4) :: seal, rech, well, weks, weke, rive, lake
    integer(I4) :: prec, evap, wtab, mass
  end type ftype_in
  type(ftype_in), public :: st_in_type

  type :: retn_in
    integer(I4) :: vana, vann, resi
  end type retn_in
  type(retn_in), public :: st_retf_type

  type :: parm_in
    integer(I4) :: pakx, paky, pakz, pass, pats
  end type parm_in
  type(parm_in), public :: st_parf_type

  type :: geog_in
    integer(I4) :: geoz, geor, geoa
  end type geog_in
  type(geog_in), public :: st_geof_type

  type :: rive_in
    integer(I4) :: wlev, wdep, blev, dept, widt, leng
  end type rive_in
  type(rive_in), public :: st_rivf_type

  type :: lake_in
    integer(I4) :: wlev, wdep, blev, area
  end type lake_in
  type(lake_in), public :: st_lakf_type

  type :: path_in
    character(:), allocatable :: grid, retn, parm, geog, init
    character(:), allocatable :: seal, rech, well, weks, weke, rive, lake
    character(:), allocatable :: prec, evap, mass
  end type path_in
  type(path_in), public :: st_in_path

  type :: retn_path
    character(:), allocatable :: vana, vann, resi
  end type retn_path
  type(retn_path), public :: st_retn_path

  type :: parm_path
    character(:), allocatable :: pakx, paky, pakz, pass, pats
  end type parm_path
  type(parm_path), public :: st_parm_path

  type :: geog_path
    character(:), allocatable :: geoz, geor, geoa
  end type geog_path
  type(geog_path), public :: st_geog_path

  type :: retn_fnum
    integer(I4) :: vana, vann, resi
  end type retn_fnum
  type(retn_fnum), public :: st_retn_fnum

  type :: parm_fnum
    integer(I4) :: pakx, paky, pakz, pass, pats
  end type parm_fnum
  type(parm_fnum), public :: st_parm_fnum

  type :: geog_fnum
    integer(I4) :: geoz, geor, geoa
  end type geog_fnum
  type(geog_fnum), public :: st_geog_fnum

  type :: unit_in
    character(:), allocatable :: init, seal, rech, well, prec, evap
  end type unit_in
  type(unit_in), public :: st_in_unit

  ! input structure
  type :: clas
    integer(I4) :: totn, fnum
    character(VARLEN), allocatable :: name(:)
    integer(I4), allocatable :: num(:)
    integer(I4), allocatable :: i(:,:), j(:,:), k(:,:)
  end type clas
  type(clas), public :: st_clas

  type :: retn
    integer(I4) :: totn, fnum
    character(VARLEN), allocatable :: name(:)
    real(SP), allocatable :: a(:), n(:), r(:)
  end type retn
  type(retn), public :: st_retn

  type :: parm
    integer(I4) :: totn, fnum
    character(VARLEN), allocatable :: name(:)
    real(SP), allocatable :: ksx(:), ksy(:), ksz(:), ss(:), ts(:)
  end type parm
  type(parm), public :: st_parm

  type :: init
    integer(I4) :: fnum
    real(SP) :: multi, depth, rest_time
  end type init
  type(init), public :: st_init

  type :: seal
    integer(I4) :: totn, fnum
    real(SP) :: etime, multi
    character(VARLEN), allocatable :: name(:)
    integer(I4), allocatable :: i(:), j(:), k(:)
    real(SP), allocatable :: value(:)
  end type seal
  type(seal), public :: st_seal

  type, public :: st_surfb
    integer(I4) :: totn, fnum
    real(SP) :: etime, uni_conv, multi
    character(VARLEN), allocatable :: name(:)
    real(SP), allocatable :: value(:)
  end type st_surfb

  type(st_surfb), public :: st_rech, st_prec, st_evap

  type :: well
    integer(I4) :: totn, fnum
    real(SP) :: etime, uni_conv, multi
    integer(I4), allocatable :: i(:), j(:), k(:), ij(:), ks(:), ke(:)
    real(SP), allocatable :: value(:)
  end type well
  type(well), public :: st_well

  type, public :: st_surfw
    integer(I4) :: totn, fnum, inttype, intfnum
    real(SP) :: etime, multi, intstep
    character(VARLEN), allocatable :: name(:)
    integer(I4), allocatable :: i(:), j(:)
    real(SP), allocatable :: value(:)
    character(:), allocatable :: intpath
  end type st_surfw

  type(st_surfw), public :: st_riwl, st_riwd, st_ribl, st_ride, st_riwi, st_rile
  type(st_surfw), public :: st_lawl, st_lawd, st_labl, st_laar

  type :: step_flag
    integer(I4) :: seal, rech, well, prec, evap, riwl, riwd, ribl, ride, riwi, rile
    integer(I4) :: lawl, lawd, labl, laar
  end type step_flag
  type(step_flag), public :: st_step_flag

  ! output file number, path, unit
  type :: ftype_out
    integer(I4) :: srat, wtab, mass, velc, rivr, lakr, sufr, dunr, seal, rech, well, calg
  end type ftype_out
  type(ftype_out), public :: st_out_type

  type :: path_out
    character(:), allocatable :: conv, head, rest, srat, wtab, mass, velx, vely, velz
    character(:), allocatable :: rivr, lakr, sufr, dunr
    character(:), allocatable :: seal, well, rech, calg
  end type path_out
  type(path_out), public :: st_out_path

  type :: unit_out
    character(:), allocatable :: head, rest, srat, wtab, mass, velc
    character(:), allocatable :: rivr, lakr, sufr, dunr
    character(:), allocatable :: seal, well, rech
  end type unit_out
  type(unit_out), public :: st_out_unit

  type :: out_time
    integer(I4) :: head, rest, srat, wtab, mass, velc, rivr, lakr, sufr, dunr
    integer(I4) :: seal, well, rech
  end type out_time
  type(out_time), public :: st_out_time

  type :: out_step
    real(SP) :: head, rest, srat, wtab, mass, velc, rivr, lakr, sufr, dunr
    real(SP) :: seal, well, rech
  end type out_step
  type(out_step), public :: st_out_step

  !input solution file
  integer(I4), public :: maxout_iter, maxinn_iter, precon_type
  integer(I4), public :: nlevel, maxvcy_iter, amg_nlevel, max_sweep
  real(DP), public :: criteria, errtol, newper, newper_inv
  real(SP), public :: jac_omega, amg_theta
  !time unit
  character(TIMELEN), allocatable, public :: unit_list(:)
  !file type
  integer(I4), allocatable, public :: in_type(:), out_type(:)

  ! -- local

  contains

  subroutine init_msg(in_stime)
  !***************************************************************************************
  ! init_msg -- Initialize msg
  !***************************************************************************************
    ! -- module
#ifdef MPI_MSG
    use mpi_initfin, only: init_mpi
#endif
    ! -- inout
    integer(I4), intent(in) :: in_stime(:)
    ! -- local

    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      ! -- Initialize log file (log)
        call init_log(in_stime)
    end if

#ifdef MPI_MSG
    ! -- Initialize mpi (mpi)
      call init_mpi(log_fnum, pro_totn, my_rank)
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
  !***************************************************************************************
  ! init_log -- Initialize log file
  !***************************************************************************************
    ! -- module
    use utility_module, only: open_new_wtxt
    ! -- inout
    integer(I4), intent(in) :: in_stime(:)
    ! -- local
    integer(I4) :: i, ierr
    character(:), allocatable :: log_file
    10 format(/"Run start date and time(yyyy/mm/dd hh:mm:ss) : ",i4,"/",i2.2,"/",i2.2,1x,&
              i2,":",i2.2,":",i2.2,/)
    !-------------------------------------------------------------------------------------
    log_file = 'msg123_log.txt'
    ! -- Open new read text file (new_rtxt)
      call open_new_wtxt(log_file, "msg123 log", log_fnum)

    rewind(log_fnum)
    write(log_fnum,10) (in_stime(i), i = 1, 3), (in_stime(i), i = 5, 7)

    ierr = 0
    do while(.true.)
      read(unit=log_fnum,fmt=*,iostat=ierr)
      if (ierr /= 0) then
        exit
      end if
    end do

    deallocate(log_file)

  end subroutine init_log

  subroutine init_var()
  !***************************************************************************************
  ! init_var -- Initialize variables
  !***************************************************************************************
    ! -- module
    use constval_module, only: SZERO, DZERO, DONE, SINFI, MACHI_EPS, INF_SPEC, INF_CLAS,&
                               INF_POIN, INF_2DTX, INF_2DBI, INF_3DTX, INF_3DBI,&
                               INF_EXTR, OUTF_TABL, OUTF_2DBI, OUTF_3DBI
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
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
    maxout_iter = 20 ; maxinn_iter = 10 ; precon_type = 0
    nlevel = 0 ; maxvcy_iter = 0 ; amg_nlevel = 0 ; max_sweep = 0
    criteria = 1.00E-03_DP ; errtol = DZERO
    newper = MACHI_EPS ; newper_inv = DONE/newper
    jac_omega = 0.67_SP ; amg_theta = 0.05_SP

    ! time unit
    unit_list = ["SEC", "MIN", "HOU", "DAY", "YEA"]

    !file type
    allocate(in_type(0:7), out_type(3))
    in_type(:) = [INF_SPEC, INF_CLAS, INF_POIN, INF_2DTX, INF_2DBI, INF_3DTX, INF_3DBI, INF_EXTR]
    out_type(:) = [OUTF_TABL, OUTF_2DBI, OUTF_3DBI]


  end subroutine init_var

  !$ subroutine init_omp()
  !***************************************************************************************
  ! init_omp -- Initialize OpenMP
  !***************************************************************************************
    ! -- module
    !$ use omp_lib
    ! -- inout

    ! -- local
    !$ integer(KIND=OMP_SCHED_KIND) :: kind
    !$ integer(I4) :: chunk_size
    !-------------------------------------------------------------------------------------
    ! set internal control variables

    !$ call OMP_GET_SCHEDULE(kind, chunk_size)
    !$ if (kind /= 1) then
    !$   call OMP_SET_SCHEDULE(1, 0)
    !$ end if

    !$ if (OMP_GET_DYNAMIC()) then
    !$   call OMP_SET_DYNAMIC(.false.)
    !$ end if

    !$ if (OMP_GET_NESTED()) then
    !$   call OMP_SET_NESTED(.false.)
    !$ end if

    !$ call OMP_SET_MAX_ACTIVE_LEVELS(1)

  !$ end subroutine init_omp

end module initial_module
