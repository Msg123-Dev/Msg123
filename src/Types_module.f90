module types_module
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: VARLEN

  implicit none
  private
  public :: mpi_set, sim_set, ctrl_set, geom_set, time_set, grid_set, ftype_in, retn_in,&
            parm_in, geog_in, rive_in, lake_in, path_in, retn_path, parm_path, geog_path,&
            retn_fnum, parm_fnum, geog_fnum, unit_in, clas_set, retn_set, parm_set, init_set,&
            seal_set, surfb_set, well_set, surfw_set, step_flag, ftype_out, path_out, unit_out,&
            out_time, out_step, msout_set, hydr_set, bcnd_set, forc_set, bflag_set, bcalc_set,&
            bcount_set, rlbc_set, bfview_set, bound_fview

  type :: mpi_set
    integer(I4) :: totn, rank, comm
  end type mpi_set

  type :: sim_set
    integer(I4) :: sim_type, res_type, reg_type, reg_neib
    integer(I4) :: sta_date(6), end_date(6)
    character(:), allocatable :: sim_name, cal_unit, reg_name, inact_name
    real(SP) :: end_time, ini_step, max_step, inc_fact, dec_fact, cal_fact
  end type sim_set

  type :: ctrl_set
    integer(I4) :: tstep_type, maxout_iter, picard_iter, maxinn_iter, precon_type
    integer(I4) :: nlevel, maxvcy_iter, amg_nlevel, max_sweep
    integer(I4) :: noclas_flag = 0
    real(DP) :: criteria, errtol, newper, newper_inv
    real(SP) :: jac_omega, amg_theta
  end type ctrl_set

  type :: geom_set
    real(DP), allocatable :: cell_vol(:)
    real(DP), allocatable :: surf_elev(:)
    real(DP), allocatable :: cell_top(:)
    real(DP), allocatable :: cell_cent(:)
    real(DP), allocatable :: cell_bot(:)
  end type geom_set

  type :: time_set
    integer(I4) :: now_date(6), out_iter, form_switch
    real(SP) :: current_t, inter_time, now_time
    real(DP) :: delt, delt_inv
    logical :: conv_flag
  end type time_set

  type :: grid_set
    integer(I4) :: nx, ny, nz, nxyz
    integer(I4) :: fnum
  end type grid_set

  type :: ftype_in
    integer(I4) :: grid, retn, parm, geog, init
    integer(I4) :: seal, rech, well, weks, weke, rive, lake
    integer(I4) :: prec, evap, wtab, mass
  end type ftype_in

  type :: retn_in
    integer(I4) :: vana, vann, resi
  end type retn_in

  type :: parm_in
    integer(I4) :: pakx, paky, pakz, pass, pats
  end type parm_in

  type :: geog_in
    integer(I4) :: geoz, geor, geoa
  end type geog_in

  type :: rive_in
    integer(I4) :: wlev, wdep, blev, dept, widt, leng
  end type rive_in

  type :: lake_in
    integer(I4) :: wlev, wdep, blev, area
  end type lake_in

  type :: path_in
    character(:), allocatable :: grid, retn, parm, geog, init
    character(:), allocatable :: seal, rech, well, weks, weke, rive, lake
    character(:), allocatable :: prec, evap, mass
  end type path_in

  type :: retn_path
    character(:), allocatable :: vana, vann, resi
  end type retn_path

  type :: parm_path
    character(:), allocatable :: pakx, paky, pakz, pass, pats
  end type parm_path

  type :: geog_path
    character(:), allocatable :: geoz, geor, geoa
  end type geog_path

  type :: retn_fnum
    integer(I4) :: vana, vann, resi
  end type retn_fnum

  type :: parm_fnum
    integer(I4) :: pakx, paky, pakz, pass, pats
  end type parm_fnum

  type :: geog_fnum
    integer(I4) :: geoz, geor, geoa
  end type geog_fnum

  type :: unit_in
    character(:), allocatable :: init, seal, rech, well, prec, evap
  end type unit_in

  type :: clas_set
    integer(I4) :: totn, fnum
    character(VARLEN), allocatable :: name(:)
    integer(I4), allocatable :: num(:)
    integer(I4), allocatable :: i(:,:), j(:,:), k(:,:)
  end type clas_set

  type :: retn_set
    integer(I4) :: totn, fnum
    character(VARLEN), allocatable :: name(:)
    real(SP), allocatable :: a(:), n(:), r(:)
  end type retn_set

  type :: parm_set
    integer(I4) :: totn, fnum
    character(VARLEN), allocatable :: name(:)
    real(SP), allocatable :: ksx(:), ksy(:), ksz(:), ss(:), ts(:)
  end type parm_set

  type :: init_set
    integer(I4) :: fnum
    real(SP) :: multi, depth, rest_time
  end type init_set

  type :: seal_set
    integer(I4) :: totn, fnum
    real(SP) :: etime, multi
    character(VARLEN), allocatable :: name(:)
    integer(I4), allocatable :: i(:), j(:), k(:)
    real(SP), allocatable :: value(:)
  end type seal_set

  type, public :: surfb_set
    integer(I4) :: totn, fnum
    real(SP) :: etime, uni_conv, multi
    character(VARLEN), allocatable :: name(:)
    real(SP), allocatable :: value(:)
  end type surfb_set

  type :: well_set
    integer(I4) :: totn, fnum
    real(SP) :: etime, uni_conv, multi
    integer(I4), allocatable :: i(:), j(:), k(:), ij(:), ks(:), ke(:)
    real(SP), allocatable :: value(:)
  end type well_set

  type, public :: surfw_set
    integer(I4) :: totn, fnum, inttype, intfnum
    real(SP) :: etime, multi, intstep
    character(VARLEN), allocatable :: name(:)
    integer(I4), allocatable :: i(:), j(:)
    real(SP), allocatable :: value(:)
    character(:), allocatable :: intpath
  end type surfw_set

  type :: step_flag
    integer(I4) :: seal, rech, well, prec, evap, riwl, riwd, ribl, ride, riwi, rile
    integer(I4) :: lawl, lawd, labl, laar
  end type step_flag

  type :: ftype_out
    integer(I4) :: srat, wtab, mass, velc, rivr, lakr, sufr, dunr, seal, rech, well, calg
  end type ftype_out

  type :: path_out
    character(:), allocatable :: conv, head, rest, srat, wtab, mass, velx, vely, velz
    character(:), allocatable :: rivr, lakr, sufr, dunr
    character(:), allocatable :: seal, well, rech, calg
  end type path_out

  type :: unit_out
    character(:), allocatable :: head, rest, srat, wtab, mass, velc
    character(:), allocatable :: rivr, lakr, sufr, dunr
    character(:), allocatable :: seal, well, rech
  end type unit_out

  type :: out_time
    integer(I4) :: head, rest, srat, wtab, mass, velc, rivr, lakr, sufr, dunr
    integer(I4) :: seal, well, rech
  end type out_time

  type :: out_step
    real(SP) :: head, rest, srat, wtab, mass, velc, rivr, lakr, sufr, dunr
    real(SP) :: seal, well, rech
  end type out_step

  type :: msout_set
    real(DP), allocatable :: sto(:), con(:), sea(:), wel(:), rec(:), sur(:), riv(:)
    real(DP), allocatable :: lak(:), tot(:)
  end type msout_set

  type :: hydr_set
    real(SP), allocatable :: read_vana(:), read_vann(:), read_resi(:)
    real(SP), allocatable :: read_hydx(:), read_hydy(:), read_hydz(:)
    real(SP), allocatable :: read_spst(:), read_pors(:)
    real(DP), allocatable :: read_init(:)
    real(DP), allocatable :: surf_bott(:), surf_reli(:), surf_parm(:)
    real(DP), allocatable :: sat_hydf(:), hydf_surf(:), abyd_surf(:), abyd_well(:)
    real(DP), allocatable :: area_dis(:), conn_dis(:), surf_area(:), rech_area(:)
    real(DP), allocatable :: sea_hydf(:), sea_dis(:), sea_abyd(:)
    real(DP), allocatable :: surf_top(:)
    real(DP), allocatable :: abyd_conn(:), hydf_conn(:), inv_dis(:)
    real(DP), allocatable :: abyd_seal(:), hydf_seal(:), dis_seal(:)
  end type hydr_set

  type :: bcnd_set
    integer(I4) :: rech_num = 0, well_num = 0, prec_num = 0
    integer(I4) :: evap_num = 0, rive_num = 0, lake_num = 0
    integer(I4) :: seal_num
    integer(I4), allocatable :: rech2cals(:), rive2cals(:), lake2cals(:)
    integer(I4), allocatable :: rech_cflag(:), prec_cflag(:), evap_cflag(:)
    integer(I4), allocatable :: well_index(:), well_conn(:)
    integer(I4), allocatable :: sea2cal(:), sea2sea(:)
    integer(I4), allocatable :: seal2calc(:), seal2seal(:)
  end type bcnd_set

  type :: forc_set
    real(DP), allocatable :: read_head(:), abyd_rive(:), abyd_lake(:)
    real(DP), allocatable :: calc_rech(:)
    real(DP), allocatable :: rive_head(:), rive_bott(:), rive_area(:)
    real(DP), allocatable :: lake_head(:), lake_bott(:), lake_area(:)
    real(SP), allocatable :: read_seal(:), read_rech(:), read_well(:)
    real(SP), allocatable :: read_prec(:), read_evap(:)
    real(DP), allocatable :: well_top(:), well_bott(:), calc_well(:)
  end type forc_set

  type :: bflag_set
    integer(I4), allocatable :: wl(:), wd(:), bl(:), de(:), wi(:), le(:), ar(:)
  end type bflag_set

  type :: bcalc_set
    real(SP), allocatable :: wl(:), wd(:), bl(:), de(:), wi(:), le(:), ar(:)
  end type bcalc_set

  type :: bcount_set
    integer(I4) :: wl = 0, wd = 0, bl = 0, de = 0, wi = 0, le = 0, ar = 0
  end type bcount_set

  type :: rlbc_set
    type(bflag_set) :: cflag
    type(bcalc_set) :: calc
    type(bcount_set) :: num
  end type rlbc_set

  type :: bfview_set
    integer(I4) :: wl = 0, wd = 0, bl = 0, de = 0, wi = 0, le = 0, ar = 0
  end type bfview_set

  type :: bound_fview
    integer(I4) :: seal = 0, rech = 0, well = 0, prec = 0, evap = 0
  end type bound_fview

end module types_module
