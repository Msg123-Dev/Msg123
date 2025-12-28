module mpi_set
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SZERO, DZERO, FACE
  use utility_module, only: write_err_stop
  use initial_module, only: my_rank, st_sim, st_grid, st_in_type, st_out_type,&
                            precon_type
  use mpi_initfin, only: my_comm
  use mpi

  implicit none
  private
  public :: bcast_ici_set, bcast_sim_flag, bcast_xyz_num, bcast_clas_set, bcast_glob_flag
  public :: set_calc_view, set_seal_view, set_cell_view, set_rest_view, set_write_fview
  public :: senrec_reg_info, senrec_grid_num
  public :: bcast_calc_ftype, bcast_sim_val, bcast_glob_xyzv
  public :: bcast_retn_clas, bcast_parm_clas, bcast_init_dep
  public :: bcast_bound_ftype, bcast_out_type, bcast_solval
  public :: senrec_neibval, senrec_faceval
  public :: scatter_xyval, scatter_xyzval
  public :: bcast_clas_val, bcast_2dpoint, bcast_3dpoint, bcast_wellpoint
  integer(I4), public :: cals_i4view, calc_i4view
  integer(I4), public :: cals_r4view, cals_r4hview, calc_r4view, calc_r4hview
  integer(I4), public :: surf_r4view, surf_r4hview, cell_r4view, cell_r4hview, rest_view
  integer(I4), public :: write_2dview, write_3dview, write_no2dview, write_no3dview

  interface scatter_xyval
    module procedure scatter_i4xy
    module procedure scatter_r4xy
    module procedure scatter_r8xy
  end interface

  interface scatter_xyzval
    module procedure scatter_i4xyz
    module procedure scatter_r4xyz
    module procedure scatter_r8xyz
  end interface

  interface senrec_neibval
    module procedure senrec_neib_i4
    module procedure senrec_neib_r4
    module procedure senrec_neib_r8
  end interface

  interface senrec_faceval
    module procedure senrec_face_i4
    module procedure senrec_face_r4
    module procedure senrec_face_r8
  end interface

  ! -- local

  contains

  subroutine bcast_ici_set(mcomp, mgrid, ici_file, get_mat, put_mat, get_cama, put_cama)
  !***************************************************************************************
  ! bcast_ici_set -- Bcast ici set
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(inout) :: mcomp, mgrid, ici_file
    logical, intent(inout) :: get_mat, put_mat, get_cama, put_cama
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: length
    !-------------------------------------------------------------------------------------
    ierr = 0
    length = len_trim(mcomp)
    call MPI_BCAST(mcomp, length, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the component name for ILS.")
      end if
    end if

    length = len_trim(mgrid)
    call MPI_BCAST(mgrid, length, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the grid name for ILS.")
      end if
    end if

    length = len_trim(ici_file)
    call MPI_BCAST(ici_file, length, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the coupling conf file for ILS.")
      end if
    end if

    call MPI_BCAST(get_mat, 1, MPI_LOGICAL, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the MATSIRO get flag for ILS.")
      end if
    end if

    call MPI_BCAST(put_mat, 1, MPI_LOGICAL, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the MATSIRO put flag for ILS.")
      end if
    end if

    call MPI_BCAST(get_cama, 1, MPI_LOGICAL, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the CaMa-Flood get flag for ILS.")
      end if
    end if

    call MPI_BCAST(put_cama, 1, MPI_LOGICAL, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the CaMa-Flood put flag for ILS.")
      end if
    end if

  end subroutine bcast_ici_set

  subroutine bcast_sim_flag()
  !***************************************************************************************
  ! bcast_sim_flag -- Bcast simulation flag
  !***************************************************************************************
    ! -- module
    use initial_module, only: noclas_flag
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(st_sim%sim_type, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the simulation type.")
      end if
    end if

    call MPI_BCAST(st_sim%res_type, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the calculation time type.")
      end if
    end if

    call MPI_BCAST(st_sim%reg_type, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the calculation region type.")
      end if
    end if

    call MPI_BCAST(st_sim%reg_neib, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the calculation neighbor region.")
      end if
    end if

    call MPI_BCAST(noclas_flag, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the classification flag.")
      end if
    end if

    call MPI_BCAST(st_out_type%calg, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the calculation grid flag for output.")
      end if
    end if

    call MPI_BCAST(precon_type, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the preconditoner type.")
      end if
    end if

  end subroutine bcast_sim_flag

  subroutine bcast_xyz_num()
  !***************************************************************************************
  ! bcast_xyz_num -- Bcast xyz number
  !***************************************************************************************
    ! -- module

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(st_grid%nx, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the grid size in x direcition.")
      end if
    end if

    call MPI_BCAST(st_grid%ny, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the grid size in y direcition.")
      end if
    end if

    call MPI_BCAST(st_grid%nz, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the grid size in z direcition.")
      end if
    end if

    call MPI_BCAST(st_grid%nxyz, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the grid size in all direcition.")
      end if
    end if

  end subroutine bcast_xyz_num

  subroutine bcast_clas_set()
  !***************************************************************************************
  ! bcast_clas_set -- Bcast classification setting
  !***************************************************************************************
    ! -- module
    use constval_module, only: VARLEN
    use initial_module, only: st_clas
    ! -- inout

    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: char_leng
    integer(I4) :: clas_totn, max_clas, clas_len
    integer(I4), allocatable :: clas_num(:)
    integer(I4), allocatable :: clas_i(:,:), clas_j(:,:), clas_k(:,:)
    character(VARLEN), allocatable :: clas_name(:)
    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      clas_totn = st_clas%totn
    end if

    ierr = 0
    call MPI_BCAST(clas_totn, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the classification size.")
      end if
    end if

    char_leng = len_trim(st_sim%inact_name)
    call MPI_BCAST(st_sim%inact_name, char_leng, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the inactive region name.")
      end if
    end if

    if (my_rank /= 0) then
      st_clas%totn = clas_totn
    end if

    allocate(clas_name(clas_totn))
    !$omp parallel
    !$omp workshare
    clas_name(:) = ""
    !$omp end workshare

    if (my_rank == 0) then
      !$omp do private(i)
      do i = 1, clas_totn
        clas_name(i) = st_clas%name(i)
      end do
      !$omp end do
    end if
    !$omp end parallel

    do i = 1, clas_totn
      char_leng = len_trim(clas_name(i))
      call MPI_BCAST(clas_name(i), char_leng, MPI_CHARACTER, 0, my_comm, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (my_rank == 0) then
          call write_err_stop("Broadcast the classification name "//clas_name(i)//".")
        end if
      end if
    end do

    if (my_rank /= 0) then
      allocate(st_clas%name(clas_totn))
      !$omp parallel workshare
      st_clas%name(:) = "" ; st_clas%name(:) = clas_name(:)
      !$omp end parallel workshare
    end if

    allocate(clas_num(clas_totn))
    !$omp parallel workshare
    clas_num(:) = 0
    !$omp end parallel workshare

    if (my_rank == 0) then
      !$omp parallel do private(i)
      do i = 1, clas_totn
        clas_num(i) = st_clas%num(i)
      end do
      !$omp end parallel do
    end if

    call MPI_BCAST(clas_num, clas_totn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the size in each classification name.")
      end if
    end if

    if (my_rank /= 0) then
      allocate(st_clas%num(clas_totn))
      !$omp parallel workshare
      st_clas%num(:) = 0 ; st_clas%num(:) = clas_num(:)
      !$omp end parallel workshare
    end if

    max_clas = 0
    !$omp parallel do private(i) reduction(max:max_clas)
    do i = 1, clas_totn
      if (clas_num(i) > max_clas) then
        max_clas = clas_num(i)
      end if
    end do
    !$omp end parallel do

    allocate(clas_i(max_clas,clas_totn), clas_j(max_clas,clas_totn))
    allocate(clas_k(max_clas,clas_totn))
    !$omp parallel workshare
    clas_i(:,:) = 0 ; clas_j(:,:) = 0 ; clas_k(:,:) = 0
    !$omp end parallel workshare

    if (my_rank == 0) then
      !$omp parallel workshare
      clas_i(:,:) = st_clas%i(:,:) ; clas_j(:,:) = st_clas%j(:,:)
      clas_k(:,:) = st_clas%k(:,:)
      !$omp end parallel workshare
    end if

    clas_len = max_clas*clas_totn

    call MPI_BCAST(clas_i, clas_len, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the x number in classification file.")
      end if
    end if

    call MPI_BCAST(clas_j, clas_len, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the y number in classification file.")
      end if
    end if

    call MPI_BCAST(clas_k, clas_len, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the z number in classification file.")
      end if
    end if

    if (my_rank /= 0) then
      allocate(st_clas%i(max_clas,clas_totn), st_clas%j(max_clas,clas_totn))
      allocate(st_clas%k(max_clas,clas_totn))
      !$omp parallel workshare
      st_clas%i(:,:) = 0 ; st_clas%j(:,:) = 0 ; st_clas%k(:,:) = 0
      st_clas%i(:,:) = clas_i(:,:) ; st_clas%j(:,:) = clas_j(:,:)
      st_clas%k(:,:) = clas_k(:,:)
      !$omp end parallel workshare
    end if

    deallocate(clas_name, clas_num, clas_i, clas_j, clas_k)

  end subroutine bcast_clas_set

  subroutine bcast_glob_flag(totreg_num, reg_flag, mpi_flag)
  !***************************************************************************************
  ! bcast_glob_flag -- Bcast global flag value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: totreg_num
    integer(I4), intent(inout) :: reg_flag(:), mpi_flag(:)
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(totreg_num, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast the total region size.")
      end if
    end if

    call MPI_BCAST(reg_flag(1), st_grid%nxyz, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast global region flag.")
      end if
    end if

    call MPI_BCAST(mpi_flag(1), st_grid%nxyz, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast global mpi flag.")
      end if
    end if

  end subroutine bcast_glob_flag

  subroutine set_calc_view(loc_ncals, loc_ncalc, l2g_ij, l2g_ijk)
  !***************************************************************************************
  ! set_calc_view -- Set calculation view
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_ncals, loc_ncalc
    integer(I4), intent(in) :: l2g_ij(:), l2g_ijk(:)
    ! -- local
    integer(I4) :: i, ij, ijk, ierr
    integer(I4), allocatable :: xyblock(:), xyzblock(:), xydis(:), xyzdis(:)
    integer(I4), allocatable :: xytype(:), xyztype(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(xyblock(loc_ncals+2), xydis(loc_ncals+2), xytype(loc_ncals+2))
    !$omp parallel
    !$omp workshare
    xyblock(:) = 1 ; xydis(:) = 0 ; xytype(:) = MPI_INTEGER
    !$omp end workshare

    xytype(1) = MPI_LB ; xytype(loc_ncals+2) = MPI_UB
    xydis(loc_ncals+2) = st_grid%nx*st_grid%ny*I4

    !$omp do private(i, ij)
    do i = 1, loc_ncals
      ij = i + 1 ; xydis(ij) = (l2g_ij(i)-1)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncals+2, xyblock, xydis, xytype, cals_i4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype without header for within range in xy direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(cals_i4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype without header for within range in xy direction.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)

    allocate(xyblock(loc_ncals+2), xydis(loc_ncals+2), xytype(loc_ncals+2))
    !$omp parallel
    !$omp workshare
    xyblock(:) = 1 ; xydis(:) = 0 ; xytype(:) = MPI_REAL4
    !$omp end workshare

    xytype(1) = MPI_LB ; xytype(loc_ncals+2) = MPI_UB
    xydis(loc_ncals+2) = st_grid%nx*st_grid%ny*I4

    !$omp do private(i, ij)
    do i = 1, loc_ncals
      ij = i + 1 ; xydis(ij) = (l2g_ij(i)-1)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncals+2, xyblock, xydis, xytype, cals_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype without header for within range in xy direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(cals_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype without header for within range in xy direction.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)

    allocate(xyblock(loc_ncals+3), xydis(loc_ncals+3), xytype(loc_ncals+3))
    !$omp parallel
    !$omp workshare
    xyblock(:) = 1 ; xydis(:) = 0 ; xytype(:) = MPI_REAL4
    !$omp end workshare

    xytype(1) = MPI_LB ; xytype(loc_ncals+3) = MPI_UB
    xydis(loc_ncals+3) = (st_grid%nx*st_grid%ny+1)*I4

    !$omp do private(i, ij)
    do i = 1, loc_ncals
      ij = i + 2 ; xydis(ij) = l2g_ij(i)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncals+3, xyblock, xydis, xytype, cals_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype with header for within range in xy direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(cals_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype with header for within range in xy direction.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)

    allocate(xyzblock(loc_ncalc+2), xyzdis(loc_ncalc+2), xyztype(loc_ncalc+2))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_INTEGER
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_ncalc+2) = MPI_UB
    xyzdis(loc_ncalc+2) = st_grid%nxyz*I4

    !$omp do private(i, ijk)
    do i = 1, loc_ncalc
      ijk = i + 1 ; xyzdis(ijk) = (l2g_ijk(i)-1)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncalc+2, xyzblock, xyzdis, xyztype, calc_i4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype without header for in range in xyz direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(calc_i4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype without header for in range in xyz direction.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

    allocate(xyzblock(loc_ncalc+2), xyzdis(loc_ncalc+2), xyztype(loc_ncalc+2))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_REAL4
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_ncalc+2) = MPI_UB
    xyzdis(loc_ncalc+2) = st_grid%nxyz*I4

    !$omp do private(i, ijk)
    do i = 1, loc_ncalc
      ijk = i + 1 ; xyzdis(ijk) = (l2g_ijk(i)-1)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncalc+2, xyzblock, xyzdis, xyztype, calc_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype without header for in range in xyz direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(calc_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype without header for in range in xyz direction.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

    allocate(xyzblock(loc_ncalc+3), xyzdis(loc_ncalc+3), xyztype(loc_ncalc+3))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_REAL4
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_ncalc+3) = MPI_UB
    xyzdis(loc_ncalc+3) = (st_grid%nxyz+1)*I4

    !$omp do private(i, ijk)
    do i = 1, loc_ncalc
      ijk = i + 2 ; xyzdis(ijk) = l2g_ijk(i)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncalc+3, xyzblock, xyzdis, xyztype, calc_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype with header for within range in xyz direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(calc_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype with header for within range in xyz direction.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

  end subroutine set_calc_view

  subroutine set_seal_view(loc_seas, loc_seac, l2g_ij, l2g_ijk)
  !***************************************************************************************
  ! set_seal_view -- Set seal view
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_seas, loc_seac
    integer(I4), intent(in) :: l2g_ij(:), l2g_ijk(:)
    ! -- local
    integer(I4) :: i, ij, ijk, ierr
    integer(I4), allocatable :: xyblock(:), xyzblock(:), xydis(:), xyzdis(:)
    integer(I4), allocatable :: xytype(:), xyztype(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(xyblock(loc_seas+2), xydis(loc_seas+2), xytype(loc_seas+2))
    !$omp parallel
    !$omp workshare
    xyblock(:) = 1 ; xydis(:) = 0 ; xytype(:) = MPI_REAL4
    !$omp end workshare

    xytype(1) = MPI_LB ; xytype(loc_seas+2) = MPI_UB
    xydis(loc_seas+2) = st_grid%nx*st_grid%ny*I4

    !$omp do private(i, ij)
    do i = 1, loc_seas
      ij = i + 1 ; xydis(ij) = (l2g_ij(i)-1)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_seas+2, xyblock, xydis, xytype, surf_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype without header for out of range in xy direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(surf_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype without header for out of range in xy direction.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)

    allocate(xyblock(loc_seas+3), xydis(loc_seas+3), xytype(loc_seas+3))
    !$omp parallel
    !$omp workshare
    xyblock(:) = 1 ; xydis(:) = 0 ; xytype(:) = MPI_REAL4
    !$omp end workshare

    xytype(1) = MPI_LB ; xytype(loc_seas+3) = MPI_UB
    xydis(loc_seas+3) = (st_grid%nx*st_grid%ny+1)*I4

    !$omp do private(i, ij)
    do i = 1, loc_seas
      ij = i + 2 ; xydis(ij) = l2g_ij(i)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_seas+3, xyblock, xydis, xytype, surf_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype with header for out of range in xy direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(surf_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype with header for out of range in xy direction.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)


    allocate(xyzblock(loc_seac+2), xyzdis(loc_seac+2), xyztype(loc_seac+2))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_REAL4
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_seac+2) = MPI_UB
    xyzdis(loc_seac+2) = st_grid%nxyz*I4

    !$omp do private(i, ijk)
    do i = 1, loc_seac
      ijk = i + 1 ; xyzdis(ijk) = (l2g_ijk(i)-1)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_seac+2, xyzblock, xyzdis, xyztype, cell_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype without header for out of range in xyz direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(cell_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype without header for out of range in xyz direction.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

    allocate(xyzblock(loc_seac+3), xyzdis(loc_seac+3), xyztype(loc_seac+3))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_REAL4
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_seac+3) = MPI_UB
    xyzdis(loc_seac+3) = (st_grid%nxyz+1)*I4

    !$omp do private(i, ijk)
    do i = 1, loc_seac
      ijk = i + 2 ; xyzdis(ijk) = l2g_ijk(i)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_seac+3, xyzblock, xyzdis, xyztype, cell_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype with header for out of range in xyz direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(cell_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype with header for out of range in xyz direction.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

  end subroutine set_seal_view

  subroutine set_cell_view(loc_nsurf, loc_ncell, l2g_ij, l2g_ijk)
  !***************************************************************************************
  ! set_cell_view -- Set cell view
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_nsurf, loc_ncell
    integer(I4), intent(in) :: l2g_ij(:), l2g_ijk(:)
    ! -- local
    integer(I4) :: i, ij, ijk, ierr
    integer(I4), allocatable :: xyblock(:), xyzblock(:), xydis(:), xyzdis(:)
    integer(I4), allocatable :: xytype(:), xyztype(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(xyblock(loc_nsurf+2), xydis(loc_nsurf+2), xytype(loc_nsurf+2))
    !$omp parallel
    !$omp workshare
    xyblock(:) = 1 ; xydis(:) = 0 ; xytype(:) = MPI_REAL4
    !$omp end workshare

    xytype(1) = MPI_LB ; xytype(loc_nsurf+2) = MPI_UB
    xydis(loc_nsurf+2) = st_grid%nx*st_grid%ny*I4

    !$omp do private(i, ij)
    do i = 1, loc_nsurf
      ij = i + 1 ; xydis(ij) = (l2g_ij(i)-1)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_nsurf+2, xyblock, xydis, xytype, surf_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype without header for within range and sea in xy direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(surf_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype without header for within range and sea in xy direction.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)

    allocate(xyblock(loc_nsurf+3), xydis(loc_nsurf+3), xytype(loc_nsurf+3))
    !$omp parallel
    !$omp workshare
    xyblock(:) = 1 ; xydis(:) = 0 ; xytype(:) = MPI_REAL4
    !$omp end workshare

    xytype(1) = MPI_LB ; xytype(loc_nsurf+3) = MPI_UB
    xydis(loc_nsurf+3) = (st_grid%nx*st_grid%ny+1)*I4

    !$omp do private(i, ij)
    do i = 1, loc_nsurf
      ij = i + 2 ; xydis(ij) = l2g_ij(i)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_nsurf+3, xyblock, xydis, xytype, surf_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype with header for within range and sea in xy direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(surf_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype with header for within range and sea in xy direction.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)

    allocate(xyzblock(loc_ncell+2), xyzdis(loc_ncell+2), xyztype(loc_ncell+2))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_REAL4
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_ncell+2) = MPI_UB
    xyzdis(loc_ncell+2) = st_grid%nxyz*I4

    !$omp do private(i, ijk)
    do i = 1, loc_ncell
      ijk = i + 1 ; xyzdis(ijk) = (l2g_ijk(i)-1)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncell+2, xyzblock, xyzdis, xyztype, cell_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype without header for within range and sea in xyz direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(cell_r4view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype without header for within range and sea in xyz direction.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

    allocate(xyzblock(loc_ncell+3), xyzdis(loc_ncell+3), xyztype(loc_ncell+3))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_REAL4
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_ncell+3) = MPI_UB
    xyzdis(loc_ncell+3) = (st_grid%nxyz+1)*I4

    !$omp do private(i, ijk)
    do i = 1, loc_ncell
      ijk = i + 2 ; xyzdis(ijk) = l2g_ijk(i)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncell+3, xyzblock, xyzdis, xyztype, cell_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype with header for within range and sea in xyz direction.")
      end if
    end if

    call MPI_TYPE_COMMIT(cell_r4hview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype with header for within range and sea in xyz direction.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

  end subroutine set_cell_view

  subroutine set_rest_view(loc_ncalc, glob_ncalc, l2g_ijk)
  !***************************************************************************************
  ! set_rest_view -- Set restart view
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_ncalc, glob_ncalc
    integer(I4), intent(in) :: l2g_ijk(:)
    ! -- local
    integer(I4) :: i, ijk, ierr
    integer(I4), allocatable :: xyzblock(:), xyzdis(:), xyztype(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(xyzblock(loc_ncalc+3), xyzdis(loc_ncalc+3), xyztype(loc_ncalc+3))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_REAL8
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_ncalc+3) = MPI_UB
    xyzdis(loc_ncalc+3) = (glob_ncalc+1)*8

    !$omp do private(i, ijk)
    do i = 1, loc_ncalc
      ijk = i + 2 ; xyzdis(ijk) = l2g_ijk(i)*8
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_ncalc+3, xyzblock, xyzdis, xyztype, rest_view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype for restart.")
      end if
    end if

    call MPI_TYPE_COMMIT(rest_view, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype for restart.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

  end subroutine set_rest_view

  subroutine set_write_fview(cals, calc, l2g_ij, l2g_ijk, l2g_nos, l2g_noc)
  !***************************************************************************************
  ! set_write_fview -- Set write file view
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: cals, calc
    integer(I4), intent(in) :: l2g_ij(:), l2g_ijk(:), l2g_nos(:), l2g_noc(:)
    ! -- local
    integer(I4) :: i, ij, ijk, ierr
    integer(I4) :: no_cals, no_calc, loc_cals, loc_calc
    integer(I4), allocatable :: xyblock(:), xydis(:), xytype(:)
    integer(I4), allocatable :: xyzblock(:), xyzdis(:), xyztype(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    no_cals = size(l2g_nos(:)) ; no_calc = size(l2g_noc(:))

    loc_cals = cals + no_cals + 3
    allocate(xyblock(loc_cals), xydis(loc_cals), xytype(loc_cals))
    !$omp parallel
    !$omp workshare
    xyblock(:) = 1 ; xydis(:) = 0 ; xytype(:) = MPI_REAL4
    !$omp end workshare

    xytype(1) = MPI_LB ; xytype(loc_cals) = MPI_UB
    xydis(loc_cals) = (st_grid%nx*st_grid%ny+1)*I4

    !$omp do private(i, ij)
    do i = 1, cals
      ij = i + 2 ; xydis(ij) = l2g_ij(i)*I4
    end do
    !$omp end do
    !$omp do private(i, ij)
    do i = 1, no_cals
      ij = i + cals + 2 ; xydis(ij) = l2g_nos(i)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_cals, xyblock, xydis, xytype, write_2dview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype for 2d output.")
      end if
    end if

    call MPI_TYPE_COMMIT(write_2dview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype for 2d output.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)

    loc_calc = calc + no_calc + 3
    allocate(xyzblock(loc_calc), xyzdis(loc_calc), xyztype(loc_calc))
    !$omp parallel
    !$omp workshare
    xyzblock(:) = 1 ; xyzdis(:) = 0 ; xyztype(:) = MPI_REAL4
    !$omp end workshare

    xyztype(1) = MPI_LB ; xyztype(loc_calc) = MPI_UB
    xyzdis(loc_calc) = (st_grid%nxyz+1)*I4

    !$omp do private(i, ijk)
    do i = 1, calc
      ijk = i + 2 ; xyzdis(ijk) = l2g_ijk(i)*I4
    end do
    !$omp end do
    !$omp do private(i, ijk)
    do i = 1, no_calc
      ijk = i + calc + 2 ; xyzdis(ijk) = l2g_noc(i)*I4
    end do
    !$omp end do
    !$omp end parallel

    call MPI_TYPE_STRUCT(loc_calc, xyzblock, xyzdis, xyztype, write_3dview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Create struct datatype for 3d output.")
      end if
    end if

    call MPI_TYPE_COMMIT(write_3dview, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Commit struct datatype for 3d output.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

  end subroutine set_write_fview

  subroutine senrec_reg_info(pron, s_nreg, r_nreg)
  !***************************************************************************************
  ! senrec_reg_info -- Send and Receive region information
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: pron
    integer(I4), intent(in) :: s_nreg(:)
    integer(I4), allocatable, intent(out) :: r_nreg(:)
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: nreg
    integer(I4), allocatable :: istat(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(istat(MPI_STATUS_SIZE))
    !$omp parallel workshare
    istat(:) = 0
    !$omp end parallel workshare

    if (pron == my_rank+1) then
      call MPI_SEND(size(s_nreg), 1, MPI_INTEGER, 0, 0, my_comm, ierr)
    else if (my_rank == 0) then
      call MPI_RECV(nreg, 1, MPI_INTEGER, pron-1, 0, my_comm, istat, ierr)
    end if

    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send&Receive local region size.")
      end if
    end if

    if (pron == my_rank+1) then
      call MPI_SEND(s_nreg(1), size(s_nreg), MPI_INTEGER, 0, 0, my_comm, ierr)
    else if (my_rank == 0) then
      allocate(r_nreg(nreg))
      !$omp parallel workshare
      r_nreg(:) = 0
      !$omp end parallel workshare
      call MPI_RECV(r_nreg(1), nreg, MPI_INTEGER, pron-1, 0, my_comm, istat, ierr)
    end if

    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send&Receive local region flag.")
      end if
    end if

    deallocate(istat)

  end subroutine senrec_reg_info

  subroutine senrec_grid_num(pron, x_num, y_num, z_num)
  !***************************************************************************************
  ! senrec_grid_num -- Send and Receive grid number
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: pron
    integer(I4), intent(inout) :: x_num, y_num, z_num
    ! -- local
    integer(I4) :: ierr
    integer(I4), allocatable :: istat(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(istat(MPI_STATUS_SIZE))
    !$omp parallel workshare
    istat(:) = 0
    !$omp end parallel workshare

    if (pron == my_rank+1) then
      call MPI_SEND(x_num, 1, MPI_INTEGER, 0, 1, my_comm, ierr)
      call MPI_SEND(y_num, 1, MPI_INTEGER, 0, 2, my_comm, ierr)
      call MPI_SEND(z_num, 1, MPI_INTEGER, 0, 3, my_comm, ierr)
    else if (my_rank == 0) then
      call MPI_RECV(x_num, 1, MPI_INTEGER, pron-1, 1, my_comm, istat, ierr)
      call MPI_RECV(y_num, 1, MPI_INTEGER, pron-1, 2, my_comm, istat, ierr)
      call MPI_RECV(z_num, 1, MPI_INTEGER, pron-1, 3, my_comm, istat, ierr)
    end if

    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send&Receive local grid number.")
      end if
    end if

    deallocate(istat)

  end subroutine senrec_grid_num

  subroutine bcast_calc_ftype()
  !***************************************************************************************
  ! bcast_calc_ftype -- Bcast calculation file type
  !***************************************************************************************
    ! -- module

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(st_in_type%retn, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast retention file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%parm, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast parameter file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%geog, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast geography file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%init, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast initial file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%wtab, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast water table file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%mass, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast massbalance file type.")
      end if
    end if

  end subroutine bcast_calc_ftype

  subroutine bcast_sim_val()
  !***************************************************************************************
  ! bcast_sim_val -- Bcast simulation value
  !***************************************************************************************
    ! -- module
    use initial_module, only: nlevel
    use read_input, only: len_scal, len_scal_inv
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    integer(I4) :: char_leng
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(st_sim%sta_date, 6, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast simulation date.")
      end if
    end if

    call MPI_BCAST(nlevel, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast multigrid level.")
      end if
    end if

    call MPI_BCAST(st_sim%end_time, 1, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast simulation end time.")
      end if
    end if

    call MPI_BCAST(st_sim%cal_fact, 1, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast calculation factor for end time.")
      end if
    end if

    call MPI_BCAST(len_scal, 1, MPI_REAL8, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast length scale.")
      end if
    end if

    call MPI_BCAST(len_scal_inv, 1, MPI_REAL8, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast inverse length scale.")
      end if
    end if

    char_leng = len_trim(st_sim%cal_unit)
    call MPI_BCAST(st_sim%cal_unit, char_leng, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast calculation unit.")
      end if
    end if

  end subroutine bcast_sim_val

  subroutine bcast_glob_xyzv()
  !***************************************************************************************
  ! bcast_glob_xyzv -- Bcast global xyz value
  !***************************************************************************************
    ! -- module
    use read_input, only: glob_x, glob_y, glob_z
    ! -- inout

    ! -- local
    integer(I4) :: ierr, nxy, nxyz
    !-------------------------------------------------------------------------------------
    ierr = 0
    nxy = (st_grid%nx+1)*(st_grid%ny+1) ; nxyz = nxy*(st_grid%nz+1)

    call MPI_BCAST(glob_x(1,1), nxy, MPI_REAL8, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast global grid value in x direcition.")
      end if
    end if

    call MPI_BCAST(glob_y(1,1), nxy, MPI_REAL8, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast global grid value in y direcition.")
      end if
    end if

    call MPI_BCAST(glob_z(1,1,1), nxyz, MPI_REAL8, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast global grid value in z direcition.")
      end if
    end if

  end subroutine bcast_glob_xyzv

  subroutine bcast_bound_ftype()
  !***************************************************************************************
  ! bcast_bound_ftype -- Bcast boundary file type
  !***************************************************************************************
    ! -- module

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(st_in_type%seal, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast input sea level file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%rech, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast input recharge file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%well, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast input well file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%prec, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast input precipitation file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%evap, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast input evapotranspiration file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%rive, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast input river file type.")
      end if
    end if

    call MPI_BCAST(st_in_type%lake, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast input lake file type.")
      end if
    end if

  end subroutine bcast_bound_ftype

  subroutine bcast_out_type()
  !***************************************************************************************
  ! bcast_out_type -- Bcast output file number
  !***************************************************************************************
    ! -- module

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(st_out_type%wtab, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output water table file type.")
      end if
    end if

    call MPI_BCAST(st_out_type%mass, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output massbalance file type.")
      end if
    end if

    call MPI_BCAST(st_out_type%velc, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output velocity file type.")
      end if
    end if

    call MPI_BCAST(st_out_type%rivr, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output river runoff file type.")
      end if
    end if

    call MPI_BCAST(st_out_type%lakr, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output lake runoff file type.")
      end if
    end if

    call MPI_BCAST(st_out_type%sufr, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output surface runoff file type.")
      end if
    end if

    call MPI_BCAST(st_out_type%dunr, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output dunne overland runoff file type.")
      end if
    end if

    call MPI_BCAST(st_out_type%well, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output well results file type.")
      end if
    end if

    call MPI_BCAST(st_out_type%rech, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast output recharge results file type.")
      end if
    end if

  end subroutine bcast_out_type

  subroutine senrec_neib_i4(nbtot, nbnum, sind, rind, sitem, ritem, inv, outv)
  !***************************************************************************************
  ! senrec_neib_i4 -- Send and Receive neighbor integer value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: nbtot
    integer(I4), intent(in) :: nbnum(:), sind(0:), rind(0:), sitem(:), ritem(:), inv(:)
    integer(I4), intent(out) :: outv(:)
    ! -- local
    integer(I4) :: i, j, k, jj, kk, ierr
    integer(I4) :: sbuflen, rbuflen, send_count, recv_count
    integer(I4) :: is_sta, is_end, ir_sta, ir_end
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_send(:,:), stat_recv(:,:)
    integer(I4), allocatable :: sbufint(:), rbufint(:)
    !-------------------------------------------------------------------------------------
    allocate(requ_send(nbtot), requ_recv(nbtot))
    allocate(stat_send(MPI_STATUS_SIZE,nbtot), stat_recv(MPI_STATUS_SIZE,nbtot))
    !$omp parallel
    !$omp workshare
    requ_send(:) = 0 ; requ_recv(:) = 0
    stat_send(:,:) = 0 ; stat_recv(:,:) = 0
    !$omp end workshare

    send_count = 0 ; recv_count = 0
    !$omp do private(i, j, k) reduction(+:send_count, recv_count)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        send_count = send_count + 1
      end do
      do k = rind(i-1)+1, rind(i)
        recv_count = recv_count + 1
      end do
    end do
    !$omp end do
    !$omp end parallel
    allocate(sbufint(send_count), rbufint(recv_count))
    !$omp parallel
    !$omp workshare
    sbufint(:) = 0 ; rbufint(:) = 0
    !$omp end workshare

    !$omp do private(i, j, jj)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        jj = sitem(j)
        sbufint(j) = inv(jj)
      end do
    end do
    !$omp end do

    ierr = 0
    !$omp do private(i, is_sta, is_end, sbuflen)
    do i = 1, nbtot
      is_sta = sind(i-1)+1 ; is_end = sind(i)
      sbuflen = is_end - is_sta + 1
      if (sbuflen /= 0) then
        call MPI_ISEND(sbufint(is_sta), sbuflen, MPI_REAL8, nbnum(i), 0, my_comm,&
                       requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, ir_sta, ir_end, rbuflen)
    do i = 1, nbtot
      ir_sta = rind(i-1)+1 ; ir_end = rind(i)
      rbuflen = ir_end - ir_sta + 1
      if (rbuflen /= 0) then
        call MPI_IRECV(rbufint(ir_sta), rbuflen, MPI_REAL8, nbnum(i), 0, my_comm,&
                       requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(nbtot, requ_recv, stat_recv, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Receive neighbor integer values.")
      end if
    end if

    call MPI_WAITALL(nbtot, requ_send, stat_send, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send neighbor integer values.")
      end if
    end if

    !$omp do private(i, k, kk)
    do i = 1, nbtot
      do k = rind(i-1)+1, rind(i)
        kk = ritem(k)
        outv(kk) = rbufint(k)
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_send, stat_recv, sbufint, rbufint)

  end subroutine senrec_neib_i4

  subroutine senrec_neib_r4(nbtot, nbnum, sind, rind, sitem, ritem, inv, outv)
  !***************************************************************************************
  ! senrec_neib_r4 -- Send and Receive neighbor real4 value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: nbtot
    integer(I4), intent(in) :: nbnum(:), sind(0:), rind(0:), sitem(:), ritem(:)
    real(SP), intent(in) :: inv(:)
    real(SP), intent(out) :: outv(:)
    ! -- local
    integer(I4) :: i, j, k, jj, kk, ierr
    integer(I4) :: sbuflen, rbuflen, send_count, recv_count
    integer(I4) :: is_sta, is_end, ir_sta, ir_end
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_send(:,:), stat_recv(:,:)
    real(SP), allocatable :: sbufreal(:), rbufreal(:)
    !-------------------------------------------------------------------------------------
    allocate(requ_send(nbtot), requ_recv(nbtot))
    allocate(stat_send(MPI_STATUS_SIZE,nbtot), stat_recv(MPI_STATUS_SIZE,nbtot))
    !$omp parallel
    !$omp workshare
    requ_send(:) = 0 ; requ_recv(:) = 0
    stat_send(:,:) = 0 ; stat_recv(:,:) = 0
    !$omp end workshare

    send_count = 0 ; recv_count = 0
    !$omp do private(i, j, k) reduction(+:send_count, recv_count)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        send_count = send_count + 1
      end do
      do k = rind(i-1)+1, rind(i)
        recv_count = recv_count + 1
      end do
    end do
    !$omp end do
    !$omp end parallel
    allocate(sbufreal(send_count), rbufreal(recv_count))
    !$omp parallel
    !$omp workshare
    sbufreal(:) = SZERO ; rbufreal(:) = SZERO
    !$omp end workshare

    !$omp do private(i, j, jj)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        jj = sitem(j)
        sbufreal(j) = inv(jj)
      end do
    end do
    !$omp end do

    ierr = 0
    !$omp do private(i, is_sta, is_end, sbuflen)
    do i = 1, nbtot
      is_sta = sind(i-1)+1 ; is_end = sind(i)
      sbuflen = is_end - is_sta + 1
      if (sbuflen /= 0) then
        call MPI_ISEND(sbufreal(is_sta), sbuflen, MPI_REAL4, nbnum(i), 0, my_comm,&
                       requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, ir_sta, ir_end, rbuflen)
    do i = 1, nbtot
      ir_sta = rind(i-1)+1 ; ir_end = rind(i)
      rbuflen = ir_end - ir_sta + 1
      if (rbuflen /= 0) then
        call MPI_IRECV(rbufreal(ir_sta), rbuflen, MPI_REAL4, nbnum(i), 0, my_comm,&
                       requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(nbtot, requ_recv, stat_recv, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Receive neighbor real4 values.")
      end if
    end if

    call MPI_WAITALL(nbtot, requ_send, stat_send, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send neighbor real4 values.")
      end if
    end if

    !$omp do private(i, k, kk)
    do i = 1, nbtot
      do k = rind(i-1)+1, rind(i)
        kk = ritem(k)
        outv(kk) = rbufreal(k)
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_send, stat_recv, sbufreal, rbufreal)

  end subroutine senrec_neib_r4

  subroutine senrec_neib_r8(nbtot, nbnum, sind, rind, sitem, ritem, inv, outv)
  !***************************************************************************************
  ! senrec_neib_r8 -- Send and Receive neighbor real8 value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: nbtot
    integer(I4), intent(in) :: nbnum(:), sind(0:), rind(0:), sitem(:), ritem(:)
    real(DP), intent(in) :: inv(:)
    real(DP), intent(out) :: outv(:)
    ! -- local
    integer(I4) :: i, j, k, jj, kk, ierr
    integer(I4) :: sbuflen, rbuflen, send_count, recv_count
    integer(I4) :: is_sta, is_end, ir_sta, ir_end
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_send(:,:), stat_recv(:,:)
    real(DP), allocatable :: sbufreal(:), rbufreal(:)
    !-------------------------------------------------------------------------------------
    allocate(requ_send(nbtot), requ_recv(nbtot))
    allocate(stat_send(MPI_STATUS_SIZE,nbtot), stat_recv(MPI_STATUS_SIZE,nbtot))
    !$omp parallel
    !$omp workshare
    requ_send(:) = 0 ; requ_recv(:) = 0
    stat_send(:,:) = 0 ; stat_recv(:,:) = 0
    !$omp end workshare

    send_count = 0 ; recv_count = 0
    !$omp do private(i, j, k) reduction(+:send_count, recv_count)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        send_count = send_count + 1
      end do
      do k = rind(i-1)+1, rind(i)
        recv_count = recv_count + 1
      end do
    end do
    !$omp end do
    !$omp end parallel
    allocate(sbufreal(send_count), rbufreal(recv_count))
    !$omp parallel
    !$omp workshare
    sbufreal(:) = DZERO ; rbufreal(:) = DZERO
    !$omp end workshare

    !$omp do private(i, j, jj)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        jj = sitem(j)
        sbufreal(j) = inv(jj)
      end do
    end do
    !$omp end do

    ierr = 0
    !$omp do private(i, is_sta, is_end, sbuflen)
    do i = 1, nbtot
      is_sta = sind(i-1)+1 ; is_end = sind(i)
      sbuflen = is_end - is_sta + 1
      if (sbuflen /= 0) then
        call MPI_ISEND(sbufreal(is_sta), sbuflen, MPI_REAL8, nbnum(i), 0, my_comm,&
                       requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, ir_sta, ir_end, rbuflen)
    do i = 1, nbtot
      ir_sta = rind(i-1)+1 ; ir_end = rind(i)
      rbuflen = ir_end - ir_sta + 1
      if (rbuflen /= 0) then
        call MPI_IRECV(rbufreal(ir_sta), rbuflen, MPI_REAL8, nbnum(i), 0, my_comm,&
                       requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(nbtot, requ_recv, stat_recv, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Receive neighbor real8 values.")
      end if
    end if

    call MPI_WAITALL(nbtot, requ_send, stat_send, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send neighbor real8 values.")
      end if
    end if

    !$omp do private(i, k, kk)
    do i = 1, nbtot
      do k = rind(i-1)+1, rind(i)
        kk = ritem(k)
        outv(kk) = rbufreal(k)
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_send, stat_recv, sbufreal, rbufreal)

  end subroutine senrec_neib_r8

  subroutine senrec_face_i4(nbtot, nbnum, sind, rind, sitem, ritem, ind, outd)
  !***************************************************************************************
  ! senrec_face_i4 -- Send and Receive face integer value
  !***************************************************************************************
    ! -- module
    use constval_module, only: FACE
    ! -- inout
    integer(I4), intent(in) :: nbtot
    integer(I4), intent(in) :: nbnum(:), sind(0:), rind(0:), sitem(:), ritem(:)
    integer(I4), intent(in) :: ind(:,:)
    integer(I4), intent(out) :: outd(:,:)
    ! -- local
    integer(I4) :: i, j, k, jj, kk, ierr
    integer(I4) :: sbuflen, rbuflen, send_count, recv_count
    integer(I4) :: is_sta, is_end, ir_sta, ir_end
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_send(:,:), stat_recv(:,:)
    integer(I4), allocatable :: sbufreal(:), rbufreal(:)
    !-------------------------------------------------------------------------------------
    allocate(requ_send(nbtot), requ_recv(nbtot))
    allocate(stat_send(MPI_STATUS_SIZE,nbtot), stat_recv(MPI_STATUS_SIZE,nbtot))
    !$omp parallel
    !$omp workshare
    stat_send(:,:) = 0 ; stat_recv(:,:) = 0
    requ_send(:) = 0 ; requ_recv(:) = 0
    !$omp end workshare

    send_count = 0 ; recv_count = 0
    !$omp do private(i, j, k) reduction(+:send_count, recv_count)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        send_count = send_count + 1
      end do
      do k = rind(i-1)+1, rind(i)
        recv_count = recv_count + 1
      end do
    end do
    !$omp end do
    !$omp end parallel
    allocate(sbufreal(send_count*FACE), rbufreal(recv_count*FACE))
    !$omp parallel
    !$omp workshare
    sbufreal(:) = 0 ; rbufreal(:) = 0
    !$omp end workshare

    !$omp do private(i, j, jj, k, kk)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        jj = sitem(j)
        do k = 1, FACE
          kk = (j-1)*FACE + k
          sbufreal(kk) = ind(jj,k)
        end do
      end do
    end do
    !$omp end do

    ierr = 0
    !$omp do private(i, is_sta, is_end, sbuflen)
    do i = 1, nbtot
      is_sta = sind(i-1)*FACE+1 ; is_end = sind(i)*FACE
      sbuflen = is_end - is_sta + 1
      if (sbuflen /= 0) then
        call MPI_ISEND(sbufreal(is_sta), sbuflen, MPI_INTEGER, nbnum(i), 0, my_comm,&
                       requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, ir_sta, ir_end, rbuflen)
    do i = 1, nbtot
      ir_sta = rind(i-1)*FACE+1 ; ir_end = rind(i)*FACE
      rbuflen = ir_end - ir_sta + 1
      if (rbuflen /= 0) then
        call MPI_IRECV(rbufreal(ir_sta), rbuflen, MPI_INTEGER, nbnum(i), 0, my_comm,&
                       requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(nbtot, requ_recv, stat_recv, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Receive face integer values.")
      end if
    end if

    call MPI_WAITALL(nbtot, requ_send, stat_send, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send face integer values.")
      end if
    end if

    !$omp do private(i, j, jj, k, kk)
    do i = 1, nbtot
      do j = rind(i-1)+1, rind(i)
        jj = ritem(j)
        do k = 1, FACE
          kk = (j-1)*FACE + k
          outd(jj,k) = rbufreal(kk)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_send, stat_recv, sbufreal, rbufreal)

  end subroutine senrec_face_i4

  subroutine senrec_face_r4(nbtot, nbnum, sind, rind, sitem, ritem, ind, outd)
  !***************************************************************************************
  ! senrec_face_r4 -- Send and Receive face real4 value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: nbtot
    integer(I4), intent(in) :: nbnum(:), sind(0:), rind(0:), sitem(:), ritem(:)
    real(SP), intent(in) :: ind(:,:)
    real(SP), intent(out) :: outd(:,:)
    ! -- local
    integer(I4) :: i, j, k, jj, kk, ierr
    integer(I4) :: sbuflen, rbuflen, send_count, recv_count
    integer(I4) :: is_sta, is_end, ir_sta, ir_end
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_send(:,:), stat_recv(:,:)
    real(SP), allocatable :: sbufreal(:), rbufreal(:)
    !-------------------------------------------------------------------------------------
    allocate(requ_send(nbtot), requ_recv(nbtot))
    allocate(stat_send(MPI_STATUS_SIZE,nbtot), stat_recv(MPI_STATUS_SIZE,nbtot))
    !$omp parallel
    !$omp workshare
    stat_send(:,:) = 0 ; stat_recv(:,:) = 0
    requ_send(:) = 0 ; requ_recv(:) = 0
    !$omp end workshare

    send_count = 0 ; recv_count = 0
    !$omp do private(i, j, k) reduction(+:send_count, recv_count)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        send_count = send_count + 1
      end do
      do k = rind(i-1)+1, rind(i)
        recv_count = recv_count + 1
      end do
    end do
    !$omp end do
    !$omp end parallel
    allocate(sbufreal(send_count*FACE), rbufreal(recv_count*FACE))
    !$omp parallel
    !$omp workshare
    sbufreal(:) = SZERO ; rbufreal(:) = SZERO
    !$omp end workshare

    !$omp do private(i, j, jj, k, kk)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        jj = sitem(j)
        do k = 1, FACE
          kk = (j-1)*FACE + k
          sbufreal(kk) = ind(jj,k)
        end do
      end do
    end do
    !$omp end do

    ierr = 0
    !$omp do private(i, is_sta, is_end, sbuflen)
    do i = 1, nbtot
      is_sta = sind(i-1)*FACE+1 ; is_end = sind(i)*FACE
      sbuflen = is_end - is_sta + 1
      if (sbuflen /= 0) then
        call MPI_ISEND(sbufreal(is_sta), sbuflen, MPI_REAL4, nbnum(i), 0, my_comm,&
                       requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, ir_sta, ir_end, rbuflen)
    do i = 1, nbtot
      ir_sta = rind(i-1)*FACE+1 ; ir_end = rind(i)*FACE
      rbuflen = ir_end - ir_sta + 1
      if (rbuflen /= 0) then
        call MPI_IRECV(rbufreal(ir_sta), rbuflen, MPI_REAL4, nbnum(i), 0, my_comm,&
                       requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(nbtot, requ_recv, stat_recv, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Receive face real4 values.")
      end if
    end if

    call MPI_WAITALL(nbtot, requ_send, stat_send, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send face real4 values.")
      end if
    end if

    !$omp do private(i, j, jj, k, kk)
    do i = 1, nbtot
      do j = rind(i-1)+1, rind(i)
        jj = ritem(j)
        do k = 1, FACE
          kk = (j-1)*FACE + k
          outd(jj,k) = rbufreal(kk)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_send, stat_recv, sbufreal, rbufreal)

  end subroutine senrec_face_r4

  subroutine senrec_face_r8(nbtot, nbnum, sind, rind, sitem, ritem, ind, outd)
  !***************************************************************************************
  ! senrec_face_r8 -- Send and Receive face real8 value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: nbtot
    integer(I4), intent(in) :: nbnum(:), sind(0:), rind(0:), sitem(:), ritem(:)
    real(DP), intent(in) :: ind(:,:)
    real(DP), intent(out) :: outd(:,:)
    ! -- local
    integer(I4) :: i, j, k, jj, kk, ierr
    integer(I4) :: sbuflen, rbuflen, send_count, recv_count
    integer(I4) :: is_sta, is_end, ir_sta, ir_end
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_send(:,:), stat_recv(:,:)
    real(DP), allocatable :: sbufreal(:), rbufreal(:)
    !-------------------------------------------------------------------------------------
    allocate(requ_send(nbtot), requ_recv(nbtot))
    allocate(stat_send(MPI_STATUS_SIZE,nbtot), stat_recv(MPI_STATUS_SIZE,nbtot))
    !$omp parallel
    !$omp workshare
    stat_send(:,:) = 0 ; stat_recv(:,:) = 0
    requ_send(:) = 0 ; requ_recv(:) = 0
    !$omp end workshare

    send_count = 0 ; recv_count = 0
    !$omp do private(i, j, k) reduction(+:send_count, recv_count)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        send_count = send_count + 1
      end do
      do k = rind(i-1)+1, rind(i)
        recv_count = recv_count + 1
      end do
    end do
    !$omp end do
    !$omp end parallel
    allocate(sbufreal(send_count*FACE), rbufreal(recv_count*FACE))
    !$omp parallel
    !$omp workshare
    sbufreal(:) = DZERO ; rbufreal(:) = DZERO
    !$omp end workshare

    !$omp do private(i, j, jj, k, kk)
    do i = 1, nbtot
      do j = sind(i-1)+1, sind(i)
        jj = sitem(j)
        do k = 1, FACE
          kk = (j-1)*FACE + k
          sbufreal(kk) = ind(jj,k)
        end do
      end do
    end do
    !$omp end do

    ierr = 0
    !$omp do private(i, is_sta, is_end, sbuflen)
    do i = 1, nbtot
      is_sta = sind(i-1)*FACE+1 ; is_end = sind(i)*FACE
      sbuflen = is_end - is_sta + 1
      if (sbuflen /= 0) then
        call MPI_ISEND(sbufreal(is_sta), sbuflen, MPI_REAL8, nbnum(i), 0, my_comm,&
                       requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, ir_sta, ir_end, rbuflen)
    do i = 1, nbtot
      ir_sta = rind(i-1)*FACE+1 ; ir_end = rind(i)*FACE
      rbuflen = ir_end - ir_sta + 1
      if (rbuflen /= 0) then
        call MPI_IRECV(rbufreal(ir_sta), rbuflen, MPI_REAL8, nbnum(i), 0, my_comm,&
                       requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(nbtot, requ_recv, stat_recv, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Receive face real8 values.")
      end if
    end if

    call MPI_WAITALL(nbtot, requ_send, stat_send, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Send face real8 values.")
      end if
    end if

    !$omp do private(i, j, jj, k, kk)
    do i = 1, nbtot
      do j = rind(i-1)+1, rind(i)
        jj = ritem(j)
        do k = 1, FACE
          kk = (j-1)*FACE + k
          outd(jj,k) = rbufreal(kk)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_send, stat_recv, sbufreal, rbufreal)

  end subroutine senrec_face_r8

  subroutine bcast_retn_clas(clasn, retn_name, reta, retn, res)
  !***************************************************************************************
  ! bcast_retn_clas -- Bcast retention classification value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: clasn
    character(*), intent(inout) :: retn_name(:)
    real(SP), intent(inout) :: reta(:), retn(:), res(:)
    ! -- local
    integer(I4) :: i, ierr, char_leng
    !-------------------------------------------------------------------------------------
    ierr = 0
    do i = 1, clasn
      char_leng = len_trim(retn_name(i))
      call MPI_BCAST(retn_name(i), char_leng, MPI_CHARACTER, 0, my_comm, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (my_rank == 0) then
          call write_err_stop("Broadcast retention name "//retn_name(i)//".")
        end if
      end if
    end do

    call MPI_BCAST(reta, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast van genuchten parameter alpha.")
      end if
    end if

    call MPI_BCAST(retn, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast van genuchten parameter n.")
      end if
    end if

    call MPI_BCAST(res, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast residual water content.")
      end if
    end if

  end subroutine bcast_retn_clas

  subroutine bcast_parm_clas(clasn, parm_name, ksx, ksy, ksz, ss, ts)
  !***************************************************************************************
  ! bcast_parm_clas -- Bcast parameter classification value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: clasn
    character(*), intent(inout) :: parm_name(:)
    real(SP), intent(inout) :: ksx(:), ksy(:), ksz(:), ss(:), ts(:)
    ! -- local
    integer(I4) :: i, ierr, char_leng
    !-------------------------------------------------------------------------------------
    ierr = 0
    do i = 1, clasn
      char_leng = len_trim(parm_name(i))
      call MPI_BCAST(parm_name(i), char_leng, MPI_CHARACTER, 0, my_comm, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (my_rank == 0) then
          call write_err_stop("Broadcast parameter name "//parm_name(i)//".")
        end if
      end if
    end do

    call MPI_BCAST(ksx, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast saturated hydraulic conductivity value in x direction.")
      end if
    end if

    call MPI_BCAST(ksy, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast saturated hydraulic conductivity value in y direction.")
      end if
    end if

    call MPI_BCAST(ksz, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast saturated hydraulic conductivity value in z direction.")
      end if
    end if

    call MPI_BCAST(ss, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast specific storage value.")
      end if
    end if

    call MPI_BCAST(ts, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast porosity value.")
      end if
    end if

  end subroutine bcast_parm_clas

  subroutine bcast_init_dep()
  !***************************************************************************************
  ! bcast_init_dep -- Bcast initial depth
  !***************************************************************************************
    ! -- module
    use initial_module, only: st_init
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(st_init%depth, 1, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast initial depth value.")
      end if
    end if

  end subroutine bcast_init_dep

  subroutine bcast_clas_val(clasn, cname, cval)
  !***************************************************************************************
  ! bcast_clas_val -- Bcast classification value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: clasn
    character(*), intent(inout) :: cname(:)
    real(SP), intent(inout) :: cval(:)
    ! -- local
    integer(I4) :: i, ierr, char_leng
    !-------------------------------------------------------------------------------------
    ierr = 0
    do i = 1, clasn
      char_leng = len_trim(cname(i))
      call MPI_BCAST(cname(i), char_leng, MPI_CHARACTER, 0, my_comm, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (my_rank == 0) then
          call write_err_stop("Broadcast classification name "//cname(i)//".")
        end if
      end if
    end do

    call MPI_BCAST(cval, clasn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast classification value.")
      end if
    end if

  end subroutine bcast_clas_val

  subroutine bcast_2dpoint(pointn, pi, pj, pval)
  !***************************************************************************************
  ! bcast_2dpoint -- Bcast 2d point value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: pointn
    integer(I4), intent(inout) :: pi(:), pj(:)
    real(SP), intent(inout) :: pval(:)
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(pi, pointn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 2d point position in x direction.")
      end if
    end if

    call MPI_BCAST(pj, pointn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 2d point position in y direction.")
      end if
    end if

    call MPI_BCAST(pval, pointn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 2d point value.")
      end if
    end if


  end subroutine bcast_2dpoint

  subroutine bcast_3dpoint(pointn, pi, pj, pk, pval)
  !***************************************************************************************
  ! bcast_3dpoint -- Bcast 3d point value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: pointn
    integer(I4), intent(inout) :: pi(:), pj(:), pk(:)
    real(SP), intent(inout) :: pval(:)
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(pi, pointn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 3d point position in x direction.")
      end if
    end if

    call MPI_BCAST(pj, pointn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 3d point position in y direction.")
      end if
    end if

    call MPI_BCAST(pk, pointn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 3d point position in z direction.")
      end if
    end if

    call MPI_BCAST(pval, pointn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 3d point value.")
      end if
    end if

  end subroutine bcast_3dpoint

  subroutine bcast_wellpoint(wpn, wid, wi, wj, wks, wke, wval)
  !***************************************************************************************
  ! bcast_wellpoint -- Bcast well point value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: wpn
    integer(I4), intent(inout) :: wid(:), wi(:), wj(:), wks(:), wke(:)
    real(SP), intent(inout) :: wval(:)
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(wid, wpn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast well point id.")
      end if
    end if

    call MPI_BCAST(wi, wpn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast well point position in x direction.")
      end if
    end if

    call MPI_BCAST(wj, wpn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast well point position in y direction.")
      end if
    end if

    call MPI_BCAST(wks, wpn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast well point start position in z direction.")
      end if
    end if

    call MPI_BCAST(wke, wpn, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast well point end position in z direction.")
      end if
    end if

    call MPI_BCAST(wval, wpn, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast well point value.")
      end if
    end if

  end subroutine bcast_wellpoint

  subroutine scatter_i4xy(loc_num, l2g_ij, xycell_in, xycalc_out)
  !***************************************************************************************
  ! scatter_i4xy -- Scatter integer xy value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: l2g_ij(:)
    integer(I4), intent(in) :: xycell_in(:)
    integer(I4), intent(out) :: xycalc_out(:)
    ! -- local
    integer(I4) :: i, s_num, ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(xycell_in, st_grid%nx*st_grid%ny, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 2d integer value.")
      end if
    end if

    !$omp parallel do private(i, s_num)
    do i = 1, loc_num
      s_num = l2g_ij(i)
      xycalc_out(i) = xycell_in(s_num)
    end do
    !$omp end parallel do

  end subroutine scatter_i4xy

  subroutine scatter_r4xy(loc_num, l2g_ij, xycell_in, xycalc_out)
  !***************************************************************************************
  ! scatter_r4xy -- Scatter real4 xy value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: l2g_ij(:)
    real(SP), intent(in) :: xycell_in(:)
    real(SP), intent(out) :: xycalc_out(:)
    ! -- local
    integer(I4) :: i, s_num, ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(xycell_in, st_grid%nx*st_grid%ny, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 2d real4 value.")
      end if
    end if

    !$omp parallel do private(i, s_num)
    do i = 1, loc_num
      s_num = l2g_ij(i)
      xycalc_out(i) = xycell_in(s_num)
    end do
    !$omp end parallel do

  end subroutine scatter_r4xy

  subroutine scatter_r8xy(loc_num, l2g_ij, xycell_in, xycalc_out)
  !***************************************************************************************
  ! scatter_r8xy -- Scatter real8 xy value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: l2g_ij(:)
    real(DP), intent(in) :: xycell_in(:)
    real(DP), intent(out) :: xycalc_out(:)
    ! -- local
    integer(I4) :: i, s_num, ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(xycell_in, st_grid%nx*st_grid%ny, MPI_REAL8, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 2d real8 value.")
      end if
    end if

    !$omp parallel do private(i, s_num)
    do i = 1, loc_num
      s_num = l2g_ij(i)
      xycalc_out(i) = xycell_in(s_num)
    end do
    !$omp end parallel do

  end subroutine scatter_r8xy

  subroutine scatter_i4xyz(loc_num, l2g_ijk, xyzcalc_in, xyzcalc_out)
  !***************************************************************************************
  ! scatter_i4xyz -- Scatter integer xyz value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: l2g_ijk(:)
    integer(I4), intent(in) :: xyzcalc_in(:)
    integer(I4), intent(out) :: xyzcalc_out(:)
    ! -- local
    integer(I4) :: i, c_num, ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(xyzcalc_in, st_grid%nxyz, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 3d integer value.")
      end if
    end if

    !$omp parallel do private(i, c_num)
    do i = 1, loc_num
      c_num = l2g_ijk(i)
      xyzcalc_out(i) = xyzcalc_in(c_num)
    end do
    !$omp end parallel do

  end subroutine scatter_i4xyz

  subroutine scatter_r4xyz(loc_num, l2g_ijk, xyzcalc_in, xyzcalc_out)
  !***************************************************************************************
  ! scatter_r4xyz -- Scatter real4 xyz value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: l2g_ijk(:)
    real(SP), intent(in) :: xyzcalc_in(:)
    real(SP), intent(out) :: xyzcalc_out(:)
    ! -- local
    integer(I4) :: i, c_num, ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(xyzcalc_in, st_grid%nxyz, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 3d real4 value.")
      end if
    end if

    !$omp parallel do private(i, c_num)
    do i = 1, loc_num
      c_num = l2g_ijk(i)
      xyzcalc_out(i) = xyzcalc_in(c_num)
    end do
    !$omp end parallel do

  end subroutine scatter_r4xyz

  subroutine scatter_r8xyz(loc_num, l2g_ijk, xyzcalc_in, xyzcalc_out)
  !***************************************************************************************
  ! scatter_r8xyz -- Scatter real8 xyz value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: l2g_ijk(:)
    real(DP), intent(in) :: xyzcalc_in(:)
    real(DP), intent(out) :: xyzcalc_out(:)
    ! -- local
    integer(I4) :: i, c_num, ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(xyzcalc_in, st_grid%nxyz, MPI_REAL8, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast 3d real8 value.")
      end if
    end if

    !$omp parallel do private(i, c_num)
    do i = 1, loc_num
      c_num = l2g_ijk(i)
      xyzcalc_out(i) = xyzcalc_in(c_num)
    end do
    !$omp end parallel do

  end subroutine scatter_r8xyz

  subroutine bcast_solval()
  !***************************************************************************************
  ! bcast_solval -- Bcast solution value
  !***************************************************************************************
    ! -- module
    use mpi_utility, only: bcast_val
    use initial_module, only: maxout_iter, maxinn_iter, amg_nlevel, maxvcy_iter,&
                              max_sweep, criteria, jac_omega, amg_theta
    use read_input, only: temp_maxinn_iter
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    ! -- Bcast scalar value (val)
      call bcast_val(maxout_iter, " maximum outer iteration number")

    if (st_sim%sim_type /= -1) then
      ! -- Bcast scalar value (val)
        call bcast_val(maxinn_iter, " maximum inner iteration number")
    end if

    ! -- Bcast scalar value (val)
      call bcast_val(amg_nlevel, " multigrid level number")

    ! -- Bcast scalar value (val)
      call bcast_val(st_sim%ini_step, " initial time step value")
    ! -- Bcast scalar value (val)
      call bcast_val(st_sim%max_step, " maximun time step value")
    ! -- Bcast scalar value (val)
      call bcast_val(st_sim%inc_fact, " increment multiplier value")
    ! -- Bcast scalar value (val)
      call bcast_val(criteria, " outer criteria (max norm) value")

    if (precon_type == 1) then
      ! -- Bcast scalar value (val)
        call bcast_val(maxvcy_iter, " maximum v-cycle number")
      ! -- Bcast scalar value (val)
        call bcast_val(max_sweep, " maximum sweep number")
      ! -- Bcast scalar value (val)
        call bcast_val(jac_omega, " jacobian omega value")
      ! -- Bcast scalar value (val)
        call bcast_val(amg_theta, " amg theta value")
    end if

    ! -- Bcast scalar value (val)
      call bcast_val(temp_maxinn_iter, " maximum inner iteration number")

  end subroutine bcast_solval

end module mpi_set
