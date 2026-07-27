module set_cell
  ! -- modules
  use kind_module, only: I4, SP
  use constval_module, only: SNOVAL, CHALEN
  use types_module, only: conn_set, gmap_set
  use utility_module, only: st_mpi, iquick_sort, iquick_sort2, write_err_stop
  use utility_module, only: gmap_init, gmap_put, gmap_get, gmap_free
  use initial_module, only: st_sim, st_grid, st_clas

  implicit none
  private
  public :: set_cell_info, get_cals_grid, get_calc_grid
  integer(I4), public :: amg_setflag
  integer(I4), public :: ncalc, ncals, ncell, nsurf, no_ncalc, no_ncals
  integer(I4), public :: seal_snum, seal_cnum, neib_mpi_totn, neib_ncals, neib_ncalc
#ifdef MPI_MSG
  integer(I4), allocatable, public :: send_cind(:), recv_cind(:)
  integer(I4), allocatable, public :: send_citem(:), recv_citem(:)
  integer(I4), allocatable, public :: neib_num(:), send2recv(:), calc2recv(:)
  integer(I4), allocatable, public :: loc2glo_nos(:), loc2glo_noc(:)
  integer(I4), allocatable, public :: loc2unk_ij(:)
#endif
  type(conn_set), public :: st_conn

  ! -- local
  integer(I4) :: nc_unknow
  integer(I4) :: totnreg, loc_regn, neib_hash
  logical :: div_reg_3d
  integer(I4), allocatable :: glob_reg_flag(:), glob_mpi_flag(:)
  integer(I4), allocatable :: calc_end(:), loc_nreg(:)
  integer(I4), allocatable :: l2g_ij(:), l2g_ijk(:)
  integer(I4), allocatable :: neib_hash_key(:), neib_hash_mpi(:), neib_hash_reg(:)
  integer(I4), allocatable :: own_calc_glo(:), own_calc_reg(:), own_nocal_glo(:)
  integer(I4) :: read_sta, read_num
  integer(I4), allocatable :: read_reg(:), read_mpi(:)
  type(gmap_set) :: surf_map

  contains

  subroutine set_cell_info()
  !*********************************************************************************************
  ! set_cell_info -- Set cell information
  !*********************************************************************************************
    ! -- modules
    use utility_module, only: open_new_wtxt, close_file
    use initial_module, only: out_type, st_ctrl, st_out_type, st_out_path
    use check_condition, only: check_calc_region, read_sea_points
#ifdef MPI_MSG
    use initial_module, only: st_in_path
    use mpi_read, only: read_dist_seaval
#endif
#ifdef MPI_MSG
    use mpi_utility, only: barrier_proc, mpisum_val, bcast_val, bcast_file
    use mpi_set, only: bcast_sim_flag, bcast_xyz_num, bcast_clas_set, bcast_glob_flag,&
                       set_calc_view, set_seal_view, set_rest_view, set_write_fview,&
                       senrec_reg_info, senrec_grid_num
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i, j, n, nxy
    integer(I4) :: sea_num, sea_pos, sea_mode, sea_ptn, do_fallback
    integer(I4), allocatable :: sea_glo(:), sp_i(:), sp_j(:), sp_k(:)
#ifdef MPI_MSG
    integer(I4) :: sea_li
    real(SP), allocatable :: read_seaval(:)
    character(CHALEN) :: seal_path_bc
#endif
    integer(I4) :: mpi_ncals, mpi_ncalc
    integer(I4) :: sta_calc, end_calc, i_num, j_num, k_num, pro_nreg
    integer(I4) :: calg_num
    integer(I4), allocatable :: cur_nreg(:)
#ifdef MPI_MSG
    integer(I4), allocatable :: sort_sglo(:), sort_cglo(:)
    integer(I4), allocatable :: loc2unk_ijk(:)
#endif
    !-------------------------------------------------------------------------------------------
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Bcast simulation flag (sim_flag)
        call bcast_sim_flag()
      ! -- Bcast xyz number (xyz_num)
        call bcast_xyz_num()
    end if
#endif

    if (st_ctrl%precon_type == 1) then
      amg_setflag = 0
    else
      st_ctrl%amg_nlevel = 1
    end if

    if (st_ctrl%noclas_flag /= 1) then
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Bcast classification setting (clas_set)
          call bcast_clas_set()
      end if
#endif
    else
      st_clas%totn = 1
    end if

    allocate(glob_reg_flag(st_grid%nxyz), glob_mpi_flag(st_grid%nxyz))
    !$omp parallel do private(i)
    do i = 1, st_grid%nxyz
      glob_reg_flag(i) = 0 ; glob_mpi_flag(i) = 0
    end do
    !$omp end parallel do

    nxy = st_grid%nx*st_grid%ny

    ! -- Partition dimension (2d columns / 3d cells). 3d is selected by:
    !      div_reg_3d = st_sim%reg_neib /= 1 .and. &
    !                   (st_sim%reg_type == in_type(5) .or. st_sim%reg_type == in_type(6) .or. &
    !                    st_mpi%totn >= nxy)
    !    Forced 2d until the 3d downstream (column maps, surface comm, seepage) is ready.
    div_reg_3d = .false.

    ! -- Set global region (glob_reg)
      call set_glob_reg()

#ifdef MPI_MSG
    ! -- Region flag must be on all ranks before the distributed partition
    if (st_mpi%totn /= 1) then
      call bcast_glob_flag(totnreg, glob_reg_flag, glob_mpi_flag)
    end if
    ! -- Divide calculation region, distributed (calc_reg_dist)
      call div_calc_reg_dist()
    ! -- Divide no calculation flag, distributed (nocalc_dist)
      call div_nocalc_dist()
#else
    ! -- Divide calculation region (calc_reg)
      call div_calc_reg()
    ! -- Divide no calculation flag (nocalc_flag)
      call div_nocalc_flag()
#endif

    ! -- Build this rank's read range (PRE-sea) via the unit->cell exchange
      call build_read_range()

    sea_num = 0 ; sea_mode = 0 ; do_fallback = 0
    if (st_mpi%rank == 0) then
      call read_sea_points(sea_mode, sea_ptn, sp_i, sp_j, sp_k)
      do_fallback = 0
      if (sea_mode == 9) then
        do_fallback = 1
      end if
#ifndef MPI_MSG
      if (sea_mode == 2 .or. sea_mode == 3) then
        do_fallback = 1
      end if
#endif
      ! -- Check the calculation regin (calc_region)
        call check_calc_region(glob_reg_flag)
        where (glob_reg_flag(:) == 0 .and. glob_mpi_flag(:) > 0) glob_mpi_flag(:) = 0
      if (do_fallback == 1) then
        sea_num = count(glob_mpi_flag(:) == 0)
        allocate(sea_glo(sea_num))
        sea_pos = 0
        do i = 1, st_grid%nxyz
          if (glob_mpi_flag(i) == 0) then
            sea_pos = sea_pos + 1 ; sea_glo(sea_pos) = i
          end if
        end do
      end if
    end if

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      call bcast_val(sea_mode, "sea mode")
      call bcast_val(do_fallback, "sea fallback")
      if (sea_mode == 1) then
        call bcast_val(sea_ptn, "sea point count")
        if (st_mpi%rank /= 0) then
          allocate(sp_i(sea_ptn+1), sp_j(sea_ptn+1), sp_k(sea_ptn+1))
        end if
        if (sea_ptn > 0) then
          call bcast_val(sp_i, "sea pi") ; call bcast_val(sp_j, "sea pj")
          call bcast_val(sp_k, "sea pk")
        end if
      else if (do_fallback == 1) then
        call bcast_val(sea_num, "sea count")
        if (st_mpi%rank /= 0) then
          allocate(sea_glo(sea_num))
        end if
        if (sea_num > 0) then
          call bcast_val(sea_glo, "sea list")
        end if
      end if
      ! -- Bcast global flag (glob_flag) (scaffold for the shadow check)
        call bcast_glob_flag(totnreg, glob_reg_flag, glob_mpi_flag)
    end if
#endif

    ! -- Apply sea to the read range, then shadow-verify against the post-sea global slice
    if (sea_mode == 1) then
      call punch_sea_points(sea_ptn, sp_i, sp_j, sp_k)
    else if (do_fallback == 1) then
      call punch_sea_read_range(sea_num, sea_glo)
#ifdef MPI_MSG
    else if (sea_mode == 2 .or. sea_mode == 3) then
      ! -- the sea path is read only on rank 0 (read_main_file) and bcast later in set_boundary,
      !    so broadcast it here before the collective MPI-IO read
      seal_path_bc = ""
      if (st_mpi%rank == 0) then
        seal_path_bc = st_in_path%seal
      end if
      if (st_mpi%totn /= 1) then
        call bcast_file(seal_path_bc, "sea level path")
      end if
      allocate(read_seaval(read_num+1))
      call read_dist_seaval(sea_mode == 3, trim(seal_path_bc), read_sta, read_num, read_seaval)
      do sea_li = 1, read_num
        if (read_seaval(sea_li) /= SNOVAL .and. read_reg(sea_li) > 0) then
          read_reg(sea_li) = 0 ; read_mpi(sea_li) = 0
        end if
      end do
      deallocate(read_seaval)
#endif
    end if
      call verify_read_range()
    if (allocated(sea_glo)) then
      deallocate(sea_glo)
    end if
    if (allocated(sp_i)) then
      deallocate(sp_i, sp_j, sp_k)
    end if

    ! -- Redistribute owned cells to their owner ranks (owned2owner)
      call redist_owned2owener()

    ncalc = count(glob_mpi_flag(:) == st_mpi%rank+1)
    ncals = count(glob_mpi_flag(1:nxy) == st_mpi%rank+1)

    allocate(st_conn%glo2loc_ijk(st_grid%nxyz), l2g_ijk(ncalc))
    allocate(l2g_ij(ncals))
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_grid%nxyz
      st_conn%glo2loc_ijk(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncalc
      l2g_ijk(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncals
      l2g_ij(i) = 0
    end do
    !$omp end do
    !$omp end parallel

    ! -- Set relationship between global&local (rel_gloloc)
      call set_rel_gloloc()

    neib_ncals = 0 ; neib_ncalc = 0 ; neib_mpi_totn = 0

    ! -- Set hash neighbor (hash_neib)
      call set_hash_neib()

#ifdef MPI_MSG
    allocate(calc2recv(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      calc2recv(i) = 0
    end do
    !$omp end parallel do

    ! -- Set mpi relationship (mpi_rel)
      call set_mpi_rel()
#endif

    ! -- Set relationship of sea region (rel_seareg)
      call set_rel_seareg()

    deallocate(neib_hash_key, neib_hash_mpi, neib_hash_reg)

    deallocate(glob_reg_flag)
    deallocate(read_reg, read_mpi)

    mpi_ncals = ncals + neib_ncals ; mpi_ncalc = ncalc + neib_ncalc
    nsurf = mpi_ncals + seal_snum ; ncell = mpi_ncalc + seal_cnum

    allocate(st_conn%loc2glo_ij(nsurf), st_conn%loc2glo_ijk(ncell))
    !$omp parallel
    !$omp do private(i)
    do i = 1, nsurf
      st_conn%loc2glo_ij(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncell
      st_conn%loc2glo_ijk(i) = 0
    end do
    !$omp end do
    !$omp end parallel

    st_conn%loc2glo_ij(:mpi_ncals) = l2g_ij(:)
    st_conn%loc2glo_ijk(:mpi_ncalc) = l2g_ijk(:)

    do i = 1, surf_map%table_num
      if (surf_map%table_key(i) /= 0 .and. surf_map%table_val(i) > mpi_ncals) then
        st_conn%loc2glo_ij(surf_map%table_val(i)) = surf_map%table_key(i)
      end if
    end do
    do i = 1, st_grid%nxyz
      if (st_conn%glo2loc_ijk(i) > mpi_ncalc) then
        st_conn%loc2glo_ijk(st_conn%glo2loc_ijk(i)) = i
      end if
    end do

#ifdef MPI_MSG
    ! -- add output-only orphan sea cells (punctured sea with no active neighbour, in no
    !    calc/nocalc/seal list) to the nocalc list so the mpi view-write fills SNOVAL there
    !    like the serial full-grid write
      call add_orphan_nocal()
#endif

    no_ncalc = size(own_nocal_glo)
    no_ncals = 0
    do i = 1, no_ncalc
      if (own_nocal_glo(i) <= nxy) then
        no_ncals = no_ncals + 1
      end if
    end do

#ifdef MPI_MSG
    allocate(loc2glo_nos(no_ncals), loc2glo_noc(no_ncalc))
    do i = 1, no_ncalc
      loc2glo_noc(i) = own_nocal_glo(i)
    end do
    do i = 1, no_ncals
      loc2glo_nos(i) = own_nocal_glo(i)
    end do
#endif

    deallocate(l2g_ij, l2g_ijk, glob_mpi_flag)
    deallocate(own_calc_glo, own_calc_reg, own_nocal_glo)

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(ncalc, "calculation number", nc_unknow)
    else
      nc_unknow = ncalc
    end if
#else
    nc_unknow = ncalc
#endif

#ifdef MPI_MSG
    allocate(loc2unk_ij(ncals), loc2unk_ijk(ncalc))

    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      loc2unk_ij(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncalc
      loc2unk_ijk(i) = 0
    end do
    !$omp end do
    !$omp end parallel

    ! -- the unknown number of each own cell, ranked across ranks without any global gather
    call set_dist_unknum(st_conn%loc2glo_ij(1:ncals), ncals, nxy, loc2unk_ij)
    call set_dist_unknum(st_conn%loc2glo_ijk(1:ncalc), ncalc, st_grid%nxyz, loc2unk_ijk)

    allocate(sort_sglo(ncals+no_ncals+seal_snum), sort_cglo(ncalc+no_ncalc+seal_cnum))
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      sort_sglo(i) = st_conn%loc2glo_ij(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, no_ncals
      sort_sglo(ncals+i) = loc2glo_nos(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, seal_snum
      sort_sglo(ncals+no_ncals+i) = st_conn%loc2glo_ij(mpi_ncals+i)
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, ncalc
      sort_cglo(i) = st_conn%loc2glo_ijk(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, no_ncalc
      sort_cglo(ncalc+i) = loc2glo_noc(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, seal_cnum
      sort_cglo(ncalc+no_ncalc+i) = st_conn%loc2glo_ijk(mpi_ncalc+i)
    end do
    !$omp end do
    !$omp end parallel

    call iquick_sort(sort_sglo, 1, ncals+no_ncals+seal_snum)
    call iquick_sort(sort_cglo, 1, ncalc+no_ncalc+seal_cnum)

    ! -- Set calculation view (calc_view)
      call set_calc_view(ncals, ncalc, st_conn%loc2glo_ij, st_conn%loc2glo_ijk)
    ! -- Set seal view (seal_view)
      call set_seal_view(seal_snum, seal_cnum, st_conn%loc2glo_ij(mpi_ncals+1:),&
                         st_conn%loc2glo_ijk(mpi_ncalc+1:))
    ! -- Set restart view (rest_view)
      call set_rest_view(ncalc, nc_unknow, loc2unk_ijk)
    ! -- Set write file view (write_fview)
      call set_write_fview(ncals, ncalc, st_conn%loc2glo_ij, st_conn%loc2glo_ijk,&
                           sort_sglo, sort_cglo)
    deallocate(sort_sglo, sort_cglo)

#ifdef ICI
    deallocate(loc2unk_ijk)
#else
    deallocate(loc2unk_ij, loc2unk_ijk)
#endif
#endif

    call gmap_free(surf_map)

    if (st_ctrl%noclas_flag /= 1) then
      ! -- Set local cell classification (loc_cell_clas)
        call set_loc_cell_clas()
    end if

    if (st_out_type%calg == out_type(1)) then
      if (st_mpi%rank == 0) then
        call open_new_wtxt(st_out_path%calg, "output calculation grid number", calg_num)
        write(calg_num,'(a)') "Calclation_No,I,J,K"
      end if
      do n = 1, st_mpi%totn
#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
          ! -- Send and Receive region information (reg_info)
            call senrec_reg_info(n, loc_nreg, cur_nreg)
          pro_nreg = size(cur_nreg)
        else
          pro_nreg = loc_regn
          allocate(cur_nreg(pro_nreg))
          !$omp parallel do private(i)
          do i = 1, pro_nreg
            cur_nreg(i) = loc_nreg(i)
          end do
          !$omp end parallel do
        end if
#else
        pro_nreg = loc_regn
        allocate(cur_nreg(pro_nreg))
        !$omp parallel do private(i)
        do i = 1, pro_nreg
          cur_nreg(i) = loc_nreg(i)
        end do
        !$omp end parallel do
#endif
        if (st_mpi%rank == 0) then
          write(calg_num,'(a,i0)') "Rank", n-1
        end if
        do i = 1, pro_nreg
          if (st_mpi%rank == 0) then
            write(calg_num,'(a,i0)') "Region", cur_nreg(i)
          end if
          sta_calc = calc_end(i-1) ; end_calc = calc_end(i)
          do j = sta_calc+1, end_calc
            if (n == st_mpi%rank+1) then
              ! -- Get calculation number from grid number (calc_grid)
                call get_calc_grid(j, i_num, j_num, k_num)
            end if
#ifdef MPI_MSG
            if (st_mpi%totn /= 1) then
              ! -- Send and Receive grid number (grid_num)
                call senrec_grid_num(n, i_num, j_num, k_num)
            end if
#endif
            if (st_mpi%rank == 0) then
              write(calg_num,'(*(i0:,","))') j, i_num, j_num, k_num
            end if
          end do
        end do
#ifdef MPI_MSG
        ! -- Barrier process (proc)
          call barrier_proc()
#endif
      end do
      deallocate(cur_nreg)
      if (st_mpi%rank == 0) then
        call close_file(calg_num)
      end if
    end if

    deallocate(calc_end, loc_nreg)

  end subroutine set_cell_info

  subroutine build_clas_flag(clas_k, flag)
  !*********************************************************************************************
  ! build_clas_flag -- Build the 0/1 cell flag of one class from its specs (i,j,k with -1
  !   wildcards). Shared expansion so set_glob_reg gets the same cells the (removed) global
  !   clas array held. (段5-4b)
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: clas_k
    integer(I4), intent(out) :: flag(:)
    ! -- local
    integer(I4) :: i, ii, jj, kk, c_num
    integer(I4) :: clasi, clasj, clask
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, st_grid%nxyz
      flag(i) = 0
    end do
    !$omp end parallel do

    do i = 1, st_clas%num(clas_k)
      clasi = st_clas%i(i,clas_k) ; clasj = st_clas%j(i,clas_k) ; clask = st_clas%k(i,clas_k)
      if (clasi == -1 .and. clasj == -1 .and. clask == -1) then !all cells
        flag(1:st_grid%nxyz) = 1
      else if (clasi == -1 .and. clasj == -1) then !i,j cells
        do jj = 1, st_grid%ny
          do ii = 1, st_grid%nx
            c_num = st_grid%nx*(st_grid%ny*(clask-1) + (jj-1)) + ii
            flag(c_num) = 1
          end do
        end do
      else if (clasi == -1 .and. clask == -1) then !i,k cells
        do kk = 1, st_grid%nz
          do ii = 1, st_grid%nx
            c_num = st_grid%nx*(st_grid%ny*(kk-1) + (clasj-1)) + ii
            flag(c_num) = 1
          end do
        end do
      else if (clasj == -1 .and. clask == -1) then !j,k cells
        do kk = 1, st_grid%nz
          do jj = 1, st_grid%ny
            c_num = st_grid%nx*(st_grid%ny*(kk-1) + (jj-1)) + clasi
            flag(c_num) = 1
          end do
        end do
      else if (clasi == -1) then !only i cell
        do ii = 1, st_grid%nx
          c_num = st_grid%nx*(st_grid%ny*(clask-1) + (clasj-1)) + ii
          flag(c_num) = 1
        end do
      else if (clasj == -1) then !only j cell
        do jj = 1, st_grid%ny
          c_num = st_grid%nx*(st_grid%ny*(clask-1) + (jj-1)) + clasi
          flag(c_num) = 1
        end do
      else if (clask == -1) then !only k cell
        do kk = 1, st_grid%nz
          c_num = st_grid%nx*(st_grid%ny*(kk-1) + (clasj-1)) + clasi
          flag(c_num) = 1
        end do
      else !others
        c_num = st_grid%nx*(st_grid%ny*(clask-1) + (clasj-1)) + clasi
        flag(c_num) = 1
      end if
    end do

  end subroutine build_clas_flag

  subroutine set_glob_reg()
  !*********************************************************************************************
  ! set_glob_reg -- Set global region
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: VARLEN
    use utility_module, only: open_new_rtxt, open_new_rbin
    use initial_module, only: in_type
    use read_module, only: read_2d_calcreg, read_3d_calcreg
#ifdef MPI_MSG
    use mpi_read, only: read_dist_calcreg
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: reg_type, glob_ncalc, calreg_fnum
    integer(I4) :: len_char, first_pos, sta_pos, end_pos, temp_num
    integer(I4), allocatable :: temp_reg(:), type_2d(:), type_3d(:), cflag(:)
    character(1), allocatable :: temp_char(:)
    character(VARLEN), allocatable :: reg_name(:)
    logical, allocatable :: mask(:)
    !-------------------------------------------------------------------------------------------
    reg_type = st_sim%reg_type
    ! -- Only the binary file read (in_type 4/6) is collective; other region types are
    !    built on rank 0 and broadcast, so non-root ranks return early for those
    if (st_mpi%rank /= 0 .and. reg_type /= in_type(4) .and. reg_type /= in_type(6)) then
      return
    end if
    allocate(type_2d(2), type_3d(2))
    type_2d(:) = [in_type(3:4)] ; type_3d(:) = [in_type(5:6)]

    if (reg_type == in_type(0)) then
      glob_ncalc = st_grid%nxyz ; totnreg = 1 ; glob_reg_flag(:) = 1

    else if (reg_type == in_type(1)) then
      glob_ncalc = 0 ; len_char = len(trim(adjustl(st_sim%reg_name)))
      allocate(temp_char(len_char))
      temp_char = transfer(st_sim%reg_name, ' ', size = len_char)
      totnreg = count(temp_char(:) == ",") + 1
      allocate(reg_name(totnreg))
      first_pos = index(st_sim%reg_name, ",")
      if (first_pos == 1) then
        call write_err_stop("check the calculation region name.")
      end if
      sta_pos = 1 ; end_pos = first_pos
      do i = 1, totnreg-1
        reg_name(i) = st_sim%reg_name(sta_pos:end_pos-1) ; sta_pos = end_pos + 1
        end_pos = index(st_sim%reg_name(sta_pos:), ",") + sta_pos - 1
      end do
      reg_name(totnreg) = st_sim%reg_name(sta_pos:)
      deallocate(temp_char)
      allocate(mask(st_grid%nxyz), cflag(st_grid%nxyz))
      do j = 1, st_clas%totn
        call build_clas_flag(j, cflag)
        if (trim(adjustl(st_sim%inact_name)) == trim(adjustl(st_clas%name(j)))) then
          temp_num = sum(cflag(:))
          allocate(temp_reg(temp_num))
          !$omp parallel
          !$omp do private(i)
          do i = 1, temp_num
            temp_reg(i) = -1
          end do
          !$omp end do
          !$omp do private(i)
          do i = 1, st_grid%nxyz
            mask(i) = (cflag(i) == 1)
          end do
          !$omp end do
          !$omp end parallel
          glob_reg_flag(:) = unpack(temp_reg(:), mask(:), glob_reg_flag(:))
          deallocate(temp_reg)
        end if
        temp_num = 0
        do i = 1, totnreg
          if (trim(adjustl(reg_name(i))) == trim(adjustl(st_clas%name(j)))) then
            temp_num = sum(cflag(:))
            allocate(temp_reg(temp_num))
            !$omp parallel
            !$omp do private(k)
            do k = 1, temp_num
              temp_reg(k) = i
            end do
            !$omp end do
            !$omp do private(k)
            do k = 1, st_grid%nxyz
              mask(k) = (cflag(k) == 1)
            end do
            !$omp end do
            !$omp end parallel
            glob_reg_flag(:) = unpack(temp_reg(:), mask(:), glob_reg_flag(:))
            deallocate(temp_reg)
          end if
        end do
        glob_ncalc = glob_ncalc + temp_num
      end do
      deallocate(mask, cflag, reg_name)
      if (.not.(any(glob_reg_flag(:) == -1))) then
        call write_err_stop("check the no calculation region name.")
      end if

    else if (any(reg_type == type_2d(:)) .or. any(reg_type == type_3d(:))) then
      if (reg_type == type_2d(2) .or. reg_type == type_3d(2)) then
        ! -- Binary calc region file
#ifdef MPI_MSG
        ! -- Read calc region file for distributed column (dist_calcreg)
          call read_dist_calcreg(reg_type, glob_reg_flag, glob_ncalc, totnreg)
#else
          call open_new_rbin(1, 1, st_sim%reg_name, "calculation reigion", calreg_fnum)
          if (reg_type == type_2d(2)) then
            call read_2d_calcreg(calreg_fnum, reg_type, glob_reg_flag, glob_ncalc)
          else
            call read_3d_calcreg(calreg_fnum, reg_type, glob_reg_flag, glob_ncalc)
          end if
          totnreg = maxval(glob_reg_flag(:))
#endif
      else
        ! -- Text calc region file
          call open_new_rtxt(1, 1, st_sim%reg_name, "calculation reigion", calreg_fnum)
          if (reg_type == type_2d(1)) then
            call read_2d_calcreg(calreg_fnum, reg_type, glob_reg_flag, glob_ncalc)
          else
            call read_3d_calcreg(calreg_fnum, reg_type, glob_reg_flag, glob_ncalc)
          end if
          totnreg = maxval(glob_reg_flag(:))
      end if
      if (.not.(any(glob_reg_flag(:) == -1))) then
        call write_err_stop("check the no calculation region name.")
      end if
    end if

    deallocate(type_2d, type_3d)

  end subroutine set_glob_reg

  subroutine div_calc_reg()
  !*********************************************************************************************
  ! div_calc_reg -- Divide calculation region (serial; div_reg_3d selects 2d/3d)
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, n, nxy, ntot, ij, reg_id
    integer(I4) :: rough_divn, sum_cals, sum_mpi, pro_num, quot, max_reg, pre_cals
    integer(I4) :: reg_mpi_ncals, reg_mpi_remain, sta_mpi, end_mpi, mpi
    integer(I4), allocatable :: reg_ncals(:), reg_remain(:)
    integer(I4), allocatable :: reg_mpi_num(:), reg_mpi_end(:)
    integer(I4), allocatable :: reg_off(:), all_num(:), reg_pos(:)
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    if (div_reg_3d) then
      ntot = st_grid%nxyz
    else
      ntot = nxy
    end if
    allocate(reg_ncals(totnreg), reg_remain(totnreg))
    allocate(reg_mpi_num(totnreg), reg_mpi_end(0:totnreg))
    reg_mpi_end(0) = 0
    !$omp parallel do private(i)
    do i = 1, totnreg
      reg_ncals(i) = 0 ; reg_remain(i) = 0 ; reg_mpi_num(i) = 0 ; reg_mpi_end(i) = 0
    end do
    !$omp end parallel do
    do ij = 1, ntot
      reg_id = glob_reg_flag(ij)
      if (reg_id >= 1) reg_ncals(reg_id) = reg_ncals(reg_id) + 1
    end do

    rough_divn = sum(reg_ncals(:))/st_mpi%totn
    if (rough_divn == 0) then
      rough_divn = 1
    end if

    sum_cals = 0 ; sum_mpi = 0
    do i = 1, totnreg
      sum_cals = sum_cals + reg_ncals(i)
      quot = sum_cals/rough_divn ; reg_remain(i) = mod(sum_cals, rough_divn)
      if (quot == 0) then
        reg_mpi_num(i) = 1 ; reg_remain(i) = 0 ; reg_mpi_end(i) = sum_mpi + 1
      else if (sum_cals > reg_ncals(i) .and. i /= totnreg) then
        reg_mpi_end(i) = sum_mpi + 1 ; sum_cals = 0
      else
        reg_mpi_num(i) = quot ; sum_mpi = sum_mpi + quot ; reg_mpi_end(i) = sum_mpi
        sum_cals = 0
      end if
    end do

    pro_num = 0
    if (reg_mpi_end(totnreg) < st_mpi%totn) then
      do j = 1, st_mpi%totn-reg_mpi_end(totnreg)
        max_reg = maxloc(reg_remain(:),1) ;
        reg_mpi_num(max_reg) = reg_mpi_num(max_reg) + 1
        if (reg_remain(max_reg) < reg_ncals(max_reg)/2) then
          pro_num = max_reg + 1
        else
          pro_num = max_reg
        end if
        reg_remain(max_reg) = 0
        do k = pro_num, totnreg
          reg_mpi_end(k) = reg_mpi_end(k) + 1
        end do
      end do
    end if
    deallocate(reg_remain)

    allocate(reg_off(0:totnreg))
    reg_off(0) = 0
    do i = 1, totnreg
      reg_off(i) = reg_off(i-1) + reg_ncals(i)
    end do
    allocate(all_num(reg_off(totnreg)), reg_pos(totnreg))
    !$omp parallel do private(i)
    do i = 1, totnreg
      reg_pos(i) = 0
    end do
    !$omp end parallel do
    do ij = 1, ntot
      reg_id = glob_reg_flag(ij)
      if (reg_id >= 1) then
        reg_pos(reg_id) = reg_pos(reg_id) + 1
        all_num(reg_off(reg_id-1) + reg_pos(reg_id)) = ij
      end if
    end do
    do i = 1, totnreg
      sum_cals = 0 ; pre_cals = 0
      sta_mpi = reg_mpi_end(i-1) ; end_mpi = reg_mpi_end(i)
      if (end_mpi == sta_mpi) then
        sta_mpi = sta_mpi - 1
      end if
      reg_mpi_ncals = reg_ncals(i)/(end_mpi-sta_mpi)
      reg_mpi_remain = mod(reg_ncals(i), end_mpi-sta_mpi)
      do n = 1, end_mpi-sta_mpi
        if (n <= reg_mpi_remain) then
          sum_cals = pre_cals + reg_mpi_ncals + 1
        else
          sum_cals = pre_cals + reg_mpi_ncals
        end if
        mpi = sta_mpi + n
        glob_mpi_flag(all_num(reg_off(i-1)+pre_cals+1 : reg_off(i-1)+sum_cals)) = mpi
        pre_cals = sum_cals
      end do
    end do

    deallocate(reg_ncals, reg_mpi_num, reg_mpi_end, reg_off, all_num, reg_pos)

    ! -- 2d: replicate the column ownership to every layer
    if (.not. div_reg_3d) then
      !$omp parallel do private(i)
      do i = 1, st_grid%nz-1
        glob_mpi_flag(nxy*i+1:nxy*(i+1)) = glob_mpi_flag(1:nxy)
      end do
      !$omp end parallel do
    end if

  end subroutine div_calc_reg

#ifdef MPI_MSG
  subroutine div_calc_reg_dist()
  !*********************************************************************************************
  ! div_calc_reg_dist -- Divide calculation region, distributed
  !*********************************************************************************************
    ! -- modules
    use mpi_utility, only: mpisum_val, mpiexscan_val, gather_val
    ! -- inout

    ! -- local
    integer(I4) :: i, k, e, idx, nxy, ntot, base, rem, rang_sta, rang_num
    integer(I4) :: reg_id, pos_in_reg, sta_mpi, end_mpi, reg_nproc, proc_cals, remain_cals
    integer(I4) :: proc_pos, rough_divn, sum_cals, sum_mpi, pro_num, quot, max_reg
    integer(I4), allocatable :: reg_ncals(:), reg_remain(:), reg_mpi_num(:), reg_mpi_end(:)
    integer(I4), allocatable :: my_reg_ncals(:), cals_end(:), seen(:), rang_mpi(:)
    integer(I4), allocatable :: glo_mpi(:)
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    if (div_reg_3d) then
      ntot = st_grid%nxyz
    else
      ntot = nxy
    end if

    ! -- reading range over 1..ntot (contiguous, decided identically on all ranks)
    base = ntot/st_mpi%totn ; rem = mod(ntot, st_mpi%totn)
    rang_num = base
    if (st_mpi%rank < rem) then
      rang_num = base + 1
    end if
    rang_sta = st_mpi%rank*base + min(st_mpi%rank, rem) + 1

    ! -- step 1: count active elements per region in my range
    allocate(my_reg_ncals(totnreg))
    !$omp parallel do private(i)
    do i = 1, totnreg
      my_reg_ncals(i) = 0
    end do
    !$omp end parallel do
    do idx = 1, rang_num
      e = rang_sta + idx - 1
      reg_id = glob_reg_flag(e)
      if (reg_id >= 1) then
        my_reg_ncals(reg_id) = my_reg_ncals(reg_id) + 1
      end if
    end do

    ! -- step 2: global per-region counts (= div_calc_reg's reg_ncals)
    allocate(reg_ncals(totnreg))
    call mpisum_val(my_reg_ncals, "region column count", reg_ncals)

    ! -- step 3: reproduce the assignment rule (deterministic from reg_ncals)
    allocate(reg_remain(totnreg), reg_mpi_num(totnreg), reg_mpi_end(0:totnreg))
    reg_mpi_end(0) = 0
    !$omp parallel do private(i)
    do i = 1, totnreg
      reg_remain(i) = 0 ; reg_mpi_num(i) = 0 ; reg_mpi_end(i) = 0
    end do
    !$omp end parallel do
    rough_divn = sum(reg_ncals(:))/st_mpi%totn
    if (rough_divn == 0) then
      rough_divn = 1
    end if
    sum_cals = 0 ; sum_mpi = 0
    do i = 1, totnreg
      sum_cals = sum_cals + reg_ncals(i)
      quot = sum_cals/rough_divn ; reg_remain(i) = mod(sum_cals, rough_divn)
      if (quot == 0) then
        reg_mpi_num(i) = 1 ; reg_remain(i) = 0 ; reg_mpi_end(i) = sum_mpi + 1
      else if (sum_cals > reg_ncals(i) .and. i /= totnreg) then
        reg_mpi_end(i) = sum_mpi + 1 ; sum_cals = 0
      else
        reg_mpi_num(i) = quot ; sum_mpi = sum_mpi + quot ; reg_mpi_end(i) = sum_mpi
        sum_cals = 0
      end if
    end do
    pro_num = 0
    if (reg_mpi_end(totnreg) < st_mpi%totn) then
      do k = 1, st_mpi%totn-reg_mpi_end(totnreg)
        max_reg = maxloc(reg_remain(:),1)
        reg_mpi_num(max_reg) = reg_mpi_num(max_reg) + 1
        if (reg_remain(max_reg) < reg_ncals(max_reg)/2) then
          pro_num = max_reg + 1
        else
          pro_num = max_reg
        end if
        reg_remain(max_reg) = 0
        do i = pro_num, totnreg
          reg_mpi_end(i) = reg_mpi_end(i) + 1
        end do
      end do
    end if

    ! -- step 4: per-region counts held by earlier ranks (exclusive prefix sum)
    allocate(cals_end(totnreg))
    call mpiexscan_val(my_reg_ncals, "region column prefix", cals_end)

    ! -- step 5: owner rank of each element in my range
    allocate(rang_mpi(rang_num), seen(totnreg))
    !$omp parallel do private(i)
    do i = 1, rang_num
      rang_mpi(i) = 0
    end do
    !$omp end parallel do
    !$omp parallel do private(i)
    do i = 1, totnreg
      seen(i) = 0
    end do
    !$omp end parallel do
    do idx = 1, rang_num
      e = rang_sta + idx - 1
      reg_id = glob_reg_flag(e)
      if (reg_id >= 1) then
        seen(reg_id) = seen(reg_id) + 1
        pos_in_reg = cals_end(reg_id) + seen(reg_id)
        sta_mpi = reg_mpi_end(reg_id-1) ; end_mpi = reg_mpi_end(reg_id)
        if (end_mpi == sta_mpi) then
          sta_mpi = sta_mpi - 1
        end if
        reg_nproc = end_mpi - sta_mpi
        proc_cals = reg_ncals(reg_id)/reg_nproc
        remain_cals = mod(reg_ncals(reg_id), reg_nproc)
        if (pos_in_reg <= remain_cals*(proc_cals+1)) then
          proc_pos = (pos_in_reg-1)/(proc_cals+1) + 1
        else
          proc_pos = remain_cals + (pos_in_reg-1 - remain_cals*(proc_cals+1))/proc_cals + 1
        end if
        rang_mpi(idx) = sta_mpi + proc_pos
      end if
    end do

    ! -- step 6: gather owners into a full temp, then copy into the global flag
    allocate(glo_mpi(ntot))
    call gather_val(st_mpi%totn, rang_num, rang_mpi, glo_mpi, "mpi flag")
    !$omp parallel do private(i)
    do i = 1, ntot
      glob_mpi_flag(i) = glo_mpi(i)
    end do
    !$omp end parallel do
    deallocate(glo_mpi)

    ! -- 2d: replicate the column ownership to every layer
    if (.not. div_reg_3d) then
      !$omp parallel do private(i)
      do i = 1, st_grid%nz-1
        glob_mpi_flag(nxy*i+1:nxy*(i+1)) = glob_mpi_flag(1:nxy)
      end do
      !$omp end parallel do
    end if

    deallocate(my_reg_ncals, reg_ncals, reg_remain, reg_mpi_num, reg_mpi_end)
    deallocate(cals_end, seen, rang_mpi)

  end subroutine div_calc_reg_dist
#endif

  subroutine div_nocalc_flag()
  !*********************************************************************************************
  ! div_nocalc_flag -- Divide no calculation flag (serial; div_reg_3d: 2d per-layer / 3d whole)
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, k, nxy, nloop, nsize, nxyz0, nxyz1
    integer(I4) :: glo_nocals, rough_divn, nocals_remain
    integer(I4) :: nocals_sta, nocals_end
    integer(I4), allocatable :: nocals_mpi_num(:), grid_num(:), nocals_glo_num(:)
    integer(I4), allocatable :: nocals_mpi_glo(:)
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    if (div_reg_3d) then
      nloop = 1 ; nsize = st_grid%nxyz
    else
      nloop = st_grid%nz ; nsize = nxy
    end if
    do k = 1, nloop
      nxyz0 = nsize*(k-1)+1 ; nxyz1 = nsize*k
      glo_nocals = count(glob_mpi_flag(nxyz0:nxyz1) == 0)
      rough_divn = glo_nocals/st_mpi%totn ; nocals_remain = mod(glo_nocals, st_mpi%totn)
      allocate(nocals_mpi_num(st_mpi%totn), grid_num(nsize), nocals_glo_num(glo_nocals))
      !$omp parallel
      !$omp do private(i)
      do i = 1, st_mpi%totn
        nocals_mpi_num(i) = rough_divn
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, st_mpi%totn
        if (i <= nocals_remain) then
          nocals_mpi_num(i) = nocals_mpi_num(i) + 1
        end if
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, nsize
        grid_num(i) = nsize*(k-1) + i
      end do
      !$omp end do
      !$omp end parallel
      nocals_glo_num(:) = pack(grid_num(:), glob_mpi_flag(nxyz0:nxyz1) == 0)

      deallocate(grid_num)

      nocals_sta = 0 ; nocals_end = 0
      do i = 1, st_mpi%totn
        nocals_end = nocals_sta + nocals_mpi_num(i)
        allocate(nocals_mpi_glo(nocals_end-nocals_sta))
        nocals_mpi_glo(:) = nocals_glo_num(nocals_sta+1:nocals_end)
        glob_mpi_flag((/ nocals_mpi_glo /)) = -i
        nocals_sta = nocals_end
        deallocate(nocals_mpi_glo)
      end do

      deallocate(nocals_mpi_num, nocals_glo_num)

    end do

  end subroutine div_nocalc_flag

#ifdef MPI_MSG
  subroutine div_nocalc_dist()
  !*********************************************************************************************
  ! div_nocalc_dist -- Divide no calculation flag, distributed (reproduces div_nocalc_flag)
  !*********************************************************************************************
    ! -- modules
    use mpi_utility, only: mpisum_val, mpiexscan_val, gather_val
    ! -- inout

    ! -- local
    integer(I4) :: i, e, idx, nxy, ntot, base, rem, rang_sta, rang_num
    integer(I4) :: glo_nocals, rough_divn, nocals_remain, pos_nocal, proc_pos, seen
    integer(I4) :: my_nocal(1), all_nocal(1), nocal_end(1)
    integer(I4), allocatable :: rang_mpi(:), glo_mpi(:)
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    if (div_reg_3d) then
      ntot = st_grid%nxyz
    else
      ntot = nxy
    end if

    ! -- reading range over 1..ntot (contiguous, decided identically on all ranks)
    base = ntot/st_mpi%totn ; rem = mod(ntot, st_mpi%totn)
    rang_num = base
    if (st_mpi%rank < rem) then
      rang_num = base + 1
    end if
    rang_sta = st_mpi%rank*base + min(st_mpi%rank, rem) + 1

    ! -- step 1: count inactive cells (mpi flag == 0) in my range
    my_nocal(1) = 0
    do idx = 1, rang_num
      e = rang_sta + idx - 1
      if (glob_mpi_flag(e) == 0) then
        my_nocal(1) = my_nocal(1) + 1
      end if
    end do

    ! -- step 2: global inactive count, then round-robin sizes (deterministic)
    call mpisum_val(my_nocal, "nocalc count", all_nocal)
    glo_nocals = all_nocal(1)
    rough_divn = glo_nocals/st_mpi%totn
    nocals_remain = mod(glo_nocals, st_mpi%totn)

    ! -- step 4: inactive cells held by earlier ranks (exclusive prefix sum)
    call mpiexscan_val(my_nocal, "nocalc prefix", nocal_end)

    ! -- step 5: owner rank of each element (active keeps its owner; inactive gets -rank)
    allocate(rang_mpi(rang_num))
    seen = 0
    do idx = 1, rang_num
      e = rang_sta + idx - 1
      if (glob_mpi_flag(e) == 0) then
        seen = seen + 1
        pos_nocal = nocal_end(1) + seen
        if (pos_nocal <= nocals_remain*(rough_divn+1)) then
          proc_pos = (pos_nocal-1)/(rough_divn+1) + 1
        else
          proc_pos = nocals_remain + (pos_nocal-1 - nocals_remain*(rough_divn+1))/rough_divn + 1
        end if
        rang_mpi(idx) = -proc_pos
      else
        rang_mpi(idx) = glob_mpi_flag(e)
      end if
    end do

    ! -- step 6: gather into a full temp, then copy into the global flag
    allocate(glo_mpi(ntot))
    call gather_val(st_mpi%totn, rang_num, rang_mpi, glo_mpi, "nocalc flag")
    !$omp parallel do private(i)
    do i = 1, ntot
      glob_mpi_flag(i) = glo_mpi(i)
    end do
    !$omp end parallel do
    deallocate(glo_mpi)

    ! -- 2d: replicate the column flag to every layer
    if (.not. div_reg_3d) then
      !$omp parallel do private(i)
      do i = 1, st_grid%nz-1
        glob_mpi_flag(nxy*i+1:nxy*(i+1)) = glob_mpi_flag(1:nxy)
      end do
      !$omp end parallel do
    end if

    deallocate(rang_mpi)

  end subroutine div_nocalc_dist
#endif

  subroutine build_read_range()
  !*********************************************************************************************
  ! build_read_range -- Build this rank's read range read_reg/read_mpi (PRE-sea) over 1..nxyz
  !   via the unit->cell exchange. 2d: unit=column(nxy), key=mod(g-1,nxy)+1 ; 3d: unit=cell,
  !   key=g (handles z-varying regions). Sea is applied per-cell later (punch_sea_read_range).
  !   Runs before the sea puncture, where the global flags are z-uniform (2d). Scaffold slices
  !   the unit charge from the global flags; 段5-5 sources it from the distributed routines.
  !*********************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_utility, only: alltoall_val, alltoallv_val
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i, base_num, rem_num
#ifdef MPI_MSG
    integer(I4) :: g, k, key, nxy, unit_size, ubase, urem, unum, usta
    integer(I4) :: cand_raw, nkey, proc_id, dest_rank, send_tot, recv_tot, pos
    integer(I4), allocatable :: my_unit_reg(:), my_unit_mpi(:)
    integer(I4), allocatable :: cand_list(:), key_list(:), key_home(:)
    integer(I4), allocatable :: send_count(:), recv_count(:), fill_pos(:)
    integer(I4), allocatable :: q_key(:), recv_key(:), reply_reg(:), reply_owner(:)
    integer(I4), allocatable :: ans_reg(:), ans_owner(:)
    type(gmap_set) :: key_map
#endif
    !-------------------------------------------------------------------------------------------
    base_num = st_grid%nxyz/st_mpi%totn ; rem_num = mod(st_grid%nxyz, st_mpi%totn)
    read_num = base_num
    if (st_mpi%rank < rem_num) then
      read_num = base_num + 1
    end if
    read_sta = st_mpi%rank*base_num + min(st_mpi%rank, rem_num) + 1
    allocate(read_reg(read_num), read_mpi(read_num))

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      nxy = st_grid%nx*st_grid%ny
      if (div_reg_3d) then
        unit_size = st_grid%nxyz
      else
        unit_size = nxy
      end if
      ubase = unit_size/st_mpi%totn ; urem = mod(unit_size, st_mpi%totn)
      unum = ubase
      if (st_mpi%rank < urem) then
        unum = ubase + 1
      end if
      usta = st_mpi%rank*ubase + min(st_mpi%rank, urem) + 1
      allocate(my_unit_reg(unum), my_unit_mpi(unum))
      do i = 1, unum
        my_unit_reg(i) = glob_reg_flag(usta + i - 1)
        my_unit_mpi(i) = glob_mpi_flag(usta + i - 1)
      end do

      ! -- distinct keys needed by my read cells (2d: column, 3d: the cell itself)
      allocate(cand_list(read_num + 1))
      cand_raw = 0
      do i = 1, read_num
        g = read_sta + i - 1
        if (div_reg_3d) then
          key = g
        else
          key = mod(g-1, nxy) + 1
        end if
        cand_raw = cand_raw + 1 ; cand_list(cand_raw) = key
      end do
      if (cand_raw > 1) then
        call iquick_sort(cand_list, 1, cand_raw)
      end if
      allocate(key_list(cand_raw + 1))
      nkey = 0
      do k = 1, cand_raw
        if (k == 1) then
          nkey = 1 ; key_list(1) = cand_list(1)
        else if (cand_list(k) /= cand_list(k-1)) then
          nkey = nkey + 1 ; key_list(nkey) = cand_list(k)
        end if
      end do

      ! -- bucket keys to their home rank (inverse of the unit split), then exchange
      allocate(key_home(nkey + 1), send_count(st_mpi%totn), fill_pos(st_mpi%totn))
      do k = 1, nkey
        key = key_list(k)
        if (key <= urem*(ubase+1)) then
          key_home(k) = (key-1)/(ubase+1)
        else
          key_home(k) = urem + (key - 1 - urem*(ubase+1))/ubase
        end if
      end do
      do proc_id = 1, st_mpi%totn
        send_count(proc_id) = 0
      end do
      do k = 1, nkey
        dest_rank = key_home(k) + 1 ; send_count(dest_rank) = send_count(dest_rank) + 1
      end do
      send_tot = 0
      do proc_id = 1, st_mpi%totn
        fill_pos(proc_id) = send_tot ; send_tot = send_tot + send_count(proc_id)
      end do
      allocate(q_key(send_tot + 1))
      do k = 1, nkey
        dest_rank = key_home(k) + 1
        fill_pos(dest_rank) = fill_pos(dest_rank) + 1
        q_key(fill_pos(dest_rank)) = key_list(k)
      end do
      allocate(recv_count(st_mpi%totn))
      call alltoall_val(send_count, "read range count", recv_count)
      recv_tot = 0
      do proc_id = 1, st_mpi%totn
        recv_tot = recv_tot + recv_count(proc_id)
      end do
      allocate(recv_key(recv_tot + 1), reply_reg(recv_tot + 1), reply_owner(recv_tot + 1))
      call alltoallv_val(q_key, send_count, recv_count, "read range key", recv_key)
      do i = 1, recv_tot
        pos = recv_key(i) - usta + 1
        reply_reg(i) = my_unit_reg(pos) ; reply_owner(i) = my_unit_mpi(pos)
      end do
      allocate(ans_reg(send_tot + 1), ans_owner(send_tot + 1))
      call alltoallv_val(reply_reg, recv_count, send_count, "read range reg", ans_reg)
      call alltoallv_val(reply_owner, recv_count, send_count, "read range owner", ans_owner)

      ! -- key -> position in q_key, then fill each read cell from the exchanged value
      call gmap_init(key_map, send_tot)
      do k = 1, send_tot
        call gmap_put(key_map, q_key(k), k)
      end do
      do i = 1, read_num
        g = read_sta + i - 1
        if (div_reg_3d) then
          key = g
        else
          key = mod(g-1, nxy) + 1
        end if
        pos = gmap_get(key_map, key)
        read_reg(i) = ans_reg(pos) ; read_mpi(i) = ans_owner(pos)
      end do
      call gmap_free(key_map)
      deallocate(my_unit_reg, my_unit_mpi, cand_list, key_list, key_home)
      deallocate(send_count, fill_pos, recv_count, q_key, recv_key)
      deallocate(reply_reg, reply_owner, ans_reg, ans_owner)
    else
      do i = 1, read_num
        read_reg(i) = glob_reg_flag(read_sta + i - 1)
        read_mpi(i) = glob_mpi_flag(read_sta + i - 1)
      end do
    end if
#else
    do i = 1, read_num
      read_reg(i) = glob_reg_flag(read_sta + i - 1)
      read_mpi(i) = glob_mpi_flag(read_sta + i - 1)
    end do
#endif

  end subroutine build_read_range

  subroutine punch_sea_read_range(sea_num, sea_glo)
  !*********************************************************************************************
  ! punch_sea_read_range -- Zero read_reg/read_mpi for the sea cells listed in sea_glo that fall
  !   in this rank's read range (rank0-global fallback path for non-point sea formats).
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: sea_num
    integer(I4), intent(in) :: sea_glo(:)
    ! -- local
    integer(I4) :: i, g, li
    !-------------------------------------------------------------------------------------------
    do i = 1, sea_num
      g = sea_glo(i)
      if (g >= read_sta .and. g <= read_sta + read_num - 1) then
        li = g - read_sta + 1
        read_reg(li) = 0 ; read_mpi(li) = 0
      end if
    end do

  end subroutine punch_sea_read_range

  subroutine punch_sea_points(sea_ptn, sp_i, sp_j, sp_k)
  !*********************************************************************************************
  ! punch_sea_points -- Distributed sea determination: for each cell in this rank's read
  !   range, if it is a calc cell (read_reg>0) matched by any sea point spec (i,j,k with -1
  !   wildcards), zero read_reg/read_mpi. Per-cell, so z-varying sea is handled. Point specs are
  !   the same for all ranks (broadcast); no global array is used.
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: sea_ptn
    integer(I4), intent(in) :: sp_i(:), sp_j(:), sp_k(:)
    ! -- local
    integer(I4) :: li, g, p, nxy, i_num, j_num, k_num, ij
    logical :: is_sea
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    do li = 1, read_num
      if (read_reg(li) <= 0) then
        cycle
      end if
      g = read_sta + li - 1
      k_num = (g-1)/nxy + 1
      ij = g - nxy*(k_num-1)
      j_num = (ij-1)/st_grid%nx + 1
      i_num = ij - (j_num-1)*st_grid%nx
      is_sea = .false.
      do p = 1, sea_ptn
        if ((sp_i(p) == -1 .or. sp_i(p) == i_num) .and. &
            (sp_j(p) == -1 .or. sp_j(p) == j_num) .and. &
            (sp_k(p) == -1 .or. sp_k(p) == k_num)) then
          is_sea = .true.
        end if
      end do
      if (is_sea) then
        read_reg(li) = 0 ; read_mpi(li) = 0
      end if
    end do

  end subroutine punch_sea_points

  subroutine verify_read_range()
  !*********************************************************************************************
  ! verify_read_range -- Shadow-verify the distributed read range against the post-sea global
  !   slice (段5-2/5-3 scaffold, removed with the global flags in 段5-5).
  !*********************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_utility, only: mpimax_val
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i, g, mismatch, ndiag
#ifdef MPI_MSG
    integer(I4) :: mism_loc(1), mism_all(1)
#endif
    !-------------------------------------------------------------------------------------------
    mismatch = 0 ; ndiag = 0
    do i = 1, read_num
      g = read_sta + i - 1
      if (read_reg(i) /= glob_reg_flag(g) .or. read_mpi(i) /= glob_mpi_flag(g)) then
        mismatch = 1
        if (ndiag < 10) then
          write(*,'(a,6i9)') "S52-READ-MISS:", st_mpi%rank, g, glob_reg_flag(g),&
            glob_mpi_flag(g), read_reg(i), read_mpi(i)
          flush(6)
          ndiag = ndiag + 1
        end if
      end if
    end do
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      mism_loc(1) = mismatch
      call mpimax_val(mism_loc, "read range verification", mism_all)
      mismatch = mism_all(1)
    end if
#endif
    if (mismatch /= 0) then
      call write_err_stop("distributed read range vs global mismatch (see S52-READ-MISS).")
    end if

  end subroutine verify_read_range

  subroutine redist_owned2owener()
  !*********************************************************************************************
  ! redist_owned2owener -- Redistribute owned cells to their owner ranks
  !*********************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_utility, only: alltoall_val, alltoallv_val, mpimax_val
#endif
    ! -- inout

    ! -- local
    integer(I4) :: cell_id, ref_pos, mismatch
#ifdef MPI_MSG
    integer(I4) :: cell_sta, cell_num, base_num, rem_num, list_pos
    integer(I4) :: owner_code, dest_rank, proc_id, buf_pos
    integer(I4) :: send_tot, recv_tot, calc_cnt, nocal_cnt
    integer(I4) :: mism_loc(1), mism_all(1)
    integer(I4), allocatable :: send_count(:), recv_count(:), fill_pos(:)
    integer(I4), allocatable :: send_gidx(:), send_reg(:), recv_gidx(:), recv_reg(:)
#endif
    !-------------------------------------------------------------------------------------------
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      base_num = st_grid%nxyz/st_mpi%totn ; rem_num = mod(st_grid%nxyz, st_mpi%totn)
      cell_num = base_num
      if (st_mpi%rank < rem_num) then
        cell_num = base_num + 1
      end if
      cell_sta = st_mpi%rank*base_num + min(st_mpi%rank, rem_num) + 1

      allocate(send_count(st_mpi%totn))
      !$omp parallel do private(proc_id)
      do proc_id = 1, st_mpi%totn
        send_count(proc_id) = 0
      end do
      !$omp end parallel do
      do list_pos = 1, cell_num
        owner_code = read_mpi(list_pos)
        if (owner_code /= 0) then
          dest_rank = abs(owner_code)
          send_count(dest_rank) = send_count(dest_rank) + 1
        end if
      end do

      allocate(fill_pos(st_mpi%totn))
      send_tot = 0
      do proc_id = 1, st_mpi%totn
        fill_pos(proc_id) = send_tot
        send_tot = send_tot + send_count(proc_id)
      end do
      allocate(send_gidx(send_tot), send_reg(send_tot))
      do list_pos = 1, cell_num
        cell_id = cell_sta + list_pos - 1
        owner_code = read_mpi(list_pos)
        if (owner_code /= 0) then
          dest_rank = abs(owner_code)
          fill_pos(dest_rank) = fill_pos(dest_rank) + 1
          buf_pos = fill_pos(dest_rank)
          if (owner_code > 0) then
            send_gidx(buf_pos) = cell_id
          else
            send_gidx(buf_pos) = -cell_id
          end if
          send_reg(buf_pos) = read_reg(list_pos)
        end if
      end do

      allocate(recv_count(st_mpi%totn))
      call alltoall_val(send_count, "redist count", recv_count)
      recv_tot = 0
      do proc_id = 1, st_mpi%totn
        recv_tot = recv_tot + recv_count(proc_id)
      end do
      allocate(recv_gidx(recv_tot), recv_reg(recv_tot))
      call alltoallv_val(send_gidx, send_count, recv_count, "redist gidx", recv_gidx)
      call alltoallv_val(send_reg, send_count, recv_count, "redist reg", recv_reg)

      calc_cnt = 0 ; nocal_cnt = 0
      do list_pos = 1, recv_tot
        if (recv_gidx(list_pos) > 0) then
          calc_cnt = calc_cnt + 1
        else
          nocal_cnt = nocal_cnt + 1
        end if
      end do
      allocate(own_calc_glo(calc_cnt), own_calc_reg(calc_cnt), own_nocal_glo(nocal_cnt))
      calc_cnt = 0 ; nocal_cnt = 0
      do list_pos = 1, recv_tot
        if (recv_gidx(list_pos) > 0) then
          calc_cnt = calc_cnt + 1
          own_calc_glo(calc_cnt) = recv_gidx(list_pos)
          own_calc_reg(calc_cnt) = recv_reg(list_pos)
        else
          nocal_cnt = nocal_cnt + 1
          own_nocal_glo(nocal_cnt) = -recv_gidx(list_pos)
        end if
      end do

      if (calc_cnt > 1) then
        call iquick_sort2(own_calc_glo, own_calc_reg, 1, calc_cnt)
      end if
      if (nocal_cnt > 1) then
        call iquick_sort(own_nocal_glo, 1, nocal_cnt)
      end if

      deallocate(send_count, recv_count, fill_pos)
      deallocate(send_gidx, send_reg, recv_gidx, recv_reg)
    else
      call pack_owned_local()
    end if
#else
    call pack_owned_local()
#endif

    ! -- verification against the full-array scaffold (A3-a only)
    mismatch = 0 ; ref_pos = 0
    do cell_id = 1, st_grid%nxyz
      if (glob_mpi_flag(cell_id) == st_mpi%rank+1) then
        ref_pos = ref_pos + 1
        if (ref_pos > size(own_calc_glo)) then
          mismatch = 1
        else if (own_calc_glo(ref_pos) /= cell_id) then
          mismatch = 1
        else if (own_calc_reg(ref_pos) /= glob_reg_flag(cell_id)) then
          mismatch = 1
        end if
      end if
    end do
    if (ref_pos /= size(own_calc_glo)) then
      mismatch = 1
    end if
    ref_pos = 0
    do cell_id = 1, st_grid%nxyz
      if (glob_mpi_flag(cell_id) == -(st_mpi%rank+1)) then
        ref_pos = ref_pos + 1
        if (ref_pos > size(own_nocal_glo)) then
          mismatch = 1
        else if (own_nocal_glo(ref_pos) /= cell_id) then
          mismatch = 1
        end if
      end if
    end do
    if (ref_pos /= size(own_nocal_glo)) then
      mismatch = 1
    end if

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      mism_loc(1) = mismatch
      call mpimax_val(mism_loc, "redistribute verification", mism_all)
      mismatch = mism_all(1)
    end if
#endif
    if (mismatch /= 0) then
      call write_err_stop("redistribute owned verification.")
    end if

  end subroutine redist_owned2owener

  subroutine pack_owned_local()
  !*********************************************************************************************
  ! pack_owned_local -- Pack owned local info
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: cell_id, calc_cnt, nocal_cnt
    !-------------------------------------------------------------------------------------------
    calc_cnt = count(glob_mpi_flag(:) == st_mpi%rank+1)
    nocal_cnt = count(glob_mpi_flag(:) == -(st_mpi%rank+1))
    allocate(own_calc_glo(calc_cnt), own_calc_reg(calc_cnt), own_nocal_glo(nocal_cnt))
    calc_cnt = 0 ; nocal_cnt = 0
    do cell_id = 1, st_grid%nxyz
      if (glob_mpi_flag(cell_id) == st_mpi%rank+1) then
        calc_cnt = calc_cnt + 1
        own_calc_glo(calc_cnt) = cell_id
        own_calc_reg(calc_cnt) = glob_reg_flag(cell_id)
      else if (glob_mpi_flag(cell_id) == -(st_mpi%rank+1)) then
        nocal_cnt = nocal_cnt + 1
        own_nocal_glo(nocal_cnt) = cell_id
      end if
    end do

  end subroutine pack_owned_local

  subroutine set_rel_gloloc()
  !*********************************************************************************************
  ! set_rel_gloloc -- Set relationship global&local
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, own_pos, ij, nxy, count_reg, reg_flag, temp_nreg
    integer(I4) :: count_calc, count_cals
    integer(I4), allocatable :: temp_cend(:)
    logical, allocatable :: reg_seen(:)
    !-------------------------------------------------------------------------------------------
    allocate(reg_seen(totnreg))
    !$omp parallel do private(i)
    do i = 1, totnreg
      reg_seen(i) = .false.
    end do
    !$omp end parallel do
    count_reg = 0
    do own_pos = 1, ncalc
      reg_flag = own_calc_reg(own_pos)
      if (.not. reg_seen(reg_flag)) then
        reg_seen(reg_flag) = .true. ; count_reg = count_reg + 1
      end if
    end do
    loc_regn = count_reg

    deallocate(reg_seen)

    allocate(loc_nreg(loc_regn), temp_cend(0:loc_regn))
    allocate(st_conn%calc2reg(ncalc))
    temp_cend(0) = 0
    !$omp parallel
    !$omp do private(i)
    do i = 1, loc_regn
      loc_nreg(i) = 0 ; temp_cend(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncalc
      st_conn%calc2reg(i) = 0
    end do
    !$omp end do
    !$omp end parallel

    nxy = st_grid%nx*st_grid%ny
    count_reg = 0 ; count_calc = 0 ; count_cals = 0
    do own_pos = 1, ncalc
      ij = own_calc_glo(own_pos)
      reg_flag = 0
      if (count_reg /= 0) then
        do j = 1, count_reg
          if (loc_nreg(j) == own_calc_reg(own_pos)) then
            reg_flag = 1 ; temp_nreg = j
          end if
        end do
        if (reg_flag == 0) then
          count_reg = count_reg + 1 ; loc_nreg(count_reg) = own_calc_reg(own_pos)
          temp_nreg = count_reg
        end if
      else
        count_reg = 1 ; loc_nreg(1) = own_calc_reg(own_pos) ; temp_nreg = count_reg
      end if
      if (ij <= nxy) then
        count_cals = count_cals + 1
        l2g_ij(count_cals) = ij
      end if
      count_calc = count_calc + 1 ; temp_cend(temp_nreg) = temp_cend(temp_nreg) + 1
      l2g_ijk(count_calc) = ij ; st_conn%glo2loc_ijk(ij) = count_calc
      st_conn%calc2reg(count_calc) = temp_nreg
    end do

    allocate(calc_end(0:loc_regn))
    calc_end(0) = 0
    !$omp parallel do private(i)
    do i = 1, loc_regn
      calc_end(i) = calc_end(i-1) + temp_cend(i)
    end do
    !$omp end parallel do

    deallocate(temp_cend)

  end subroutine set_rel_gloloc

#ifdef MPI_MSG
  subroutine set_mpi_rel()
  !*********************************************************************************************
  ! set_mpi_rel -- Set mpi relationship
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, ij
    integer(I4) :: i_num, j_num, k_num, g_num, n_num, w_num, e_num, s_num, u_num, d_num
    integer(I4) :: neib_mpi_num, sta_send, end_send, sendrecv_num, send_loc, recv_loc
    integer(I4) :: cals_niebn, calc_niebn, reg_num, own_reg
    logical :: find_flag
    integer(I4), allocatable :: temp_neib_num(:), temp_neib_flag(:), temp_mpi_num(:)
    integer(I4), allocatable :: neib_glos(:), neib_gloc(:), neib_locc(:)
    integer(I4), allocatable :: sort_recv_num(:), temp_calc_reg(:)
    integer(I4), allocatable :: recv_num(:), send_num(:)
    integer(I4), allocatable :: temp_sort(:), sort_mpi_num(:), loc_send_num(:)
    integer(I4), allocatable :: mpi_l2g_ij(:), mpi_l2g_ijk(:), mpi_calc2reg(:)
    !-------------------------------------------------------------------------------------------
    if (st_mpi%totn /= 1) then
      allocate(neib_glos(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        neib_glos(i) = 0
      end do
      !$omp end parallel do
      do i = 1, ncals
        g_num = l2g_ij(i)
        own_reg = own_calc_reg(st_conn%glo2loc_ijk(g_num))
        k_num = (g_num-1)/(st_grid%nx*st_grid%ny) + 1
        ij = g_num - st_grid%nx*st_grid%ny*(k_num-1)
        j_num = (ij-1)/st_grid%nx + 1
        i_num = ij - (j_num-1)*st_grid%nx
        ! north direction
        if (j_num /= 1) then
          n_num = g_num-st_grid%nx
          call get_hash_info(n_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncals = neib_ncals + 1 ; neib_glos(neib_ncals) = n_num
            end if
          end if
        end if
        ! west direction
        if (i_num /= 1) then
          w_num = g_num-1
          call get_hash_info(w_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncals = neib_ncals + 1 ; neib_glos(neib_ncals) = w_num
            end if
          end if
        end if
        ! east direction
        if (i_num /= st_grid%nx) then
          e_num = g_num+1
          call get_hash_info(e_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncals = neib_ncals + 1 ; neib_glos(neib_ncals) = e_num
            end if
          end if
        end if
        ! south direction
        if (j_num /= st_grid%ny) then
          s_num = g_num+st_grid%nx
          call get_hash_info(s_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncals = neib_ncals + 1 ; neib_glos(neib_ncals) = s_num
            end if
          end if
        end if
      end do

      allocate(temp_neib_num(st_mpi%totn), temp_neib_flag(st_mpi%totn))
      allocate(temp_mpi_num(ncalc), neib_locc(ncalc))
      allocate(neib_gloc(ncalc), temp_calc_reg(ncalc))
      !$omp parallel
      !$omp do private(i)
      do i = 1, st_mpi%totn
        temp_neib_num(i) = 0 ; temp_neib_flag(i) = 0
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, ncalc
        temp_mpi_num(i) = 0 ; neib_locc(i) = 0 ; neib_gloc(i) = 0 ; temp_calc_reg(i) = 0
      end do
      !$omp end do
      !$omp end parallel

      do i = 1, ncalc
        g_num = l2g_ijk(i)
        own_reg = own_calc_reg(st_conn%glo2loc_ijk(g_num))
        k_num = (g_num-1)/(st_grid%nx*st_grid%ny) + 1
        ij = g_num - st_grid%nx*st_grid%ny*(k_num-1)
        j_num = (ij-1)/st_grid%nx + 1
        i_num = ij - (j_num-1)*st_grid%nx
        ! up direction
        if (k_num /= 1) then
          u_num = g_num-st_grid%nx*st_grid%ny
          call get_hash_info(u_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1
              temp_calc_reg(neib_ncalc) = reg_num ; temp_mpi_num(neib_ncalc) = neib_mpi_num
              neib_locc(neib_ncalc) = i ; neib_gloc(neib_ncalc) = u_num
              if (temp_neib_flag(neib_mpi_num) == 0) then
                neib_mpi_totn = neib_mpi_totn + 1 ; temp_neib_flag(neib_mpi_num) = 1
                temp_neib_num(neib_mpi_totn) = neib_mpi_num
              end if
            end if
          end if
        end if
        ! north direction
        if (j_num /= 1) then
          n_num = g_num-st_grid%nx
          call get_hash_info(n_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1
              temp_calc_reg(neib_ncalc) = reg_num ; temp_mpi_num(neib_ncalc) = neib_mpi_num
              neib_locc(neib_ncalc) = i ; neib_gloc(neib_ncalc) = n_num
              if (temp_neib_flag(neib_mpi_num) == 0) then
                neib_mpi_totn = neib_mpi_totn + 1 ; temp_neib_flag(neib_mpi_num) = 1
                temp_neib_num(neib_mpi_totn) = neib_mpi_num
              end if
            end if
          end if
        end if
        ! west direction
        if (i_num /= 1) then
          w_num = g_num-1
          call get_hash_info(w_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1
              temp_calc_reg(neib_ncalc) = reg_num ; temp_mpi_num(neib_ncalc) = neib_mpi_num
              neib_locc(neib_ncalc) = i ; neib_gloc(neib_ncalc) = w_num
              if (temp_neib_flag(neib_mpi_num) == 0) then
                neib_mpi_totn = neib_mpi_totn + 1 ; temp_neib_flag(neib_mpi_num) = 1
                temp_neib_num(neib_mpi_totn) = neib_mpi_num
              end if
            end if
          end if
        end if
        ! east direction
        if (i_num /= st_grid%nx) then
          e_num = g_num+1
          call get_hash_info(e_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1
              temp_calc_reg(neib_ncalc) = reg_num ; temp_mpi_num(neib_ncalc) = neib_mpi_num
              neib_locc(neib_ncalc) = i ; neib_gloc(neib_ncalc) = e_num
              if (temp_neib_flag(neib_mpi_num) == 0) then
                neib_mpi_totn = neib_mpi_totn + 1 ; temp_neib_flag(neib_mpi_num) = 1
                temp_neib_num(neib_mpi_totn) = neib_mpi_num
              end if
            end if
          end if
        end if
        ! south direction
        if (j_num /= st_grid%ny) then
          s_num = g_num+st_grid%nx
          call get_hash_info(s_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1
              temp_calc_reg(neib_ncalc) = reg_num ; temp_mpi_num(neib_ncalc) = neib_mpi_num
              neib_locc(neib_ncalc) = i ; neib_gloc(neib_ncalc) = s_num
              if (temp_neib_flag(neib_mpi_num) == 0) then
                neib_mpi_totn = neib_mpi_totn + 1 ; temp_neib_flag(neib_mpi_num) = 1
                temp_neib_num(neib_mpi_totn) = neib_mpi_num
              end if
            end if
          end if
        end if
        ! down direction
        if (k_num /= st_grid%nz) then
          d_num = g_num+st_grid%nx*st_grid%ny
          call get_hash_info(d_num, find_flag, neib_mpi_num, reg_num)
          if (reg_num > 0 .and. neib_mpi_num /= st_mpi%rank+1) then
            if (reg_num == own_reg .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1
              temp_calc_reg(neib_ncalc) = reg_num ; temp_mpi_num(neib_ncalc) = neib_mpi_num
              neib_locc(neib_ncalc) = i ; neib_gloc(neib_ncalc) = d_num
              if (temp_neib_flag(neib_mpi_num) == 0) then
                neib_mpi_totn = neib_mpi_totn + 1 ; temp_neib_flag(neib_mpi_num) = 1
                temp_neib_num(neib_mpi_totn) = neib_mpi_num
              end if
            end if
          end if
        end if
      end do

      deallocate(temp_neib_flag)

      allocate(mpi_l2g_ij(ncals+neib_ncals), mpi_l2g_ijk(ncalc+neib_ncalc))
      allocate(mpi_calc2reg(ncalc+neib_ncalc))
      !$omp parallel
      !$omp do private(i)
      do i = 1, ncals
        mpi_l2g_ij(i) = l2g_ij(i)
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, ncalc
        mpi_l2g_ijk(i) = l2g_ijk(i)
        mpi_calc2reg(i) = st_conn%calc2reg(i)
      end do
      !$omp end do
      !$omp end parallel
      deallocate(l2g_ij, l2g_ijk, st_conn%calc2reg)

      !$omp parallel
      !$omp do private(i, cals_niebn)
      do i = 1, neib_ncals
        cals_niebn = ncals + i
        mpi_l2g_ij(cals_niebn) = neib_glos(i)
      end do
      !$omp end do
      !$omp do private(i, calc_niebn)
      do i = 1, neib_ncalc
        calc_niebn = ncalc + i
        mpi_l2g_ijk(calc_niebn) = neib_gloc(i)
        mpi_calc2reg(calc_niebn) = temp_calc_reg(i)
        st_conn%glo2loc_ijk(neib_gloc(i)) = calc_niebn
      end do
      !$omp end do
      !$omp end parallel

      deallocate(neib_glos, temp_calc_reg)

      allocate(l2g_ij(ncals+neib_ncals), l2g_ijk(ncalc+neib_ncalc))
      allocate(st_conn%calc2reg(ncalc+neib_ncalc))
      !$omp parallel
      !$omp do private(i)
      do i = 1, ncals+neib_ncals
        l2g_ij(i) = mpi_l2g_ij(i)
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, ncalc+neib_ncalc
        l2g_ijk(i) = mpi_l2g_ijk(i)
        st_conn%calc2reg(i) = mpi_calc2reg(i)
      end do
      !$omp end do
      !$omp end parallel

      deallocate(mpi_l2g_ij, mpi_l2g_ijk, mpi_calc2reg)

      sendrecv_num = count(temp_mpi_num(:) /= 0)
      allocate(neib_num(neib_mpi_totn))
      allocate(send_cind(0:neib_mpi_totn), recv_cind(0:neib_mpi_totn))
      allocate(send_citem(sendrecv_num), recv_citem(sendrecv_num))
      allocate(temp_sort(neib_ncalc), sort_recv_num(neib_ncalc), sort_mpi_num(neib_ncalc))
      allocate(send2recv(neib_ncalc), loc_send_num(neib_ncalc))
      send_cind(0) = 0 ; recv_cind(0) = 0
      !$omp parallel
      !$omp do private(i)
      do i = 1, neib_mpi_totn
        send_cind(i) = 0 ; recv_cind(i) = 0
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, sendrecv_num
        send_citem(i) = 0 ; recv_citem(i) = 0
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, neib_ncalc
        temp_sort(i) = neib_gloc(i)
        sort_recv_num(i) = 0 ; send2recv(i) = 0
        loc_send_num(i) = neib_locc(i)
      end do
      !$omp end do
      !$omp end parallel

      call iquick_sort(temp_sort, 1, neib_ncalc)
      !$omp parallel
      !$omp do private(i, recv_loc)
      do i = 1, neib_ncalc
        recv_loc = findloc(neib_gloc(1:neib_ncalc), value = temp_sort(i), dim = 1)
        neib_gloc(recv_loc) = 0
        sort_recv_num(i) = neib_locc(recv_loc) ; sort_mpi_num(i) = temp_mpi_num(recv_loc)
      end do
      !$omp end do
      !$omp end parallel
      deallocate(neib_gloc)

      do i = 1, neib_mpi_totn
        neib_mpi_num = temp_neib_num(i) ; neib_num(i) = neib_mpi_num - 1
        send_cind(i) = count(sort_mpi_num(:) == neib_mpi_num) + send_cind(i-1)
        recv_cind(i) = send_cind(i)
        sta_send = send_cind(i-1) ; end_send = send_cind(i)
        allocate(recv_num(count(sort_mpi_num(:) == neib_mpi_num)))
        allocate(send_num(count(sort_mpi_num(:) == neib_mpi_num)))
        recv_num(:) = pack(sort_recv_num(:), sort_mpi_num(:) == neib_mpi_num)
        send_num(:) = pack(neib_locc(:), temp_mpi_num(:) == neib_mpi_num)
        do j = 1, end_send-sta_send
          k = sta_send + j ; recv_citem(k) = ncalc + k
          send_loc = findloc(loc_send_num(:), value = recv_num(j), dim = 1)
          send2recv(send_loc) = k ; send_citem(k) = send_num(j)
          loc_send_num(send_loc) = 0
          if (calc2recv(recv_num(j)) == 0) then
            calc2recv(recv_num(j)) = ncalc + k
          end if
        end do
        deallocate(recv_num, send_num)
      end do

      deallocate(temp_neib_num, temp_mpi_num, sort_recv_num, sort_mpi_num)
      deallocate(neib_locc, loc_send_num)

    end if

  end subroutine set_mpi_rel
#endif

  subroutine set_rel_seareg()
  !*********************************************************************************************
  ! set_rel_seareg -- Set relationship of sea region
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, ij, nxy
    integer(I4) :: i_num, j_num, k_num, g_num, n_num, w_num, e_num, s_num, u_num, d_num
    integer(I4) :: cand_num, surf_num, sdir_num, neib_snum(4)
    !-------------------------------------------------------------------------------------------
    seal_snum = 0
    nxy = st_grid%nx*st_grid%ny

    ! -- (A3-d2) pass 1: count seal candidates so the surface map is sized by the actual amount
    !    (duplicates are counted too, which only makes the table slightly larger = safe).
    !    A provisional map holding own + halo is needed to answer "already mapped?" here.
    surf_num = ncals + neib_ncals
    call gmap_init(surf_map, surf_num)
    do i = 1, surf_num
      call gmap_put(surf_map, l2g_ij(i), i)
    end do

    cand_num = 0
    do i = 1, ncals
      call get_neib_snum(l2g_ij(i), neib_snum, sdir_num)
      do j = 1, sdir_num
        if (glob_reg_flag(neib_snum(j)) == 0 .and. gmap_get(surf_map, neib_snum(j)) == 0) then
          cand_num = cand_num + 1
        end if
      end do
    end do

    ! -- rebuild at the final capacity, then register own + halo again
    call gmap_init(surf_map, surf_num + cand_num)
    do i = 1, surf_num
      call gmap_put(surf_map, l2g_ij(i), i)
    end do

    do i = 1, ncals
      g_num = l2g_ij(i)
      k_num = (g_num-1)/(nxy) + 1
      ij = g_num - nxy*(k_num-1)
      j_num = (ij-1)/st_grid%nx + 1
      i_num = ij - (j_num-1)*st_grid%nx
      ! north direction
      if (j_num /= 1) then
        n_num = g_num-st_grid%nx
        if (glob_reg_flag(n_num) == 0 .and. gmap_get(surf_map, n_num) == 0) then
          seal_snum = seal_snum + 1
          call gmap_put(surf_map, n_num, ncals + neib_ncals + seal_snum)
        end if
      end if
      ! west direction
      if (i_num /= 1) then
        w_num = g_num-1
        if (glob_reg_flag(w_num) == 0 .and. gmap_get(surf_map, w_num) == 0) then
          seal_snum = seal_snum + 1
          call gmap_put(surf_map, w_num, ncals + neib_ncals + seal_snum)
        end if
      end if
      ! east direction
      if (i_num /= st_grid%nx) then
        e_num = g_num+1
        if (glob_reg_flag(e_num) == 0 .and. gmap_get(surf_map, e_num) == 0) then
          seal_snum = seal_snum + 1
          call gmap_put(surf_map, e_num, ncals + neib_ncals + seal_snum)
        end if
      end if
      ! south direction
      if (j_num /= st_grid%ny) then
        s_num = g_num+st_grid%nx
        if (glob_reg_flag(s_num) == 0 .and. gmap_get(surf_map, s_num) == 0) then
          seal_snum = seal_snum + 1
          call gmap_put(surf_map, s_num, ncals + neib_ncals + seal_snum)
        end if
      end if
    end do

    seal_cnum = 0

    do i = 1, ncalc
      g_num = l2g_ijk(i)
      k_num = (g_num-1)/(nxy) + 1
      ij = g_num - nxy*(k_num-1)
      j_num = (ij-1)/st_grid%nx + 1
      i_num = ij - (j_num-1)*st_grid%nx
      ! up direction
      if (k_num /= 1) then
        u_num = g_num-nxy
        if (glob_reg_flag(u_num) == 0 .and. st_conn%glo2loc_ijk(u_num) == 0) then
          seal_cnum = seal_cnum + 1
          st_conn%glo2loc_ijk(u_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! north direction
      if (j_num /= 1) then
        n_num = g_num-st_grid%nx
        if (glob_reg_flag(n_num) == 0 .and. st_conn%glo2loc_ijk(n_num) == 0) then
          seal_cnum = seal_cnum + 1
          st_conn%glo2loc_ijk(n_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! west direction
      if (i_num /= 1) then
        w_num = g_num-1
        if (glob_reg_flag(w_num) == 0 .and. st_conn%glo2loc_ijk(w_num) == 0) then
          seal_cnum = seal_cnum + 1
          st_conn%glo2loc_ijk(w_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! east direction
      if (i_num /= st_grid%nx) then
        e_num = g_num+1
        if (glob_reg_flag(e_num) == 0 .and. st_conn%glo2loc_ijk(e_num) == 0) then
          seal_cnum = seal_cnum + 1
          st_conn%glo2loc_ijk(e_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! south direction
      if (j_num /= st_grid%ny) then
        s_num = g_num+st_grid%nx
        if (glob_reg_flag(s_num) == 0 .and. st_conn%glo2loc_ijk(s_num) == 0) then
          seal_cnum = seal_cnum + 1
          st_conn%glo2loc_ijk(s_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! down direction
      if (k_num /= st_grid%nz) then
        d_num = g_num+nxy
        if (glob_reg_flag(d_num) == 0 .and. st_conn%glo2loc_ijk(d_num) == 0) then
          seal_cnum = seal_cnum + 1
          st_conn%glo2loc_ijk(d_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
    end do

  end subroutine set_rel_seareg

  subroutine set_loc_cell_clas()
  !*********************************************************************************************
  ! set_loc_cell_clas -- Set local cell classification. (段5-4) When classes exist, each cell's
  !   flags are built from the class specs (i,j,k with -1 wildcards). The noclas case leaves the
  !   flags zero (unused downstream). Class specs are freed here after use.
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_ctrl
    ! -- inout

    ! -- local
    integer(I4) :: i, j, p, g, nxy, ij, ii, jj, kk, ci, cj, ck, cflag
    !-------------------------------------------------------------------------------------------
    allocate(st_conn%clas_flag(ncell,st_clas%totn))
    !$omp parallel do private(i)
    do i = 1, ncell
      st_conn%clas_flag(i,:) = 0
    end do
    !$omp end parallel do

    if (st_ctrl%noclas_flag == 1) then
      ! -- no classification: specs are not available; clas_flag stays zero (unused downstream)
      return
    end if

    nxy = st_grid%nx*st_grid%ny
    !$omp parallel do private(i, j, p, g, ij, ii, jj, kk, ci, cj, ck, cflag)
    do i = 1, ncell
      g = st_conn%loc2glo_ijk(i)
      kk = (g-1)/nxy + 1
      ij = g - nxy*(kk-1)
      jj = (ij-1)/st_grid%nx + 1
      ii = ij - (jj-1)*st_grid%nx
      do j = 1, st_clas%totn
        cflag = 0
        do p = 1, st_clas%num(j)
          ci = st_clas%i(p,j) ; cj = st_clas%j(p,j) ; ck = st_clas%k(p,j)
          if ((ci == -1 .or. ci == ii) .and. (cj == -1 .or. cj == jj) .and. &
              (ck == -1 .or. ck == kk)) then
            cflag = 1
          end if
        end do
        st_conn%clas_flag(i,j) = cflag
      end do
    end do
    !$omp end parallel do

    deallocate(st_clas%i, st_clas%j, st_clas%k, st_clas%num)

  end subroutine set_loc_cell_clas

  subroutine get_cals_grid(cal_num, x_num, y_num)
  !*********************************************************************************************
  ! get_cals_grid -- Get surface number from grid number
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: cal_num
    integer(I4), intent(out) :: x_num, y_num
    ! -- local
    integer(I4) :: s_num
    !-------------------------------------------------------------------------------------------
    s_num = st_conn%loc2glo_ij(cal_num)
    y_num = (s_num-1)/st_grid%nx + 1
    x_num = s_num - (y_num-1)*st_grid%nx

  end subroutine get_cals_grid

  subroutine get_calc_grid(cal_num, x_num, y_num, z_num)
  !*********************************************************************************************
  ! get_calc_grid -- Get calculation number from grid number
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: cal_num
    integer(I4), intent(out) :: x_num, y_num, z_num
    ! -- local
    integer(I4) :: c_num, xy_num
    !-------------------------------------------------------------------------------------------
    c_num = st_conn%loc2glo_ijk(cal_num)
    z_num = (c_num-1)/(st_grid%nx*st_grid%ny) + 1
    xy_num = c_num - st_grid%nx*st_grid%ny*(z_num-1)
    y_num = (xy_num-1)/st_grid%nx + 1
    x_num = xy_num - (y_num-1)*st_grid%nx

  end subroutine get_calc_grid

  subroutine get_neib_gnum(g_num, neib_gnum, dir_num)
  !*********************************************************************************************
  ! get_neib_gnum -- Get in-bounds 6-neighbor global numbers of a cell
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: g_num
    integer(I4), intent(out) :: neib_gnum(6), dir_num
    ! -- local
    integer(I4) :: i_num, j_num, k_num, ij, nxy
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    k_num = (g_num-1)/nxy + 1
    ij = g_num - nxy*(k_num-1)
    j_num = (ij-1)/st_grid%nx + 1
    i_num = ij - (j_num-1)*st_grid%nx
    dir_num = 0
    if (k_num /= 1) then
      dir_num = dir_num + 1 ; neib_gnum(dir_num) = g_num-nxy
    end if
    if (j_num /= 1) then
      dir_num = dir_num + 1 ; neib_gnum(dir_num) = g_num-st_grid%nx
    end if
    if (i_num /= 1) then
      dir_num = dir_num + 1 ; neib_gnum(dir_num) = g_num-1
    end if
    if (i_num /= st_grid%nx) then
      dir_num = dir_num + 1 ; neib_gnum(dir_num) = g_num+1
    end if
    if (j_num /= st_grid%ny) then
      dir_num = dir_num + 1 ; neib_gnum(dir_num) = g_num+st_grid%nx
    end if
    if (k_num /= st_grid%nz) then
      dir_num = dir_num + 1 ; neib_gnum(dir_num) = g_num+nxy
    end if

  end subroutine get_neib_gnum

  subroutine add_orphan_nocal()
  !*********************************************************************************************
  ! add_orphan_nocal -- Punctured sea cells with no active neighbour (orphans) belong to no
  !   calc/nocalc/seal list, so the mpi view-write would leave them 0 while the serial full-grid
  !   write fills SNOVAL. Append them to the nocalc list (output-only SNOVAL, no computational
  !   role) so both outputs match. Each rank claims the orphans in its read range; adjacency is
  !   judged with the global flag. (段5-5: replace glob_mpi_flag with distributed range owners.)
  !*********************************************************************************************
    ! -- modules

    ! -- local
    integer(I4) :: g, n, j, norphan, dir_num, neib_gnum(6)
    integer(I4), allocatable :: temp_nocal(:)
    logical :: is_orphan
    !-------------------------------------------------------------------------------------------
    norphan = 0
    do g = read_sta, read_sta + read_num - 1
      if (glob_mpi_flag(g) /= 0) then
        cycle
      end if
      call get_neib_gnum(g, neib_gnum, dir_num)
      is_orphan = .true.
      do n = 1, dir_num
        if (glob_mpi_flag(neib_gnum(n)) > 0) then
          is_orphan = .false. ; exit
        end if
      end do
      if (is_orphan) then
        norphan = norphan + 1
      end if
    end do

    if (norphan == 0) then
      return
    end if

    allocate(temp_nocal(size(own_nocal_glo) + norphan))
    do j = 1, size(own_nocal_glo)
      temp_nocal(j) = own_nocal_glo(j)
    end do
    j = size(own_nocal_glo)
    do g = read_sta, read_sta + read_num - 1
      if (glob_mpi_flag(g) /= 0) then
        cycle
      end if
      call get_neib_gnum(g, neib_gnum, dir_num)
      is_orphan = .true.
      do n = 1, dir_num
        if (glob_mpi_flag(neib_gnum(n)) > 0) then
          is_orphan = .false. ; exit
        end if
      end do
      if (is_orphan) then
        j = j + 1 ; temp_nocal(j) = g
      end if
    end do
    call move_alloc(temp_nocal, own_nocal_glo)
    call iquick_sort(own_nocal_glo, 1, size(own_nocal_glo))

  end subroutine add_orphan_nocal

  subroutine get_neib_snum(g_num, neib_snum, sdir_num)
  !*********************************************************************************************
  ! get_neib_snum -- Get in-bounds 4-neighbor surface global numbers (n, w, e, s order)
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: g_num
    integer(I4), intent(out) :: neib_snum(4), sdir_num
    ! -- local
    integer(I4) :: i_num, j_num, k_num, ij, nxy
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    k_num = (g_num-1)/nxy + 1
    ij = g_num - nxy*(k_num-1)
    j_num = (ij-1)/st_grid%nx + 1
    i_num = ij - (j_num-1)*st_grid%nx
    sdir_num = 0
    if (j_num /= 1) then
      sdir_num = sdir_num + 1 ; neib_snum(sdir_num) = g_num-st_grid%nx
    end if
    if (i_num /= 1) then
      sdir_num = sdir_num + 1 ; neib_snum(sdir_num) = g_num-1
    end if
    if (i_num /= st_grid%nx) then
      sdir_num = sdir_num + 1 ; neib_snum(sdir_num) = g_num+1
    end if
    if (j_num /= st_grid%ny) then
      sdir_num = sdir_num + 1 ; neib_snum(sdir_num) = g_num+st_grid%nx
    end if

  end subroutine get_neib_snum

#ifdef MPI_MSG
  subroutine set_dist_unknum(own_gnum, own_num, gnum_range, unk_num)
  !*********************************************************************************************
  ! set_dist_unknum -- Rank of each own global number within the ascending order of every rank's
  !                    own cells, obtained without gathering the whole numbering (A3-f).
  !                    Buckets split [1..gnum_range] into rank-ordered intervals, so the global
  !                    ascending order is the buckets concatenated in rank order; the rank of a
  !                    cell is therefore (entries owned by lower buckets) + (position inside).
  !*********************************************************************************************
    ! -- modules
    use mpi_utility, only: alltoall_val, alltoallv_val, mpiexscan_val

    ! -- inout
    integer(I4), intent(in) :: own_num, gnum_range
    integer(I4), intent(in) :: own_gnum(own_num)
    integer(I4), intent(out) :: unk_num(own_num)
    ! -- local
    integer(I4) :: i, proc_id, dest_rank, send_tot, recv_tot, list_pos, sort_pos
    integer(I4) :: base_num, rem_num, gnum
    integer(I4) :: my_count(1), off_count(1)
    integer(I4), allocatable :: send_count(:), recv_count(:), fill_pos(:)
    integer(I4), allocatable :: q_gnum(:), q_orig(:), recv_gnum(:), recv_orig(:)
    integer(I4), allocatable :: reply_num(:), ans_num(:)
    !-------------------------------------------------------------------------------------------
    if (st_mpi%totn == 1) then
      do i = 1, own_num
        unk_num(i) = i
      end do
      return
    end if

    base_num = gnum_range/st_mpi%totn ; rem_num = mod(gnum_range, st_mpi%totn)

    ! -- bucket the own numbers by the rank in charge of their interval
    allocate(send_count(st_mpi%totn), recv_count(st_mpi%totn), fill_pos(st_mpi%totn))
    do proc_id = 1, st_mpi%totn
      send_count(proc_id) = 0
    end do
    do i = 1, own_num
      gnum = own_gnum(i)
      if (gnum <= rem_num*(base_num+1)) then
        dest_rank = (gnum-1)/(base_num+1) + 1
      else
        dest_rank = rem_num + (gnum - 1 - rem_num*(base_num+1))/base_num + 1
      end if
      send_count(dest_rank) = send_count(dest_rank) + 1
    end do
    send_tot = 0
    do proc_id = 1, st_mpi%totn
      fill_pos(proc_id) = send_tot ; send_tot = send_tot + send_count(proc_id)
    end do

    allocate(q_gnum(send_tot+1), q_orig(send_tot+1))
    do i = 1, own_num
      gnum = own_gnum(i)
      if (gnum <= rem_num*(base_num+1)) then
        dest_rank = (gnum-1)/(base_num+1) + 1
      else
        dest_rank = rem_num + (gnum - 1 - rem_num*(base_num+1))/base_num + 1
      end if
      fill_pos(dest_rank) = fill_pos(dest_rank) + 1
      q_gnum(fill_pos(dest_rank)) = gnum ; q_orig(fill_pos(dest_rank)) = i
    end do

    call alltoall_val(send_count, "unknown count", recv_count)
    recv_tot = 0
    do proc_id = 1, st_mpi%totn
      recv_tot = recv_tot + recv_count(proc_id)
    end do

    allocate(recv_gnum(recv_tot+1), recv_orig(recv_tot+1), reply_num(recv_tot+1))
    call alltoallv_val(q_gnum, send_count, recv_count, "unknown gnum", recv_gnum)

    ! -- order this bucket, remembering where each entry came from
    do list_pos = 1, recv_tot
      recv_orig(list_pos) = list_pos
    end do
    if (recv_tot > 1) then
      call iquick_sort2(recv_gnum, recv_orig, 1, recv_tot)
    end if

    ! -- entries held by lower-numbered buckets (MPI_EXSCAN leaves rank 0 undefined)
    my_count(1) = recv_tot
    call mpiexscan_val(my_count, "unknown offset", off_count)
    if (st_mpi%rank == 0) then
      off_count(1) = 0
    end if

    do sort_pos = 1, recv_tot
      reply_num(recv_orig(sort_pos)) = off_count(1) + sort_pos
    end do

    allocate(ans_num(send_tot+1))
    call alltoallv_val(reply_num, recv_count, send_count, "unknown number", ans_num)
    do list_pos = 1, send_tot
      unk_num(q_orig(list_pos)) = ans_num(list_pos)
    end do

    deallocate(send_count, recv_count, fill_pos, q_gnum, q_orig)
    deallocate(recv_gnum, recv_orig, reply_num, ans_num)

  end subroutine set_dist_unknum

#endif

  subroutine set_hash_neib()
  !*********************************************************************************************
  ! set_hash_neib -- Set hash neighbor
  !*********************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_utility, only: alltoall_val, alltoallv_val, mpimax_val
#endif
    ! -- inout

    ! -- local
    integer(I4) :: own_pos, dir_pos, dir_num, neib_gnum(6), nb_gnum
    integer(I4) :: k, cand_raw, nhalo, active_num, mismatch, bad, ndiag, rf
    integer(I4) :: base_num, rem_num, hmpi, hreg
    logical :: find_flag
    integer(I4), allocatable :: cand_list(:), halo_gnum(:)
    integer(I4), allocatable :: q_gnum(:), ans_reg(:), ans_owner(:)
#ifdef MPI_MSG
    integer(I4) :: gnum, local_idx, proc_id, dest_rank, send_tot, recv_tot, query_pos
    integer(I4) :: mism_loc(1), mism_all(1)
    integer(I4), allocatable :: halo_home(:), send_count(:), recv_count(:), fill_pos(:)
    integer(I4), allocatable :: recv_gnum(:), reply_reg(:), reply_owner(:)
#endif
    !-------------------------------------------------------------------------------------------
    base_num = st_grid%nxyz/st_mpi%totn ; rem_num = mod(st_grid%nxyz, st_mpi%totn)

    allocate(cand_list(6*ncalc + 1))
    cand_raw = 0
    do own_pos = 1, ncalc
      call get_neib_gnum(l2g_ijk(own_pos), neib_gnum, dir_num)
      do dir_pos = 1, dir_num
        nb_gnum = neib_gnum(dir_pos)
        if (st_conn%glo2loc_ijk(nb_gnum) == 0) then
          cand_raw = cand_raw + 1 ; cand_list(cand_raw) = nb_gnum
        end if
      end do
    end do
    if (cand_raw > 1) then
      call iquick_sort(cand_list, 1, cand_raw)
    end if
    allocate(halo_gnum(cand_raw + 1))
    nhalo = 0
    do k = 1, cand_raw
      if (k == 1) then
        nhalo = 1 ; halo_gnum(1) = cand_list(1)
      else if (cand_list(k) /= cand_list(k-1)) then
        nhalo = nhalo + 1 ; halo_gnum(nhalo) = cand_list(k)
      end if
    end do

    allocate(q_gnum(nhalo + 1), ans_reg(nhalo + 1), ans_owner(nhalo + 1))
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      allocate(halo_home(nhalo + 1))
      do k = 1, nhalo
        gnum = halo_gnum(k)
        if (gnum <= rem_num*(base_num+1)) then
          halo_home(k) = (gnum-1)/(base_num+1)
        else
          halo_home(k) = rem_num + (gnum - 1 - rem_num*(base_num+1))/base_num
        end if
      end do
      allocate(send_count(st_mpi%totn), fill_pos(st_mpi%totn))
      do proc_id = 1, st_mpi%totn
        send_count(proc_id) = 0
      end do
      do k = 1, nhalo
        dest_rank = halo_home(k) + 1 ; send_count(dest_rank) = send_count(dest_rank) + 1
      end do
      send_tot = 0
      do proc_id = 1, st_mpi%totn
        fill_pos(proc_id) = send_tot ; send_tot = send_tot + send_count(proc_id)
      end do
      do k = 1, nhalo
        dest_rank = halo_home(k) + 1
        fill_pos(dest_rank) = fill_pos(dest_rank) + 1
        q_gnum(fill_pos(dest_rank)) = halo_gnum(k)
      end do
      allocate(recv_count(st_mpi%totn))
      call alltoall_val(send_count, "hash count", recv_count)
      recv_tot = 0
      do proc_id = 1, st_mpi%totn
        recv_tot = recv_tot + recv_count(proc_id)
      end do
      allocate(recv_gnum(recv_tot + 1))
      call alltoallv_val(q_gnum, send_count, recv_count, "hash gnum", recv_gnum)
      allocate(reply_reg(recv_tot + 1), reply_owner(recv_tot + 1))
      do query_pos = 1, recv_tot
        local_idx = recv_gnum(query_pos) - read_sta + 1
        reply_reg(query_pos) = read_reg(local_idx)
        reply_owner(query_pos) = read_mpi(local_idx)
      end do
      call alltoallv_val(reply_reg, recv_count, send_count, "hash reg", ans_reg)
      call alltoallv_val(reply_owner, recv_count, send_count, "hash owner", ans_owner)
      deallocate(halo_home, send_count, fill_pos, recv_count)
      deallocate(recv_gnum, reply_reg, reply_owner)
    else
      do k = 1, nhalo
        q_gnum(k) = halo_gnum(k)
        ans_reg(k) = read_reg(halo_gnum(k) - read_sta + 1)
        ans_owner(k) = read_mpi(halo_gnum(k) - read_sta + 1)
      end do
    end if
#else
    do k = 1, nhalo
      q_gnum(k) = halo_gnum(k)
      ans_reg(k) = read_reg(halo_gnum(k) - read_sta + 1)
      ans_owner(k) = read_mpi(halo_gnum(k) - read_sta + 1)
    end do
#endif

    active_num = 0
    do k = 1, nhalo
      if (ans_reg(k) > 0) then
        active_num = active_num + 1
      end if
    end do
    neib_hash = get_hash_size(active_num)
    allocate(neib_hash_key(neib_hash), neib_hash_mpi(neib_hash), neib_hash_reg(neib_hash))
    !$omp parallel do private(k)
    do k = 1, neib_hash
      neib_hash_key(k) = 0 ; neib_hash_mpi(k) = 0 ; neib_hash_reg(k) = 0
    end do
    !$omp end parallel do
    do k = 1, nhalo
      if (ans_reg(k) > 0) then
        call set_hash_info(q_gnum(k), ans_owner(k), ans_reg(k))
      end if
    end do

    mismatch = 0 ; ndiag = 0
    do own_pos = 1, ncalc
      call get_neib_gnum(l2g_ijk(own_pos), neib_gnum, dir_num)
      do dir_pos = 1, dir_num
        nb_gnum = neib_gnum(dir_pos)
        call get_hash_info(nb_gnum, find_flag, hmpi, hreg)
        bad = 0
        if (glob_reg_flag(nb_gnum) > 0 .and. glob_mpi_flag(nb_gnum) /= st_mpi%rank+1) then
          if (.not. find_flag) then
            bad = 1
          else if (hmpi /= glob_mpi_flag(nb_gnum) .or. hreg /= glob_reg_flag(nb_gnum)) then
            bad = 1
          end if
        else
          if (find_flag) then
            bad = 1
          end if
        end if
        if (bad == 1) then
          mismatch = 1
          if (ndiag < 10) then
            rf = 0
            if (find_flag) then
              rf = 1
            end if
            write(*,'(a,7i9)') "A3C-HASH-MISS:", st_mpi%rank, nb_gnum, &
              glob_mpi_flag(nb_gnum), glob_reg_flag(nb_gnum), rf, hmpi, hreg
            flush(6)
            ndiag = ndiag + 1
          end if
        end if
      end do
    end do
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      mism_loc(1) = mismatch
      call mpimax_val(mism_loc, "hash exchange verification", mism_all)
      mismatch = mism_all(1)
    end if
#endif
    if (mismatch /= 0) then
      call write_err_stop("set hash neighbor exchange vs scaffold mismatch.")
    end if

    deallocate(cand_list, halo_gnum, q_gnum, ans_reg, ans_owner)

  end subroutine set_hash_neib

  subroutine set_hash_info(hash_key, hash_mpi, hash_reg)
  !*********************************************************************************************
  ! set_hash_info -- Set hash infomation
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: hash_key, hash_mpi, hash_reg
    ! -- local
    integer(I4) :: slot
    !-------------------------------------------------------------------------------------------
    slot = mod(hash_key - 1, neib_hash) + 1

    do while (neib_hash_key(slot) /= 0 .and. neib_hash_key(slot) /= hash_key)
      slot = mod(slot, neib_hash) + 1
    end do

    neib_hash_key(slot) = hash_key
    neib_hash_mpi(slot) = hash_mpi
    neib_hash_reg(slot) = hash_reg

  end subroutine set_hash_info

  subroutine get_hash_info(hash_key, find_flag, neib_mpi, neib_reg)
  !*********************************************************************************************
  ! get_hash_info -- Get hash infomation
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: hash_key
    logical, intent(out) :: find_flag
    integer(I4), intent(out) :: neib_mpi, neib_reg
    ! -- local
    integer(I4) :: slot
    !-------------------------------------------------------------------------------------------
    slot = mod(hash_key - 1, neib_hash) + 1

    do while (neib_hash_key(slot) /= 0 .and. neib_hash_key(slot) /= hash_key)
      slot = mod(slot, neib_hash) + 1
    end do

    if (neib_hash_key(slot) == hash_key) then
      find_flag = .true. ; neib_mpi = neib_hash_mpi(slot) ; neib_reg = neib_hash_reg(slot)
    else
      find_flag = .false. ; neib_mpi = 0 ; neib_reg = 0
    end if

  end subroutine get_hash_info

  function get_hash_size(table_num) result(hash_size)
  !*********************************************************************************************
  ! get_hash_size -- Get hash table size
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: table_num
    ! -- local
    integer(I4) :: hash_size
    !-------------------------------------------------------------------------------------------
    hash_size = table_num*2 + (table_num/2) + 3

    if (mod(table_num, 2) == 0) then
      hash_size = hash_size + 1
    end if

  end function get_hash_size

end module set_cell
