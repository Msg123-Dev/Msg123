module set_cell
  ! -- modules
  use kind_module, only: I4, SP
  use constval_module, only: SNOVAL, CHALEN, VARLEN
  use types_module, only: conn_set, gmap_set
  use utility_module, only: st_mpi, iquick_sort, iquick_sort2, open_new_rtxt, write_err_stop,&
                            gmap_init, gmap_put, gmap_get, gmap_free
  use initial_module, only: st_sim, st_grid, st_clas, in_type, st_ctrl
  use read_module, only: read_2d_calcreg, read_3d_calcreg
#ifdef MPI_MSG
  use mpi_utility, only: mpisum_val, mpimax_val, mpiexscan_val, bcast_val, scatterv_val
  use mpi_utility, only: alltoall_val, alltoallv_val
#endif

  implicit none
  private
  public :: set_cell_info, get_cals_grid, get_calc_grid
  integer(I4), public :: amg_setflag
  integer(I4), public :: ncalc, ncals, ncell, nsurf, no_ncalc, no_ncals
  integer(I4), public :: seal_snum, seal_cnum, neib_mpi_totn, neib_ncals, neib_ncalc
  type(conn_set), public :: st_conn
#ifdef MPI_MSG
  integer(I4), allocatable, public :: send_cind(:), recv_cind(:)
  integer(I4), allocatable, public :: send_citem(:), recv_citem(:)
  integer(I4), allocatable, public :: neib_num(:), send2recv(:), calc2recv(:)
  integer(I4), allocatable, public :: loc2glo_nos(:), loc2glo_noc(:)
  integer(I4), allocatable, public :: loc2unk_ij(:)
#endif

  ! -- local
  integer(I4) :: totnreg, loc_regn, neib_hash, read_sta, read_num
  integer(I4), allocatable :: calc_end(:), loc_nreg(:)
  integer(I4), allocatable :: l2g_ij(:), l2g_ijk(:)
  integer(I4), allocatable :: neib_hash_key(:), neib_hash_mpi(:), neib_hash_reg(:)
  integer(I4), allocatable :: own_calc_glo(:), own_calc_reg(:), own_nocal_glo(:)
  integer(I4), allocatable :: read_reg(:), read_mpi(:)
  logical :: div_reg_3d
  type(gmap_set) :: surf_map
#ifndef MPI_MSG
  integer(I4), allocatable :: glob_reg_flag(:), glob_mpi_flag(:)
#else
  integer(I4) :: part_sta, part_num
  integer(I4), allocatable :: own_part_reg(:), own_part_mpi(:)
#endif

  contains

  subroutine set_cell_info()
  !*********************************************************************************************
  ! set_cell_info -- Set cell information
  !*********************************************************************************************
    ! -- modules
    use utility_module, only: open_new_wtxt, close_file
    use initial_module, only: out_type, st_out_type, st_out_path
    use check_condition, only: read_seal_set, read_seal_clasf, read_seal_point, read_sea_allv
#ifdef MPI_MSG
    use mpi_utility, only: barrier_proc, bcast_file
    use mpi_read, only: read_dist_seaval
    use mpi_set, only: bcast_sim_flag, bcast_xyz_num, bcast_clas_set, bcast_clas_val,&
                       set_calc_view, set_seal_view, set_rest_view, set_write_fview,&
                       senrec_reg_info, senrec_grid_num
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i, j, n, nxy, nc_unknow, mpi_ncals, mpi_ncalc, calg_num
    integer(I4) :: sta_calc, end_calc, i_num, j_num, k_num, pro_nreg, loc_num
    integer(I4) :: sea_ftype, sea_ptn, sea_nclas, sea_null
    integer(I4), allocatable :: cur_nreg(:), sp_i(:), sp_j(:), sp_k(:), sea_hit(:), pre_sea(:)
    real(SP), allocatable :: sp_v(:), sea_clas_val(:), read_seaval(:)
    character(CHALEN) :: sea_path
    character(VARLEN), allocatable :: sea_clas_name(:)
#ifdef MPI_MSG
    integer(I4) :: sea_chk(1), sea_all(1)
    integer(I4), allocatable :: sort_sglo(:), sort_cglo(:), loc2unk_ijk(:)
    real(SP), allocatable :: glo_seaval(:)
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

#ifndef MPI_MSG
    allocate(glob_reg_flag(st_grid%nxyz), glob_mpi_flag(st_grid%nxyz))
    !$omp parallel do private(i)
    do i = 1, st_grid%nxyz
      glob_reg_flag(i) = 0 ; glob_mpi_flag(i) = 0
    end do
    !$omp end parallel do
#endif

    nxy = st_grid%nx*st_grid%ny

    ! -- Partition dimension (2d columns / 3d cells). 3d is selected by:
    !      div_reg_3d = st_sim%reg_neib /= 1 .and. &
    !                   (st_sim%reg_type == in_type(5) .or. st_sim%reg_type == in_type(6) .or. &
    !                    st_mpi%totn >= nxy)
    !    Forced 2d until the 3d downstream (column maps, surface comm, seepage) is ready.
    div_reg_3d = .false.

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      call bcast_file(st_sim%reg_name, "calculation region name")
    end if
#endif

#ifdef MPI_MSG
    ! -- Set own rank's region (own_reg)
      call set_own_reg()
    ! -- Divide calculation region, distributed (calc_reg_dist)
      call div_calc_reg_dist()
    ! -- Divide no calculation flag, distributed (nocalc_dist)
      call div_nocalc_dist()
#else
    ! -- Set global region (glob_reg)
      call set_glob_reg()
    ! -- Divide calculation region (calc_reg)
      call div_calc_reg()
    ! -- Divide no calculation flag (nocalc_flag)
      call div_nocalc_flag()
#endif

    ! -- Set read range (read_range)
      call set_read_range()

    sea_ftype = 0 ; sea_ptn = 0 ; sea_nclas = 0 ; sea_path = repeat(' ', CHALEN)
    if (st_mpi%rank == 0) then
      ! -- Read sea level file set (seal_set)
        call read_seal_set(sea_ftype, sea_path)
      if (sea_ftype == in_type(1)) then
        ! -- Read sea level class file (seal_clasf)
          call read_seal_clasf(sea_nclas, sea_clas_name, sea_clas_val)
      else if (sea_ftype == in_type(2)) then
        ! -- Read sea level point file (seal_point)
          call read_seal_point(sea_ptn, sp_i, sp_j, sp_k, sp_v)
      end if
    end if

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      call bcast_val(sea_ftype, "sea input type")
      call bcast_file(sea_path, "sea input path")
      if (sea_ftype == in_type(1)) then
        call bcast_val(sea_nclas, "sea class count")
        if (st_mpi%rank /= 0) then
          allocate(sea_clas_name(max(sea_nclas,1)), sea_clas_val(max(sea_nclas,1)))
          do i = 1, max(sea_nclas,1)
            sea_clas_name(i) = "" ; sea_clas_val(i) = SNOVAL
          end do
        end if
        if (sea_nclas > 0) then
          call bcast_clas_val(sea_nclas, sea_clas_name, sea_clas_val)
        end if
      else if (sea_ftype == in_type(2)) then
        call bcast_val(sea_ptn, "sea point count")
        if (st_mpi%rank /= 0) then
          allocate(sp_i(sea_ptn+1), sp_j(sea_ptn+1), sp_k(sea_ptn+1))
        end if
        if (st_mpi%rank /= 0) then
          allocate(sp_v(sea_ptn+1))
        end if
        if (sea_ptn > 0) then
          call bcast_val(sp_i, "sea pi") ; call bcast_val(sp_j, "sea pj")
          call bcast_val(sp_k, "sea pk") ; call bcast_val(sp_v, "sea pv")
        end if
      end if
    end if
#endif

    allocate(sea_hit(max(read_num,1)), pre_sea(max(read_num,1)))
    do i = 1, max(read_num,1)
      sea_hit(i) = 0 ; pre_sea(i) = 0
    end do
    do i = 1, read_num
      if (read_reg(i) == 0) then
        pre_sea(i) = 1
      end if
    end do

    if (sea_ftype /= 0) then
      allocate(read_seaval(read_num+1))
      do i = 1, read_num+1
        read_seaval(i) = SNOVAL
      end do
      if (sea_ftype == in_type(1)) then
        ! -- Set sea level value from class file (seaval_clasf)
          call set_seaval_clasf(sea_nclas, sea_clas_name, sea_clas_val, read_seaval)
      else if (sea_ftype == in_type(2)) then
        ! -- Set sea level value from point file (seaval_point)
          call set_seaval_point(sea_ptn, sp_i, sp_j, sp_k, sp_v, read_seaval)

#ifdef MPI_MSG
      else if (sea_ftype == in_type(4) .or. sea_ftype == in_type(6)) then
        call read_dist_seaval(sea_ftype, trim(sea_path), read_sta, read_num, read_seaval)
      else
        if (st_mpi%rank == 0) then
          allocate(glo_seaval(st_grid%nxyz))
          do i = 1, st_grid%nxyz
            glo_seaval(i) = SNOVAL
          end do
          call read_sea_allv(glo_seaval)
        else
          allocate(glo_seaval(1))
          glo_seaval(1) = SNOVAL
        end if
        if (st_mpi%totn /= 1) then
          call scatterv_val(st_mpi%totn, read_num, glo_seaval, read_seaval, "sea all value")
        else
          do i = 1, read_num
            read_seaval(i) = glo_seaval(read_sta + i - 1)
          end do
        end if
        deallocate(glo_seaval)
#else
      else
        call read_sea_allv(read_seaval)
#endif
      end if
      call set_own_seaval(read_seaval, sea_hit)
      deallocate(read_seaval)
    end if

    if (sea_ftype /= 0) then
      sea_null = 0
      do i = 1, read_num
        if (pre_sea(i) == 1 .and. sea_hit(i) == 0) then
          sea_null = 1
          exit
        end if
      end do
#ifdef MPI_MSG
      sea_chk(1) = sea_null
      call mpimax_val(sea_chk, "sea region check", sea_all)
      sea_null = sea_all(1)
#endif
      if (sea_null /= 0) then
        call write_err_stop("Null value in sea region.")
      end if
    end if
    deallocate(sea_hit, pre_sea)
    if (allocated(sp_i)) then
      deallocate(sp_i, sp_j, sp_k)
    end if
    if (allocated(sp_v)) then
      deallocate(sp_v)
    end if
    if (allocated(sea_clas_name)) then
      deallocate(sea_clas_name, sea_clas_val)
    end if

    ! -- Redistribute owned cells to their owner ranks (owned2owner)
      call redist_owned2owner()

    ncalc = size(own_calc_glo)
    ncals = 0
    do i = 1, ncalc
      if (own_calc_glo(i) <= nxy) then
        ncals = ncals + 1
      end if
    end do

    allocate(l2g_ijk(ncalc), l2g_ij(ncals))
    call gmap_init(st_conn%glo2loc_map, ncalc)
    !$omp parallel
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

#ifndef MPI_MSG
    deallocate(glob_reg_flag)
#endif

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
    do i = 1, st_conn%glo2loc_map%table_num
      loc_num = st_conn%glo2loc_map%table_val(i)
      if (st_conn%glo2loc_map%table_key(i) /= 0 .and. loc_num > mpi_ncalc) then
        st_conn%loc2glo_ijk(loc_num) = st_conn%glo2loc_map%table_key(i)
      end if
    end do

#ifdef MPI_MSG
    ! -- Add orphan sea cells to nocalc (orphan nocalc)
      call add_orphan_nocalc()
#endif
    deallocate(read_reg, read_mpi)

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

    deallocate(l2g_ij, l2g_ijk)

#ifdef MPI_MSG
    deallocate(own_part_reg, own_part_mpi)
#else
    deallocate(glob_mpi_flag)
#endif

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

    ! -- Set distributed unknown numbering (dist_unknum)
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

#ifdef MPI_MSG
  subroutine set_own_reg()
  !*********************************************************************************************
  ! set_own_reg -- Set own rank's region
  !*********************************************************************************************
    ! -- modules
    use mpi_read, only: read_mpi_calcreg
    ! -- inout

    ! -- local
    integer(I4) :: i, j, g, nxy, part_size, pbase, prem, reg_type, calreg_fnum
    integer(I4) :: len_char, first_pos, sta_pos, end_pos, glob_ncalc
    integer(I4) :: tot_chk
    integer(I4) :: mism_loc(1), mism_all(1)
    integer(I4), allocatable :: cflag(:), type_2d(:), type_3d(:), glo_tmp(:)
    character(1), allocatable :: temp_char(:)
    character(VARLEN), allocatable :: reg_name(:)
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    if (div_reg_3d) then
      part_size = st_grid%nxyz
    else
      part_size = nxy
    end if
    pbase = part_size/st_mpi%totn ; prem = mod(part_size, st_mpi%totn)
    part_num = pbase
    if (st_mpi%rank < prem) then
      part_num = pbase + 1
    end if
    part_sta = st_mpi%rank*pbase + min(st_mpi%rank, prem) + 1

    allocate(own_part_reg(max(part_num,1)), own_part_mpi(max(part_num,1)))
    do i = 1, max(part_num,1)
      own_part_reg(i) = 0 ; own_part_mpi(i) = 0
    end do

    reg_type = st_sim%reg_type
    allocate(type_2d(2), type_3d(2))
    type_2d(:) = [in_type(3:4)] ; type_3d(:) = [in_type(5:6)]

    if (reg_type == in_type(0)) then
      do i = 1, part_num
        own_part_reg(i) = 1
      end do
      totnreg = 1

    else if (reg_type == in_type(1)) then
      len_char = len(trim(adjustl(st_sim%reg_name)))
      allocate(temp_char(len_char))
      temp_char = transfer(st_sim%reg_name, ' ', size = len_char)
      allocate(reg_name(count(temp_char(:) == ",") + 1))
      first_pos = index(st_sim%reg_name, ",")
      sta_pos = 1 ; end_pos = first_pos
      do i = 1, size(reg_name) - 1
        reg_name(i) = st_sim%reg_name(sta_pos:end_pos-1) ; sta_pos = end_pos + 1
        end_pos = index(st_sim%reg_name(sta_pos:), ",") + sta_pos - 1
      end do
      reg_name(size(reg_name)) = st_sim%reg_name(sta_pos:)
      deallocate(temp_char)

      allocate(cflag(max(part_num,1)))
      do j = 1, st_clas%totn
        call set_clas_flag(j, part_sta, part_num, cflag)
        if (trim(adjustl(st_sim%inact_name)) == trim(adjustl(st_clas%name(j)))) then
          do i = 1, part_num
            if (cflag(i) == 1) then
              own_part_reg(i) = -1
            end if
          end do
        end if
        do i = 1, size(reg_name)
          if (trim(adjustl(reg_name(i))) == trim(adjustl(st_clas%name(j)))) then
            do g = 1, part_num
              if (cflag(g) == 1) then
                own_part_reg(g) = i
              end if
            end do
          end if
        end do
      end do
      totnreg = size(reg_name)
      deallocate(cflag, reg_name)

    else if (reg_type == type_2d(2) .or. reg_type == type_3d(2)) then
      call read_mpi_calcreg(reg_type, div_reg_3d, part_sta, part_num, own_part_reg, totnreg)

    else if (reg_type == type_2d(1) .or. reg_type == type_3d(1)) then
      tot_chk = 0
      if (st_mpi%rank == 0) then
        allocate(glo_tmp(st_grid%nxyz))
        do i = 1, st_grid%nxyz
          glo_tmp(i) = 0
        end do
        call open_new_rtxt(1, 1, st_sim%reg_name, "calculation reigion", calreg_fnum)
        if (reg_type == type_2d(1)) then
          call read_2d_calcreg(calreg_fnum, reg_type, glo_tmp, glob_ncalc)
        else
          call read_3d_calcreg(calreg_fnum, reg_type, glo_tmp, glob_ncalc)
        end if
        tot_chk = maxval(glo_tmp(:))
      else
        allocate(glo_tmp(1))
        glo_tmp(1) = 0
      end if
      call scatterv_val(st_mpi%totn, part_num, glo_tmp, own_part_reg, "part region")
      call bcast_val(tot_chk, "text region count")
      totnreg = tot_chk
      deallocate(glo_tmp)
    end if
    deallocate(type_2d, type_3d)

    if (reg_type /= in_type(0)) then
      mism_loc(1) = 0
      do i = 1, part_num
        if (own_part_reg(i) == -1) then
          mism_loc(1) = 1
          exit
        end if
      end do
      call mpimax_val(mism_loc, "inactive region check", mism_all)
      if (mism_all(1) == 0) then
        call write_err_stop("check the no calculation region name.")
      end if
    end if

  end subroutine set_own_reg
#endif

#ifdef MPI_MSG
  subroutine div_calc_reg_dist()
  !*********************************************************************************************
  ! div_calc_reg_dist -- Divide calculation region, distributed
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, k, part_pos, nxy, ntot
    integer(I4) :: reg_id, pos_in_reg, sta_mpi, end_mpi, reg_nproc, proc_cals, remain_cals
    integer(I4) :: proc_pos, rough_divn, sum_cals, sum_mpi, pro_num, quot, max_reg
    integer(I4), allocatable :: reg_ncals(:), reg_remain(:), reg_mpi_num(:), reg_mpi_end(:)
    integer(I4), allocatable :: my_reg_ncals(:), cals_end(:), seen(:), rang_mpi(:)
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    if (div_reg_3d) then
      ntot = st_grid%nxyz
    else
      ntot = nxy
    end if

    allocate(my_reg_ncals(totnreg))
    !$omp parallel do private(i)
    do i = 1, totnreg
      my_reg_ncals(i) = 0
    end do
    !$omp end parallel do
    do part_pos = 1, part_num
      reg_id = own_part_reg(part_pos)
      if (reg_id >= 1) then
        my_reg_ncals(reg_id) = my_reg_ncals(reg_id) + 1
      end if
    end do

    allocate(reg_ncals(totnreg))
    call mpisum_val(my_reg_ncals, "region column count", reg_ncals)

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

    allocate(cals_end(totnreg))
    call mpiexscan_val(my_reg_ncals, "region column prefix", cals_end)

    allocate(rang_mpi(part_num), seen(totnreg))
    !$omp parallel
    !$omp do private(i)
    do i = 1, part_num
      rang_mpi(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, totnreg
      seen(i) = 0
    end do
    !$omp end do
    !$omp end parallel
    do part_pos = 1, part_num
      reg_id = own_part_reg(part_pos)
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
        rang_mpi(part_pos) = sta_mpi + proc_pos
      end if
    end do

    do part_pos = 1, part_num
      own_part_mpi(part_pos) = rang_mpi(part_pos)
    end do

    deallocate(my_reg_ncals, reg_ncals, reg_remain, reg_mpi_num, reg_mpi_end)
    deallocate(cals_end, seen, rang_mpi)

  end subroutine div_calc_reg_dist
#endif

#ifdef MPI_MSG
  subroutine div_nocalc_dist()
  !*********************************************************************************************
  ! div_nocalc_dist -- Divide no calculation flag, distributed
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: part_pos, nxy, ntot
    integer(I4) :: glo_nocals, rough_divn, nocals_remain, pos_nocal, proc_pos, seen
    integer(I4) :: my_nocal(1), all_nocal(1), nocal_end(1)
    integer(I4), allocatable :: rang_mpi(:)
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    if (div_reg_3d) then
      ntot = st_grid%nxyz
    else
      ntot = nxy
    end if

    my_nocal(1) = 0
    do part_pos = 1, part_num
      if (own_part_mpi(part_pos) == 0) then
        my_nocal(1) = my_nocal(1) + 1
      end if
    end do

    call mpisum_val(my_nocal, "nocalc count", all_nocal)
    glo_nocals = all_nocal(1)
    rough_divn = glo_nocals/st_mpi%totn
    nocals_remain = mod(glo_nocals, st_mpi%totn)

    call mpiexscan_val(my_nocal, "nocalc prefix", nocal_end)

    allocate(rang_mpi(part_num))
    seen = 0
    do part_pos = 1, part_num
      if (own_part_mpi(part_pos) == 0) then
        seen = seen + 1
        pos_nocal = nocal_end(1) + seen
        if (pos_nocal <= nocals_remain*(rough_divn+1)) then
          proc_pos = (pos_nocal-1)/(rough_divn+1) + 1
        else
          proc_pos = nocals_remain + (pos_nocal-1 - nocals_remain*(rough_divn+1))/rough_divn + 1
        end if
        rang_mpi(part_pos) = -proc_pos
      else
        rang_mpi(part_pos) = own_part_mpi(part_pos)
      end if
    end do

    do part_pos = 1, part_num
      own_part_mpi(part_pos) = rang_mpi(part_pos)
    end do

    deallocate(rang_mpi)

  end subroutine div_nocalc_dist
#endif

#ifndef MPI_MSG
  subroutine set_glob_reg()
  !*********************************************************************************************
  ! set_glob_reg -- Set global region
  !*********************************************************************************************
    ! -- modules
    use utility_module, only: open_new_rbin
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
        call set_clas_flag(j, 1, st_grid%nxyz, cflag)
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
          call open_new_rbin(1, 1, st_sim%reg_name, "calculation reigion", calreg_fnum)
          if (reg_type == type_2d(2)) then
            call read_2d_calcreg(calreg_fnum, reg_type, glob_reg_flag, glob_ncalc)
          else
            call read_3d_calcreg(calreg_fnum, reg_type, glob_reg_flag, glob_ncalc)
          end if
          totnreg = maxval(glob_reg_flag(:))
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
#endif

#ifndef MPI_MSG
  subroutine div_calc_reg()
  !*********************************************************************************************
  ! div_calc_reg -- Divide calculation region
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

    if (.not. div_reg_3d) then
      !$omp parallel do private(i)
      do i = 1, st_grid%nz-1
        glob_mpi_flag(nxy*i+1:nxy*(i+1)) = glob_mpi_flag(1:nxy)
      end do
      !$omp end parallel do
    end if

  end subroutine div_calc_reg
#endif

#ifndef MPI_MSG
  subroutine div_nocalc_flag()
  !*********************************************************************************************
  ! div_nocalc_flag -- Divide no calculation flag
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
#endif

  subroutine set_read_range()
  !*********************************************************************************************
  ! set_read_range -- Set read range
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, base_num, rem_num
#ifdef MPI_MSG
    integer(I4) :: g, k, key, nxy, part_size, pbase, prem, psta
    integer(I4) :: cand_raw, nkey, proc_id, dest_rank, send_tot, recv_tot, pos
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
        part_size = st_grid%nxyz
      else
        part_size = nxy
      end if
      pbase = part_size/st_mpi%totn ; prem = mod(part_size, st_mpi%totn)
      psta = part_sta

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

      allocate(key_home(nkey + 1), send_count(st_mpi%totn), fill_pos(st_mpi%totn))
      do k = 1, nkey
        key = key_list(k)
        if (key <= prem*(pbase+1)) then
          key_home(k) = (key-1)/(pbase+1)
        else
          key_home(k) = prem + (key - 1 - prem*(pbase+1))/pbase
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
        pos = recv_key(i) - psta + 1
        reply_reg(i) = own_part_reg(pos) ; reply_owner(i) = own_part_mpi(pos)
      end do
      allocate(ans_reg(send_tot + 1), ans_owner(send_tot + 1))
      call alltoallv_val(reply_reg, recv_count, send_count, "read range reg", ans_reg)
      call alltoallv_val(reply_owner, recv_count, send_count, "read range owner", ans_owner)

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
      deallocate(cand_list, key_list, key_home)
      deallocate(send_count, fill_pos, recv_count, q_key, recv_key)
      deallocate(reply_reg, reply_owner, ans_reg, ans_owner)

    else
      nxy = st_grid%nx*st_grid%ny
      do i = 1, read_num
        g = read_sta + i - 1
        if (div_reg_3d) then
          key = g
        else
          key = mod(g-1, nxy) + 1
        end if
        read_reg(i) = own_part_reg(key) ; read_mpi(i) = own_part_mpi(key)
      end do

    end if
#else
    do i = 1, read_num
      read_reg(i) = glob_reg_flag(read_sta + i - 1)
      read_mpi(i) = glob_mpi_flag(read_sta + i - 1)
    end do
#endif

  end subroutine set_read_range

  subroutine set_seaval_clasf(nclas, clas_name, clas_val, seaval)
  !*********************************************************************************************
  ! set_seaval_clasf -- Set sea level value from class file
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: nclas
    character(VARLEN), intent(in) :: clas_name(:)
    real(SP), intent(in) :: clas_val(:)
    real(SP), intent(inout) :: seaval(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4), allocatable :: cflag(:)
    !-------------------------------------------------------------------------------------------
    allocate(cflag(max(read_num,1)))
    do i = 1, nclas
      if (clas_val(i) == SNOVAL) then
        cycle
      end if
      do k = 1, st_clas%totn
        if (trim(adjustl(clas_name(i))) /= trim(adjustl(st_clas%name(k)))) then
          cycle
        end if
        call set_clas_flag(k, read_sta, read_num, cflag)
        do j = 1, read_num
          if (cflag(j) == 1) then
            seaval(j) = clas_val(i)
          end if
        end do
      end do
    end do
    deallocate(cflag)

  end subroutine set_seaval_clasf

  subroutine set_seaval_point(sea_ptn, sp_i, sp_j, sp_k, sp_v, seaval)
  !*********************************************************************************************
  ! set_seaval_point -- Set sea level value from point file
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: sea_ptn
    integer(I4), intent(in) :: sp_i(:), sp_j(:), sp_k(:)
    real(SP), intent(in) :: sp_v(:)
    real(SP), intent(inout) :: seaval(:)
    ! -- local
    integer(I4) :: li, g, p, nxy, i_num, j_num, k_num, ij
    !-------------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    do li = 1, read_num
      g = read_sta + li - 1
      k_num = (g-1)/nxy + 1
      ij = g - nxy*(k_num-1)
      j_num = (ij-1)/st_grid%nx + 1
      i_num = ij - (j_num-1)*st_grid%nx
      do p = 1, sea_ptn
        if ((sp_i(p) == -1 .or. sp_i(p) == i_num) .and. &
            (sp_j(p) == -1 .or. sp_j(p) == j_num) .and. &
            (sp_k(p) == -1 .or. sp_k(p) == k_num)) then
          seaval(li) = sp_v(p)
        end if
      end do
    end do

  end subroutine set_seaval_point

  subroutine set_own_seaval(seaval, sea_hit)
  !*********************************************************************************************
  ! set_own_seaval -- Set owned sea cell values
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(SP), intent(in) :: seaval(:)
    integer(I4), intent(inout) :: sea_hit(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    do i = 1, read_num
      if (seaval(i) /= SNOVAL) then
        sea_hit(i) = 1
        if (read_reg(i) > 0) then
          read_reg(i) = 0 ; read_mpi(i) = 0
        end if
      end if
    end do

  end subroutine set_own_seaval

  subroutine redist_owned2owner()
  !*********************************************************************************************
  ! redist_owned2owner -- Redistribute owned cells to their owner ranks
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
#ifdef MPI_MSG
    integer(I4) :: cell_id
    integer(I4) :: cell_sta, cell_num, base_num, rem_num, list_pos
    integer(I4) :: owner_code, dest_rank, proc_id, buf_pos
    integer(I4) :: send_tot, recv_tot, calc_cnt, nocal_cnt
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

  end subroutine redist_owned2owner

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
      l2g_ijk(count_calc) = ij
      call gmap_put(st_conn%glo2loc_map, ij, count_calc)
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
    integer(I4), allocatable :: temp_neib_num(:), temp_neib_flag(:), temp_mpi_num(:)
    integer(I4), allocatable :: neib_glos(:), neib_gloc(:), neib_locc(:)
    integer(I4), allocatable :: sort_recv_num(:), temp_calc_reg(:)
    integer(I4), allocatable :: recv_num(:), send_num(:)
    integer(I4), allocatable :: temp_sort(:), sort_mpi_num(:), loc_send_num(:)
    integer(I4), allocatable :: mpi_l2g_ij(:), mpi_l2g_ijk(:), mpi_calc2reg(:)
    logical :: find_flag
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
        own_reg = own_calc_reg(i)
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
        own_reg = own_calc_reg(i)
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
    integer(I4) :: cand_num, surf_num, volm_num, sdir_num, cdir_num, neib_mpi_num, reg_num
    integer(I4) :: neib_snum(4), neib_gnum(6)
    logical :: find_flag, is_sea
    !-------------------------------------------------------------------------------------------
    seal_snum = 0
    nxy = st_grid%nx*st_grid%ny

    surf_num = ncals + neib_ncals
    call gmap_init(surf_map, surf_num)
    do i = 1, surf_num
      call gmap_put(surf_map, l2g_ij(i), i)
    end do

    cand_num = 0
    do i = 1, ncals
      call get_neib_snum(l2g_ij(i), neib_snum, sdir_num)
      do j = 1, sdir_num
        call get_hash_info(neib_snum(j), find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(surf_map, neib_snum(j)) == 0) then
          cand_num = cand_num + 1
        end if
      end do
    end do

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
        call get_hash_info(n_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(surf_map, n_num) == 0) then
          seal_snum = seal_snum + 1
          call gmap_put(surf_map, n_num, ncals + neib_ncals + seal_snum)
        end if
      end if
      ! west direction
      if (i_num /= 1) then
        w_num = g_num-1
        call get_hash_info(w_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(surf_map, w_num) == 0) then
          seal_snum = seal_snum + 1
          call gmap_put(surf_map, w_num, ncals + neib_ncals + seal_snum)
        end if
      end if
      ! east direction
      if (i_num /= st_grid%nx) then
        e_num = g_num+1
        call get_hash_info(e_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(surf_map, e_num) == 0) then
          seal_snum = seal_snum + 1
          call gmap_put(surf_map, e_num, ncals + neib_ncals + seal_snum)
        end if
      end if
      ! south direction
      if (j_num /= st_grid%ny) then
        s_num = g_num+st_grid%nx
        call get_hash_info(s_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(surf_map, s_num) == 0) then
          seal_snum = seal_snum + 1
          call gmap_put(surf_map, s_num, ncals + neib_ncals + seal_snum)
        end if
      end if
    end do

    seal_cnum = 0

    volm_num = ncalc + neib_ncalc
    call gmap_init(st_conn%glo2loc_map, volm_num)
    do i = 1, volm_num
      call gmap_put(st_conn%glo2loc_map, l2g_ijk(i), i)
    end do

    cand_num = 0
    do i = 1, ncalc
      call get_neib_cnum(l2g_ijk(i), neib_gnum, cdir_num)
      do j = 1, cdir_num
        call get_hash_info(neib_gnum(j), find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(st_conn%glo2loc_map, neib_gnum(j)) == 0) then
          cand_num = cand_num + 1
        end if
      end do
    end do

    call gmap_init(st_conn%glo2loc_map, volm_num + cand_num)
    do i = 1, volm_num
      call gmap_put(st_conn%glo2loc_map, l2g_ijk(i), i)
    end do

    do i = 1, ncalc
      g_num = l2g_ijk(i)
      k_num = (g_num-1)/(nxy) + 1
      ij = g_num - nxy*(k_num-1)
      j_num = (ij-1)/st_grid%nx + 1
      i_num = ij - (j_num-1)*st_grid%nx
      ! up direction
      if (k_num /= 1) then
        u_num = g_num-nxy
        call get_hash_info(u_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(st_conn%glo2loc_map, u_num) == 0) then
          seal_cnum = seal_cnum + 1
          call gmap_put(st_conn%glo2loc_map, u_num, ncalc + neib_ncalc + seal_cnum)
        end if
      end if
      ! north direction
      if (j_num /= 1) then
        n_num = g_num-st_grid%nx
        call get_hash_info(n_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(st_conn%glo2loc_map, n_num) == 0) then
          seal_cnum = seal_cnum + 1
          call gmap_put(st_conn%glo2loc_map, n_num, ncalc + neib_ncalc + seal_cnum)
        end if
      end if
      ! west direction
      if (i_num /= 1) then
        w_num = g_num-1
        call get_hash_info(w_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(st_conn%glo2loc_map, w_num) == 0) then
          seal_cnum = seal_cnum + 1
          call gmap_put(st_conn%glo2loc_map, w_num, ncalc + neib_ncalc + seal_cnum)
        end if
      end if
      ! east direction
      if (i_num /= st_grid%nx) then
        e_num = g_num+1
        call get_hash_info(e_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(st_conn%glo2loc_map, e_num) == 0) then
          seal_cnum = seal_cnum + 1
          call gmap_put(st_conn%glo2loc_map, e_num, ncalc + neib_ncalc + seal_cnum)
        end if
      end if
      ! south direction
      if (j_num /= st_grid%ny) then
        s_num = g_num+st_grid%nx
        call get_hash_info(s_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(st_conn%glo2loc_map, s_num) == 0) then
          seal_cnum = seal_cnum + 1
          call gmap_put(st_conn%glo2loc_map, s_num, ncalc + neib_ncalc + seal_cnum)
        end if
      end if
      ! down direction
      if (k_num /= st_grid%nz) then
        d_num = g_num+nxy
        call get_hash_info(d_num, find_flag, neib_mpi_num, reg_num)
        is_sea = find_flag .and. reg_num == 0
        if (is_sea .and. gmap_get(st_conn%glo2loc_map, d_num) == 0) then
          seal_cnum = seal_cnum + 1
          call gmap_put(st_conn%glo2loc_map, d_num, ncalc + neib_ncalc + seal_cnum)
        end if
      end if
    end do

  end subroutine set_rel_seareg

#ifdef MPI_MSG
  subroutine add_orphan_nocalc()
  !*********************************************************************************************
  ! add_orphan_nocalc -- Add orphan sea cells to nocalc
  !*********************************************************************************************
    ! -- modules

    ! -- local
    integer(I4) :: g, n, j, k, norphan, dir_num, read_end, nb_gnum, out_raw, nout
    integer(I4) :: base_num, rem_num, proc_id, dest_rank, send_tot, recv_tot, pos
    integer(I4) :: neib_gnum(6)
    integer(I4), allocatable :: temp_nocal(:), cand_list(:), out_key(:), out_home(:)
    integer(I4), allocatable :: send_count(:), recv_count(:), fill_pos(:)
    integer(I4), allocatable :: q_gnum(:), recv_gnum(:), reply_own(:), ans_own(:)
    logical :: is_orphan
    type(gmap_set) :: out_map
    !-------------------------------------------------------------------------------------------
    read_end = read_sta + read_num - 1
    base_num = st_grid%nxyz/st_mpi%totn ; rem_num = mod(st_grid%nxyz, st_mpi%totn)

    allocate(cand_list(6*read_num + 1))
    out_raw = 0
    do g = read_sta, read_end
      if (read_mpi(g-read_sta+1) /= 0) then
        cycle
      end if
      call get_neib_cnum(g, neib_gnum, dir_num)
      do n = 1, dir_num
        nb_gnum = neib_gnum(n)
        if (nb_gnum < read_sta .or. nb_gnum > read_end) then
          out_raw = out_raw + 1 ; cand_list(out_raw) = nb_gnum
        end if
      end do
    end do
    if (out_raw > 1) then
      call iquick_sort(cand_list, 1, out_raw)
    end if
    allocate(out_key(out_raw + 1))
    nout = 0
    do k = 1, out_raw
      if (k == 1) then
        nout = 1 ; out_key(1) = cand_list(1)
      else if (cand_list(k) /= cand_list(k-1)) then
        nout = nout + 1 ; out_key(nout) = cand_list(k)
      end if
    end do

    allocate(ans_own(nout + 1))
    do k = 1, nout + 1
      ans_own(k) = 0
    end do
    if (st_mpi%totn /= 1) then
      allocate(out_home(nout + 1), send_count(st_mpi%totn), fill_pos(st_mpi%totn))
      do k = 1, nout
        if (out_key(k) <= rem_num*(base_num+1)) then
          out_home(k) = (out_key(k)-1)/(base_num+1)
        else
          out_home(k) = rem_num + (out_key(k) - 1 - rem_num*(base_num+1))/base_num
        end if
      end do
      do proc_id = 1, st_mpi%totn
        send_count(proc_id) = 0
      end do
      do k = 1, nout
        dest_rank = out_home(k) + 1 ; send_count(dest_rank) = send_count(dest_rank) + 1
      end do
      send_tot = 0
      do proc_id = 1, st_mpi%totn
        fill_pos(proc_id) = send_tot ; send_tot = send_tot + send_count(proc_id)
      end do
      allocate(q_gnum(send_tot + 1))
      do k = 1, nout
        dest_rank = out_home(k) + 1
        fill_pos(dest_rank) = fill_pos(dest_rank) + 1
        q_gnum(fill_pos(dest_rank)) = out_key(k)
      end do
      allocate(recv_count(st_mpi%totn))
      call alltoall_val(send_count, "orphan neighbour count", recv_count)
      recv_tot = 0
      do proc_id = 1, st_mpi%totn
        recv_tot = recv_tot + recv_count(proc_id)
      end do
      allocate(recv_gnum(recv_tot + 1), reply_own(recv_tot + 1))
      call alltoallv_val(q_gnum, send_count, recv_count, "orphan neighbour key", recv_gnum)
      do k = 1, recv_tot
        reply_own(k) = read_mpi(recv_gnum(k) - read_sta + 1)
      end do
      call alltoallv_val(reply_own, recv_count, send_count, "orphan neighbour owner", ans_own)
      call gmap_init(out_map, send_tot)
      do k = 1, send_tot
        call gmap_put(out_map, q_gnum(k), k)
      end do
      deallocate(out_home, send_count, fill_pos, recv_count, recv_gnum, reply_own)
    else
      call gmap_init(out_map, max(nout,1))
    end if

    norphan = 0
    do j = 1, 2
      if (j == 2) then
        if (norphan == 0) then
          exit
        end if
        allocate(temp_nocal(size(own_nocal_glo) + norphan))
        do k = 1, size(own_nocal_glo)
          temp_nocal(k) = own_nocal_glo(k)
        end do
        norphan = size(own_nocal_glo)
      end if
      do g = read_sta, read_end
        if (read_mpi(g-read_sta+1) /= 0) then
          cycle
        end if
        call get_neib_cnum(g, neib_gnum, dir_num)
        is_orphan = .true.
        do n = 1, dir_num
          nb_gnum = neib_gnum(n)
          if (nb_gnum >= read_sta .and. nb_gnum <= read_end) then
            if (read_mpi(nb_gnum-read_sta+1) > 0) then
              is_orphan = .false. ; exit
            end if
          else
            pos = gmap_get(out_map, nb_gnum)
            if (pos /= 0) then
              if (ans_own(pos) > 0) then
                is_orphan = .false. ; exit
              end if
            end if
          end if
        end do
        if (is_orphan) then
          norphan = norphan + 1
          if (j == 2) then
            temp_nocal(norphan) = g
          end if
        end if
      end do
    end do

    call gmap_free(out_map)

    deallocate(cand_list, out_key, ans_own)

    if (allocated(q_gnum)) then
      deallocate(q_gnum)
    end if

    if (allocated(temp_nocal)) then
      call move_alloc(temp_nocal, own_nocal_glo)
      call iquick_sort(own_nocal_glo, 1, size(own_nocal_glo))
    end if

  end subroutine add_orphan_nocalc
#endif

#ifdef MPI_MSG
  subroutine set_dist_unknum(own_gnum, own_num, gnum_range, unk_num)
  !*********************************************************************************************
  ! set_dist_unknum -- Set distributed unknown numbering
  !*********************************************************************************************
    ! -- modules

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

    do list_pos = 1, recv_tot
      recv_orig(list_pos) = list_pos
    end do
    if (recv_tot > 1) then
      call iquick_sort2(recv_gnum, recv_orig, 1, recv_tot)
    end if

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

  subroutine set_loc_cell_clas()
  !*********************************************************************************************
  ! set_loc_cell_clas -- Set local cell classification
  !*********************************************************************************************
    ! -- modules

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

  subroutine pack_owned_local()
  !*********************************************************************************************
  ! pack_owned_local -- Pack owned local info
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: cell_id, calc_cnt, nocal_cnt
    !-------------------------------------------------------------------------------------------
    calc_cnt = 0 ; nocal_cnt = 0
    do cell_id = 1, read_num
      if (read_mpi(cell_id) == st_mpi%rank+1) then
        calc_cnt = calc_cnt + 1
      else if (read_mpi(cell_id) == -(st_mpi%rank+1)) then
        nocal_cnt = nocal_cnt + 1
      end if
    end do
    allocate(own_calc_glo(calc_cnt), own_calc_reg(calc_cnt), own_nocal_glo(nocal_cnt))
    calc_cnt = 0 ; nocal_cnt = 0
    do cell_id = 1, read_num
      if (read_mpi(cell_id) == st_mpi%rank+1) then
        calc_cnt = calc_cnt + 1
        own_calc_glo(calc_cnt) = read_sta + cell_id - 1
        own_calc_reg(calc_cnt) = read_reg(cell_id)
      else if (read_mpi(cell_id) == -(st_mpi%rank+1)) then
        nocal_cnt = nocal_cnt + 1
        own_nocal_glo(nocal_cnt) = read_sta + cell_id - 1
      end if
    end do

  end subroutine pack_owned_local

  subroutine set_hash_neib()
  !*********************************************************************************************
  ! set_hash_neib -- Set hash neighbor
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: own_pos, dir_pos, dir_num, nb_gnum, k, cand_raw, nhalo, active_num
    integer(I4) :: base_num, rem_num
    integer(I4) :: neib_gnum(6)
    integer(I4), allocatable :: cand_list(:), halo_gnum(:)
    integer(I4), allocatable :: q_gnum(:), ans_reg(:), ans_owner(:)
#ifdef MPI_MSG
    integer(I4) :: gnum, local_idx, proc_id, dest_rank, send_tot, recv_tot, query_pos
    integer(I4), allocatable :: halo_home(:), send_count(:), recv_count(:), fill_pos(:)
    integer(I4), allocatable :: recv_gnum(:), reply_reg(:), reply_owner(:)
#endif
    !-------------------------------------------------------------------------------------------
    base_num = st_grid%nxyz/st_mpi%totn ; rem_num = mod(st_grid%nxyz, st_mpi%totn)

    allocate(cand_list(6*ncalc + 1))
    cand_raw = 0
    do own_pos = 1, ncalc
      call get_neib_cnum(l2g_ijk(own_pos), neib_gnum, dir_num)
      do dir_pos = 1, dir_num
        nb_gnum = neib_gnum(dir_pos)
        if (gmap_get(st_conn%glo2loc_map, nb_gnum) == 0) then
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

    active_num = nhalo
    neib_hash = get_hash_size(active_num)

    allocate(neib_hash_key(neib_hash), neib_hash_mpi(neib_hash), neib_hash_reg(neib_hash))
    !$omp parallel do private(k)
    do k = 1, neib_hash
      neib_hash_key(k) = 0 ; neib_hash_mpi(k) = 0 ; neib_hash_reg(k) = 0
    end do
    !$omp end parallel do
    do k = 1, nhalo
      call set_hash_info(q_gnum(k), ans_owner(k), ans_reg(k))
    end do

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

  subroutine set_clas_flag(clas_k, sta_num, ran_num, flag)
  !*********************************************************************************************
  ! set_clas_flag -- Set clas flag
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: clas_k, sta_num, ran_num
    integer(I4), intent(out) :: flag(:)
    ! -- local
    integer(I4) :: i, ii, jj, kk, c_num, pos, end_num
    integer(I4) :: clasi, clasj, clask
    !-------------------------------------------------------------------------------------------
    end_num = sta_num + ran_num - 1
    !$omp parallel do private(i)
    do i = 1, ran_num
      flag(i) = 0
    end do
    !$omp end parallel do

    do i = 1, st_clas%num(clas_k)
      clasi = st_clas%i(i,clas_k) ; clasj = st_clas%j(i,clas_k) ; clask = st_clas%k(i,clas_k)
      if (clasi == -1 .and. clasj == -1 .and. clask == -1) then !all cells
        !$omp parallel do private(pos)
        do pos = 1, ran_num
          flag(pos) = 1
        end do
        !$omp end parallel do
      else if (clasi == -1 .and. clasj == -1) then !i,j cells
        do jj = 1, st_grid%ny
          do ii = 1, st_grid%nx
            c_num = st_grid%nx*(st_grid%ny*(clask-1) + (jj-1)) + ii
            if (c_num >= sta_num .and. c_num <= end_num) then
              flag(c_num-sta_num+1) = 1
            end if
          end do
        end do
      else if (clasi == -1 .and. clask == -1) then !i,k cells
        do kk = 1, st_grid%nz
          do ii = 1, st_grid%nx
            c_num = st_grid%nx*(st_grid%ny*(kk-1) + (clasj-1)) + ii
            if (c_num >= sta_num .and. c_num <= end_num) then
              flag(c_num-sta_num+1) = 1
            end if
          end do
        end do
      else if (clasj == -1 .and. clask == -1) then !j,k cells
        do kk = 1, st_grid%nz
          do jj = 1, st_grid%ny
            c_num = st_grid%nx*(st_grid%ny*(kk-1) + (jj-1)) + clasi
            if (c_num >= sta_num .and. c_num <= end_num) then
              flag(c_num-sta_num+1) = 1
            end if
          end do
        end do
      else if (clasi == -1) then !only i cell
        do ii = 1, st_grid%nx
          c_num = st_grid%nx*(st_grid%ny*(clask-1) + (clasj-1)) + ii
          if (c_num >= sta_num .and. c_num <= end_num) then
            flag(c_num-sta_num+1) = 1
          end if
        end do
      else if (clasj == -1) then !only j cell
        do jj = 1, st_grid%ny
          c_num = st_grid%nx*(st_grid%ny*(clask-1) + (jj-1)) + clasi
          if (c_num >= sta_num .and. c_num <= end_num) then
            flag(c_num-sta_num+1) = 1
          end if
        end do
      else if (clask == -1) then !only k cell
        do kk = 1, st_grid%nz
          c_num = st_grid%nx*(st_grid%ny*(kk-1) + (clasj-1)) + clasi
          if (c_num >= sta_num .and. c_num <= end_num) then
            flag(c_num-sta_num+1) = 1
          end if
        end do
      else !others
        c_num = st_grid%nx*(st_grid%ny*(clask-1) + (clasj-1)) + clasi
        if (c_num >= sta_num .and. c_num <= end_num) then
          flag(c_num-sta_num+1) = 1
        end if
      end if
    end do

  end subroutine set_clas_flag

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

  subroutine get_neib_snum(g_num, neib_snum, sdir_num)
  !*********************************************************************************************
  ! get_neib_snum -- Get neighbor surface numbers
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

  subroutine get_neib_cnum(g_num, neib_gnum, dir_num)
  !*********************************************************************************************
  ! get_neib_cnum -- Get neighbor cell numbers
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

  end subroutine get_neib_cnum

end module set_cell
