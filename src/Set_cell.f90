module set_cell
  ! -- modules
  use kind_module, only: I4
  use utility_module, only: iquick_sort
  use initial_module, only: pro_totn, my_rank, st_sim, st_grid, st_clas

  implicit none
  private
  public :: set_cell_info, get_cals_grid, get_calc_grid
  integer(I4), public :: amg_setflag
  integer(I4), public :: ncalc, ncals, ncell, nsurf, no_ncalc, no_ncals
  integer(I4), public :: seal_cnum, neib_mpi_totn, neib_ncalc
  integer(I4), allocatable, public :: clas_flag(:,:)
  integer(I4), allocatable, public :: calc2reg(:)
  integer(I4), allocatable, public :: glo2loc_ijk(:), loc2glo_ijk(:), loc2glo_ij(:)
#ifdef MPI_MSG
  integer(I4), allocatable, public :: send_cind(:), recv_cind(:)
  integer(I4), allocatable, public :: send_citem(:), recv_citem(:)
  integer(I4), allocatable, public :: neib_num(:), send2recv(:), calc2recv(:)
  integer(I4), allocatable, public :: loc2glo_nos(:), loc2glo_noc(:)
  integer(I4), allocatable, public :: loc2unk_ij(:)
#endif

  ! -- local
  integer(I4) :: ns_unknow, nc_unknow, seal_snum, neib_ncals
  integer(I4) :: totnreg, loc_regn
  integer(I4), allocatable :: glob_reg_flag(:), glob_mpi_flag(:)
  integer(I4), allocatable :: glob_clas_flag(:,:)
  integer(I4), allocatable :: calc_end(:), loc_nreg(:)
  integer(I4), allocatable :: l2g_ij(:), l2g_ijk(:)
  integer(I4), allocatable :: glo2loc_ij(:)

  contains

  subroutine set_cell_info()
  !***************************************************************************************
  ! set_cell_info -- Set cell information
  !***************************************************************************************
    ! -- modules
    use utility_module, only: open_new_wtxt, close_file
    use initial_module, only: precon_type, amg_nlevel, st_out_type, st_out_path,&
                              noclas_flag, out_type
    use check_condition, only: check_calc_region
#ifdef MPI_MSG
    use mpi_utility, only: barrier_proc, mpisum_val, gather_val
    use mpi_set, only: bcast_sim_flag, bcast_xyz_num, bcast_clas_set, bcast_glob_flag,&
                       set_calc_view, set_seal_view, set_rest_view, set_write_fview,&
                       senrec_reg_info, senrec_grid_num
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i, j, n, nxy
    integer(I4) :: mpi_ncals, mpi_ncalc
    integer(I4) :: sta_calc, end_calc, i_num, j_num, k_num, pro_nreg
    integer(I4) :: calg_num
    integer(I4), allocatable :: cur_nreg(:), glob_num(:), temp_locs(:), temp_locc(:)
#ifdef MPI_MSG
    integer(I4), allocatable :: cals_glob(:), calc_glob(:), sort_sglo(:), sort_cglo(:)
    integer(I4), allocatable :: nsun_num(:), ncun_num(:)
    integer(I4), allocatable :: glo2unk_ij(:)
    integer(I4), allocatable :: loc2unk_ijk(:), glo2unk_ijk(:)
#endif
    !-------------------------------------------------------------------------------------
#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast simulation flag (sim_flag)
        call bcast_sim_flag()
      ! -- Bcast xyz number (xyz_num)
        call bcast_xyz_num()
    end if
#endif

    if (precon_type == 1) then
      amg_setflag = 0
    else
      amg_nlevel = 1
    end if

    if (noclas_flag /= 1) then
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast classification setting (clas_set)
          call bcast_clas_set()
      end if
#endif
     allocate(glob_clas_flag(st_grid%nxyz,st_clas%totn))
      ! -- Set global cell classification (glob_cell_clas)
        call set_glob_cell_clas()
    else
      st_clas%totn = 1
      allocate(glob_clas_flag(st_grid%nxyz,st_clas%totn))
    end if

    allocate(glob_reg_flag(st_grid%nxyz), glob_mpi_flag(st_grid%nxyz))
    !$omp parallel workshare
    glob_reg_flag(:) = 0 ; glob_mpi_flag(:) = 0
    !$omp end parallel workshare

    nxy = st_grid%nx*st_grid%ny

    if (my_rank == 0) then
      ! -- Set global region (glob_reg)
        call set_glob_reg()
      ! -- Divide calculation region for 2d (calc_reg_2d)
        call div_calc_reg_2d()
      ! -- Divide no calculation flag for 2d (nocalc_flag_2d)
        call div_nocalc_flag_2d()
      ! -- Check the calculation regin (calc_region)
        call check_calc_region(glob_clas_flag, glob_reg_flag)
        allocate(glob_num(count(glob_reg_flag(:) == 0)))
        glob_num(:) = 0
        glob_mpi_flag(:) = unpack(glob_num(:), glob_reg_flag(:) == 0, glob_mpi_flag(:))
        deallocate(glob_num)
!      if (st_sim%reg_neib == 1) then
!        ! -- Divide calculation region for 2d (calc_reg_2d)
!          call div_calc_reg_2d()
!        ! -- Divide no calculation flag for 2d (nocalc_flag_2d)
!          call div_nocalc_flag_2d()
!      else if (st_sim%reg_type /= in_type(5) .and. st_sim%reg_type /= in_type(6) .and. &
!               pro_totn < nxy) then
!        ! -- Divide calculation region for 2d (calc_reg_2d)
!          call div_calc_reg_2d()
!        ! -- Divide no calculation flag for 2d (nocalc_flag_2d)
!          call div_nocalc_flag_2d()
!      else
!        ! -- Divide calculation region for 3d (calc_reg_3d)
!          call div_calc_reg_3d()
!        ! -- Divide no calculation flag for 3d (nocalc_flag_3d)
!          call div_nocalc_flag_3d()
!      end if
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast global flag (glob_flag)
        call bcast_glob_flag(totnreg, glob_reg_flag, glob_mpi_flag)
    end if
#endif

    ncalc = count(glob_mpi_flag(:) == my_rank+1)
    ncals = count(glob_mpi_flag(1:nxy) == my_rank+1)

    allocate(glo2loc_ijk(st_grid%nxyz), l2g_ijk(ncalc))
    allocate(glo2loc_ij(nxy), l2g_ij(ncals))
    !$omp parallel workshare
    glo2loc_ijk(:) = 0 ; l2g_ijk(:) = 0 ; glo2loc_ij(:) = 0 ; l2g_ij(:) = 0
    !$omp end parallel workshare

    ! -- Set relationship between global&local (rel_gloloc)
      call set_rel_gloloc()

    neib_ncals = 0 ; neib_ncalc = 0 ; neib_mpi_totn = 0

#ifdef MPI_MSG
    allocate(calc2recv(ncalc))
    !$omp parallel workshare
    calc2recv(:) = 0
    !$omp end parallel workshare

    ! -- Set mpi relationship (mpi_rel)
      call set_mpi_rel()
#endif

    ! -- Set relationship of sea region (rel_seareg)
      call set_rel_seareg()

    deallocate(glob_reg_flag)

    mpi_ncals = ncals + neib_ncals ; mpi_ncalc = ncalc + neib_ncalc
    nsurf = mpi_ncals + seal_snum ; ncell = mpi_ncalc + seal_cnum

    allocate(loc2glo_ij(nsurf), loc2glo_ijk(ncell))
    allocate(glob_num(st_grid%nxyz))
    allocate(temp_locs(seal_snum), temp_locc(seal_cnum))
    !$omp parallel
    !$omp workshare
    loc2glo_ij(:) = 0 ; loc2glo_ijk(:) = 0 ; temp_locs(:) = 0 ; temp_locc(:) = 0
    !$omp end workshare
    !$omp do private(i)
    do i = 1, st_grid%nxyz
      glob_num(i) = i
    end do
    !$omp end do
    !$omp workshare
    loc2glo_ij(:mpi_ncals) = l2g_ij(:)
    loc2glo_ijk(:mpi_ncalc) = l2g_ijk(:)

    temp_locs(:) = pack(glo2loc_ij(:), glo2loc_ij(:) > mpi_ncals)
    loc2glo_ij((/ temp_locs /)) = pack(glob_num(1:nxy), glo2loc_ij(:) > mpi_ncals)
    temp_locc(:) = pack(glo2loc_ijk(:), glo2loc_ijk(:) > mpi_ncalc)
    loc2glo_ijk((/ temp_locc /)) = pack(glob_num(:), glo2loc_ijk(:) > mpi_ncalc)
    !$omp end workshare
    !$omp end parallel

    no_ncals = count(glob_mpi_flag(1:nxy) == -(my_rank+1))
    no_ncalc = count(glob_mpi_flag(:) == -(my_rank+1))

#ifdef MPI_MSG
    allocate(loc2glo_nos(no_ncals), loc2glo_noc(no_ncalc))
    loc2glo_nos(:) = pack(glob_num(1:nxy), (glob_mpi_flag(1:nxy) == -(my_rank+1)))
    loc2glo_noc(:) = pack(glob_num(:), (glob_mpi_flag(:) == -(my_rank+1)))
#endif

    deallocate(l2g_ij, l2g_ijk, temp_locs, temp_locc, glob_num, glob_mpi_flag, glo2loc_ij)

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(ncals, "calculation surface number", ns_unknow)
      ! -- Sum value for MPI (val)
        call mpisum_val(ncalc, "calculation number", nc_unknow)
    else
      ns_unknow = ncals ; nc_unknow = ncalc
    end if
#else
    ns_unknow = ncals ; nc_unknow = ncalc
#endif

#ifdef MPI_MSG
    allocate(cals_glob(ns_unknow), calc_glob(nc_unknow))
    !$omp parallel workshare
    cals_glob(:) = 0 ; calc_glob(:) = 0
    !$omp end parallel workshare
    if (pro_totn /= 1) then
      ! -- Gather array (val)
        call gather_val(pro_totn, ncals, loc2glo_ij(1:ncals), cals_glob, "calculation number")
      ! -- Gather array (val)
        call gather_val(pro_totn, ncalc, loc2glo_ijk(1:ncalc), calc_glob, "calculation number")
    else
      !$omp parallel workshare
      cals_glob(:) = loc2glo_ij(1:ncals) ; calc_glob(:) = loc2glo_ijk(1:ncalc)
      !$omp end parallel workshare
    end if

    allocate(sort_sglo(ns_unknow), sort_cglo(nc_unknow))
    !$omp parallel workshare
    sort_sglo(:) = cals_glob(:) ; sort_cglo(:) = calc_glob(:)
    !$omp end parallel workshare

    call iquick_sort(sort_sglo, 1, ns_unknow)
    call iquick_sort(sort_cglo, 1, nc_unknow)
    deallocate(cals_glob, calc_glob)
    allocate(nsun_num(ns_unknow), loc2unk_ij(ncals), glo2unk_ij(nxy))
    allocate(ncun_num(nc_unknow), loc2unk_ijk(ncalc), glo2unk_ijk(st_grid%nxyz))

    !$omp parallel
    !$omp do private(i)
    do i = 1, ns_unknow
      nsun_num(i) = i
    end do
    do i = 1, nc_unknow
      ncun_num(i) = i
    end do
    !$omp end do
    !$omp workshare
    loc2unk_ij(:) = 0 ; glo2unk_ij(:) = 0 ; loc2unk_ijk(:) = 0 ; glo2unk_ijk(:) = 0
    glo2unk_ij((/ sort_sglo /)) = nsun_num(:)
    glo2unk_ijk((/ sort_cglo /)) = ncun_num(:)
    loc2unk_ij(:) = pack(glo2unk_ij(:), (glo2loc_ij(:) > 0 .and. glo2loc_ij(:) <= ncals))
    loc2unk_ijk(:) = pack(glo2unk_ijk(:), (glo2loc_ijk(:) > 0 .and. glo2loc_ijk(:) <= ncalc))
    !$omp end workshare
    !$omp end parallel
    deallocate(sort_sglo, nsun_num, glo2unk_ij)
    deallocate(sort_cglo, ncun_num, glo2unk_ijk)

    ! -- Set calculation view (calc_view)
      call set_calc_view(ncals, ncalc, loc2glo_ij, loc2glo_ijk)
    ! -- Set seal view (seal_view)
      call set_seal_view(no_ncals, no_ncalc, loc2glo_nos, loc2glo_noc)
    ! -- Set restart view (rest_view)
      call set_rest_view(ncalc, nc_unknow, loc2unk_ijk)
    ! -- Set write file view (write_fview)
      call set_write_fview(ncals, ncalc, loc2glo_ij, loc2glo_ijk, loc2glo_nos, loc2glo_noc)
#ifdef ICI
    deallocate(loc2unk_ijk)
#else
    deallocate(loc2unk_ij, loc2unk_ijk)
#endif
#endif

    if (noclas_flag /= 1) then
      ! -- Set local cell classification (loc_cell_clas)
        call set_loc_cell_clas()
    end if

    deallocate(glob_clas_flag)

    if (st_out_type%calg == out_type(1)) then
      if (my_rank == 0) then
        call open_new_wtxt(st_out_path%calg, "output calculation grid number", calg_num)
        write(calg_num,'(a)') "Calclation_No,I,J,K"
      end if
      do n = 1, pro_totn
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Send and Receive region information (reg_info)
            call senrec_reg_info(n, loc_nreg, cur_nreg)
          pro_nreg = size(cur_nreg)
        else
          pro_nreg = loc_regn
          allocate(cur_nreg(pro_nreg))
          !$omp parallel workshare
          cur_nreg(:) = loc_nreg(:)
          !$omp end parallel workshare
        end if
#else
        pro_nreg = loc_regn
        allocate(cur_nreg(pro_nreg))
        !$omp parallel workshare
        cur_nreg(:) = loc_nreg(:)
        !$omp end parallel workshare
#endif
        if (my_rank == 0) then
          write(calg_num,'(a,i0)') "Rank", n-1
        end if
        do i = 1, pro_nreg
          if (my_rank == 0) then
            write(calg_num,'(a,i0)') "Region", cur_nreg(i)
          end if
          sta_calc = calc_end(i-1) ; end_calc = calc_end(i)
          do j = sta_calc+1, end_calc
            if (n == my_rank+1) then
              ! -- Get calculation number from grid number (calc_grid)
                call get_calc_grid(j, i_num, j_num, k_num)
            end if
#ifdef MPI_MSG
            if (pro_totn /= 1) then
              ! -- Send and Receive grid number (grid_num)
                call senrec_grid_num(n, i_num, j_num, k_num)
            end if
#endif
            if (my_rank == 0) then
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
      if (my_rank == 0) then
        call close_file(calg_num)
      end if
    end if

    deallocate(calc_end, loc_nreg)

  end subroutine set_cell_info

  subroutine set_glob_cell_clas()
  !***************************************************************************************
  ! set_glob_cell_clas -- Set global cell classification
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, ii, jj, kk, c_num
    integer(I4) :: clasi, clasj, clask
    !-------------------------------------------------------------------------------------
    !$omp parallel workshare
    glob_clas_flag(:,:) = 0
    !$omp end parallel workshare

    do j = 1, st_clas%totn
      do i = 1, st_clas%num(j)
        clasi = st_clas%i(i,j) ; clasj = st_clas%j(i,j) ; clask = st_clas%k(i,j)
        if (clasi == -1 .and. clasj == -1 .and. clask == -1) then !all cells
          glob_clas_flag(1:st_grid%nxyz,j) = 1
        else if (clasi == -1 .and. clasj == -1) then !i,j cells
          do jj = 1, st_grid%ny
            do ii = 1, st_grid%nx
              c_num = st_grid%nx*(st_grid%ny*(clask-1) + (jj-1)) + ii
              glob_clas_flag(c_num,j) = 1
            end do
          end do
        else if (clasi == -1 .and. clask == -1) then !i,k cells
          do kk = 1, st_grid%nz
            do ii = 1, st_grid%nx
              c_num = st_grid%nx*(st_grid%ny*(kk-1) + (clasj-1)) + ii
              glob_clas_flag(c_num,j) = 1
            end do
          end do
        else if (clasj == -1 .and. clask == -1) then !j,k cells
          do kk = 1, st_grid%nz
            do jj = 1, st_grid%ny
              c_num = st_grid%nx*(st_grid%ny*(kk-1) + (jj-1)) + clasi
              glob_clas_flag(c_num,j) = 1
            end do
          end do
        else if (clasi == -1) then !only i cell
          do ii = 1, st_grid%nx
            c_num = st_grid%nx*(st_grid%ny*(clask-1) + (clasj-1)) + ii
            glob_clas_flag(c_num,j) = 1
          end do
        else if (clasj == -1) then !only j cell
          do jj = 1, st_grid%ny
            c_num = st_grid%nx*(st_grid%ny*(clask-1) + (jj-1)) + clasi
            glob_clas_flag(c_num,j) = 1
          end do
        else if (clask == -1) then !only k cell
          do kk = 1, st_grid%nz
            c_num = st_grid%nx*(st_grid%ny*(kk-1) + (clasj-1)) + clasi
            glob_clas_flag(c_num,j) = 1
          end do
        else !others
          c_num = st_grid%nx*(st_grid%ny*(clask-1) + (clasj-1)) + clasi
          glob_clas_flag(c_num,j) = 1
        end if
      end do
    end do

    deallocate(st_clas%i, st_clas%j, st_clas%k, st_clas%num)

  end subroutine set_glob_cell_clas

  subroutine set_glob_reg()
  !***************************************************************************************
  ! set_glob_reg -- Set global region
  !***************************************************************************************
    ! -- modules
    use constval_module, only: VARLEN
    use utility_module, only: open_new_rtxt, open_new_rbin, write_err_stop
    use initial_module, only: in_type
    use read_module, only: read_2d_calcreg, read_3d_calcreg
    ! -- inout

    ! -- local
    integer(I4) :: i, j
    integer(I4) :: reg_type, glob_ncalc, calreg_fnum
    integer(I4) :: len_char, first_pos, sta_pos, end_pos, temp_num
    integer(I4), allocatable :: temp_reg(:), type_2d(:), type_3d(:)
    character(1), allocatable :: temp_char(:)
    character(VARLEN), allocatable :: reg_name(:)
    logical, allocatable :: mask(:)
    !-------------------------------------------------------------------------------------
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
      allocate(mask(st_grid%nxyz))
      do j = 1, st_clas%totn
        if (trim(adjustl(st_sim%inact_name)) == trim(adjustl(st_clas%name(j)))) then
          temp_num = sum(glob_clas_flag(:,j))
          allocate(temp_reg(temp_num))
          !$omp parallel workshare
          temp_reg(:) = -1 ; mask(:) = (glob_clas_flag(:,j) == 1)
          !$omp end parallel workshare
          glob_reg_flag(:) = unpack(temp_reg(:), mask(:), glob_reg_flag(:))
          deallocate(temp_reg)
        end if
        temp_num = 0
        do i = 1, totnreg
          if (trim(adjustl(reg_name(i))) == trim(adjustl(st_clas%name(j)))) then
            temp_num = sum(glob_clas_flag(:,j))
            allocate(temp_reg(temp_num))
            !$omp parallel workshare
            temp_reg(:) = i ; mask(:) = (glob_clas_flag(:,j) == 1)
            !$omp end parallel workshare
            glob_reg_flag(:) = unpack(temp_reg(:), mask(:), glob_reg_flag(:))
            deallocate(temp_reg)
          end if
        end do
        glob_ncalc = glob_ncalc + temp_num
      end do
      deallocate(mask, reg_name)
      if (.not.(any(glob_reg_flag(:) == -1))) then
        call write_err_stop("check the no calculation region name.")
      end if

    else if (any(reg_type == type_2d(:)) .or. any(reg_type == type_3d(:))) then
      if (reg_type == type_2d(1) .or. reg_type == type_3d(1)) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, st_sim%reg_name, "calculation reigion", calreg_fnum)
      else if (reg_type == type_2d(2) .or. reg_type == type_3d(2)) then
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, st_sim%reg_name, "calculation reigion", calreg_fnum)
      end if
      if (any(reg_type == type_2d(:))) then
        ! -- Read calcluation region file from 2d file (2d_calcreg)
          call read_2d_calcreg(calreg_fnum, reg_type, glob_reg_flag, glob_ncalc)
      else if (any(reg_type == type_3d(:))) then
        ! -- Read calcluation region file from 3d file (3d_calcreg)
          call read_3d_calcreg(calreg_fnum, reg_type, glob_reg_flag, glob_ncalc)
      end if
      totnreg = maxval(glob_reg_flag(:))
      if (.not.(any(glob_reg_flag(:) == -1))) then
        call write_err_stop("check the no calculation region name.")
      end if
    end if

    deallocate(type_2d, type_3d)

  end subroutine set_glob_reg

  subroutine div_calc_reg_2d()
  !***************************************************************************************
  ! div_calc_reg_2d -- Divide calculation region for 2d
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, n, nxy
    integer(I4) :: rough_divn, sum_cals, sum_mpi, pro_num, quot, max_reg, pre_cals
    integer(I4) :: reg_mpi_ncals, reg_mpi_remain, sta_mpi, end_mpi, mpi
    integer(I4), allocatable :: reg_ncals(:), reg_remain(:)
    integer(I4), allocatable :: reg_mpi_num(:), reg_mpi_end(:)
    integer(I4), allocatable :: reg_num(:), grid_num(:), reg_num_mpi(:)
!    logical, allocatable :: mask(:)
    !-------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    allocate(reg_ncals(totnreg), reg_remain(totnreg))
    allocate(reg_mpi_num(totnreg), reg_mpi_end(0:totnreg))
    !$omp parallel
    !$omp workshare
    reg_ncals(:) = 0 ; reg_remain(:) = 0 ; reg_mpi_num(:) = 0 ; reg_mpi_end(:) = 0
    !$omp end workshare
    !$omp do private(i)
    do i = 1, totnreg
      reg_ncals(i) = count(glob_reg_flag(1:nxy) == i)
    end do
    !$omp end do
    !$omp end parallel

    rough_divn = sum(reg_ncals(:))/pro_totn
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
    if (reg_mpi_end(totnreg) < pro_totn) then
      do j = 1, pro_totn-reg_mpi_end(totnreg)
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

    allocate(grid_num(nxy))
    !$omp parallel do private(i)
    do i = 1, nxy
      grid_num(i) = i
    end do
    !$omp end parallel do
    do i = 1, totnreg
      reg_mpi_ncals = reg_ncals(i)/reg_mpi_num(i) ; sum_cals = 0 ; pre_cals = 0
      reg_mpi_remain = mod(reg_ncals(i), reg_mpi_num(i))
      sta_mpi = reg_mpi_end(i-1) ; end_mpi = reg_mpi_end(i)
      if (end_mpi == sta_mpi) then
        sta_mpi = sta_mpi - 1
      end if
      allocate(reg_num(reg_ncals(i)))
      reg_num(:) = pack(grid_num(:), glob_reg_flag(1:nxy) == i)
      do n = 1, end_mpi-sta_mpi
        if (n <= reg_mpi_remain) then
          sum_cals = pre_cals + reg_mpi_ncals + 1
        else
          sum_cals = pre_cals + reg_mpi_ncals
        end if
        mpi = sta_mpi + n
        allocate(reg_num_mpi(sum_cals-pre_cals))
        reg_num_mpi(:) = reg_num(pre_cals+1:sum_cals)
        glob_mpi_flag((/ reg_num_mpi /)) = mpi
        pre_cals = sum_cals
        deallocate(reg_num_mpi)
      end do
      deallocate(reg_num)
    end do

    deallocate(reg_ncals, reg_mpi_num, reg_mpi_end, grid_num)

    !$omp parallel do private(i)
    do i = 1, st_grid%nz-1
      glob_mpi_flag(nxy*i+1:nxy*(i+1)) = glob_mpi_flag(1:nxy)
    end do
    !$omp end parallel do

  end subroutine div_calc_reg_2d

!  subroutine div_calc_reg_3d()
!  !***************************************************************************************
!  ! div_calc_reg_3d -- Divide calculation region for 3d
!  !***************************************************************************************
!    ! -- modules
!
!    ! -- inout
!
!    ! -- local
!    integer(I4) :: i, j, k, n
!    integer(I4) :: rough_divn, sum_calc, sum_mpi, pro_num, quot, max_reg, pre_calc
!    integer(I4) :: reg_mpi_ncalc, reg_mpi_remain, sta_mpi, end_mpi, mpi
!    integer(I4), allocatable :: reg_ncalc(:), reg_remain(:)
!    integer(I4), allocatable :: reg_mpi_num(:), reg_mpi_end(:)
!    integer(I4), allocatable :: reg_num(:), grid_num(:), reg_num_mpi(:)
!    !-------------------------------------------------------------------------------------
!    allocate(reg_ncalc(totnreg), reg_remain(totnreg))
!    allocate(reg_mpi_num(totnreg), reg_mpi_end(0:totnreg))
!    !$omp parallel
!    !$omp workshare
!    reg_ncalc(:) = 0 ; reg_remain(:) = 0 ; reg_mpi_num(:) = 0 ; reg_mpi_end(:) = 0
!    !$omp end workshare
!    !$omp do private(i)
!    do i = 1, totnreg
!      reg_ncalc(i) = count(glob_reg_flag(:) == i)
!    end do
!    !$omp end do
!    !$omp end parallel
!
!    rough_divn = sum(reg_ncalc(:))/pro_totn
!    if (rough_divn == 0) then
!      rough_divn = 1
!    end if
!
!    sum_calc = 0 ; sum_mpi = 0
!    do i = 1, totnreg
!      sum_calc = sum_calc + reg_ncalc(i)
!      quot = sum_calc/rough_divn ; reg_remain(i) = mod(sum_calc, rough_divn)
!      if (quot == 0) then
!        reg_mpi_num(i) = 1 ; reg_remain(i) = 0 ; reg_mpi_end(i) = sum_mpi + 1
!      else if (sum_calc > reg_ncalc(i) .and. i /= totnreg) then
!        reg_mpi_end(i) = sum_mpi + 1 ; sum_calc = 0
!      else
!        reg_mpi_num(i) = quot ; sum_mpi = sum_mpi + quot ; reg_mpi_end(i) = sum_mpi
!        sum_calc = 0
!      end if
!    end do
!
!    pro_num = 0
!    if (reg_mpi_end(totnreg) < pro_totn) then
!      do j = 1, pro_totn-reg_mpi_end(totnreg)
!        max_reg = maxloc(reg_remain(:),1) ;
!        reg_mpi_num(max_reg) = reg_mpi_num(max_reg) + 1
!        if (reg_remain(max_reg) < reg_ncalc(max_reg)/2) then
!          pro_num = max_reg + 1
!        else
!          pro_num = max_reg
!        end if
!        reg_remain(max_reg) = 0
!        do k = pro_num, totnreg
!          reg_mpi_end(k) = reg_mpi_end(k) + 1
!        end do
!      end do
!    end if
!    deallocate(reg_remain)
!
!    allocate(grid_num(st_grid%nxyz))
!    !$omp parallel do private(i)
!    do i = 1, st_grid%nxyz
!      grid_num(i) = i
!    end do
!    !$omp end parallel do
!    do i = 1, totnreg
!      reg_mpi_ncalc = reg_ncalc(i)/reg_mpi_num(i) ; sum_calc = 0 ; pre_calc = 0
!      reg_mpi_remain = mod(reg_ncalc(i), reg_mpi_num(i))
!      sta_mpi = reg_mpi_end(i-1) ; end_mpi = reg_mpi_end(i)
!      allocate(reg_num(count(glob_reg_flag(:) == i)))
!      reg_num(:) = pack(grid_num(:), glob_reg_flag(:) == i)
!      if (end_mpi == sta_mpi) then
!        sta_mpi = sta_mpi - 1
!      end if
!      do n = 1, end_mpi-sta_mpi
!        if (n <= reg_mpi_remain) then
!          sum_calc = pre_calc + reg_mpi_ncalc + 1
!        else
!          sum_calc = pre_calc + reg_mpi_ncalc
!        end if
!        mpi = sta_mpi + n
!        allocate(reg_num_mpi(sum_calc-pre_calc))
!        reg_num_mpi(:) = reg_num(pre_calc+1:sum_calc)
!        glob_mpi_flag((/ reg_num_mpi /)) = mpi
!        pre_calc = sum_calc
!        deallocate(reg_num_mpi)
!      end do
!      deallocate(reg_num)
!    end do
!
!    deallocate(reg_ncalc, reg_mpi_num, reg_mpi_end, grid_num)
!
!  end subroutine div_calc_reg_3d

  subroutine div_nocalc_flag_2d()
  !***************************************************************************************
  ! div_nocalc_flag_2d -- Divide no calculation flag for 2d
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, k, nxy, nxyz0, nxyz1
    integer(I4) :: glo_nocals, rough_divn, nocals_remain
    integer(I4) :: nocals_sta, nocals_end
    integer(I4), allocatable :: nocals_mpi_num(:), grid_num(:), nocals_glo_num(:)
    integer(I4), allocatable :: nocals_mpi_glo(:)
    !-------------------------------------------------------------------------------------
    nxy = st_grid%nx*st_grid%ny
    do k = 1, st_grid%nz
      nxyz0 = nxy*(k-1)+1 ; nxyz1 = nxy*k
      glo_nocals = count(glob_mpi_flag(nxyz0:nxyz1) == 0)
      rough_divn = glo_nocals/pro_totn ; nocals_remain = mod(glo_nocals, pro_totn)
      allocate(nocals_mpi_num(pro_totn), grid_num(nxy), nocals_glo_num(glo_nocals))
      !$omp parallel
      !$omp workshare
      nocals_mpi_num(:) = rough_divn
      !$omp end workshare
      !$omp do private(i)
      do i = 1, pro_totn
        if (i <= nocals_remain) then
          nocals_mpi_num(i) = nocals_mpi_num(i) + 1
        end if
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, nxy
        grid_num(i) = nxy*(k-1) + i
      end do
      !$omp end do
      !$omp workshare
      nocals_glo_num(:) = pack(grid_num(:), glob_mpi_flag(nxyz0:nxyz1) == 0)
      !$omp end workshare
      !$omp end parallel
      deallocate(grid_num)

      nocals_sta = 0 ; nocals_end = 0
      do i = 1, pro_totn
        nocals_end = nocals_sta + nocals_mpi_num(i)
        allocate(nocals_mpi_glo(nocals_end-nocals_sta))
        nocals_mpi_glo(:) = nocals_glo_num(nocals_sta+1:nocals_end)
        glob_mpi_flag((/ nocals_mpi_glo /)) = -i
        nocals_sta = nocals_end
        deallocate(nocals_mpi_glo)
      end do

      deallocate(nocals_mpi_num, nocals_glo_num)

    end do

  end subroutine div_nocalc_flag_2d

!  subroutine div_nocalc_flag_3d()
!  !***************************************************************************************
!  ! div_nocalc_flag_3d -- Divide no calculation flag for 3d
!  !***************************************************************************************
!    ! -- modules
!
!    ! -- inout
!
!    ! -- local
!    integer(I4) :: i, nxyz
!    integer(I4) :: glo_nocalc, rough_divn, nocalc_remain
!    integer(I4) :: nocalc_sta, nocalc_end
!    integer(I4), allocatable :: nocalc_mpi_num(:), grid_num(:), nocalc_glo_num(:)
!    integer(I4), allocatable :: nocalc_mpi_glo(:)
!    !-------------------------------------------------------------------------------------
!    nxyz = st_grid%nxyz
!    glo_nocalc = count(glob_mpi_flag(:) == 0)
!    rough_divn = glo_nocalc/pro_totn ; nocalc_remain = mod(glo_nocalc, pro_totn)
!    allocate(nocalc_mpi_num(pro_totn), grid_num(nxyz), nocalc_glo_num(glo_nocalc))
!    !$omp parallel
!    !$omp workshare
!    nocalc_mpi_num(:) = rough_divn
!    !$omp end workshare
!    !$omp do private(i)
!    do i = 1, pro_totn
!      if (i <= nocalc_remain) then
!        nocalc_mpi_num(i) = nocalc_mpi_num(i) + 1
!      end if
!    end do
!    !$omp end do
!    !$omp do private(i)
!    do i = 1, nxyz
!      grid_num(i) = i
!    end do
!    !$omp end do
!    !$omp workshare
!    nocalc_glo_num(:) = pack(grid_num(:), glob_mpi_flag(:) == 0)
!    !$omp end workshare
!    !$omp end parallel
!
!    deallocate(grid_num)
!
!    nocalc_sta = 0 ; nocalc_end = 0
!    do i = 1, pro_totn
!      nocalc_end = nocalc_sta + nocalc_mpi_num(i)
!      allocate(nocalc_mpi_glo(nocalc_end-nocalc_sta))
!      nocalc_mpi_glo(:) = nocalc_glo_num(nocalc_sta+1:nocalc_end)
!      glob_mpi_flag((/ nocalc_mpi_glo /)) = -i
!      nocalc_sta = nocalc_end
!      deallocate(nocalc_mpi_glo)
!    end do
!
!    deallocate(nocalc_mpi_num, nocalc_glo_num)
!
!  end subroutine div_nocalc_flag_3d

  subroutine set_rel_gloloc()
  !***************************************************************************************
  ! set_rel_gloloc -- Set relationship global&local
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, ij, count_reg, reg_flag, temp_nreg
    integer(I4) :: count_calc, count_cals
    integer(I4), allocatable :: temp_mpi_reg(:), temp_cend(:)
    logical, allocatable :: mask(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_mpi_reg(ncalc), mask(ncalc))
    temp_mpi_reg(:) = pack(glob_reg_flag(:), glob_mpi_flag(:) == my_rank+1)
    count_reg = 0
    !$omp parallel do private(i), reduction(+:count_reg)
    do i = 1, totnreg
      mask(:) = (temp_mpi_reg(:) == i)
      if (any(mask)) then
        count_reg = count_reg + 1
      end if
    end do
    !$omp end parallel do
    loc_regn = count_reg

    deallocate(temp_mpi_reg, mask)

    allocate(loc_nreg(loc_regn), temp_cend(0:loc_regn))
    allocate(calc2reg(ncalc))
    !$omp parallel workshare
    loc_nreg(:) = 0 ; temp_cend(:) = 0 ; calc2reg(:) = 0
    !$omp end parallel workshare

    count_reg = 0 ; count_calc = 0 ; count_cals = 0
    do k = 1, st_grid%nz
      do i = 1, st_grid%nx*st_grid%ny
        ij = (k-1)*st_grid%nx*st_grid%ny + i
        if (glob_mpi_flag(ij) == my_rank+1) then
          reg_flag = 0
          if (count_reg /= 0) then
            do j = 1, count_reg
              if (loc_nreg(j) == glob_reg_flag(ij)) then
                reg_flag = 1 ; temp_nreg = j
              end if
            end do
            if (reg_flag == 0) then
              count_reg = count_reg + 1 ; loc_nreg(count_reg) = glob_reg_flag(ij)
              temp_nreg = count_reg
            end if
          else
            count_reg = 1 ; loc_nreg(1) = glob_reg_flag(ij) ; temp_nreg = count_reg
          end if
          if (k == 1) then
            count_cals = count_cals + 1
            l2g_ij(count_cals) = i ; glo2loc_ij(i) = count_cals
          end if
          count_calc = count_calc + 1 ; temp_cend(temp_nreg) = temp_cend(temp_nreg) + 1
          l2g_ijk(count_calc) = ij ; glo2loc_ijk(ij) = count_calc
          calc2reg(count_calc) = temp_nreg
        end if
      end do
    end do

    allocate(calc_end(0:loc_regn))
    !$omp parallel
    !$omp workshare
    calc_end(:) = 0
    !$omp end workshare
    !$omp do private(i)
    do i = 1, loc_regn
      calc_end(i) = calc_end(i-1) + temp_cend(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(temp_cend)

  end subroutine set_rel_gloloc

#ifdef MPI_MSG
  subroutine set_mpi_rel()
  !***************************************************************************************
  ! set_mpi_rel -- Set mpi relationship
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, ij
    integer(I4) :: i_num, j_num, k_num, g_num, n_num, w_num, e_num, s_num, u_num, d_num
    integer(I4) :: neib_mpi_num, sta_send, end_send, sendrecv_num, send_loc, recv_loc
    integer(I4) :: cals_niebn, calc_niebn, reg_num
    integer(I4), allocatable :: temp_neib_num(:), temp_neib_flag(:), temp_mpi_num(:)
    integer(I4), allocatable :: neib_glos(:), neib_gloc(:), neib_locc(:)
    integer(I4), allocatable :: sort_recv_num(:), temp_calc_reg(:)
    integer(I4), allocatable :: recv_num(:), send_num(:)
    integer(I4), allocatable :: temp_sort(:), sort_mpi_num(:), loc_send_num(:)
    integer(I4), allocatable :: mpi_l2g_ij(:), mpi_l2g_ijk(:), mpi_calc2reg(:)
    !-------------------------------------------------------------------------------------
    if (pro_totn /= 1) then
      allocate(neib_glos(ncals))
      !$omp parallel workshare
      neib_glos(:) = 0
      !$omp end parallel workshare
      do i = 1, ncals
        g_num = l2g_ij(i)
        k_num = (g_num-1)/(st_grid%nx*st_grid%ny) + 1
        ij = g_num - st_grid%nx*st_grid%ny*(k_num-1)
        j_num = (ij-1)/st_grid%nx + 1
        i_num = ij - (j_num-1)*st_grid%nx
        ! north direction
        if (j_num /= 1) then
          n_num = g_num-st_grid%nx ; reg_num = glob_reg_flag(n_num)
          if (reg_num > 0 .and. glob_mpi_flag(n_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncals = neib_ncals + 1 ; neib_glos(neib_ncals) = n_num
            end if
          end if
        end if
        ! west direction
        if (i_num /= 1) then
          w_num = g_num-1 ; reg_num = glob_reg_flag(w_num)
          if (reg_num > 0 .and. glob_mpi_flag(w_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncals = neib_ncals + 1 ; neib_glos(neib_ncals) = w_num
            end if
          end if
        end if
        ! east direction
        if (i_num /= st_grid%nx) then
          e_num = g_num+1 ; reg_num = glob_reg_flag(e_num)
          if (reg_num > 0 .and. glob_mpi_flag(e_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncals = neib_ncals + 1 ; neib_glos(neib_ncals) = e_num
            end if
          end if
        end if
        ! south direction
        if (j_num /= st_grid%ny) then
          s_num = g_num+st_grid%nx ; reg_num = glob_reg_flag(s_num)
          if (reg_num > 0 .and. glob_mpi_flag(s_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncals = neib_ncals + 1 ; neib_glos(neib_ncals) = s_num
            end if
          end if
        end if
      end do

      allocate(temp_neib_num(pro_totn), temp_neib_flag(pro_totn))
      allocate(temp_mpi_num(ncalc), neib_locc(ncalc))
      allocate(neib_gloc(ncalc), temp_calc_reg(ncalc))
      !$omp parallel workshare
      temp_neib_num(:) = 0 ; temp_neib_flag(:) = 0 ; temp_mpi_num(:) = 0
      neib_locc(:) = 0 ; neib_gloc(:) = 0 ; temp_calc_reg(:) = 0
      !$omp end parallel workshare

      do i = 1, ncalc
        g_num = l2g_ijk(i)
        k_num = (g_num-1)/(st_grid%nx*st_grid%ny) + 1
        ij = g_num - st_grid%nx*st_grid%ny*(k_num-1)
        j_num = (ij-1)/st_grid%nx + 1
        i_num = ij - (j_num-1)*st_grid%nx
        ! up direction
        if (k_num /= 1) then
          u_num = g_num-st_grid%nx*st_grid%ny ; reg_num = glob_reg_flag(u_num)
          if (reg_num > 0 .and. glob_mpi_flag(u_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1 ; neib_mpi_num = glob_mpi_flag(u_num)
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
          n_num = g_num-st_grid%nx ; reg_num = glob_reg_flag(n_num)
          if (reg_num > 0 .and. glob_mpi_flag(n_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1 ; neib_mpi_num = glob_mpi_flag(n_num)
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
          w_num = g_num-1 ; reg_num = glob_reg_flag(w_num)
          if (reg_num > 0 .and. glob_mpi_flag(w_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1 ; neib_mpi_num = glob_mpi_flag(w_num)
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
          e_num = g_num+1 ; reg_num = glob_reg_flag(e_num)
          if (reg_num > 0 .and. glob_mpi_flag(e_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1 ; neib_mpi_num = glob_mpi_flag(e_num)
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
          s_num = g_num+st_grid%nx ; reg_num = glob_reg_flag(s_num)
          if (reg_num > 0 .and. glob_mpi_flag(s_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1 ; neib_mpi_num = glob_mpi_flag(s_num)
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
          d_num = g_num+st_grid%nx*st_grid%ny ; reg_num = glob_reg_flag(d_num)
          if (reg_num > 0 .and. glob_mpi_flag(d_num) /= my_rank+1) then
            if (reg_num == glob_reg_flag(g_num) .or. st_sim%reg_neib == 1) then
              neib_ncalc = neib_ncalc + 1 ; neib_mpi_num = glob_mpi_flag(d_num)
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
      !$omp parallel workshare
      mpi_l2g_ij(1:ncals) = l2g_ij(:) ; mpi_l2g_ijk(1:ncalc) = l2g_ijk(:)
      mpi_calc2reg(1:ncalc) = calc2reg(:)
      !$omp end parallel workshare
      deallocate(l2g_ij, l2g_ijk, calc2reg)

      !$omp parallel
      !$omp do private(i, cals_niebn)
      do i = 1, neib_ncals
        cals_niebn = ncals + i
        mpi_l2g_ij(cals_niebn) = neib_glos(i)
        glo2loc_ij(neib_glos(i)) = cals_niebn
      end do
      !$omp end do
      !$omp do private(i, calc_niebn)
      do i = 1, neib_ncalc
        calc_niebn = ncalc + i
        mpi_l2g_ijk(calc_niebn) = neib_gloc(i)
        mpi_calc2reg(calc_niebn) = temp_calc_reg(i)
        glo2loc_ijk(neib_gloc(i)) = calc_niebn
      end do
      !$omp end do
      !$omp end parallel

      deallocate(neib_glos, temp_calc_reg)

      allocate(l2g_ij(ncals+neib_ncals), l2g_ijk(ncalc+neib_ncalc))
      allocate(calc2reg(ncalc+neib_ncalc))
      !$omp parallel workshare
      l2g_ij(:) = mpi_l2g_ij(:) ; l2g_ijk(:) = mpi_l2g_ijk(:)
      calc2reg(:) = mpi_calc2reg(:)
      !$omp end parallel workshare

      deallocate(mpi_l2g_ij, mpi_l2g_ijk, mpi_calc2reg)

      sendrecv_num = count(temp_mpi_num(:) /= 0)
      allocate(neib_num(neib_mpi_totn))
      allocate(send_cind(0:neib_mpi_totn), recv_cind(0:neib_mpi_totn))
      allocate(send_citem(sendrecv_num), recv_citem(sendrecv_num))
      allocate(temp_sort(neib_ncalc), sort_recv_num(neib_ncalc), sort_mpi_num(neib_ncalc))
      allocate(send2recv(neib_ncalc), loc_send_num(neib_ncalc))
      !$omp parallel workshare
      send_cind(:) = 0 ; recv_cind(:) = 0
      send_citem(:) = 0 ; recv_citem(:) = 0
      temp_sort(:) = neib_gloc(1:neib_ncalc) ; sort_recv_num(:) = 0 ; sort_mpi_num = 0
      send2recv(:) = 0 ; loc_send_num(:) = neib_locc(1:neib_ncalc)
      !$omp end parallel workshare

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
  !***************************************************************************************
  ! set_rel_seareg -- Set relationship of sea region
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, ij, nxy
    integer(I4) :: i_num, j_num, k_num, g_num, n_num, w_num, e_num, s_num, u_num, d_num
    !-------------------------------------------------------------------------------------
    seal_snum = 0
    nxy = st_grid%nx*st_grid%ny

    do i = 1, ncals
      g_num = l2g_ij(i)
      k_num = (g_num-1)/(nxy) + 1
      ij = g_num - nxy*(k_num-1)
      j_num = (ij-1)/st_grid%nx + 1
      i_num = ij - (j_num-1)*st_grid%nx
      ! north direction
      if (j_num /= 1) then
        n_num = g_num-st_grid%nx
        if (glob_reg_flag(n_num) == 0 .and. glo2loc_ij(n_num) == 0) then
          seal_snum = seal_snum + 1 ; glo2loc_ij(n_num) = ncals + neib_ncals + seal_snum
        end if
      end if
      ! west direction
      if (i_num /= 1) then
        w_num = g_num-1
        if (glob_reg_flag(w_num) == 0 .and. glo2loc_ij(w_num) == 0) then
          seal_snum = seal_snum + 1 ; glo2loc_ij(w_num) = ncals + neib_ncals + seal_snum
        end if
      end if
      ! east direction
      if (i_num /= st_grid%nx) then
        e_num = g_num+1
        if (glob_reg_flag(e_num) == 0 .and. glo2loc_ij(e_num) == 0) then
          seal_snum = seal_snum + 1 ; glo2loc_ij(e_num) = ncals + neib_ncals + seal_snum
        end if
      end if
      ! south direction
      if (j_num /= st_grid%ny) then
        s_num = g_num+st_grid%nx
        if (glob_reg_flag(s_num) == 0 .and. glo2loc_ij(s_num) == 0) then
          seal_snum = seal_snum + 1 ; glo2loc_ij(s_num) = ncals + neib_ncals + seal_snum
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
        if (glob_reg_flag(u_num) == 0 .and. glo2loc_ijk(u_num) == 0) then
          seal_cnum = seal_cnum + 1 ; glo2loc_ijk(u_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! north direction
      if (j_num /= 1) then
        n_num = g_num-st_grid%nx
        if (glob_reg_flag(n_num) == 0 .and. glo2loc_ijk(n_num) == 0) then
          seal_cnum = seal_cnum + 1 ; glo2loc_ijk(n_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! west direction
      if (i_num /= 1) then
        w_num = g_num-1
        if (glob_reg_flag(w_num) == 0 .and. glo2loc_ijk(w_num) == 0) then
          seal_cnum = seal_cnum + 1 ; glo2loc_ijk(w_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! east direction
      if (i_num /= st_grid%nx) then
        e_num = g_num+1
        if (glob_reg_flag(e_num) == 0 .and. glo2loc_ijk(e_num) == 0) then
          seal_cnum = seal_cnum + 1 ; glo2loc_ijk(e_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! south direction
      if (j_num /= st_grid%ny) then
        s_num = g_num+st_grid%nx
        if (glob_reg_flag(s_num) == 0 .and. glo2loc_ijk(s_num) == 0) then
          seal_cnum = seal_cnum + 1 ; glo2loc_ijk(s_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
      ! down direction
      if (k_num /= st_grid%nz) then
        d_num = g_num+nxy
        if (glob_reg_flag(d_num) == 0 .and. glo2loc_ijk(d_num) == 0) then
          seal_cnum = seal_cnum + 1 ; glo2loc_ijk(d_num) = ncalc + neib_ncalc + seal_cnum
        end if
      end if
    end do

  end subroutine set_rel_seareg

  subroutine set_loc_cell_clas()
  !***************************************************************************************
  ! set_loc_cell_clas -- Set local cell classification
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    allocate(clas_flag(ncell,st_clas%totn))
    !$omp parallel
    !$omp workshare
    clas_flag(:,:) = 0
    !$omp end workshare
    !$omp do private(i)
    do i = 1, st_grid%nxyz
      if (glo2loc_ijk(i) > 0) then
        clas_flag(glo2loc_ijk(i),:) = glob_clas_flag(i,:)
      end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine set_loc_cell_clas

  subroutine get_cals_grid(cal_num, x_num, y_num)
  !***************************************************************************************
  ! get_cals_grid -- Get surface number from grid number
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: cal_num
    integer(I4), intent(out) :: x_num, y_num
    ! -- local
    integer(I4) :: s_num
    !-------------------------------------------------------------------------------------
    s_num = loc2glo_ij(cal_num)
    y_num = (s_num-1)/st_grid%nx + 1
    x_num = s_num - (y_num-1)*st_grid%nx

  end subroutine get_cals_grid

  subroutine get_calc_grid(cal_num, x_num, y_num, z_num)
  !***************************************************************************************
  ! get_calc_grid -- Get calculation number from grid number
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: cal_num
    integer(I4), intent(out) :: x_num, y_num, z_num
    ! -- local
    integer(I4) :: c_num, xy_num
    !-------------------------------------------------------------------------------------
    c_num = loc2glo_ijk(cal_num)
    z_num = (c_num-1)/(st_grid%nx*st_grid%ny) + 1
    xy_num = c_num - st_grid%nx*st_grid%ny*(z_num-1)
    y_num = (xy_num-1)/st_grid%nx + 1
    x_num = xy_num - (y_num-1)*st_grid%nx

  end subroutine get_calc_grid

end module set_cell
