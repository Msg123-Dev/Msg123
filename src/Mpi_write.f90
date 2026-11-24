module mpi_write
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SNOVAL, DZERO
  use utility_module, only: st_mpi, write_err_write, write_err_stop
  use mpi_utility, only: mpimax_val
  use initial_module, only: st_grid
  use set_cell, only: get_calc_grid, ncalc, ncals, neib_mpi_totn, st_conn
  use make_cell, only: st_geom
  use mpi

  implicit none
  private
  public :: write_mpi_2dbin, write_mpi_3dbin, write_mpi_rest
  public :: redu_mpi_mass, set_senrec_wtab, calc_mpi_wtable

  ! -- local
  integer(I4) :: neib_wtab_stotn, neib_wtab_rtotn
  integer(I4), allocatable :: neib_wtab_snum(:), neib_wtab_rnum(:)
  integer(I4), allocatable :: send_wtab_cind(:), recv_wtab_cind(:)
  integer(I4), allocatable :: send_wtab_citem(:), recv_wtab_glo(:)

  contains

  subroutine write_mpi_2dbin(out_fh, out_totn, calc_num, out_unit, out_val, ntime, out_base)
  !*********************************************************************************************
  ! write_mpi_2dbin -- Write MPI 2d binary
  !*********************************************************************************************
    ! -- module
    use mpi_set, only: write_2d_ind, write_2d_num
    ! -- inout
    integer(I4), intent(in) :: out_fh, out_totn
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit, ntime
    real(DP), intent(in) :: out_val(:)
    real(DP), intent(in), optional :: out_base
    ! -- local
    integer(I4) :: i, s, ierr, max_ierr
    integer(I4) :: all_ncals
    integer(I4), allocatable :: istat(:)
    real(SP), allocatable :: vari_sp(:)
    real(DP) :: base_val
    !-------------------------------------------------------------------------------------------
    base_val = DZERO
    if (present(out_base)) then
      base_val = out_base
    end if
    all_ncals = write_2d_num
    allocate(istat(MPI_STATUS_SIZE))
    allocate(vari_sp(all_ncals))
    !$omp parallel
    !$omp do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, all_ncals
      vari_sp(i) = SNOVAL
    end do
    !$omp end do
    !$omp do private(i, s)
    do i = 1, out_totn
      s = write_2d_ind(calc_num(i)) + 1
      vari_sp(s) = real(out_val(i)*out_unit+base_val, kind=SP)
    end do
    !$omp end do
    !$omp end parallel
    vari_sp(1) = ntime

    ierr = 0
    call MPI_FILE_WRITE(out_fh, vari_sp, all_ncals, MPI_REAL4, istat, ierr)
    ! -- Max value for MPI (val)
      call mpimax_val(ierr, "write mpi file error", max_ierr)
    ierr = max_ierr
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_write(out_fh)
      end if
    end if

    deallocate(istat)
    deallocate(vari_sp)

  end subroutine write_mpi_2dbin

  subroutine write_mpi_3dbin(out_fh, out_totn, calc_num, out_unit, out_val, ntime, out_base)
  !*********************************************************************************************
  ! write_mpi_3dbin -- Write MPI 3d binary
  !*********************************************************************************************
    ! -- module
    use mpi_set, only: write_3d_ind
    use set_cell, only: no_ncalc, seal_cnum
    ! -- inout
    integer(I4), intent(in) :: out_fh, out_totn
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit, ntime
    real(DP), intent(in) :: out_val(:)
    real(DP), intent(in), optional :: out_base
    ! -- local
    integer(I4) :: i, c, ierr, max_ierr
    integer(I4) :: all_ncalc
    integer(I4), allocatable :: istat(:)
    real(SP), allocatable :: vari_sp(:)
    real(DP) :: base_val
    !-------------------------------------------------------------------------------------------
    base_val = DZERO
    if (present(out_base)) then
      base_val = out_base
    end if
    all_ncalc = ncalc + no_ncalc + seal_cnum + 1
    allocate(istat(MPI_STATUS_SIZE))
    allocate(vari_sp(all_ncalc))
    !$omp parallel
    !$omp do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, all_ncalc
      vari_sp(i) = SNOVAL
    end do
    !$omp end do
    !$omp do private(i, c)
    do i = 1, out_totn
      c = write_3d_ind(calc_num(i)) + 1
      vari_sp(c) = real(out_val(i)*out_unit+base_val, kind=SP)
    end do
    !$omp end do
    !$omp end parallel
    vari_sp(1) = ntime

    ierr = 0
    call MPI_FILE_WRITE(out_fh, vari_sp, all_ncalc, MPI_REAL4, istat, ierr)
    ! -- Max value for MPI (val)
      call mpimax_val(ierr, "write mpi file error", max_ierr)
    ierr = max_ierr
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_write(out_fh)
      end if
    end if

    deallocate(istat)
    deallocate(vari_sp)

  end subroutine write_mpi_3dbin

  subroutine write_mpi_rest(out_fh, out_time, out_delt, out_unit, out_val)
  !*********************************************************************************************
  ! write_mpi_rest -- Write mpi restart value
  !*********************************************************************************************
    ! -- modules
    use read_input, only: z_base
    ! -- inout
    integer(I4), intent(in) :: out_fh
    real(DP), intent(in) :: out_time, out_delt
    real(SP), intent(in) :: out_unit
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, ierr, max_ierr
    integer(I4), allocatable :: istat(:)
    real(DP) :: sync_time
    real(DP), allocatable :: out_rest(:)
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------------
    sync_time = out_time
    call MPI_BCAST(sync_time, 1, MPI_REAL8, 0, st_mpi%comm, ierr)
    allocate(istat(MPI_STATUS_SIZE), out_rest(ncalc+2))
    !$omp parallel
    !$omp do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncalc
      out_rest(i+1) = out_val(i)*real(out_unit, kind=DP) + z_base
    end do
    !$omp end do
    !$omp end parallel
    out_rest(1) = sync_time ; out_rest(ncalc+2) = out_delt

    ierr = 0 ; head_dis = 0
    call MPI_FILE_WRITE_AT_ALL(out_fh, head_dis, out_rest, ncalc+2, MPI_REAL8, istat, ierr)
    ! -- Max value for MPI (val)
      call mpimax_val(ierr, "write restart file error", max_ierr)
    ierr = max_ierr
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_write(out_fh)
      end if
    end if

    deallocate(istat, out_rest)

  end subroutine write_mpi_rest

  subroutine redu_mpi_mass(num_mass, inout_st)
  !*********************************************************************************************
  ! redu_mpi_mass -- Calculate output massbalance for MPI
  !*********************************************************************************************
    ! -- module
    use assign_calc, only: msout_tnum
    use types_module, only: msout_set
    ! -- inout
    integer(I4), intent(in) :: num_mass
    type(msout_set), intent(inout) :: inout_st
    ! -- local
    integer(I4) :: i, ierr
    real(DP), allocatable :: mpi_sto(:), mpi_con(:), mpi_sea(:), mpi_wel(:)
    real(DP), allocatable :: mpi_rec(:), mpi_sur(:), mpi_riv(:), mpi_lak(:), mpi_tot(:)
    !-------------------------------------------------------------------------------------------
    allocate(mpi_sto(num_mass), mpi_con(num_mass), mpi_sea(num_mass), mpi_wel(num_mass))
    allocate(mpi_rec(num_mass), mpi_sur(num_mass), mpi_riv(num_mass), mpi_lak(num_mass))
    allocate(mpi_tot(num_mass))
    !$omp parallel do private(i)
    do i = 1, num_mass
      mpi_sto(i) = DZERO ; mpi_con(i) = DZERO ; mpi_sea(i) = DZERO ; mpi_wel(i) = DZERO
      mpi_rec(i) = DZERO ; mpi_sur(i) = DZERO ; mpi_riv(i) = DZERO ; mpi_lak(i) = DZERO
      mpi_tot(i) = DZERO
    end do
    !$omp end parallel do

    ierr = 0
    call MPI_REDUCE(inout_st%sto(1), mpi_sto(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum storage for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%con(1), mpi_con(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum connect flow for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%sea(1), mpi_sea(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum sea discharge for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%wel(1), mpi_wel(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum well pumping for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%rec(1), mpi_rec(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum recharge for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%sur(1), mpi_sur(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum surface runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%riv(1), mpi_riv(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum river runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%lak(1), mpi_lak(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum lake runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%tot(1), mpi_tot(1), num_mass, MPI_REAL8, MPI_SUM, 0, st_mpi%comm,&
                    ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum total volume for massbalance.")
      end if
    end if

    if (st_mpi%rank == 0) then
      !$omp parallel do private(i)
      do i = 1, msout_tnum
        inout_st%sto(i) = mpi_sto(i) ; inout_st%con(i) = mpi_con(i)
        inout_st%sea(i) = mpi_sea(i) ; inout_st%wel(i) = mpi_wel(i)
        inout_st%rec(i) = mpi_rec(i) ; inout_st%sur(i) = mpi_sur(i)
        inout_st%riv(i) = mpi_riv(i) ; inout_st%lak(i) = mpi_lak(i)
        inout_st%tot(i) = mpi_tot(i)
      end do
      !$omp end parallel do
    end if

    deallocate(mpi_sto, mpi_con, mpi_sea, mpi_wel)
    deallocate(mpi_rec, mpi_sur, mpi_riv, mpi_lak)
    deallocate(mpi_tot)

  end subroutine redu_mpi_mass

  subroutine set_senrec_wtab(st_matrix)
  !*********************************************************************************************
  ! set_senrec_wtab -- Set send and receive for water table
  !*********************************************************************************************
    ! -- module
    use types_module, only: matrix_set
    use utility_module, only: iquick_sort, iquick_sort2
    use set_cell, only: recv_cind, neib_num
    ! -- inout
    type(matrix_set), intent(in) :: st_matrix
    ! -- local
    integer(I4) :: i, k, ind, neib_ind, temp_neib_num
    integer(I4) :: temp_snum, temp_rnum, sta_ind, end_ind
    integer(I4), allocatable :: temp_send_cind(:), temp_recv_cind(:), temp_send_glo(:)
    !-------------------------------------------------------------------------------------------
    allocate(temp_send_cind(0:neib_mpi_totn), temp_recv_cind(0:neib_mpi_totn))
    !$omp parallel do private(i)
    do i = 0, neib_mpi_totn
      temp_send_cind(i) = 0 ; temp_recv_cind(i) = 0
    end do
    !$omp end parallel do

    do i = 1, ncalc
      do ind = st_matrix%crs_index(1)%offind(i-1)+1, st_matrix%crs_index(1)%offind(i)
        temp_neib_num = st_matrix%crs_index(1)%offrow(ind)
        if (temp_neib_num > ncalc) then
          neib_ind = count(recv_cind(1:neib_mpi_totn) < temp_neib_num-ncalc) + 1
          if (st_matrix%dir_conn(ind) == 1) then
            temp_send_cind(neib_ind) = temp_send_cind(neib_ind) + 1
          else if (st_matrix%dir_conn(ind) == 6) then
            temp_recv_cind(neib_ind) = temp_recv_cind(neib_ind) + 1
          end if
        end if
      end do
    end do

    temp_snum = 0 ; temp_rnum = 0
    do i = 1, neib_mpi_totn
      k = temp_send_cind(i) ; temp_send_cind(i) = temp_snum ; temp_snum = temp_snum + k
      k = temp_recv_cind(i) ; temp_recv_cind(i) = temp_rnum ; temp_rnum = temp_rnum + k
    end do
    allocate(send_wtab_citem(temp_snum), recv_wtab_glo(temp_rnum), temp_send_glo(temp_snum))

    do i = 1, ncalc
      do ind = st_matrix%crs_index(1)%offind(i-1)+1, st_matrix%crs_index(1)%offind(i)
        temp_neib_num = st_matrix%crs_index(1)%offrow(ind)
        if (temp_neib_num > ncalc) then
          neib_ind = count(recv_cind(1:neib_mpi_totn) < temp_neib_num-ncalc) + 1
          if (st_matrix%dir_conn(ind) == 1) then
            temp_send_cind(neib_ind) = temp_send_cind(neib_ind) + 1
            send_wtab_citem(temp_send_cind(neib_ind)) = i
            temp_send_glo(temp_send_cind(neib_ind)) = st_conn%loc2glo_ijk(i)
          else if (st_matrix%dir_conn(ind) == 6) then
            temp_recv_cind(neib_ind) = temp_recv_cind(neib_ind) + 1
            recv_wtab_glo(temp_recv_cind(neib_ind)) = st_conn%loc2glo_ijk(i)&
                                                      + st_grid%nx*st_grid%ny
          end if
        end if
      end do
    end do

    do i = 1, neib_mpi_totn
      sta_ind = temp_send_cind(i-1) + 1 ; end_ind = temp_send_cind(i)
      call iquick_sort2(temp_send_glo, send_wtab_citem, sta_ind, end_ind)
      sta_ind = temp_recv_cind(i-1) + 1 ; end_ind = temp_recv_cind(i)
      call iquick_sort(recv_wtab_glo, sta_ind, end_ind)
    end do

    neib_wtab_stotn = count(temp_send_cind(1:) /= temp_send_cind(:neib_mpi_totn-1))
    neib_wtab_rtotn = count(temp_recv_cind(1:) /= temp_recv_cind(:neib_mpi_totn-1))
    allocate(neib_wtab_snum(neib_wtab_stotn), neib_wtab_rnum(neib_wtab_rtotn))
    allocate(send_wtab_cind(0:neib_wtab_stotn), recv_wtab_cind(0:neib_wtab_rtotn))
    send_wtab_cind(0) = 0 ; recv_wtab_cind(0) = 0

    temp_snum = 0 ; temp_rnum = 0
    do i = 1, neib_mpi_totn
      if (temp_send_cind(i) /= temp_send_cind(i-1)) then
        temp_snum = temp_snum + 1
        neib_wtab_snum(temp_snum) = neib_num(i) ; send_wtab_cind(temp_snum) = temp_send_cind(i)
      end if
      if (temp_recv_cind(i) /= temp_recv_cind(i-1)) then
        temp_rnum = temp_rnum + 1
        neib_wtab_rnum(temp_rnum) = neib_num(i) ; recv_wtab_cind(temp_rnum) = temp_recv_cind(i)
      end if
    end do

  end subroutine set_senrec_wtab

  subroutine calc_mpi_wtable(hnew)
  !*********************************************************************************************
  ! calc_mpi_wtable -- Calculate water table for MPI
  !*********************************************************************************************
    ! -- module
    use mpi_utility, only: mpisum_val
    use set_cell, only: get_cals_grid
    use allocate_output, only: wtable
    ! -- inout
    real(DP), intent(in) :: hnew(:)
    ! -- local
    integer(I4) :: i, j, k, ierr, nxyz, loc_n, loc_r, i_num, j_num
    integer(I4) :: wtab_sendn, wtab_recvn, sum_sendn, sum_recvn
    integer(I4) :: isend_sta, isend_end, irecv_sta, irecv_end
    integer(I4) :: send_len, recv_len
    integer(I4) :: rank_flag, allp_flag
    integer(I4), allocatable :: send_flag(:), recv_flag(:)
    integer, allocatable :: flag_send(:), req_flag_recv(:)
    integer, allocatable :: level_send(:), req_level_recv(:), psi_send(:), req_psi_recv(:)
    integer, allocatable :: ztop_send(:), req_ztop_recv(:)
    integer, allocatable :: stat_s(:,:), stat_r(:,:)
    real(DP) :: top_psi, below_psi, below_top
    real(DP), allocatable :: send_level(:), send_psi(:), send_ztop(:)
    real(DP), allocatable :: recv_level(:), recv_psi(:), recv_ztop(:)
    logical :: below_flag
    !-------------------------------------------------------------------------------------------
    wtab_sendn = send_wtab_cind(neib_wtab_stotn)
    wtab_recvn = recv_wtab_cind(neib_wtab_rtotn)
    allocate(send_flag(wtab_sendn), recv_flag(wtab_recvn))
    allocate(flag_send(neib_mpi_totn), req_flag_recv(neib_mpi_totn))
    allocate(level_send(neib_mpi_totn), req_level_recv(neib_mpi_totn))
    allocate(psi_send(neib_mpi_totn), req_psi_recv(neib_mpi_totn))
    allocate(ztop_send(neib_mpi_totn), req_ztop_recv(neib_mpi_totn))
    allocate(stat_s(MPI_STATUS_SIZE,neib_mpi_totn), stat_r(MPI_STATUS_SIZE,neib_mpi_totn))
    allocate(send_level(wtab_sendn), send_psi(wtab_sendn), send_ztop(wtab_sendn))
    allocate(recv_level(wtab_recvn), recv_psi(wtab_recvn), recv_ztop(wtab_recvn))
    !$omp parallel
    !$omp do private(i)
    do i = 1, wtab_sendn
      send_flag(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, wtab_recvn
      recv_flag(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, neib_mpi_totn
      flag_send(i) = MPI_REQUEST_NULL ; req_flag_recv(i) = MPI_REQUEST_NULL
      level_send(i) = MPI_REQUEST_NULL ; req_level_recv(i) = MPI_REQUEST_NULL
      psi_send(i) = MPI_REQUEST_NULL ; req_psi_recv(i) = MPI_REQUEST_NULL
      ztop_send(i) = MPI_REQUEST_NULL ; req_ztop_recv(i) = MPI_REQUEST_NULL
    end do
    !$omp end do
    !$omp do private(j)
    do j = 1, neib_mpi_totn
      stat_s(:,j) = 0 ; stat_r(:,j) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, wtab_sendn
      send_level(i) = DZERO ; send_psi(i) = DZERO ; send_ztop(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, wtab_recvn
      recv_level(i) = DZERO ; recv_psi(i) = DZERO ; recv_ztop(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    sum_sendn = 0 ; sum_recvn = 0
    !$omp parallel
    !$omp do private(i) reduction(+:sum_sendn)
    do i = 1, wtab_sendn
      sum_sendn = sum_sendn + send_flag(i)
    end do
    !$omp end do

    !$omp do private(i) reduction(+:sum_recvn)
    do i = 1, wtab_recvn
      sum_recvn = sum_recvn + recv_flag(i)
    end do
    !$omp end do
    !$omp end parallel

    rank_flag = 0 ; allp_flag = 0
    ! -- Set send flag (send_flag)
      call set_send_flag(send_flag)
    wtab_fix_loop: do while (allp_flag /= st_mpi%totn)
      ! -- Set send variable (send_vari)
        call set_send_vari(recv_flag, send_flag, hnew, recv_level, recv_psi, recv_ztop,&
                           send_level, send_psi, send_ztop)

      do i = 1, neib_wtab_stotn
        isend_sta = send_wtab_cind(i-1)+1 ; isend_end = send_wtab_cind(i)
        send_len = isend_end - isend_sta + 1
        if (send_len /= 0) then
          call MPI_ISEND(send_flag(isend_sta), send_len, MPI_INTEGER, neib_wtab_snum(i), 0,&
                         st_mpi%comm, flag_send(i), ierr)
          call MPI_ISEND(send_level(isend_sta), send_len, MPI_REAL8, neib_wtab_snum(i), 1,&
                         st_mpi%comm, level_send(i), ierr)
          call MPI_ISEND(send_psi(isend_sta), send_len, MPI_REAL8, neib_wtab_snum(i), 2,&
                         st_mpi%comm, psi_send(i), ierr)
          call MPI_ISEND(send_ztop(isend_sta), send_len, MPI_REAL8, neib_wtab_snum(i), 3,&
                         st_mpi%comm, ztop_send(i), ierr)
        end if
      end do

      do i = 1, neib_wtab_rtotn
        irecv_sta = recv_wtab_cind(i-1)+1 ; irecv_end = recv_wtab_cind(i)
        recv_len = irecv_end - irecv_sta + 1
        if (irecv_end /= 0) then
          call MPI_IRECV(recv_flag(irecv_sta), recv_len, MPI_INTEGER, neib_wtab_rnum(i), 0,&
                         st_mpi%comm, req_flag_recv(i), ierr)
          call MPI_IRECV(recv_level(irecv_sta), recv_len, MPI_REAL8, neib_wtab_rnum(i), 1,&
                         st_mpi%comm, req_level_recv(i), ierr)
          call MPI_IRECV(recv_psi(irecv_sta), recv_len, MPI_REAL8, neib_wtab_rnum(i), 2,&
                         st_mpi%comm, req_psi_recv(i), ierr)
          call MPI_IRECV(recv_ztop(irecv_sta), recv_len, MPI_REAL8, neib_wtab_rnum(i), 3,&
                         st_mpi%comm, req_ztop_recv(i), ierr)
        end if
      end do

      call MPI_WAITALL(neib_wtab_rtotn, req_flag_recv, stat_r, ierr)
      call MPI_WAITALL(neib_wtab_rtotn, req_level_recv, stat_r, ierr)
      call MPI_WAITALL(neib_wtab_rtotn, req_psi_recv, stat_r, ierr)
      call MPI_WAITALL(neib_wtab_rtotn, req_ztop_recv, stat_r, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (st_mpi%rank == 0) then
          call write_err_stop("Receive water table information.")
        end if
      end if

      call MPI_WAITALL(neib_wtab_stotn, flag_send, stat_s, ierr)
      call MPI_WAITALL(neib_wtab_stotn, level_send, stat_s, ierr)
      call MPI_WAITALL(neib_wtab_stotn, psi_send, stat_s, ierr)
      call MPI_WAITALL(neib_wtab_stotn, ztop_send, stat_s, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (st_mpi%rank == 0) then
          call write_err_stop("Send water table information.")
        end if
      end if

      sum_sendn = 0 ; sum_recvn = 0
      !$omp parallel
      !$omp do private(i) reduction(+:sum_sendn)
      do i = 1, wtab_sendn
        sum_sendn = sum_sendn + send_flag(i)
      end do
      !$omp end do
      !$omp do private(i) reduction(+:sum_recvn)
      do i = 1, wtab_recvn
        sum_recvn = sum_recvn + recv_flag(i)
      end do
      !$omp end do
      !$omp end parallel

      if (sum_sendn == wtab_sendn .and. sum_recvn == wtab_recvn) then
        rank_flag = 1
      end if
      ! -- Sum value for MPI (val)
        call mpisum_val(rank_flag, "water table", allp_flag)
    end do wtab_fix_loop

    !$omp parallel do private(i, k, nxyz, loc_n, loc_r, i_num, j_num, top_psi, below_psi,&
    !$omp&                    below_top, below_flag)
    do i = 1, ncals
      call get_cals_grid(i, i_num, j_num)
      below_flag = .false. ; below_psi = DZERO ; below_top = DZERO
      wtab: do k = st_grid%nz, 1, -1
        nxyz = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
        loc_r = 0 ; loc_r = findloc(recv_wtab_glo(:), value = nxyz, dim = 1)
        if (loc_r /= 0) then
          wtable(i) = recv_level(loc_r)
          if (recv_psi(loc_r) < DZERO) then
            exit wtab
          end if
          below_flag = .true. ; below_psi = recv_psi(loc_r) ; below_top = recv_ztop(loc_r)
        else
          loc_n = 0 ; loc_n = findloc(st_conn%loc2glo_ijk(:ncalc), value = nxyz, dim = 1)
          if (loc_n /= 0) then
            top_psi = hnew(loc_n) - st_geom%cell_top(loc_n)
            if (top_psi < DZERO .and. below_flag) then
              wtable(i) = below_top + below_psi*(st_geom%cell_top(loc_n) - below_top)&
                          /(below_psi - top_psi)
              exit wtab
            end if
            wtable(i) = hnew(loc_n)
            if (top_psi < DZERO) then
              exit wtab
            end if
            below_flag = .true. ; below_psi = top_psi ; below_top = st_geom%cell_top(loc_n)
          end if
        end if
      end do wtab
    end do
    !$omp end parallel do

    deallocate(send_flag, recv_flag, flag_send, req_flag_recv, level_send, req_level_recv)
    deallocate(psi_send, req_psi_recv, ztop_send, req_ztop_recv, stat_s, stat_r)
    deallocate(send_level, send_psi, send_ztop, recv_level, recv_psi, recv_ztop)

  end subroutine calc_mpi_wtable

  subroutine set_send_flag(flag_send)
  !*********************************************************************************************
  ! set_send_flag -- Set send flag
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(out) :: flag_send(:)
    ! -- local
    integer(I4) :: i, j, k, c_num, isend, nxyz, loc_n, loc_r
    integer(I4) :: is_sta, is_end
    integer(I4) :: i_num, j_num, k_num
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num, isend, nxyz, loc_n, loc_r, is_sta, is_end,&
    !$omp&                    i_num, j_num, k_num)
    do i = 1, neib_wtab_stotn
      is_sta = send_wtab_cind(i-1) ; is_end = send_wtab_cind(i)
      do isend = 1, is_end-is_sta
        j = is_sta + isend ; c_num = send_wtab_citem(j)
        call get_calc_grid(c_num, i_num, j_num, k_num)
        sr_flag: do k = st_grid%nz, k_num, -1
          nxyz = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
          loc_r = 0 ; loc_r = findloc(recv_wtab_glo(:), value = nxyz, dim = 1)
          if (loc_r /= 0) then
            flag_send(j) = 0
            exit sr_flag
          end if
          loc_n = 0 ; loc_n = findloc(st_conn%loc2glo_ijk(:ncalc), value = nxyz, dim = 1)
          if (loc_n /= 0) then
            flag_send(j) = 1
            exit sr_flag
          end if
        end do sr_flag
      end do
    end do
    !$omp end parallel do

  end subroutine set_send_flag

  subroutine set_send_vari(r_flag, s_flag, c_head, r_level, r_psi, r_ztop, s_level, s_psi,&
                           s_ztop)
  !*********************************************************************************************
  ! set_send_vari -- Set send variable
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: r_flag(:)
    integer(I4), intent(inout) :: s_flag(:)
    real(DP), intent(in) :: c_head(:)
    real(DP), intent(in) :: r_level(:), r_psi(:), r_ztop(:)
    real(DP), intent(out) :: s_level(:), s_psi(:), s_ztop(:)
    ! -- local
    integer(I4) :: i, j, k, c_num, isend, nxyz, loc_n, loc_r
    integer(I4) :: is_sta, is_end
    integer(I4) :: i_num, j_num, k_num
    real(DP) :: top_psi
    logical :: below_flag
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num, isend, nxyz, loc_n, loc_r, is_sta, is_end,&
    !$omp&                    i_num, j_num, k_num, top_psi, below_flag)
    do i = 1, neib_wtab_stotn
      is_sta = send_wtab_cind(i-1) ; is_end = send_wtab_cind(i)
      do isend = 1, is_end-is_sta
        j = is_sta + isend ; c_num = send_wtab_citem(j)
        call get_calc_grid(c_num, i_num, j_num, k_num)
        s_level(j) = DZERO ; s_psi(j) = DZERO ; s_ztop(j) = DZERO ; below_flag = .false.
        unsat: do k = st_grid%nz, k_num, -1
          nxyz = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
          loc_r = 0 ; loc_r = findloc(recv_wtab_glo(:), value = nxyz, dim = 1)
          if (loc_r /= 0) then
            if (r_flag(loc_r) /= 0) then
              s_flag(j) = 1
              s_level(j) = r_level(loc_r) ; s_psi(j) = r_psi(loc_r)
              s_ztop(j) = r_ztop(loc_r)
              if (s_psi(j) < DZERO) then
                exit unsat
              end if
              below_flag = .true.
            end if
          else
            loc_n = 0 ; loc_n = findloc(st_conn%loc2glo_ijk(:ncalc), value = nxyz, dim = 1)
            if (loc_n /= 0) then
              top_psi = c_head(loc_n) - st_geom%cell_top(loc_n)
              if (top_psi < DZERO .and. below_flag) then
                s_level(j) = s_ztop(j) + s_psi(j)*(st_geom%cell_top(loc_n) - s_ztop(j))&
                             /(s_psi(j) - top_psi)
              else
                s_level(j) = c_head(loc_n)
              end if
              s_psi(j) = top_psi ; s_ztop(j) = st_geom%cell_top(loc_n) ; below_flag = .true.
              if (top_psi < DZERO) then
                exit unsat
              end if
            end if
          end if
        end do unsat
      end do
    end do
    !$omp end parallel do

  end subroutine set_send_vari

end module mpi_write
