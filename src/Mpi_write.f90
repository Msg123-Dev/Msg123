module mpi_write
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SNOVAL, DZERO, DONE
  use utility_module, only: st_mpi, write_err_write, write_err_stop
  use initial_module, only: st_grid
  use set_cell, only: get_calc_grid, ncalc, ncals, neib_mpi_totn, st_conn
  use mpi

  implicit none
  private
  public :: write_mpi_2dbin, write_mpi_3dbin, write_mpi_rest
  public :: redu_mpi_mass, set_senrec_wtab, calc_mpi_wtable

  ! -- local
  integer(I4) :: neib_wtab_stotn, neib_wtab_rtotn
  integer(I4), allocatable :: neib_wtab_snum(:), neib_wtab_rnum(:)
  integer(I4), allocatable :: send_wtab_cind(:), recv_wtab_cind(:)
  integer(I4), allocatable :: send_wtab_citem(:), recv_wtab_citem(:)

  contains

  subroutine write_mpi_2dbin(out_fh, out_totn, calc_num, out_unit, out_val, ntime)
  !*********************************************************************************************
  ! write_mpi_2dbin -- Write MPI 2d binary
  !*********************************************************************************************
    ! -- module
    use mpi_set, only: write_2d_ind
    use set_cell, only: no_ncals, seal_snum
    ! -- inout
    integer(I4), intent(in) :: out_fh, out_totn
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit, ntime
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, s, ierr
    integer(I4) :: all_ncals
    integer(I4), allocatable :: istat(:)
    real(SP), allocatable :: vari_sp(:)
    !-------------------------------------------------------------------------------------------
    all_ncals = ncals + no_ncals + seal_snum + 1
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
      vari_sp(s) = real(out_val(i)*out_unit, kind=SP)
    end do
    !$omp end do
    !$omp end parallel
    vari_sp(1) = ntime

    ierr = 0
    call MPI_FILE_WRITE(out_fh, vari_sp, all_ncals, MPI_REAL4, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_write(out_fh)
      end if
    end if

    deallocate(istat)
    deallocate(vari_sp)

  end subroutine write_mpi_2dbin

  subroutine write_mpi_3dbin(out_fh, out_totn, calc_num, out_unit, out_val, ntime)
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
    ! -- local
    integer(I4) :: i, c, ierr
    integer(I4) :: all_ncalc
    integer(I4), allocatable :: istat(:)
    real(SP), allocatable :: vari_sp(:)
    !-------------------------------------------------------------------------------------------
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
      vari_sp(c) = real(out_val(i)*out_unit, kind=SP)
    end do
    !$omp end do
    !$omp end parallel
    vari_sp(1) = ntime

    ierr = 0
    call MPI_FILE_WRITE(out_fh, vari_sp, all_ncalc, MPI_REAL4, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_write(out_fh)
      end if
    end if

    deallocate(istat)
    deallocate(vari_sp)

  end subroutine write_mpi_3dbin

  subroutine write_mpi_rest(out_fh, out_time, out_unit, out_val)
  !*********************************************************************************************
  ! write_mpi_rest -- Write mpi restart value
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: out_fh
    real(SP), intent(in) :: out_time, out_unit
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, ierr
    integer(I4), allocatable :: istat(:)
    real(SP) :: sync_time
    real(DP), allocatable :: out_rest(:)
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------------
    sync_time = out_time
    call MPI_BCAST(sync_time, 1, MPI_REAL4, 0, st_mpi%comm, ierr)
    allocate(istat(MPI_STATUS_SIZE), out_rest(ncalc+1))
    !$omp parallel
    !$omp do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncalc
      out_rest(i+1) = out_val(i)*real(out_unit, kind=DP)
    end do
    !$omp end do
    !$omp end parallel
    out_rest(1) = real(sync_time, kind=DP)

    ierr = 0 ; head_dis = 0
    call MPI_FILE_WRITE_AT_ALL(out_fh, head_dis, out_rest, ncalc+1, MPI_REAL8, istat, ierr)
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
    call MPI_REDUCE(inout_st%sto(1), mpi_sto(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum storage for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%con(1), mpi_con(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum connect flow for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%sea(1), mpi_sea(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum sea discharge for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%wel(1), mpi_wel(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum well pumping for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%rec(1), mpi_rec(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum recharge for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%sur(1), mpi_sur(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum surface runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%riv(1), mpi_riv(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum river runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%lak(1), mpi_lak(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Reduce sum lake runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%tot(1), mpi_tot(1), num_mass, MPI_REAL8, MPI_SUM, 0,&
                    st_mpi%comm, ierr)
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

  subroutine set_senrec_wtab()
  !*********************************************************************************************
  ! set_senrec_wtab -- Set send and receive for water table
  !*********************************************************************************************
    ! -- module
    use set_cell, only: send_cind, send_citem, neib_num
    use allocate_solution, only: nreg_num, dir_conn, crs_index
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, jj
    integer(I4) :: temp_item_num, temp_snum, temp_rnum, temp_neib_num
    integer(I4) :: sta_ind, end_ind, ind
    integer(I4), allocatable :: temp_wtab_snum(:), temp_wtab_rnum(:)
    integer(I4), allocatable :: temp_send_cind(:), temp_recv_cind(:)
    integer(I4), allocatable :: temp_send_citem(:), temp_recv_citem(:)
    !-------------------------------------------------------------------------------------------
    temp_item_num = crs_index(1)%offind(nreg_num) - crs_index(1)%offind(ncalc)
    allocate(temp_wtab_snum(neib_mpi_totn), temp_wtab_rnum(neib_mpi_totn))
    allocate(temp_send_cind(0:neib_mpi_totn), temp_recv_cind(0:neib_mpi_totn))
    allocate(temp_send_citem(temp_item_num), temp_recv_citem(temp_item_num))
    temp_send_cind(0) = 0 ; temp_recv_cind(0) = 0
    !$omp parallel
    !$omp do private(i)
    do i = 1, neib_mpi_totn
      temp_wtab_snum(i) = -1 ; temp_wtab_rnum(i) = -1
      temp_send_cind(i) = 0 ; temp_recv_cind(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, temp_item_num
      temp_send_citem(i) = 0 ; temp_recv_citem(i) = 0
    end do
    !$omp end do
    !$omp end parallel

    temp_snum = 0 ; temp_rnum = 0
    do i = 1, neib_mpi_totn
      do j = send_cind(i-1)+1, send_cind(i)
        jj = send_citem(j)
        sta_ind = crs_index(1)%offind(jj-1) ; end_ind = crs_index(1)%offind(jj)
        do k = 1, end_ind-sta_ind
          ind = sta_ind + k ; temp_neib_num = crs_index(1)%offrow(ind)
          if (dir_conn(ind) == 1 .and. temp_neib_num > ncalc) then
            temp_snum = temp_snum + 1 ; temp_send_citem(temp_snum) = jj
            if (temp_wtab_snum(i) == -1) then
              temp_wtab_snum(i) = neib_num(i)
            end if
          else if (dir_conn(ind) == 6 .and. temp_neib_num > ncalc) then
            temp_rnum = temp_rnum + 1 ; temp_recv_citem(temp_rnum) = temp_neib_num
            if (temp_wtab_rnum(i) == -1) then
              temp_wtab_rnum(i) = neib_num(i)
            end if
          end if
        end do
      end do
      temp_send_cind(i) = temp_snum ; temp_recv_cind(i) = temp_rnum
    end do

    neib_wtab_stotn = count(temp_wtab_snum(:) /= -1)
    neib_wtab_rtotn = count(temp_wtab_rnum(:) /= -1)

    allocate(neib_wtab_snum(neib_wtab_stotn), neib_wtab_rnum(neib_wtab_rtotn))
    allocate(send_wtab_cind(0:neib_wtab_stotn), recv_wtab_cind(0:neib_wtab_rtotn))
    allocate(send_wtab_citem(temp_snum), recv_wtab_citem(temp_rnum))

    neib_wtab_snum(:) = pack(temp_wtab_snum(:), temp_wtab_snum(:) /= -1)
    neib_wtab_rnum(:) = pack(temp_wtab_rnum(:), temp_wtab_rnum(:) /= -1)
    send_wtab_citem(:) = temp_send_citem(1:temp_snum)
    recv_wtab_citem(:) = temp_recv_citem(1:temp_rnum)

    temp_snum = 0 ; temp_rnum = 0
    do i = 1, neib_mpi_totn
      if (temp_send_cind(i) /= temp_send_cind(i-1)) then
        temp_snum = temp_snum + 1
        send_wtab_cind(temp_snum) = temp_send_cind(i)
      end if
      if (temp_recv_cind(i) /= temp_recv_cind(i-1)) then
        temp_rnum = temp_rnum + 1
        recv_wtab_cind(temp_rnum) = temp_recv_cind(i)
      end if
    end do
    send_wtab_cind(0) = 0 ; recv_wtab_cind(0) = 0

  end subroutine set_senrec_wtab

  subroutine calc_mpi_wtable(hnew, snew)
  !*********************************************************************************************
  ! calc_mpi_wtable -- Calculate water table for MPI
  !*********************************************************************************************
    ! -- module
    use mpi_utility, only: mpisum_val
    use set_cell, only: get_cals_grid
    use allocate_output, only: wtable
    ! -- inout
    real(DP), intent(in) :: hnew(:), snew(:)
    ! -- local
    integer(I4) :: i, j, k, ierr, nxyz, loc_n, loc_r, i_num, j_num
    integer(I4) :: wtab_sendn, wtab_recvn, sum_sendn, sum_recvn
    integer(I4) :: isend_sta, isend_end, irecv_sta, irecv_end
    integer(I4) :: send_len, recv_len
    integer(I4) :: rank_flag, allp_flag
    integer(I4), allocatable :: send_flag(:), recv_flag(:)
    integer, allocatable :: flag_send(:), req_flag_recv(:)
    integer, allocatable :: head_send(:), req_head_recv(:), srat_send(:), req_srat_recv(:)
    integer, allocatable :: stat_s(:,:), stat_r(:,:)
    real(DP), allocatable :: send_head(:), send_srat(:), recv_head(:), recv_srat(:)
    !-------------------------------------------------------------------------------------------
    wtab_sendn = send_wtab_cind(neib_wtab_stotn)
    wtab_recvn = recv_wtab_cind(neib_wtab_rtotn)
    allocate(send_flag(wtab_sendn), recv_flag(wtab_recvn))
    allocate(flag_send(neib_mpi_totn), req_flag_recv(neib_mpi_totn))
    allocate(head_send(neib_mpi_totn), req_head_recv(neib_mpi_totn))
    allocate(srat_send(neib_mpi_totn), req_srat_recv(neib_mpi_totn))
    allocate(stat_s(MPI_STATUS_SIZE,neib_mpi_totn), stat_r(MPI_STATUS_SIZE,neib_mpi_totn))
    allocate(send_head(wtab_sendn), send_srat(wtab_sendn))
    allocate(recv_head(wtab_recvn), recv_srat(wtab_recvn))
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
      head_send(i) = MPI_REQUEST_NULL ; req_head_recv(i) = MPI_REQUEST_NULL
      srat_send(i) = MPI_REQUEST_NULL ; req_srat_recv(i) = MPI_REQUEST_NULL
    end do
    !$omp end do
    !$omp do private(j)
    do j = 1, neib_mpi_totn
      stat_s(:,j) = 0 ; stat_r(:,j) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, wtab_sendn
      send_head(i) = DZERO ; send_srat(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, wtab_recvn
      recv_head(i) = DZERO ; recv_srat(i) = DZERO
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
        call set_send_vari(recv_flag, send_flag, hnew, snew, recv_head, recv_srat, send_head,&
                           send_srat)

      do i = 1, neib_wtab_stotn
        isend_sta = send_wtab_cind(i-1)+1 ; isend_end = send_wtab_cind(i)
        send_len = isend_end - isend_sta + 1
        if (send_len /= 0) then
          call MPI_ISEND(send_flag(isend_sta), send_len, MPI_INTEGER, neib_wtab_snum(i), 0,&
                         st_mpi%comm, flag_send(i), ierr)
          call MPI_ISEND(send_head(isend_sta), send_len, MPI_REAL8, neib_wtab_snum(i), 1,&
                         st_mpi%comm, head_send(i), ierr)
          call MPI_ISEND(send_srat(isend_sta), send_len, MPI_REAL8, neib_wtab_snum(i), 2,&
                         st_mpi%comm, srat_send(i), ierr)
        end if
      end do

      do i = 1, neib_wtab_rtotn
        irecv_sta = recv_wtab_cind(i-1)+1 ; irecv_end = recv_wtab_cind(i)
        recv_len = irecv_end - irecv_sta + 1
        if (irecv_end /= 0) then
          call MPI_IRECV(recv_flag(irecv_sta), recv_len, MPI_INTEGER, neib_wtab_rnum(i), 0,&
                         st_mpi%comm, req_flag_recv(i), ierr)
          call MPI_IRECV(recv_head(irecv_sta), recv_len, MPI_REAL8, neib_wtab_rnum(i), 1,&
                         st_mpi%comm, req_head_recv(i), ierr)
          call MPI_IRECV(recv_srat(irecv_sta), recv_len, MPI_REAL8, neib_wtab_rnum(i), 2,&
                         st_mpi%comm, req_srat_recv(i), ierr)
        end if
      end do

      call MPI_WAITALL(neib_wtab_rtotn, req_flag_recv, stat_r, ierr)
      call MPI_WAITALL(neib_wtab_rtotn, req_head_recv, stat_r, ierr)
      call MPI_WAITALL(neib_wtab_rtotn, req_srat_recv, stat_r, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (st_mpi%rank == 0) then
          call write_err_stop("Receive water table information.")
        end if
      end if

      call MPI_WAITALL(neib_wtab_stotn, flag_send, stat_s, ierr)
      call MPI_WAITALL(neib_wtab_stotn, head_send, stat_s, ierr)
      call MPI_WAITALL(neib_wtab_stotn, srat_send, stat_s, ierr)
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

    !$omp parallel do private(i, k, nxyz, loc_n, loc_r, i_num, j_num)
    do i = 1, ncals
      call get_cals_grid(i, i_num, j_num)
      wtab: do k = st_grid%nz, 1, -1
        nxyz = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
        loc_n = 0 ; loc_n = findloc(st_conn%loc2glo_ijk(:), value = nxyz, dim = 1)
        if (loc_n /= 0) then
          loc_r = 0 ; loc_r = findloc(recv_wtab_citem(:), value = loc_n, dim = 1)
          if (loc_r == 0) then
            if (snew(loc_n) /= DONE) then
              wtable(i) = hnew(loc_n)
              exit wtab
            else
              wtable(i) = hnew(loc_n)
            end if
          else
            if (recv_srat(loc_r) /= DONE) then
              wtable(i) = recv_head(loc_r)
              exit wtab
            end if
          end if
        end if
      end do wtab
    end do
    !$omp end parallel do

    deallocate(send_flag, recv_flag, flag_send, req_flag_recv, head_send, req_head_recv)
    deallocate(srat_send, req_srat_recv, stat_s, stat_r)
    deallocate(send_head, send_srat, recv_head, recv_srat)

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
          loc_n = 0 ; loc_n = findloc(st_conn%loc2glo_ijk(:), value = nxyz, dim = 1)
          if (loc_n /= 0) then
            loc_r = 0 ; loc_r = findloc(recv_wtab_citem(:), value = loc_n, dim = 1)
            if (loc_r /= 0) then
              flag_send(j) = 0
            else
              flag_send(j) = 1
            end if
            exit sr_flag
          end if
        end do sr_flag
      end do
    end do
    !$omp end parallel do

  end subroutine set_send_flag

  subroutine set_send_vari(r_flag, s_flag, c_head, c_srat, r_head, r_srat, s_head, s_srat)
  !*********************************************************************************************
  ! set_send_vari -- Set send variable
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: r_flag(:)
    integer(I4), intent(inout) :: s_flag(:)
    real(DP), intent(in) :: c_head(:), c_srat(:)
    real(DP), intent(in) :: r_head(:), r_srat(:)
    real(DP), intent(out) :: s_head(:), s_srat(:)
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
        unsat: do k = st_grid%nz, k_num, -1
          nxyz = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
          loc_n = 0 ; loc_n = findloc(st_conn%loc2glo_ijk(:), value = nxyz, dim = 1)
          if (loc_n /= 0) then
            loc_r = 0 ; loc_r = findloc(recv_wtab_citem(:), value = loc_n, dim = 1)
            if (loc_r == 0) then
              s_head(j) = c_head(loc_n) ; s_srat(j) = c_srat(loc_n)
            else if (r_flag(loc_r) /= 0) then
              s_flag(j) = 1
              if (r_srat(loc_r) /= DONE) then
                s_head(j) = r_head(loc_r) ; s_srat(j) = r_srat(loc_r)
              else
                s_head(j) = c_head(loc_n) ; s_srat(j) = c_srat(loc_n)
              end if
            else
              s_head(j) = c_head(loc_n) ; s_srat(j) = c_srat(loc_n)
            end if
            if (s_srat(j) /= DONE) then
              exit unsat
            end if
          end if
        end do unsat
      end do
    end do
    !$omp end parallel do

  end subroutine set_send_vari

end module mpi_write
