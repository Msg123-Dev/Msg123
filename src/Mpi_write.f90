module mpi_write
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SNOVAL, DZERO, DONE
  use utility_module, only: write_err_write, write_err_stop
  use initial_module, only: my_rank, st_grid
  use set_cell, only: ncals, ncalc, neib_mpi_totn, loc2glo_ijk, get_calc_grid
  use mpi_initfin, only: my_comm
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
  !***************************************************************************************
  ! write_mpi_2dbin -- Write MPI 2d binary
  !***************************************************************************************
    ! -- module
    use set_cell, only: no_ncals
    ! -- inout
    integer(I4), intent(in) :: out_fh, out_totn
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit, ntime
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, s, ierr
    integer(I4), allocatable :: istat(:)
    real(SP), allocatable :: vari_sp(:)
    !-------------------------------------------------------------------------------------
    allocate(istat(MPI_STATUS_SIZE))
    allocate(vari_sp(ncals+no_ncals+1))
    !$omp parallel
    !$omp workshare
    istat(:) = 0 ; vari_sp(:) = SNOVAL ; vari_sp(1) = ntime
    !$omp end workshare

    !$omp do private(i, s)
    do i = 1, out_totn
      s = calc_num(i) + 1
      vari_sp(s) = real(out_val(i)*out_unit, kind=SP)
    end do
    !$omp end do
    !$omp end parallel

    ierr = 0
    call MPI_FILE_WRITE(out_fh, vari_sp, ncals+no_ncals+1, MPI_REAL4, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_write(out_fh)
      end if
    end if

    deallocate(istat)
    deallocate(vari_sp)

  end subroutine write_mpi_2dbin

  subroutine write_mpi_3dbin(out_fh, out_totn, calc_num, out_unit, out_val, ntime)
  !***************************************************************************************
  ! write_mpi_3dbin -- Write MPI 3d binary
  !***************************************************************************************
    ! -- module
    use set_cell, only: no_ncalc
    ! -- inout
    integer(I4), intent(in) :: out_fh, out_totn
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit, ntime
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, c, ierr
    integer(I4), allocatable :: istat(:)
    real(SP), allocatable :: vari_sp(:)
    !-------------------------------------------------------------------------------------
    allocate(istat(MPI_STATUS_SIZE))
    allocate(vari_sp(ncalc+no_ncalc+1))
    !$omp parallel
    !$omp workshare
    istat(:) = 0 ; vari_sp(:) = SNOVAL ; vari_sp(1) = ntime
    !$omp end workshare

    !$omp do private(i, c)
    do i = 1, out_totn
      c = calc_num(i) + 1
      vari_sp(c) = real(out_val(i)*out_unit, kind=SP)
    end do
    !$omp end do
    !$omp end parallel

    ierr = 0
    call MPI_FILE_WRITE(out_fh, vari_sp, ncalc+no_ncalc+1, MPI_REAL4, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_write(out_fh)
      end if
    end if

    deallocate(istat)
    deallocate(vari_sp)

  end subroutine write_mpi_3dbin

  subroutine write_mpi_rest(out_fh, out_time, out_unit, out_val)
  !***************************************************************************************
  ! write_mpi_rest -- Write mpi restart value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: out_fh
    real(SP), intent(in) :: out_time, out_unit
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: ierr
    integer(I4), allocatable :: istat(:)
    real(DP), allocatable :: out_rest(:)
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------
    allocate(istat(MPI_STATUS_SIZE), out_rest(ncalc+1))
    !$omp parallel workshare
    istat(:) = 0 ; out_rest(:) = DZERO
    out_rest(1) = real(out_time, kind=DP)
    out_rest(2:) = out_val(1:ncalc)*real(out_unit, kind=DP)
    !$omp end parallel workshare

    ierr = 0 ; head_dis = 0
    call MPI_FILE_WRITE_AT_ALL(out_fh, head_dis, out_rest, ncalc+1, MPI_REAL8, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_write(out_fh)
      end if
    end if

    deallocate(istat, out_rest)

  end subroutine write_mpi_rest

  subroutine redu_mpi_mass(num_mass, inout_st)
  !***************************************************************************************
  ! redu_mpi_mass -- Calculate output massbalance for MPI
  !***************************************************************************************
    ! -- module
    use allocate_output, only: st_msout
    ! -- inout
    integer(I4), intent(in) :: num_mass
    type(st_msout), intent(inout) :: inout_st
    ! -- local
    integer(I4) :: ierr
    real(DP), allocatable :: mpi_sto(:), mpi_con(:), mpi_sea(:), mpi_wel(:)
    real(DP), allocatable :: mpi_rec(:), mpi_sur(:), mpi_riv(:), mpi_lak(:), mpi_tot(:)
    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      allocate(mpi_sto(num_mass), mpi_con(num_mass), mpi_sea(num_mass), mpi_wel(num_mass))
      allocate(mpi_rec(num_mass), mpi_sur(num_mass), mpi_riv(num_mass), mpi_lak(num_mass))
      allocate(mpi_tot(num_mass))
      !$omp parallel workshare
      mpi_sto(:) = DZERO ; mpi_con(:) = DZERO ; mpi_sea(:) = DZERO ; mpi_wel(:) = DZERO
      mpi_rec(:) = DZERO ; mpi_sur(:) = DZERO ; mpi_riv(:) = DZERO ; mpi_lak(:) = DZERO
      mpi_tot(:) = DZERO
      !$omp end parallel workshare
    end if

    ierr = 0
    call MPI_REDUCE(inout_st%sto(1), mpi_sto(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum storage for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%con(1), mpi_con(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum connect flow for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%sea(1), mpi_sea(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum sea discharge for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%wel(1), mpi_wel(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum well pumping for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%rec(1), mpi_rec(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum recharge for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%sur(1), mpi_sur(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum surface runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%riv(1), mpi_riv(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum river runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%lak(1), mpi_lak(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum lake runoff for massbalance.")
      end if
    end if

    call MPI_REDUCE(inout_st%tot(1), mpi_tot(1), num_mass, MPI_REAL8, MPI_SUM, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Reduce sum total volume for massbalance.")
      end if
    end if

    if (my_rank == 0) then
      !$omp parallel workshare
      inout_st%sto(:) = mpi_sto(:) ; inout_st%con(:) = mpi_con(:)
      inout_st%sea(:) = mpi_sea(:) ; inout_st%wel(:) = mpi_wel(:)
      inout_st%rec(:) = mpi_rec(:) ; inout_st%sur(:) = mpi_sur(:)
      inout_st%riv(:) = mpi_riv(:) ; inout_st%lak(:) = mpi_lak(:)
      inout_st%tot(:) = mpi_tot(:)
      !$omp end parallel workshare
      deallocate(mpi_sto, mpi_con, mpi_sea, mpi_wel)
      deallocate(mpi_rec, mpi_sur, mpi_riv, mpi_lak)
      deallocate(mpi_tot)
    end if

  end subroutine redu_mpi_mass

  subroutine set_senrec_wtab()
  !***************************************************************************************
  ! set_senrec_wtab -- Set send and receive for water table
  !***************************************************************************************
    ! -- module
    use set_cell, only: neib_num, send_cind, send_citem
    use allocate_solution, only: nreg_num, crs_index, dir_conn
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, jj
    integer(I4) :: temp_item_num, temp_snum, temp_rnum, temp_neib_num
    integer(I4) :: sta_ind, end_ind, ind
    integer(I4), allocatable :: temp_wtab_snum(:), temp_wtab_rnum(:)
    integer(I4), allocatable :: temp_send_cind(:), temp_recv_cind(:)
    integer(I4), allocatable :: temp_send_citem(:), temp_recv_citem(:)
    !-------------------------------------------------------------------------------------
    temp_item_num = crs_index(1)%offind(nreg_num) - crs_index(1)%offind(ncalc)
    allocate(temp_wtab_snum(neib_mpi_totn), temp_wtab_rnum(neib_mpi_totn))
    allocate(temp_send_cind(0:neib_mpi_totn), temp_recv_cind(0:neib_mpi_totn))
    allocate(temp_send_citem(temp_item_num), temp_recv_citem(temp_item_num))
    !$omp parallel workshare
    temp_wtab_snum(:) = -1 ; temp_wtab_rnum(:) = -1
    temp_send_cind(:) = 0 ; temp_recv_cind(:) = 0
    temp_send_citem(:) = 0 ; temp_recv_citem(:) = 0
    !$omp end parallel workshare

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

    !$omp parallel workshare
    neib_wtab_snum(:) = pack(temp_wtab_snum(:), temp_wtab_snum(:) /= -1)
    neib_wtab_rnum(:) = pack(temp_wtab_rnum(:), temp_wtab_rnum(:) /= -1)
    send_wtab_citem(:) = temp_send_citem(1:temp_snum)
    recv_wtab_citem(:) = temp_recv_citem(1:temp_rnum)
    !$omp end parallel workshare

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
  !***************************************************************************************
  ! calc_mpi_wtable -- Calculate water table for MPI
  !***************************************************************************************
    ! -- module
    use initial_module, only: pro_totn
    use set_cell, only: get_cals_grid
    use allocate_output, only: wtable
    use mpi_utility, only: mpisum_val
    ! -- inout
    real(DP), intent(in) :: hnew(:), snew(:)
    ! -- local
    integer(I4) :: i, k, ierr, nxyz, loc_n, loc_r, i_num, j_num
    integer(I4) :: wtab_sendn, wtab_recvn, sum_sendn, sum_recvn
    integer(I4) :: isend_sta, isend_end, irecv_sta, irecv_end
    integer(I4) :: send_len, recv_len
    integer(I4) :: rank_flag, allp_flag
    integer(I4), allocatable :: send_flag(:), recv_flag(:)
    integer(I4), allocatable :: flag_send(:), flag_recv(:)
    integer(I4), allocatable :: head_send(:), head_recv(:), srat_send(:), srat_recv(:)
    integer(I4), allocatable :: stat_s(:,:), stat_r(:,:)
    real(DP), allocatable :: send_head(:), send_srat(:), recv_head(:), recv_srat(:)
    !-------------------------------------------------------------------------------------
    wtab_sendn = send_wtab_cind(neib_wtab_stotn)
    wtab_recvn = recv_wtab_cind(neib_wtab_rtotn)
    allocate(send_flag(wtab_sendn), recv_flag(wtab_recvn))
    allocate(flag_send(neib_mpi_totn), flag_recv(neib_mpi_totn))
    allocate(head_send(neib_mpi_totn), head_recv(neib_mpi_totn))
    allocate(srat_send(neib_mpi_totn), srat_recv(neib_mpi_totn))
    allocate(stat_s(MPI_STATUS_SIZE,neib_mpi_totn), stat_r(MPI_STATUS_SIZE,neib_mpi_totn))
    allocate(send_head(wtab_sendn), send_srat(wtab_sendn))
    allocate(recv_head(wtab_recvn), recv_srat(wtab_recvn))
    !$omp parallel workshare
    send_flag(:) = 0 ; recv_flag(:) = 0 ; flag_send(:) = 0 ; flag_recv(:) = 0
    head_send(:) = 0 ; head_recv(:) = 0 ; srat_send(:) = 0 ; srat_recv(:) = 0
    stat_s(:,:) = 0 ; stat_r(:,:) = 0
    send_head(:) = DZERO ; send_srat(:) = DZERO
    recv_head(:) = DZERO ; recv_srat(:) = DZERO
    sum_sendn = sum(send_flag(:)) ; sum_recvn = sum(recv_flag(:))
    !$omp end parallel workshare

    rank_flag = 0 ; allp_flag = 0
    ! -- Set send flag (send_flag)
      call set_send_flag(send_flag)
    wtab_fix_loop: do while (allp_flag /= pro_totn)
      ! -- Set send variable (send_vari)
        call set_send_vari(recv_flag, send_flag, hnew, snew, recv_head, recv_srat,&
                           send_head, send_srat)

      do i = 1, neib_wtab_stotn
        isend_sta = send_wtab_cind(i-1)+1 ; isend_end = send_wtab_cind(i)
        send_len = isend_end - isend_sta + 1
        if (send_len /= 0) then
          call MPI_ISEND(send_flag(isend_sta), send_len, MPI_INTEGER, neib_wtab_snum(i),&
                         0, my_comm, flag_send(i), ierr)
          call MPI_ISEND(send_head(isend_sta), send_len, MPI_REAL8, neib_wtab_snum(i), 0,&
                         my_comm, head_send(i), ierr)
          call MPI_ISEND(send_srat(isend_sta), send_len, MPI_REAL8, neib_wtab_snum(i), 0,&
                         my_comm, srat_send(i), ierr)
        end if
      end do

      do i = 1, neib_wtab_rtotn
        irecv_sta = recv_wtab_cind(i-1)+1 ; irecv_end = recv_wtab_cind(i)
        recv_len = irecv_end - irecv_sta + 1
        if (irecv_end /= 0) then
          call MPI_IRECV(recv_flag(irecv_sta), recv_len, MPI_INTEGER, neib_wtab_rnum(i),&
                         0, my_comm, flag_recv(i), ierr)
          call MPI_IRECV(recv_head(irecv_sta), recv_len, MPI_REAL8, neib_wtab_rnum(i), 0,&
                         my_comm, head_recv(i), ierr)
          call MPI_IRECV(recv_srat(irecv_sta), recv_len, MPI_REAL8, neib_wtab_rnum(i), 0,&
                         my_comm, srat_recv(i), ierr)
        end if
      end do

      call MPI_WAITALL(neib_wtab_rtotn, flag_recv, stat_r, ierr)
      call MPI_WAITALL(neib_wtab_rtotn, head_recv, stat_r, ierr)
      call MPI_WAITALL(neib_wtab_rtotn, srat_recv, stat_r, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (my_rank == 0) then
          call write_err_stop("Receive water table information.")
        end if
      end if

      call MPI_WAITALL(neib_wtab_stotn, flag_send, stat_s, ierr)
      call MPI_WAITALL(neib_wtab_stotn, head_send, stat_s, ierr)
      call MPI_WAITALL(neib_wtab_stotn, srat_send, stat_s, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (my_rank == 0) then
          call write_err_stop("Send water table information.")
        end if
      end if

      !$omp parallel workshare
      sum_sendn = sum(send_flag(:)) ; sum_recvn = sum(recv_flag(:))
      !$omp end parallel workshare
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
        loc_n = 0 ; loc_n = findloc(loc2glo_ijk(:), value = nxyz, dim = 1)
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

    deallocate(send_flag, recv_flag, flag_send, flag_recv, head_send, head_recv)
    deallocate(srat_send, srat_recv, stat_s, stat_r)
    deallocate(send_head, send_srat, recv_head, recv_srat)

  end subroutine calc_mpi_wtable

  subroutine set_send_flag(flag_send)
  !***************************************************************************************
  ! set_send_flag -- Set send flag
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(out) :: flag_send(:)
    ! -- local
    integer(I4) :: i, j, k, c_num, isend, nxyz, loc_n, loc_r
    integer(I4) :: is_sta, is_end
    integer(I4) :: i_num, j_num, k_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num, isend, nxyz, loc_n, loc_r, is_sta, is_end, i_num, j_num, k_num)
    do i = 1, neib_wtab_stotn
      is_sta = send_wtab_cind(i-1) ; is_end = send_wtab_cind(i)
      do isend = 1, is_end-is_sta
        j = is_sta + isend ; c_num = send_wtab_citem(j)
        call get_calc_grid(c_num, i_num, j_num, k_num)
        sr_flag: do k = st_grid%nz, k_num, -1
          nxyz = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
          loc_n = 0 ; loc_n = findloc(loc2glo_ijk(:), value = nxyz, dim = 1)
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
  !***************************************************************************************
  ! set_send_vari -- Set send variable
  !***************************************************************************************
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
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num, isend, nxyz, loc_n, loc_r, is_sta, is_end, i_num, j_num, k_num)
    do i = 1, neib_wtab_stotn
      is_sta = send_wtab_cind(i-1) ; is_end = send_wtab_cind(i)
      do isend = 1, is_end-is_sta
        j = is_sta + isend ; c_num = send_wtab_citem(j)
        call get_calc_grid(c_num, i_num, j_num, k_num)
        unsat: do k = st_grid%nz, k_num, -1
          nxyz = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
          loc_n = 0 ; loc_n = findloc(loc2glo_ijk(:), value = nxyz, dim = 1)
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

  subroutine set_recv_vari(r_flag, s_flag, c_head, c_srat, r_head, r_srat, s_head, s_srat)
  !***************************************************************************************
  ! set_recv_vari -- Set receive variable
  !***************************************************************************************
    ! -- module
    use constval_module, only: DONE
    ! -- inout
    integer(I4), intent(in) :: r_flag(:)
    integer(I4), intent(inout) :: s_flag(:)
    real(DP), intent(in) :: c_head(:), c_srat(:)
    real(DP), intent(in) :: r_head(:), r_srat(:)
    real(DP), intent(out) :: s_head(:), s_srat(:)
    ! -- local
    integer(I4) :: i, j, k, c_num, irecv, nxyz, loc_n, loc_r
    integer(I4) :: ir_sta, ir_end
    integer(I4) :: i_num, j_num, k_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num, irecv, nxyz, loc_n, loc_r, ir_sta, ir_end, i_num, j_num, k_num)
    do i = 1, neib_wtab_rtotn
      ir_sta = recv_wtab_cind(i-1) ; ir_end = recv_wtab_cind(i)
      do irecv = 1, ir_end-ir_sta
        j = ir_sta + irecv ; c_num = recv_wtab_citem(j)
        call get_calc_grid(c_num, i_num, j_num, k_num)
        unsat: do k = st_grid%nz, k_num, -1
          nxyz = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
          loc_n = 0 ; loc_n = findloc(loc2glo_ijk(:), value = nxyz, dim = 1)
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
            end if
            if (s_srat(j) /= DONE) then
              exit unsat
            end if
          end if
        end do unsat
      end do
    end do
    !$omp end parallel do

  end subroutine set_recv_vari

end module mpi_write
