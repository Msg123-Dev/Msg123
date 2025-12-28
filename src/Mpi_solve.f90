module mpi_solve
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO
  use mpi_initfin, only: my_comm
  use utility_module, only: write_err_stop
  use initial_module, only: my_rank, pro_totn
  use set_cell, only: ncalc, neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                      recv_citem
  use allocate_solution, only: nreg_num, crs_index, dir_conn, left_offr
  use mpi_utility, only: mpisum_val
  use mpi

  implicit none
  private
  public :: senrec_ivectv, senrec_rvectv
  public :: precon_mpi_dilu, solve_mpi_ilu
  public :: check_mpimaxerr, bcast_convinfo

  ! -- local

  contains

  subroutine senrec_ivectv(ivector)
  !***************************************************************************************
  ! senrec_ivectv -- Send and Recieve integer vector value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: ivector(:)
    ! -- local
    integer(I4) :: i, j, k, jj, kk, ierr
    integer(I4) :: nsenrev, isend_sta, isend_end, irecv_sta, irecv_end
    integer(I4) :: buflen_send, buflen_recv
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_s(:,:), stat_r(:,:)
    integer(I4), allocatable :: sbufint(:), rbufint(:)
    !-------------------------------------------------------------------------------------
    allocate(requ_send(neib_mpi_totn), requ_recv(neib_mpi_totn))
    allocate(stat_s(MPI_STATUS_SIZE,neib_mpi_totn), stat_r(MPI_STATUS_SIZE,neib_mpi_totn))
    !$omp parallel
    !$omp workshare
    requ_send(:) = 0 ; requ_recv(:) = 0
    stat_s(:,:) = 0 ; stat_r(:,:) = 0
    !$omp end workshare

    nsenrev = max(send_cind(neib_mpi_totn), recv_cind(neib_mpi_totn))

    allocate(sbufint(nsenrev), rbufint(nsenrev))
    !$omp workshare
    sbufint(:) = 0 ; rbufint(:) = 0
    !$omp end workshare

    !$omp do private(i, j, jj)
    do i = 1, neib_mpi_totn
      do j = send_cind(i-1)+1, send_cind(i)
        jj = send_citem(j)
        sbufint(j) = ivector(jj)
      end do
    end do
    !$omp end do

    !$omp do private(i, isend_sta, isend_end, buflen_send)
    do i = 1, neib_mpi_totn
      isend_sta = send_cind(i-1)+1
      isend_end = send_cind(i)
      buflen_send = isend_end - isend_sta + 1
      if (buflen_send /= 0) then
        call MPI_ISEND(sbufint(isend_sta), buflen_send, MPI_INTEGER, neib_num(i), 0,&
                       my_comm, requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, irecv_sta, irecv_end, buflen_recv)
    do i = 1, neib_mpi_totn
      irecv_sta = recv_cind(i-1)+1
      irecv_end = recv_cind(i)
      buflen_recv = irecv_end - irecv_sta + 1
      if (buflen_recv /= 0) then
        call MPI_IRECV(rbufint(irecv_sta), buflen_recv, MPI_INTEGER, neib_num(i), 0,&
                       my_comm, requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(neib_mpi_totn, requ_recv, stat_r, ierr)
    call MPI_WAITALL(neib_mpi_totn, requ_send, stat_s, ierr)

    !$omp do private(i, k, kk)
    do i = 1, neib_mpi_totn
      do k = recv_cind(i-1)+1, recv_cind(i)
        kk = recv_citem(k)
        ivector(kk) = rbufint(k)
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_s, stat_r, sbufint, rbufint)

  end subroutine senrec_ivectv

  subroutine senrec_rvectv(rvector)
  !***************************************************************************************
  ! senrec_rvectv -- Send and Recieve real vector value
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(inout) :: rvector(:)
    ! -- local
    integer(I4) :: i, j, k, jj, kk, ierr
    integer(I4) :: nsenrev, isend_sta, isend_end, irecv_sta, irecv_end
    integer(I4) :: buflen_send, buflen_recv
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_s(:,:), stat_r(:,:)
    real(DP), allocatable :: sbufreal(:), rbufreal(:)
    !-------------------------------------------------------------------------------------
    allocate(requ_send(neib_mpi_totn), requ_recv(neib_mpi_totn))
    allocate(stat_s(MPI_STATUS_SIZE,neib_mpi_totn), stat_r(MPI_STATUS_SIZE,neib_mpi_totn))
    !$omp parallel
    !$omp workshare
    requ_send(:) = 0 ; requ_recv(:) = 0
    stat_s(:,:) = 0 ; stat_r(:,:) = 0
    !$omp end workshare

    nsenrev = max(send_cind(neib_mpi_totn), recv_cind(neib_mpi_totn))

    allocate(sbufreal(nsenrev), rbufreal(nsenrev))
    !$omp workshare
    sbufreal(:) = DZERO ; rbufreal(:) = DZERO
    !$omp end workshare

    !$omp do private(i, j, jj)
    do i = 1, neib_mpi_totn
      do j = send_cind(i-1)+1, send_cind(i)
        jj = send_citem(j)
        sbufreal(j) = rvector(jj)
      end do
    end do
    !$omp end do

    !$omp do private(i, isend_sta, isend_end, buflen_send)
    do i = 1, neib_mpi_totn
      isend_sta = send_cind(i-1)+1
      isend_end = send_cind(i)
      buflen_send = isend_end - isend_sta + 1
      if (buflen_send /= 0) then
        call MPI_ISEND(sbufreal(isend_sta), buflen_send, MPI_REAL8, neib_num(i), 0,&
                       my_comm, requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, irecv_sta, irecv_end, buflen_recv)
    do i = 1, neib_mpi_totn
      irecv_sta = recv_cind(i-1)+1
      irecv_end = recv_cind(i)
      buflen_recv = irecv_end - irecv_sta + 1
      if (buflen_recv /= 0) then
        call MPI_IRECV(rbufreal(irecv_sta), buflen_recv, MPI_REAL8, neib_num(i), 0,&
                       my_comm, requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(neib_mpi_totn, requ_recv, stat_r, ierr)
    call MPI_WAITALL(neib_mpi_totn, requ_send, stat_s, ierr)

    !$omp do private(i, k, kk)
    do i = 1, neib_mpi_totn
      do k = recv_cind(i-1)+1, recv_cind(i)
        kk = recv_citem(k)
        rvector(kk) = rbufreal(k)
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_s, stat_r, sbufreal, rbufreal)

  end subroutine senrec_rvectv

  subroutine senrec_lumatv(lu_mat)
  !***************************************************************************************
  ! senrec_lumatv -- Send and Recieve lu matrix value
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(inout) :: lu_mat(:)
    ! -- local
    integer(I4) :: i, k, kk, ierr
    integer(I4) :: nsenrev, isend_sta, isend_end, irecv_sta, irecv_end
    integer(I4) :: buflen_send, buflen_recv
    integer(I4) :: off_row_num, count_offrow, sta_off, end_off
    integer(I4), allocatable :: requ_send(:), requ_recv(:)
    integer(I4), allocatable :: stat_s(:,:), stat_r(:,:)
    integer(I4), allocatable :: temp_offrow(:)
    real(DP), allocatable :: sbufreal(:), rbufreal(:)
    !-------------------------------------------------------------------------------------
    off_row_num = crs_index(1)%offind(ncalc)
    allocate(temp_offrow(off_row_num))
    !$omp parallel workshare
    temp_offrow(:) = 0
    !$omp end parallel workshare
    count_offrow = 0
    do i = 1, ncalc
      sta_off = crs_index(1)%offind(i-1) + 1
      end_off = crs_index(1)%offind(i)
      do k = sta_off, end_off
        if (crs_index(1)%offrow(k) > ncalc) then
          count_offrow = count_offrow + 1
          temp_offrow(count_offrow) = lu_mat(k)
        end if
      end do
    end do

    allocate(requ_send(neib_mpi_totn), requ_recv(neib_mpi_totn))
    allocate(stat_s(MPI_STATUS_SIZE,neib_mpi_totn), stat_r(MPI_STATUS_SIZE,neib_mpi_totn))
    !$omp parallel workshare
    requ_send(:) = 0 ; requ_recv(:) = 0
    stat_s(:,:) = 0 ; stat_r(:,:) = 0
    !$omp end parallel workshare

    nsenrev = max(send_cind(neib_mpi_totn), recv_cind(neib_mpi_totn))

    allocate(sbufreal(nsenrev), rbufreal(nsenrev))

    !$omp parallel
    !$omp workshare
    sbufreal(:) = DZERO ; rbufreal(:) = DZERO
    !$omp end workshare

    !$omp do private(i, k)
    do i = 1, neib_mpi_totn
      do k = send_cind(i-1)+1, send_cind(i)
        sbufreal(k) = temp_offrow(k)
      end do
    end do
    !$omp end do

    !$omp do private(i, isend_sta, isend_end, buflen_send)
    do i = 1, neib_mpi_totn
      isend_sta = send_cind(i-1)+1
      isend_end = send_cind(i)
      buflen_send = isend_end - isend_sta + 1
      if (buflen_send /= 0) then
        call MPI_ISEND(sbufreal(isend_sta), buflen_send, MPI_REAL8, neib_num(i), 0,&
                       my_comm, requ_send(i), ierr)
      end if
    end do
    !$omp end do

    !$omp do private(i, irecv_sta, irecv_end, buflen_recv)
    do i = 1, neib_mpi_totn
      irecv_sta = recv_cind(i-1)+1
      irecv_end = recv_cind(i)
      buflen_recv = irecv_end - irecv_sta + 1
      if (buflen_recv /= 0) then
        call MPI_IRECV(rbufreal(irecv_sta), buflen_recv, MPI_REAL8, neib_num(i), 0,&
                       my_comm, requ_recv(i), ierr)
      end if
    end do
    !$omp end do

    call MPI_WAITALL(neib_mpi_totn, requ_recv, stat_r, ierr)
    call MPI_WAITALL(neib_mpi_totn, requ_send, stat_s, ierr)

    !$omp do private(i, k, kk)
    do i = 1, neib_mpi_totn
      do k = recv_cind(i-1)+1, recv_cind(i)
        kk = crs_index(1)%offind(recv_citem(k))
        lu_mat(kk) = rbufreal(k)
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(requ_send, requ_recv, stat_s, stat_r, sbufreal, rbufreal, temp_offrow)

  end subroutine senrec_lumatv

  subroutine precon_mpi_dilu(pre_ind, pre_inlu, pre_d)
  !***************************************************************************************
  ! precon_mpi_dilu -- Preconditon mpi incomplete lu diagonal
  !***************************************************************************************
    ! -- module
    use constval_module, only: DONE
    ! -- inout
    real(DP), intent(in) :: pre_ind(:), pre_inlu(:)
    real(DP), intent(inout) :: pre_d(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: off_sta, off_end, off_sta2, off_end2
    integer(I4) :: offr, offr2, off_left, off_row_num
    integer(I4) :: rank_flag, allp_flag
    real(DP) :: d_invk
    integer(I4), allocatable :: fix_flag(:), offr_flag(:)
    real(DP), allocatable :: temp_pred(:)
    !-------------------------------------------------------------------------------------
    off_row_num = crs_index(1)%offind(nreg_num)
    allocate(fix_flag(nreg_num), offr_flag(off_row_num))
    allocate(temp_pred(nreg_num))
    !$omp parallel workshare
    fix_flag(:) = 0 ; offr_flag(:) = 0
    pre_d(:) = pre_ind(:) ; temp_pred(:) = pre_ind(:)
    !$omp end parallel workshare
    rank_flag = 0 ; allp_flag = 0
    prefix_loop: do while (allp_flag /= pro_totn)
      do i = 1, ncalc
        if (fix_flag(i) == 0) then
          off_sta = crs_index(1)%offind(i-1) + 1
          off_end = crs_index(1)%offind(i)
          do k = off_sta, off_end
            offr = crs_index(1)%offrow(k) ; off_left = 0
            off_left = left_offr(k)
            if (off_left /= 0 .and. fix_flag(offr) == 1) then
              offr_flag(k) = 1
              d_invk = DONE/pre_d(offr)
              off_sta2 = crs_index(1)%offind(offr-1) + 1
              off_end2 = crs_index(1)%offind(offr)
              do j = off_sta2, off_end2
                offr2 = crs_index(1)%offrow(j)
                if (3 < dir_conn(j)) then
                  if (offr2 == i) then
                    temp_pred(i) = temp_pred(i) - pre_inlu(k)*d_invk*pre_inlu(j)
                  end if
                end if
              end do
            else if (off_left == 0 .and. offr_flag(k) == 0) then
              offr_flag(k) = 1
            end if
          end do
          if (sum(offr_flag(off_sta:off_end)) == off_end-off_sta+1) then
            fix_flag(i) = 1 ; pre_d(i) = temp_pred(i)
          end if
        end if
      end do
      ! -- Send and Receive integer vector value (ivectv)
        call senrec_ivectv(fix_flag)
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(pre_d)
      if (sum(fix_flag(:)) == nreg_num) then
        rank_flag = 1
      end if
      ! -- Sum value for MPI (val)
        call mpisum_val(rank_flag, "check precondtion", allp_flag)
    end do prefix_loop

    deallocate(fix_flag, offr_flag)

  end subroutine precon_mpi_dilu

  subroutine solve_mpi_ilu(inrhs, indmat, inlumat, outx)
  !***************************************************************************************
  ! solve_mpi_ilu -- Solve mpi ilu factorization
  !***************************************************************************************
    ! -- module
    use allocate_solution, only: right_offr
    ! -- inout
    real(DP), intent(in) :: inrhs(:), indmat(:), inlumat(:)
    real(DP), intent(inout) :: outx(:)
    ! -- local
    integer(I4) :: i, k
    integer(I4) :: off_sta, off_end, offr
    integer(I4) :: off_row_num, off_left, off_right
    integer(I4) :: rank_flag, allp_flag
    integer(I4), allocatable :: fix_flag(:), offr_flag(:)
    real(DP), allocatable :: temp_outx(:)
    !-------------------------------------------------------------------------------------
    off_row_num = crs_index(1)%offind(nreg_num)
    allocate(fix_flag(nreg_num), offr_flag(off_row_num))
    allocate(temp_outx(nreg_num))
    !$omp parallel workshare
    fix_flag(:) = 0 ; offr_flag(:) = 0 ; outx(:) = inrhs(:)
    temp_outx(:) = inrhs(:)
    !$omp end parallel workshare
    rank_flag = 0 ; allp_flag = 0
    forwfix_loop: do while (allp_flag /= pro_totn)
      do i = 1, ncalc
        if (fix_flag(i) == 0) then
          off_sta = crs_index(1)%offind(i-1) + 1
          off_end = crs_index(1)%offind(i)
          do k = off_sta, off_end
            offr = crs_index(1)%offrow(k) ; off_left = 0
            off_left = left_offr(k)
            if (off_left /= 0 .and. fix_flag(offr) == 1) then
              offr_flag(k) = 1 ; temp_outx(i) = temp_outx(i) - inlumat(k)*outx(offr)
            else if (off_left == 0 .and. offr_flag(k) == 0) then
              offr_flag(k) = 1
            end if
          end do
          if (sum(offr_flag(off_sta:off_end)) == off_end-off_sta+1) then
            fix_flag(i) = 1 ; outx(i) = temp_outx(i)/indmat(i)
          end if
        end if
      end do
      ! -- Send and Receive integer vector value (ivectv)
        call senrec_ivectv(fix_flag)
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(outx)
      if (sum(fix_flag(:)) == nreg_num) then
        rank_flag = 1
      end if
      ! -- Sum value for MPI (val)
        call mpisum_val(rank_flag, "forward substitution", allp_flag)
    end do forwfix_loop

    !$omp parallel workshare
    fix_flag(:) = 0 ; offr_flag(:) = 0 ; temp_outx(:) = DZERO
    !$omp end parallel workshare
    rank_flag = 0 ; allp_flag = 0
    backfix_loop: do while (allp_flag /= pro_totn)
      do i = ncalc, 1, -1
        if (fix_flag(i) == 0) then
          off_sta = crs_index(1)%offind(i-1) + 1
          off_end = crs_index(1)%offind(i)
          do k = off_sta, off_end
            offr = crs_index(1)%offrow(k) ; off_right = 0
            off_right = right_offr(k)
            if (off_right /= 0 .and. fix_flag(offr) == 1) then
              offr_flag(k) = 1 ; temp_outx(i) = temp_outx(i) + inlumat(k)*outx(offr)
            else if (off_right == 0 .and. offr_flag(k) == 0) then
              offr_flag(k) = 1
            end if
          end do
          if (sum(offr_flag(off_sta:off_end)) == off_end-off_sta+1) then
            fix_flag(i) = 1 ; outx(i) = outx(i) - temp_outx(i)/indmat(i)
          end if
        end if
      end do
      ! -- Send and Receive integer vector value (ivectv)
        call senrec_ivectv(fix_flag)
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(outx)
      if (sum(fix_flag(:)) == nreg_num) then
        rank_flag = 1
      end if
      ! -- Sum value for MPI (val)
        call mpisum_val(rank_flag, "backward substitution", allp_flag)
    end do backfix_loop

    deallocate(fix_flag, offr_flag)
    deallocate(temp_outx)

  end subroutine solve_mpi_ilu

  subroutine check_mpimaxerr(change, unknow, abs_max_ch, max_un, max_ch)
  !***************************************************************************************
  ! check_mpimaxerr -- Check mpi max error
  !***************************************************************************************
    ! -- module
    use constval_module, only: DNOVAL
    use mpi_utility, only: mpimax_val
    ! -- inout
    real(DP), intent(in) :: change, unknow
    real(DP), intent(out) :: abs_max_ch, max_un, max_ch
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: pnum, wrank
    real(DP) :: mpimax_err
    !-------------------------------------------------------------------------------------
    ! -- MAX value for MPI (val)
      call mpimax_val(abs(change), "absolute change", mpimax_err)
    abs_max_ch = mpimax_err
    ! -- MAX value for MPI (val)
      call mpimax_val(abs(unknow), "absolute unknown", mpimax_err)
    max_un = mpimax_err

    pnum = 0 ; wrank = 0
    if (abs(change) == abs_max_ch) then
      pnum = my_rank
      max_ch = change
    else
      max_ch = DNOVAL
    end if

    ierr = 0
    call MPI_ALLREDUCE(pnum, wrank, 1, MPI_INTEGER, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Allreduce absolute change process.")
      end if
    end if

    call MPI_BCAST(max_ch, 1, MPI_REAL8, wrank, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast absolute change value.")
      end if
    end if

  end subroutine check_mpimaxerr

  subroutine bcast_convinfo(cval, cdmat, crhs, chead, vmax)
  !***************************************************************************************
  ! bcast_convinfo -- Bcast converge information
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(inout) :: cval
    real(DP), intent(inout) :: cdmat, crhs, chead, vmax
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: pnum, wrank
    integer(I4) :: str_len
    !-------------------------------------------------------------------------------------
    pnum = 0 ; wrank = 0
    if (len_trim(adjustl(cval)) /= 0) then
      pnum = my_rank
    end if

    ierr = 0
    call MPI_ALLREDUCE(pnum, wrank, 1, MPI_INTEGER, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Allreduce max absolute change process.")
      end if
    end if

    str_len = len_trim(cval)
    call MPI_BCAST(cval, str_len, MPI_CHARACTER, wrank, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast absolute change grid information.")
      end if
    end if

    call MPI_BCAST(cdmat, 1, MPI_REAL8, wrank, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast absolute change diagonal information.")
      end if
    end if

    call MPI_BCAST(crhs, 1, MPI_REAL8, wrank, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast absolute change right hand information.")
      end if
    end if

    call MPI_BCAST(chead, 1, MPI_REAL8, wrank, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast absolute change maximun value information.")
      end if
    end if

    call MPI_BCAST(vmax, 1, MPI_REAL8, wrank, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        call write_err_stop("Broadcast absolute change unknown value information.")
      end if
    end if

  end subroutine bcast_convinfo

end module mpi_solve
