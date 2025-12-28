module make_amg_matrix
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE
  use allocate_solution, only: crs_index, array_var, pro_var, res_var

  implicit none
  private
  public :: make_amgmat

  ! -- local
  integer(I4) :: amglev, nfine, que_size, ncoase
  real(DP), allocatable :: temp_dmat(:)
  integer(I4), allocatable :: aggr_luflag(:), nonaggr_lu(:), nonaggr_index(:)
  integer(I4), allocatable :: aggr_inresult(:), aggr_belongnum(:)
  integer(I4), allocatable :: que_flag(:), que2glob(:), glob_num(:), lay2_num(:)

  contains

  subroutine make_amgmat()
  !***************************************************************************************
  ! make_amgmat -- Make amg matrix
  !***************************************************************************************
    ! -- modules
    use initial_module, only: nlevel, amg_nlevel
    use set_cell, only: amg_setflag
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    amg_setflag = 0

    amg_level: do amglev = 2, amg_nlevel
      nfine = crs_index(amglev-1)%unknow
      if (nfine == 0) then
        nlevel = amglev - 1
        exit amg_level
      else
        nlevel = amglev
      end if
      ! -- Make graph for next coase level matrix by filtering (grapfilt)
        call make_grapfilt()
      ! -- Make aggregation (aggregation)
        call make_aggregation()

      crs_index(amglev)%unknow = ncoase

      ! -- Make prolongation matrix by soomthed aggregation method (prolomat)
        call make_prolomat()
      ! -- Transpose prolongation matrix to restriction matrix (inter2restr)
        call trans_prolo2restr()
      ! -- Make galerkin coase operater (galerkin)
        call make_galerkin()

      allocate(array_var(amglev)%x(ncoase), array_var(amglev)%rhs(ncoase))
      !$omp parallel workshare
      array_var(amglev)%x(:) = DZERO ; array_var(amglev)%rhs(:) = DZERO
      !$omp end parallel workshare

    end do amg_level

    amg_setflag = 1

  end subroutine make_amgmat

  subroutine make_grapfilt()
  !***************************************************************************************
  ! make_grapfilt -- Make graph for next coase level matrix by filtering
  !***************************************************************************************
    ! -- modules
    use constval_module, only: DTWO, FACE
    use initial_module, only: amg_theta
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, st_off, en_off, nonaggr
    real(DP) :: ddmat, lulumat
    !-------------------------------------------------------------------------------------
    allocate(temp_dmat(nfine))
    allocate(aggr_luflag(nfine), nonaggr_lu(nfine*FACE), nonaggr_index(0:nfine))
    !$omp parallel workshare
    temp_dmat(:) = DZERO
    aggr_luflag(:) = 0 ; nonaggr_lu(:) = 0 ; nonaggr_index(:) = 0
    !$omp end parallel workshare

    nonaggr = 0

    do i = 1, nfine
      temp_dmat(i) = array_var(amglev-1)%dmat(i)
      aggr_luflag(i) = 1
      st_off = crs_index(amglev-1)%offind(i-1) + 1
      en_off = crs_index(amglev-1)%offind(i)
      do k = st_off, en_off
        j = crs_index(amglev-1)%offrow(k)
        ddmat = array_var(amglev-1)%dmat(i)*array_var(amglev-1)%dmat(j)
        lulumat = array_var(amglev-1)%lumat(k)
        if (lulumat**DTWO > abs(ddmat*amg_theta**DTWO) .and. ddmat*lulumat > DZERO) then
          nonaggr = nonaggr + 1
          nonaggr_lu(nonaggr) = k
        else
          temp_dmat(i) = temp_dmat(i) + lulumat
        end if
      end do
      nonaggr_index(i) = nonaggr
      if (nonaggr_index(i) == nonaggr_index(i-1)) then
        aggr_luflag(i) = 0
      end if
    end do

  end subroutine make_grapfilt

  subroutine make_aggregation()
  !***************************************************************************************
  ! make_aggregation -- Make aggregation
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k
    integer(I4), allocatable :: aggr_num(:), aggr_index(:)
    !-------------------------------------------------------------------------------------
    allocate(aggr_num(nfine), aggr_index(0:nfine), que_flag(nfine))
    allocate(que2glob(nfine), glob_num(nfine), lay2_num(nfine))
    !$omp parallel workshare
    aggr_num(:) = 0 ; aggr_index(:) = 0 ; que_flag(:) = 0
    que2glob(:) = 0 ; glob_num(:) = 0 ; lay2_num(:) = 0
    !$omp end parallel workshare

    que_size = 0

    do i = 1, nfine
      if (aggr_luflag(i) == 1) then
        que_size = que_size + 1
        aggr_luflag(i) = que_size
        que_flag(que_size) = 1
        que2glob(que_size) = i
      end if
    end do

    ncoase = 0

    que_loop: do
      if (que_flag(1) == 0) then
        exit que_loop
      end if
      i = que2glob(1)
      ncoase = ncoase + 1
      aggr_num(i) = ncoase
      do k = nonaggr_index(i-1)+1, nonaggr_index(i)
        j = crs_index(amglev-1)%offrow(nonaggr_lu(k))
        aggr_num(j) = ncoase
        if (aggr_luflag(j) > 0) then
          ! -- Remove queue (que)
            call remove_que(j)
        end if
        aggr_luflag(j) = 0
        glob_num(j) = i
      end do
      ! -- Remove queue (que)
        call remove_que(i)
      aggr_luflag(i) = 0
      glob_num(i) = i
      ! -- Remove and Enter queue deep layer (quedeep)
        call remove_enter_quedeep(i)
    end do que_loop

    deallocate(que_flag, que2glob, glob_num, lay2_num)

    ! aggregate check
    aggre_check: do i = 1, nfine
      if (aggr_luflag(i) == 0) then
        do k = nonaggr_index(i-1)+1, nonaggr_index(i)
          j = crs_index(amglev-1)%offrow(nonaggr_lu(k))
          if (aggr_num(j) > 0) then
            aggr_num(i) = aggr_num(j)
            exit aggre_check
          end if
        end do
!        if (aggr_num(i) == 0) then
!          write(*,*) "independent node in aggregate()", aggr_luflag(i), i
!        end if
      end if
    end do aggre_check

    allocate(aggr_inresult(0:ncoase), aggr_belongnum(nfine))
    !$omp parallel workshare
    aggr_inresult(:) = 0 ; aggr_belongnum(:) = 0
    !$omp end parallel workshare

    !$omp parallel do private(i, j)
    do i = 1, nfine
      j = aggr_num(i)
      if (j > 0) then
        aggr_index(j) = aggr_index(j) + 1
      end if
    end do
    !$omp end parallel do

    do i = 1, ncoase
      aggr_index(i) = aggr_index(i) + aggr_index(i-1)
      aggr_inresult(i) = aggr_index(i-1)
    end do

    !$omp parallel do private(i, j, k)
    do i = 1, nfine
      j = aggr_num(i)
      if (j > 0) then
        k = aggr_inresult(j) + 1
        aggr_inresult(j) = k
        aggr_belongnum(k) = i
      end if
    end do
    !$omp end parallel do

!    write(*,*) ncoase, " aggregates are selected from ", aggr_inresult(ncoase)

    deallocate(aggr_num, aggr_index, aggr_luflag)

  end subroutine make_aggregation

  subroutine make_prolomat()
  !***************************************************************************************
  ! make_prolomat -- Make prolongation matrix by soomthed aggregation method
  !***************************************************************************************
    ! -- modules
    use initial_module, only: jac_omega
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, kk, tsapce, fine_row
    real(DP), allocatable :: temp_prov(:,:)
    integer(I4), allocatable :: temp_proag(:,:), temp_prorow(:), glob2aggrn(:)
    !-------------------------------------------------------------------------------------
    tsapce = 30 + 80*(amglev-2)
    if (tsapce > ncoase) then
      tsapce = ncoase
    end if

    allocate(temp_prov(tsapce, nfine))
    allocate(temp_proag(tsapce, nfine), temp_prorow(nfine), glob2aggrn(nfine))
    !$omp parallel workshare
    temp_prov(:,:) = DZERO
    temp_proag(:,:) = 0 ; temp_prorow(:) = 0 ; glob2aggrn(:) = 0
    !$omp end parallel workshare
    !$omp parallel do private(i, k, fine_row)
    do k = 1, ncoase
      do i = aggr_inresult(k-1)+1, aggr_inresult(k)
        fine_row = aggr_belongnum(i)
        glob2aggrn(fine_row) = k
      end do
    end do
    !$omp end parallel do

    do i = 1, nfine
      do kk = nonaggr_index(i-1)+1, nonaggr_index(i)
        j = crs_index(amglev-1)%offrow(nonaggr_lu(kk))
        if (glob2aggrn(j) == 0) then
          cycle
        end if
        prolo_1: do k = 1, temp_prorow(i)
          if (temp_proag(k,i) == glob2aggrn(j)) then
            exit prolo_1
          end if
        end do prolo_1
        if (k > temp_prorow(i)) then
          if (k <= tsapce) then
            temp_prorow(i) = k
            temp_proag(k,i) = glob2aggrn(j)
            temp_prov(k,i) = DZERO
          else
            stop 'Prolongation size should be enlarged. Stop in off-diagonal'
          end if
        end if
      end do

      if (glob2aggrn(i) == 0) then
        cycle
      end if
      prolo_2: do k = 1, temp_prorow(i)
        if (temp_proag(k,i) == glob2aggrn(i)) then
          exit prolo_2
        end if
      end do prolo_2
      if (k > temp_prorow(i)) then
        if (k <= tsapce) then
          temp_prorow(i) = k
          temp_proag(k,i) = glob2aggrn(i)
          temp_prov(k,i) = DZERO
        else
          stop 'Prolongation size should be enlarged. Stop in diagonal'
        end if
      end if
    end do

    do i = 1, nfine
      if (temp_prorow(i) == 0) then
        cycle
      end if
      do kk = nonaggr_index(i-1)+1, nonaggr_index(i)
        j = crs_index(amglev-1)%offrow(nonaggr_lu(kk))
        if (glob2aggrn(j) == 0) then
          cycle
        end if
        prolo_3: do k = 1, temp_prorow(i)
          if (temp_proag(k,i) == glob2aggrn(j)) then
            exit prolo_3
          end if
        end do prolo_3
        if (k > temp_prorow(i)) then
          stop 'Unexpected error in smooth aggregat'
        end if
        temp_prov(k,i) = temp_prov(k,i) - array_var(amglev-1)%lumat(nonaggr_lu(kk))
      end do
    end do
    !$omp parallel do private(i, k)
    do i = 1, nfine
      do k = 1, temp_prorow(i)
        temp_prov(k,i) = temp_prov(k,i)*jac_omega/temp_dmat(i)
      end do
    end do
    !$omp end parallel do

    do i = 1, nfine
      if (glob2aggrn(i) == 0) then
        cycle
      end if
      prolo_4: do k = 1, temp_prorow(i)
        if (temp_proag(k,i) == glob2aggrn(i)) then
          exit prolo_4
        end if
      end do prolo_4
      if (k > temp_prorow(i)) then
        stop 'Unexpected error in smooth aggregat'
      end if
      temp_prov(k,i) = temp_prov(k,i) + DONE - jac_omega
    end do

    deallocate(glob2aggrn)
    deallocate(aggr_inresult, aggr_belongnum)
    deallocate(nonaggr_lu, nonaggr_index, temp_dmat)

    allocate(pro_var(amglev)%pindex(0:nfine))
    !$omp parallel workshare
    pro_var(amglev)%pindex(:) = 0
    !$omp end parallel workshare

    do i = 1, nfine
      pro_var(amglev)%pindex(i) = pro_var(amglev)%pindex(i-1) + temp_prorow(i)
    end do
    k = pro_var(amglev)%pindex(nfine)

    allocate(pro_var(amglev)%poffrow(k))
    allocate(pro_var(amglev)%pval(k))
    !$omp parallel workshare
    pro_var(amglev)%poffrow(:) = 0 ; pro_var(amglev)%pval(:) = DZERO
    !$omp end parallel workshare
    j = 0
    do i = 1, nfine
      do k = 1, temp_prorow(i)
        j = j + 1
        pro_var(amglev)%poffrow(j) = temp_proag(k,i)
        pro_var(amglev)%pval(j) = temp_prov(k,i)
      end do
    end do

    deallocate(temp_prorow, temp_proag, temp_prov)

  end subroutine make_prolomat

  subroutine trans_prolo2restr()
  !***************************************************************************************
  ! trans_prolo2restr -- Transpose prolongation matrix to restriction matrix
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, kk, psta, pend
    !-------------------------------------------------------------------------------------
    allocate(res_var(amglev)%rindex(0:ncoase))
    !$omp parallel workshare
    res_var(amglev)%rindex(:) = 0
    !$omp end parallel workshare
    kk = 0
    !$omp parallel do private(i, j, k, psta, pend) reduction(+:kk)
    do i = 1, nfine
      psta = pro_var(amglev)%pindex(i-1) + 1
      pend = pro_var(amglev)%pindex(i)
      do j = psta, pend
        k = pro_var(amglev)%poffrow(j) + 1
        kk = kk + 1
        if (k <= ncoase) then
          res_var(amglev)%rindex(k) = res_var(amglev)%rindex(k) + 1
        end if
      end do
    end do
    !$omp end parallel do

    allocate(res_var(amglev)%roffrow(kk))
    allocate(res_var(amglev)%rval(kk))

    !$omp parallel workshare
    res_var(amglev)%roffrow(:) = 0 ; res_var(amglev)%rval(:) = DZERO
    !$omp end parallel workshare
    do i = 1, ncoase
      res_var(amglev)%rindex(i) = res_var(amglev)%rindex(i) + res_var(amglev)%rindex(i-1)
    end do

    do i = 1, nfine
      psta = pro_var(amglev)%pindex(i-1) + 1
      pend = pro_var(amglev)%pindex(i)
      do j = psta, pend
        kk = pro_var(amglev)%poffrow(j)
        k = res_var(amglev)%rindex(kk) + 1
        res_var(amglev)%rindex(kk) = k
        res_var(amglev)%roffrow(k) = i
        res_var(amglev)%rval(k) = pro_var(amglev)%pval(j)
      end do
    end do

  end subroutine trans_prolo2restr

  subroutine make_galerkin()
  !***************************************************************************************
  ! make_galerkin -- Make galerkin coase operater
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, kk, st_off, en_off, rsta, rend
    integer(I4) :: galfin_flag, cancel_flag, fin_num, fine_lu, size_flag
    integer(I4) :: node_size, coa_num, ns, new_tconn_num, coase_index, row
    integer(I4), allocatable :: temp_size(:), temp_f2c(:,:), temp_node(:)
    integer(I4), allocatable :: new_row(:), new_col(:), temp_aggr(:)
    real(DP) :: rest_v
    real(DP), allocatable :: temp_val(:,:), temp_lumat(:), new_lumat(:)
    !-------------------------------------------------------------------------------------
    galfin_flag = 0
    allocate(temp_size(nfine), temp_f2c(7,nfine), temp_node(ncoase))
    allocate(crs_index(amglev)%offind(0:ncoase))
    allocate(new_row(nfine*7), new_col(nfine*7))
    allocate(temp_aggr(ncoase))
    allocate(temp_val(7,nfine), new_lumat(nfine*7), temp_lumat(ncoase))
    !$omp parallel workshare
    temp_size(:) = 0 ; temp_f2c(:,:) = 0 ; temp_node(:) = 0
    crs_index(amglev)%offind(:) = 0
    new_row(:) = 0 ; new_col(:) = 0
    temp_aggr(:) = 0
    temp_val(:,:) =  DZERO ; new_lumat(:) = DZERO ; temp_lumat(:) = DZERO
    !$omp end parallel workshare
    new_tconn_num = 0 ; coase_index = 0

    gal_loop: do
      if (galfin_flag == 1) then
        exit gal_loop
      end if
      galfin_flag = 1
      temp_size(:) = 0
      do i = 1, ncoase
        if (temp_aggr(i) > 0) then
          cycle
        end if
        cancel_flag = 0
        rsta = res_var(amglev)%rindex(i-1) + 1
        rend = res_var(amglev)%rindex(i)
        rest_1: do j = rsta, rend
          fin_num = res_var(amglev)%roffrow(j)
          if (fin_num <= nfine .and. temp_size(fin_num) >= 7) then
            cancel_flag = 1
            exit rest_1
          end if
        end do rest_1
        if (cancel_flag == 0) then
          rsta = res_var(amglev)%rindex(i-1) + 1
          rend = res_var(amglev)%rindex(i)
          rest_2: do j = rsta, rend
            fin_num = res_var(amglev)%roffrow(j)
            st_off = crs_index(amglev-1)%offind(fin_num-1) + 1
            en_off = crs_index(amglev-1)%offind(fin_num)
            rest_3: do k = st_off, en_off
              fine_lu = crs_index(amglev-1)%offrow(k)
              if (fine_lu <= nfine .and. temp_size(fine_lu) >= 7) then
                cancel_flag = 1
                exit rest_3
              end if
            end do rest_3
            if (cancel_flag == 1) then
              exit rest_2
            end if
          end do rest_2
        end if

        if (galfin_flag == 0 .or. cancel_flag == 1) then
          galfin_flag = 0
        end if
        if (cancel_flag == 1) then
          cycle
        end if

        temp_aggr(i) = 1

        rsta = res_var(amglev)%rindex(i-1) + 1
        rend = res_var(amglev)%rindex(i)
        do j = rsta, rend
          fin_num = res_var(amglev)%roffrow(j)
          rest_v = res_var(amglev)%rval(j)
          if (fin_num <= nfine) then
            kk = temp_size(fin_num)
            if (kk == 0) then
              size_flag = 1
            else
              size_flag = 0
              if (temp_f2c(kk,fin_num) /= i) then
                size_flag = 1
              end if
            end if

            if (size_flag == 1) then
              kk = kk + 1
              temp_size(fin_num) = kk
              temp_f2c(kk,fin_num) = i
              temp_val(kk,fin_num) = rest_v*array_var(amglev-1)%dmat(fin_num)
            else
              temp_val(kk,fin_num) = temp_val(kk,fin_num) + rest_v*array_var(amglev-1)%dmat(fin_num)
            end if
          end if

          st_off = crs_index(amglev-1)%offind(fin_num-1) + 1
          en_off = crs_index(amglev-1)%offind(fin_num)
          do k = st_off, en_off
            fine_lu = crs_index(amglev-1)%offrow(k)
            if (fine_lu > nfine) then
              cycle
            end if
            kk = temp_size(fine_lu)
            if (kk == 0) then
              size_flag = 1
            else
              size_flag = 0
              if (temp_f2c(kk,fine_lu) /= i) then
                size_flag = 1
              end if
            end if

            if (size_flag == 1) then
              kk = kk + 1
              temp_size(fine_lu) = kk
              temp_f2c(kk,fine_lu) = i
              temp_val(kk,fine_lu) = rest_v*array_var(amglev-1)%lumat(k)
            else
              temp_val(kk,fine_lu) = temp_val(kk,fine_lu) + rest_v*array_var(amglev-1)%lumat(k)
            end if
          end do
        end do
      end do

      do i = 1, ncoase
        node_size = 0
        rsta = res_var(amglev)%rindex(i-1) + 1
        rend = res_var(amglev)%rindex(i)
        do j = rsta, rend
          fin_num = res_var(amglev)%roffrow(j)
          rest_v = res_var(amglev)%rval(j)
          do k = 1, temp_size(fin_num)
            coa_num = temp_f2c(k,fin_num)
            rest_4: do ns = 1, node_size
              if (temp_node(ns) == coa_num) then
                exit rest_4
              end if
            end do rest_4
            if (ns > node_size) then
              node_size = node_size + 1
              temp_node(node_size) = coa_num
            end if
            temp_lumat(coa_num) = temp_lumat(coa_num) + rest_v*temp_val(k,fin_num)
          end do
        end do

        do j = 1, node_size
          if (i /= temp_node(j)) then
            new_tconn_num = new_tconn_num + 1
            k = temp_node(j) + 1
            if (k <= ncoase) then
              crs_index(amglev)%offind(k) = crs_index(amglev)%offind(k) + 1
            end if
          end if
          coase_index = coase_index + 1
          new_row(coase_index) = temp_node(j)
          new_col(coase_index) = i
          new_lumat(coase_index) = temp_lumat(temp_node(j))
          temp_lumat(temp_node(j)) = DZERO
        end do
      end do
    end do gal_loop

    deallocate(temp_size, temp_f2c, temp_val, temp_node, temp_lumat, temp_aggr)

    do i = 1, ncoase
      crs_index(amglev)%offind(i) = crs_index(amglev)%offind(i) + crs_index(amglev)%offind(i-1)
    end do

    allocate(crs_index(amglev)%offrow(new_tconn_num))
    allocate(array_var(amglev)%lumat(new_tconn_num), array_var(amglev)%dmat(ncoase))
    !$omp parallel workshare
    crs_index(amglev)%offrow(:) = 0
    array_var(amglev)%lumat(:) = DZERO ; array_var(amglev)%dmat(:) = DONE
    !$omp end parallel workshare
    crs_index(amglev)%lunum = new_tconn_num
    !$omp parallel do private(j, k, row)
    do j = 1, coase_index
      row = new_row(j)
      if (row /= new_col(j)) then
        k = crs_index(amglev)%offind(row) + 1
        crs_index(amglev)%offind(row) = k
        crs_index(amglev)%offrow(k) = new_col(j)
        array_var(amglev)%lumat(k) = new_lumat(j)
      else
        array_var(amglev)%dmat(row) = new_lumat(j)
      end if
    end do
    !$omp end parallel do

    deallocate(new_row, new_col, new_lumat)

  end subroutine make_galerkin

  subroutine remove_enter_quedeep(n)
  !***************************************************************************************
  ! remove_enter_quedeep -- Remove and Enter queue deep layer
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: n
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: lay2_count, lay2, lay3
    !-------------------------------------------------------------------------------------
    lay2_count = 0

    ! lay2
    do k = nonaggr_index(n-1)+1, nonaggr_index(n)
      j = crs_index(amglev-1)%offrow(nonaggr_lu(k))
      do i = nonaggr_index(j-1)+1, nonaggr_index(j)
        lay2 = crs_index(amglev-1)%offrow(nonaggr_lu(i))
        if (glob_num(lay2) /= n) then
          glob_num(lay2) = n
          lay2_count = lay2_count + 1
          lay2_num(lay2_count) = lay2
          if (aggr_luflag(lay2) > 0) then
            ! -- Remove queue (que)
              call remove_que(lay2)
            aggr_luflag(lay2) = 0
          end if
        end if
      end do
    end do

    ! lay3
    do k = 1, lay2_count
      j = lay2_num(k)
      do i = nonaggr_index(j-1)+1, nonaggr_index(j)
        lay3 = crs_index(amglev-1)%offrow(nonaggr_lu(i))
        if (glob_num(lay3) /= n) then
          glob_num(lay3) = n
          if (aggr_luflag(lay3) > 0) then
            que_flag(aggr_luflag(lay3)) = que_flag(aggr_luflag(lay3)) + 1
            ! -- Enter queue (que)
              call enter_que(lay3)
          end if
        end if
      end do
    end do

  end subroutine remove_enter_quedeep

  subroutine remove_que(n)
  !***************************************************************************************
  ! remove_que -- Remove queue
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: n
    ! -- local
    integer(I4) :: qnum, temq, temg
    !-------------------------------------------------------------------------------------
    qnum = aggr_luflag(n)

    que_flag(qnum) = 0
    que2glob(qnum) = 0

    temq = qnum*2

    rmque_loop: do
      if (temq > que_size) then
        exit rmque_loop
      end if
      if (temq+1 <= que_size .and. que_flag(temq+1) > que_flag(temq)) then
        temq = temq + 1
      end if
      if (que_flag(temq) == 0) then
        exit rmque_loop
      end if

      temg = que2glob(temq)
      if (temg == 0) then
!        write(*,*) qnum, que_flag(qnum), que2glob(qnum), temg, que_flag(temq), que2glob(temq)
      end if

      que_flag(qnum) = que_flag(temq)
      que2glob(qnum) = temg
      aggr_luflag(temg) = qnum

      que_flag(temq) = 0
      que2glob(temq) = 0

      qnum = temq
      temq = qnum*2
    end do rmque_loop

  end subroutine remove_que

  subroutine enter_que(n)
  !***************************************************************************************
  ! enter_que -- Enter queue
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: n
    ! -- local
    integer(I4) :: qnum, temq, temp, temg1, temg2
    !-------------------------------------------------------------------------------------
    qnum = aggr_luflag(n)
    temq = qnum/2

    inque_loop: do
      if (qnum <= 1) then
        exit inque_loop
      end if
      if (que_flag(temq) >= que_flag(qnum)) then
        exit inque_loop
      end if

      temp = que_flag(temq)
      que_flag(temq) = que_flag(qnum)
      que_flag(qnum) = temp

      temg1 = que2glob(qnum)
      temg2 = que2glob(temq)
      que2glob(temq) = temg1
      que2glob(qnum) = temg2

      aggr_luflag(temg2) = qnum
      aggr_luflag(temg1) = temq

      qnum = temq
      temq = qnum/2
    end do inque_loop

  end subroutine enter_que

end module make_amg_matrix
