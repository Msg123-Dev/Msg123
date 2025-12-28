module make_linearsystem
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE
  use initial_module, only: newper, newper_inv
  use set_cell, only: ncalc, ncals
  use set_condition, only: hydf_surf
  use allocate_solution, only: head_new, srat_new, rel_perm, nreg_num, crs_index
  use prep_calculation, only: form_switch
#ifdef MPI_MSG
  use initial_module, only: pro_totn
  use mpi_solve, only: senrec_rvectv
#endif

  implicit none
  private
  public :: make_matvec

  ! -- local
  real(DP), allocatable :: per_srat(:), per_relp(:)

  contains

  subroutine make_matvec()
  !***************************************************************************************
  ! make_matvec -- Make matrix and vector
  !***************************************************************************************
    ! -- modules
    use initial_module, only: precon_type
    use allocate_solution, only: array_var
    use calc_function, only: calc_func
    use make_amg_matrix, only: make_amgmat
    ! -- inout

    ! -- local
    real(DP), allocatable :: temp_rhs(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_rhs(nreg_num))
    !$omp parallel workshare
    temp_rhs(:) = DZERO
    !$omp end parallel workshare

    ! -- Calculate function value (func)
      call calc_func(head_new, temp_rhs)

    !$omp parallel workshare
    array_var(1)%rhs(:) = -temp_rhs(:)
    !$omp end parallel workshare

    ! -- Make matrix (matrix)
      call make_matrix(array_var(1)%dmat, array_var(1)%lumat)

#ifdef MPI_MSG
    if (pro_totn /= 1) then
    ! -- Send and Receive real vector value (rvectv)
      call senrec_rvectv(array_var(1)%dmat)
      call senrec_rvectv(array_var(1)%rhs)
    end if
#endif

    if (precon_type == 1) then
      ! -- Make amg matrix (amgmat)
        call make_amgmat()
    end if

    deallocate(temp_rhs)

  end subroutine make_matvec

  subroutine make_matrix(diamat, lumat)
  !***************************************************************************************
  ! make_matrix -- Make matrix
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_sim
    use calc_parameter, only: calc_srat_rperm
    ! -- inout
    real(DP), intent(out) :: diamat(:), lumat(:)
    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: ss_alp(:)
    real(DP), allocatable :: stod(:), cond(:), sead(:), dmats(:)
    real(DP), allocatable :: rivd(:), lakd(:), surd(:)
    !-------------------------------------------------------------------------------------
    allocate(per_srat(ncalc), per_relp(nreg_num), ss_alp(ncalc))
    allocate(stod(ncalc), cond(nreg_num), sead(ncalc), dmats(ncalc))
    allocate(rivd(ncals), lakd(ncals), surd(ncals))
    !$omp parallel workshare
    per_srat(:) = DZERO ; per_relp(:) = DZERO ; ss_alp(:) = DZERO
    stod(:) = DZERO ; cond(:) = DZERO ; sead(:) = DZERO ; dmats(:) = DZERO
    rivd(:) = DZERO ; lakd(:) = DZERO ; surd(:) = DZERO
    !$omp end parallel workshare

    ! -- Calculate saturation and relative permeability (srat_rperm)
      call calc_srat_rperm(ncalc, newper, head_new, per_srat, per_relp, ss_alp)
    ! -- Calculate saturation and relative permeability (srat_rperm)
      call calc_srat_rperm(ncalc, DZERO, head_new, srat_new, rel_perm, ss_alp)

    if (st_sim%sim_type >= 0) then
      ! -- Form storage change (stochn)
        call form_stochn(ss_alp, stod)
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
    ! -- Send and Receive real vector value (rvectv)
      call senrec_rvectv(head_new)
      call senrec_rvectv(rel_perm)
      call senrec_rvectv(per_relp)
    end if
#endif

    ! -- Form connect flow from adjacent cells (connflow)
      call form_connflow(cond, lumat)

    ! -- Set river boundary dmat (rivebound)
      call set_rivebound(rivd)

    ! -- Set lake boundary dmat (lakebound)
      call set_lakebound(lakd)

    ! -- Set surface boundary dmat (surfbound)
      call set_surfbound(surd)

    ! -- Set sea boundary dmat (seabound)
      call set_seabound(sead)

    !$omp parallel
    !$omp do private(s)
    do s = 1, ncals
      dmats(s) = surd(s) + rivd(s) + lakd(s)
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, ncalc
      diamat(i) = dmats(i) + stod(i) + cond(i) + sead(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(per_srat, per_relp, ss_alp)
    deallocate(stod, cond, sead, dmats)
    deallocate(rivd, lakd, surd)

  end subroutine make_matrix

  subroutine form_stochn(alp, dmat_sto)
  !***************************************************************************************
  ! form_stochn -- Form storage change
  !***************************************************************************************
    ! -- modules
    use make_cell, only: cell_vol
    use assign_calc, only: read_poro
    use prep_calculation, only: delt_inv
    ! -- inout
    real(DP), intent(in) :: alp(:)
    real(DP), intent(out) :: dmat_sto(:)
    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: deri_srat(:), deri_stor(:)
    !-------------------------------------------------------------------------------------
    allocate(deri_srat(ncalc), deri_stor(ncalc))
    !$omp parallel
    !$omp workshare
    deri_srat(:) = DZERO ; deri_stor(:) = DZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncalc
      deri_srat(i) = (per_srat(i)-srat_new(i))*newper_inv
    end do
    !$omp end do

    if (form_switch == 1) then
      !$omp do private(i)
      do i = 1, ncalc
        deri_stor(i) = alp(i)*deri_srat(i)*head_new(i)
      end do
      !$omp end do
    end if

    !$omp do private(i)
    do i = 1, ncalc
      dmat_sto(i) = -(read_poro(i)*deri_srat(i)+alp(i)*srat_new(i)+deri_stor(i))*&
                    delt_inv*cell_vol(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(deri_srat, deri_stor)

  end subroutine form_stochn

  subroutine form_connflow(dmat_con, lumat_con)
  !***************************************************************************************
  ! form_connflow -- Form connect flow from adjacent cells
  !***************************************************************************************
    ! -- modules
    use calc_parameter, only: calc_hyd_upwind
    use allocate_solution, only: abyd_conn, hydf_conn
    ! -- inout
    real(DP), intent(out) :: dmat_con(:)
    real(DP), intent(out) :: lumat_con(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: tot_ind, sta_ind, end_ind, ind
    real(DP) :: delhead, relp1, relp2, per_head1, per_head2, relat, deri_hyd1, deri_hyd2
    real(DP), allocatable :: deri_dcon(:), rel_hyd(:)
    real(DP), allocatable :: deri_lucon(:), deri_con1(:), deri_con2(:)
    !-------------------------------------------------------------------------------------
    tot_ind = crs_index(1)%offind(nreg_num)
    allocate(deri_dcon(tot_ind), rel_hyd(tot_ind))
    allocate(deri_lucon(tot_ind), deri_con1(tot_ind), deri_con2(tot_ind))
    !$omp parallel
    !$omp workshare
    deri_dcon(:) = DZERO ; rel_hyd(:) = DZERO
    deri_lucon(:) = DZERO ; deri_con1(:) = DZERO ; deri_con2(:) = DZERO
    !$omp end workshare

    !$omp do private(i, j, k, sta_ind, end_ind, ind, relat, delhead, relp1, relp2)
    do i = 1, nreg_num
      sta_ind = crs_index(1)%offind(i-1) ; end_ind = crs_index(1)%offind(i)
      do k = 1, end_ind-sta_ind
        ind = sta_ind + k ; j = crs_index(1)%offrow(ind)
        delhead = head_new(j) - head_new(i)
        relp1 = rel_perm(i) ; relp2 = rel_perm(j)

        ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
          call calc_hyd_upwind(-delhead, relp1, relp2, relat)

        rel_hyd(ind) = hydf_conn(ind)*relat
        deri_dcon(ind) = -rel_hyd(ind)*abyd_conn(ind)
        deri_lucon(ind) = rel_hyd(ind)*abyd_conn(ind)
      end do
    end do
    !$omp end do

    if (form_switch == 1) then
      !$omp do private(i, j, k, sta_ind, end_ind, ind, relat, delhead, relp1, relp2, per_head1, per_head2, deri_hyd1, deri_hyd2)
      do i = 1, nreg_num
        sta_ind = crs_index(1)%offind(i-1) ; end_ind = crs_index(1)%offind(i)
        do k = 1, end_ind-sta_ind
          ind = sta_ind + k
          j = crs_index(1)%offrow(ind)

          delhead = head_new(j) - head_new(i)
          relp1 = per_relp(i) ; relp2 = rel_perm(j)

          per_head1 = -delhead + newper
          per_head2 = -delhead - newper

          ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
            call calc_hyd_upwind(per_head1, relp1, relp2, relat)
          deri_hyd1 = hydf_conn(ind)*relat

          relp1 = rel_perm(i) ; relp2 = per_relp(j)

          ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
            call calc_hyd_upwind(per_head2, relp1, relp2, relat)
          deri_hyd2 = hydf_conn(ind)*relat

          deri_con1(ind) = (deri_hyd1-rel_hyd(ind))*newper_inv*delhead*abyd_conn(ind)
          deri_con2(ind) = (deri_hyd2-rel_hyd(ind))*newper_inv*delhead*abyd_conn(ind)
        end do
      end do
      !$omp end do
    end if

    !$omp do private(i, k, sta_ind, end_ind, ind)
    do i = 1, nreg_num
      sta_ind = crs_index(1)%offind(i-1) ; end_ind = crs_index(1)%offind(i)
      do k = 1, end_ind-sta_ind
        ind = sta_ind + k
        lumat_con(ind) = deri_lucon(ind) + deri_con2(ind)
        dmat_con(i) = dmat_con(i) + deri_dcon(ind) + deri_con1(ind)
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(deri_dcon, rel_hyd)
    deallocate(deri_lucon, deri_con1, deri_con2)

  end subroutine form_connflow

  subroutine set_rivebound(dmat_riv)
  !***************************************************************************************
  ! set_rivebound -- Set river boundary to dmat
  !***************************************************************************************
    ! -- modules
    use calc_boundary, only: rive2cals, rive_head, rive_bott
    use set_boundary, only: rive_num, abyd_rive
    ! -- inout
    real(DP), intent(out) :: dmat_riv(:)
    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: over_riv(:), deri_r(:), deri_ks(:), delh_r(:)
    real(DP), allocatable :: per_riv(:), rel_riv(:), tran_riv(:)
    !-------------------------------------------------------------------------------------
    allocate(over_riv(rive_num), deri_r(rive_num), deri_ks(rive_num), delh_r(rive_num))
    allocate(per_riv(rive_num), rel_riv(rive_num), tran_riv(rive_num))
    !$omp parallel
    !$omp workshare
    over_riv(:) = DZERO ; deri_r(:) = DZERO ; deri_ks(:) = DZERO ; delh_r(:) = DZERO
    per_riv(:) = DZERO ; rel_riv(:) = DZERO ; tran_riv(:) = DZERO
    !$omp end workshare

    !$omp do private(i, s)
    do i = 1, rive_num
      s = rive2cals(i)
      per_riv(i) = per_relp(s) ; rel_riv(i) = rel_perm(s)
      if (head_new(s) >= rive_bott(i)) then
        delh_r(i) = rive_head(i) - head_new(s)
        over_riv(i) = DONE
      else if (rive_head(i) > rive_bott(i)) then
        delh_r(i) = rive_head(i) - rive_bott(i)
        over_riv(i) = DZERO
      else
!        delh_r(i) = rive_head(i) - head_new(s)
        delh_r(i) = DZERO
        over_riv(i) = DONE
      end if
      tran_riv(i) = hydf_surf(s)*abyd_rive(i)
      deri_r(i) = -tran_riv(i)*rel_riv(i)
    end do
    !$omp end do

    if (form_switch == 1) then
      !$omp do private(i)
      do i = 1, rive_num
        deri_ks(i) = (per_riv(i)-rel_riv(i))*newper_inv*delh_r(i)*tran_riv(i)*over_riv(i)
      end do
      !$omp end do
    end if

    !$omp do private(i, s)
    do i = 1, rive_num
      s = rive2cals(i)
      dmat_riv(s) = deri_r(i) + deri_ks(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(over_riv, deri_r, deri_ks, delh_r, per_riv, rel_riv, tran_riv)

  end subroutine set_rivebound

  subroutine set_lakebound(dmat_lak)
  !***************************************************************************************
  ! set_lakebound -- Set lake boundary to dmat
  !***************************************************************************************
    ! -- modules
    use calc_boundary, only: lake2cals, lake_head, lake_bott
    use set_boundary, only: lake_num, abyd_lake
    ! -- inout
    real(DP), intent(out) :: dmat_lak(:)
    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: over_lak(:), deri_l(:), deri_ks(:), delh_l(:)
    real(DP), allocatable :: per_lak(:), rel_lak(:), tran_lak(:)
    !-------------------------------------------------------------------------------------
    allocate(over_lak(lake_num), deri_l(lake_num), deri_ks(lake_num), delh_l(lake_num))
    allocate(per_lak(lake_num), rel_lak(lake_num), tran_lak(lake_num))
    !$omp parallel
    !$omp workshare
    over_lak(:) = DZERO ; deri_l(:) = DZERO ; deri_ks(:) = DZERO ; delh_l(:) = DZERO
    per_lak(:) = DZERO ; rel_lak(:) =  DZERO ; tran_lak(:) = DZERO
    !$omp end workshare
    !$omp do private(i, s)
    do i = 1, lake_num
      s = lake2cals(i)
      per_lak(i) = per_relp(s) ; rel_lak(i) = rel_perm(s)
      if (head_new(s) >= lake_bott(i)) then
        delh_l(i) = lake_head(i) - head_new(s)
        over_lak(i) = DONE
      else if (lake_head(i) > lake_bott(i)) then
        delh_l(i) = lake_head(i) - lake_bott(i)
        over_lak(i) = DZERO
      else
!        delh_l(i) = lake_head(i) - head_new(s)
        delh_l(i) = DZERO
        over_lak(i) = DONE
      end if
      tran_lak(i) = hydf_surf(s)*abyd_lake(i)
      deri_l(i) = -tran_lak(i)*rel_lak(i)
    end do
    !$omp end do

    if (form_switch == 1) then
      !$omp do private(i)
      do i = 1, lake_num
        deri_ks(i) = (per_lak(i)-rel_lak(i))*newper_inv*delh_l(i)*tran_lak(i)*over_lak(i)
      end do
      !$omp end do
    end if

    !$omp do private(i, s)
    do i = 1, lake_num
      s = lake2cals(i)
      dmat_lak(s) = deri_l(i) + deri_ks(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(over_lak, deri_l, deri_ks, delh_l, per_lak, rel_lak, tran_lak)

  end subroutine set_lakebound

  subroutine set_surfbound(dmat_sur)
  !***************************************************************************************
  ! set_surfbound -- Set surface boundary to dmat
  !***************************************************************************************
    ! -- modules
!    use make_cell, only: surf_elev
    use assign_calc, only: surf_bott
    use prep_calculation, only: surf_top
    use set_condition, only: abyd_surf
    use allocate_solution, only: surf_head
    ! -- inout
    real(DP), intent(out) :: dmat_sur(:)
    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: over_sur(:), deri_s(:), deri_ks(:), delh_s(:), tran_sur(:)
    !-------------------------------------------------------------------------------------
    allocate(over_sur(ncals), deri_s(ncals), deri_ks(ncals), delh_s(ncals), tran_sur(ncals))
    !$omp parallel
    !$omp workshare
    over_sur(:) = DZERO ; deri_s(:) = DZERO ; deri_ks(:) = DZERO ; delh_s(:) = DZERO
    tran_sur(:) = DZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
!      if (head_new(i) >= surf_elev(i)) then
!        if (surf_head(i) > surf_bott(i)) then
!          delh_s(i) = surf_head(i) - head_new(i)
!        else
!          delh_s(i) = surf_elev(i) - head_new(i)
!        end if
!        over_sur(i) = DONE
!      else if (surf_head(i) >= surf_elev(i)) then
!        delh_s(i) = surf_head(i) - surf_elev(i)
!        over_sur(i) = DZERO
!      else
!!        delh_s(i) = surf_head(i) - head_new(i)
!        delh_s(i) = DZERO
!        over_sur(i) = DONE
!      end if
!      if (surf_reli(i) == DZERO) then
!        delh_s(i) = DZERO
!      else
        if (head_new(i) >= surf_top(i)) then
          if (surf_head(i) > surf_top(i)) then
            delh_s(i) = surf_head(i) - head_new(i)
          else
            delh_s(i) = surf_top(i) - head_new(i)
          end if
          over_sur(i) = DONE
        else if (head_new(i) >= surf_bott(i)) then
          if (surf_head(i) > surf_bott(i)) then
            delh_s(i) = surf_head(i) - head_new(i)
          else
            delh_s(i) = surf_bott(i) - head_new(i)
          end if
          over_sur(i) = DONE
        else if(surf_head(i) > surf_bott(i)) then
          delh_s(i) = surf_head(i) - surf_bott(i)
          over_sur(i) = DZERO
        else
          delh_s(i) = DZERO
          over_sur(i) = DONE
        end if
        tran_sur(i) = hydf_surf(i)*abyd_surf(i)
        deri_s(i) = -tran_sur(i)*rel_perm(i)
!      end if
    end do
    !$omp end do

    if (form_switch == 1) then
      !$omp do private(i)
      do i = 1, ncals
        deri_ks(i) = (per_relp(i)-rel_perm(i))*newper_inv*delh_s(i)*tran_sur(i)*over_sur(i)
      end do
      !$omp end do
    end if

    !$omp do private(i)
    do i = 1, ncals
      dmat_sur(i) = deri_s(i) + deri_ks(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(over_sur, deri_s, deri_ks, delh_s, tran_sur)

  end subroutine set_surfbound

  subroutine set_seabound(dmat_sea)
  !***************************************************************************************
  ! set_seabound -- Set sea boundary to dmat
  !***************************************************************************************
    ! -- modules
    use set_condition, only: nseal
    use assign_boundary, only: read_seal
    use allocate_solution, only: seal2calc, seal2seal, hydf_seal, abyd_seal
    ! -- inout
    real(DP), intent(inout) :: dmat_sea(:)
    ! -- local
    integer(I4) :: i, c, s
    real(DP), allocatable :: deri_sea(:), deri_ks(:), delh_sea(:)
    real(DP), allocatable :: per_sea(:), rel_sea(:), tran_sea(:)
    !-------------------------------------------------------------------------------------
    allocate(deri_sea(nseal), deri_ks(nseal), delh_sea(nseal))
    allocate(per_sea(nseal), rel_sea(nseal), tran_sea(nseal))
    !$omp parallel
    !$omp workshare
    deri_sea(:) = DZERO ; deri_ks(:) = DZERO ; delh_sea(:) = DZERO
    per_sea(:) = DZERO ; rel_sea(:) = DZERO ; tran_sea(:) = DZERO
    !$omp end workshare

    !$omp do private(i, c, s)
    do i = 1, nseal
      c = seal2calc(i) ; s = seal2seal(i)
      delh_sea(i) = read_seal(s) - head_new(c)
      per_sea(i) = per_relp(c)*hydf_seal(i) ; rel_sea(i) = rel_perm(c)*hydf_seal(i)
      deri_sea(i) = -rel_sea(i)*abyd_seal(i)
    end do
    !$omp end do

    if (form_switch == 1) then
      !$omp do private(i)
      do i = 1, nseal
        deri_ks(i) = (per_sea(i)-rel_sea(i))*newper_inv*delh_sea(i)*abyd_seal(i)
      end do
      !$omp end do
    end if
    !$omp end parallel

    do i = 1, nseal
      c = seal2calc(i)
      dmat_sea(c) = dmat_sea(c) + deri_sea(i) + deri_ks(i)
    end do

    deallocate(deri_sea, deri_ks, delh_sea, per_sea, rel_sea, tran_sea)

  end subroutine set_seabound

end module make_linearsystem
