module make_linearsystem
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE
  use types_module, only: coef_set
  use initial_module, only: st_ctrl
  use set_cell, only: ncalc, ncals
  use set_condition, only: st_hydr, st_bcnd
  use prep_calculation, only: st_time
  use assign_boundary, only: st_forc
  use allocate_solution, only: nreg_num, st_sol, crs_index
#ifdef MPI_MSG
  use utility_module, only: st_mpi
  use mpi_solve, only: senrec_rvectv
#endif

  implicit none
  private
  public :: allocate_matvec, make_matvec

  ! -- local

  contains

  subroutine allocate_matvec(st_coef)
  !*********************************************************************************************
  ! allocate_matvec -- Allocate for matrix and vector
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    type(coef_set), intent(inout) :: st_coef
    ! -- local
    integer(I4) :: tot_ind
    !-------------------------------------------------------------------------------------------
    tot_ind = crs_index(1)%offind(nreg_num)

    allocate(st_coef%per_srat(ncalc), st_coef%per_relp(nreg_num), st_coef%temp_rhs(nreg_num))
    allocate(st_coef%stod(ncalc), st_coef%cond(nreg_num), st_coef%sead(ncalc))
    allocate(st_coef%dmats(ncalc))
    allocate(st_coef%rivd(ncals), st_coef%lakd(ncals), st_coef%surd(ncals))
    allocate(st_coef%deri_srat(ncalc), st_coef%deri_stor(ncalc))
    allocate(st_coef%deri_dcon(tot_ind), st_coef%rel_hyd(tot_ind), st_coef%deri_lucon(tot_ind))
    allocate(st_coef%deri_con1(tot_ind), st_coef%deri_con2(tot_ind))
    allocate(st_coef%over_riv(st_bcnd%rive_num), st_coef%deri_r(st_bcnd%rive_num))
    allocate(st_coef%deri_ks_riv(st_bcnd%rive_num), st_coef%delh_r(st_bcnd%rive_num))
    allocate(st_coef%per_riv(st_bcnd%rive_num), st_coef%rel_riv(st_bcnd%rive_num))
    allocate(st_coef%tran_riv(st_bcnd%rive_num))
    allocate(st_coef%over_lak(st_bcnd%lake_num), st_coef%deri_l(st_bcnd%lake_num))
    allocate(st_coef%deri_ks_lak(st_bcnd%lake_num), st_coef%delh_l(st_bcnd%lake_num))
    allocate(st_coef%per_lak(st_bcnd%lake_num), st_coef%rel_lak(st_bcnd%lake_num))
    allocate(st_coef%tran_lak(st_bcnd%lake_num))
    allocate(st_coef%over_sur(ncals), st_coef%deri_s(ncals), st_coef%deri_ks_sur(ncals))
    allocate(st_coef%delh_s(ncals), st_coef%tran_sur(ncals))
    allocate(st_coef%deri_sea(st_bcnd%seal_num), st_coef%deri_ks_sea(st_bcnd%seal_num))
    allocate(st_coef%delh_sea(st_bcnd%seal_num), st_coef%per_sea(st_bcnd%seal_num))
    allocate(st_coef%rel_sea(st_bcnd%seal_num), st_coef%tran_sea(st_bcnd%seal_num))

  end subroutine allocate_matvec

  subroutine make_matvec(st_coef)
  !*********************************************************************************************
  ! make_matvec -- Make matrix and vector
  !*********************************************************************************************
    ! -- modules
    use allocate_solution, only: array_var
    use calc_function, only: calc_func
    use make_amg_matrix, only: make_amgmat
    ! -- inout
    type(coef_set), intent(inout) :: st_coef
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, nreg_num
      st_coef%temp_rhs(i) = DZERO
    end do
    !$omp end parallel do

    ! -- Calculate function value (func)
      call calc_func(st_sol%head_new, st_coef%temp_rhs)

    array_var(1)%rhs(:) = -st_coef%temp_rhs(:)

    ! -- Make matrix (matrix)
      call make_matrix(array_var(1)%dmat, array_var(1)%lumat, st_coef)

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
    ! -- Send and Receive real vector value (rvectv)
      call senrec_rvectv(array_var(1)%dmat)
      call senrec_rvectv(array_var(1)%rhs)
    end if
#endif

    if (st_ctrl%precon_type == 1) then
      ! -- Make amg matrix (amgmat)
        call make_amgmat()
    end if

  end subroutine make_matvec

  subroutine make_matrix(diamat, lumat, st_coef)
  !*********************************************************************************************
  ! make_matrix -- Make matrix
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_sim
    use calc_parameter, only: calc_srat_rperm
    use calc_function, only: alp_ss
    ! -- inout
    real(DP), intent(out) :: diamat(:), lumat(:)
    type(coef_set), intent(inout) :: st_coef
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, ncalc
      st_coef%per_srat(i) = DZERO
      st_coef%stod(i) = DZERO ; st_coef%sead(i) = DZERO ; st_coef%dmats(i) = DZERO
    end do
    !$omp end parallel do
    !$omp parallel do private(i)
    do i = 1, nreg_num
      st_coef%per_relp(i) = DZERO
      st_coef%cond(i) = DZERO
    end do
    !$omp end parallel do
    !$omp parallel do private(i)
    do i = 1, ncals
      st_coef%rivd(i) = DZERO ; st_coef%lakd(i) = DZERO ; st_coef%surd(i) = DZERO
    end do
    !$omp end parallel do

    ! -- Calculate saturation and relative permeability (srat_rperm)
      call calc_srat_rperm(ncalc, st_ctrl%newper, st_sol%head_new, st_coef%per_srat, st_coef%per_relp)

    if (st_sim%sim_type >= 0) then
      ! -- Form storage change (stochn)
        call form_stochn(alp_ss, st_coef%stod, st_coef%per_srat, st_coef%deri_srat,&
                         st_coef%deri_stor)
    end if

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
    ! -- Send and Receive real vector value (rvectv)
      call senrec_rvectv(st_coef%per_relp)
    end if
#endif

    ! -- Form connect flow from adjacent cells (connflow)
      call form_connflow(st_coef%cond, lumat, st_coef%per_relp, st_coef%deri_dcon,&
                         st_coef%rel_hyd, st_coef%deri_lucon, st_coef%deri_con1,&
                         st_coef%deri_con2)

    ! -- Set river boundary dmat (rivebound)
      call set_rivebound(st_coef%rivd, st_coef%per_relp, st_coef%over_riv, st_coef%deri_r,&
                         st_coef%deri_ks_riv, st_coef%delh_r, st_coef%per_riv, st_coef%rel_riv,&
                         st_coef%tran_riv)

    ! -- Set lake boundary dmat (lakebound)
      call set_lakebound(st_coef%lakd, st_coef%per_relp, st_coef%over_lak, st_coef%deri_l,&
                         st_coef%deri_ks_lak, st_coef%delh_l, st_coef%per_lak, st_coef%rel_lak,&
                         st_coef%tran_lak)

    ! -- Set surface boundary dmat (surfbound)
      call set_surfbound(st_coef%surd, st_coef%per_relp, st_coef%over_sur, st_coef%deri_s,&
                         st_coef%deri_ks_sur, st_coef%delh_s, st_coef%tran_sur)

    ! -- Set sea boundary dmat (seabound)
      call set_seabound(st_coef%sead, st_coef%per_relp, st_coef%deri_sea, st_coef%deri_ks_sea,&
                        st_coef%delh_sea, st_coef%per_sea, st_coef%rel_sea, st_coef%tran_sea)

    !$omp parallel
    !$omp do private(s)
    do s = 1, ncals
      st_coef%dmats(s) = st_coef%surd(s) + st_coef%rivd(s) + st_coef%lakd(s)
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, ncalc
      diamat(i) = st_coef%dmats(i) + st_coef%stod(i) + st_coef%cond(i) + st_coef%sead(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine make_matrix

  subroutine form_stochn(alp, dmat_sto, per_srat, deri_srat, deri_stor)
  !*********************************************************************************************
  ! form_stochn -- Form storage change
  !*********************************************************************************************
    ! -- modules
    use make_cell, only: st_geom
    ! -- inout
    real(DP), intent(in) :: alp(:)
    real(DP), intent(out) :: dmat_sto(:)
    real(DP), intent(in) :: per_srat(:)
    real(DP), intent(out) :: deri_srat(:), deri_stor(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      deri_srat(i) = DZERO ; deri_stor(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncalc
      deri_srat(i) = (per_srat(i)-st_sol%srat_new(i))*st_ctrl%newper_inv
    end do
    !$omp end do

    if (st_time%form_switch == 1) then
      !$omp do private(i)
      do i = 1, ncalc
        deri_stor(i) = alp(i)*deri_srat(i)*st_sol%head_new(i)
      end do
      !$omp end do
    end if

    !$omp do private(i)
    do i = 1, ncalc
      dmat_sto(i) = -(st_hydr%read_pors(i)*deri_srat(i)+alp(i)*st_sol%srat_new(i)+deri_stor(i))&
                    *st_time%delt_inv*st_geom%cell_vol(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine form_stochn

  subroutine form_connflow(dmat_con, lumat_con, per_relp, deri_dcon, rel_hyd, deri_lucon,&
                           deri_con1, deri_con2)
  !*********************************************************************************************
  ! form_connflow -- Form connect flow from adjacent cells
  !*********************************************************************************************
    ! -- modules
    use calc_parameter, only: calc_hyd_upwind
    ! -- inout
    real(DP), intent(out) :: dmat_con(:)
    real(DP), intent(out) :: lumat_con(:)
    real(DP), intent(in) :: per_relp(:)
    real(DP), intent(out) :: deri_dcon(:), rel_hyd(:), deri_lucon(:), deri_con1(:), deri_con2(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: tot_ind, sta_ind, end_ind, ind
    real(DP) :: delhead, relp1, relp2, per_head1, per_head2, relat, deri_hyd1, deri_hyd2
    !-------------------------------------------------------------------------------------------
    tot_ind = crs_index(1)%offind(nreg_num)
    !$omp parallel
    !$omp do private(i)
    do i = 1, tot_ind
      deri_dcon(i) = DZERO ; rel_hyd(i) = DZERO ; deri_lucon(i) = DZERO
      deri_con1(i) = DZERO ; deri_con2(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i, j, k, sta_ind, end_ind, ind, relat, delhead, relp1, relp2)
    do i = 1, nreg_num
      sta_ind = crs_index(1)%offind(i-1) ; end_ind = crs_index(1)%offind(i)
      do k = 1, end_ind-sta_ind
        ind = sta_ind + k ; j = crs_index(1)%offrow(ind)
        delhead = st_sol%head_new(j) - st_sol%head_new(i)
        relp1 = st_sol%rel_perm(i) ; relp2 = st_sol%rel_perm(j)

        ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
          call calc_hyd_upwind(-delhead, relp1, relp2, relat)

        rel_hyd(ind) = st_hydr%hydf_conn(ind)*relat
        deri_dcon(ind) = -rel_hyd(ind)*st_hydr%abyd_conn(ind)
        deri_lucon(ind) = rel_hyd(ind)*st_hydr%abyd_conn(ind)
      end do
    end do
    !$omp end do

    if (st_time%form_switch == 1) then
      !$omp do private(i, j, k, sta_ind, end_ind, ind, relat, delhead, relp1, relp2,&
      !$omp&           per_head1, per_head2, deri_hyd1, deri_hyd2)
      do i = 1, nreg_num
        sta_ind = crs_index(1)%offind(i-1) ; end_ind = crs_index(1)%offind(i)
        do k = 1, end_ind-sta_ind
          ind = sta_ind + k
          j = crs_index(1)%offrow(ind)

          delhead = st_sol%head_new(j) - st_sol%head_new(i)
          relp1 = per_relp(i) ; relp2 = st_sol%rel_perm(j)

          per_head1 = -delhead + st_ctrl%newper
          per_head2 = -delhead - st_ctrl%newper

          ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
            call calc_hyd_upwind(per_head1, relp1, relp2, relat)
          deri_hyd1 = st_hydr%hydf_conn(ind)*relat

          relp1 = st_sol%rel_perm(i) ; relp2 = per_relp(j)

          ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
            call calc_hyd_upwind(per_head2, relp1, relp2, relat)
          deri_hyd2 = st_hydr%hydf_conn(ind)*relat

          deri_con1(ind) = (deri_hyd1-rel_hyd(ind))*st_ctrl%newper_inv*delhead*&
                           st_hydr%abyd_conn(ind)
          deri_con2(ind) = (deri_hyd2-rel_hyd(ind))*st_ctrl%newper_inv*delhead*&
                           st_hydr%abyd_conn(ind)
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

  end subroutine form_connflow

  subroutine set_rivebound(dmat_riv, per_relp, over_riv, deri_r, deri_ks_riv, delh_r, per_riv,&
                           rel_riv, tran_riv)
  !*********************************************************************************************
  ! set_rivebound -- Set river boundary to dmat
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    real(DP), intent(inout) :: dmat_riv(:)
    real(DP), intent(in) :: per_relp(:)
    real(DP), intent(out) :: over_riv(:), deri_r(:), deri_ks_riv(:), delh_r(:), per_riv(:)
    real(DP), intent(out) :: rel_riv(:), tran_riv(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%rive_num
      over_riv(i) = DZERO ; deri_r(i) = DZERO ; deri_ks_riv(i) = DZERO
      delh_r(i) = DZERO ; per_riv(i) = DZERO ; rel_riv(i) = DZERO
      tran_riv(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i, s)
    do i = 1, st_bcnd%rive_num
      s = st_bcnd%rive2cals(i)
      per_riv(i) = per_relp(s) ; rel_riv(i) = st_sol%rel_perm(s)
      if (st_sol%head_new(s) >= st_forc%rive_bott(i)) then
        delh_r(i) = st_forc%rive_head(i) - st_sol%head_new(s)
        over_riv(i) = DONE
      else if (st_forc%rive_head(i) > st_forc%rive_bott(i)) then
        delh_r(i) = st_forc%rive_head(i) - st_forc%rive_bott(i)
        over_riv(i) = DZERO
      else
        delh_r(i) = DZERO
        over_riv(i) = DONE
      end if
      tran_riv(i) = st_hydr%hydf_surf(s)*st_forc%abyd_rive(i)
      deri_r(i) = -tran_riv(i)*rel_riv(i)
    end do
    !$omp end do

    if (st_time%form_switch == 1) then
      !$omp do private(i)
      do i = 1, st_bcnd%rive_num
        deri_ks_riv(i) = (per_riv(i)-rel_riv(i))*st_ctrl%newper_inv*delh_r(i)*tran_riv(i)&
                         *over_riv(i)
      end do
      !$omp end do
    end if

    !$omp do private(i, s)
    do i = 1, st_bcnd%rive_num
      s = st_bcnd%rive2cals(i)
      dmat_riv(s) = deri_r(i) + deri_ks_riv(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine set_rivebound

  subroutine set_lakebound(dmat_lak, per_relp, over_lak, deri_l, deri_ks_lak, delh_l, per_lak,&
                           rel_lak, tran_lak)
  !*********************************************************************************************
  ! set_lakebound -- Set lake boundary to dmat
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    real(DP), intent(inout) :: dmat_lak(:)
    real(DP), intent(in) :: per_relp(:)
    real(DP), intent(out) :: over_lak(:), deri_l(:), deri_ks_lak(:), delh_l(:), per_lak(:)
    real(DP), intent(out) :: rel_lak(:), tran_lak(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%lake_num
      over_lak(i) = DZERO ; deri_l(i) = DZERO ; deri_ks_lak(i) = DZERO
      delh_l(i) = DZERO ; per_lak(i) = DZERO ; rel_lak(i) =  DZERO
      tran_lak(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i, s)
    do i = 1, st_bcnd%lake_num
      s = st_bcnd%lake2cals(i)
      per_lak(i) = per_relp(s) ; rel_lak(i) = st_sol%rel_perm(s)
      if (st_sol%head_new(s) >= st_forc%lake_bott(i)) then
        delh_l(i) = st_forc%lake_head(i) - st_sol%head_new(s)
        over_lak(i) = DONE
      else if (st_forc%lake_head(i) > st_forc%lake_bott(i)) then
        delh_l(i) = st_forc%lake_head(i) - st_forc%lake_bott(i)
        over_lak(i) = DZERO
      else
        delh_l(i) = DZERO
        over_lak(i) = DONE
      end if
      tran_lak(i) = st_hydr%hydf_surf(s)*st_forc%abyd_lake(i)
      deri_l(i) = -tran_lak(i)*rel_lak(i)
    end do
    !$omp end do

    if (st_time%form_switch == 1) then
      !$omp do private(i)
      do i = 1, st_bcnd%lake_num
        deri_ks_lak(i) = (per_lak(i)-rel_lak(i))*st_ctrl%newper_inv*delh_l(i)*tran_lak(i)&
                         *over_lak(i)
      end do
      !$omp end do
    end if

    !$omp do private(i, s)
    do i = 1, st_bcnd%lake_num
      s = st_bcnd%lake2cals(i)
      dmat_lak(s) = deri_l(i) + deri_ks_lak(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine set_lakebound

  subroutine set_surfbound(dmat_sur, per_relp, over_sur, deri_s, deri_ks_sur, delh_s, tran_sur)
  !*********************************************************************************************
  ! set_surfbound -- Set surface boundary to dmat
  !*********************************************************************************************
    ! -- modules
!    use make_cell, only: surf_elev
    ! -- inout
    real(DP), intent(out) :: dmat_sur(:)
    real(DP), intent(in) :: per_relp(:)
    real(DP), intent(out) :: over_sur(:), deri_s(:), deri_ks_sur(:), delh_s(:), tran_sur(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      over_sur(i) = DZERO ; deri_s(i) = DZERO ; deri_ks_sur(i) = DZERO
      delh_s(i) = DZERO ; tran_sur(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncals
!      if (st_sol%head_new(i) >= surf_elev(i)) then
!        if (st_sol%surf_head(i) > surf_bott(i)) then
!          delh_s(i) = st_sol%surf_head(i) - st_sol%head_new(i)
!        else
!          delh_s(i) = surf_elev(i) - st_sol%head_new(i)
!        end if
!        over_sur(i) = DONE
!      else if (st_sol%surf_head(i) >= surf_elev(i)) then
!        delh_s(i) = st_sol%surf_head(i) - surf_elev(i)
!        over_sur(i) = DZERO
!      else
!!        delh_s(i) = st_sol%surf_head(i) - st_sol%head_new(i)
!        delh_s(i) = DZERO
!        over_sur(i) = DONE
!      end if
!      if (surf_reli(i) == DZERO) then
!        delh_s(i) = DZERO
!      else
        if (st_sol%head_new(i) >= st_hydr%surf_top(i)) then
          if (st_sol%surf_head(i) > st_hydr%surf_top(i)) then
            delh_s(i) = st_sol%surf_head(i) - st_sol%head_new(i)
          else
            delh_s(i) = st_hydr%surf_top(i) - st_sol%head_new(i)
          end if
          over_sur(i) = DONE
        else if (st_sol%head_new(i) >= st_hydr%surf_bott(i)) then
          if (st_sol%surf_head(i) > st_hydr%surf_bott(i)) then
            delh_s(i) = st_sol%surf_head(i) - st_sol%head_new(i)
          else
            delh_s(i) = st_hydr%surf_bott(i) - st_sol%head_new(i)
          end if
          over_sur(i) = DONE
        else if(st_sol%surf_head(i) > st_hydr%surf_bott(i)) then
          delh_s(i) = st_sol%surf_head(i) - st_hydr%surf_bott(i)
          over_sur(i) = DZERO
        else
          delh_s(i) = DZERO
          over_sur(i) = DONE
        end if
        tran_sur(i) = st_hydr%hydf_surf(i)*st_hydr%abyd_surf(i)
        deri_s(i) = -tran_sur(i)*st_sol%rel_perm(i)
!      end if
    end do
    !$omp end do

    if (st_time%form_switch == 1) then
      !$omp do private(i)
      do i = 1, ncals
        deri_ks_sur(i) = (per_relp(i)-st_sol%rel_perm(i))*st_ctrl%newper_inv*delh_s(i)*tran_sur(i)&
                         *over_sur(i)
      end do
      !$omp end do
    end if

    !$omp do private(i)
    do i = 1, ncals
      dmat_sur(i) = deri_s(i) + deri_ks_sur(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine set_surfbound

  subroutine set_seabound(dmat_sea, per_relp, deri_sea, deri_ks_sea, delh_sea, per_sea,&
                          rel_sea, tran_sea)
  !*********************************************************************************************
  ! set_seabound -- Set sea boundary to dmat
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    real(DP), intent(inout) :: dmat_sea(:)
    real(DP), intent(in) :: per_relp(:)
    real(DP), intent(out) :: deri_sea(:), deri_ks_sea(:), delh_sea(:), per_sea(:), rel_sea(:)
    real(DP), intent(out) :: tran_sea(:)
    ! -- local
    integer(I4) :: i, c, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%seal_num
      deri_sea(i) = DZERO ; deri_ks_sea(i) = DZERO ; delh_sea(i) = DZERO
      per_sea(i) = DZERO ; rel_sea(i) = DZERO ; tran_sea(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i, c, s)
    do i = 1, st_bcnd%seal_num
      c = st_bcnd%seal2calc(i) ; s = st_bcnd%seal2seal(i)
      delh_sea(i) = st_forc%read_seal(s) - st_sol%head_new(c)
      per_sea(i) = per_relp(c)*st_hydr%hydf_seal(i)
      rel_sea(i) = st_sol%rel_perm(c)*st_hydr%hydf_seal(i)
      deri_sea(i) = -rel_sea(i)*st_hydr%abyd_seal(i)
    end do
    !$omp end do

    if (st_time%form_switch == 1) then
      !$omp do private(i)
      do i = 1, st_bcnd%seal_num
        deri_ks_sea(i) = (per_sea(i)-rel_sea(i))*st_ctrl%newper_inv*delh_sea(i)&
                         *st_hydr%abyd_seal(i)
      end do
      !$omp end do
    end if

    !$omp do private(i, c)
    do i = 1, st_bcnd%seal_num
      c = st_bcnd%seal2calc(i)
      !$omp atomic
      dmat_sea(c) = dmat_sea(c) + deri_sea(i) + deri_ks_sea(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine set_seabound

end module make_linearsystem
