module calc_function
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE
  use initial_module, only: st_sim
  use set_cell, only: ncalc, ncals
  use set_condition, only: nseal, hydf_surf
  use set_boundary, only: rive_num, lake_num
  use calc_parameter, only: calc_srat_rperm
  use allocate_solution, only: srat_new, rel_perm, crs_index
#ifdef MPI_MSG
  use initial_module, only: pro_totn
  use mpi_solve, only: senrec_rvectv
#endif

  implicit none
  private
  public :: allocate_calfun, calc_func, calc_mass, calc_vecjacf
  public :: func_rechterm, func_wellterm, func_surfterm, func_riveterm
  public :: func_laketerm, func_sealterm, alp_ss

  ! -- local
  real(DP), allocatable, private :: stof(:), conf(:), welf(:), seaf(:)
  real(DP), allocatable, private :: funcvs(:)
  real(DP), allocatable          :: alp_ss(:)
  real(DP), allocatable, private :: recf(:), surf(:), rivf(:), lakf(:)
  real(DP), allocatable, private :: ch_stor1(:), ch_stor2(:), ch_stor(:)
  real(DP), allocatable, private :: conn_flow(:)
  real(DP), allocatable, private :: delh_s(:), elev_rati(:)
  real(DP), allocatable, private :: delh_r(:), delh_l(:), seal_flow(:)

  contains

  subroutine allocate_calfun()
  !*********************************************************************************************
  ! allocate_calfun -- Allocate for calculate function value
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------------
    allocate(stof(ncalc), conf(ncalc), welf(ncalc), seaf(ncalc))
    allocate(funcvs(ncalc), alp_ss(ncalc))
    allocate(recf(ncals), surf(ncals), rivf(ncals), lakf(ncals))
    allocate(ch_stor1(ncalc), ch_stor2(ncalc), ch_stor(ncalc))
    allocate(conn_flow(ncalc))
    allocate(delh_s(ncals), elev_rati(ncals))
    allocate(delh_r(rive_num), delh_l(lake_num), seal_flow(nseal))

  end subroutine allocate_calfun

  subroutine calc_func(infx, funcv)
  !*********************************************************************************************
  ! calc_func -- Calculate function value
  !*********************************************************************************************
    ! -- modules
    use allocate_solution, only: surf_head
    ! -- inout
    real(DP), intent(inout) :: infx(:)
    real(DP), intent(out) :: funcv(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      stof(i) = DZERO ; conf(i) = DZERO ; welf(i) = DZERO ; seaf(i) = DZERO
      funcvs(i) = DZERO ; alp_ss(i) = DZERO ; funcv(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncals
      recf(i) = DZERO ; surf(i) = DZERO ; rivf(i) = DZERO ; lakf(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    ! -- Calculate saturation and relative permeability (srat_rperm)
      call calc_srat_rperm(ncalc, DZERO, infx, srat_new, rel_perm, alp_ss)

    if (st_sim%sim_type >= 0) then
      ! -- Function storage change (stochn)
        call func_stochn(infx, alp_ss, alp_ss, stof)
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(infx)
        call senrec_rvectv(rel_perm)
    end if
#endif

    ! -- Function connect flow from adjacent cells (connflow)
      call func_connflow(infx, conf)

    ! -- Function recharge term (rechterm)
      call func_rechterm(recf)

    ! -- Function well term (wellterm)
      call func_wellterm(welf)

    ! -- Function surface term (surfterm)
      call func_surfterm(infx, surf_head, surf)

    ! -- Function river term (riveterm)
      call func_riveterm(infx, rivf)

    ! -- Function lake term (laketerm)
      call func_laketerm(infx, lakf)

    ! -- Function sea level term (sealterm)
      call func_sealterm(infx, seaf)

    !$omp parallel do private(s)
    do s = 1, ncals
      funcvs(s) = recf(s) + surf(s) + rivf(s) + lakf(s)
    end do
    !$omp end parallel do

    !$omp parallel do private(i)
    do i = 1, ncalc
      funcv(i) = funcvs(i) + stof(i) + conf(i) + welf(i) + seaf(i)
    end do
    !$omp end parallel do

  end subroutine calc_func

  subroutine calc_mass(sfla, inmxn, inmxo, stom, conm, seam, welm, recm, surm, rivm, lakm)
  !*********************************************************************************************
  ! calc_mass -- Calculate function value for massbalance
  !*********************************************************************************************
    ! -- modules
    use prep_calculation, only: delt
    use allocate_solution, only: surf_old
    ! -- inout
    integer(I4), intent(in) :: sfla
    real(DP), intent(in) :: inmxo(:)
    real(DP), intent(inout) :: inmxn(:)
    real(DP), intent(inout) :: stom(:), conm(:), seam(:), welm(:)
    real(DP), intent(inout) :: recm(:), surm(:), rivm(:), lakm(:)
    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: alp_ss_new(:), alp_ss_old(:)
    !-------------------------------------------------------------------------------------------
    allocate(alp_ss_new(ncalc), alp_ss_old(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      alp_ss_new(i) = DZERO ; alp_ss_old(i) = DZERO
    end do
    !$omp end parallel do

    ! -- Calculate saturation and relative permeability (srat_rperm)
      call calc_srat_rperm(ncalc, DZERO, inmxo, srat_new, rel_perm, alp_ss_old)
    ! -- Calculate saturation and relative permeability (srat_rperm)
      call calc_srat_rperm(ncalc, DZERO, inmxn, srat_new, rel_perm, alp_ss_new)

    if (st_sim%sim_type >= 0) then
      ! -- Function storage change (stochn)
        call func_stochn(inmxn, alp_ss_new, alp_ss_old, stom)
    end if
#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(inmxn)
        call senrec_rvectv(rel_perm)
    end if
#endif

    ! -- Function connect flow from adjacent cells (connflow)
      call func_connflow(inmxn, conm)

    ! -- Function recharge term (rechterm)
      call func_rechterm(recm)

    ! -- Function well term (wellterm)
      call func_wellterm(welm)

    if (sfla == 0) then
      ! -- Function surface term (surfterm)
        call func_surfterm(inmxn, surf_old, surm)
    end if

    ! -- Function river term (riveterm)
      call func_riveterm(inmxn, rivm)

    ! -- Function lake term (laketerm)
      call func_laketerm(inmxn, lakm)

    ! -- Function sea level term (sealterm)
      call func_sealterm(inmxn, seam)

    if (st_sim%sim_type >= 0) then
      !$omp parallel
      !$omp do private(i)
      do i = 1, ncalc
        stom(i) = stom(i)*delt ; conm(i) = conm(i)*delt
        seam(i) = seam(i)*delt ; welm(i) = welm(i)*delt
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, ncals
        recm(i) = recm(i)*delt ; surm(i) = surm(i)*delt
        rivm(i) = rivm(i)*delt ; lakm(i) = lakm(i)*delt
      end do
      !$omp end do
      !$omp end parallel
    end if

    deallocate(alp_ss_new, alp_ss_old)

  end subroutine calc_mass

  subroutine func_stochn(infstoc, alp_new, alp_old, stofunc)
  !*********************************************************************************************
  ! func_stochn -- Function storage change
  !*********************************************************************************************
    ! -- modules
    use make_cell, only: cell_vol
    use assign_calc, only: read_poro
    use prep_calculation, only: delt_inv
    use allocate_solution, only: srat_old, head_old
    ! -- inout
    real(DP), intent(in) :: infstoc(:), alp_new(:), alp_old(:)
    real(DP), intent(out) :: stofunc(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      ch_stor1(i) = DZERO ; ch_stor2(i) = DZERO ; ch_stor(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, ncalc
      ch_stor1(i) = alp_old(i)*srat_old(i)*head_old(i)
      ch_stor2(i) = alp_new(i)*srat_new(i)*infstoc(i)
      ch_stor(i) = (ch_stor2(i)-ch_stor1(i)) + read_poro(i)*(srat_new(i)-srat_old(i))
      stofunc(i) = -ch_stor(i)*delt_inv*cell_vol(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_stochn

  subroutine func_connflow(infconn, confunc)
  !*********************************************************************************************
  ! func_connflow -- Function connect flow from adjacent cells
  !*********************************************************************************************
    ! -- modules
    use calc_parameter, only: calc_hyd_upwind
    use allocate_solution, only: abyd_conn, hydf_conn
    ! -- inout
    real(DP), intent(in) :: infconn(:)
    real(DP), intent(out) :: confunc(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: sta_ind, end_ind, ind
    real(DP) :: relat, delhead, relp1, relp2
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      conn_flow(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i, j, k, sta_ind, end_ind, ind, delhead, relp1, relp2, relat)
    do i = 1, ncalc
      sta_ind = crs_index(1)%offind(i-1) ; end_ind = crs_index(1)%offind(i)
      do k = 1, end_ind-sta_ind
        ind = sta_ind + k ; j = crs_index(1)%offrow(ind)
        delhead = infconn(j) - infconn(i)
        relp1 = rel_perm(i) ; relp2 = rel_perm(j)

        ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
          call calc_hyd_upwind(-delhead, relp1, relp2, relat)

        conn_flow(i) = conn_flow(i) + hydf_conn(ind)*relat*abyd_conn(ind)*delhead
      end do
      confunc(i) = conn_flow(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_connflow

  subroutine func_rechterm(recfunc)
  !*********************************************************************************************
  ! func_rechterm -- Function recharge term
  !*********************************************************************************************
    ! -- modules
    use calc_boundary, only: rech2cals, calc_rech
    use set_boundary, only: rech_num
    ! -- inout
    real(DP), intent(inout) :: recfunc(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, s)
    do i = 1, rech_num
      s = rech2cals(i)
      recfunc(s) = calc_rech(i)
    end do
    !$omp end parallel do

  end subroutine func_rechterm

  subroutine func_wellterm(welfunc)
  !*********************************************************************************************
  ! func_wellterm -- Function well term
  !*********************************************************************************************
    ! -- modules
    use set_condition, only: well_index, well_conn
    use assign_boundary, only: calc_well
    use set_boundary, only: well_num
    ! -- inout
    real(DP), intent(inout) :: welfunc(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: sta_wind, end_wind, wind
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, sta_wind, end_wind, wind)
    do i = 1, well_num
      sta_wind = well_index(i-1) ; end_wind = well_index(i)
      do k = 1, end_wind-sta_wind
        wind = sta_wind + k
        j = well_conn(wind)
        welfunc(j) = calc_well(j)
      end do
    end do
    !$omp end parallel do

  end subroutine func_wellterm

  subroutine func_surfterm(infsurf, inshead, surfunc)
  !*********************************************************************************************
  ! func_surfterm -- Function surface term
  !*********************************************************************************************
    ! -- modules
!    use make_cell, only: surf_elev
    use set_condition, only: abyd_surf
    use assign_calc, only: surf_bott, surf_reli, surf_parm
    use prep_calculation, only: surf_top
    use allocate_solution, only: surf_rati
    ! -- inout
    real(DP), intent(in) :: infsurf(:), inshead(:)
    real(DP), intent(inout) :: surfunc(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      delh_s(i) = DZERO ; elev_rati(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, ncals
!      if (infsurf(i) >= surf_elev(i)) then
!        if (inshead(i) > surf_bott(i)) then
!          delh_s(i) = inshead(i) - infsurf(i)
!        else
!          delh_s(i) = surf_bott(i) - infsurf(i)
!        end if
!        surf_rati(i) = DONE
!      else if (inshead(i) >= surf_elev(i)) then
!        delh_s(i) = inshead(i) - surf_elev(i)
!        surf_rati(i) = DONE
!      else
!!        delh_s(i) = inshead(i) - infsurf(i)
!        delh_s(i) = DZERO
!        surf_rati(i) = DONE
!      end if
!      if (surf_reli(i) == DZERO) then
!        surf_rati(i) = DONE ; delh_s(i) = DZERO
!      else if (infsurf(i) >= surf_top(i)) then
      if (infsurf(i) >= surf_top(i)) then
        if (inshead(i) > surf_top(i)) then
          delh_s(i) = inshead(i) - infsurf(i)
        else
          delh_s(i) = surf_top(i) - infsurf(i)
        end if
        surf_rati(i) = DONE
      else if (infsurf(i) >= surf_bott(i)) then
        if (inshead(i) > surf_bott(i)) then
          delh_s(i) = inshead(i) - infsurf(i)
        else
          delh_s(i) = surf_bott(i) - infsurf(i)
        end if
        if (surf_reli(i) /= DZERO) then
          elev_rati(i) = (infsurf(i)-surf_bott(i))/surf_reli(i)
          if (elev_rati(i) > DONE) then
            elev_rati(i) = DONE
          end if
          surf_rati(i) = elev_rati(i)**surf_parm(i)
        else
          surf_rati(i) = DONE
        end if
      else if (inshead(i) > surf_bott(i)) then
        delh_s(i) = inshead(i) - surf_bott(i)
        if (surf_reli(i) /= DZERO) then
          elev_rati(i) = (inshead(i)-surf_bott(i))/surf_reli(i)
          if (elev_rati(i) > DONE) then
            elev_rati(i) = DONE
          end if
          surf_rati(i) = elev_rati(i)**surf_parm(i)
        else
          surf_rati(i) = DONE
        end if
      else
        delh_s(i) = DZERO
      end if

      surfunc(i) = hydf_surf(i)*abyd_surf(i)*delh_s(i)*rel_perm(i)*surf_rati(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_surfterm

  subroutine func_riveterm(infrive, rivfunc)
  !*********************************************************************************************
  ! func_riveterm -- Function river term
  !*********************************************************************************************
    ! -- modules
    use calc_boundary, only: rive2cals, rive_head, rive_bott
    use set_boundary, only: abyd_rive
    ! -- inout
    real(DP), intent(in) :: infrive(:)
    real(DP), intent(inout) :: rivfunc(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, rive_num
      delh_r(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i, s)
    do i = 1, rive_num
      s = rive2cals(i)
      if (infrive(s) >= rive_bott(i)) then
        delh_r(i) = rive_head(i) - infrive(s)
      else if (rive_head(i) > rive_bott(i)) then
        delh_r(i) = rive_head(i) - rive_bott(i)
      else
        delh_r(i) = DZERO
      end if

      rivfunc(s) = hydf_surf(s)*abyd_rive(i)*delh_r(i)*rel_perm(s)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_riveterm

  subroutine func_laketerm(inflake, lakfunc)
  !*********************************************************************************************
  ! func_laketerm -- Function lake term
  !*********************************************************************************************
    ! -- modules
    use calc_boundary, only: lake2cals, lake_head, lake_bott
    use set_boundary, only: abyd_lake
    ! -- inout
    real(DP), intent(in) :: inflake(:)
    real(DP), intent(inout) :: lakfunc(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, lake_num
      delh_l(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i, s)
    do i = 1, lake_num
      s = lake2cals(i)
      if (inflake(s) >= lake_bott(i)) then
        delh_l(i) = lake_head(i) - inflake(s)
      else if (lake_head(i) > lake_bott(i)) then
        delh_l(i) = lake_head(i) - lake_bott(i)
      else
        delh_l(i) = DZERO
      end if

      lakfunc(s) = hydf_surf(s)*abyd_lake(i)*delh_l(i)*rel_perm(s)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_laketerm

  subroutine func_sealterm(infseal, seafunc)
  !*********************************************************************************************
  ! func_sealterm -- Function sea level term
  !*********************************************************************************************
    ! -- modules
    use assign_boundary, only: read_seal
    use allocate_solution, only: seal2calc, seal2seal, hydf_seal, abyd_seal
    ! -- inout
    real(DP), intent(in) :: infseal(:)
    real(DP), intent(inout) :: seafunc(:)
    ! -- local
    integer(I4) :: i, c, s
    real(DP) :: delhead
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, nseal
      seal_flow(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i, c, s, delhead)
    do i = 1, nseal
      c = seal2calc(i) ; s = seal2seal(i)
      delhead = read_seal(s) - infseal(c)
      seal_flow(i) = hydf_seal(i)*rel_perm(c)*abyd_seal(i)*delhead
    end do
    !$omp end do

    !$omp do private(i, c)
    do i = 1, nseal
      c = seal2calc(i)
      !$omp atomic
      seafunc(c) = seafunc(c) + seal_flow(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_sealterm

  subroutine calc_vecjacf(vjlevel, injx, injvec, outjvec)
  !*********************************************************************************************
  ! calc_vecjacf -- Calculate vector by jacobi-free
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: MACHI_EPS
#ifdef MPI_MSG
    use mpi_utility, only: mpisum_val
#endif
    ! -- inout
    integer(I4), intent(in) :: vjlevel
    real(DP), intent(in) :: injvec(:)
    real(DP), intent(inout) :: injx(:)
    real(DP), intent(out) :: outjvec(:)
    ! -- local
    integer(I4) :: i
    integer(I4) :: vj_num, vj_regnum
    real(DP) :: eps, eps_inv, l2_x, l2_v, l1_v, sign
    real(DP), allocatable :: jcvec(:), tempf1(:), tempf2(:)
#ifdef MPI_MSG
    real(DP) :: sum_l2
#endif
    !-------------------------------------------------------------------------------------------
    vj_num = crs_index(vjlevel)%unknow
    vj_regnum = size(injx)

    allocate(jcvec(vj_regnum), tempf1(vj_num), tempf2(vj_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, vj_regnum
      jcvec(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, vj_num
      tempf1(i) = DZERO ; tempf2(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    ! -- Calculate function value (func)
      call calc_func(injx, tempf1)

    ! Brown and Saad version
    l2_v = DZERO ; l2_x = DZERO ; l1_v = DZERO
    !$omp parallel
    !$omp do private(i) reduction(+:l2_v, l2_x)
    do i = 1, vj_num
      l2_v = l2_v + injvec(i)*injvec(i)
      l2_x = l2_x + injx(i)*injvec(i)
    end do
    !$omp end do
    !$omp do private(i) reduction(+:l1_v)
    do i = 1, vj_num
      l1_v = l1_v + abs(injvec(i))
    end do
    !$omp end do
    !$omp end parallel

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_v, "jacobian-free v l2-norm", sum_l2)
      l2_v = sum_l2
      ! -- Sum value for MPI (val)
        call mpisum_val(l2_x, "jacobian-free x l2-norm", sum_l2)
      l2_x = sum_l2
      ! -- Sum value for MPI (val)
        call mpisum_val(l1_v, "jacobian-free v l1-norm", sum_l2)
      l1_v = sum_l2
    end if
#endif

    if (l2_x < DZERO) then
      sign = -DONE
    else
      sign = DONE
    end if

    eps = sign*sqrt(MACHI_EPS)*max(abs(l2_x),l1_v)/l2_v
    eps_inv = DONE/eps

    !$omp parallel do private(i)
    do i = 1, vj_num
      jcvec(i) = injx(i) + eps*injvec(i)
    end do
    !$omp end parallel do

    ! -- Calculate function value (func)
      call calc_func(jcvec, tempf2)

    !$omp parallel do private(i)
    do i = 1, vj_num
      outjvec(i) = -(tempf2(i)-tempf1(i))*eps_inv
    end do
    !$omp end parallel do

    deallocate(jcvec, tempf1, tempf2)

  end subroutine calc_vecjacf

end module calc_function
