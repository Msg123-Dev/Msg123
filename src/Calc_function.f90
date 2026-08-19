module calc_function
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE
  use initial_module, only: st_sim
  use set_cell, only: ncalc, ncals
  use set_condition, only: st_hydr, st_bcnd
  use prep_calculation, only: st_time
  use assign_boundary, only: st_forc
  use calc_parameter, only: calc_srat_rperm
  use allocate_solution, only: crs_index
#ifdef MPI_MSG
  use utility_module, only: st_mpi
  use mpi_solve, only: senrec_rvectv
#endif

  implicit none
  private
  public :: allocate_calfun, calc_func, calc_mass, calc_vecjacf
  public :: func_rechterm, func_wellterm, func_surfterm, func_riveterm
  public :: func_laketerm, func_sealterm
  ! -- local
  real(DP), allocatable :: stof(:), conf(:), welf(:), seaf(:)
  real(DP), allocatable :: funcvs(:)
  real(DP), allocatable :: recf(:), surf(:), rivf(:), lakf(:)
  real(DP), allocatable :: wk_stor(:)
  real(DP), allocatable :: conn_flow(:)
  real(DP), allocatable :: delh_s(:), elev_rati(:)
  real(DP), allocatable :: delh_r(:), delh_l(:), seal_flow(:)
  real(DP), allocatable :: jcvec(:), tempf1(:), tempf2(:)

  contains

  subroutine allocate_calfun()
  !*********************************************************************************************
  ! allocate_calfun -- Allocate for calculate function value
  !*********************************************************************************************
    ! -- modules
    use allocate_solution, only: nreg_num

    ! -- inout

    ! -- local
    !-------------------------------------------------------------------------------------------
    allocate(stof(ncalc), conf(ncalc), welf(ncalc), seaf(ncalc))
    allocate(funcvs(ncalc))
    allocate(recf(ncals), surf(ncals), rivf(ncals), lakf(ncals))
    allocate(wk_stor(ncalc))
    allocate(conn_flow(ncalc))
    allocate(delh_s(ncals), elev_rati(ncals))
    allocate(delh_r(st_bcnd%rive_num), delh_l(st_bcnd%lake_num), seal_flow(st_bcnd%seal_num))
    allocate(jcvec(nreg_num), tempf1(ncalc), tempf2(ncalc))

  end subroutine allocate_calfun

  subroutine calc_func(stold, stnew, surfh, infx, snew, rperm, surfr, funcv)
  !*********************************************************************************************
  ! calc_func -- Calculate function value
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    real(DP), intent(in) :: stold(:), surfh(:)
    real(DP), intent(inout) :: stnew(:), infx(:), snew(:), rperm(:), surfr(:)
    real(DP), intent(out) :: funcv(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      stof(i) = DZERO ; conf(i) = DZERO ; welf(i) = DZERO ; seaf(i) = DZERO
      funcvs(i) = DZERO ; funcv(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncals
      recf(i) = DZERO ; surf(i) = DZERO ; rivf(i) = DZERO ; lakf(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    ! -- Calculate saturation and relative permeability (srat_rperm)
      call calc_srat_rperm(ncalc, DZERO, infx, snew, rperm, stnew)

    if (st_sim%sim_type >= 0) then
      ! -- Function storage change (stochn)
        call func_stochn(stnew, stold, stof)
    end if

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(infx)
        call senrec_rvectv(rperm)
    end if
#endif

    ! -- Function connect flow from adjacent cells (connflow)
      call func_connflow(infx, rperm, conf)

    ! -- Function recharge term (rechterm)
      call func_rechterm(recf)

    ! -- Function well term (wellterm)
      call func_wellterm(welf)

    ! -- Function surface term (surfterm)
      call func_surfterm(infx, surfh, rperm, surf, surfr)

    ! -- Function river term (riveterm)
      call func_riveterm(infx, rperm, rivf)

    ! -- Function lake term (laketerm)
      call func_laketerm(infx, rperm, lakf)

    ! -- Function sea level term (sealterm)
      call func_sealterm(infx, rperm, seaf)

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

  subroutine calc_mass(sfla, stold, surf_old, hnew, stom, conm, seam, welm, recm,&
                       surm, rivm, lakm, snew, rperm, surfr)
  !*********************************************************************************************
  ! calc_mass -- Calculate function value for massbalance
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    integer(I4), intent(in) :: sfla
    real(DP), intent(in) :: stold(:), surf_old(:)
    real(DP), intent(inout) :: hnew(:)
    real(DP), intent(inout) :: stom(:), conm(:), seam(:), welm(:)
    real(DP), intent(inout) :: recm(:), surm(:), rivm(:), lakm(:)
    real(DP), intent(inout) :: snew(:), rperm(:), surfr(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    ! -- Calculate saturation and relative permeability (srat_rperm)
      call calc_srat_rperm(ncalc, DZERO, hnew, snew, rperm, wk_stor)

    if (st_sim%sim_type >= 0) then
      ! -- Function storage change (stochn)
        call func_stochn(wk_stor, stold, stom)
    end if
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(hnew)
        call senrec_rvectv(rperm)
    end if
#endif

    ! -- Function connect flow from adjacent cells (connflow)
      call func_connflow(hnew, rperm, conm)

    ! -- Function recharge term (rechterm)
      call func_rechterm(recm)

    ! -- Function well term (wellterm)
      call func_wellterm(welm)

    if (sfla == 0) then
      ! -- Function surface term (surfterm)
        call func_surfterm(hnew, surf_old, rperm, surm, surfr)
    end if

    ! -- Function river term (riveterm)
      call func_riveterm(hnew, rperm, rivm)

    ! -- Function lake term (laketerm)
      call func_laketerm(hnew, rperm, lakm)

    ! -- Function sea level term (sealterm)
      call func_sealterm(hnew, rperm, seam)

    if (st_sim%sim_type >= 0) then
      !$omp parallel
      !$omp do private(i)
      do i = 1, ncalc
        stom(i) = stom(i)*st_time%delt ; conm(i) = conm(i)*st_time%delt
        seam(i) = seam(i)*st_time%delt ; welm(i) = welm(i)*st_time%delt
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, ncals
        recm(i) = recm(i)*st_time%delt ; surm(i) = surm(i)*st_time%delt
        rivm(i) = rivm(i)*st_time%delt ; lakm(i) = lakm(i)*st_time%delt
      end do
      !$omp end do
      !$omp end parallel
    end if


  end subroutine calc_mass

  subroutine func_stochn(stnew, stold, stofunc)
  !*********************************************************************************************
  ! func_stochn -- Function storage change
  !*********************************************************************************************
    ! -- modules
    use make_cell, only: st_geom
    ! -- inout
    real(DP), intent(in) :: stnew(:), stold(:)
    real(DP), intent(out) :: stofunc(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, ncalc
      stofunc(i) = -(stnew(i)-stold(i))*st_time%delt_inv*st_geom%cell_vol(i)
    end do
    !$omp end parallel do

  end subroutine func_stochn

  subroutine func_connflow(infconn, rperm, confunc)
  !*********************************************************************************************
  ! func_connflow -- Function connect flow from adjacent cells
  !*********************************************************************************************
    ! -- modules
    use calc_parameter, only: calc_hyd_upwind
    ! -- inout
    real(DP), intent(in) :: infconn(:), rperm(:)
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
        relp1 = rperm(i) ; relp2 = rperm(j)

        ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
          call calc_hyd_upwind(-delhead, relp1, relp2, relat)

        conn_flow(i) = conn_flow(i) + st_hydr%hydf_conn(ind)*relat*st_hydr%abyd_conn(ind)*&
                       delhead
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
    ! -- inout
    real(DP), intent(inout) :: recfunc(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, s)
    do i = 1, st_bcnd%rech_num
      s = st_bcnd%rech2cals(i)
      recfunc(s) = st_forc%calc_rech(i)
    end do
    !$omp end parallel do

  end subroutine func_rechterm

  subroutine func_wellterm(welfunc)
  !*********************************************************************************************
  ! func_wellterm -- Function well term
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    real(DP), intent(inout) :: welfunc(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: sta_wind, end_wind, wind
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, sta_wind, end_wind, wind)
    do i = 1, st_bcnd%well_num
      sta_wind = st_bcnd%well_index(i-1) ; end_wind = st_bcnd%well_index(i)
      do k = 1, end_wind-sta_wind
        wind = sta_wind + k
        j = st_bcnd%well_conn(wind)
        welfunc(j) = st_forc%calc_well(j)
      end do
    end do
    !$omp end parallel do

  end subroutine func_wellterm

  subroutine func_surfterm(infsurf, inshead, rperm, surfunc, surfr)
  !*********************************************************************************************
  ! func_surfterm -- Function surface term
  !*********************************************************************************************
    ! -- modules
!    use make_cell, only: surf_elev
    ! -- inout
    real(DP), intent(in) :: infsurf(:), inshead(:), rperm(:)
    real(DP), intent(inout) :: surfunc(:), surfr(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      delh_s(i) = DZERO ; elev_rati(i) = DZERO ; surfr(i) = DONE
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
!        surfr(i) = DONE
!      else if (inshead(i) >= surf_elev(i)) then
!        delh_s(i) = inshead(i) - surf_elev(i)
!        surfr(i) = DONE
!      else
!!        delh_s(i) = inshead(i) - infsurf(i)
!        delh_s(i) = DZERO
!        surfr(i) = DONE
!      end if
!      if (surf_reli(i) == DZERO) then
!        surfr(i) = DONE ; delh_s(i) = DZERO
!      else if (infsurf(i) >= surf_top(i)) then
      if (infsurf(i) >= st_hydr%surf_top(i)) then
        if (inshead(i) > st_hydr%surf_top(i)) then
          delh_s(i) = inshead(i) - infsurf(i)
        else
          delh_s(i) = st_hydr%surf_top(i) - infsurf(i)
        end if
        surfr(i) = DONE
      else if (infsurf(i) >= st_hydr%surf_bott(i)) then
        if (inshead(i) > st_hydr%surf_bott(i)) then
          delh_s(i) = inshead(i) - infsurf(i)
        else
          delh_s(i) = st_hydr%surf_bott(i) - infsurf(i)
        end if
        if (st_hydr%surf_reli(i) /= DZERO) then
          elev_rati(i) = (infsurf(i)-st_hydr%surf_bott(i))/st_hydr%surf_reli(i)
          if (elev_rati(i) > DONE) then
            elev_rati(i) = DONE
          end if
          surfr(i) = elev_rati(i)**st_hydr%surf_parm(i)
        else
          surfr(i) = DONE
        end if
      else if (inshead(i) > st_hydr%surf_bott(i)) then
        delh_s(i) = inshead(i) - st_hydr%surf_bott(i)
        if (st_hydr%surf_reli(i) /= DZERO) then
          elev_rati(i) = (inshead(i)-st_hydr%surf_bott(i))/st_hydr%surf_reli(i)
          if (elev_rati(i) > DONE) then
            elev_rati(i) = DONE
          end if
          surfr(i) = elev_rati(i)**st_hydr%surf_parm(i)
        else
          surfr(i) = DONE
        end if
      else
        delh_s(i) = DZERO
      end if

      surfunc(i) = st_hydr%hydf_surf(i)*st_hydr%abyd_surf(i)*delh_s(i)*rperm(i)*surfr(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_surfterm

  subroutine func_riveterm(infrive, rperm, rivfunc)
  !*********************************************************************************************
  ! func_riveterm -- Function river term
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    real(DP), intent(in) :: infrive(:), rperm(:)
    real(DP), intent(inout) :: rivfunc(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%rive_num
      delh_r(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i, s)
    do i = 1, st_bcnd%rive_num
      s = st_bcnd%rive2cals(i)
      if (infrive(s) >= st_forc%rive_bott(i)) then
        delh_r(i) = st_forc%rive_head(i) - infrive(s)
      else if (st_forc%rive_head(i) > st_forc%rive_bott(i)) then
        delh_r(i) = st_forc%rive_head(i) - st_forc%rive_bott(i)
      else
        delh_r(i) = DZERO
      end if

      rivfunc(s) = st_hydr%hydf_surf(s)*st_forc%abyd_rive(i)*delh_r(i)*rperm(s)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_riveterm

  subroutine func_laketerm(inflake, rperm, lakfunc)
  !*********************************************************************************************
  ! func_laketerm -- Function lake term
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    real(DP), intent(in) :: inflake(:), rperm(:)
    real(DP), intent(inout) :: lakfunc(:)
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%lake_num
      delh_l(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i, s)
    do i = 1, st_bcnd%lake_num
      s = st_bcnd%lake2cals(i)
      if (inflake(s) >= st_forc%lake_bott(i)) then
        delh_l(i) = st_forc%lake_head(i) - inflake(s)
      else if (st_forc%lake_head(i) > st_forc%lake_bott(i)) then
        delh_l(i) = st_forc%lake_head(i) - st_forc%lake_bott(i)
      else
        delh_l(i) = DZERO
      end if

      lakfunc(s) = st_hydr%hydf_surf(s)*st_forc%abyd_lake(i)*delh_l(i)*rperm(s)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_laketerm

  subroutine func_sealterm(infseal, rperm, seafunc)
  !*********************************************************************************************
  ! func_sealterm -- Function sea level term
  !*********************************************************************************************
    ! -- modules
    ! -- inout
    real(DP), intent(in) :: infseal(:), rperm(:)
    real(DP), intent(inout) :: seafunc(:)
    ! -- local
    integer(I4) :: i, c, s
    real(DP) :: delhead
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%seal_num
      seal_flow(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i, c, s, delhead)
    do i = 1, st_bcnd%seal_num
      c = st_bcnd%seal2calc(i) ; s = st_bcnd%seal2seal(i)
      delhead = st_forc%read_seal(s) - infseal(c)
      seal_flow(i) = st_hydr%hydf_seal(i)*rperm(c)*st_hydr%abyd_seal(i)*delhead
    end do
    !$omp end do

    !$omp do private(i, c)
    do i = 1, st_bcnd%seal_num
      c = st_bcnd%seal2calc(i)
      !$omp atomic
      seafunc(c) = seafunc(c) + seal_flow(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine func_sealterm

  subroutine calc_vecjacf(vjlevel, injvec, stold, stnew, surfh, injx, snew, rperm, surfr,&
                          outjvec)
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
    real(DP), intent(in) :: injvec(:), stold(:), surfh(:)
    real(DP), intent(inout) :: stnew(:), injx(:), snew(:), rperm(:), surfr(:)
    real(DP), intent(out) :: outjvec(:)
    ! -- local
    integer(I4) :: i
    integer(I4) :: vj_num, vj_regnum
    real(DP) :: eps, eps_inv, l2_x, l2_v, l1_v, sign
#ifdef MPI_MSG
    real(DP) :: sum_l2
#endif
    !-------------------------------------------------------------------------------------------
    vj_num = crs_index(vjlevel)%unknow
    vj_regnum = size(injx)

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
      call calc_func(stold, stnew, surfh, injx, snew, rperm, surfr, tempf1)

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
    if (st_mpi%totn /= 1) then
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
      call calc_func(stold, stnew, surfh, jcvec, snew, rperm, surfr, tempf2)

    !$omp parallel do private(i)
    do i = 1, vj_num
      outjvec(i) = -(tempf2(i)-tempf1(i))*eps_inv
    end do
    !$omp end parallel do

  end subroutine calc_vecjacf

end module calc_function
