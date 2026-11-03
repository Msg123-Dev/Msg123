module calc_parameter
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: DZERO, DONE, DHALF
  use initial_module, only: st_schm
  use set_condition, only: st_hydr

  implicit none
  private
  public :: calc_srat_rperm, calc_kr_cell, calc_hyd_geoharm, calc_hyd_upwind
  public :: calc_hyd_harm, calc_hyd_geo, calc_hyd_arith

  ! -- local

  contains

  subroutine calc_srat_rperm(num, pres, srat, rperm, stor, dkr_dpsi, dstor_dpsi)
  !*********************************************************************************************
  ! calc_srat_rperm -- Calculate saturation and relative permeability
  !*********************************************************************************************
    ! -- modules
    use make_cell, only: st_geom
    ! -- inout
    integer(I4), intent(in) :: num
    real(DP), intent(in) :: pres(:)
    real(DP), intent(out), optional :: srat(:)
    real(DP), intent(out), optional :: rperm(:)
    real(DP), intent(out), optional :: stor(:)
    real(DP), intent(out), optional :: dkr_dpsi(:)
    real(DP), intent(out), optional :: dstor_dpsi(:)
    ! -- local
    integer(I4) :: i
    real(DP) :: retm, beta, theta, kr, se, stor_val, srat_val
    real(DP) :: phead, kr_phead, dkr
    real(DP) :: pors_val, resi_val, dtheta_dpsi, dstor_val
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, pors_val, resi_val, retm, phead, dkr, dstor_val) &
    !$omp             private(beta, se, theta, srat_val, stor_val, dtheta_dpsi, kr_phead, kr)
    do i = 1, num
      pors_val = st_hydr%read_pors(i) ; resi_val = st_hydr%read_resi(i)
      retm = DONE - DONE/st_hydr%read_vann(i)
      phead = pres(i) - st_geom%cell_top(i)
      ! -- Storage and saturation ratio at the cell top
      if (phead < DZERO) then
        beta = abs(phead*st_hydr%read_vana(i))**st_hydr%read_vann(i)
        se = (DONE+beta)**(-retm)
        theta = resi_val + (pors_val-resi_val)*se
        srat_val = theta/pors_val
        if (st_schm%stor_type == 1) then
          stor_val = theta + st_hydr%read_spst(i)*phead*srat_val
        else
          stor_val = theta
        end if
        if (present(dstor_dpsi)) then
          dtheta_dpsi = (pors_val-resi_val)*se*retm*st_hydr%read_vann(i)*beta&
                        /(abs(phead)*(DONE+beta))
          if (st_schm%stor_type == 1) then
            dstor_dpsi(i) = dtheta_dpsi*(DONE+st_hydr%read_spst(i)*phead/pors_val)&
                            + st_hydr%read_spst(i)*srat_val
          else
            dstor_dpsi(i) = dtheta_dpsi
          end if
        end if
      else
        srat_val = DONE
        stor_val = pors_val + st_hydr%read_spst(i)*phead
        if (present(dstor_dpsi)) then
          dstor_dpsi(i) = st_hydr%read_spst(i)
        end if
      end if

      ! -- Storage averaged over the cell thickness given by stor_avg_type
      if (st_schm%stor_avg_type == 1 .and. st_geom%cell_top(i) > st_geom%cell_bot(i)) then
        call calc_stor_avg(i, phead, st_geom%cell_top(i)-st_geom%cell_bot(i), stor_val,&
                           dstor_val)
        if (present(dstor_dpsi)) then
          dstor_dpsi(i) = dstor_val
        end if
      end if

      ! -- Relative permeability at the head position given by krpos_type
      if (st_schm%krpos_type == 1) then
        kr_phead = pres(i) - st_geom%cell_cent(i)
      else
        kr_phead = phead
      end if
      ! -- Relative permeability of the cell (kr_cell)
        call calc_kr_cell(i, kr_phead, kr, dkr)

      if (present(rperm)) then
        rperm(i) = kr
      end if
      if (present(dkr_dpsi)) then
        dkr_dpsi(i) = dkr
      end if
      if (present(srat)) then
        srat(i) = srat_val
      end if
      if (present(stor)) then
        stor(i) = stor_val
      end if
    end do
    !$omp end parallel do

  end subroutine calc_srat_rperm

  subroutine calc_stor_avg(i, ptop, thick, stor, dstor)
  !*********************************************************************************************
  ! calc_stor_avg -- Calculate storage averaged over the cell thickness
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: i
    real(DP), intent(in) :: ptop, thick
    real(DP), intent(inout) :: stor
    real(DP), intent(out) :: dstor
    ! -- local
    integer(I4), parameter :: GAUS_NUM = 16
    integer(I4) :: k
    real(DP), parameter :: GAUS_NODE(GAUS_NUM) = [&
      -9.8940093499164994E-01_DP, -9.4457502307323260E-01_DP, -8.6563120238783164E-01_DP,&
      -7.5540440835500311E-01_DP, -6.1787624440264377E-01_DP, -4.5801677765722737E-01_DP,&
      -2.8160355077925892E-01_DP, -9.5012509837637441E-02_DP, 9.5012509837637441E-02_DP,&
      2.8160355077925892E-01_DP, 4.5801677765722737E-01_DP, 6.1787624440264377E-01_DP,&
      7.5540440835500311E-01_DP, 8.6563120238783164E-01_DP, 9.4457502307323260E-01_DP,&
      9.8940093499164994E-01_DP]
    real(DP), parameter :: GAUS_WGHT(GAUS_NUM) = [&
      2.7152459411754055E-02_DP, 6.2253523938647609E-02_DP, 9.5158511682492994E-02_DP,&
      1.2462897125553395E-01_DP, 1.4959598881657665E-01_DP, 1.6915651939500262E-01_DP,&
      1.8260341504492361E-01_DP, 1.8945061045506859E-01_DP, 1.8945061045506859E-01_DP,&
      1.8260341504492361E-01_DP, 1.6915651939500262E-01_DP, 1.4959598881657665E-01_DP,&
      1.2462897125553395E-01_DP, 9.5158511682492994E-02_DP, 6.2253523938647609E-02_DP,&
      2.7152459411754055E-02_DP]
    real(SP) :: vana, vann, spst
    real(DP) :: pors_val, resi_val, retm, pbot, punsat, thick_uns, beta, se, theta, elog
    real(DP) :: stor_bot, stor_sum, slog_arg, slog_rat, slog_bot, slog_mid, slog_half, gsum
    !-------------------------------------------------------------------------------------------
    vana = st_hydr%read_vana(i) ; vann = st_hydr%read_vann(i) ; spst = st_hydr%read_spst(i)
    pors_val = st_hydr%read_pors(i) ; resi_val = st_hydr%read_resi(i)
    retm = DONE - DONE/vann
    pbot = ptop + thick
    if (ptop >= DZERO) then
      ! -- Saturated over the whole thickness
      stor = pors_val + spst*DHALF*(ptop+pbot) ; dstor = spst
    else
      ! -- Storage at the cell bottom and its difference from the top
      if (pbot >= DZERO) then
        stor_bot = pors_val + spst*pbot
      else
        beta = abs(pbot*vana)**vann
        se = (DONE+beta)**(-retm)
        theta = resi_val + (pors_val-resi_val)*se
        if (st_schm%stor_type == 1) then
          stor_bot = theta + spst*pbot*(theta/pors_val)
        else
          stor_bot = theta
        end if
      end if
      dstor = (stor_bot-stor)/thick

      ! -- Saturated part and residual part
      if (pbot > DZERO) then
        punsat = DZERO ; thick_uns = -ptop
        stor_sum = pors_val*pbot + spst*DHALF*pbot*pbot
      else
        punsat = pbot ; thick_uns = thick
        stor_sum = DZERO
      end if
      stor_sum = stor_sum + resi_val*thick_uns
      if (st_schm%stor_type == 1) then
        stor_sum = stor_sum + spst*(resi_val/pors_val)*DHALF*thick_uns*(punsat+ptop)
      end if

      ! -- Effective saturation part by Gauss-Legendre in log(1 - vana*phead)
      slog_rat = vana*thick_uns/(DONE-vana*punsat) ; slog_arg = DONE + slog_rat
      if (slog_arg > DONE) then
        slog_half = DHALF*log(slog_arg)*(slog_rat/(slog_arg-DONE))
      else
        slog_half = DHALF*slog_rat
      end if
      slog_bot = log(DONE-vana*punsat) ; slog_mid = slog_bot + slog_half
      gsum = DZERO
      do k = 1, GAUS_NUM
        elog = exp(slog_mid+slog_half*GAUS_NODE(k))
        beta = (elog-DONE)**vann
        se = (DONE+beta)**(-retm)
        if (st_schm%stor_type == 1) then
          gsum = gsum + GAUS_WGHT(k)*se*elog*(DONE+spst*(DONE-elog)/(vana*pors_val))
        else
          gsum = gsum + GAUS_WGHT(k)*se*elog
        end if
      end do
      stor_sum = stor_sum + (pors_val-resi_val)*slog_half*gsum/vana
      stor = stor_sum/thick
    end if

  end subroutine calc_stor_avg

  subroutine calc_kr_cell(i, phead, kr, dkr)
  !*********************************************************************************************
  ! calc_kr_cell -- Calculate relative permeability of a cell at a pressure head
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: i
    real(DP), intent(in) :: phead
    real(DP), intent(out) :: kr
    real(DP), intent(out), optional :: dkr
    ! -- local
    real(SP) :: vana, vann
    real(DP) :: retm, phead0, kr_lin, kr_grad, dkr_val
    !-------------------------------------------------------------------------------------------
    vana = st_hydr%read_vana(i) ; vann = st_hydr%read_vann(i)
    retm = DONE - DONE/vann
    phead0 = -st_schm%krlin_head
    if (phead >= DZERO) then
      kr = DONE ; dkr_val = DZERO
    else if (phead > phead0) then
      ! -- Relative permeability and its slope at the linear bridge head
      call calc_kr_vgm(phead0, vana, vann, retm, kr_lin, dkr_val)
      kr_grad = (DONE-kr_lin)/st_schm%krlin_head
      kr = kr_lin + (phead-phead0)*kr_grad ; dkr_val = kr_grad
    else
      call calc_kr_vgm(phead, vana, vann, retm, kr, dkr_val)
    end if
    if (present(dkr)) then
      dkr = dkr_val
    end if

  end subroutine calc_kr_cell

  subroutine calc_kr_vgm(phead, vana, vann, retm, kr, dkr)
  !*********************************************************************************************
  ! calc_kr_vgm -- Calculate relative permeability by van Genuchten and Mualem
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DTWO
    ! -- inout
    real(DP), intent(in) :: phead, retm
    real(SP), intent(in) :: vana, vann
    real(DP), intent(out) :: kr, dkr
    ! -- local
    real(DP) :: beta, se, term
    !-------------------------------------------------------------------------------------------
    beta = abs(phead*vana)**vann
    se = (DONE+beta)**(-retm)
    term = DONE-(DONE-se**(DONE/retm))**retm
    kr = se**(DHALF)*term**DTWO
    dkr = -(retm*vann/((DONE+beta)*phead))*sqrt(se)*term*(DHALF*term*beta + DTWO*(DONE-term))

  end subroutine calc_kr_vgm

  subroutine calc_hyd_geoharm(srato1, srato2, hyd_c1, hyd_c2, d1, d2, hyd_gh)
  !*********************************************************************************************
  ! calc_hyd_geoharm -- Calculate hydradulic conductivity by geometric and harmonic
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: srato1, srato2, hyd_c1, hyd_c2, d1, d2
    real(DP), intent(out) :: hyd_gh
    ! -- local
    real(DP) :: c1, c2, d12, cd
    !-------------------------------------------------------------------------------------------
    d12 = d1 + d2

    if (srato1 == DONE .and. srato2 == DONE) then !harmonic mean
      c1 = hyd_c1*d12
      c2 = hyd_c2*d12
      cd = DONE/(d1*c2+d2*c1)
      hyd_gh = (c1*c2)*cd
    else !geometric mean
      d12 = DONE/d12
      c1 = d1*log(hyd_c1)
      c2 = d2*log(hyd_c2)
      hyd_gh = exp((c1+c2)*d12)
    end if

  end subroutine calc_hyd_geoharm

  subroutine calc_hyd_upwind(head_dif, hyd_c1, hyd_c2, hyd_upwind)
  !*********************************************************************************************
  ! calc_hyd_upwind -- Calculate hydradulic conductivity by upwind
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: head_dif, hyd_c1, hyd_c2
    real(DP), intent(out) :: hyd_upwind
    ! -- local

    !-------------------------------------------------------------------------------------------
    if (head_dif > DZERO) then
      hyd_upwind = hyd_c1
    else if (head_dif < DZERO) then
      hyd_upwind = hyd_c2
    else
      hyd_upwind = DHALF*(hyd_c1 + hyd_c2)
    end if

  end subroutine calc_hyd_upwind

  subroutine calc_hyd_harm(hyd_c1, hyd_c2, d1, d2, hyd_harm)
  !*********************************************************************************************
  ! calc_hyd_harm -- Calculate hydradulic conductivity by harmonic mean
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: hyd_c1, hyd_c2, d1, d2
    real(DP), intent(out) :: hyd_harm
    ! -- local
    real(DP) :: c1, c2, d12, cd, cd12
    !-------------------------------------------------------------------------------------------
    d12 = d1 + d2

    c1 = hyd_c1*d12
    c2 = hyd_c2*d12
    cd12 = d1*c2+d2*c1
    cd = DONE/cd12
    hyd_harm = c1*c2*cd

  end subroutine calc_hyd_harm

  subroutine calc_hyd_geo(hyd_c1, hyd_c2, d1, d2, hyd_geo)
  !*********************************************************************************************
  ! calc_hyd_geo -- Calculate hydradulic conductivity by geometric mean
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: hyd_c1, hyd_c2, d1, d2
    real(DP), intent(out) :: hyd_geo
    ! -- local
    real(DP) :: c1, c2, d12
    !-------------------------------------------------------------------------------------------
    d12 = d1 + d2
    d12 = DONE/d12

    c1 = d1*log(hyd_c1)
    c2 = d2*log(hyd_c2)
    hyd_geo = exp((c1+c2)*d12)

  end subroutine calc_hyd_geo

  subroutine calc_hyd_arith(hyd_c1, hyd_c2, d1, d2, hyd_arith)
  !*********************************************************************************************
  ! calc_hyd_arith -- Calculate hydradulic conductivity by arithmetic mean
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: hyd_c1, hyd_c2, d1, d2
    real(DP), intent(out) :: hyd_arith
    ! -- local
    real(DP) :: c1, c2, d12
    !-------------------------------------------------------------------------------------------
    d12 = d1 + d2
    d12 = DONE/d12

    c1 = d1*hyd_c1
    c2 = d2*hyd_c2
    hyd_arith = (c1+c2)*d12

  end subroutine calc_hyd_arith

end module calc_parameter
