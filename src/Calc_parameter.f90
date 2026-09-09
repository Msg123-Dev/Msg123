module calc_parameter
  ! -- modules
  use kind_module, only: DP
  use constval_module, only: DZERO, DONE

  implicit none
  private
  public :: calc_srat_rperm, calc_hyd_geoharm, calc_hyd_upwind
  public :: calc_hyd_harm, calc_hyd_geo, calc_hyd_arith

  ! -- local

  contains

  subroutine calc_srat_rperm(num, pertur, pres, srat, rperm, stor, dkr_dpsi, dstor_dpsi)
  !*********************************************************************************************
  ! calc_srat_rperm -- Calculate saturation and relative permeability
  !*********************************************************************************************
    ! -- modules
    use kind_module, only: I4
    use initial_module, only: st_schm
    use make_cell, only: st_geom
    use set_condition, only: st_hydr
    ! -- inout
    integer(I4), intent(in) :: num
    real(DP), intent(in) :: pertur
    real(DP), intent(in) :: pres(:)
    real(DP), intent(out) :: rperm(:)
    real(DP), intent(out), optional :: srat(:)
    real(DP), intent(out), optional :: stor(:)
    real(DP), intent(out), optional :: dkr_dpsi(:)
    real(DP), intent(out), optional :: dstor_dpsi(:)
    ! -- local
    integer(I4) :: i
    real(DP) :: retm, beta, theta, kr, se, stor_val, srat_val
    real(DP) :: per_phead, phead0, kr_phead, dkr
    real(DP) :: pors_val, resi_val, dtheta_dpsi
    real(DP) :: kr_lin, kr_grad
    !-------------------------------------------------------------------------------------------
    phead0 = -st_schm%krlin_head

    !$omp parallel do private(i, pors_val, resi_val, retm, per_phead, kr_lin, kr_grad, dkr) &
    !$omp             private(beta, se, theta, srat_val, stor_val, dtheta_dpsi, kr_phead, kr)
    do i = 1, num
      pors_val = st_hydr%read_pors(i) ; resi_val = st_hydr%read_resi(i)
      retm = DONE - DONE/st_hydr%read_vann(i)
      per_phead = pres(i) - st_geom%cell_top(i) + pertur
      kr_lin = DONE ; kr_grad = DZERO
      ! -- Relative permeability and its slope at the linear bridge head
      if (st_schm%krlin_head > DZERO) then
        call calc_kr_vgm(phead0, st_hydr%read_vana(i), st_hydr%read_vann(i), retm, kr_lin, dkr)
        kr_grad = (DONE-kr_lin)/st_schm%krlin_head
      end if

      ! -- Storage and saturation ratio at the cell top
      if (per_phead < DZERO) then
        beta = abs(per_phead*st_hydr%read_vana(i))**st_hydr%read_vann(i)
        se = (DONE+beta)**(-retm)
        theta = resi_val + (pors_val-resi_val)*se
        srat_val = theta/pors_val
        if (st_schm%stor_type == 1) then
          stor_val = theta + st_hydr%read_spst(i)*per_phead*srat_val
        else
          stor_val = theta
        end if
        if (present(dstor_dpsi)) then
          dtheta_dpsi = (pors_val-resi_val)*se*retm*st_hydr%read_vann(i)*beta&
                        /(abs(per_phead)*(DONE+beta))
          if (st_schm%stor_type == 1) then
            dstor_dpsi(i) = dtheta_dpsi*(DONE+st_hydr%read_spst(i)*per_phead/pors_val)&
                            + st_hydr%read_spst(i)*srat_val
          else
            dstor_dpsi(i) = dtheta_dpsi
          end if
        end if
      else
        srat_val = DONE
        stor_val = pors_val + st_hydr%read_spst(i)*per_phead
        if (present(dstor_dpsi)) then
          dstor_dpsi(i) = st_hydr%read_spst(i)
        end if
      end if

      ! -- Relative permeability at the head position given by krpos_type
      if (st_schm%krpos_type == 1) then
        kr_phead = pres(i) - st_geom%cell_cent(i) + pertur
      else
        kr_phead = per_phead
      end if
      if (kr_phead >= DZERO) then
        kr = DONE ; dkr = DZERO
      else if (kr_phead > phead0) then
        kr = kr_lin + (kr_phead-phead0)*kr_grad ; dkr = kr_grad
      else
        call calc_kr_vgm(kr_phead, st_hydr%read_vana(i), st_hydr%read_vann(i), retm, kr, dkr)
      end if

      rperm(i) = kr
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

  subroutine calc_kr_vgm(phead, vana, vann, retm, kr, dkr)
  !*********************************************************************************************
  ! calc_kr_vgm -- Calculate relative permeability by van Genuchten and Mualem
  !*********************************************************************************************
    ! -- modules
    use kind_module, only: SP
    use constval_module, only: DHALF, DTWO
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
    if (head_dif >= DZERO) then
      hyd_upwind = hyd_c1
    else
      hyd_upwind = hyd_c2
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
