module calc_parameter
  ! -- modules
  use kind_module, only: DP
  use constval_module, only: DONE, DZERO

  implicit none
  private
  public :: calc_srat_rperm, calc_hyd_geoharm, calc_hyd_upwind
  public :: calc_hyd_harm, calc_hyd_geo, calc_hyd_arith

  ! -- local

  contains

  subroutine calc_srat_rperm(num, pertur, pres, srat, rperm, sstor)
  !***************************************************************************************
  ! calc_srat_rperm -- Calculate saturation and relative permeability
  !***************************************************************************************
    ! -- modules
    use kind_module, only: I4
    use constval_module, only: DTWO, DHALF
    use make_cell, only: cell_top
    use assign_calc, only: read_retn, read_reta, read_poro, read_resi, read_ss
    ! -- inout
    integer(I4), intent(in) :: num
    real(DP), intent(in) :: pertur
    real(DP), intent(in) :: pres(:)
    real(DP), intent(out) :: srat(:), rperm(:)
    real(DP), intent(out), optional :: sstor(:)
    ! -- local
    integer(I4) :: i
    real(DP) :: retm, beta, theta, kr, se, ss
    real(DP) :: per_phead, phead0
    real(DP), allocatable :: temp_ss(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_ss(num))
    !$omp parallel workshare
    temp_ss(:) = DZERO
    !$omp end parallel workshare

    phead0 = DZERO

    !$omp parallel do private(i, per_phead, retm, beta, theta, kr, se, ss)
    do i = 1, num
      per_phead = pres(i) - cell_top(i) + pertur
      if (per_phead < DZERO) then
        retm = DONE - DONE/read_retn(i)
        if (per_phead <= phead0) then
          beta = abs(per_phead*read_reta(i))**read_retn(i)
          theta = read_resi(i) + (read_poro(i)-read_resi(i))*((DONE+beta))**(-retm)
        else
          beta = DZERO
          theta = read_poro(i)
        end if
        srat(i) = theta/read_poro(i)
        se = (DONE+beta)**(-retm)
        kr = se**(DHALF)*(DONE-(DONE-se**(DONE/retm))**retm)**DTWO
        ss = DZERO
      else
        srat(i) = DONE
        kr = DONE
        ss = read_ss(i)
      end if
      rperm(i) = kr
      temp_ss(i) = ss
    end do
    !$omp end parallel do

    if (present(sstor)) then
      !$omp parallel workshare
      sstor(:) = temp_ss(:)
      !$omp end parallel workshare
    end if

    deallocate(temp_ss)

  end subroutine calc_srat_rperm

  subroutine calc_hyd_geoharm(srato1, srato2, hyd_c1, hyd_c2, d1, d2, hyd_gh)
  !***************************************************************************************
  ! calc_hyd_geoharm -- Calculate hydradulic conductivity by geometric and harmonic
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: srato1, srato2, hyd_c1, hyd_c2, d1, d2
    real(DP), intent(out) :: hyd_gh
    ! -- local
    real(DP) :: c1, c2, d12, cd
    !-------------------------------------------------------------------------------------
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
  !***************************************************************************************
  ! calc_hyd_upwind -- Calculate hydradulic conductivity by upwind
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: head_dif, hyd_c1, hyd_c2
    real(DP), intent(out) :: hyd_upwind
    ! -- local

    !-------------------------------------------------------------------------------------
    if (head_dif >= DZERO) then
      hyd_upwind = hyd_c1
    else
      hyd_upwind = hyd_c2
    end if

  end subroutine calc_hyd_upwind

  subroutine calc_hyd_harm(hyd_c1, hyd_c2, d1, d2, hyd_harm)
  !***************************************************************************************
  ! calc_hyd_harm -- Calculate hydradulic conductivity by harmonic mean
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: hyd_c1, hyd_c2, d1, d2
    real(DP), intent(out) :: hyd_harm
    ! -- local
    real(DP) :: c1, c2, d12, cd, cd12
    !-------------------------------------------------------------------------------------
    d12 = d1 + d2

    c1 = hyd_c1*d12
    c2 = hyd_c2*d12
    cd12 = d1*c2+d2*c1
    cd = DONE/cd12
    hyd_harm = c1*c2*cd

  end subroutine calc_hyd_harm

  subroutine calc_hyd_geo(hyd_c1, hyd_c2, d1, d2, hyd_geo)
  !***************************************************************************************
  ! calc_hyd_geo -- Calculate hydradulic conductivity by geometric mean
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: hyd_c1, hyd_c2, d1, d2
    real(DP), intent(out) :: hyd_geo
    ! -- local
    real(DP) :: c1, c2, d12
    !-------------------------------------------------------------------------------------
    d12 = d1 + d2
    d12 = DONE/d12

    c1 = d1*log(hyd_c1)
    c2 = d2*log(hyd_c2)
    hyd_geo = exp((c1+c2)*d12)

  end subroutine calc_hyd_geo

  subroutine calc_hyd_arith(hyd_c1, hyd_c2, d1, d2, hyd_arith)
  !***************************************************************************************
  ! calc_hyd_arith -- Calculate hydradulic conductivity by arithmetic mean
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: hyd_c1, hyd_c2, d1, d2
    real(DP), intent(out) :: hyd_arith
    ! -- local
    real(DP) :: c1, c2, d12
    !-------------------------------------------------------------------------------------
    d12 = d1 + d2
    d12 = DONE/d12

    c1 = d1*hyd_c1
    c2 = d2*hyd_c2
    hyd_arith = (c1+c2)*d12

  end subroutine calc_hyd_arith

end module calc_parameter
