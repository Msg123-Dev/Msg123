module calc_output
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE
  use set_cell, only: ncalc, ncals
  use prep_calculation, only: delt
  use allocate_solution, only: head_new
  use allocate_output, only: st_msloc

  implicit none
  private
  public :: calc_wtable, calc_cell_mas, calc_out_mass, calc_outvelc
  public :: calc_rivr_off, calc_lakr_off, calc_sufr_off, calc_dunr_off
  public :: calc_seal_res, calc_rech_res, calc_well_res

  ! -- local

  contains

  subroutine calc_wtable(hyd_head, deg_satu)
  !***************************************************************************************
  ! calc_wtable -- Calculate water table
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_grid
    use set_cell, only: loc2glo_ijk, get_cals_grid
    use allocate_output, only: wtable
    ! -- inout
    real(DP), intent(in) :: hyd_head(:), deg_satu(:)
    ! -- local
    integer(I4) :: i, k, nijk, c_num
    integer(I4) :: i_num, j_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, k, i_num, j_num, c_num, nijk)
    do i = 1, ncals
      call get_cals_grid(i, i_num, j_num)
      unsat: do k = st_grid%nz, 1, -1
        nijk = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
        c_num = 0 ; c_num = findloc(loc2glo_ijk(:), value = nijk, dim = 1)
        if (c_num /= 0 .and. c_num <= ncalc) then
          if (deg_satu(c_num) /= DONE) then
            wtable(i) = hyd_head(c_num)
            exit unsat
          else
            wtable(i) = hyd_head(c_num)
          end if
        end if
      end do unsat
    end do
    !$omp end parallel do

  end subroutine calc_wtable

  subroutine calc_cell_mas()
  !***************************************************************************************
  ! calc_cell_mas -- Calculate cell massbalance
  !***************************************************************************************
    ! -- modules
    use calc_function, only: calc_mass
    use allocate_solution, only: head_old
    ! -- inout

    ! -- local
    real(DP), allocatable :: ms_st(:), ms_co(:), ms_se(:), ms_we(:)
    real(DP), allocatable :: ms_re(:), ms_su(:), ms_ri(:), ms_la(:)
    !-------------------------------------------------------------------------------------
    allocate(ms_st(ncalc), ms_co(ncalc), ms_se(ncalc), ms_we(ncalc))
    allocate(ms_re(ncals), ms_su(ncals), ms_ri(ncals), ms_la(ncals))
    !$omp parallel workshare
    ms_st(:) = DZERO ; ms_co(:) = DZERO ; ms_se(:) = DZERO ; ms_we(:) = DZERO
    ms_re(:) = DZERO ; ms_su(:) = DZERO ; ms_ri(:) = DZERO ; ms_la(:) = DZERO
    !$omp end parallel workshare

    ! -- Calculate massbalance (mass)
      call calc_mass(0, head_new, head_old, ms_st, ms_co, ms_se, ms_we, ms_re, ms_su,&
                     ms_ri, ms_la)

    !$omp parallel workshare
    st_msloc%sto(:) = st_msloc%sto(:) + ms_st(:)
    st_msloc%con(:) = st_msloc%con(:) + ms_co(:)
    st_msloc%sea(:) = st_msloc%sea(:) + ms_se(:)
    st_msloc%wel(:) = st_msloc%wel(:) + ms_we(:)
    st_msloc%rec(:) = st_msloc%rec(:) + ms_re(:)
    st_msloc%sur(:) = st_msloc%sur(:) + ms_su(:)
    st_msloc%riv(:) = st_msloc%riv(:) + ms_ri(:)
    st_msloc%lak(:) = st_msloc%lak(:) + ms_la(:)
    !$omp end parallel workshare

    deallocate(ms_re, ms_we, ms_st, ms_co, ms_su, ms_ri, ms_la)

  end subroutine calc_cell_mas

  subroutine calc_out_mass()
  !***************************************************************************************
  ! calc_out_mass -- Calculate output massbalance
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_in_type, in_type
    use assign_calc, only: mass_num, msout_tnum, mass2calc, int_mass
    use allocate_output, only: st_msglo
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k
    real(DP), allocatable :: ms_sto(:), ms_con(:), ms_sea(:), ms_wel(:)
    real(DP), allocatable :: ms_rec(:), ms_sur(:), ms_riv(:), ms_lak(:)
    !-------------------------------------------------------------------------------------
    allocate(ms_sto(msout_tnum), ms_con(msout_tnum), ms_sea(msout_tnum), ms_wel(msout_tnum))
    allocate(ms_rec(msout_tnum), ms_sur(msout_tnum), ms_riv(msout_tnum), ms_lak(msout_tnum))
    !$omp parallel
    !$omp workshare
    st_msglo%sto(:) = DZERO ; st_msglo%con(:) = DZERO ; st_msglo%sea(:) = DZERO
    st_msglo%wel(:) = DZERO ; st_msglo%rec(:) = DZERO ; st_msglo%sur(:) = DZERO
    st_msglo%riv(:) = DZERO ; st_msglo%lak(:) = DZERO ; st_msglo%tot(:) = DZERO
    ms_sto(:) = DZERO ; ms_con(:) = DZERO ; ms_sea(:) = DZERO ; ms_wel(:) = DZERO
    ms_rec(:) = DZERO ; ms_sur(:) = DZERO ; ms_riv(:) = DZERO ; ms_lak(:) = DZERO
    !$omp end workshare

    if (st_in_type%mass /= in_type(7)) then
      !$omp do private(i, j, k)
      do i = 1, mass_num
        j = mass2calc(i) ; k = int_mass(i)
        if (j <= ncals) then
          ms_rec(k) = ms_rec(k) + st_msloc%rec(j)
          ms_sur(k) = ms_sur(k) + st_msloc%sur(j)
          ms_riv(k) = ms_riv(k) + st_msloc%riv(j)
          ms_lak(k) = ms_lak(k) + st_msloc%lak(j)
        end if
        ms_sto(k) = ms_sto(k) + st_msloc%sto(j)
        ms_con(k) = ms_con(k) + st_msloc%con(j)
        ms_sea(k) = ms_sea(k) + st_msloc%sea(j)
        ms_wel(k) = ms_wel(k) + st_msloc%wel(j)
      end do
      !$omp end do

    else if (st_in_type%mass == in_type(7)) then
      !$omp workshare
      ms_sto(:) = sum(st_msloc%sto) ; ms_con(:) = sum(st_msloc%con)
      ms_sea(:) = sum(st_msloc%sea) ; ms_wel(:) = sum(st_msloc%wel)
      ms_rec(:) = sum(st_msloc%rec) ; ms_sur(:) = sum(st_msloc%sur)
      ms_riv(:) = sum(st_msloc%riv) ; ms_lak(:) = sum(st_msloc%lak)
      !$omp end workshare
    end if

    !$omp workshare
    st_msglo%sto(:) = ms_sto(:) ; st_msglo%con(:) = ms_con(:)
    st_msglo%sea(:) = ms_sea(:) ; st_msglo%wel(:) = ms_wel(:)
    st_msglo%rec(:) = ms_rec(:) ; st_msglo%sur(:) = ms_sur(:)
    st_msglo%riv(:) = ms_riv(:) ; st_msglo%lak(:) = ms_lak(:)
    !$omp end workshare

    !$omp workshare
    st_msglo%tot(:) = ms_sto(:) + ms_con(:) + ms_sea(:) + ms_wel(:) + ms_rec(:) +&
                      ms_sur(:) + ms_riv(:) + ms_lak(:)
    !$omp end workshare

    !$omp workshare
    st_msloc%sto(:) = DZERO ; st_msloc%con(:) = DZERO ; st_msloc%sea(:) = DZERO
    st_msloc%wel(:) = DZERO ; st_msloc%rec(:) = DZERO ; st_msloc%sur(:) = DZERO
    st_msloc%riv(:) = DZERO ; st_msloc%lak(:) = DZERO
    !$omp end workshare
    !$omp end parallel

    deallocate(ms_sto, ms_con, ms_sea, ms_wel, ms_rec, ms_sur, ms_riv, ms_lak)

  end subroutine calc_out_mass

  subroutine calc_outvelc()
  !***************************************************************************************
  ! calc_outvelc -- Calculate output velocity
  !***************************************************************************************
    ! -- modules
    use constval_module, only: DHALF
    use make_cell, only: surf_elev, cell_cent
    use set_condition, only: nseal, hydf_surf
    use assign_boundary, only: read_seal
    use calc_parameter, only: calc_hyd_upwind
    use allocate_solution, only: crs_index, hydf_conn, inv_dis, surf_head, rel_perm,&
                                 dir_conn, hydf_seal, dis_seal, seal2calc, seal2seal,&
                                 dir_seal
    use allocate_output, only: facev, pointv
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, c, s
    integer(I4) :: sta_ind, end_ind, ind
    integer(I4) :: dir, d
    real(DP) :: delhead, relat, invdis, relp1, relp2
    !-------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i, j, k, sta_ind, end_ind, ind, dir, delhead, relat, relp1, relp2)
    do i = 1, ncalc
      sta_ind = crs_index(1)%offind(i-1) ; end_ind = crs_index(1)%offind(i)
      do k = 1, end_ind-sta_ind
        ind = sta_ind + k ; j = crs_index(1)%offrow(ind)
        delhead = head_new(j) - head_new(i)
        relp1 = rel_perm(i) ; relp2 = rel_perm(j)

        ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
          call calc_hyd_upwind(-delhead, relp1, relp2, relat)

        dir = dir_conn(ind) ; facev(i,dir) = hydf_conn(ind)*relat*delhead*inv_dis(ind)
      end do
    end do
    !$omp end do

    !$omp do private(i, invdis, delhead)
    do i = 1, ncals
      invdis = DONE/(surf_elev(i)-cell_cent(i))
      delhead = surf_head(i) - head_new(i)
      facev(i,1) = hydf_surf(i)*delhead*invdis*rel_perm(i)
    end do
    !$omp end do

    !$omp do private(i, c, s, dir, delhead)
    do i = 1, nseal
      c = seal2calc(i) ; dir = dir_seal(i) ; s = seal2seal(i)
      delhead = read_seal(s) - head_new(c)
      facev(c,dir) = hydf_seal(i)*delhead*dis_seal(i)*rel_perm(c)
    end do
    !$omp end do

    !$omp do private(i, d)
    do i = 1, ncalc
      do d = 1, 3
        if (facev(i,d) > DZERO .and. facev(i,7-d) < DZERO) then
          pointv(i,4-d) = (facev(i,d) - facev(i,7-d))*DHALF
        else if (facev(i,d) < DZERO .and. facev(i,7-d) > DZERO) then
          pointv(i,4-d) = (facev(i,d) - facev(i,7-d))*DHALF
        else
          pointv(i,4-d) = facev(i,d) - facev(i,7-d)
        end if
        if (d == 1) then
          pointv(i,4-d) = -pointv(i,4-d)
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine calc_outvelc

  subroutine calc_rivr_off()
  !***************************************************************************************
  ! calc_rivr_off -- Calculate river runoff
  !***************************************************************************************
    ! -- modules
    use calc_boundary, only: rive2cals
    use set_boundary, only: rive_num
    use calc_function, only: func_riveterm
    use allocate_output, only: roff_rive, rive_sumtime
    ! -- inout

    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: rives(:), temp_rive(:)
    !-------------------------------------------------------------------------------------
    allocate(rives(ncals), temp_rive(rive_num))
    !$omp parallel workshare
    rives(:) = DZERO ; temp_rive(:) = roff_rive(:)
    !$omp end parallel workshare

    ! -- Function river term (riveterm)
      call func_riveterm(head_new, rives)

    !$omp parallel do private(i, s)
    do i = 1, rive_num
      s = rive2cals(i)
      roff_rive(i) = temp_rive(i) - rives(s)*delt
    end do
    !$omp end parallel do

    deallocate(rives, temp_rive)

    rive_sumtime = rive_sumtime + delt

  end subroutine calc_rivr_off

  subroutine calc_lakr_off()
  !***************************************************************************************
  ! calc_lakr_off -- Calculate lake runoff
  !***************************************************************************************
    ! -- modules
    use calc_boundary, only: lake2cals
    use set_boundary, only: lake_num
    use calc_function, only: func_laketerm
    use allocate_output, only: roff_lake, lake_sumtime
    ! -- inout

    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: lakes(:), temp_lake(:)
    !-------------------------------------------------------------------------------------
    allocate(lakes(ncals), temp_lake(lake_num))
    !$omp parallel workshare
    lakes(:) = DZERO ; temp_lake(:) = roff_lake(:)
    !$omp end parallel workshare

    ! -- Function lake term (laketerm)
      call func_laketerm(head_new, lakes)

    !$omp parallel do private(i, s)
    do i = 1, lake_num
      s = lake2cals(i)
      roff_lake(i) = temp_lake(i) - lakes(s)*delt
    end do
    !$omp end parallel do

    deallocate(lakes, temp_lake)

    lake_sumtime = lake_sumtime + delt

  end subroutine calc_lakr_off

  subroutine calc_sufr_off()
  !***************************************************************************************
  ! calc_sufr_off -- Calculate surface runoff
  !***************************************************************************************
    ! -- modules
    use calc_function, only: func_surfterm
    use allocate_solution, only: surf_old
    use allocate_output, only: roff_surf, surf_sumtime
    ! -- inout

    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: surfs(:), temp_surf(:)
    !-------------------------------------------------------------------------------------
    allocate(surfs(ncals), temp_surf(ncals))
    !$omp parallel workshare
    surfs(:) = DZERO ; temp_surf(:) = roff_surf(:)
    !$omp end parallel workshare

    ! -- Function surface term (surfterm)
      call func_surfterm(head_new, surf_old, surfs)

    !$omp parallel do private(i)
    do i = 1, ncals
      roff_surf(i) = temp_surf(i) - surfs(i)*delt
    end do
    !$omp end parallel do

    deallocate(surfs, temp_surf)

    surf_sumtime = surf_sumtime + delt

  end subroutine calc_sufr_off

  subroutine calc_dunr_off()
  !***************************************************************************************
  ! calc_dunr_off -- Calculate dunne runoff
  !***************************************************************************************
    ! -- modules
    use set_condition, only: rech_area
    use assign_boundary, only: read_rech
    use calc_boundary, only: calc_rech, rech2cals
    use set_boundary, only: rech_num
    use allocate_output, only: roff_dunn, dunn_sumtime
    ! -- inout

    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: dunns(:), temp_dunn(:)
    !-------------------------------------------------------------------------------------
    allocate(dunns(rech_num), temp_dunn(rech_num))
    !$omp parallel
    !$omp workshare
    dunns(:) = DZERO ; temp_dunn(:) = roff_dunn(:)
    !$omp end workshare

    !$omp do private(i, s)
    do i = 1, rech_num
      s = rech2cals(i)
      dunns(i) = read_rech(i) - calc_rech(i)/rech_area(s)
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, rech_num
      roff_dunn(i) = temp_dunn(i) + dunns(i)*delt
    end do
    !$omp end do
    !$omp end parallel

    deallocate(dunns, temp_dunn)

    dunn_sumtime = dunn_sumtime + delt

  end subroutine calc_dunr_off

  subroutine calc_seal_res()
  !***************************************************************************************
  ! calc_seal_res -- Calculate sea level results
  !***************************************************************************************
    ! -- modules
    use calc_function, only: func_sealterm
    use allocate_output, only: res_seal, res_snum
    ! -- inout

    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: sealr(:), temp_seal(:)
    !-------------------------------------------------------------------------------------
    allocate(sealr(ncalc), temp_seal(ncalc))
    !$omp parallel workshare
    sealr(:) = DZERO ; temp_seal(:) = res_seal(:)
    !$omp end parallel workshare

    ! -- Function sea level term (sealterm)
      call func_sealterm(head_new, sealr)

    !$omp parallel do private(i)
    do i = 1, ncalc
      res_seal(i) = temp_seal(i) + sealr(i)*delt
      res_snum(i) = i
    end do
    !$omp end parallel do

    deallocate(sealr, temp_seal)

  end subroutine calc_seal_res

  subroutine calc_rech_res()
  !***************************************************************************************
  ! calc_rech_res -- Calculate recharge results
  !***************************************************************************************
    ! -- modules
    use calc_function, only: func_rechterm
    use allocate_output, only: res_rech, res_rnum
    ! -- inout

    ! -- local
    integer(I4) :: s
    real(DP), allocatable :: rechr(:), temp_rech(:)
    !-------------------------------------------------------------------------------------
    allocate(rechr(ncals), temp_rech(ncals))
    !$omp parallel workshare
    rechr(:) = DZERO ; temp_rech(:) = res_rech(:)
    !$omp end parallel workshare

    ! -- Function recharge term (rechterm)
      call func_rechterm(rechr)

    !$omp parallel do private(s)
    do s = 1, ncals
      res_rech(s) = temp_rech(s) + rechr(s)*delt
      res_rnum(s) = s
    end do
    !$omp end parallel do

    deallocate(rechr, temp_rech)

  end subroutine calc_rech_res

  subroutine calc_well_res()
  !***************************************************************************************
  ! calc_well_res -- Calculate well results
  !***************************************************************************************
    ! -- modules
    use calc_function, only: func_wellterm
    use allocate_output, only: res_well, res_wnum
    ! -- inout

    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: wellr(:), temp_well(:)
    !-------------------------------------------------------------------------------------
    allocate(wellr(ncalc), temp_well(ncalc))
    !$omp parallel workshare
    wellr(:) = DZERO ; temp_well(:) = res_well(:)
    !$omp end parallel workshare

    ! -- Function well term (wellterm)
      call func_wellterm(wellr)

    !$omp parallel do private(i)
    do i = 1, ncalc
      res_well(i) = temp_well(i) + wellr(i)*delt
      res_wnum(i) = i
    end do
    !$omp end parallel do

    deallocate(wellr, temp_well)

  end subroutine calc_well_res

end module calc_output
