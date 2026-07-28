module calc_output
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE
  use types_module, only: sol_set
  use set_cell, only: ncalc, ncals, st_conn
  use set_condition, only: st_hydr, st_bcnd
  use prep_calculation, only: st_time
  use assign_boundary, only: st_forc
  use allocate_output, only: st_msloc

  implicit none
  private
  public :: calc_wtable, calc_cell_mas, calc_out_mass, calc_outvelc
  public :: calc_rivr_off, calc_lakr_off, calc_sufr_off, calc_dunr_off
  public :: calc_seal_res, calc_rech_res, calc_well_res

  ! -- local

  contains

  subroutine calc_wtable(hyd_head, deg_satu)
  !*********************************************************************************************
  ! calc_wtable -- Calculate water table
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_grid
    use set_cell, only: get_cals_grid
    use allocate_output, only: wtable
    ! -- inout
    real(DP), intent(in) :: hyd_head(:), deg_satu(:)
    ! -- local
    integer(I4) :: i, k, nijk, c_num
    integer(I4) :: i_num, j_num
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, k, i_num, j_num, c_num, nijk)
    do i = 1, ncals
      call get_cals_grid(i, i_num, j_num)
      unsat: do k = st_grid%nz, 1, -1
        nijk = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(j_num-1) + i_num
        c_num = 0 ; c_num = findloc(st_conn%loc2glo_ijk(:), value = nijk, dim = 1)
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

  subroutine calc_cell_mas(st_sol)
  !*********************************************************************************************
  ! calc_cell_mas -- Calculate cell massbalance
  !*********************************************************************************************
    ! -- modules
    use calc_function, only: calc_mass
    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: ms_st(:), ms_co(:), ms_se(:), ms_we(:)
    real(DP), allocatable :: ms_re(:), ms_su(:), ms_ri(:), ms_la(:)
    !-------------------------------------------------------------------------------------------
    allocate(ms_st(ncalc), ms_co(ncalc), ms_se(ncalc), ms_we(ncalc))
    allocate(ms_re(ncals), ms_su(ncals), ms_ri(ncals), ms_la(ncals))
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      ms_st(i) = DZERO ; ms_co(i) = DZERO ; ms_se(i) = DZERO ; ms_we(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncals
      ms_re(i) = DZERO ; ms_su(i) = DZERO ; ms_ri(i) = DZERO ; ms_la(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    ! -- Calculate massbalance (mass)
      call calc_mass(0, st_sol%head_old, st_sol%srat_old, st_sol%surf_old, st_sol%head_new,&
                     ms_st, ms_co, ms_se, ms_we, ms_re, ms_su, ms_ri, ms_la, st_sol%srat_new,&
                     st_sol%rel_perm, st_sol%surf_rati)
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      st_msloc%sto(i) = st_msloc%sto(i) + ms_st(i)
      st_msloc%con(i) = st_msloc%con(i) + ms_co(i)
      st_msloc%sea(i) = st_msloc%sea(i) + ms_se(i)
      st_msloc%wel(i) = st_msloc%wel(i) + ms_we(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncals
      st_msloc%rec(i) = st_msloc%rec(i) + ms_re(i)
      st_msloc%sur(i) = st_msloc%sur(i) + ms_su(i)
      st_msloc%riv(i) = st_msloc%riv(i) + ms_ri(i)
      st_msloc%lak(i) = st_msloc%lak(i) + ms_la(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(ms_re, ms_we, ms_st, ms_co, ms_su, ms_ri, ms_la)

  end subroutine calc_cell_mas

  subroutine calc_out_mass()
  !*********************************************************************************************
  ! calc_out_mass -- Calculate output massbalance
  !*********************************************************************************************
    ! -- modules
    use assign_calc, only: mass_num, msout_tnum, int_mass, mass2calc
    use allocate_output, only: st_msglo
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k
    real(DP), allocatable :: ms_sto(:), ms_con(:), ms_sea(:), ms_wel(:)
    real(DP), allocatable :: ms_rec(:), ms_sur(:), ms_riv(:), ms_lak(:)
    !-------------------------------------------------------------------------------------------
    allocate(ms_sto(msout_tnum), ms_con(msout_tnum), ms_sea(msout_tnum), ms_wel(msout_tnum))
    allocate(ms_rec(msout_tnum), ms_sur(msout_tnum), ms_riv(msout_tnum), ms_lak(msout_tnum))
    !$omp parallel
    !$omp do private(i)
    do i = 1, msout_tnum
      st_msglo%sto(i) = DZERO ; st_msglo%con(i) = DZERO
      st_msglo%sea(i) = DZERO ; st_msglo%wel(i) = DZERO
      st_msglo%rec(i) = DZERO ; st_msglo%sur(i) = DZERO
      st_msglo%riv(i) = DZERO ; st_msglo%lak(i) = DZERO ; st_msglo%tot(i) = DZERO
      ms_sto(i) = DZERO ; ms_con(i) = DZERO ; ms_sea(i) = DZERO ; ms_wel(i) = DZERO
      ms_rec(i) = DZERO ; ms_sur(i) = DZERO ; ms_riv(i) = DZERO ; ms_lak(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

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

    !$omp parallel
    !$omp do private(i)
    do i = 1, msout_tnum
      st_msglo%sto(i) = ms_sto(i) ; st_msglo%con(i) = ms_con(i)
      st_msglo%sea(i) = ms_sea(i) ; st_msglo%wel(i) = ms_wel(i)
      st_msglo%rec(i) = ms_rec(i) ; st_msglo%sur(i) = ms_sur(i)
      st_msglo%riv(i) = ms_riv(i) ; st_msglo%lak(i) = ms_lak(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, msout_tnum
      st_msglo%tot(i) = ms_sto(i) + ms_con(i) + ms_sea(i) + ms_wel(i) + ms_rec(i) +&
                        ms_sur(i) + ms_riv(i) + ms_lak(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncalc
      st_msloc%sto(i) = DZERO ; st_msloc%con(i) = DZERO
      st_msloc%sea(i) = DZERO ; st_msloc%wel(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncals
      st_msloc%rec(i) = DZERO ; st_msloc%sur(i) = DZERO
      st_msloc%riv(i) = DZERO ; st_msloc%lak(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    deallocate(ms_sto, ms_con, ms_sea, ms_wel, ms_rec, ms_sur, ms_riv, ms_lak)

  end subroutine calc_out_mass

  subroutine calc_outvelc(st_sol)
  !*********************************************************************************************
  ! calc_outvelc -- Calculate output velocity
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DHALF
    use make_cell, only: st_geom
    use calc_parameter, only: calc_hyd_upwind
    use allocate_solution, only: dir_conn, dir_seal, crs_index
    use allocate_output, only: pointv, facev
    ! -- inout

    type(sol_set), intent(in) :: st_sol
    ! -- local
    integer(I4) :: i, j, k, c, s
    integer(I4) :: sta_ind, end_ind, ind
    integer(I4) :: dir, d
    real(DP) :: delhead, relat, invdis, relp1, relp2
    !-------------------------------------------------------------------------------------------
    !$omp parallel
    !$omp do private(i, j, k, sta_ind, end_ind, ind, dir, delhead, relat, relp1, relp2)
    do i = 1, ncalc
      sta_ind = crs_index(1)%offind(i-1) ; end_ind = crs_index(1)%offind(i)
      do k = 1, end_ind-sta_ind
        ind = sta_ind + k ; j = crs_index(1)%offrow(ind)
        delhead = st_sol%head_new(j) - st_sol%head_new(i)
        relp1 = st_sol%rel_perm(i) ; relp2 = st_sol%rel_perm(j)

        ! -- Calculate hydradulic conductivity by upwind (hyd_upwind)
          call calc_hyd_upwind(-delhead, relp1, relp2, relat)

        dir = dir_conn(ind)
        facev(i,dir) = st_hydr%hydf_conn(ind)*relat*delhead*st_hydr%inv_dis(ind)
      end do
    end do
    !$omp end do

    !$omp do private(i, invdis, delhead)
    do i = 1, ncals
      invdis = DONE/(st_geom%surf_elev(i)-st_geom%cell_cent(i))
      delhead = st_sol%surf_head(i) - st_sol%head_new(i)
      facev(i,1) = st_hydr%hydf_surf(i)*delhead*invdis*st_sol%rel_perm(i)
    end do
    !$omp end do

    !$omp do private(i, c, s, dir, delhead)
    do i = 1, st_bcnd%seal_num
      c = st_bcnd%seal2calc(i) ; dir = dir_seal(i) ; s = st_bcnd%seal2seal(i)
      delhead = st_forc%read_seal(s) - st_sol%head_new(c)
      facev(c,dir) = st_hydr%hydf_seal(i)*delhead*st_hydr%dis_seal(i)*st_sol%rel_perm(c)
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

  subroutine calc_rivr_off(st_sol)
  !*********************************************************************************************
  ! calc_rivr_off -- Calculate river runoff
  !*********************************************************************************************
    ! -- modules
    use calc_function, only: func_riveterm
    use allocate_output, only: rive_sumtime, roff_rive
    ! -- inout

    type(sol_set), intent(in) :: st_sol
    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: rives(:), temp_rive(:)
    !-------------------------------------------------------------------------------------------
    allocate(rives(ncals), temp_rive(st_bcnd%rive_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      rives(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, st_bcnd%rive_num
      temp_rive(i) = roff_rive(i)
    end do
    !$omp end do
    !$omp end parallel

    ! -- Function river term (riveterm)
      call func_riveterm(st_sol%head_new, st_sol%rel_perm, rives)

    !$omp parallel do private(i, s)
    do i = 1, st_bcnd%rive_num
      s = st_bcnd%rive2cals(i)
      roff_rive(i) = temp_rive(i) - rives(s)*st_time%delt
    end do
    !$omp end parallel do

    deallocate(rives, temp_rive)

    rive_sumtime = rive_sumtime + st_time%delt

  end subroutine calc_rivr_off

  subroutine calc_lakr_off(st_sol)
  !*********************************************************************************************
  ! calc_lakr_off -- Calculate lake runoff
  !*********************************************************************************************
    ! -- modules
    use calc_function, only: func_laketerm
    use allocate_output, only: lake_sumtime, roff_lake
    ! -- inout

    type(sol_set), intent(in) :: st_sol
    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: lakes(:), temp_lake(:)
    !-------------------------------------------------------------------------------------------
    allocate(lakes(ncals), temp_lake(st_bcnd%lake_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      lakes(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, st_bcnd%lake_num
      temp_lake(i) = roff_lake(i)
    end do
    !$omp end do
    !$omp end parallel

    ! -- Function lake term (laketerm)
      call func_laketerm(st_sol%head_new, st_sol%rel_perm, lakes)

    !$omp parallel do private(i, s)
    do i = 1, st_bcnd%lake_num
      s = st_bcnd%lake2cals(i)
      roff_lake(i) = temp_lake(i) - lakes(s)*st_time%delt
    end do
    !$omp end parallel do

    deallocate(lakes, temp_lake)

    lake_sumtime = lake_sumtime + st_time%delt

  end subroutine calc_lakr_off

  subroutine calc_sufr_off(st_sol)
  !*********************************************************************************************
  ! calc_sufr_off -- Calculate surface runoff
  !*********************************************************************************************
    ! -- modules
    use calc_function, only: func_surfterm
    use allocate_output, only: surf_sumtime, roff_surf
    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: surfs(:), temp_surf(:)
    !-------------------------------------------------------------------------------------------
    allocate(surfs(ncals), temp_surf(ncals))
    !$omp parallel do private(i)
    do i = 1, ncals
      surfs(i) = DZERO
      temp_surf(i) = roff_surf(i)
    end do
    !$omp end parallel do

    ! -- Function surface term (surfterm)
      call func_surfterm(st_sol%head_new, st_sol%surf_old, st_sol%rel_perm, surfs,&
                         st_sol%surf_rati)

    !$omp parallel do private(i)
    do i = 1, ncals
      roff_surf(i) = temp_surf(i) - surfs(i)*st_time%delt
    end do
    !$omp end parallel do

    deallocate(surfs, temp_surf)

    surf_sumtime = surf_sumtime + st_time%delt

  end subroutine calc_sufr_off

  subroutine calc_dunr_off()
  !*********************************************************************************************
  ! calc_dunr_off -- Calculate dunne runoff
  !*********************************************************************************************
    ! -- modules
    use allocate_output, only: dunn_sumtime, roff_dunn
    ! -- inout

    ! -- local
    integer(I4) :: i, s
    real(DP), allocatable :: dunns(:), temp_dunn(:)
    !-------------------------------------------------------------------------------------------
    allocate(dunns(st_bcnd%rech_num), temp_dunn(st_bcnd%rech_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%rech_num
      dunns(i) = DZERO
      temp_dunn(i) = roff_dunn(i)
    end do
    !$omp end do

    !$omp do private(i, s)
    do i = 1, st_bcnd%rech_num
      s = st_bcnd%rech2cals(i)
      dunns(i) = st_forc%read_rech(i) - st_forc%calc_rech(i)/st_hydr%rech_area(s)
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, st_bcnd%rech_num
      roff_dunn(i) = temp_dunn(i) + dunns(i)*st_time%delt
    end do
    !$omp end do
    !$omp end parallel

    deallocate(dunns, temp_dunn)

    dunn_sumtime = dunn_sumtime + st_time%delt

  end subroutine calc_dunr_off

  subroutine calc_seal_res(st_sol)
  !*********************************************************************************************
  ! calc_seal_res -- Calculate sea level results
  !*********************************************************************************************
    ! -- modules
    use calc_function, only: func_sealterm
    use allocate_output, only: res_snum, res_seal
    ! -- inout

    type(sol_set), intent(in) :: st_sol
    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: sealr(:), temp_seal(:)
    !-------------------------------------------------------------------------------------------
    allocate(sealr(ncalc), temp_seal(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      sealr(i) = DZERO
      temp_seal(i) = res_seal(i)
    end do
    !$omp end parallel do

    ! -- Function sea level term (sealterm)
      call func_sealterm(st_sol%head_new, st_sol%rel_perm, sealr)

    !$omp parallel do private(i)
    do i = 1, ncalc
      res_seal(i) = temp_seal(i) + sealr(i)*st_time%delt
      res_snum(i) = i
    end do
    !$omp end parallel do

    deallocate(sealr, temp_seal)

  end subroutine calc_seal_res

  subroutine calc_rech_res()
  !*********************************************************************************************
  ! calc_rech_res -- Calculate recharge results
  !*********************************************************************************************
    ! -- modules
    use calc_function, only: func_rechterm
    use allocate_output, only: res_rnum, res_rech
    ! -- inout

    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: rechr(:), temp_rech(:)
    !-------------------------------------------------------------------------------------------
    allocate(rechr(ncals), temp_rech(ncals))
    !$omp parallel do private(i)
    do i = 1, ncals
      rechr(i) = DZERO
      temp_rech(i) = res_rech(i)
    end do
    !$omp end parallel do

    ! -- Function recharge term (rechterm)
      call func_rechterm(rechr)

    !$omp parallel do private(i)
    do i = 1, ncals
      res_rech(i) = temp_rech(i) + rechr(i)*st_time%delt
      res_rnum(i) = i
    end do
    !$omp end parallel do

    deallocate(rechr, temp_rech)

  end subroutine calc_rech_res

  subroutine calc_well_res()
  !*********************************************************************************************
  ! calc_well_res -- Calculate well results
  !*********************************************************************************************
    ! -- modules
    use calc_function, only: func_wellterm
    use allocate_output, only: res_wnum, res_well
    ! -- inout

    ! -- local
    integer(I4) :: i
    real(DP), allocatable :: wellr(:), temp_well(:)
    !-------------------------------------------------------------------------------------------
    allocate(wellr(ncalc), temp_well(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      wellr(i) = DZERO
      temp_well(i) = res_well(i)
    end do
    !$omp end parallel do

    ! -- Function well term (wellterm)
      call func_wellterm(wellr)

    !$omp parallel do private(i)
    do i = 1, ncalc
      res_well(i) = temp_well(i) + wellr(i)*st_time%delt
      res_wnum(i) = i
    end do
    !$omp end parallel do

    deallocate(wellr, temp_well)

  end subroutine calc_well_res

end module calc_output
