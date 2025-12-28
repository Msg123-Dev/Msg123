module time_module
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SZERO, DZERO, DONE, SNOVAL, DHALF
  use initial_module, only: st_sim, st_in_type, my_rank, st_rivf_type, st_lakf_type,&
                            st_step_flag, st_seal, st_rech, st_well, st_prec, st_evap,&
                            st_riwl, st_riwd, st_ribl, st_ride, st_riwi, st_rile,&
                            st_lawl, st_lawd, st_labl, st_laar
  use open_file, only: st_intre, st_intpr, st_intev
  use check_condition, only: st_out_fnum
  use set_cell, only: ncalc, ncals
  use prep_calculation, only: current_t, delt, conv_flag
  use set_condition, only: rech_area, well_index, well_conn, abyd_well
  use assign_boundary, only: read_rech, read_well, read_prec, read_evap, well_top,&
                             well_bott, calc_well
  use calc_boundary, only: rech2cals, calc_rech
  use set_boundary, only: rech_num, well_num, cflag_riv, cflag_lak, criv, clak, prec_num,&
                          evap_num
  use allocate_solution, only: head_new, surf_head, rel_perm, head_old

  implicit none
  private
  public :: update_tstep
  real(SP), public :: now_time

  ! -- local
  integer(I4) :: timestep_num, boundstep
  real(SP) :: next_time
  real(DP) :: delt_old1, delt_old2
  real(DP), allocatable :: whead_new(:), head_old2(:)

  contains

  subroutine update_tstep()
  !***************************************************************************************
  ! update_tstep -- Update time step
  !***************************************************************************************
    ! -- module
    use constval_module, only: DSMAL
    use utility_module, only: write_err_stop
    use assign_calc, only: geog_num
    use prep_calculation, only: delt_inv
#ifdef ICI
    use ici_module, only: get_var
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------
    ! -- Calculate next time step (nextst)
      call calc_nextst()

    if (conv_flag == 1 .and. st_sim%sim_type == 1) then
      ! -- Set next end time (nextet)
        call set_nextet()
      ! -- Set delta time (delt)
        call set_delt()

      if (boundstep >= 1) then
#ifdef ICI
        ! -- Get variables (var)
          call get_var()
#endif
      end if
      ! -- Set next variable (nextvar)
        call set_nextvar()
    else if (next_time > st_sim%end_time) then
      delt = real(st_sim%end_time - current_t, kind=DP)
    end if

    if (delt < DSMAL .and. st_sim%sim_type /= -1 .and. my_rank == 0) then
      call write_err_stop("Time Step is too small.")
    else if (st_sim%sim_type /= -1) then
      current_t = current_t + real(delt, kind=SP)
    else if (st_sim%sim_type == -1) then
      current_t = SZERO
    end if

    now_time = current_t/st_sim%cal_fact
    delt_inv = DONE/delt

    if (conv_flag == 1 .or. timestep_num == 0) then
      if (well_num /= 0) then
        ! -- Set virtual well head (vwell_head)
          call set_vwell_head()
      end if
      if (rech_num /= 0) then
        !$omp parallel do private(i, s)
        do i = 1, rech_num
          s = rech2cals(i)
          calc_rech(i) = read_rech(i)*rech_area(s)
        end do
        !$omp end parallel do
        if (geog_num /= 0) then
          ! -- Change the recharge volume
            call change_recharge()
        end if
      end if
    end if

!    if (rech_num /= 0) then
!      ! -- Change the volume rate using time step and cell volume (volrate)
!        call change_volrate()
!    end if

!    if (timestep_num /= 0 .and. st_sim%sim_type /= -1) then
!      ! -- Set estimate head (est_head)
!        call set_est_head(head_old2, head_old, head_new)
!    end if

  end subroutine update_tstep

  subroutine calc_nextst()
  !***************************************************************************************
  ! calc_nextst -- Calculate next time step
  !***************************************************************************************
    ! -- modules
    use constval_module, only: DNOVAL
    use utility_module, only: conv_unit
    use initial_module, only: maxinn_iter, st_init
    use read_input, only: temp_maxinn_iter, len_scal
    use prep_calculation, only: now_date, out_iter, inter_time
    use set_boundary, only: read_head
    use calc_parameter, only: calc_srat_rperm
    use allocate_solution, only: srat_old, srat_new, surf_old
#ifdef MPI_MSG
    use mpi_utility, only: bcast_val
    use initial_module, only: pro_totn
    use mpi_write, only: write_mpi_3dbin
#else
    use write_module, only: write_header_bin, write_3dbin
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: calc2calc(:)
    real(SP), save :: resi_time
    real(DP), allocatable :: cell_srat(:)
    !-------------------------------------------------------------------------------------
    if (current_t == DZERO) then
      timestep_num = 0 ; boundstep = 0
      delt_old1 = DZERO ; resi_time = SZERO
      if (st_sim%sim_type == -1) then
        delt = DNOVAL
        maxinn_iter = 100
      else
        allocate(head_old2(ncalc))
        !$omp parallel workshare
        head_old2(:) = read_head(:)
        !$omp end parallel workshare
        delt = st_sim%ini_step
        maxinn_iter = temp_maxinn_iter
      end if

      allocate(cell_srat(ncalc), calc2calc(ncalc))
      !$omp parallel workshare
      cell_srat(:) = DZERO ; calc2calc(:) = 0
      !$omp end parallel workshare
      ! -- Calculate saturation and relative permeability (srat_rperm)
        call calc_srat_rperm(ncalc, DZERO, read_head, cell_srat, rel_perm)
      !$omp parallel
      !$omp do private(i)
      do i = 1, ncalc
        head_old(i) = read_head(i) ; head_new(i) = read_head(i)
        srat_old(i) = cell_srat(i) ; srat_new(i) = cell_srat(i)
        calc2calc(i) = i
      end do
      !$omp end do

      !$omp do private(i)
      do i = 1, ncals
        surf_head(i) = read_head(i)
      end do
      !$omp end do
      !$omp end parallel

      deallocate(cell_srat)

      if (st_sim%sim_type /= -1) then
        deallocate(read_head)
      end if

      !$omp parallel do private(i)
      do i = 1, ncals
        surf_old(i) = surf_head(i)
      end do
      !$omp end parallel do

      if (st_sim%res_type == 1) then
        current_t = st_init%rest_time
      end if

      now_time = current_t/st_sim%cal_fact

#ifdef MPI_MSG
      ! -- Write MPI 3D binary file (mpi_3dbin)
        call write_mpi_3dbin(st_out_fnum%head, ncalc, calc2calc, len_scal, head_new, now_time)
#else
      ! -- Write header binary file (header_bin)
        call write_header_bin(st_out_fnum%head, now_time)
      ! -- Write 3D binary file (3dbin)
        call write_3dbin(st_out_fnum%head, ncalc, calc2calc, len_scal, head_new)
#endif
      deallocate(calc2calc)

    else if (st_sim%sim_type >= 0) then
      if (conv_flag == 1) then
        if (timestep_num == 0) then
          delt_old1 = delt
        end if
        timestep_num = timestep_num + 1
        delt_old2 = delt_old1
        delt_old1 = delt

        if (my_rank == 0) then
          ! -- Write boundary change information
            call write_bound_change(boundstep)
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(boundstep, "boundary change number")
        end if
#endif
        if (boundstep >= 1) then
          inter_time = current_t - resi_time
          ! -- Set next date (date)
            call set_date(inter_time, now_date, resi_time)
          resi_time = current_t + resi_time
          if (trim(adjustl(st_sim%cal_unit)) == "YEA") then
            call conv_unit(my_rank, st_sim%cal_unit, "main file", now_date, st_sim%cal_fact)
          end if
!          delt = delt*DHALF
          delt = delt
        else if (st_sim%sim_type == 0) then
          delt = delt*st_sim%inc_fact
!          ! -- Apply heuristic time stepping (heuri)
!            call apply_heuri(out_iter)
        else if (st_sim%sim_type == 1) then
          ! -- Apply heuristic time stepping (heuri)
            call apply_heuri(out_iter)
!          ! -- Apply adaptive time stepping (adapt)
!            call apply_adapt(head_old, head_old2, head_new)
        end if

        ! -- Set value exchange (valexc)
          call set_valexc(ncalc, head_old, head_old2)
          call set_valexc(ncalc, head_new, head_old)
          call set_valexc(ncalc, srat_new, srat_old)
          call set_valexc(ncals, surf_head, surf_old)
      else
        current_t = current_t - real(delt, kind=SP)
        delt = delt*DHALF
        ! -- Set value exchange (valexc)
          call set_valexc(ncalc, head_old, head_new)
          call set_valexc(ncalc, srat_old, srat_new)
          call set_valexc(ncals, surf_old, surf_head)
        if (st_sim%sim_type == 1) then
          ! -- Reset stepflag (stepf)
            call reset_stepf()
        end if
      end if
    else if (st_sim%sim_type == -1) then
      ! -- Set value exchange (valexc)
        call set_valexc(ncalc, head_new, head_old)
        call set_valexc(ncalc, srat_new, srat_old)
    end if

    next_time = current_t + real(delt, kind=SP)

  end subroutine calc_nextst

  subroutine set_nextet()
  !***************************************************************************************
  ! set_nextet -- Set next end time
  !***************************************************************************************
    ! -- modules
    use utility_module, only: write_logf
    use initial_module, only: in_type
    use open_file, only: st_intse, st_intwe
    use read_module, only: read_next, read_intn
    use assign_boundary, only: read_seal
#ifdef MPI_MSG
    use mpi_read, only: set_real4_fview
    use set_boundary, only: bfview, rfview, lfview
#endif
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    character(:), allocatable :: bound_name, err_mes, str_time
    !-------------------------------------------------------------------------------------
    ierr = 0 ; str_time = ""
    if (st_step_flag%seal == 1) then
      bound_name = "sea level"
      if (st_in_type%seal /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_in_type%seal, st_seal%fnum, st_seal%multi, bound_name,&
                         st_seal%totn, st_step_flag%seal, ierr, st_seal%etime)
      else if (st_in_type%seal == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_intse%type, st_intse%fnum, st_seal%multi, st_intse%step,&
                         bound_name, st_seal%fnum, st_step_flag%seal, ierr, st_seal%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_seal%fnum, bfview%seal, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%seal = 0
        st_seal%etime = st_sim%end_time
      else
        deallocate(read_seal)
      end if
    end if

    if (st_step_flag%rech == 1) then
      bound_name = "recharge"
      if (st_in_type%rech /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_in_type%rech, st_rech%fnum, st_rech%multi, bound_name,&
                         st_rech%totn, st_step_flag%rech, ierr, st_rech%etime)
      else if (st_in_type%rech == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_intre%type, st_intre%fnum, st_rech%multi, st_intre%step,&
                         bound_name, st_rech%fnum, st_step_flag%rech, ierr, st_rech%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_rech%fnum, bfview%rech, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%rech = 0
        st_rech%etime = st_sim%end_time
        ! -- Reset cell value (value)
          call reset_value(rech_num, read_rech)
      else
        deallocate(read_rech, rech2cals, calc_rech)
      end if
    end if

    if (st_step_flag%well == 1) then
      bound_name = "well"
      if (st_in_type%well /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_in_type%well, st_well%fnum, st_well%multi, bound_name,&
                         st_well%totn, st_step_flag%well, ierr, st_well%etime)
      else if (st_in_type%well == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_intwe%type, st_intwe%fnum, st_well%multi, st_intwe%step,&
                         bound_name, st_well%fnum, st_step_flag%well, ierr, st_well%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_well%fnum, bfview%well, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%well = 0
        st_well%etime = st_sim%end_time
        ! -- Reset cell value (value)
          call reset_value(well_num, read_well)
      else
        deallocate(read_well, calc_well, well_top, well_bott, well_index, well_conn)
        if (well_num /= 0) then
          deallocate(abyd_well)
        end if
      end if
    end if

    if (st_step_flag%prec == 1) then
      bound_name = "precipitation"
      if (st_in_type%prec /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_in_type%prec, st_prec%fnum, st_prec%multi, bound_name,&
                         st_prec%totn, st_step_flag%prec, ierr, st_prec%etime)
      else if (st_in_type%prec == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_intpr%type, st_intpr%fnum, st_prec%multi, st_intpr%step,&
                         bound_name, st_prec%fnum, st_step_flag%prec, ierr, st_prec%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_prec%fnum, bfview%prec, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%prec = 0
        st_prec%etime = st_sim%end_time
        ! -- Reset cell value (value)
          call reset_value(prec_num, read_prec)
      else
        if (prec_num > 0) then
          deallocate(read_prec)
        end if
      end if
    end if

    if (st_step_flag%evap == 1) then
      bound_name = "evapotranspiration"
      if (st_in_type%evap /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_in_type%evap, st_evap%fnum, st_evap%multi, bound_name,&
                         st_evap%totn, st_step_flag%evap, ierr, st_evap%etime)
      else if (st_in_type%evap == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_intev%type, st_intev%fnum, st_evap%multi, st_intev%step,&
                         bound_name, st_evap%fnum, st_step_flag%evap, ierr, st_evap%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_evap%fnum, bfview%evap, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%evap = 0
        st_evap%etime = st_sim%end_time
        ! -- Reset cell value (value)
          call reset_value(evap_num, read_evap)
      else
        if (evap_num > 0) then
          deallocate(read_evap)
        end if
      end if
    end if

    if (st_step_flag%riwl == 1) then
      bound_name = "river water level"
      if (st_rivf_type%wlev /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_rivf_type%wlev, st_riwl%fnum, st_riwl%multi, bound_name,&
                         st_riwl%totn, st_step_flag%riwl, ierr, st_riwl%etime)
      else if (st_rivf_type%wlev == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_riwl%inttype, st_riwl%intfnum, st_riwl%multi,&
                         st_riwl%intstep, bound_name, st_riwl%fnum, st_step_flag%riwl,&
                         ierr, st_riwl%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_riwl%fnum, rfview%wl, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%riwl = 0
        st_riwl%etime = st_sim%end_time
      else
        deallocate(cflag_riv%wl, criv%wl)
      end if
    end if

    if (st_step_flag%riwd == 1) then
      bound_name = "river water depth"
      if (st_rivf_type%wdep /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_rivf_type%wdep, st_riwd%fnum, st_riwd%multi, bound_name,&
                         st_riwd%totn, st_step_flag%riwd, ierr, st_riwd%etime)
      else if (st_rivf_type%wdep == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_riwd%inttype, st_riwd%intfnum, st_riwd%multi,&
                         st_riwd%intstep, bound_name, st_riwd%fnum, st_step_flag%riwd,&
                         ierr, st_riwd%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_riwd%fnum, rfview%wd, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%riwd = 0
        st_riwd%etime = st_sim%end_time
      end if
    end if

    if (st_step_flag%ribl == 1) then
      bound_name = "river bottom level"
      if (st_rivf_type%blev /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_rivf_type%blev, st_ribl%fnum, st_ribl%multi, bound_name,&
                         st_ribl%totn, st_step_flag%ribl, ierr, st_ribl%etime)
      else if (st_rivf_type%blev == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_ribl%inttype, st_ribl%intfnum, st_ribl%multi,&
                         st_ribl%intstep, bound_name, st_ribl%fnum, st_step_flag%ribl,&
                         ierr, st_ribl%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_ribl%fnum, rfview%bl, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%ribl = 0
        st_ribl%etime = st_sim%end_time
      else
        deallocate(cflag_riv%bl, criv%bl)
      end if
    end if

    if (st_step_flag%ride == 1) then
      bound_name = "river depth"
      if (st_rivf_type%dept /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_rivf_type%dept, st_ride%fnum, st_ride%multi, bound_name,&
                         st_ride%totn, st_step_flag%ride, ierr, st_ride%etime)
      else if (st_rivf_type%dept == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_ride%inttype, st_ride%intfnum, st_ride%multi,&
                         st_ride%intstep, bound_name, st_ride%fnum, st_step_flag%ride,&
                         ierr, st_ride%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_ride%fnum, rfview%de, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%ride = 0
        st_ride%etime = st_sim%end_time
      else
        if (st_ride%totn > 0) then
          deallocate(cflag_riv%de, criv%de)
        end if
      end if
    end if

    if (st_step_flag%riwi == 1) then
      bound_name = "river width"
      if (st_rivf_type%widt /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_rivf_type%widt, st_riwi%fnum, st_riwi%multi, bound_name,&
                         st_riwi%totn, st_step_flag%riwi, ierr, st_riwi%etime)
      else if (st_rivf_type%widt == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_riwi%inttype, st_riwi%intfnum, st_riwi%multi,&
                         st_riwi%intstep, bound_name, st_riwi%fnum, st_step_flag%riwi,&
                         ierr, st_riwi%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_riwi%fnum, rfview%wi, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%riwi = 0
        st_riwi%etime = st_sim%end_time
      else
        if (st_riwi%totn > 0) then
          deallocate(cflag_riv%wi, criv%wi)
        end if
      end if
    end if

    if (st_step_flag%rile == 1) then
      bound_name = "river length"
      if (st_rivf_type%leng /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_rivf_type%leng, st_rile%fnum, st_rile%multi, bound_name,&
                         st_rile%totn, st_step_flag%rile, ierr, st_rile%etime)
      else if (st_rivf_type%leng == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_rile%inttype, st_rile%intfnum, st_rile%multi,&
                         st_rile%intstep, bound_name, st_rile%fnum, st_step_flag%rile,&
                         ierr, st_rile%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_rile%fnum, rfview%le, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%rile = 0
        st_rile%etime = st_sim%end_time
      else
        if (st_rile%totn > 0) then
          deallocate(cflag_riv%le, criv%le)
        end if
      end if
    end if

    if (st_step_flag%lawl == 1) then
      bound_name = "lake water level"
      if (st_lakf_type%wlev /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_lakf_type%wlev, st_lawl%fnum, st_lawl%multi, bound_name,&
                         st_lawl%totn, st_step_flag%lawl, ierr, st_lawl%etime)
      else if (st_lakf_type%wlev == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_lawl%inttype, st_lawl%intfnum, st_lawl%multi,&
                         st_lawl%intstep, bound_name, st_lawl%fnum, st_step_flag%lawl,&
                         ierr, st_lawl%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_lawl%fnum, lfview%wl, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%lawl = 0
        st_lawl%etime = st_sim%end_time
      else
        deallocate(cflag_lak%wl, clak%wl)
      end if
    end if

    if (st_step_flag%lawd == 1) then
      bound_name = "lake water depth"
      if (st_lakf_type%wdep /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_lakf_type%wdep, st_lawd%fnum, st_lawd%multi, bound_name,&
                         st_lawd%totn, st_step_flag%lawd, ierr, st_lawd%etime)
      else if (st_lakf_type%wdep == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_lawd%inttype, st_lawd%intfnum, st_lawd%multi,&
                         st_lawd%intstep, bound_name, st_lawd%fnum, st_step_flag%lawd,&
                         ierr, st_lawd%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_lawd%fnum, lfview%wd, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%lawd = 0
        st_lawd%etime = st_sim%end_time
      else
        if (st_lawd%totn > 0) then
          deallocate(cflag_lak%wd, clak%wd)
        end if
      end if
    end if

    if (st_step_flag%labl == 1) then
      bound_name = "lake bottom level"
      if (st_lakf_type%blev /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_lakf_type%blev, st_labl%fnum, st_labl%multi, bound_name,&
                         st_labl%totn, st_step_flag%labl, ierr, st_labl%etime)
      else if (st_lakf_type%blev == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_labl%inttype, st_labl%intfnum, st_labl%multi,&
                         st_labl%intstep, bound_name, st_labl%fnum, st_step_flag%labl,&
                         ierr, st_labl%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_labl%fnum, lfview%bl, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%labl = 0
        st_labl%etime = st_sim%end_time
      else
        deallocate(cflag_lak%bl, clak%bl)
      end if
    end if

    if (st_step_flag%laar == 1) then
      bound_name = "lake area"
      if (st_lakf_type%area /= in_type(7)) then
        ! -- Read next time (next)
          call read_next(st_lakf_type%area, st_laar%fnum, st_laar%multi, bound_name,&
                         st_laar%totn, st_step_flag%laar, ierr, st_laar%etime)
      else if (st_lakf_type%area == in_type(7)) then
        ! -- Read next time for time interval file (intn)
          call read_intn(st_laar%inttype, st_laar%intfnum, st_laar%multi,&
                         st_laar%intstep, bound_name, st_laar%fnum, st_step_flag%laar,&
                         ierr, st_laar%etime)
#ifdef MPI_MSG
        ! -- Set real4 file view (real4_fview)
          call set_real4_fview(st_laar%fnum, lfview%ar, bound_name)
#endif
      end if

      if (ierr /= 0) then
        if (my_rank == 0) then
          write(str_time,'(f0.3)') now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%laar = 0
        st_laar%etime = st_sim%end_time
      else
        if (st_laar%totn > 0) then
          deallocate(cflag_lak%ar, clak%ar)
        end if
      end if
    end if

    if (allocated(bound_name)) then
      deallocate(bound_name)
    end if

  end subroutine set_nextet

  subroutine set_delt()
  !***************************************************************************************
  ! set_delt -- Set delta time
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    real(SP) :: min_step
    !-------------------------------------------------------------------------------------
    min_step = min(st_rech%etime, st_well%etime, st_seal%etime, st_prec%etime,&
                   st_evap%etime, st_riwl%etime, st_riwd%etime, st_ribl%etime,&
                   st_ride%etime, st_riwi%etime, st_rile%etime, st_lawl%etime,&
                   st_lawd%etime, st_labl%etime, st_laar%etime)

    if (min_step == current_t) then
      delt = delt
    else if (min_step < next_time .and. next_time < st_sim%end_time) then
      delt = real(min_step - current_t, kind=DP)
    else if (next_time > st_sim%end_time) then
      delt = real(st_sim%end_time - current_t, kind=DP)
    end if

    if (delt > st_sim%max_step) then
      delt = real(st_sim%max_step, kind=DP)
    end if

    next_time = current_t + real(delt, kind=SP)

  end subroutine set_delt

  subroutine set_nextvar()
  !***************************************************************************************
  ! set_nextvar -- Set next variable
  !***************************************************************************************
    ! -- modules
    use assign_calc, only: read_ksx, read_ksy
    use set_condition, only: set_srabyd, set_chabyd, set_wellconn
    use assign_boundary, only: assign_sealv, assign_surfbv, assign_wellv, assign_rilav,&
                               rech_cflag, prec_cflag, evap_cflag
    use calc_boundary, only: calc_reprev, calc_rivea, calc_wlbd, conv_rech2calc,&
                             count_rivecalc, count_lakecalc, rive2cals, lake2cals,&
                             rive_head, lake_head, rive_bott, rive_area, lake_bott,&
                             lake_area
    use set_boundary, only: rive_num, lake_num, rivnum, laknum, abyd_rive, abyd_lake
#ifdef MPI_MSG
    use mpi_utility, only: mpisum_val
#endif
    ! -- inout

    ! -- local

    integer(I4) :: prec_stepflag, evap_stepflag, rive_stepflag, lake_stepflag
    integer(I4) :: rive_aflag, sum_ribln, sum_riwln, sum_lawln, sum_labln
    real(DP), allocatable :: temp_area(:)
    !-------------------------------------------------------------------------------------
    if (st_step_flag%seal == 1) then
      ! -- Assign sea level value (sealv)
        call assign_sealv(st_in_type%seal)

      st_step_flag%seal = 0

    else if (st_seal%etime == next_time) then
      st_step_flag%seal = 1
    end if

    if (st_step_flag%rech == 1) then
      allocate(rech_cflag(ncals), read_rech(ncals))
      !$omp parallel workshare
      rech_cflag(:) = 0 ; read_rech(:) = SNOVAL
      !$omp end parallel workshare
      ! -- Assign recharge value
        call assign_surfbv(st_in_type%rech, st_intre%type, st_rech, rech_num, rech_cflag,&
                           read_rech)

      call conv_rech2calc(rech_num)

      st_step_flag%rech = 0
    else if (st_rech%etime == next_time) then
      st_step_flag%rech = 1
    end if

    if (st_step_flag%well == 1) then
      ! -- Assign well value (wellv)
        call assign_wellv(st_in_type%well, st_in_type%weks, st_in_type%weke, well_num)

      st_step_flag%well = 0
      if (well_num /= 0) then
        ! -- Set well connectivity (wellconn)
          call set_wellconn(well_num, read_ksx, read_ksy)
      end if

    else if (st_well%etime == next_time) then
      st_step_flag%well = 1
    end if

    prec_stepflag = 0

    if (st_step_flag%prec == 1) then
      allocate(prec_cflag(ncals), read_prec(ncals))
      !$omp parallel workshare
      prec_cflag(:) = 0 ; read_prec(:) = SNOVAL
      !$omp end parallel workshare
      ! -- Assign precipitation value
        call assign_surfbv(st_in_type%prec, st_intpr%type, st_prec, prec_num, prec_cflag,&
                           read_prec)
      prec_stepflag = prec_stepflag + 1
      deallocate(prec_cflag)

      st_step_flag%prec = 0
    else if (st_prec%etime == next_time) then
      st_step_flag%prec = 1
    end if

    evap_stepflag = 0

    if (st_step_flag%evap == 1) then
      allocate(evap_cflag(ncals), read_evap(ncals))
      !$omp parallel workshare
      evap_cflag(:) = 0 ; read_evap(:) = SNOVAL
      !$omp end parallel workshare
      ! -- Assign evapotranspiration value
        call assign_surfbv(st_in_type%evap, st_intev%type, st_evap, evap_num, evap_cflag,&
                           read_evap)
      evap_stepflag = evap_stepflag + 1
      deallocate(evap_cflag)

      st_step_flag%evap = 0
    else if (st_evap%etime == next_time) then
      st_step_flag%evap = 1
    end if

    if (prec_stepflag /= 0 .or. evap_stepflag /= 0) then
      deallocate(read_rech, rech2cals, calc_rech)
      ! -- Calculate recharge from precipitation and evapotranspiration (reprev)
        call calc_reprev(rech_num)
      call conv_rech2calc(rech_num)
    end if

    rive_stepflag = 0 ; rive_aflag = 0

    if (st_step_flag%riwl == 1) then
      allocate(cflag_riv%wl(ncals), criv%wl(ncals))
      !$omp parallel workshare
      cflag_riv%wl(:) = 0 ; criv%wl(:) = SNOVAL
      !$omp end parallel workshare
      ! -- Assign river water level value
        call assign_rilav(st_rivf_type%wlev, 0, st_riwl, rivnum%wl, cflag_riv%wl, criv%wl)

      st_step_flag%riwl = 0 ; rive_stepflag = rive_stepflag + 1
    else if (st_riwl%etime == next_time) then
      st_step_flag%riwl = 1
    end if

    if (st_step_flag%ribl == 1) then
      allocate(cflag_riv%bl(ncals), criv%bl(ncals))
      !$omp parallel workshare
      cflag_riv%bl(:) = 0 ; criv%bl(:) = SNOVAL
      !$omp end parallel workshare
      ! -- Assign river bottom level value
        call assign_rilav(st_rivf_type%blev, 0, st_ribl, rivnum%bl, cflag_riv%bl, criv%bl)

      st_step_flag%ribl = 0 ; rive_stepflag = rive_stepflag + 1
    else if (st_ribl%etime == next_time) then
      st_step_flag%ribl = 1
    end if

    if (st_step_flag%riwd == 1) then
      if (st_riwd%totn > 0) then
        allocate(cflag_riv%wd(ncals), criv%wd(ncals))
        !$omp parallel workshare
        cflag_riv%wd(:) = 0 ; criv%wd(:) = SNOVAL
        !$omp end parallel workshare
      end if
      ! -- Assign river water depth value
        call assign_rilav(st_rivf_type%wdep, 0, st_riwd, rivnum%wd, cflag_riv%wd, criv%wd)

#ifdef MPI_MSG
      ! -- Sum value for MPI (val)
        call mpisum_val(rivnum%wl, "river water level", sum_riwln)
        call mpisum_val(rivnum%bl, "river bottom level", sum_ribln)
#else
      sum_riwln = rivnum%wl ; sum_ribln = rivnum%bl
#endif

      if (sum_riwln == 0 .and. sum_ribln /= 0) then
        ! -- Calculate water level from river bottom level (wlrb)
          call calc_wlbd(cflag_riv%bl, criv%bl, cflag_riv%wd, criv%wd, cflag_riv%wl,&
                         criv%wl, rivnum%wl)
        rivnum%wl = 0
        deallocate(cflag_riv%wd, criv%wd)
      end if

      st_step_flag%riwd = 0 ; rive_stepflag = rive_stepflag + 1
    else if (st_riwd%etime == next_time) then
      st_step_flag%riwd = 1
    end if

    if (st_step_flag%riwi == 1) then
      if (st_riwi%totn > 0) then
        allocate(cflag_riv%wi(ncals), criv%wi(ncals))
        !$omp parallel workshare
        cflag_riv%wi(:) = 0 ; criv%wi(:) = SNOVAL
        !$omp end parallel workshare
      end if
      ! -- Assign river width value
        call assign_rilav(st_rivf_type%widt, 0, st_riwi, rivnum%wi, cflag_riv%wi, criv%wi)

      st_step_flag%riwi = 0 ; rive_aflag = rive_aflag + 1
    else if (st_riwi%etime == next_time) then
      st_step_flag%riwi = 1
    end if

    if (st_step_flag%rile == 1) then
      if (st_rile%totn > 0) then
        allocate(cflag_riv%le(ncals), criv%le(ncals))
        !$omp parallel workshare
        cflag_riv%le(:) = 0 ; criv%le(:) = SNOVAL
        !$omp end parallel workshare
      end if
      ! -- Assign river length value
        call assign_rilav(st_rivf_type%leng, 0, st_rile, rivnum%le, cflag_riv%le, criv%le)

      st_step_flag%rile = 0 ; rive_aflag = rive_aflag + 1
    else if (st_rile%etime == next_time) then
      st_step_flag%rile = 1
    end if

    if (rive_aflag > 0) then
      !$omp parallel workshare
      cflag_riv%ar(:) = 0 ; criv%ar(:) = SNOVAL
      !$omp end parallel workshare
      ! -- Calculate river area (rivea)
        call calc_rivea(cflag_riv%wi, cflag_riv%le, criv%wi, criv%le, cflag_riv%ar,&
                        criv%ar, rivnum%ar)
    end if

    if (rive_stepflag > 0) then
      if (rive_num /= 0) then
        allocate(temp_area(rive_num))
        !$omp parallel workshare
        temp_area(:) = -rive_area(:)
        !$omp end parallel workshare
        ! -- Set surface&recharge area and area by distance (srabyd)
          call set_srabyd(rive_num, rive_bott, temp_area, rive2cals, abyd_rive)
        deallocate(rive_head, rive_bott, rive_area, temp_area)
      end if
      deallocate(rive2cals)
      ! -- Count river calculation (rivecalc)
        call count_rivecalc(cflag_riv%wl, cflag_riv%bl, cflag_riv%ar, criv%wl, criv%bl,&
                            criv%ar, rive_num)
      if (rive_num /= 0) then
        deallocate(abyd_rive)
        allocate(abyd_rive(rive_num))
        !$omp parallel workshare
        abyd_rive(:) = DZERO
        !$omp end parallel workshare
        ! -- Set surface&recharge area and area by distance (srabyd)
          call set_srabyd(rive_num, rive_bott, rive_area, rive2cals, abyd_rive)
      end if
    end if

    lake_stepflag = 0

    if (st_step_flag%lawl == 1) then
      allocate(cflag_lak%wl(ncals), clak%wl(ncals))
      cflag_lak%wl(:) = 0 ; clak%wl(:) = SNOVAL
      ! -- Assign lake water level value
        call assign_rilav(st_lakf_type%wlev, 0, st_lawl, laknum%wl, cflag_lak%wl, clak%wl)

      st_step_flag%lawl = 0 ; lake_stepflag = lake_stepflag + 1
    else if (st_lawl%etime == next_time) then
      st_step_flag%lawl = 1
    end if

    if (st_step_flag%labl == 1) then
      allocate(cflag_lak%bl(ncals), clak%bl(ncals))
      !$omp parallel workshare
      cflag_lak%bl(:) = 0 ; clak%bl(:) = SNOVAL
      !$omp end parallel workshare
      ! -- Assign lake bottom level value
        call assign_rilav(st_lakf_type%blev, 0, st_labl, laknum%bl, cflag_lak%bl, clak%bl)

      st_step_flag%labl = 0 ; lake_stepflag = lake_stepflag + 1
    else if (st_labl%etime == next_time) then
      st_step_flag%labl = 1
    end if

    if (st_step_flag%lawd == 1) then
      if (st_lawd%totn > 0) then
        allocate(cflag_lak%wd(ncals), clak%wd(ncals))
        !$omp parallel workshare
        cflag_lak%wd(:) = 0 ; clak%wd(:) = SNOVAL
        !$omp end parallel workshare
      end if
      ! -- Assign lake water depth value
        call assign_rilav(st_lakf_type%wdep, 0, st_lawd, laknum%wd, cflag_lak%wd, clak%wd)

#ifdef MPI_MSG
      ! -- Sum value for MPI (val)
        call mpisum_val(laknum%wl, "lake water level", sum_lawln)
        call mpisum_val(laknum%bl, "lake bottom level", sum_labln)
#else
      sum_lawln = laknum%wl ; sum_labln = laknum%bl
#endif

      if (sum_lawln == 0 .and. sum_labln /= 0) then
        ! -- Calculate water level from bottom level and water depth (wlbd)
          call calc_wlbd(cflag_lak%bl, clak%bl, cflag_lak%wd, clak%wd, cflag_lak%wl,&
                         clak%wl, laknum%wl)
        laknum%wl = 0
        deallocate(cflag_lak%wd, clak%wd)
      end if

      st_step_flag%lawd = 0 ; lake_stepflag = lake_stepflag + 1
    else if (st_lawd%etime == next_time) then
      st_step_flag%lawd = 1
    end if

    if (st_step_flag%laar == 1) then
      if (st_laar%totn > 0) then
        allocate(cflag_lak%ar(ncals), clak%ar(ncals))
        !$omp parallel workshare
        cflag_lak%ar(:) = 0 ; clak%ar(:) = SNOVAL
        !$omp end parallel workshare
      end if
      ! -- Assign lake area value
        call assign_rilav(st_lakf_type%area, 1, st_laar, laknum%ar, cflag_lak%ar, clak%ar)

      st_step_flag%laar = 0 ; lake_stepflag = lake_stepflag + 1
    else if (st_laar%etime == next_time) then
      st_step_flag%laar = 1
    end if

    if (lake_stepflag > 0) then
      if (lake_num /= 0) then
        allocate(temp_area(lake_num))
        !$omp parallel workshare
        temp_area(:) = -lake_area(:)
        !$omp end parallel workshare
        ! -- Set surface&recharge area and area by distance (srabyd)
          call set_srabyd(lake_num, lake_bott, temp_area, lake2cals, abyd_lake)
        deallocate(lake_head, lake_bott, lake_area, temp_area)
      end if
      deallocate(lake2cals)
      ! -- Count lake calculation cell (lakecalc)
        call count_lakecalc(cflag_lak%wl, cflag_lak%bl, cflag_lak%ar, clak%wl, clak%bl,&
                            clak%ar, lake_num)
      if (lake_num /= 0) then
        deallocate(abyd_lake)
        allocate(abyd_lake(lake_num))
        !$omp parallel workshare
        abyd_lake(:) = DZERO
        !$omp end parallel workshare
        ! -- Set surface&recharge area and area by distance (srabyd)
          call set_srabyd(lake_num, lake_bott, lake_area, lake2cals, abyd_lake)
      end if
    end if

    if (rive_stepflag > 0 .or. lake_stepflag > 0) then
      ! -- Set charge area by distance (chabyd)
        call set_chabyd()
    end if

  end subroutine set_nextvar

  subroutine set_vwell_head()
  !***************************************************************************************
  ! set_vwell_head -- Set virtual well head
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    ! -- Calculate virtual well head without well pumping (vheadout)
      call calc_vheadout()

    !$omp parallel workshare
    calc_well(:) = DZERO
    !$omp end parallel workshare

    ! -- Calculate virtual well head with well puming (vheadin)
      call calc_vheadin()

  end subroutine set_vwell_head

  subroutine change_recharge()
  !***************************************************************************************
  ! change_recharge -- Change the recharge volume
  !***************************************************************************************
    ! -- modules
    use assign_calc, only: surf_bott, surf_parm, surf_reli
    ! -- inout

    ! -- local
    integer(I4) :: i, s
    real(DP) :: norm_elev
    real(DP), allocatable :: water_dep(:), rech_rati(:)
    !-------------------------------------------------------------------------------------
    allocate(water_dep(ncals), rech_rati(rech_num))
    !$omp parallel
    !$omp workshare
    water_dep(:) = DZERO ; rech_rati(:) = DONE
    !$omp end workshare

    !$omp do private(i, s, norm_elev)
    do i = 1, rech_num
      s = rech2cals(i)
      water_dep(s) = head_new(s) - surf_bott(s)
      if (read_rech(i) < SZERO .or. water_dep(s) < DZERO) then
        rech_rati(i) = DONE
      else if (water_dep(s) < surf_reli(s) .and. surf_reli(s) /= DZERO) then
        norm_elev = water_dep(s)/surf_reli(s)
        rech_rati(i) = norm_elev**surf_parm(s)
        rech_rati(i) = DONE - rech_rati(i)
      else if (water_dep(s) >= surf_reli(s) .and. surf_reli(s) /= DZERO) then
        rech_rati(i) = DZERO
      else
        rech_rati(i) = DONE
      end if
      calc_rech(i) = read_rech(i)*rech_area(s)*rech_rati(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(water_dep, rech_rati)

  end subroutine change_recharge

!  subroutine change_volrate()
!  !***************************************************************************************
!  ! change_volrate -- Change the volume rate using time step and cell volume
!  !***************************************************************************************
!    ! -- modules
!    use make_cell, only: surf_elev
!    ! -- inout
!
!    ! -- local
!    integer(I4) :: i, s
!    !-------------------------------------------------------------------------------------
!    !$omp parallel do private(i, s)
!    do i = 1, rech_num
!      s = rech2cals(i)
!      if (surf_head(s) >= surf_elev(s)) then !exist surface water
!        surf_head(s) = surf_head(s) + read_rech(i)*delt
!        if (surf_head(s) < surf_elev(s)) then !below groundlevel
!          calc_rech(i) = (surf_head(s)-surf_elev(s))/delt*rech_area(s)
!          surf_head(s) = surf_elev(s)
!        else
!          calc_rech(i) = DZERO
!        end if
!      else
!        calc_rech(i) = read_rech(i)*rech_area(s)
!      end if
!    end do
!    !$omp end parallel do
!
!  end subroutine change_volrate

!  subroutine set_est_head(old2_head, old1_head, new_head)
!  !***************************************************************************************
!  ! set_est_head -- Set estimate head
!  !***************************************************************************************
!    ! -- modules
!
!    ! -- inout
!    real(DP), intent(in) :: old1_head(:), old2_head(:)
!    real(DP), intent(out) :: new_head(:)
!    ! -- local
!    integer(I4) :: i
!    real(DP) :: delt2_inv
!    !-------------------------------------------------------------------------------------
!    delt2_inv = DONE/delt_old2
!    !$omp parallel do private(i)
!    do i = 1, ncalc
!      new_head(i) = old1_head(i) + delt_old1*((old1_head(i)-old2_head(i))*delt2_inv)
!    end do
!    !$omp end parallel do
!
!  end subroutine set_est_head

  subroutine calc_vheadout()
  !***************************************************************************************
  ! calc_vheadout -- Calculate virtual well head without well pumping
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k
    real(DP) :: tot_cond, tot_flux
    !-------------------------------------------------------------------------------------
    allocate(whead_new(well_num))
    !$omp parallel
    !$omp workshare
    whead_new(:) = DZERO
    !$omp end workshare

    !$omp do private(i, j, k, tot_cond, tot_flux)
    do i = 1, well_num
      tot_cond = DZERO ; tot_flux = DZERO
      do k = well_index(i-1)+1, well_index(i)
        j = well_conn(k)
        tot_cond = tot_cond + rel_perm(j)*abyd_well(k)
        tot_flux = tot_flux + rel_perm(j)*abyd_well(k)*head_new(j)
      end do
      if (tot_cond /= DZERO) then
        whead_new(i) = tot_flux/tot_cond
      end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine calc_vheadout

  subroutine calc_vheadin()
  !***************************************************************************************
  ! calc_vheadin -- Calculate virtual well head with well pumping
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, j, k
    real(DP) :: tot_cond
    real(DP), allocatable :: temp_whead(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_whead(well_num))
    !$omp parallel
    !$omp workshare
    temp_whead(:) = whead_new(:)
    !$omp end workshare

    !$omp do private(i, j, k, tot_cond)
    do i = 1, well_num
      tot_cond = DZERO
      do k = well_index(i-1)+1, well_index(i)
        j = well_conn(k)
        tot_cond = tot_cond + rel_perm(j)*abyd_well(k)
      end do
      if (tot_cond /= DZERO) then
        whead_new(i) = temp_whead(i) + read_well(i)/tot_cond

!        if (whead_new(i) > well_top(i)) then
!          whead_new(i) = well_top(i)
!        else if (whead_new(i) < well_bott(i)) then
!          whead_new(i) = well_bott(i)
!        end if

        do k = well_index(i-1)+1, well_index(i)
          j = well_conn(k)
          if (whead_new(i) > well_bott(i) .or. head_new(j) > well_bott(i)) then
            calc_well(j) = calc_well(j) + rel_perm(j)*abyd_well(k)*(whead_new(i)-head_new(j))
          end if
        end do
      end if
    end do
    !$omp end do
    !$omp end parallel

    deallocate(whead_new, temp_whead)

  end subroutine calc_vheadin

  subroutine set_valexc(excn, exc_inval, exc_outval)
  !***************************************************************************************
  ! set_valexc -- Set value exchange
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: excn
    real(DP), intent(in) :: exc_inval(:)
    real(DP), intent(out) :: exc_outval(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, excn
      exc_outval(i) = exc_inval(i)
    end do
    !$omp end parallel do

  end subroutine set_valexc

  subroutine write_bound_change(bchange)
  !***************************************************************************************
  ! write_bound_change -- Write boundary change information
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(out) :: bchange
    ! -- local
    integer(I4) :: conv_fnum
    character(9) :: cond_format
    !-------------------------------------------------------------------------------------
    bchange = st_step_flag%rech + st_step_flag%well + st_step_flag%seal +&
              st_step_flag%prec + st_step_flag%evap + st_step_flag%riwl +&
              st_step_flag%riwd + st_step_flag%ribl + st_step_flag%ride +&
              st_step_flag%riwi + st_step_flag%lawl + st_step_flag%lawd +&
              st_step_flag%labl + st_step_flag%laar

    conv_fnum = st_out_fnum%conv ; cond_format = "(a,f10.3)"

    !$omp parallel shared(conv_fnum, cond_format)
    !$omp single
    if (st_step_flag%rech == 1) then
      write(conv_fnum,cond_format) "Changed recharge condition at ", now_time
    end if
    if (st_step_flag%well == 1) then
      write(conv_fnum,cond_format) "Changed well condition at ", now_time
    end if
    if (st_step_flag%seal == 1) then
      write(conv_fnum,cond_format) "Changed sea condition at ", now_time
    end if
    if (st_step_flag%prec == 1) then
      write(conv_fnum,cond_format) "Changed precipitation condition at ", now_time
    end if
    if (st_step_flag%evap == 1) then
      write(conv_fnum,cond_format) "Changed evapotranspiration condition at ", now_time
    end if
    if (st_step_flag%riwl == 1) then
      write(conv_fnum,cond_format) "Changed river water level at ", now_time
    end if
    if (st_step_flag%riwd == 1) then
      write(conv_fnum,cond_format) "Changed river water depth at ", now_time
    end if
    if (st_step_flag%ribl == 1) then
      write(conv_fnum,cond_format) "Changed river bottom level at ", now_time
    end if
    if (st_step_flag%ride == 1) then
      write(conv_fnum,cond_format) "Changed river depth at ", now_time
    end if
    if (st_step_flag%riwi == 1) then
      write(conv_fnum,cond_format) "Changed river width at ", now_time
    end if
    if (st_step_flag%lawl == 1) then
      write(conv_fnum,cond_format) "Changed lake water level at ", now_time
    end if
    if (st_step_flag%lawd == 1) then
      write(conv_fnum,cond_format) "Changed lake water depth at ", now_time
    end if
    if (st_step_flag%labl == 1) then
      write(conv_fnum,cond_format) "Changed lake bottom level at ", now_time
    end if
    if (st_step_flag%laar == 1) then
      write(conv_fnum,cond_format) "Changed lake area at ", now_time
    end if
    !$omp end single
    !$omp end parallel

  end subroutine write_bound_change

  subroutine reset_stepf()
  !***************************************************************************************
  ! reset_stepf -- Reset stepflag
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    !$omp parallel sections
    !$omp section
    if (st_step_flag%rech == 1) then
      st_step_flag%rech = 0
    end if
    !$omp section
    if (st_step_flag%well == 1) then
      st_step_flag%well = 0
    end if
    !$omp section
    if (st_step_flag%seal == 1) then
      st_step_flag%seal = 0
    end if
    !$omp section
    if (st_step_flag%prec == 1) then
      st_step_flag%prec = 0
    end if
    !$omp section
    if (st_step_flag%evap == 1) then
      st_step_flag%evap = 0
    end if
    !$omp section
    if (st_step_flag%riwl == 1) then
      st_step_flag%riwl = 0
    end if
    !$omp section
    if (st_step_flag%riwd == 1) then
      st_step_flag%riwd = 0
    end if
    !$omp section
    if (st_step_flag%ribl == 1) then
      st_step_flag%ribl = 0
    end if
    !$omp section
    if (st_step_flag%ride == 1) then
      st_step_flag%ride = 0
    end if
    !$omp section
    if (st_step_flag%riwi == 1) then
      st_step_flag%riwi = 0
    end if
    !$omp section
    if (st_step_flag%lawl == 1) then
      st_step_flag%lawl = 0
    end if
    !$omp section
    if (st_step_flag%lawd == 1) then
      st_step_flag%lawd = 0
    end if
    !$omp section
    if (st_step_flag%labl == 1) then
      st_step_flag%labl = 0
    end if
    !$omp section
    if (st_step_flag%laar == 1) then
      st_step_flag%laar = 0
    end if
    !$omp end parallel sections

  end subroutine reset_stepf

  subroutine reset_value(targn, targ_value)
  !***************************************************************************************
  ! reset_value -- Reset cell value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: targn
    real(SP), intent(out) :: targ_value(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, targn
      targ_value(i) = SZERO
    end do
    !$omp end parallel do

  end subroutine reset_value

  subroutine apply_heuri(out_num)
  !***************************************************************************************
  ! apply_heuri -- Apply heuristic time stepping
  !***************************************************************************************
    ! -- modules
    use initial_module, only: maxout_iter
    ! -- inout
    integer(I4), intent(in) :: out_num
    ! -- local
    integer(I4) :: incr_num, decr_num
    real(DP) :: incr_fac = 1.2E+00_DP, decr_fac = 0.8E+00_DP
    !-------------------------------------------------------------------------------------
    incr_num = int(maxout_iter*0.4) ; decr_num = int(maxout_iter*0.8)
    if (out_num <= incr_num) then
      delt = delt_old1*incr_fac
    else if (out_num <= decr_num) then
      delt = delt_old1
    else
      delt = delt_old1*decr_fac
    end if

  end subroutine apply_heuri

!  subroutine apply_adapt(old1_head, old2_head, new_head)
!  !***************************************************************************************
!  ! apply_adapt -- Apply adaptive time stepping
!  !***************************************************************************************
!    ! -- modules
!    use constval_module, only: MACHI_EPS
!    use initial_module, only: criteria
!    use read_input, only: len_scal_inv
!    ! -- inout
!    real(DP), intent(in) :: old1_head(:), old2_head(:), new_head(:)
!    ! -- local
!    integer(I4) :: loop_max
!    real(DP) :: tru_err, max_err
!    real(DP) :: temp_step, rmin = 0.1_DP
!    real(DP) :: rmax, max_tau
!    !-------------------------------------------------------------------------------------
!    rmax = st_sim%inc_fact ; max_tau = criteria*10_DP*len_scal_inv
!    ! -- Calculate truncation error (trunerr)
!      call calc_trunerr(max_tau, new_head, old1_head, old2_head, tru_err, max_err)
!
!    loop_max = 1 ; temp_step = delt_old1
!    trun_loop: do while (tru_err >= DZERO)
!      if (loop_max == 10) then
!        exit trun_loop
!      end if
!      delt_old1 = delt_old1*max(sqrt(max_tau/max(max_err, MACHI_EPS)), rmin)
!      ! -- Calculate truncation error (trunerr)
!        call calc_trunerr(max_tau, new_head, old1_head, old2_head, tru_err, max_err)
!      loop_max = loop_max + 1
!    end do trun_loop
!
!    if (loop_max == 10) then
!      delt = temp_step*rmin
!    else
!      delt = delt_old1*min(sqrt(max_tau/max(max_err, MACHI_EPS)), rmax)
!    end if
!
!  end subroutine apply_adapt

!  subroutine calc_trunerr(abs_err, newh, old1h, old2h, maxterr, maxerr)
!  !***************************************************************************************
!  ! calc_trunerr -- Calculate truncation error
!  !***************************************************************************************
!    ! -- modules
!#ifdef MPI_MSG
!    use mpi_utility, only: mpimax_val
!#endif
!    ! -- inout
!    real(DP), intent(in) :: abs_err
!    real(DP), intent(in) :: newh(:), old1h(:), old2h(:)
!    real(DP), intent(out) :: maxterr, maxerr
!    ! -- local
!    integer(I4) :: i
!    real(DP) :: delt1_inv, delt2_inv, trun_crit, err_max
!    real(DP), allocatable :: deri_time1(:), deri_time2(:), truc_error(:)
!    !-------------------------------------------------------------------------------------
!    allocate(deri_time1(ncalc), deri_time2(ncalc), truc_error(ncalc))
!    !$omp parallel workshare
!    deri_time1(:) = DZERO ; deri_time2(:) = DZERO ; truc_error(:) = DZERO
!    !$omp end parallel workshare
!
!    if (delt_old2 /= DZERO) then
!      delt2_inv = DONE/delt_old2
!      !$omp parallel do private(i)
!      do i = 1, ncalc
!        deri_time2(i) = (old1h(i) - old2h(i))*delt2_inv
!      end do
!      !$omp end parallel do
!    end if
!
!    delt1_inv = DONE/delt_old1
!    !$omp parallel do private(i)
!    do i = 1, ncalc
!      deri_time1(i) = (newh(i) - old1h(i))*delt1_inv
!      truc_error(i) = DHALF*delt_old1*abs(deri_time1(i)-deri_time2(i))
!    end do
!    !$omp end parallel do
!
!    trun_crit = truc_error(1) - abs_err ; err_max =  truc_error(1)
!    !$omp parallel do private(i) reduction(max:trun_crit, err_max)
!    do i = 1, ncalc
!      if ((truc_error(i)-abs_err) > trun_crit) then
!        trun_crit = truc_error(i) - abs_err ; err_max =  truc_error(i)
!      end if
!    end do
!    !$omp end parallel do
!
!#ifdef MPI_MSG
!    ! -- MAX value for MPI (val)
!      call mpimax_val(trun_crit, "truncation criteria", maxterr)
!      call mpimax_val(err_max, "truncation error", maxerr)
!#else
!    maxterr = trun_crit ; maxerr = err_max
!#endif
!
!    deallocate(deri_time1, deri_time2, truc_error)
!
!  end subroutine calc_trunerr

  subroutine set_date(inttime, ndate, restime)
  !***************************************************************************************
  ! set_date -- Set next date
  !***************************************************************************************
    ! -- modules
    use constval_module, only: MINSEC, HOURSEC, DAYSEC
    use utility_module, only: get_days
    ! -- inout
    real(SP), intent(in) :: inttime
    integer(I4), intent(inout) :: ndate(:)
    real(SP), intent(out) :: restime
    ! -- local
    integer(I4) :: i, monday
    integer(I4) :: isec, imin, ihour, iday
    real(SP) :: rsec
    !-------------------------------------------------------------------------------------
    iday = int(inttime/DAYSEC) ; rsec = inttime - iday*DAYSEC
    ihour = int(rsec/HOURSEC) ; rsec = rsec - ihour*HOURSEC
    imin = int(rsec/MINSEC) ; rsec = rsec - imin*MINSEC
    isec = int(rsec) ; restime = rsec - isec

    ndate(4) = ndate(4) + ihour ; ndate(5) = ndate(5) + imin ; ndate(6) = ndate(6) + isec
    if (ndate(6) >= 60) then
      ndate(5) = ndate(5) + 1 ; ndate(6) = ndate(6) - 60
    end if
    if (ndate(5) >= 60) then
      ndate(4) = ndate(4) + 1 ; ndate(5) = ndate(5) - 60
    end if
    if (ndate(4) >= 24) then
      ndate(3) = ndate(3) + 1 ; ndate(4) = ndate(4) - 24
    end if

    monday = get_days(ndate(1), ndate(2))
    do i = 1, iday
      ndate(3) = ndate(3) + 1
      if (ndate(3) > monday) then
        ndate(2) = ndate(2) + 1 ; ndate(3) = 1
        if (ndate(2) > 12) then
          ndate(1) = ndate(1) + 1 ; ndate(2) = 1
        end if
        monday = get_days(ndate(1), ndate(2))
      end if
    end do

  end subroutine set_date

end module time_module
