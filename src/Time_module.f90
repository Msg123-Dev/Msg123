module time_module
  ! -- modules
  use kind_module, only: I4, SP, DP
  use types_module, only: sol_set
  use constval_module, only: SZERO, SNOVAL, DZERO, DONE, DHALF
  use utility_module, only: st_mpi
  use initial_module, only: st_sim, st_ctrl, st_in_type, st_rivf_type, st_lakf_type, st_seal,&
                            st_rech, st_prec, st_evap, st_well, st_riwl, st_riwd, st_ribl,&
                            st_ride, st_riwi, st_rile, st_lawl, st_lawd, st_labl, st_laar,&
                            st_step_flag
  use open_file, only: st_intre, st_intpr, st_intev
  use check_condition, only: st_out_fnum
  use set_cell, only: ncalc, ncals
  use set_condition, only: st_hydr, st_bcnd
  use prep_calculation, only: st_time
  use assign_boundary, only: st_forc
  use set_boundary, only: st_rive, st_lake

  implicit none
  private
  public :: update_tstep

  ! -- local
  integer(I4) :: timestep_num, boundstep
  real(SP) :: next_time
  real(DP) :: delt_old1
  real(DP), allocatable :: whead_new(:)

  contains

  subroutine update_tstep(st_sol)
  !*********************************************************************************************
  ! update_tstep -- Update time step
  !*********************************************************************************************
    ! -- module
    use constval_module, only: DSMAL
    use utility_module, only: write_err_stop
    use assign_calc, only: geog_num
#ifdef ICI
    use ici_module, only: get_var
#endif
    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i, s
    !-------------------------------------------------------------------------------------------
    ! -- Calculate next time step (nextst)
      call calc_nextst(st_sol)

    if (st_time%conv_flag .and. st_sim%sim_type == 1) then
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
      st_time%delt = real(st_sim%end_time - st_time%current_t, kind=DP)
    end if

    if (st_time%delt < DSMAL .and. st_sim%sim_type /= -1 .and. st_mpi%rank == 0) then
      call write_err_stop("Time Step is too small.")
    else if (st_sim%sim_type /= -1) then
      st_time%current_t = st_time%current_t + real(st_time%delt, kind=SP)
    else if (st_sim%sim_type == -1) then
      st_time%current_t = SZERO
    end if

    st_time%now_time = st_time%current_t/st_sim%cal_fact
    st_time%delt_inv = DONE/st_time%delt

    if (st_time%conv_flag .or. timestep_num == 0) then
      if (st_bcnd%well_num /= 0) then
        ! -- Set virtual well head (vwell_head)
          call set_vwell_head(st_sol)
      end if
      if (st_bcnd%rech_num /= 0) then
        !$omp parallel do private(i, s)
        do i = 1, st_bcnd%rech_num
          s = st_bcnd%rech2cals(i)
          st_forc%calc_rech(i) = st_forc%read_rech(i)*st_hydr%rech_area(s)
        end do
        !$omp end parallel do
        if (geog_num /= 0) then
          ! -- Change the recharge volume
            call change_recharge(st_sol)
        end if
      end if
    end if

  end subroutine update_tstep

  subroutine calc_nextst(st_sol)
  !*********************************************************************************************
  ! calc_nextst -- Calculate next time step
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DNOVAL
    use utility_module, only: conv_unit
    use initial_module, only: st_init
    use read_input, only: len_scal
    use calc_parameter, only: calc_srat_rperm
#ifdef MPI_MSG
    use mpi_utility, only: bcast_val
    use mpi_write, only: write_mpi_3dbin
#else
    use write_module, only: write_3dbin, write_header_bin
#endif
    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: calc2calc(:)
    real(SP), save :: resi_time
    real(DP), allocatable :: cell_srat(:)
    !-------------------------------------------------------------------------------------------
    if (st_time%current_t == DZERO) then
      timestep_num = 0 ; boundstep = 0
      delt_old1 = DZERO ; resi_time = SZERO
      if (st_sim%sim_type == -1) then
        st_time%delt = DNOVAL
      else
        st_time%delt = st_sim%ini_step
      end if

      allocate(cell_srat(ncalc), calc2calc(ncalc))
      !$omp parallel do private(i)
      do i = 1, ncalc
        cell_srat(i) = DZERO
        calc2calc(i) = 0
      end do
      !$omp end parallel do
      ! -- Calculate saturation and relative permeability (srat_rperm)
        call calc_srat_rperm(ncalc, DZERO, st_forc%read_head, cell_srat, st_sol%rel_perm)
      !$omp parallel
      !$omp do private(i)
      do i = 1, ncalc
        st_sol%head_old(i) = st_forc%read_head(i) ; st_sol%head_new(i) = st_forc%read_head(i)
        st_sol%srat_old(i) = cell_srat(i) ; st_sol%srat_new(i) = cell_srat(i)
        calc2calc(i) = i
      end do
      !$omp end do

      !$omp do private(i)
      do i = 1, ncals
        st_sol%surf_head(i) = st_forc%read_head(i)
      end do
      !$omp end do
      !$omp end parallel

      deallocate(cell_srat)

      if (st_sim%sim_type /= -1) then
        deallocate(st_forc%read_head)
      end if

      !$omp parallel do private(i)
      do i = 1, ncals
        st_sol%surf_old(i) = st_sol%surf_head(i)
      end do
      !$omp end parallel do

      if (st_sim%res_type == 1) then
        st_time%current_t = st_init%rest_time
      end if

      st_time%now_time = st_time%current_t/st_sim%cal_fact

#ifdef MPI_MSG
      ! -- Write MPI 3D binary file (mpi_3dbin)
        call write_mpi_3dbin(st_out_fnum%head, ncalc, calc2calc, len_scal, st_sol%head_new,&
                             st_time%now_time)
#else
      ! -- Write header binary file (header_bin)
        call write_header_bin(st_out_fnum%head, st_time%now_time)
      ! -- Write 3D binary file (3dbin)
        call write_3dbin(st_out_fnum%head, ncalc, calc2calc, len_scal, st_sol%head_new)
#endif
      deallocate(calc2calc)

    else if (st_sim%sim_type >= 0) then
      if (st_time%conv_flag) then
        if (timestep_num == 0) then
          delt_old1 = st_time%delt
        end if
        timestep_num = timestep_num + 1
        delt_old1 = st_time%delt

        if (st_mpi%rank == 0) then
          ! -- Write boundary change information
            call write_bound_change(boundstep)
        end if
#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(boundstep, "boundary change number")
        end if
#endif
        if (boundstep >= 1) then
          st_time%inter_time = st_time%current_t - resi_time
          ! -- Set next date (date)
            call set_date(st_time%inter_time, st_time%now_date, resi_time)
          resi_time = st_time%current_t + resi_time
          if (trim(adjustl(st_sim%cal_unit)) == "YEA") then
            call conv_unit(st_mpi%rank, st_sim%cal_unit, "main file", st_time%now_date,&
                           st_sim%cal_fact)
          end if
          st_time%delt = st_time%delt
        else
          select case (st_ctrl%tstep_type)
          case (0)
            st_time%delt = st_time%delt*st_sim%inc_fact
          case (1)
            ! -- Apply heuristic time stepping (heuri)
              call apply_heuri(st_time%out_iter)
          end select
        end if

        ! -- Set value exchange (valexc)
          call set_valexc(ncalc, st_sol%head_new, st_sol%head_old)
          call set_valexc(ncalc, st_sol%srat_new, st_sol%srat_old)
          call set_valexc(ncals, st_sol%surf_head, st_sol%surf_old)
      else
        st_time%current_t = st_time%current_t - real(st_time%delt, kind=SP)
        st_time%delt = st_time%delt*st_sim%dec_fact
        ! -- Set value exchange (valexc)
          call set_valexc(ncalc, st_sol%head_old, st_sol%head_new)
          call set_valexc(ncalc, st_sol%srat_old, st_sol%srat_new)
          call set_valexc(ncals, st_sol%surf_old, st_sol%surf_head)
        if (st_sim%sim_type == 1) then
          ! -- Reset stepflag (stepf)
            call reset_stepf()
        end if
      end if
    else if (st_sim%sim_type == -1) then
      ! -- Set value exchange (valexc)
        call set_valexc(ncalc, st_sol%head_new, st_sol%head_old)
        call set_valexc(ncalc, st_sol%srat_new, st_sol%srat_old)
    end if

    next_time = st_time%current_t + real(st_time%delt, kind=SP)

  end subroutine calc_nextst

  subroutine set_nextet()
  !*********************************************************************************************
  ! set_nextet -- Set next end time
  !*********************************************************************************************
    ! -- modules
    use utility_module, only: write_logf
    use initial_module, only: in_type
    use read_module, only: read_next, read_intn
    use open_file, only: st_intse, st_intwe
#ifdef MPI_MSG
    use mpi_read, only: set_real4_fview
    use set_boundary, only: bfview, rfview, lfview
#endif
    ! -- inout

    ! -- local
    integer(I4) :: ierr
    character(:), allocatable :: bound_name, err_mes
    character(32) :: str_time
    !-------------------------------------------------------------------------------------------
    allocate(character(0) :: err_mes, bound_name)
    ierr = 0 ; str_time = "" ; err_mes = "" ; bound_name = ""
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%seal = 0
        st_seal%etime = st_sim%end_time
      else
        deallocate(st_forc%read_seal)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%rech = 0
        st_rech%etime = st_sim%end_time
        ! -- Reset cell value (value)
          call reset_value(st_bcnd%rech_num, st_forc%read_rech)
      else
        deallocate(st_forc%read_rech, st_bcnd%rech2cals, st_forc%calc_rech)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%well = 0
        st_well%etime = st_sim%end_time
        ! -- Reset cell value (value)
          call reset_value(st_bcnd%well_num, st_forc%read_well)
      else
        deallocate(st_forc%read_well, st_forc%calc_well, st_forc%well_top, st_forc%well_bott,&
                   st_bcnd%well_index, st_bcnd%well_conn)
        if (st_bcnd%well_num /= 0) then
          deallocate(st_hydr%abyd_well)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%prec = 0
        st_prec%etime = st_sim%end_time
        ! -- Reset cell value (value)
          call reset_value(st_bcnd%prec_num, st_forc%read_prec)
      else
        if (st_bcnd%prec_num > 0) then
          deallocate(st_forc%read_prec)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%evap = 0
        st_evap%etime = st_sim%end_time
        ! -- Reset cell value (value)
          call reset_value(st_bcnd%evap_num, st_forc%read_evap)
      else
        if (st_bcnd%evap_num > 0) then
          deallocate(st_forc%read_evap)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%riwl = 0
        st_riwl%etime = st_sim%end_time
      else
        deallocate(st_rive%cflag%wl, st_rive%calc%wl)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%ribl = 0
        st_ribl%etime = st_sim%end_time
      else
        deallocate(st_rive%cflag%bl, st_rive%calc%bl)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%ride = 0
        st_ride%etime = st_sim%end_time
      else
        if (st_ride%totn > 0) then
          deallocate(st_rive%cflag%de, st_rive%calc%de)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%riwi = 0
        st_riwi%etime = st_sim%end_time
      else
        if (st_riwi%totn > 0) then
          deallocate(st_rive%cflag%wi, st_rive%calc%wi)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%rile = 0
        st_rile%etime = st_sim%end_time
      else
        if (st_rile%totn > 0) then
          deallocate(st_rive%cflag%le, st_rive%calc%le)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%lawl = 0
        st_lawl%etime = st_sim%end_time
      else
        deallocate(st_lake%cflag%wl, st_lake%calc%wl)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%lawd = 0
        st_lawd%etime = st_sim%end_time
      else
        if (st_lawd%totn > 0) then
          deallocate(st_lake%cflag%wd, st_lake%calc%wd)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%labl = 0
        st_labl%etime = st_sim%end_time
      else
        deallocate(st_lake%cflag%bl, st_lake%calc%bl)
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
        if (st_mpi%rank == 0) then
          write(str_time,'(f0.3)') st_time%now_time
          err_mes = "Read final time step "//bound_name//" file at "//trim(str_time)//&
                    trim(st_sim%cal_unit)
          call write_logf(err_mes)
        end if
        st_step_flag%laar = 0
        st_laar%etime = st_sim%end_time
      else
        if (st_laar%totn > 0) then
          deallocate(st_lake%cflag%ar, st_lake%calc%ar)
        end if
      end if
    end if

    if (allocated(bound_name)) then
      deallocate(bound_name)
    end if

  end subroutine set_nextet

  subroutine set_delt()
  !*********************************************************************************************
  ! set_delt -- Set delta time
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    real(SP) :: min_step
    !-------------------------------------------------------------------------------------------
    min_step = min(st_rech%etime, st_well%etime, st_seal%etime, st_prec%etime,&
                   st_evap%etime, st_riwl%etime, st_riwd%etime, st_ribl%etime,&
                   st_ride%etime, st_riwi%etime, st_rile%etime, st_lawl%etime,&
                   st_lawd%etime, st_labl%etime, st_laar%etime)

    if (min_step == st_time%current_t) then
      st_time%delt = st_time%delt
    else if (min_step < next_time .and. next_time < st_sim%end_time) then
      st_time%delt = real(min_step - st_time%current_t, kind=DP)
    else if (next_time > st_sim%end_time) then
      st_time%delt = real(st_sim%end_time - st_time%current_t, kind=DP)
    end if

    if (st_time%delt > st_sim%max_step) then
      st_time%delt = real(st_sim%max_step, kind=DP)
    end if

    next_time = st_time%current_t + real(st_time%delt, kind=SP)

  end subroutine set_delt

  subroutine set_nextvar()
  !*********************************************************************************************
  ! set_nextvar -- Set next variable
  !*********************************************************************************************
    ! -- modules
    use set_condition, only: set_srabyd, set_chabyd, set_wellconn
    use assign_boundary, only: assign_sealv, assign_surfbv, assign_wellv, assign_rilav
    use calc_boundary, only: calc_reprev, conv_rech2calc, count_rivecalc, count_lakecalc,&
                             calc_wlbd, calc_rivea
#ifdef MPI_MSG
    use mpi_utility, only: mpisum_val
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4) :: prec_stepflag, evap_stepflag, rive_stepflag, lake_stepflag
    integer(I4) :: rive_aflag, sum_ribln, sum_riwln, sum_lawln, sum_labln
    real(DP), allocatable :: temp_area(:)
    !-------------------------------------------------------------------------------------------
    if (st_step_flag%seal == 1) then
      ! -- Assign sea level value (sealv)
        call assign_sealv(st_in_type%seal)

      st_step_flag%seal = 0

    else if (st_seal%etime == next_time) then
      st_step_flag%seal = 1
    end if

    if (st_step_flag%rech == 1) then
      allocate(st_bcnd%rech_cflag(ncals), st_forc%read_rech(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_bcnd%rech_cflag(i) = 0
        st_forc%read_rech(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign recharge value
        call assign_surfbv(st_in_type%rech, st_intre%type, st_rech, st_bcnd%rech_num,&
                           st_bcnd%rech_cflag,&
                           st_forc%read_rech)

      call conv_rech2calc(st_bcnd%rech_num)

      st_step_flag%rech = 0
    else if (st_rech%etime == next_time) then
      st_step_flag%rech = 1
    end if

    if (st_step_flag%well == 1) then
      ! -- Assign well value (wellv)
        call assign_wellv(st_in_type%well, st_in_type%weks, st_in_type%weke, st_bcnd%well_num)

      st_step_flag%well = 0
      if (st_bcnd%well_num /= 0) then
        ! -- Set well connectivity (wellconn)
          call set_wellconn(st_bcnd%well_num, st_hydr%read_hydx, st_hydr%read_hydy)
      end if

    else if (st_well%etime == next_time) then
      st_step_flag%well = 1
    end if

    prec_stepflag = 0

    if (st_step_flag%prec == 1) then
      allocate(st_bcnd%prec_cflag(ncals), st_forc%read_prec(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_bcnd%prec_cflag(i) = 0
        st_forc%read_prec(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign precipitation value
        call assign_surfbv(st_in_type%prec, st_intpr%type, st_prec, st_bcnd%prec_num,&
                           st_bcnd%prec_cflag,&
                           st_forc%read_prec)
      prec_stepflag = prec_stepflag + 1
      deallocate(st_bcnd%prec_cflag)

      st_step_flag%prec = 0
    else if (st_prec%etime == next_time) then
      st_step_flag%prec = 1
    end if

    evap_stepflag = 0

    if (st_step_flag%evap == 1) then
      allocate(st_bcnd%evap_cflag(ncals), st_forc%read_evap(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_bcnd%evap_cflag(i) = 0
        st_forc%read_evap(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign evapotranspiration value
        call assign_surfbv(st_in_type%evap, st_intev%type, st_evap, st_bcnd%evap_num,&
                           st_bcnd%evap_cflag,&
                           st_forc%read_evap)
      evap_stepflag = evap_stepflag + 1
      deallocate(st_bcnd%evap_cflag)

      st_step_flag%evap = 0
    else if (st_evap%etime == next_time) then
      st_step_flag%evap = 1
    end if

    if (prec_stepflag /= 0 .or. evap_stepflag /= 0) then
      deallocate(st_forc%read_rech, st_bcnd%rech2cals, st_forc%calc_rech)
      ! -- Calculate recharge from precipitation and evapotranspiration (reprev)
        call calc_reprev(st_bcnd%rech_num)
      call conv_rech2calc(st_bcnd%rech_num)
    end if

    rive_stepflag = 0 ; rive_aflag = 0

    if (st_step_flag%riwl == 1) then
      allocate(st_rive%cflag%wl(ncals), st_rive%calc%wl(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_rive%cflag%wl(i) = 0
        st_rive%calc%wl(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign river water level value
        call assign_rilav(st_rivf_type%wlev, 0, st_riwl, st_rive%num%wl, st_rive%cflag%wl,&
                          st_rive%calc%wl)

      st_step_flag%riwl = 0 ; rive_stepflag = rive_stepflag + 1
    else if (st_riwl%etime == next_time) then
      st_step_flag%riwl = 1
    end if

    if (st_step_flag%ribl == 1) then
      allocate(st_rive%cflag%bl(ncals), st_rive%calc%bl(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_rive%cflag%bl(i) = 0
        st_rive%calc%bl(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign river bottom level value
        call assign_rilav(st_rivf_type%blev, 0, st_ribl, st_rive%num%bl, st_rive%cflag%bl,&
                          st_rive%calc%bl)

      st_step_flag%ribl = 0 ; rive_stepflag = rive_stepflag + 1
    else if (st_ribl%etime == next_time) then
      st_step_flag%ribl = 1
    end if

    if (st_step_flag%riwd == 1) then
      if (st_riwd%totn > 0) then
        allocate(st_rive%cflag%wd(ncals), st_rive%calc%wd(ncals))
        !$omp parallel do private(i)
        do i = 1, ncals
          st_rive%cflag%wd(i) = 0
          st_rive%calc%wd(i) = SNOVAL
        end do
        !$omp end parallel do
      end if
      ! -- Assign river water depth value
        call assign_rilav(st_rivf_type%wdep, 0, st_riwd, st_rive%num%wd, st_rive%cflag%wd,&
                          st_rive%calc%wd)

#ifdef MPI_MSG
      ! -- Sum value for MPI (val)
        call mpisum_val(st_rive%num%wl, "river water level", sum_riwln)
        call mpisum_val(st_rive%num%bl, "river bottom level", sum_ribln)
#else
      sum_riwln = st_rive%num%wl ; sum_ribln = st_rive%num%bl
#endif

      if (sum_riwln == 0 .and. sum_ribln /= 0) then
        ! -- Calculate water level from river bottom level (wlrb)
          call calc_wlbd(st_rive%cflag%bl, st_rive%calc%bl, st_rive%cflag%wd, st_rive%calc%wd,&
                         st_rive%cflag%wl, st_rive%calc%wl, st_rive%num%wl)
        st_rive%num%wl = 0
        deallocate(st_rive%cflag%wd, st_rive%calc%wd)
      end if

      st_step_flag%riwd = 0 ; rive_stepflag = rive_stepflag + 1
    else if (st_riwd%etime == next_time) then
      st_step_flag%riwd = 1
    end if

    if (st_step_flag%riwi == 1) then
      if (st_riwi%totn > 0) then
        allocate(st_rive%cflag%wi(ncals), st_rive%calc%wi(ncals))
        !$omp parallel do private(i)
        do i = 1, ncals
          st_rive%cflag%wi(i) = 0
          st_rive%calc%wi(i) = SNOVAL
        end do
        !$omp end parallel do
      end if
      ! -- Assign river width value
        call assign_rilav(st_rivf_type%widt, 0, st_riwi, st_rive%num%wi, st_rive%cflag%wi,&
                          st_rive%calc%wi)

      st_step_flag%riwi = 0 ; rive_aflag = rive_aflag + 1
    else if (st_riwi%etime == next_time) then
      st_step_flag%riwi = 1
    end if

    if (st_step_flag%rile == 1) then
      if (st_rile%totn > 0) then
        allocate(st_rive%cflag%le(ncals), st_rive%calc%le(ncals))
        !$omp parallel do private(i)
        do i = 1, ncals
          st_rive%cflag%le(i) = 0
          st_rive%calc%le(i) = SNOVAL
        end do
        !$omp end parallel do
      end if
      ! -- Assign river length value
        call assign_rilav(st_rivf_type%leng, 0, st_rile, st_rive%num%le, st_rive%cflag%le,&
                          st_rive%calc%le)

      st_step_flag%rile = 0 ; rive_aflag = rive_aflag + 1
    else if (st_rile%etime == next_time) then
      st_step_flag%rile = 1
    end if

    if (rive_aflag > 0) then
      !$omp parallel do private(i)
      do i = 1, ncals
        st_rive%cflag%ar(i) = 0
        st_rive%calc%ar(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Calculate river area (rivea)
        call calc_rivea(st_rive%cflag%wi, st_rive%cflag%le, st_rive%calc%wi, st_rive%calc%le,&
                        st_rive%cflag%ar, st_rive%calc%ar, st_rive%num%ar)
    end if

    if (rive_stepflag > 0) then
      if (st_bcnd%rive_num /= 0) then
        allocate(temp_area(st_bcnd%rive_num))
        !$omp parallel do private(i)
        do i = 1, st_bcnd%rive_num
          temp_area(i) = -st_forc%rive_area(i)
        end do
        !$omp end parallel do
        ! -- Set surface&recharge area and area by distance (srabyd)
          call set_srabyd(st_bcnd%rive_num, st_forc%rive_bott, temp_area,&
                          st_bcnd%rive2cals, st_forc%abyd_rive)
        deallocate(st_forc%rive_head, st_forc%rive_bott, st_forc%rive_area, temp_area)
      end if
      deallocate(st_bcnd%rive2cals)
      ! -- Count river calculation (rivecalc)
        call count_rivecalc(st_rive%cflag%wl, st_rive%cflag%bl, st_rive%cflag%ar,&
                            st_rive%calc%wl, st_rive%calc%bl, st_rive%calc%ar,&
                            st_bcnd%rive_num)
      if (st_bcnd%rive_num /= 0) then
        deallocate(st_forc%abyd_rive)
        allocate(st_forc%abyd_rive(st_bcnd%rive_num))
        !$omp parallel do private(i)
        do i = 1, st_bcnd%rive_num
          st_forc%abyd_rive(i) = DZERO
        end do
        !$omp end parallel do
        ! -- Set surface&recharge area and area by distance (srabyd)
          call set_srabyd(st_bcnd%rive_num, st_forc%rive_bott, st_forc%rive_area,&
                        st_bcnd%rive2cals, st_forc%abyd_rive)
      end if
    end if

    lake_stepflag = 0

    if (st_step_flag%lawl == 1) then
      allocate(st_lake%cflag%wl(ncals), st_lake%calc%wl(ncals))
      st_lake%cflag%wl(:) = 0 ; st_lake%calc%wl(:) = SNOVAL
      ! -- Assign lake water level value
        call assign_rilav(st_lakf_type%wlev, 0, st_lawl, st_lake%num%wl, st_lake%cflag%wl,&
                          st_lake%calc%wl)

      st_step_flag%lawl = 0 ; lake_stepflag = lake_stepflag + 1
    else if (st_lawl%etime == next_time) then
      st_step_flag%lawl = 1
    end if

    if (st_step_flag%labl == 1) then
      allocate(st_lake%cflag%bl(ncals), st_lake%calc%bl(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_lake%cflag%bl(i) = 0
        st_lake%calc%bl(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign lake bottom level value
        call assign_rilav(st_lakf_type%blev, 0, st_labl, st_lake%num%bl, st_lake%cflag%bl,&
                          st_lake%calc%bl)

      st_step_flag%labl = 0 ; lake_stepflag = lake_stepflag + 1
    else if (st_labl%etime == next_time) then
      st_step_flag%labl = 1
    end if

    if (st_step_flag%lawd == 1) then
      if (st_lawd%totn > 0) then
        allocate(st_lake%cflag%wd(ncals), st_lake%calc%wd(ncals))
        !$omp parallel do private(i)
        do i = 1, ncals
          st_lake%cflag%wd(i) = 0
          st_lake%calc%wd(i) = SNOVAL
        end do
        !$omp end parallel do
      end if
      ! -- Assign lake water depth value
        call assign_rilav(st_lakf_type%wdep, 0, st_lawd, st_lake%num%wd, st_lake%cflag%wd,&
                          st_lake%calc%wd)

#ifdef MPI_MSG
      ! -- Sum value for MPI (val)
        call mpisum_val(st_lake%num%wl, "lake water level", sum_lawln)
        call mpisum_val(st_lake%num%bl, "lake bottom level", sum_labln)
#else
      sum_lawln = st_lake%num%wl ; sum_labln = st_lake%num%bl
#endif

      if (sum_lawln == 0 .and. sum_labln /= 0) then
        ! -- Calculate water level from bottom level and water depth (wlbd)
          call calc_wlbd(st_lake%cflag%bl, st_lake%calc%bl, st_lake%cflag%wd, st_lake%calc%wd,&
                         st_lake%cflag%wl, st_lake%calc%wl, st_lake%num%wl)
        st_lake%num%wl = 0
        deallocate(st_lake%cflag%wd, st_lake%calc%wd)
      end if

      st_step_flag%lawd = 0 ; lake_stepflag = lake_stepflag + 1
    else if (st_lawd%etime == next_time) then
      st_step_flag%lawd = 1
    end if

    if (st_step_flag%laar == 1) then
      if (st_laar%totn > 0) then
        allocate(st_lake%cflag%ar(ncals), st_lake%calc%ar(ncals))
        !$omp parallel do private(i)
        do i = 1, ncals
          st_lake%cflag%ar(i) = 0
          st_lake%calc%ar(i) = SNOVAL
        end do
        !$omp end parallel do
      end if
      ! -- Assign lake area value
        call assign_rilav(st_lakf_type%area, 1, st_laar, st_lake%num%ar, st_lake%cflag%ar,&
                          st_lake%calc%ar)

      st_step_flag%laar = 0 ; lake_stepflag = lake_stepflag + 1
    else if (st_laar%etime == next_time) then
      st_step_flag%laar = 1
    end if

    if (lake_stepflag > 0) then
      if (st_bcnd%lake_num /= 0) then
        allocate(temp_area(st_bcnd%lake_num))
        !$omp parallel do private(i)
        do i = 1, st_bcnd%lake_num
          temp_area(i) = -st_forc%lake_area(i)
        end do
        !$omp end parallel do
        ! -- Set surface&recharge area and area by distance (srabyd)
          call set_srabyd(st_bcnd%lake_num, st_forc%lake_bott, temp_area,&
                          st_bcnd%lake2cals, st_forc%abyd_lake)
        deallocate(st_forc%lake_head, st_forc%lake_bott, st_forc%lake_area, temp_area)
      end if
      deallocate(st_bcnd%lake2cals)
      ! -- Count lake calculation cell (lakecalc)
        call count_lakecalc(st_lake%cflag%wl, st_lake%cflag%bl, st_lake%cflag%ar,&
                            st_lake%calc%wl, st_lake%calc%bl, st_lake%calc%ar,&
                            st_bcnd%lake_num)
      if (st_bcnd%lake_num /= 0) then
        deallocate(st_forc%abyd_lake)
        allocate(st_forc%abyd_lake(st_bcnd%lake_num))
        !$omp parallel do private(i)
        do i = 1, st_bcnd%lake_num
          st_forc%abyd_lake(i) = DZERO
        end do
        !$omp end parallel do
        ! -- Set surface&recharge area and area by distance (srabyd)
          call set_srabyd(st_bcnd%lake_num, st_forc%lake_bott, st_forc%lake_area,&
                        st_bcnd%lake2cals, st_forc%abyd_lake)
      end if
    end if

    if (rive_stepflag > 0 .or. lake_stepflag > 0) then
      ! -- Set charge area by distance (chabyd)
        call set_chabyd()
    end if

  end subroutine set_nextvar

  subroutine set_vwell_head(st_sol)
  !*********************************************************************************************
  ! set_vwell_head -- Set virtual well head
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    ! -- Calculate virtual well head without well pumping (vheadout)
      call calc_vheadout(st_sol)

    !$omp parallel do private(i)
    do i = 1, ncalc
      st_forc%calc_well(i) = DZERO
    end do
    !$omp end parallel do

    ! -- Calculate virtual well head with well puming (vheadin)
      call calc_vheadin(st_sol)

  end subroutine set_vwell_head

  subroutine change_recharge(st_sol)
  !*********************************************************************************************
  ! change_recharge -- Change the recharge volume
  !*********************************************************************************************
    ! -- modules
    ! -- inout

    type(sol_set), intent(in) :: st_sol
    ! -- local
    integer(I4) :: i, s
    real(DP) :: norm_elev
    real(DP), allocatable :: water_dep(:), rech_rati(:)
    !-------------------------------------------------------------------------------------------
    allocate(water_dep(ncals), rech_rati(st_bcnd%rech_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncals
      water_dep(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, st_bcnd%rech_num
      rech_rati(i) = DONE
    end do
    !$omp end do

    !$omp do private(i, s, norm_elev)
    do i = 1, st_bcnd%rech_num
      s = st_bcnd%rech2cals(i)
      water_dep(s) = st_sol%head_new(s) - st_hydr%surf_bott(s)
      if (st_forc%read_rech(i) < SZERO .or. water_dep(s) < DZERO) then
        rech_rati(i) = DONE
      else if (water_dep(s) < st_hydr%surf_reli(s) .and. st_hydr%surf_reli(s) /= DZERO) then
        norm_elev = water_dep(s)/st_hydr%surf_reli(s)
        rech_rati(i) = norm_elev**st_hydr%surf_parm(s)
        rech_rati(i) = DONE - rech_rati(i)
      else if (water_dep(s) >= st_hydr%surf_reli(s) .and. st_hydr%surf_reli(s) /= DZERO) then
        rech_rati(i) = DZERO
      else
        rech_rati(i) = DONE
      end if
      st_forc%calc_rech(i) = st_forc%read_rech(i)*st_hydr%rech_area(s)*rech_rati(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(water_dep, rech_rati)

  end subroutine change_recharge

  subroutine calc_vheadout(st_sol)
  !*********************************************************************************************
  ! calc_vheadout -- Calculate virtual well head without well pumping
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    type(sol_set), intent(in) :: st_sol
    ! -- local
    integer(I4) :: i, j, k
    real(DP) :: tot_cond, tot_flux
    !-------------------------------------------------------------------------------------------
    allocate(whead_new(st_bcnd%well_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%well_num
      whead_new(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i, j, k, tot_cond, tot_flux)
    do i = 1, st_bcnd%well_num
      tot_cond = DZERO ; tot_flux = DZERO
      do k = st_bcnd%well_index(i-1)+1, st_bcnd%well_index(i)
        j = st_bcnd%well_conn(k)
        tot_cond = tot_cond + st_sol%rel_perm(j)*st_hydr%abyd_well(k)
        tot_flux = tot_flux + st_sol%rel_perm(j)*st_hydr%abyd_well(k)*st_sol%head_new(j)
      end do
      if (tot_cond /= DZERO) then
        whead_new(i) = tot_flux/tot_cond
      end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine calc_vheadout

  subroutine calc_vheadin(st_sol)
  !*********************************************************************************************
  ! calc_vheadin -- Calculate virtual well head with well pumping
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    type(sol_set), intent(in) :: st_sol
    ! -- local
    integer(I4) :: i, j, k
    real(DP) :: tot_cond
    real(DP), allocatable :: temp_whead(:)
    !-------------------------------------------------------------------------------------------
    allocate(temp_whead(st_bcnd%well_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, st_bcnd%well_num
      temp_whead(i) = whead_new(i)
    end do
    !$omp end do

    !$omp do private(i, j, k, tot_cond)
    do i = 1, st_bcnd%well_num
      tot_cond = DZERO
      do k = st_bcnd%well_index(i-1)+1, st_bcnd%well_index(i)
        j = st_bcnd%well_conn(k)
        tot_cond = tot_cond + st_sol%rel_perm(j)*st_hydr%abyd_well(k)
      end do
      if (tot_cond /= DZERO) then
        whead_new(i) = temp_whead(i) + st_forc%read_well(i)/tot_cond

!        if (whead_new(i) > well_top(i)) then
!          whead_new(i) = well_top(i)
!        else if (whead_new(i) < well_bott(i)) then
!          whead_new(i) = well_bott(i)
!        end if

        do k = st_bcnd%well_index(i-1)+1, st_bcnd%well_index(i)
          j = st_bcnd%well_conn(k)
          if (whead_new(i) > st_forc%well_bott(i) .or. st_sol%head_new(j) > st_forc%well_bott(i)) then
            st_forc%calc_well(j) = st_forc%calc_well(j) + st_sol%rel_perm(j)*st_hydr%abyd_well(k)*&
                           (whead_new(i)-st_sol%head_new(j))
          end if
        end do
      end if
    end do
    !$omp end do
    !$omp end parallel

    deallocate(whead_new, temp_whead)

  end subroutine calc_vheadin

  subroutine set_valexc(excn, exc_inval, exc_outval)
  !*********************************************************************************************
  ! set_valexc -- Set value exchange
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: excn
    real(DP), intent(in) :: exc_inval(:)
    real(DP), intent(out) :: exc_outval(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, excn
      exc_outval(i) = exc_inval(i)
    end do
    !$omp end parallel do

  end subroutine set_valexc

  subroutine write_bound_change(bchange)
  !*********************************************************************************************
  ! write_bound_change -- Write boundary change information
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(out) :: bchange
    ! -- local
    integer(I4) :: conv_fnum
    character(9) :: cond_format
    !-------------------------------------------------------------------------------------------
    bchange = st_step_flag%rech + st_step_flag%well + st_step_flag%seal +&
              st_step_flag%prec + st_step_flag%evap + st_step_flag%riwl +&
              st_step_flag%riwd + st_step_flag%ribl + st_step_flag%ride +&
              st_step_flag%riwi + st_step_flag%lawl + st_step_flag%lawd +&
              st_step_flag%labl + st_step_flag%laar

    conv_fnum = st_out_fnum%conv ; cond_format = "(a,f10.3)"

    !$omp parallel shared(conv_fnum, cond_format)
    !$omp single
    if (st_step_flag%rech == 1) then
      write(conv_fnum,cond_format) "Changed recharge condition at ", st_time%now_time
    end if
    if (st_step_flag%well == 1) then
      write(conv_fnum,cond_format) "Changed well condition at ", st_time%now_time
    end if
    if (st_step_flag%seal == 1) then
      write(conv_fnum,cond_format) "Changed sea condition at ", st_time%now_time
    end if
    if (st_step_flag%prec == 1) then
      write(conv_fnum,cond_format) "Changed precipitation condition at ", st_time%now_time
    end if
    if (st_step_flag%evap == 1) then
      write(conv_fnum,cond_format) "Changed evapotranspiration condition at ", st_time%now_time
    end if
    if (st_step_flag%riwl == 1) then
      write(conv_fnum,cond_format) "Changed river water level at ", st_time%now_time
    end if
    if (st_step_flag%riwd == 1) then
      write(conv_fnum,cond_format) "Changed river water depth at ", st_time%now_time
    end if
    if (st_step_flag%ribl == 1) then
      write(conv_fnum,cond_format) "Changed river bottom level at ", st_time%now_time
    end if
    if (st_step_flag%ride == 1) then
      write(conv_fnum,cond_format) "Changed river depth at ", st_time%now_time
    end if
    if (st_step_flag%riwi == 1) then
      write(conv_fnum,cond_format) "Changed river width at ", st_time%now_time
    end if
    if (st_step_flag%lawl == 1) then
      write(conv_fnum,cond_format) "Changed lake water level at ", st_time%now_time
    end if
    if (st_step_flag%lawd == 1) then
      write(conv_fnum,cond_format) "Changed lake water depth at ", st_time%now_time
    end if
    if (st_step_flag%labl == 1) then
      write(conv_fnum,cond_format) "Changed lake bottom level at ", st_time%now_time
    end if
    if (st_step_flag%laar == 1) then
      write(conv_fnum,cond_format) "Changed lake area at ", st_time%now_time
    end if
    !$omp end single
    !$omp end parallel

  end subroutine write_bound_change

  subroutine reset_stepf()
  !*********************************************************************************************
  ! reset_stepf -- Reset stepflag
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------------
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
  !*********************************************************************************************
  ! reset_value -- Reset cell value
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: targn
    real(SP), intent(out) :: targ_value(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, targn
      targ_value(i) = SZERO
    end do
    !$omp end parallel do

  end subroutine reset_value

  subroutine apply_heuri(out_num)
  !*********************************************************************************************
  ! apply_heuri -- Apply heuristic time stepping
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: out_num
    ! -- local
    integer(I4) :: incr_num, decr_num
    !-------------------------------------------------------------------------------------------
    incr_num = int(st_ctrl%maxout_iter*0.4) ; decr_num = int(st_ctrl%maxout_iter*0.8)
    if (out_num <= incr_num) then
      st_time%delt = delt_old1*st_sim%inc_fact
    else if (out_num <= decr_num) then
      st_time%delt = delt_old1
    else
      st_time%delt = delt_old1*st_sim%dec_fact
    end if

  end subroutine apply_heuri

  subroutine set_date(inttime, ndate, restime)
  !*********************************************************************************************
  ! set_date -- Set next date
  !*********************************************************************************************
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
    !-------------------------------------------------------------------------------------------
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
