program msg123
!*****************************************************************************************
! Main program of Multi-scale groundwater 1phase 2resolution 3dimensional model (Msg123)
!*****************************************************************************************
  ! -- modules
  use kind_module, only: I4, DP
  use utility_module, only: log_fnum
  use initial_module, only: my_rank, init_msg
  use read_input, only: read_main_file
  use set_cell, only: set_cell_info
  use prep_calculation, only: prepare_calc, conv_flag
  use set_boundary, only: set_bound
  use allocate_solution, only: allocate_solvar
  use time_module, only: update_tstep, now_time
  use nonlinear_solution, only: calc_numsol
  use write_output, only: write_outf
  use check_simulation, only: check_lastts, lasttime_flag
#ifdef ICI
  use ici_module, only: set_mapt, put_initv, get_var, alloc_outvar, put_var, fin_ici
#endif
#ifdef MPI_MSG
  use mpi_initfin, only: fin_mpi
#endif

  implicit none

  ! -- local
  integer(I4) :: i
  integer(I4) :: sta_value(8), end_value(8)
  real(DP) :: tot_stime, tot_etime, loop_stime, loop_etime
  ! -- format
  11 format(/"Run end date and time(yyyy/mm/dd hh:mm:ss) : ",i4,"/",i2.2,"/",i2.2,1x,i2,":",i2.2,":",i2.2,/)
  12 format(/"Total cpu time : ", es15.6, " (sec)")
  13 format(/"Time loop cpu time : ", es15.6, " (sec)")
  !---------------------------------------------------------------------------------------
  if (my_rank == 0) then
    ! -- Start time
      call DATE_AND_TIME(values = sta_value)
    ! -- Calculation start time
      call CPU_TIME(tot_stime)
  end if

  ! -- Initialize msg (msg)
    call init_msg(sta_value)

  if (my_rank == 0) then
    ! -- Read main files (main_file)
      call read_main_file()
  end if

  ! -- Set cell information (cell_info)
    call set_cell_info()

  ! -- Prepare calculation (calc)
    call prepare_calc()

#ifdef ICI
  ! -- Set mapping table (mapt)
    call set_mapt()
  ! -- Put initial value (initv)
    call put_initv()
  ! -- Get variables (var)
    call get_var()
#endif

  ! -- Set boundary (bound)
    call set_bound()

#ifdef ICI
  ! -- Allocate output variables (outvar)
    call alloc_outvar()
#endif

  ! -- Allocate solution variable for time step (solvar)
    call allocate_solvar()

  if (my_rank == 0) then
    ! -- Time loop start time
      call CPU_TIME(loop_stime)
  end if

  ! start time step loop
  tstep_loop: do
    ! -- Update time step (tstep)
      call update_tstep()

    ! -- Calculate numerical solution (numsol)
      call calc_numsol()

    if (conv_flag == 1) then
      ! -- Check last time step conditions (lastts)
        call check_lastts()
      ! -- Write output file (outf)
        call write_outf(now_time)
#ifdef ICI
      ! -- Put variables (var)
        call put_var()
#endif
      if(lasttime_flag == 1) then
        exit tstep_loop
      end if
    end if

  end do tstep_loop

#ifdef ICI
  ! -- Finalize ici (ici)
    call fin_ici()
#endif

#ifdef MPI_MSG
  ! -- Finalize mpi (mpi)
    call fin_mpi(my_rank, log_fnum)
#endif

  if (my_rank == 0) then
    ! -- Time loop end time
      call CPU_TIME(loop_etime)
    ! -- Calculation end time
      call CPU_TIME(tot_etime)

    write(log_fnum,13) loop_etime - loop_stime
    write(log_fnum,12) tot_etime - tot_stime

    ! -- End time
      call DATE_AND_TIME(values = end_value)

    write(log_fnum,11) (end_value(i), i = 1, 3), (end_value(i), i = 5, 7)
  end if

end program msg123
