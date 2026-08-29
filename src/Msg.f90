program msg123
!***********************************************************************************************
! Main program of Multi-scale groundwater 1phase 2resolution 3dimensional model (Msg123)
!***********************************************************************************************
  ! -- modules
  use kind_module, only: I4, DP
  use types_module, only: kryl_set, amgt_set, coef_set, sol_set
  use utility_module, only: log_fnum, st_mpi, dilu_shift_num, slope_sign_num
  use utility_module, only: nan_recv_num, maxstep_num, satlim_num
  use initial_module, only: init_msg, st_ctrl
  use read_input, only: read_main_file
  use set_cell, only: set_cell_info
  use prep_calculation, only: prepare_calc, st_time
  use set_boundary, only: set_bound
  use allocate_solution, only: allocate_solvar
  use calc_function, only: allocate_calfun
  use make_linearsystem, only: allocate_matvec
  use check_simulation, only: check_lastts, lasttime_flag
  use linear_solution, only: allocate_amgalg, allocate_krylov
  use time_module, only: update_tstep
  use nonlinear_solution, only: allocate_nonlin, calc_numsol
  use write_output, only: write_outf
#ifdef ICI
  use ici_module, only: set_mapt, put_initv, get_var, alloc_outvar, put_var, fin_ici
#endif
#ifdef MPI_MSG
  use mpi_initfin, only: fin_mpi
  use mpi_utility, only: mpisum_val
  use mpi_solve, only: allocate_mpisolve
#endif

  implicit none

  ! -- local
  integer(I4) :: i
  integer(I4) :: sta_value(8), end_value(8)
  real(DP) :: tot_stime, tot_etime, loop_stime, loop_etime
  type(sol_set) :: st_sol
  type(kryl_set) :: st_kryl
  type(amgt_set) :: st_amgt
  type(coef_set) :: st_coef
#ifdef MPI_MSG
  integer(I4) :: sum_shift, sum_slope, sum_nanr, sum_maxs, sum_satl
#endif
  ! -- format
  11 format(/"Run end date and time(yyyy/mm/dd hh:mm:ss) : ",i4,"/",i2.2,"/",i2.2,1x,i2,":",&
            i2.2,":",i2.2,/)
  12 format(/"Total cpu time : ", es15.6, " (sec)")
  13 format(/"Time loop cpu time : ", es15.6, " (sec)")
  14 format(/a," : ", i12, " times")
  !--------------------------------------------------------------------------------------------
  if (st_mpi%rank == 0) then
    ! -- Start time
      call DATE_AND_TIME(values = sta_value)
    ! -- Calculation start time
      call CPU_TIME(tot_stime)
  end if

  ! -- Initialize msg (msg)
    call init_msg(sta_value)

  if (st_mpi%rank == 0) then
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
    call allocate_solvar(st_sol)
  ! -- Allocate for calculate function value (calfun)
    call allocate_calfun()
  ! -- Allocate for matrix and vector (matvec)
    call allocate_matvec(st_coef)
  ! -- Allocate for nonlinear solution work arrays (nonlin)
    call allocate_nonlin()
  ! -- Allocate for krylov work vectors (krylov)
    call allocate_krylov(st_kryl)
#ifdef MPI_MSG
  if (st_mpi%totn /= 1) then
    ! -- Allocate for mpi solve work buffers (mpisolve)
      call allocate_mpisolve()
  end if
#endif
  if (st_ctrl%precon_type == 1) then
    ! -- Allocate for amg algebra (amgalg)
      call allocate_amgalg(st_amgt)
  end if

  if (st_mpi%rank == 0) then
    ! -- Time loop start time
      call CPU_TIME(loop_stime)
  end if

  ! start time step loop
  tstep_loop: do
    ! -- Update time step (tstep)
      call update_tstep(st_sol)

    ! -- Calculate numerical solution (numsol)
      call calc_numsol(st_kryl, st_amgt, st_coef, st_sol)

    if (st_time%conv_flag) then
      ! -- Check last time step conditions (lastts)
        call check_lastts()
      ! -- Write output file (outf)
        call write_outf(st_time%now_time, st_sol)
#ifdef ICI
      ! -- Put variables (var)
        call put_var(st_sol)
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
  if (st_mpi%totn /= 1) then
    ! -- Sum value for MPI (val)
      call mpisum_val(dilu_shift_num, "dilu pivot shift count", sum_shift)
      call mpisum_val(slope_sign_num, "non-negative slope count", sum_slope)
      call mpisum_val(nan_recv_num, "nan step recover count", sum_nanr)
      call mpisum_val(maxstep_num, "maximum step taken count", sum_maxs)
      call mpisum_val(satlim_num, "saturation limit count", sum_satl)
    dilu_shift_num = sum_shift ; slope_sign_num = sum_slope
    nan_recv_num = sum_nanr ; maxstep_num = sum_maxs ; satlim_num = sum_satl
  end if
  ! -- Finalize mpi (mpi)
    call fin_mpi(st_mpi%rank, log_fnum)
#endif

  if (st_mpi%rank == 0) then
    if (dilu_shift_num /= 0) then
      write(log_fnum,14) "Dilu pivot shifted", dilu_shift_num
    end if
    if (slope_sign_num /= 0) then
      write(log_fnum,14) "Non-negative slope", slope_sign_num
    end if
    if (nan_recv_num /= 0) then
      write(log_fnum,14) "Nan step recovered", nan_recv_num
    end if
    if (maxstep_num /= 0) then
      write(log_fnum,14) "Maximum step taken", maxstep_num
    end if
    if (satlim_num /= 0) then
      write(log_fnum,14) "Sat limit applied ", satlim_num
    end if

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
