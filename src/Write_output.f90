module write_output
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SZERO, DZERO
  use utility_module, only: close_file
  use initial_module, only: my_rank
  use read_input, only: len_scal
  use check_condition, only: st_out_fnum
  use set_cell, only: ncalc, ncals
  use set_condition, only: rech_area
  use assign_calc, only: msout_tnum
  use allocate_solution, only: head_new, srat_new
  use check_simulation, only: lasttime_flag
  use write_module, only: write_header_bin, write_2dbin, write_3dbin
#ifdef MPI_MSG
  use initial_module, only: pro_totn
  use mpi_read, only: close_mpi_file
  use mpi_write, only: write_mpi_2dbin, write_mpi_3dbin
#endif

  implicit none
  private
  public :: write_outf

  ! -- local
  character(:), allocatable :: msformat

  contains

  subroutine write_outf(time_val)
  !***************************************************************************************
  ! write_outf -- write output file
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_out_step, st_out_type, out_type
    use prep_calculation, only: current_t
    use check_simulation, only: check_outtiming, write_flag
    use allocate_output, only: allocate_outvar
    use calc_output, only: calc_cell_mas, calc_rivr_off, calc_lakr_off, calc_sufr_off,&
                           calc_dunr_off, calc_seal_res, calc_rech_res, calc_well_res
#ifdef MPI_MSG
    use mpi_write, only: set_senrec_wtab, write_mpi_rest
#endif
    ! -- inout
    real(SP), intent(in) :: time_val
    ! -- local
    integer(I4) :: rest_fnum
    integer(I4) :: header_flag = 0, allocate_flag = 0
    !-------------------------------------------------------------------------------------
    ! -- Check output timing (outtiming)
      call check_outtiming()

    if (allocate_flag == 0) then
      allocate_flag = 1
      ! -- Allocate output variable (outvar)
        call allocate_outvar()
    end if

    if (header_flag == 0 .and. write_flag == 1) then
      header_flag = 1
      if (my_rank == 0 .and. st_out_type%mass == out_type(1)) then
        ! -- Write massbalance file header (mass_header)
          call write_mass_header()
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1 .and. st_out_type%wtab == out_type(2)) then
        ! -- Set send and receive for water table (senrec_wtab)
          call set_senrec_wtab()
      end if
#endif
    end if

    if (st_out_type%mass == out_type(1)) then
      ! -- Calculate cell massbalance (cell_mas)
        call calc_cell_mas()
    end if

    if (st_out_type%rivr == out_type(2)) then
      ! -- Calculate river runoff (rivr_off)
        call calc_rivr_off()
    end if

    if (st_out_type%lakr == out_type(2)) then
      ! -- Calculate lake runoff (lakr_off)
        call calc_lakr_off()
    end if

    if (st_out_type%sufr == out_type(2)) then
      ! -- Calculate surface runoff (sufr_off)
        call calc_sufr_off()
    end if

    if (st_out_type%dunr == out_type(2)) then
      ! -- Calculate dunne runoff (dunr_off)
        call calc_dunr_off()
    end if

    if (st_out_type%seal == out_type(3)) then
      ! -- Calculate sea results (seal_res)
        call calc_seal_res()
    end if

    if (st_out_type%rech == out_type(2)) then
      ! -- Calculate recharge results (rech_res)
        call calc_rech_res()
    end if

    if (st_out_type%well == out_type(3)) then
      ! -- Calculate well results (well_res)
        call calc_well_res()
    end if

    if (write_flag == 1) then

      if (st_out_step%head == SZERO) then
        ! -- Write output head file (out_headf)
          call write_out_headf(time_val)
      else if (mod(current_t,st_out_step%head) == 0) then
        ! -- Write output head file (out_headf)
          call write_out_headf(time_val)
      else if (lasttime_flag == 1) then
        ! -- Write output head file (out_headf)
          call write_out_headf(time_val)
      end if

      if (lasttime_flag == 1) then
        rest_fnum = st_out_fnum%rest
#ifdef MPI_MSG
        ! -- Write mpi restart file (mpi_rest)
          call write_mpi_rest(rest_fnum, time_val, len_scal, head_new)
          call close_mpi_file(rest_fnum)
#else
        ! -- Write restart file (out_restf)
          call write_out_restf(rest_fnum, time_val)
#endif
      end if

      if (st_out_type%srat == out_type(3)) then
        if (st_out_step%srat == SZERO) then
          ! -- Write output saturation file (out_sratf)
            call write_out_sratf(time_val)
        else if (mod(current_t,st_out_step%srat) == 0) then
          ! -- Write output saturation file (out_sratf)
            call write_out_sratf(time_val)
        else if (lasttime_flag == 1) then
          ! -- Write output saturation file (out_sratf)
            call write_out_sratf(time_val)
        end if
      end if

      if (st_out_type%wtab == out_type(2)) then
        if (st_out_step%wtab == SZERO) then
          ! -- Write watertable file (out_wtabf)
            call write_out_wtabf(time_val)
        else if (mod(current_t,st_out_step%wtab) == 0) then
          ! -- Write watertable file (out_wtabf)
            call write_out_wtabf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write watertable file (out_wtabf)
            call write_out_wtabf(time_val)
        end if
      end if

      if (st_out_type%mass == out_type(1)) then
        if (st_out_step%mass == SZERO) then
          ! -- Write massbalance file (out_massf)
            call write_out_massf(time_val)
        else if (mod(current_t,st_out_step%mass) == 0) then
          ! -- Write massbalance file (out_massf)
            call write_out_massf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write massbalance file (out_massf)
            call write_out_massf(time_val)
        end if
      end if

      if (st_out_type%velc == out_type(3)) then
        if (st_out_step%velc == SZERO) then
          ! -- Write velocity file (out_velcf)
            call write_out_velcf(time_val)
        else if (mod(current_t,st_out_step%velc) == 0) then
          ! -- Write velocity file (out_velcf)
            call write_out_velcf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write velocity file (out_velcf)
            call write_out_velcf(time_val)
        end if
      end if

      if (st_out_type%rivr == out_type(2)) then
        if (st_out_step%rivr == SZERO) then
          ! -- Write river runoff file (out_rivrf)
            call write_out_rivrf(time_val)
        else if (mod(current_t,st_out_step%rivr) == 0) then
          ! -- Write river runoff file (out_rivrf)
            call write_out_rivrf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write river runoff file (out_rivrf)
            call write_out_rivrf(time_val)
        end if
      end if

      if (st_out_type%lakr == out_type(2)) then
        if (st_out_step%lakr == SZERO) then
          ! -- Write lake runoff file (out_lakrf)
            call write_out_lakrf(time_val)
        else if (mod(current_t,st_out_step%lakr) == 0) then
          ! -- Write lake runoff file (out_lakrf)
            call write_out_lakrf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write lake runoff file (out_lakrf)
            call write_out_lakrf(time_val)
        end if
      end if

      if (st_out_type%sufr == out_type(2)) then
        if (st_out_step%sufr == SZERO) then
          ! -- Write surface runoff file (out_sufrf)
            call write_out_sufrf(time_val)
        else if (mod(current_t,st_out_step%sufr) == 0) then
          ! -- Write surface runoff file (out_sufrf)
            call write_out_sufrf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write surface runoff file (out_sufrf)
            call write_out_sufrf(time_val)
        end if
      end if

      if (st_out_type%dunr == out_type(2)) then
        if (st_out_step%dunr == SZERO) then
          ! -- Write dunne runoff file (out_dunrf)
            call write_out_dunrf(time_val)
        else if (mod(current_t,st_out_step%dunr) == 0) then
          ! -- Write dunne runoff file (out_dunrf)
            call write_out_dunrf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write dunne runoff file (out_dunrf)
            call write_out_dunrf(time_val)
        end if
      end if

      if (st_out_type%seal == out_type(3)) then
        if (st_out_step%seal == SZERO) then
          ! -- Write sea results file (out_sealf)
            call write_out_sealf(time_val)
        else if (mod(current_t,st_out_step%seal) == 0) then
          ! -- Write sea results file (out_sealf)
            call write_out_sealf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write sea results file (out_sealf)
            call write_out_sealf(time_val)
        end if
      end if

      if (st_out_type%rech == out_type(2)) then
        if (st_out_step%rech == SZERO) then
          ! -- Write recharge results file (out_rechf)
            call write_out_rechf(time_val)
        else if (mod(current_t,st_out_step%rech) == 0) then
          ! -- Write recharge results file (out_rechf)
            call write_out_rechf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write recharge results file (out_rechf)
            call write_out_rechf(time_val)
        end if
      end if

      if (st_out_type%well == out_type(3)) then
        if (st_out_step%well == SZERO) then
          ! -- Write well pumping results file (out_wellf)
            call write_out_wellf(time_val)
        else if (mod(current_t,st_out_step%well) == 0) then
          ! -- Write well pumping results file (out_wellf)
            call write_out_wellf(time_val)
        else if (lasttime_flag== 1) then
          ! -- Write well pumping results file (out_wellf)
            call write_out_wellf(time_val)
        end if
      end if

    end if

  end subroutine write_outf

  subroutine write_mass_header()
  !***************************************************************************************
  ! write_mass_header -- Write massbalance file header
  !***************************************************************************************
    ! -- modules
    use constval_module, only: MASSCHARA, OUTFORM
    use utility_module, only: conv_i2s
    use initial_module, only: st_sim
    use allocate_output, only: ms_head
    ! -- inout

    ! -- local
    integer(I4) :: i, mass_fnum
    character(:), allocatable :: str_mstnum
    integer(I4) :: msout_tnum_format
    !-------------------------------------------------------------------------------------
    mass_fnum = st_out_fnum%mass
    write(mass_fnum,'(3a)') trim(adjustl(st_sim%cal_unit)), ",", trim(adjustl(ms_head))

    do i = 1, msout_tnum
      if (i == 1) then
        ms_head = ","//trim(adjustl(MASSCHARA))
      else
        ms_head = trim(adjustl(ms_head))//","//trim(adjustl(MASSCHARA))
      end if
    end do

    write(mass_fnum,'(a)') (trim(adjustl(ms_head)))

    msout_tnum_format = msout_tnum*9

    str_mstnum = conv_i2s(msout_tnum_format)
    msformat = "("//OUTFORM//","//trim(adjustl(str_mstnum))//"(',',"//OUTFORM//"))"

  end subroutine write_mass_header

  subroutine write_out_headf(time_out)
  !***************************************************************************************
  ! write_out_headf -- Write output head file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, head_fnum
    integer(I4), allocatable :: calc2calc(:)
    !-------------------------------------------------------------------------------------
    allocate(calc2calc(ncalc))
    !$omp parallel
    !$omp workshare
    calc2calc(:) = 0
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncalc
      calc2calc(i) = i
    end do
    !$omp end do
    !$omp end parallel

    head_fnum = st_out_fnum%head

#ifdef MPI_MSG
    ! -- Write MPI 3D binary file (mpi_3dbin)
      call write_mpi_3dbin(head_fnum, ncalc, calc2calc, len_scal, head_new, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(head_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(head_fnum, time_out)
    ! -- Write 3D binary file (3dbin)
      call write_3dbin(head_fnum, ncalc, calc2calc, len_scal, head_new)
    if (lasttime_flag == 1) then
      call close_file(head_fnum)
    end if
#endif

    deallocate(calc2calc)

  end subroutine write_out_headf

  subroutine write_out_sratf(time_out)
  !***************************************************************************************
  ! write_out_sratf -- Write output saturation file
  !***************************************************************************************
    ! -- modules
    use constval_module, only: SONE
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, srat_fnum
    integer(I4), allocatable :: calc2calc(:)
    !-------------------------------------------------------------------------------------
    allocate(calc2calc(ncalc))
    !$omp parallel
    !$omp workshare
    calc2calc(:) = 0
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncalc
      calc2calc(i) = i
    end do
    !$omp end do
    !$omp end parallel

    srat_fnum = st_out_fnum%srat

#ifdef MPI_MSG
    ! -- Write MPI 3D binary file (mpi_3dbin)
      call write_mpi_3dbin(srat_fnum, ncalc, calc2calc, SONE, srat_new, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(srat_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(srat_fnum, time_out)
    ! -- Write 3D binary file (3dbin)
      call write_3dbin(srat_fnum, ncalc, calc2calc, SONE, srat_new)
    if (lasttime_flag == 1) then
      call close_file(srat_fnum)
    end if
#endif

    deallocate(calc2calc)

  end subroutine write_out_sratf

  subroutine write_out_restf(fnum_rest, time_out)
  !***************************************************************************************
  ! write_out_restf -- Write restart file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum_rest
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    rewind(fnum_rest)
    write(fnum_rest) real(time_out, kind=DP)
    write(fnum_rest) (head_new(i)*len_scal, i = 1, ncalc)
    call close_file(fnum_rest)

  end subroutine write_out_restf

  subroutine write_out_wtabf(time_out)
  !***************************************************************************************
  ! write_out_wtabf -- Write watertable file
  !***************************************************************************************
    ! -- modules
    use allocate_output, only: wtable
    use calc_output, only: calc_wtable
#ifdef MPI_MSG
    use mpi_write, only: calc_mpi_wtable
#endif
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, wtab_fnum
    integer(I4), allocatable :: cals2cals(:)
    !-------------------------------------------------------------------------------------
#ifdef MPI_MSG
    if (pro_totn /= 1) then
    ! -- Calculate water table for MPI (mpi_wtable)
      call calc_mpi_wtable(head_new, srat_new)
    else
    ! -- Calculate water table (wtable)
      call calc_wtable(head_new, srat_new)
    end if
#else
    ! -- Calculate water table (wtable)
      call calc_wtable(head_new, srat_new)
#endif

    allocate(cals2cals(ncals))
    !$omp parallel
    !$omp workshare
    cals2cals(:) = 0
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      cals2cals(i) = i
    end do
    !$omp end do
    !$omp end parallel

    wtab_fnum = st_out_fnum%wtab

#ifdef MPI_MSG
    ! -- Write MPI 2D binary file (mpi_2dbin)
      call write_mpi_2dbin(wtab_fnum, ncals, cals2cals, len_scal, wtable, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(wtab_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(wtab_fnum, time_out)
    ! -- Write 2D binary file (2dbin)
      call write_2dbin(wtab_fnum, ncals, cals2cals, len_scal, wtable)
    if (lasttime_flag == 1) then
      call close_file(wtab_fnum)
    end if
#endif

    deallocate(cals2cals)

  end subroutine write_out_wtabf

  subroutine write_out_massf(time_out)
  !***************************************************************************************
  ! write_out_massf -- Write massbalance file
  !***************************************************************************************
    ! -- modules
    use allocate_output, only: st_msloc, st_msglo
    use calc_output, only: calc_out_mass
#ifdef MPI_MSG
    use mpi_write, only: redu_mpi_mass
#endif
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, mass_fnum
    real(SP) :: cubic
    !-------------------------------------------------------------------------------------
    cubic = len_scal**3

    !$omp parallel workshare
    st_msloc%sto(:) = st_msloc%sto(:)*cubic ; st_msloc%con(:) = st_msloc%con(:)*cubic
    st_msloc%sea(:) = st_msloc%sea(:)*cubic ; st_msloc%wel(:) = st_msloc%wel(:)*cubic
    st_msloc%rec(:) = st_msloc%rec(:)*cubic ; st_msloc%sur(:) = st_msloc%sur(:)*cubic
    st_msloc%riv(:) = st_msloc%riv(:)*cubic ; st_msloc%lak(:) = st_msloc%lak(:)*cubic
    !$omp end parallel workshare

    ! -- Calculate output massbalance (out_mass)
      call calc_out_mass()

#ifdef MPI_MSG
    ! -- Reduce output massbalance for MPI (mpi_mass)
      call redu_mpi_mass(msout_tnum, st_msglo)
#endif

    mass_fnum = st_out_fnum%mass
    if (my_rank == 0) then
      write(mass_fnum,msformat) time_out, (st_msglo%con(i), st_msglo%sto(i),&
            st_msglo%rec(i), st_msglo%wel(i), st_msglo%sur(i), st_msglo%riv(i),&
            st_msglo%lak(i), st_msglo%sea(i), st_msglo%tot(i), i = 1, msout_tnum)
    end if

  end subroutine write_out_massf

  subroutine write_out_velcf(time_out)
  !***************************************************************************************
  ! write_out_velcf -- Write velocity file
  !***************************************************************************************
    ! -- modules
    use allocate_output, only: pointv
    use calc_output, only: calc_outvelc
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, velx_fnum, vely_fnum, velz_fnum
    integer(I4), allocatable :: calc2calc(:)
    real(DP), allocatable :: velcx(:), velcy(:), velcz(:)
    !-------------------------------------------------------------------------------------
    ! -- Calculate output velocity (outvelc)
      call calc_outvelc()

    allocate(calc2calc(ncalc))
    allocate(velcx(ncalc), velcy(ncalc), velcz(ncalc))
    !$omp parallel
    !$omp workshare
    calc2calc(:) = 0
    velcx(:) = DZERO ; velcy(:) = DZERO ; velcz(:) = DZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncalc
      calc2calc(i) = i
      velcx(i) = pointv(i,1)
      velcy(i) = pointv(i,2)
      velcz(i) = pointv(i,3)
    end do
    !$omp end do
    !$omp end parallel

    velx_fnum = st_out_fnum%velx
    vely_fnum = st_out_fnum%vely
    velz_fnum = st_out_fnum%velz

#ifdef MPI_MSG
    ! -- Write MPI 3D binary file (mpi_3dbin)
      call write_mpi_3dbin(velx_fnum, ncalc, calc2calc, len_scal, velcx, time_out)
      call write_mpi_3dbin(vely_fnum, ncalc, calc2calc, len_scal, velcy, time_out)
      call write_mpi_3dbin(velz_fnum, ncalc, calc2calc, len_scal, velcz, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(velx_fnum)
      call close_mpi_file(vely_fnum)
      call close_mpi_file(velz_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(velx_fnum, time_out)
    ! -- Write header binary file (header_bin)
      call write_header_bin(vely_fnum, time_out)
    ! -- Write header binary file (header_bin)
      call write_header_bin(velz_fnum, time_out)
    ! -- Write 3D binary file (3dbin)
      call write_3dbin(velx_fnum, ncalc, calc2calc, len_scal, velcx)
      call write_3dbin(vely_fnum, ncalc, calc2calc, len_scal, velcy)
      call write_3dbin(velz_fnum, ncalc, calc2calc, len_scal, velcz)
    if (lasttime_flag == 1) then
      call close_file(velx_fnum)
      call close_file(vely_fnum)
      call close_file(velz_fnum)
    end if
#endif

    deallocate(calc2calc)
    deallocate(velcx, velcy, velcz)

  end subroutine write_out_velcf

  subroutine write_out_rivrf(time_out)
  !***************************************************************************************
  ! write_out_rivrf -- Write river runoff file
  !***************************************************************************************
    ! -- modules
    use calc_boundary, only: rive2cals
    use set_boundary, only: rive_num
    use allocate_output, only: roff_rive, rive_sumtime
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, s, rivr_fnum
    real(DP), allocatable :: rive_flux(:)
    !-------------------------------------------------------------------------------------
    allocate(rive_flux(rive_num))
    !$omp parallel
    !$omp workshare
    rive_flux(:) = DZERO
    !$omp end workshare

    !$omp do private(i, s)
    do i = 1, rive_num
      s = rive2cals(i)
      rive_flux(i) = roff_rive(i)/rech_area(s)/rive_sumtime
    end do
    !$omp end do
    !$omp end parallel

    rivr_fnum = st_out_fnum%rivr

#ifdef MPI_MSG
    ! -- Write MPI 2D binary file (mpi_2dbin)
      call write_mpi_2dbin(rivr_fnum, rive_num, rive2cals, len_scal, rive_flux, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(rivr_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(rivr_fnum, time_out)
    ! -- Write 2D binary file (2dbin)
      call write_2dbin(rivr_fnum, rive_num, rive2cals, len_scal, rive_flux)
    if (lasttime_flag == 1) then
      call close_file(rivr_fnum)
    end if
#endif

    !$omp parallel workshare
    roff_rive(:) = DZERO
    !$omp end parallel workshare
    rive_sumtime = DZERO

    deallocate(rive_flux)

  end subroutine write_out_rivrf

  subroutine write_out_lakrf(time_out)
  !***************************************************************************************
  ! write_out_lakrf -- Write lake runoff file
  !***************************************************************************************
    ! -- modules
    use calc_boundary, only: lake2cals
    use set_boundary, only: lake_num
    use allocate_output, only: roff_lake, lake_sumtime
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, s, lakr_fnum
    real(DP), allocatable :: lake_flux(:)
    !-------------------------------------------------------------------------------------
    allocate(lake_flux(lake_num))
    !$omp parallel
    !$omp workshare
    lake_flux(:) = DZERO
    !$omp end workshare

    !$omp do private(i, s)
    do i = 1, lake_num
      s = lake2cals(i)
      lake_flux(i) = roff_lake(i)/rech_area(s)/lake_sumtime
    end do
    !$omp end do
    !$omp end parallel

    lakr_fnum = st_out_fnum%lakr

#ifdef MPI_MSG
    ! -- Write MPI 2D binary file (mpi_2dbin)
      call write_mpi_2dbin(lakr_fnum, lake_num, lake2cals, len_scal, lake_flux, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(lakr_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(lakr_fnum, time_out)
    ! -- Write 2D binary file (2dbin)
      call write_2dbin(lakr_fnum, lake_num, lake2cals, len_scal, lake_flux)
    if (lasttime_flag == 1) then
      call close_file(lakr_fnum)
    end if
#endif

    !$omp parallel workshare
    roff_lake(:) = DZERO
    !$omp end parallel workshare
    lake_sumtime = DZERO

    deallocate(lake_flux)

  end subroutine write_out_lakrf

  subroutine write_out_sufrf(time_out)
  !***************************************************************************************
  ! write_out_sufrf -- Write surface runoff file
  !***************************************************************************************
    ! -- modules
    use allocate_output, only: roff_surf, surf_sumtime
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, sufr_fnum
    integer(I4), allocatable :: cals2cals(:)
    real(DP), allocatable :: surf_flux(:)
    !-------------------------------------------------------------------------------------
    allocate(cals2cals(ncals))
    allocate(surf_flux(ncals))
    !$omp parallel
    !$omp workshare
    cals2cals(:) = 0
    surf_flux(:) = DZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, ncals
      cals2cals(i) = i
      surf_flux(i) = roff_surf(i)/rech_area(i)/surf_sumtime
    end do
    !$omp end do
    !$omp end parallel

    sufr_fnum = st_out_fnum%sufr

#ifdef MPI_MSG
    ! -- Write MPI 2D binary file (mpi_2dbin)
      call write_mpi_2dbin(sufr_fnum, ncals, cals2cals, len_scal, surf_flux, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(sufr_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(sufr_fnum, time_out)
    ! -- Write 2D binary file (2dbin)
      call write_2dbin(sufr_fnum, ncals, cals2cals, len_scal, surf_flux)
    if (lasttime_flag == 1) then
      call close_file(sufr_fnum)
    end if
#endif

    !$omp parallel workshare
    roff_surf(:) = DZERO
    !$omp end parallel workshare
    surf_sumtime = DZERO

    deallocate(cals2cals)
    deallocate(surf_flux)

  end subroutine write_out_sufrf

  subroutine write_out_dunrf(time_out)
  !***************************************************************************************
  ! write_out_dunrf -- Write dunne runoff file
  !***************************************************************************************
    ! -- modules
    use calc_boundary, only: rech2cals
    use set_boundary, only: rech_num
    use allocate_output, only: roff_dunn, dunn_sumtime
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: i, dunr_file
    real(DP), allocatable :: dunn_flux(:)
    !-------------------------------------------------------------------------------------
    allocate(dunn_flux(rech_num))
    !$omp parallel
    !$omp workshare
    dunn_flux(:) = DZERO
    !$omp end workshare

    !$omp do private(i)
    do i = 1, rech_num
      dunn_flux(i) = roff_dunn(i)/dunn_sumtime
    end do
    !$omp end do
    !$omp end parallel

    dunr_file = st_out_fnum%dunr

#ifdef MPI_MSG
    ! -- Write MPI 2D binary file (mpi_2dbin)
      call write_mpi_2dbin(dunr_file, rech_num, rech2cals, len_scal, dunn_flux, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(dunr_file)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(dunr_file, time_out)
    ! -- Write 2D binary file (2dbin)
      call write_2dbin(dunr_file, rech_num, rech2cals, len_scal, dunn_flux)
    if (lasttime_flag == 1) then
      call close_file(dunr_file)
    end if
#endif

    !$omp parallel workshare
    roff_dunn(:) = DZERO
    !$omp end parallel workshare
    dunn_sumtime = DZERO

    deallocate(dunn_flux)

  end subroutine write_out_dunrf

  subroutine write_out_sealf(time_out)
  !***************************************************************************************
  ! write_out_sealf -- Write sea results file
  !***************************************************************************************
    ! -- modules
    use allocate_output, only: res_snum, res_seal
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: seal_fnum
    real(SP) :: cubic
    !-------------------------------------------------------------------------------------
    seal_fnum = st_out_fnum%seal ; cubic = len_scal**3
#ifdef MPI_MSG
    ! -- Write MPI 3D binary file (mpi_3dbin)
      call write_mpi_3dbin(seal_fnum, ncalc, res_snum, cubic, res_seal, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(seal_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(seal_fnum, time_out)
    ! -- Write 3D binary file (3dbin)
      call write_3dbin(seal_fnum, ncalc, res_snum, cubic, res_seal)
    if (lasttime_flag == 1) then
      call close_file(seal_fnum)
    end if
#endif

    !$omp parallel workshare
    res_snum(:) = 0 ; res_seal(:) = DZERO
    !$omp end parallel workshare

  end subroutine write_out_sealf

  subroutine write_out_rechf(time_out)
  !***************************************************************************************
  ! write_out_rechf -- Write recharge results file
  !***************************************************************************************
    ! -- modules
    use allocate_output, only: res_rnum, res_rech
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: rech_fnum
    real(SP) :: cubic
    !-------------------------------------------------------------------------------------
    rech_fnum = st_out_fnum%rech ; cubic = len_scal**3
#ifdef MPI_MSG
    ! -- Write MPI 2D binary file (mpi_2dbin)
      call write_mpi_2dbin(rech_fnum, ncals, res_rnum, cubic, res_rech, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(rech_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(rech_fnum, time_out)
    ! -- Write 2D binary file (2dbin)
      call write_2dbin(rech_fnum, ncals, res_rnum, cubic, res_rech)
    if (lasttime_flag == 1) then
      call close_file(rech_fnum)
    end if
#endif

    !$omp parallel workshare
    res_rnum(:) = 0 ; res_rech(:) = DZERO
    !$omp end parallel workshare

  end subroutine write_out_rechf

  subroutine write_out_wellf(time_out)
  !***************************************************************************************
  ! write_out_wellf -- Write well pumping results file
  !***************************************************************************************
    ! -- modules
    use allocate_output, only: res_wnum, res_well
    ! -- inout
    real(SP), intent(in) :: time_out
    ! -- local
    integer(I4) :: well_fnum
    real(SP) :: cubic
    !-------------------------------------------------------------------------------------
    well_fnum = st_out_fnum%well ; cubic = len_scal**3
#ifdef MPI_MSG
    ! -- Write MPI 3D binary file (mpi_3dbin)
      call write_mpi_3dbin(well_fnum, ncalc, res_wnum, cubic, res_well, time_out)
    if (lasttime_flag == 1) then
      call close_mpi_file(well_fnum)
    end if
#else
    ! -- Write header binary file (header_bin)
      call write_header_bin(well_fnum, time_out)
    ! -- Write 3D binary file (3dbin)
      call write_3dbin(well_fnum, ncalc, res_wnum, cubic, res_well)
    if (lasttime_flag == 1) then
      call close_file(well_fnum)
    end if
#endif

    !$omp parallel workshare
    res_wnum(:) = 0 ; res_well(:) = DZERO
    !$omp end parallel workshare

  end subroutine write_out_wellf

end module write_output
