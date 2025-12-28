module assign_calc
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SZERO, SNOVAL, DNOVAL, VARLEN
  use utility_module, only: write_err_stop, close_file
  use initial_module, only: pro_totn, my_rank, st_grid, st_init, in_type
  use read_input, only: len_scal, len_scal_inv
  use set_cell, only: ncals, ncalc, glo2loc_ijk
  use set_condition, only: set_clas2calc, set_2dfile2calc, set_3dfile2calc
#ifdef MPI_MSG
  use mpi_read, only: close_mpi_file
#endif

  implicit none
  private
  public :: assign_retnv, assign_parmv, assign_geogv, assign_initv, assign_massv
  integer(I4), public :: geog_num = 0, mass_num = 0, msout_tnum
  integer(I4), allocatable, public :: int_mass(:), mass2calc(:)
  real(SP), allocatable, public :: read_reta(:), read_retn(:), read_resi(:)
  real(SP), allocatable, public :: read_ksx(:), read_ksy(:), read_ksz(:)
  real(SP), allocatable, public :: read_ss(:), read_poro(:)
  real(DP), allocatable, public :: surf_bott(:), surf_reli(:), surf_parm(:)
  real(DP), allocatable, public :: read_init(:)
  character(VARLEN), allocatable, public :: massout_name(:)

  ! -- local

  contains

  subroutine assign_retnv(retn_ftype)
  !***************************************************************************************
  ! assign_retnv -- Assign retention value
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_retf_type, st_retn, st_retn_fnum
#ifdef MPI_MSG
    use mpi_set, only: bcast_retn_clas
#endif
    ! -- inout
    integer(I4), intent(in) :: retn_ftype
    ! -- local
    integer(I4) :: i, ierr
    integer(I4), allocatable :: ret_fnum(:), ret_ftype(:)
    real(SP), allocatable :: ret_val(:)
    !-------------------------------------------------------------------------------------
    allocate(read_reta(ncalc), read_retn(ncalc), read_resi(ncalc))
    !$omp parallel workshare
    read_reta(:) = SNOVAL ; read_retn(:) = SNOVAL ; read_resi(:) = SNOVAL
    !$omp end parallel workshare

    if (retn_ftype == in_type(0)) then
      allocate(ret_fnum(3), ret_ftype(3))
      allocate(ret_val(ncalc))
      ret_fnum(1) = st_retn_fnum%vana ; ret_fnum(2) = st_retn_fnum%vann
      ret_fnum(3) = st_retn_fnum%resi
      ret_ftype(1) = st_retf_type%vana ; ret_ftype(2) = st_retf_type%vann
      ret_ftype(3) = st_retf_type%resi

      do i = 1, 3
        !$omp parallel workshare
        ret_val(:) = SNOVAL
        !$omp end parallel workshare
        if (ret_ftype(i) == in_type(3) .or. ret_ftype(i) == in_type(4)) then
          call set_2dfile2calc(ret_fnum(i), ret_ftype(i), 0, ncals, SNOVAL, ret_val)
        else if (ret_ftype(i) == in_type(5) .or. ret_ftype(i) == in_type(6)) then
          call set_3dfile2calc(ret_fnum(i), ret_ftype(i), 0, SNOVAL, ret_val)
        end if

        select case (i)
        case(1)
          !$omp parallel workshare
          read_reta(:) = ret_val(:)
          !$omp end parallel workshare
        case(2)
          !$omp parallel workshare
          read_retn(:) = ret_val(:)
          !$omp end parallel workshare
        case(3)
          !$omp parallel workshare
          read_resi(:) = ret_val(:)
          !$omp end parallel workshare
        end select
      end do
      deallocate(ret_val)

    else if (retn_ftype == in_type(1)) then
      allocate(st_retn%name(st_retn%totn))
      allocate(st_retn%a(st_retn%totn), st_retn%n(st_retn%totn))
      allocate(st_retn%r(st_retn%totn))
      !$omp parallel workshare
      st_retn%name(:) = ""
      st_retn%a(:) = SNOVAL ; st_retn%n(:) = SNOVAL ; st_retn%r(:) = SNOVAL
      !$omp end parallel workshare

      if (my_rank == 0) then
        ierr = 0
        do i = 1, st_retn%totn
          read(unit=st_retn%fnum,fmt=*,iostat=ierr) st_retn%name(i), st_retn%a(i),&
                                                    st_retn%n(i), st_retn%r(i)
          if (ierr /= 0) then
            call write_err_stop("While reading in retention file.")
          end if
        end do
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        call bcast_retn_clas(st_retn%totn, st_retn%name, st_retn%a, st_retn%n, st_retn%r)
      end if
#endif
      call set_clas2calc(st_retn%totn, st_retn%name, st_retn%a, read_reta)
      call set_clas2calc(st_retn%totn, st_retn%name, st_retn%n, read_retn)
      call set_clas2calc(st_retn%totn, st_retn%name, st_retn%r, read_resi)
      deallocate(st_retn%a, st_retn%n, st_retn%r)
    end if

    !$omp parallel workshare
    read_reta(:) = read_reta(:)*len_scal
    !$omp end parallel workshare

#ifdef MPI_MSG
    if (my_rank == 0) then
      call close_file(st_retn%fnum)
    end if
    if (retn_ftype == in_type(0)) then
      do i = 1, 3
        if (ret_ftype(i) == in_type(3) .or. ret_ftype(i) == in_type(5)) then
          if (my_rank == 0) then
            call close_file(ret_fnum(i))
          end if
        else if (ret_ftype(i) == in_type(4) .or. ret_ftype(i) == in_type(6)) then
          call close_mpi_file(ret_fnum(i))
        end if
      end do
    end if
#else
    !$omp parallel sections
    !$omp section
    call close_file(st_retn%fnum)
    !$omp section
    call close_file(st_retn_fnum%vana)
    !$omp section
    call close_file(st_retn_fnum%vann)
    !$omp section
    call close_file(st_retn_fnum%resi)
    !$omp end parallel sections
#endif

    if (retn_ftype == in_type(0)) then
      deallocate(ret_fnum, ret_ftype)
    end if

  end subroutine assign_retnv

  subroutine assign_parmv(parm_ftype)
  !***************************************************************************************
  ! assign_parmv -- Assign parameter value
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_parf_type, st_parm, st_parm_fnum
#ifdef MPI_MSG
    use mpi_set, only: bcast_parm_clas
#endif
    ! -- inout
    integer(I4), intent(in) :: parm_ftype
    ! -- local
    integer(I4) :: i, ierr
    integer(I4), allocatable :: par_ftype(:), par_fnum(:)
    real(SP), allocatable :: par_val(:)
    !-------------------------------------------------------------------------------------
    allocate(read_ksx(ncalc), read_ksy(ncalc), read_ksz(ncalc))
    allocate(read_ss(ncalc), read_poro(ncalc))
    !$omp parallel workshare
    read_ksx(:) = SNOVAL ; read_ksy(:) = SNOVAL ; read_ksz(:) = SNOVAL
    read_ss(:) = SNOVAL ; read_poro(:) = SNOVAL
    !$omp end parallel workshare

    if (parm_ftype == in_type(0)) then
      allocate(par_fnum(5), par_ftype(5))
      allocate(par_val(ncalc))
      par_fnum(1) = st_parm_fnum%pakx ; par_fnum(2) = st_parm_fnum%paky
      par_fnum(3) = st_parm_fnum%pakz ; par_fnum(4) = st_parm_fnum%pass
      par_fnum(5) = st_parm_fnum%pats
      par_ftype(1) = st_parf_type%pakx ; par_ftype(2) = st_parf_type%paky
      par_ftype(3) = st_parf_type%pakz ; par_ftype(4) = st_parf_type%pass
      par_ftype(5) = st_parf_type%pats

      do i = 1, 5
        !$omp parallel workshare
        par_val(:) = SNOVAL
        !$omp end parallel workshare
        if (par_ftype(i) == in_type(3) .or. par_ftype(i) == in_type(4)) then
          call set_2dfile2calc(par_fnum(i), par_ftype(i), 0, ncals, SNOVAL, par_val)
        else if (par_ftype(i) == in_type(5) .or. par_ftype(i) == in_type(6)) then
          call set_3dfile2calc(par_fnum(i), par_ftype(i), 0, SNOVAL, par_val)
        end if

        select case (i)
        case(1)
          !$omp parallel workshare
          read_ksx(:) = par_val(:)
          !$omp end parallel workshare
        case(2)
          !$omp parallel workshare
          read_ksy(:) = par_val(:)
          !$omp end parallel workshare
        case(3)
          !$omp parallel workshare
          read_ksz(:) = par_val(:)
          !$omp end parallel workshare
        case(4)
          !$omp parallel workshare
          read_ss(:) = par_val(:)
          !$omp end parallel workshare
        case(5)
          !$omp parallel workshare
          read_poro(:) = par_val(:)
          !$omp end parallel workshare
        end select
      end do
      deallocate(par_val)

    else if (parm_ftype == in_type(1)) then
      allocate(st_parm%name(st_parm%totn))
      allocate(st_parm%ksx(st_parm%totn), st_parm%ksy(st_parm%totn))
      allocate(st_parm%ksz(st_parm%totn))
      allocate(st_parm%ss(st_parm%totn), st_parm%ts(st_parm%totn))
      !$omp parallel workshare
      st_parm%name(:) = ""
      st_parm%ksx(:) = SNOVAL ; st_parm%ksy(:) = SNOVAL ; st_parm%ksz(:) = SNOVAL
      st_parm%ss(:) = SNOVAL ; st_parm%ts(:) = SNOVAL
      !$omp end parallel workshare

      if (my_rank == 0) then
        do i = 1, st_parm%totn
          read(unit=st_parm%fnum,fmt=*,iostat=ierr) st_parm%name(i), st_parm%ksx(i),&
                                                    st_parm%ksy(i), st_parm%ksz(i),&
                                                    st_parm%ss(i), st_parm%ts(i)
          if (ierr /= 0) then
            call write_err_stop("While reading in parameter file.")
          end if
        end do
      end if

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        call bcast_parm_clas(st_parm%totn, st_parm%name, st_parm%ksx, st_parm%ksy,&
                             st_parm%ksz, st_parm%ss, st_parm%ts)
      end if
#endif
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ksx, read_ksx)
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ksy, read_ksy)
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ksz, read_ksz)
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ss, read_ss)
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ts, read_poro)
      deallocate(st_parm%ksx, st_parm%ksy, st_parm%ksz, st_parm%ss, st_parm%ts)
    end if

    !$omp parallel workshare
    read_ksx(:) = read_ksx(:)*len_scal_inv ; read_ksy(:) = read_ksy(:)*len_scal_inv
    read_ksz(:) = read_ksz(:)*len_scal_inv ; read_ss(:) = read_ss(:)*len_scal
    !$omp end parallel workshare

#ifdef MPI_MSG
    if (my_rank == 0) then
      call close_file(st_parm%fnum)
    end if
    if (parm_ftype == in_type(0)) then
      do i = 1, 5
        if (par_ftype(i) == in_type(3) .or. par_ftype(i) == in_type(5)) then
          if (my_rank == 0) then
            call close_file(par_fnum(i))
          end if
        else if (par_ftype(i) == in_type(4) .or. par_ftype(i) == in_type(6)) then
          call close_mpi_file(par_fnum(i))
        end if
      end do
    end if
#else
    !$omp parallel sections
    !$omp section
    call close_file(st_parm%fnum)
    !$omp section
    call close_file(st_parm_fnum%pakx)
    !$omp section
    call close_file(st_parm_fnum%paky)
    !$omp section
    call close_file(st_parm_fnum%pakz)
    !$omp section
    call close_file(st_parm_fnum%pass)
    !$omp section
    call close_file(st_parm_fnum%pats)
    !$omp end parallel sections
#endif

    if (parm_ftype == in_type(0)) then
      deallocate(par_fnum, par_ftype)
    end if

  end subroutine assign_parmv

  subroutine assign_geogv()
  !***************************************************************************************
  ! assign_geogv -- Assign geography value
  !***************************************************************************************
    ! -- modules
    use constval_module, only: DZERO
    use initial_module, only: st_geof_type, st_geog_fnum
    use set_condition, only: set_2dfile2cals
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: geo_fnum(:), geo_ftype(:)
    integer(I4), allocatable :: geo_cflag(:)
    real(SP), allocatable :: geo_val(:)
    !-------------------------------------------------------------------------------------
    allocate(geo_fnum(3), geo_ftype(3))
    allocate(geo_cflag(ncals), surf_parm(ncals))
    allocate(geo_val(ncals))
    !$omp parallel workshare
    geo_cflag(:) = 0 ; surf_parm(:) = DZERO
    !$omp end parallel workshare

    geo_fnum(1) = st_geog_fnum%geoz ; geo_fnum(2) = st_geog_fnum%geor
    geo_fnum(3) = st_geog_fnum%geoa
    geo_ftype(1) = st_geof_type%geoz ; geo_ftype(2) = st_geof_type%geor
    geo_ftype(3) = st_geof_type%geoa

    do i = 1, 3
      !$omp parallel workshare
      geo_val(:) = SNOVAL
      !$omp end parallel workshare

      call set_2dfile2cals(geo_fnum(i), geo_ftype(i), 0, SNOVAL, geo_val, geo_cflag, geog_num)

      select case (i)
      case(1)
        !$omp parallel workshare
        surf_bott(:) = geo_val(:)
        !$omp end parallel workshare
      case(2)
        !$omp parallel workshare
        surf_reli(:) = geo_val(:)
        !$omp end parallel workshare
      case(3)
        !$omp parallel workshare
        surf_parm(:) = geo_val(:)
        !$omp end parallel workshare
      end select
    end do

    deallocate(geo_val, geo_cflag)

    !$omp parallel workshare
    surf_bott(:) = surf_bott(:)*len_scal_inv
    surf_reli(:) = surf_reli(:)*len_scal_inv
    !$omp end parallel workshare

#ifdef MPI_MSG
    do i = 1, 3
      if (geo_ftype(i) == in_type(3)) then
        if (my_rank == 0) then
          call close_file(geo_fnum(i))
        end if
      else if (geo_ftype(i) == in_type(4)) then
        call close_mpi_file(geo_fnum(i))
      end if
    end do
#else
    !$omp parallel sections
    !$omp section
    call close_file(st_geog_fnum%geoz)
    !$omp section
    call close_file(st_geog_fnum%geor)
    !$omp section
    call close_file(st_geog_fnum%geoa)
    !$omp end parallel sections
#endif

    deallocate(geo_fnum, geo_ftype)

  end subroutine assign_geogv

  subroutine assign_initv(init_ftype, init_unit)
  !***************************************************************************************
  ! assign_initv -- Assign initial value
  !***************************************************************************************
    ! -- modules
    use set_cell, only: get_calc_grid
#ifdef MPI_MSG
    use mpi_utility, only: bcast_val, mpisum_val
    use mpi_read, only: read_mpi_restf, read_mpi_head
    use mpi_set, only: bcast_init_dep, senrec_neibval
    use set_cell, only: neib_ncalc, neib_mpi_totn, neib_num, send_cind, recv_cind,&
                        send_citem, recv_citem, calc2recv
#endif
    ! -- inout
    integer(I4), intent(in) :: init_ftype
    character(*), intent(in) :: init_unit
    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: init_fnum
    real(DP) :: temp_end
#ifdef MPI_MSG
    integer(I4) :: j, k, xn, yn, zn, xyzn
    integer(I4) :: nov_num, sum_init
    real(DP), allocatable :: recv_init(:)
#endif
    !-------------------------------------------------------------------------------------
    allocate(read_init(ncalc))
    !$omp parallel workshare
    read_init(:) = DNOVAL
    !$omp end parallel workshare

    ierr = 0
    if (init_ftype == in_type(0)) then
#ifdef MPI_MSG
    ! -- Read mpi restart file (mpi_restf)
      call read_mpi_restf(st_init%fnum, ncalc, temp_end, read_init)
      call close_mpi_file(st_init%fnum)
#else
      read(unit=st_init%fnum,iostat=ierr) temp_end
      if (ierr /= 0) then
        call write_err_stop("While reading header time in initial file.")
      end if
      read(unit=st_init%fnum,iostat=ierr) (read_init(i), i = 1, ncalc)
      if (ierr /= 0) then
        call write_err_stop("While reading initial value in initial file.")
      end if
      call close_file(st_init%fnum)
#endif
    st_init%rest_time = real(temp_end, kind=SP)
    !$omp parallel workshare
    read_init(:) = read_init(:)*len_scal_inv
    !$omp end parallel workshare

    else if (init_ftype /= in_type(7)) then
      init_fnum = st_init%fnum
      if (len_trim(adjustl(init_unit)) /= 0) then
        if (init_ftype == in_type(3) .or. init_ftype == in_type(5)) then
          if (my_rank == 0) then
            read(unit=init_fnum,fmt=*,iostat=ierr) st_init%rest_time
            if (ierr /= 0) then
              call write_err_stop("While reading header time in initial file.")
            end if
          end if
        else if (init_ftype == in_type(4) .or. init_ftype == in_type(6)) then
#ifdef MPI_MSG
          ! -- Read mpi header (mpi_head)
            call read_mpi_head(init_fnum, ierr, temp_end)
#else
          read(unit=init_fnum,iostat=ierr) temp_end
#endif
          st_init%rest_time = real(temp_end, kind=SP)
          if (ierr /= 0 .and. my_rank == 0) then
            call write_err_stop("While reading header time in initial file.")
          end if
        end if
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(st_init%rest_time, "restart time value")
        end if
#endif
      end if
      if (init_ftype == in_type(3) .or. init_ftype == in_type(4)) then
        call set_2dfile2calc(init_fnum, init_ftype, 0, ncals, DNOVAL, read_init)
      else if (init_ftype == in_type(5) .or. init_ftype == in_type(6)) then
        call set_3dfile2calc(init_fnum, init_ftype, 0, DNOVAL, read_init)
      end if
      !$omp parallel workshare
      read_init(:) = read_init(:)*len_scal_inv
      !$omp end parallel workshare
      if (init_ftype == in_type(3) .or. init_ftype == in_type(5)) then
        if (my_rank == 0) then
          call close_file(init_fnum)
        end if
      else if (init_ftype == in_type(4) .or. init_ftype == in_type(6)) then
#ifdef MPI_MSG
        call close_mpi_file(init_fnum)
#else
        call close_file(init_fnum)
#endif
      end if
    else if (init_ftype == in_type(7)) then
      st_init%rest_time = SZERO
      st_init%depth = st_init%depth*len_scal_inv
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast initial depth (init_dep)
          call bcast_init_dep()
      end if
#endif
      ! -- Set initial waterlevel from surface elevation (init_elev)
        call set_init_elev()
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      allocate(recv_init(ncalc+neib_ncalc))
      !$omp parallel workshare
      sum_init = 1 ; recv_init(:) = DNOVAL
      !$omp end parallel workshare
      do while (sum_init /= 0)
        ! -- Send and Receive neighbor value (neibval)
          call senrec_neibval(neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                              recv_citem, read_init, recv_init)
        do i = 1, ncalc
          if (calc2recv(i) > 0 .and. read_init(i) == DNOVAL) then
            read_init(i) = recv_init(calc2recv(i))
            call get_calc_grid(i, xn, yn, zn)
            do k = zn+1, st_grid%nz
              xyzn = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(yn-1) + xn
              j = 0 ; j = glo2loc_ijk(xyzn)
              if (j /= 0 .and. j <= ncalc) then
                read_init(j) = read_init(i)
              end if
            end do
          end if
        end do
        nov_num = count(read_init(:) == DNOVAL)
        ! -- Sum value for MPI (val)
          call mpisum_val(nov_num, "non initial", sum_init)
      end do
      deallocate(recv_init)
    end if
    deallocate(calc2recv)
#endif

    st_init%rest_time = st_init%rest_time*st_init%multi

  end subroutine assign_initv

  subroutine assign_massv(mass_ftype)
  !***************************************************************************************
  ! assign_massv -- Assign massbalance value
  !***************************************************************************************
    ! -- modules
    use open_file, only: inmas_fnum
    use set_condition, only: set_mass2calc
#ifdef MPI_MSG
    use mpi_utility, only: bcast_val, mpimax_val
    use mpi_set, only: bcast_clas_val
#endif
    ! -- inout
    integer(I4), intent(in) :: mass_ftype
    ! -- local
    integer(I4) :: i, ierr
    integer(I4) :: massun
    integer(I4), allocatable :: mass_cflag(:)
    integer(I4), allocatable :: all_mass_type(:)
    integer(I4), allocatable :: mass_val(:), calc_mass(:)
    real(SP), allocatable :: clas_mass(:), real_mass(:)
    !-------------------------------------------------------------------------------------
    allocate(mass_val(ncalc), mass_cflag(ncalc))
    !$omp parallel workshare
    mass_val(:) = 0 ; mass_cflag(:) = 0
    !$omp end parallel workshare

    if (mass_ftype == in_type(1)) then
      if (my_rank == 0) then
        read(unit=inmas_fnum,fmt=*,iostat=ierr) msout_tnum
        if (ierr /= 0) then
          call write_err_stop("While reading massbalance output number in mass file.")
        else if (msout_tnum <= 0) then
          call write_err_stop("Specified wrong number for massbalance output.")
        end if
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
         call bcast_val(msout_tnum, "massbalance output number")
      end if
#endif
      allocate(massout_name(msout_tnum))
      allocate(clas_mass(msout_tnum), real_mass(ncalc))
      !$omp parallel workshare
      massout_name(:) = "" ; clas_mass(:) = SZERO ; real_mass(:) = SZERO
      !$omp end parallel workshare
      if (my_rank == 0) then
        do i = 1, msout_tnum
          read(unit=inmas_fnum,fmt=*,iostat=ierr) massout_name(i)
          if (ierr /= 0) then
            call write_err_stop("While reading massbalance output name in mass file.")
          end if
          clas_mass(i) = real(i, kind=SP)
        end do
        call close_file(inmas_fnum)
      end if

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        call bcast_clas_val(msout_tnum, massout_name, clas_mass)
      end if
#endif

      call set_clas2calc(msout_tnum, massout_name, clas_mass, real_mass, mass_cflag, mass_num)
      !$omp parallel workshare
      mass_val(:) = nint(real_mass(:))
      !$omp end parallel workshare
      deallocate(clas_mass, real_mass)

    else if (mass_ftype /= in_type(7)) then
        massun = inmas_fnum
      if (mass_ftype == in_type(3) .or. mass_ftype == in_type(4)) then
        call set_2dfile2calc(massun, mass_ftype, 0, ncals, 0, mass_val, mass_cflag, mass_num)
      else if (mass_ftype == in_type(5) .or. mass_ftype == in_type(6)) then
        call set_3dfile2calc(massun, mass_ftype, 0, 0, mass_val, mass_cflag, mass_num)
      end if

      if (mass_ftype == in_type(3) .or. mass_ftype == in_type(5)) then
        if (my_rank == 0) then
          call close_file(massun)
        end if
      else if (mass_ftype == in_type(4) .or. mass_ftype == in_type(6)) then
#ifdef MPI_MSG
        call close_mpi_file(massun)
#else
        call close_file(massun)
#endif
      end if
    else
      mass_num = ncalc
      !$omp parallel workshare
      mass_val(:) = 1 ; mass_cflag(:) = 1
      !$omp end parallel workshare
    end if

    allocate(all_mass_type(6))
    all_mass_type(:) = [in_type(1), in_type(3:7)]

    if (any(mass_ftype == all_mass_type(:))) then
      allocate(mass2calc(mass_num), calc_mass(mass_num))
      !$omp parallel workshare
      mass2calc(:) = 0 ; calc_mass(:) = 0
      !$omp end parallel workshare
      call set_mass2calc(ncalc, mass_cflag, mass_val, mass2calc, calc_mass)
      allocate(int_mass(mass_num))
      !$omp parallel workshare
      int_mass(:) = 0 ; int_mass(:) = calc_mass(:)
      !$omp end parallel workshare
      deallocate(calc_mass, mass_cflag)
#ifdef MPI_MSG
      ! -- Max value for MPI (val)
        call mpimax_val(maxval(int_mass), "massbalance number", msout_tnum)
#else
      msout_tnum = maxval(int_mass)
#endif
    end if

    deallocate(mass_val)
    deallocate(all_mass_type)

  end subroutine assign_massv

  subroutine set_init_elev()
  !***************************************************************************************
  ! set_init_elev -- Set initial waterlevel from surface elevation
  !***************************************************************************************
    ! -- modules
    use set_cell, only: get_cals_grid
    use make_cell, only: surf_elev
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, xn, yn, xyzn
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, xn, yn, xyzn)
    do i = 1, ncals
      call get_cals_grid(i, xn, yn)
      do k = 1, st_grid%nz
        xyzn = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(yn-1) + xn
        j = 0 ; j = glo2loc_ijk(xyzn)
        if (j /= 0 .and. j <= ncalc) then
          read_init(j) = surf_elev(i) - st_init%depth
        end if
      end do
    end do
    !$omp end parallel do

  end subroutine set_init_elev

end module assign_calc
