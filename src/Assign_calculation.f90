module assign_calc
  ! -- modules
  use kind_module, only: I4, SP
  use constval_module, only: VARLEN, SZERO, SNOVAL, DNOVAL
  use utility_module, only: st_mpi, close_file, write_err_stop, gmap_get
  use initial_module, only: in_type, st_grid, st_init
  use read_input, only: len_scal, len_scal_inv
  use set_cell, only: ncalc, ncals, st_conn
  use set_condition, only: set_clas2calc, set_2dfile2calc, set_3dfile2calc, st_hydr
#ifdef MPI_MSG
  use mpi_read, only: close_mpi_file
#endif

  implicit none
  private
  public :: assign_retnv, assign_parmv, assign_geogv, assign_initv, assign_massv
  integer(I4), public :: geog_num = 0, mass_num = 0, msout_tnum
  integer(I4), allocatable, public :: int_mass(:), mass2calc(:)
  character(VARLEN), allocatable, public :: massout_name(:)

  ! -- local

  contains

  subroutine assign_retnv(retn_ftype)
  !*********************************************************************************************
  ! assign_retnv -- Assign retention value
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_retf_type, st_retn_fnum, st_retn
#ifdef MPI_MSG
    use mpi_set, only: bcast_retn_clas
#endif
    ! -- inout
    integer(I4), intent(in) :: retn_ftype
    ! -- local
    integer(I4) :: i, j, ierr
    integer(I4), allocatable :: ret_fnum(:), ret_ftype(:)
    real(SP), allocatable :: ret_val(:)
    !-------------------------------------------------------------------------------------------
    allocate(st_hydr%read_vana(ncalc), st_hydr%read_vann(ncalc), st_hydr%read_resi(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      st_hydr%read_vana(i) = SNOVAL
      st_hydr%read_vann(i) = SNOVAL
      st_hydr%read_resi(i) = SNOVAL
    end do
    !$omp end parallel do

    if (retn_ftype == in_type(0)) then
      allocate(ret_fnum(3), ret_ftype(3))
      allocate(ret_val(ncalc))
      ret_fnum(1) = st_retn_fnum%vana ; ret_fnum(2) = st_retn_fnum%vann
      ret_fnum(3) = st_retn_fnum%resi
      ret_ftype(1) = st_retf_type%vana ; ret_ftype(2) = st_retf_type%vann
      ret_ftype(3) = st_retf_type%resi

      do i = 1, 3
        !$omp parallel do private(j)
        do j = 1, ncalc
          ret_val(j) = SNOVAL
        end do
        !$omp end parallel do
        if (ret_ftype(i) == in_type(3) .or. ret_ftype(i) == in_type(4)) then
          call set_2dfile2calc(ret_fnum(i), ret_ftype(i), 0, ncals, SNOVAL, ret_val)
        else if (ret_ftype(i) == in_type(5) .or. ret_ftype(i) == in_type(6)) then
          call set_3dfile2calc(ret_fnum(i), ret_ftype(i), 0, SNOVAL, ret_val)
        end if

        select case (i)
        case (1)
          !$omp parallel do private(j)
          do j = 1, ncalc
            st_hydr%read_vana(j) = ret_val(j)
          end do
          !$omp end parallel do
        case (2)
          !$omp parallel do private(j)
          do j = 1, ncalc
            st_hydr%read_vann(j) = ret_val(j)
          end do
          !$omp end parallel do
        case (3)
          !$omp parallel do private(j)
          do j = 1, ncalc
            st_hydr%read_resi(j) = ret_val(j)
          end do
          !$omp end parallel do
        end select
      end do
      deallocate(ret_val)

    else if (retn_ftype == in_type(1)) then
      allocate(st_retn%name(st_retn%totn))
      allocate(st_retn%a(st_retn%totn), st_retn%n(st_retn%totn))
      allocate(st_retn%r(st_retn%totn))
      !$omp parallel do private(i)
      do i = 1, st_retn%totn
        st_retn%name(i) = "" ; st_retn%a(i) = SNOVAL
        st_retn%n(i) = SNOVAL ; st_retn%r(i) = SNOVAL
      end do
      !$omp end parallel do

      if (st_mpi%rank == 0) then
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
      if (st_mpi%totn /= 1) then
        call bcast_retn_clas(st_retn%totn, st_retn%name, st_retn%a, st_retn%n, st_retn%r)
      end if
#endif
      call set_clas2calc(st_retn%totn, st_retn%name, st_retn%a, st_hydr%read_vana)
      call set_clas2calc(st_retn%totn, st_retn%name, st_retn%n, st_hydr%read_vann)
      call set_clas2calc(st_retn%totn, st_retn%name, st_retn%r, st_hydr%read_resi)
      deallocate(st_retn%a, st_retn%n, st_retn%r)
    end if

    !$omp parallel do private(i)
    do i = 1, ncalc
      st_hydr%read_vana(i) = st_hydr%read_vana(i)*len_scal
    end do
    !$omp end parallel do

#ifdef MPI_MSG
    if (st_mpi%rank == 0) then
      call close_file(st_retn%fnum)
    end if
    if (retn_ftype == in_type(0)) then
      do i = 1, 3
        if (ret_ftype(i) == in_type(3) .or. ret_ftype(i) == in_type(5)) then
          if (st_mpi%rank == 0) then
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
  !*********************************************************************************************
  ! assign_parmv -- Assign parameter value
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_parf_type, st_parm_fnum, st_parm
#ifdef MPI_MSG
    use mpi_set, only: bcast_parm_clas
#endif
    ! -- inout
    integer(I4), intent(in) :: parm_ftype
    ! -- local
    integer(I4) :: i, j, ierr
    integer(I4), allocatable :: par_ftype(:), par_fnum(:)
    real(SP), allocatable :: par_val(:)
    !-------------------------------------------------------------------------------------------
    allocate(st_hydr%read_hydx(ncalc), st_hydr%read_hydy(ncalc), st_hydr%read_hydz(ncalc))
    allocate(st_hydr%read_spst(ncalc), st_hydr%read_pors(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      st_hydr%read_hydx(i) = SNOVAL ; st_hydr%read_hydy(i) = SNOVAL
      st_hydr%read_hydz(i) = SNOVAL ; st_hydr%read_spst(i) = SNOVAL
      st_hydr%read_pors(i) = SNOVAL
    end do
    !$omp end parallel do

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
        !$omp parallel do private(j)
        do j = 1, ncalc
          par_val(j) = SNOVAL
        end do
        !$omp end parallel do
        if (par_ftype(i) == in_type(3) .or. par_ftype(i) == in_type(4)) then
          call set_2dfile2calc(par_fnum(i), par_ftype(i), 0, ncals, SNOVAL, par_val)
        else if (par_ftype(i) == in_type(5) .or. par_ftype(i) == in_type(6)) then
          call set_3dfile2calc(par_fnum(i), par_ftype(i), 0, SNOVAL, par_val)
        end if

        select case (i)
        case (1)
          !$omp parallel do private(j)
          do j = 1, ncalc
            st_hydr%read_hydx(j) = par_val(j)
          end do
          !$omp end parallel do
        case (2)
          !$omp parallel do private(j)
          do j = 1, ncalc
            st_hydr%read_hydy(j) = par_val(j)
          end do
          !$omp end parallel do
        case (3)
          !$omp parallel do private(j)
          do j = 1, ncalc
            st_hydr%read_hydz(j) = par_val(j)
          end do
          !$omp end parallel do
        case (4)
          !$omp parallel do private(j)
          do j = 1, ncalc
            st_hydr%read_spst(j) = par_val(j)
          end do
          !$omp end parallel do
        case (5)
          !$omp parallel do private(j)
          do j = 1, ncalc
            st_hydr%read_pors(j) = par_val(j)
          end do
          !$omp end parallel do
        end select
      end do
      deallocate(par_val)

    else if (parm_ftype == in_type(1)) then
      allocate(st_parm%name(st_parm%totn))
      allocate(st_parm%ksx(st_parm%totn), st_parm%ksy(st_parm%totn))
      allocate(st_parm%ksz(st_parm%totn))
      allocate(st_parm%ss(st_parm%totn), st_parm%ts(st_parm%totn))
      !$omp parallel do private(i)
      do i = 1, st_parm%totn
        st_parm%name(i) = "" ; st_parm%ksx(i) = SNOVAL ; st_parm%ksy(i) = SNOVAL
        st_parm%ksz(i) = SNOVAL ; st_parm%ss(i) = SNOVAL ; st_parm%ts(i) = SNOVAL
      end do
      !$omp end parallel do

      if (st_mpi%rank == 0) then
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
      if (st_mpi%totn /= 1) then
        call bcast_parm_clas(st_parm%totn, st_parm%name, st_parm%ksx, st_parm%ksy,&
                             st_parm%ksz, st_parm%ss, st_parm%ts)
      end if
#endif
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ksx, st_hydr%read_hydx)
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ksy, st_hydr%read_hydy)
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ksz, st_hydr%read_hydz)
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ss, st_hydr%read_spst)
      call set_clas2calc(st_parm%totn, st_parm%name, st_parm%ts, st_hydr%read_pors)
      deallocate(st_parm%ksx, st_parm%ksy, st_parm%ksz, st_parm%ss, st_parm%ts)
    end if

    !$omp parallel do private(i)
    do i = 1, ncalc
      st_hydr%read_hydx(i) = st_hydr%read_hydx(i)*len_scal_inv
      st_hydr%read_hydy(i) = st_hydr%read_hydy(i)*len_scal_inv
      st_hydr%read_hydz(i) = st_hydr%read_hydz(i)*len_scal_inv
      st_hydr%read_spst(i) = st_hydr%read_spst(i)*len_scal
    end do
    !$omp end parallel do

#ifdef MPI_MSG
    if (st_mpi%rank == 0) then
      call close_file(st_parm%fnum)
    end if
    if (parm_ftype == in_type(0)) then
      do i = 1, 5
        if (par_ftype(i) == in_type(3) .or. par_ftype(i) == in_type(5)) then
          if (st_mpi%rank == 0) then
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
  !*********************************************************************************************
  ! assign_geogv -- Assign geography value
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DZERO
    use initial_module, only: st_geof_type, st_geog_fnum
    use set_condition, only: set_2dfile2cals
    ! -- inout

    ! -- local
    integer(I4) :: i, j
    integer(I4), allocatable :: geo_fnum(:), geo_ftype(:)
    integer(I4), allocatable :: geo_cflag(:)
    real(SP), allocatable :: geo_val(:)
    !-------------------------------------------------------------------------------------------
    allocate(geo_fnum(3), geo_ftype(3))
    allocate(geo_cflag(ncals), st_hydr%surf_parm(ncals))
    allocate(geo_val(ncals))
    !$omp parallel do private(i)
    do i = 1, ncals
      geo_cflag(i) = 0
      st_hydr%surf_parm(i) = DZERO
    end do
    !$omp end parallel do

    geo_fnum(1) = st_geog_fnum%geoz ; geo_fnum(2) = st_geog_fnum%geor
    geo_fnum(3) = st_geog_fnum%geoa
    geo_ftype(1) = st_geof_type%geoz ; geo_ftype(2) = st_geof_type%geor
    geo_ftype(3) = st_geof_type%geoa

    do i = 1, 3
      !$omp parallel do private(j)
      do j = 1, ncals
        geo_val(j) = SNOVAL
      end do
      !$omp end parallel do

      call set_2dfile2cals(geo_fnum(i), geo_ftype(i), 0, SNOVAL, geo_val, geo_cflag, geog_num)

      select case (i)
      case (1)
        !$omp parallel do private(j)
        do j = 1, ncals
          st_hydr%surf_bott(j) = geo_val(j)
        end do
        !$omp end parallel do
      case (2)
        !$omp parallel do private(j)
        do j = 1, ncals
          st_hydr%surf_reli(j) = geo_val(j)
        end do
        !$omp end parallel do
      case (3)
        !$omp parallel do private(j)
        do j = 1, ncals
          st_hydr%surf_parm(j) = geo_val(j)
        end do
        !$omp end parallel do
      end select
    end do

    deallocate(geo_val, geo_cflag)

    !$omp parallel do private(i)
    do i = 1, ncals
      st_hydr%surf_bott(i) = st_hydr%surf_bott(i)*len_scal_inv
      st_hydr%surf_reli(i) = st_hydr%surf_reli(i)*len_scal_inv
    end do
    !$omp end parallel do

#ifdef MPI_MSG
    do i = 1, 3
      if (geo_ftype(i) == in_type(3)) then
        if (st_mpi%rank == 0) then
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
  !*********************************************************************************************
  ! assign_initv -- Assign initial value
  !*********************************************************************************************
    ! -- modules
    use kind_module, only: DP
    use set_cell, only: get_calc_grid
#ifdef MPI_MSG
    use mpi_utility, only: mpisum_val, bcast_val
    use mpi_read, only: read_mpi_restf, read_mpi_head
    use mpi_set, only: bcast_init_dep, senrec_neibval
    use set_cell, only: neib_mpi_totn, neib_ncalc, send_cind, recv_cind, send_citem,&
                        recv_citem, neib_num, calc2recv
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
    !-------------------------------------------------------------------------------------------
    allocate(st_hydr%read_init(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      st_hydr%read_init(i) = DNOVAL
    end do
    !$omp end parallel do

    ierr = 0
    if (init_ftype == in_type(0)) then
#ifdef MPI_MSG
    ! -- Read mpi restart file (mpi_restf)
      call read_mpi_restf(st_init%fnum, ncalc, temp_end, st_hydr%read_init)
      call close_mpi_file(st_init%fnum)
#else
      read(unit=st_init%fnum,iostat=ierr) temp_end
      if (ierr /= 0) then
        call write_err_stop("While reading header time in initial file.")
      end if
      read(unit=st_init%fnum,iostat=ierr) (st_hydr%read_init(i), i = 1, ncalc)
      if (ierr /= 0) then
        call write_err_stop("While reading initial value in initial file.")
      end if
      call close_file(st_init%fnum)
#endif
    st_init%rest_time = real(temp_end, kind=SP)
    !$omp parallel do private(i)
    do i = 1, ncalc
      st_hydr%read_init(i) = st_hydr%read_init(i)*len_scal_inv
    end do
    !$omp end parallel do

    else if (init_ftype /= in_type(7)) then
      init_fnum = st_init%fnum
      if (len_trim(adjustl(init_unit)) /= 0) then
        if (init_ftype == in_type(3) .or. init_ftype == in_type(5)) then
          if (st_mpi%rank == 0) then
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
          if (ierr /= 0 .and. st_mpi%rank == 0) then
            call write_err_stop("While reading header time in initial file.")
          end if
        end if
#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
          ! -- Bcast scalar value (val)
            call bcast_val(st_init%rest_time, "restart time value")
        end if
#endif
      end if
      if (init_ftype == in_type(3) .or. init_ftype == in_type(4)) then
        call set_2dfile2calc(init_fnum, init_ftype, 0, ncals, DNOVAL, st_hydr%read_init)
      else if (init_ftype == in_type(5) .or. init_ftype == in_type(6)) then
        call set_3dfile2calc(init_fnum, init_ftype, 0, DNOVAL, st_hydr%read_init)
      end if
      !$omp parallel do private(i)
      do i = 1, ncalc
        st_hydr%read_init(i) = st_hydr%read_init(i)*len_scal_inv
      end do
      !$omp end parallel do
      if (init_ftype == in_type(3) .or. init_ftype == in_type(5)) then
        if (st_mpi%rank == 0) then
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
      if (st_mpi%totn /= 1) then
        ! -- Bcast initial depth (init_dep)
          call bcast_init_dep()
      end if
#endif
      ! -- Set initial waterlevel from surface elevation (init_elev)
        call set_init_elev()
    end if

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      sum_init = 1
      allocate(recv_init(ncalc+neib_ncalc))
      !$omp parallel do private(i)
      do i = 1, ncalc+neib_ncalc
        recv_init(i) = DNOVAL
      end do
      !$omp end parallel do
      do while (sum_init /= 0)
        ! -- Send and Receive neighbor value (neibval)
          call senrec_neibval(neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                              recv_citem, st_hydr%read_init, recv_init)
        do i = 1, ncalc
          if (calc2recv(i) > 0 .and. st_hydr%read_init(i) == DNOVAL) then
            st_hydr%read_init(i) = recv_init(calc2recv(i))
            call get_calc_grid(i, xn, yn, zn)
            do k = zn+1, st_grid%nz
              xyzn = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(yn-1) + xn
              j = 0 ; j = gmap_get(st_conn%glo2loc_map, xyzn)
              if (j /= 0 .and. j <= ncalc) then
                st_hydr%read_init(j) = st_hydr%read_init(i)
              end if
            end do
          end if
        end do
        nov_num = count(st_hydr%read_init(:) == DNOVAL)
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
  !*********************************************************************************************
  ! assign_massv -- Assign massbalance value
  !*********************************************************************************************
    ! -- modules
    use open_file, only: inmas_fnum
    use set_condition, only: set_mass2calc
#ifdef MPI_MSG
    use mpi_utility, only: mpimax_val, bcast_val
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
    !-------------------------------------------------------------------------------------------
    allocate(mass_val(ncalc), mass_cflag(ncalc))
    !$omp parallel do private(i)
    do i = 1, ncalc
      mass_val(i) = 0
      mass_cflag(i) = 0
    end do
    !$omp end parallel do

    if (mass_ftype == in_type(1)) then
      if (st_mpi%rank == 0) then
        read(unit=inmas_fnum,fmt=*,iostat=ierr) msout_tnum
        if (ierr /= 0) then
          call write_err_stop("While reading massbalance output number in mass file.")
        else if (msout_tnum <= 0) then
          call write_err_stop("Specified wrong number for massbalance output.")
        end if
      end if
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
         call bcast_val(msout_tnum, "massbalance output number")
      end if
#endif
      allocate(massout_name(msout_tnum))
      allocate(clas_mass(msout_tnum), real_mass(ncalc))
      !$omp parallel
      !$omp do private(i)
      do i = 1, msout_tnum
        massout_name(i) = ""
        clas_mass(i) = SZERO
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, ncalc
        real_mass(i) = SZERO
      end do
      !$omp end do
      !$omp end parallel
      if (st_mpi%rank == 0) then
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
      if (st_mpi%totn /= 1) then
        call bcast_clas_val(msout_tnum, massout_name, clas_mass)
      end if
#endif

      call set_clas2calc(msout_tnum, massout_name, clas_mass, real_mass, mass_cflag, mass_num)
      !$omp parallel do private(i)
      do i = 1, ncalc
        mass_val(i) = nint(real_mass(i))
      end do
      !$omp end parallel do
      deallocate(clas_mass, real_mass)

    else if (mass_ftype /= in_type(7)) then
        massun = inmas_fnum
      if (mass_ftype == in_type(3) .or. mass_ftype == in_type(4)) then
        call set_2dfile2calc(massun, mass_ftype, 0, ncals, 0, mass_val, mass_cflag, mass_num)
      else if (mass_ftype == in_type(5) .or. mass_ftype == in_type(6)) then
        call set_3dfile2calc(massun, mass_ftype, 0, 0, mass_val, mass_cflag, mass_num)
      end if

      if (mass_ftype == in_type(3) .or. mass_ftype == in_type(5)) then
        if (st_mpi%rank == 0) then
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
      !$omp parallel do private(i)
      do i = 1, ncalc
        mass_val(i) = 1 ; mass_cflag(i) = 1
      end do
      !$omp end parallel do
    end if

    allocate(all_mass_type(6))
    all_mass_type(:) = [in_type(1), in_type(3:7)]

    if (any(mass_ftype == all_mass_type(:))) then
      allocate(mass2calc(mass_num), calc_mass(mass_num))
      !$omp parallel do private(i)
      do i = 1, mass_num
        mass2calc(i) = 0 ; calc_mass(i) = 0
      end do
      !$omp end parallel do
      call set_mass2calc(ncalc, mass_cflag, mass_val, mass2calc, calc_mass)
      allocate(int_mass(mass_num))
      !$omp parallel do private(i)
      do i = 1, mass_num
        int_mass(i) = 0
        int_mass(i) = calc_mass(i)
      end do
      !$omp end parallel do
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
  !*********************************************************************************************
  ! set_init_elev -- Set initial waterlevel from surface elevation
  !*********************************************************************************************
    ! -- modules
    use set_cell, only: get_cals_grid
    use make_cell, only: st_geom
    ! -- inout

    ! -- local
    integer(I4) :: i, j, k, xn, yn, xyzn
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, xn, yn, xyzn)
    do i = 1, ncals
      call get_cals_grid(i, xn, yn)
      do k = 1, st_grid%nz
        xyzn = (st_grid%nx*st_grid%ny)*(k-1) + st_grid%nx*(yn-1) + xn
        j = 0 ; j = gmap_get(st_conn%glo2loc_map, xyzn)
        if (j /= 0 .and. j <= ncalc) then
          st_hydr%read_init(j) = st_geom%surf_elev(i) - st_init%depth
        end if
      end do
    end do
    !$omp end parallel do

  end subroutine set_init_elev

end module assign_calc
