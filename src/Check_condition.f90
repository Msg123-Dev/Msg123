module check_condition
  ! -- modules
  use kind_module, only: I4, SP
  use constval_module, only: SONE, SNOVAL
  use utility_module, only: write_err_stop, conv_i2s
  use initial_module, only: my_rank, st_grid, st_in_type, st_seal, in_type, out_type,&
                            st_out_type, st_out_step, st_out_path, st_out_time,&
                            st_out_unit
  use open_file, only: open_out_binf
  use read_input, only: len_scal, len_scal_inv
#ifdef MPI_MSG
  use mpi_utility, only: mpisum_val
  use mpi_read, only: set_real4_fview
  use mpi_set, only: write_2dview, write_3dview
#endif

  implicit none
  private
  public :: check_calc_region, check_calc_retn, check_calc_parm, check_calc_init
  public :: check_outf_cond

  type :: fnum_out
    integer(I4) :: conv, head, rest, srat, wtab, mass, velx, vely, velz
    integer(I4) :: rivr, lakr, sufr, dunr, seal, well, rech, calg
  end type fnum_out
  type(fnum_out), public :: st_out_fnum

  ! -- local
  integer(I4) :: intse_type

  contains

  subroutine check_calc_region(gclas_flag, greg_flag)
  !***************************************************************************************
  ! check_calc_region -- Check calculation region
  !***************************************************************************************
    ! -- modules
    use constval_module, only: SZERO
    use initial_module, only: st_in_path, st_in_unit
    ! -- inout
    integer(I4), intent(in) :: gclas_flag(:,:)
    integer(I4), intent(inout) :: greg_flag(:)
    ! -- local
    integer(I4) :: num_seal, num_calc
    integer(I4), allocatable :: temp_greg(:)
    real(SP), allocatable :: check_seal(:), temp_seal(:)
    logical, allocatable :: mask(:)
    !-------------------------------------------------------------------------------------
    if (st_in_type%seal > 0) then
      ! -- open input sea level file for check (sealf)
        call open_sealf(st_in_type%seal, st_in_path%seal, st_in_unit%seal)
    end if

    if (st_seal%totn > 0) then
      allocate(check_seal(st_grid%nxyz))
      !$omp parallel workshare
      check_seal(:) = SNOVAL
      !$omp end parallel workshare
      call read_sealf(st_in_type%seal, st_seal%totn, gclas_flag, check_seal)

      num_seal = count(greg_flag(:) == 0)
      allocate(temp_seal(num_seal), mask(st_grid%nxyz))
      !$omp parallel workshare
      temp_seal(:) = SZERO ; mask(:) = (greg_flag(:) == 0)
      temp_seal(:) = pack(check_seal(:), mask(:))
      !$omp end parallel workshare
      if (any((temp_seal(:)==SNOVAL))) then
        call write_err_stop("Null value in sea region.")
      end if

      deallocate(temp_seal)

      num_calc = count(greg_flag(:) > 0)
      allocate(temp_greg(num_calc))
      !$omp parallel workshare
      temp_greg(:) = 0
      mask(:) = (greg_flag(:) > 0) .and. (check_seal(:) /= SNOVAL)
      greg_flag(:) = unpack(temp_greg(:), mask(:), greg_flag(:))
      !$omp end parallel workshare

      deallocate(temp_greg, check_seal, mask)

    end if

  end subroutine check_calc_region

  subroutine check_calc_retn(calcn, calc_reta, calc_retn, calc_resi)
  !***************************************************************************************
  ! check_calc_retn -- Check calculation retention value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: calcn
    real(SP), intent(in) :: calc_reta(:), calc_retn(:), calc_resi(:)
    ! -- local
    integer(I4) :: check_flag, sum_checkf
    character(:), allocatable :: str_num, err_mes
    !-------------------------------------------------------------------------------------
    check_flag = 0 ; sum_checkf = 0

    ! -- Check read value (read_val)
      call check_read_val(calcn, calc_reta, len_scal_inv, check_flag)
      call check_read_val(calcn, calc_retn, SONE, check_flag)
      call check_read_val(calcn, calc_resi, SONE, check_flag)
#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(check_flag, "check retention", sum_checkf)
#else
    sum_checkf = check_flag
#endif
    if (sum_checkf > 0 .and. my_rank == 0) then
      str_num = conv_i2s(sum_checkf)
      err_mes = "There are "//str_num//" cells which have no retention value."
      call write_err_stop(err_mes)
      deallocate(err_mes)
    end if

  end subroutine check_calc_retn

  subroutine check_calc_parm(calcn, calc_ksx, calc_ksy, calc_ksz, calc_ss, calc_por)
  !***************************************************************************************
  ! check_calc_parm -- Check calculation parameter value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: calcn
    real(SP), intent(in) :: calc_ksx(:), calc_ksy(:), calc_ksz(:), calc_ss(:), calc_por(:)
    ! -- local
    integer(I4) :: check_flag, sum_checkf
    character(:), allocatable :: str_num, err_mes
    !-------------------------------------------------------------------------------------
    check_flag = 0 ; sum_checkf = 0
    ! -- Check read value (read_val)
      call check_read_val(calcn, calc_ksx, len_scal, check_flag)
      call check_read_val(calcn, calc_ksy, len_scal, check_flag)
      call check_read_val(calcn, calc_ksz, len_scal, check_flag)
      call check_read_val(calcn, calc_ss, len_scal_inv, check_flag)
      call check_read_val(calcn, calc_por, SONE, check_flag)
#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(check_flag, "check parameter", sum_checkf)
#else
    sum_checkf = check_flag
#endif
    if (sum_checkf > 0 .and. my_rank == 0) then
      str_num = conv_i2s(sum_checkf)
      err_mes = "There are "//str_num//" cells which have no parameter value."
      call write_err_stop(err_mes)
      deallocate(err_mes)
    end if

  end subroutine check_calc_parm

  subroutine check_calc_init(calcn, calc_init)
  !***************************************************************************************
  ! check_calc_init -- Check calculation initial value
  !***************************************************************************************
    ! -- modules
    use kind_module, only: DP
    ! -- inout
    integer(I4), intent(in) :: calcn
    real(DP), intent(in) :: calc_init(:)
    ! -- local
    integer(I4) :: check_flag, sum_checkf
    character(:), allocatable :: str_num, err_mes
    !-------------------------------------------------------------------------------------
    check_flag = 0 ; sum_checkf = 0
    ! -- Check read value (read_val)
      call check_read_val(calcn, real(calc_init, kind=SP), len_scal, check_flag)
#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(check_flag, "check initial", sum_checkf)
#else
    sum_checkf = check_flag
#endif
    if (sum_checkf > 0 .and. my_rank == 0) then
      str_num = conv_i2s(sum_checkf)
      err_mes = "There are "//str_num//" cells which have no initial value."
      call write_err_stop(err_mes)
      deallocate(err_mes)
    end if

  end subroutine check_calc_init

  subroutine open_sealf(seal_ftype, seal_path, seal_unit)
  !***************************************************************************************
  ! open_sealf -- Open input sea level file for check
  !***************************************************************************************
    ! -- modules
    use kind_module, only: SP
    use utility_module, only: open_new_rtxt, open_new_rbin, write_err_read, write_success
    ! -- inout
    integer(I4), intent(in) :: seal_ftype
    character(*), intent(in) :: seal_path, seal_unit
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: seal_fnum, intse_num
    integer(I4), allocatable :: type_txt(:), type_bin(:)
    real(SP) :: intse_step, intse_end
    character(:), allocatable :: err_mes, intsep
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(type_txt(4), type_bin(2))
    type_txt(:) = [in_type(1:3), in_type(5)] ; type_bin(:) = [in_type(4), in_type(6)]
    err_mes = "input sea level for check" ; intsep = ""
    if (any(seal_ftype == type_txt(:))) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, seal_path, err_mes, seal_fnum)
      st_seal%fnum = seal_fnum
    else if (any(seal_ftype == type_bin(:))) then
      ! -- Open new read binary file (new_rbin)
        call open_new_rbin(1, 1, seal_path, err_mes, seal_fnum)
      st_seal%fnum = seal_fnum
    else if (seal_ftype == in_type(7)) then
      ! -- Open new read text file (new_rtxt)
        call open_new_rtxt(1, 1, seal_path, err_mes, seal_fnum)
      read(seal_fnum,*) intse_type, intse_step, intse_end
      read(seal_fnum,'(a)') intsep
      err_mes = "input sea level time interval for check"
      if (intse_type == in_type(3) .or. intse_type == in_type(5)) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(1, 1, seal_path, err_mes, intse_num)
        st_seal%fnum = intse_num
      else if (any(intse_type == type_bin(:))) then
        ! -- Open new read binary file (new_rbin)
          call open_new_rbin(1, 1, seal_path, err_mes, intse_num)
        st_seal%fnum = intse_num
      else
        call write_err_stop("Specified wrong sea level time interval number.")
      end if
      deallocate(intsep)
    end if

    deallocate(err_mes)

    st_seal%totn = st_grid%nxyz

    if (len_trim(adjustl(seal_unit)) /= 0) then
      if (seal_ftype == in_type(1) .or. seal_ftype == in_type(2)) then
        read(unit=st_seal%fnum,fmt=*,iostat=ierr) st_seal%totn, st_seal%etime
      else if (seal_ftype == in_type(3) .or. seal_ftype == in_type(5)) then
        read(unit=st_seal%fnum,fmt=*,iostat=ierr) st_seal%etime
      else if (any(seal_ftype == type_bin(:))) then
        read(unit=st_seal%fnum,fmt=*,iostat=ierr) st_seal%etime
      end if
    else if (seal_ftype == in_type(1) .or. seal_ftype == in_type(2)) then
      read(unit=st_seal%fnum,fmt=*,iostat=ierr) st_seal%totn
    end if

    if (ierr /= 0) then
      call write_err_read(st_seal%fnum)
    end if

    deallocate(type_txt, type_bin)

  end subroutine open_sealf

  subroutine read_sealf(seal_ftype, seal_num, seal_clas, seal_val)
  !***************************************************************************************
  ! read_sealf -- Read sea level value for check
  !***************************************************************************************
    ! -- modules
    use constval_module, only: VARLEN
    use utility_module, only: close_file
    use read_module, only: read_clasf, read_3dpointf
    ! -- inout
    integer(I4), intent(in) :: seal_ftype, seal_num
    integer(I4), intent(in) :: seal_clas(:,:)
    real(SP), intent(out) :: seal_val(:)
    ! -- local
    integer(I4) :: seal_fnum
    integer(I4), allocatable :: seal_i(:), seal_j(:), seal_k(:)
    integer(I4), allocatable :: type_2d(:), type_3d(:)
    real(SP), allocatable :: seal_read(:)
    character(VARLEN), allocatable :: seal_name(:)
    !-------------------------------------------------------------------------------------
    allocate(type_2d(2), type_3d(2))
    type_2d(:) = [in_type(3:4)] ; type_3d(:) = [in_type(5:6)]
    seal_fnum = st_seal%fnum

    if (seal_ftype == in_type(1)) then
      allocate(seal_read(seal_num), seal_name(seal_num))
      !$omp parallel workshare
      seal_read(:) = SNOVAL ; seal_name(:) = ""
      !$omp end parallel workshare
      call read_clasf(seal_fnum, seal_num, seal_name, seal_read)
      call close_file(seal_fnum)

      call set_clas2seal(seal_fnum, seal_name, seal_clas, seal_read, seal_val)
      deallocate(seal_read, seal_name)
    else if (seal_ftype == in_type(2)) then
      allocate(seal_read(seal_num))
      allocate(seal_i(seal_num), seal_j(seal_num), seal_k(seal_num))
      !$omp parallel workshare
      seal_read(:) = SNOVAL
      seal_i(:) = 0 ; seal_j(:) = 0 ; seal_k(:) = 0
      !$omp end parallel workshare
      call read_3dpointf(seal_fnum, seal_num, seal_i, seal_j, seal_k, seal_read)
      call close_file(seal_fnum)

      call set_point2seal(seal_num, seal_i, seal_j, seal_k, seal_read, seal_val)
      deallocate(seal_i, seal_j, seal_k, seal_read)
    else if (any(seal_ftype == type_2d(:))) then
      call set_2dfile2seal(seal_fnum, seal_ftype, SNOVAL, seal_val)
      call close_file(seal_fnum)
    else if (any(seal_ftype == type_3d(:))) then
      call set_3dfile2seal(seal_fnum, seal_ftype, SNOVAL, seal_val)
      call close_file(seal_fnum)
    else if (seal_ftype == in_type(7)) then
      if (any(seal_ftype == type_2d(:))) then
        call set_2dfile2seal(seal_fnum, intse_type, SNOVAL, seal_val)
      else if (any(seal_ftype == type_3d(:))) then
        call set_3dfile2seal(seal_fnum, intse_type, SNOVAL, seal_val)
      end if
      call close_file(seal_fnum)
    end if

  end subroutine read_sealf

  subroutine check_outf_cond()
  !***************************************************************************************
  ! check_outf_cond -- Check output file condition
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    ! -- Check convergence file (convf)
      call check_convf()

    ! -- Check head file (headf)
      call check_headf()

    ! -- Check restart file (restf)
      call check_restf()

    ! -- Check saturation file (sratf)
      call check_sratf()

    ! -- Check water table file (wtabf)
      call check_wtabf()

    ! -- Check massbalance file (massf)
      call check_massf()

    ! -- Check velocity file (velcf)
      call check_velcf()

    ! -- Check river runoff file (rivrf)
      call check_rivrf()

    ! -- Check lake runoff file (lakrf)
      call check_lakrf()

    ! -- Check surface runoff file (sufrf)
      call check_sufrf()

    ! -- Check dunne runoff file (dunrf)
      call check_dunrf()

    ! -- Check output sea results file (out_sealf)
      call check_out_sealf()

    ! -- Check output recharge results file (out_rechf)
      call check_out_rechf()

    ! -- Check output well pumping results file (out_wellf)
      call check_out_wellf()

  end subroutine check_outf_cond

  subroutine check_read_val(readn, readv, cunit, cflag)
  !***************************************************************************************
  ! check_read_val -- Check read value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: readn
    real(SP), intent(in) :: readv(:)
    real(SP), intent(in) :: cunit
    integer(I4), intent(inout) :: cflag
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i) reduction(+:cflag)
    do i = 1, readn
      if (readv(i)*cunit == SNOVAL) then
        cflag = cflag + 1
      end if
    end do
    !$omp end parallel do

  end subroutine check_read_val

  subroutine set_clas2seal(clasn, clas_name, clas_set, clas_val, tg_val)
  !***************************************************************************************
  ! set_clas2seal -- Set sea level from classification
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_clas
    ! -- inout
    integer(I4), intent(in) :: clasn
    character(*), intent(in) :: clas_name(:)
    integer(I4), intent(in) :: clas_set(:,:)
    real(SP), intent(in) :: clas_val(:)
    real(SP), intent(out) :: tg_val(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4), allocatable :: temp_clas(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_clas(clasn))
    !$omp parallel workshare
    temp_clas(:) = 0
    !$omp end parallel workshare

    do i = 1, clasn
      do k = 1, st_clas%totn
        if (clas_name(i) == st_clas%name(k)) then
          temp_clas(i) = k
        end if
      end do
    end do

    do i = 1, clasn
      k = temp_clas(i)
      !$omp parallel do private(j)
      do j = 1, st_grid%nxyz
        if (clas_set(j,k) == 1) then
          tg_val(j) = clas_val(i)
        end if
      end do
      !$omp end parallel do
    end do

    deallocate(temp_clas)

  end subroutine set_clas2seal

  subroutine set_point2seal(pointn, p_i, p_j, p_k, poin_val, tg_val)
  !***************************************************************************************
  ! set_point2seal -- Set sea level from point
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: pointn
    integer(I4), intent(in) :: p_i(:), p_j(:), p_k(:)
    real(SP), intent(in) :: poin_val(:)
    real(SP), intent(out) :: tg_val(:)
    ! -- local
    integer(I4) :: i, c_num, ii, jj, kk
    integer(I4) :: pi, pj, pk, ist, ien, jst, jen, kst, ken
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, ii, jj, kk, pi, pj, pk, ist, ien, jst, jen, kst, ken, c_num)
    do i = 1, pointn
      pi = p_i(i) ; pj = p_j(i) ; pk = p_k(i)
      if (pi == -1 .and. pj == -1 .and. pk == -1) then !all cells
        ist = 1 ; ien = st_grid%nx ; jst = 1 ; jen = st_grid%ny ; kst = 1 ; ken = st_grid%nz
      else if (pi == -1 .and. pj == -1) then !i,j cells
        ist = 1 ; ien = st_grid%nx ; jst = 1 ; jen = st_grid%ny ; kst = pk ; ken = pk
      else if (pi == -1 .and. pk == -1) then !i,k cells
        ist = 1 ; ien = st_grid%nx ; jst = pj ; jen = pj ; kst = 1 ; ken = st_grid%nz
      else if (pj == -1 .and. pk == -1) then !j,k cells
        ist = pi ; ien = pi ; jst = 1 ; jen = st_grid%ny ; kst = 1 ; ken = st_grid%nz
      else if (pi == -1) then !only i cell
        ist = 1 ; ien = st_grid%nx ; jst = pj ; jen = pj ; kst = pk ; ken = pk
      else if (pj == -1) then !only j cell
        ist = pi ; ien = pi ; jst = 1 ; jen = st_grid%ny ; kst = pk ; ken = pk
      else if (pk == -1) then !only k cell
        ist = pi ; ien = pi ; jst = pj ; jen = pj ; kst = 1 ; ken = st_grid%nz
      else !others
        ist = pi ; ien = pi ; jst = pj ; jen = pj ; kst = pk ; ken = pk
      end if

      do ii = ist, ien
        do jj = jst, jen
          do kk = kst, ken
            c_num = st_grid%nx*(st_grid%ny*(kk-1) + (jj-1)) + ii
            tg_val(c_num) = poin_val(i)
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine set_point2seal

  subroutine set_2dfile2seal(fnum, file_type, no_val, tg_val)
  !***************************************************************************************
  ! set_2dfile2seal -- Set sea level value from 2d file
  !***************************************************************************************
    ! -- modules
    use read_module, only: read_2dtxt, read_2dbin, flat_2dto3d
    ! -- inout
    integer(I4), intent(in) :: fnum, file_type
    real(SP), intent(in) :: no_val
    real(SP), intent(out) :: tg_val(:)
    ! -- local
    real(SP), allocatable :: array_read(:,:)
    !-------------------------------------------------------------------------------------
    allocate(array_read(st_grid%nx,st_grid%ny))
    !$omp parallel workshare
    array_read(:,:) = no_val
    !$omp end parallel workshare

    if (file_type == in_type(3)) then
      call read_2dtxt(fnum, st_grid%nx, st_grid%ny, array_read)
    else if (file_type == in_type(4)) then
      call read_2dbin(fnum, st_grid%nx, st_grid%ny, no_val, array_read)
    end if

    call flat_2dto3d(st_grid%nx, st_grid%ny, st_grid%nz, array_read, tg_val)
    deallocate(array_read)

  end subroutine set_2dfile2seal

  subroutine set_3dfile2seal(fnum, file_type, no_val, tg_val)
  !***************************************************************************************
  ! set_3dfile2seal -- Set sea level value from 3d file
  !***************************************************************************************
    ! -- modules
    use read_module, only: read_3dtxt, read_3dbin, flat_3dto3d
    ! -- inout
    integer(I4), intent(in) :: fnum, file_type
    real(SP), intent(in) :: no_val
    real(SP), intent(out) :: tg_val(:)
    ! -- local
    real(SP), allocatable :: array_read(:,:,:)
    !-------------------------------------------------------------------------------------
    allocate(array_read(st_grid%nx,st_grid%ny,st_grid%nz))
    !$omp parallel workshare
    array_read(:,:,:) = no_val
    !$omp end parallel workshare

    if (file_type == in_type(5)) then
      call read_3dtxt(fnum, st_grid%nx, st_grid%ny, st_grid%nz, array_read)
    else if (file_type == in_type(6)) then
      call read_3dbin(fnum, st_grid%nx, st_grid%ny, st_grid%nz, no_val, array_read)
    end if
    call flat_3dto3d(st_grid%nx, st_grid%ny, st_grid%nz, array_read, tg_val)

    deallocate(array_read)

  end subroutine set_3dfile2seal

  subroutine check_convf()
  !***************************************************************************************
  ! check_convf -- Check convergence file
  !***************************************************************************************
    ! -- modules
    use open_file, only: open_out_convf
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    if (my_rank == 0) then
      ! -- Open output convergence file (out_convf)
        call open_out_convf(st_out_path%conv, st_out_fnum%conv)
    end if

  end subroutine check_convf

  subroutine check_headf()
  !***************************************************************************************
  ! check_headf -- Check head file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: head_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    head_file = 0 ; out_mess = "output head"
    ! -- Open output binary file (out_binf)
      call open_out_binf(head_file, st_out_time%head, st_out_path%head, st_out_unit%head,&
                         out_mess, st_out_step%head)

#ifdef MPI_MSG
    ! -- Set real4 file view (real4_fview)
      call set_real4_fview(head_file, write_3dview, out_mess)
#endif

    st_out_fnum%head = head_file

    deallocate(out_mess)

  end subroutine check_headf

  subroutine check_restf()
  !***************************************************************************************
  ! check_restf -- Check restart file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_read, only: set_real8_fview
    use mpi_set, only: rest_view
#endif
    ! -- inout

    ! -- local
    integer(I4) :: rest_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    rest_file = 0 ; out_mess = "output restart"
    ! -- Open output binary file (out_binf)
      call open_out_binf(rest_file, st_out_time%rest, st_out_path%rest, st_out_unit%rest,&
                         out_mess, st_out_step%rest)

#ifdef MPI_MSG
    ! -- Set real8 file view (real8_fview)
      call set_real8_fview(rest_file, rest_view, out_mess)
#endif

    st_out_fnum%rest = rest_file

    deallocate(out_mess)

  end subroutine check_restf

  subroutine check_sratf()
  !***************************************************************************************
  ! check_sratf -- Check saturation file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: srat_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%srat == out_type(3)) then
      srat_file = 0 ; out_mess = "output saturation"
      ! -- Open output binary file (out_binf)
        call open_out_binf(srat_file, st_out_time%srat, st_out_path%srat,&
                           st_out_unit%srat, out_mess, st_out_step%srat)

#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(srat_file, write_3dview, out_mess)
#endif

      st_out_fnum%srat = srat_file

      deallocate(out_mess)
    end if

  end subroutine check_sratf

  subroutine check_wtabf()
  !***************************************************************************************
  ! check_wtabf -- Check water table file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: wtab_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%wtab == out_type(2) .and. st_in_type%wtab == in_type(7)) then
      wtab_file = 0 ; out_mess = "output water table"
      ! -- Open output binary file (out_binf)
        call open_out_binf(wtab_file, st_out_time%wtab, st_out_path%wtab,&
                           st_out_unit%wtab, out_mess, st_out_step%wtab)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(wtab_file, write_2dview, out_mess)
#endif
      st_out_fnum%wtab = wtab_file
      deallocate(out_mess)
    end if

  end subroutine check_wtabf

  subroutine check_massf()
  !***************************************************************************************
  ! check_massf -- Check massbalance file
  !***************************************************************************************
    ! -- modules
    use utility_module, only: write_logf
    use open_file, only: open_out_massf
    ! -- inout

    ! -- local
    integer(I4) :: mass_file
    integer(I4), allocatable :: all_mass(:)
    !-------------------------------------------------------------------------------------
    allocate(all_mass(6))
    all_mass(:) = [in_type(1), in_type(3:7)]

    if (any(st_in_type%mass == all_mass(:)) .and. st_out_type%mass == out_type(1)) then
      mass_file = 0
      ! -- Open output massbalance file (out_massf)
        call open_out_massf(st_out_time%mass, st_out_path%mass, st_out_unit%mass, mass_file)
        st_out_fnum%mass = mass_file
    else if (my_rank == 0) then
      call write_logf("Set not to use massbalance output function.")
    end if

    deallocate(all_mass)

  end subroutine check_massf

  subroutine check_velcf()
  !***************************************************************************************
  ! check_velcf -- Check velocity file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: velc_file
    character(:), allocatable :: out_mess, new_velc
    !-------------------------------------------------------------------------------------
    if (st_out_type%velc == out_type(3)) then
      out_mess = "output velocity x direction" ; new_velc = st_out_path%velx
      ! -- Open output binary file (out_binf)
        call open_out_binf(velc_file, st_out_time%velc, new_velc, st_out_unit%velc,&
                           out_mess, st_out_step%velc)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(velc_file, write_3dview, out_mess)
#endif
      st_out_fnum%velx = velc_file

      out_mess = "output velocity y direction" ; new_velc = st_out_path%vely
      ! -- Open output binary file (out_binf)
        call open_out_binf(velc_file, st_out_time%velc, new_velc, st_out_unit%velc,&
                           out_mess, st_out_step%velc)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(velc_file, write_3dview, out_mess)
#endif
      st_out_fnum%vely = velc_file

      out_mess = "output velocity z direction" ; new_velc = st_out_path%velz
      ! -- Open output binary file (out_binf)
        call open_out_binf(velc_file, st_out_time%velc, new_velc, st_out_unit%velc,&
                           out_mess, st_out_step%velc)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(velc_file, write_3dview, out_mess)
#endif

      st_out_fnum%velz = velc_file

      deallocate(out_mess, new_velc)
    end if

  end subroutine check_velcf

  subroutine check_rivrf()
  !***************************************************************************************
  ! check_rivrf -- Check river runoff file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: rivr_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%rivr == out_type(2)) then
      rivr_file = 0 ; out_mess = "output river runoff"
      ! -- Open output binary file (out_binf)
        call open_out_binf(rivr_file, st_out_time%rivr, st_out_path%rivr,&
                           st_out_unit%rivr, out_mess, st_out_step%rivr)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(rivr_file, write_2dview, out_mess)
#endif
      st_out_fnum%rivr = rivr_file

      deallocate(out_mess)
    end if

  end subroutine check_rivrf

  subroutine check_lakrf()
  !***************************************************************************************
  ! check_lakrf -- Check lake runoff file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: lakr_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%lakr == out_type(2)) then
      lakr_file = 0 ; out_mess = "output lake runoff"
      ! -- Open output binary file (out_binf)
        call open_out_binf(lakr_file, st_out_time%lakr, st_out_path%lakr,&
                           st_out_unit%lakr, out_mess, st_out_step%lakr)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(lakr_file, write_2dview, out_mess)
#endif
      st_out_fnum%lakr = lakr_file

      deallocate(out_mess)
    end if

  end subroutine check_lakrf

  subroutine check_sufrf()
  !***************************************************************************************
  ! check_sufrf -- Check surface runoff file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: sufr_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%sufr == out_type(2)) then
      sufr_file = 0 ; out_mess = "output surface runoff"
      ! -- Open output binary file (out_binf)
        call open_out_binf(sufr_file, st_out_time%sufr, st_out_path%sufr,&
                           st_out_unit%sufr, out_mess, st_out_step%sufr)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(sufr_file, write_2dview, out_mess)
#endif
      st_out_fnum%sufr = sufr_file

      deallocate(out_mess)
    end if

  end subroutine check_sufrf

  subroutine check_dunrf()
  !***************************************************************************************
  ! check_dunrf -- Check dunne runoff file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: dunr_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%dunr == out_type(3)) then
      dunr_file = 0 ; out_mess = "output dunne runoff"
      ! -- Open output binary file (out_binf)
        call open_out_binf(dunr_file, st_out_time%dunr, st_out_path%dunr,&
                           st_out_unit%dunr, out_mess, st_out_step%dunr)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(dunr_file, write_2dview, out_mess)
#endif
      st_out_fnum%dunr = dunr_file

      deallocate(out_mess)
    end if

  end subroutine check_dunrf

  subroutine check_out_sealf()
  !***************************************************************************************
  ! check_out_sealf -- Check sea results file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: seal_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%seal == out_type(3)) then
      seal_file = 0 ; out_mess = "output sea"
      ! -- Open output binary file (out_binf)
        call open_out_binf(seal_file, st_out_time%seal, st_out_path%seal,&
                           st_out_unit%seal, out_mess, st_out_step%seal)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(seal_file, write_3dview, out_mess)
#endif
      st_out_fnum%seal = seal_file

      deallocate(out_mess)
    end if

  end subroutine check_out_sealf

  subroutine check_out_rechf()
  !***************************************************************************************
  ! check_out_rechf -- Check output recharge results file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: rech_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%rech == out_type(2)) then
      rech_file = 0 ; out_mess = "output recharge"
      ! -- Open output binary file (out_binf)
        call open_out_binf(rech_file, st_out_time%rech, st_out_path%rech,&
                           st_out_unit%rech, out_mess, st_out_step%rech)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(rech_file, write_2dview, out_mess)
#endif
      st_out_fnum%rech = rech_file

      deallocate(out_mess)
    end if

  end subroutine check_out_rechf

  subroutine check_out_wellf()
  !***************************************************************************************
  ! check_out_wellf -- Check output well pumping results file
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: well_file
    character(:), allocatable :: out_mess
    !-------------------------------------------------------------------------------------
    if (st_out_type%well == out_type(3)) then
      well_file = 0 ; out_mess = "output well"
      ! -- Open output binary file (out_binf)
        call open_out_binf(well_file, st_out_time%well, st_out_path%well,&
                           st_out_unit%well, out_mess, st_out_step%well)
#ifdef MPI_MSG
      ! -- Set real4 file view (real4_fview)
        call set_real4_fview(well_file, write_3dview, out_mess)
#endif
      st_out_fnum%well = well_file

      deallocate(out_mess)
    end if

  end subroutine check_out_wellf

end module check_condition
