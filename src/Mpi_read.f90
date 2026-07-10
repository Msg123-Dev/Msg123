module mpi_read
  ! -- modules
  use kind_module, only: I4, SP, DP
  use utility_module, only: st_mpi, write_err_read, write_err_stop
  use initial_module, only: st_sim, st_grid, st_init, in_type
  use mpi

  implicit none
  private
  public :: open_mpi_read_file, open_mpi_write_file
  public :: set_int4_fview, set_real4_fview, set_real8_fview
  public :: read_mpi_restf, read_mpi_head, read_mpi_file
  public :: skip_mpi_file, skip_mpi_file_int, close_mpi_file

  interface read_mpi_head
    module procedure read_mpi_i4head
    module procedure read_mpi_r4head
    module procedure read_mpi_r8head
  end interface

  interface read_mpi_file
    module procedure read_mpi_i4
    module procedure read_mpi_r4
    module procedure read_mpi_r8
  end interface

  ! -- local
  integer(I4) :: gview2d, gview3d

  contains

  subroutine open_mpi_read_file(stop_flag, write_flag, mpi_path, err_mes, mpi_fh, mpi_ier)
  !*********************************************************************************************
  ! open_mpi_read_file -- Open mpi read file
  !*********************************************************************************************
    ! -- module
    use utility_module, only: write_success
    ! -- inout
    integer(I4), intent(in) :: stop_flag, write_flag
    character(*), intent(in) :: mpi_path, err_mes
    integer(I4), intent(out) :: mpi_fh
    integer(I4), intent(out), optional :: mpi_ier
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_FILE_OPEN(st_mpi%comm, mpi_path, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_fh, ierr)

    if (ierr == MPI_SUCCESS) then
      if (st_mpi%rank == 0 .and. write_flag == 1) then
        call write_success("Open "//err_mes//" file", mpi_fh)
      end if
    else if (stop_flag == 1) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Open "//err_mes//" file.")
      end if
    end if

    if (present(mpi_ier)) then
      mpi_ier = ierr
    end if

  end subroutine open_mpi_read_file

  subroutine open_mpi_write_file(mpi_path, err_mes, mpi_fh)
  !*********************************************************************************************
  ! open_mpi_write_file -- Open mpi write file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: mpi_path, err_mes
    integer(I4), intent(out) :: mpi_fh
    ! -- local
    integer(I4) :: ierr
    integer(KIND=MPI_OFFSET_KIND) :: offset
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_FILE_OPEN(st_mpi%comm, mpi_path, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL,&
                       mpi_fh, ierr)

    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Open "//err_mes//" file.")
      end if
    end if

    offset = 0
    call MPI_FILE_SET_SIZE(mpi_fh, offset, ierr)

  end subroutine open_mpi_write_file

  subroutine set_int4_fview(fileh, file_view, err_mes)
  !*********************************************************************************************
  ! set_int4_fview -- Set int4 file view
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: fileh, file_view
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; head_dis = 0
    call MPI_FILE_SET_VIEW(fileh, head_dis, MPI_INTEGER, file_view, "native", MPI_INFO_NULL,&
                           ierr)

    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Set view "//err_mes//" file.")
      end if
    end if

  end subroutine set_int4_fview

  subroutine set_real4_fview(fileh, file_view, err_mes)
  !*********************************************************************************************
  ! set_real4_fview -- Set real4 file view
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: fileh, file_view
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; head_dis = 0
    call MPI_FILE_SET_VIEW(fileh, head_dis, MPI_REAL4, file_view, "native", MPI_INFO_NULL,&
                           ierr)

    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Set view "//err_mes//" file.")
      end if
    end if

  end subroutine set_real4_fview

  subroutine set_real8_fview(fileh, file_view, err_mes)
  !*********************************************************************************************
  ! set_real8_fview -- Set real8 file view
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: fileh, file_view
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; head_dis = 0
    call MPI_FILE_SET_VIEW(fileh, head_dis, MPI_REAL8, file_view, "native", MPI_INFO_NULL,&
                           ierr)

    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Set view "//err_mes//" file.")
      end if
    end if

  end subroutine set_real8_fview

  subroutine read_mpi_restf(fileh, calc_num, rest_time, calc_init)
  !*********************************************************************************************
  ! read_mpi_restf -- Read mpi restart file
  !*********************************************************************************************
    ! -- module
    use constval_module, only: DZERO
    ! -- inout
    integer(I4), intent(in) :: fileh, calc_num
    real(DP), intent(out) :: rest_time
    real(DP), intent(out) :: calc_init(:)
    ! -- local
    integer(I4) :: i, ierr
    integer(I4), allocatable :: istat(:)
    real(DP), allocatable :: read_rest(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(istat(MPI_STATUS_SIZE), read_rest(calc_num+1))
    !$omp parallel
    !$omp do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, calc_num+1
      read_rest(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    call MPI_FILE_READ_ALL(fileh, read_rest, calc_num+1, MPI_REAL8, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Read restart file in MPI program.")
      end if
    end if

    rest_time = read_rest(1)

    !$omp parallel do private(i)
    do i = 1, calc_num
      calc_init(i) = read_rest(i+1)
    end do
    !$omp end parallel do

    deallocate(istat, read_rest)

  end subroutine read_mpi_restf

  subroutine read_mpi_i4head(fileh, ierr, head_out)
  !*********************************************************************************************
  ! read_mpi_i4head -- Read mpi integer header
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: fileh
    integer(I4), intent(out) :: ierr
    integer(I4), intent(out) :: head_out
    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: istat(:)
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------------
    allocate(istat(MPI_STATUS_SIZE))
    !$omp parallel do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end parallel do

    call MPI_FILE_GET_POSITION(fileh, head_dis, ierr)

    call MPI_FILE_READ_AT_ALL(fileh, head_dis, head_out, 1, MPI_INTEGER, istat, ierr)

    deallocate(istat)

  end subroutine read_mpi_i4head

  subroutine read_mpi_r4head(fileh, ierr, head_out)
  !*********************************************************************************************
  ! read_mpi_r4head -- Read mpi real4 header
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: fileh
    integer(I4), intent(out) :: ierr
    real(SP), intent(out) :: head_out
    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: istat(:)
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------------
    allocate(istat(MPI_STATUS_SIZE))
    !$omp parallel do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end parallel do

    call MPI_FILE_GET_POSITION(fileh, head_dis, ierr)

    call MPI_FILE_READ_AT_ALL(fileh, head_dis, head_out, 1, MPI_REAL4, istat, ierr)

    deallocate(istat)

  end subroutine read_mpi_r4head

  subroutine read_mpi_r8head(fileh, ierr, head_out)
  !*********************************************************************************************
  ! read_mpi_r8head -- Read mpi real8 header
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: fileh
    integer(I4), intent(out) :: ierr
    real(DP), intent(out) :: head_out
    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: istat(:)
    integer(KIND=MPI_OFFSET_KIND) :: head_dis
    !-------------------------------------------------------------------------------------------
    allocate(istat(MPI_STATUS_SIZE))
    !$omp parallel do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end parallel do

    call MPI_FILE_GET_POSITION(fileh, head_dis, ierr)

    call MPI_FILE_READ_AT_ALL(fileh, head_dis, head_out, 1, MPI_REAL8, istat, ierr)

    deallocate(istat)

  end subroutine read_mpi_r8head

  subroutine read_mpi_i4(ftype, int_ftype, fileh, read_num, read_out)
  !*********************************************************************************************
  ! read_mpi_i4 -- Read integer mpi file
  !*********************************************************************************************
    ! -- module
    use constval_module, only: INOVAL
    ! -- inout
    integer(I4), intent(in) :: ftype, int_ftype, fileh, read_num
    integer(I4), intent(out) :: read_out(:)
    ! -- local
    integer(I4) :: i, ierr, mpi_rnum
    integer(I4), allocatable :: istat(:)
    integer(I4), allocatable :: read_val(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(istat(MPI_STATUS_SIZE))

    if (ftype == in_type(7) .or. int_ftype == 0) then
      mpi_rnum = read_num
      allocate(read_val(mpi_rnum))
    else
      mpi_rnum = read_num + 1
      allocate(read_val(mpi_rnum))
    end if

    !$omp parallel
    !$omp do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, mpi_rnum
      read_val(i) = INOVAL
    end do
    !$omp end do
    !$omp end parallel

    call MPI_FILE_READ_ALL(fileh, read_val, mpi_rnum, MPI_INTEGER, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_read(fileh)
      end if
    end if

    if (ftype == in_type(7) .or. int_ftype == 0) then
      !$omp parallel do private(i)
      do i = 1, mpi_rnum
        read_out(i) = read_val(i)
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(i)
      do i = 1, mpi_rnum
        read_out(i) = read_val(i+1)
      end do
      !$omp end parallel do
    end if

    deallocate(istat, read_val)

  end subroutine read_mpi_i4

  subroutine read_mpi_r4(ftype, int_ftype, fileh, read_num, read_out)
  !*********************************************************************************************
  ! read_mpi_r4 -- Read real4 mpi file
  !*********************************************************************************************
    ! -- module
    use constval_module, only: SNOVAL
    ! -- inout
    integer(I4), intent(in) :: ftype, int_ftype, fileh, read_num
    real(SP), intent(out) :: read_out(:)
    ! -- local
    integer(I4) :: i, ierr, mpi_rnum
    integer(I4), allocatable :: istat(:)
    real(SP), allocatable :: read_val(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(istat(MPI_STATUS_SIZE))

    if (ftype == in_type(7) .or. int_ftype == 0) then
      mpi_rnum = read_num
      allocate(read_val(mpi_rnum))
    else
      mpi_rnum = read_num + 1
      allocate(read_val(mpi_rnum))
    end if

    !$omp parallel
    !$omp do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, mpi_rnum
      read_val(i) = SNOVAL
    end do
    !$omp end do
    !$omp end parallel

    call MPI_FILE_READ_ALL(fileh, read_val, mpi_rnum, MPI_REAL4, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_read(fileh)
      end if
    end if

    if (ftype == in_type(7) .or. int_ftype == 0) then
      !$omp parallel do private(i)
      do i = 1, mpi_rnum
        read_out(i) = read_val(i)
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(i)
      do i = 1, mpi_rnum
        read_out(i) = read_val(i+1)
      end do
      !$omp end parallel do
    end if

    deallocate(istat, read_val)

  end subroutine read_mpi_r4

  subroutine read_mpi_r8(ftype, int_ftype, fileh, read_num, read_out)
  !*********************************************************************************************
  ! read_mpi_r8 -- Read real8 mpi file
  !*********************************************************************************************
    ! -- module
    use constval_module, only: DNOVAL
    ! -- inout
    integer(I4), intent(in) :: ftype, int_ftype, fileh, read_num
    real(DP), intent(out) :: read_out(:)
    ! -- local
    integer(I4) :: i, ierr, mpi_rnum
    integer(I4), allocatable :: istat(:)
    real(DP), allocatable :: read_val(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(istat(MPI_STATUS_SIZE))

    if (ftype == in_type(7) .or. int_ftype == 0) then
      mpi_rnum = read_num
      allocate(read_val(mpi_rnum))
    else
      mpi_rnum = read_num + 1
      allocate(read_val(mpi_rnum))
    end if

    !$omp parallel
    !$omp do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, mpi_rnum
      read_val(i) = DNOVAL
    end do
    !$omp end do
    !$omp end parallel

    call MPI_FILE_READ_ALL(fileh, read_val, mpi_rnum, MPI_REAL8, istat, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_read(fileh)
      end if
    end if

    if (ftype == in_type(7) .or. int_ftype == 0) then
      !$omp parallel do private(i)
      do i = 1, mpi_rnum
        read_out(i) = read_val(i)
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(i)
      do i = 1, mpi_rnum
        read_out(i) = read_val(i+1)
      end do
      !$omp end parallel do
    end if

    deallocate(istat, read_val)

  end subroutine read_mpi_r8

  subroutine skip_mpi_file(ftype, fnum, fview, err_mes, fmulti, fetime)
  !*********************************************************************************************
  ! skip_mpi_file -- Skip mpi file
  !*********************************************************************************************
    ! -- module
    use utility_module, only: get_ilen, conv_i2s
    ! -- inout
    integer(I4), intent(in) :: ftype, fnum, fview
    character(*), intent(in) :: err_mes
    real(SP), intent(in) :: fmulti
    real(SP), intent(out) :: fetime
    ! -- local
    integer(I4) :: i, ierr, time_flag, read_count
    integer(I4), allocatable :: istat(:)
    integer(KIND=MPI_OFFSET_KIND) :: head_dis, read_head
    real(SP) :: read_etime
    character(:), allocatable :: str_fnum
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; head_dis = 0
    allocate(istat(MPI_STATUS_SIZE))
    !$omp parallel do private(i)
    do i = 1, MPI_STATUS_SIZE
      istat(i) = 0
    end do
    !$omp end parallel do

    allocate(character(get_ilen(fnum)) :: str_fnum)
    call conv_i2s(fnum, str_fnum)

    if (ftype == in_type(4)) then
      call set_2dgrid_view(err_mes)
      call MPI_FILE_SET_VIEW(fnum, head_dis, MPI_REAL4, gview2d, "native", MPI_INFO_NULL, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (st_mpi%rank == 0) then
          call write_err_stop("Set view "//err_mes//" in MPI program.")
        end if
      end if
      read_count = st_grid%nx*st_grid%ny+1
    else if (ftype == in_type(6)) then
      call set_3dgrid_view(err_mes)
      call MPI_FILE_SET_VIEW(fnum, head_dis, MPI_REAL4, gview3d, "native", MPI_INFO_NULL, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (st_mpi%rank == 0) then
          call write_err_stop("Set view "//err_mes//" in MPI program.")
        end if
      end if
      read_count= st_grid%nxyz+1
    else
      read_count = 0
    end if

    time_flag = 0 ; read_head = 0
    do while (time_flag == 0)
      call MPI_FILE_READ_ALL(fnum, read_etime, 1, MPI_REAL4, istat, ierr)
      if (ierr /= MPI_SUCCESS) then
        if (st_mpi%rank == 0) then
          call write_err_read(fnum)
        end if
      end if
      fetime = read_etime*fmulti
      if (st_sim%res_type == 0) then
        time_flag = 1
        call MPI_FILE_SET_VIEW(fnum, read_head, MPI_REAL4, fview, "native", MPI_INFO_NULL,ierr)
        if (ierr /= MPI_SUCCESS) then
          if (st_mpi%rank == 0) then
            call write_err_stop("Set view "//err_mes//" in MPI program.")
          end if
        end if
      else if (fetime > st_init%rest_time .and. st_sim%res_type == 1) then
        time_flag = 1 ; read_head = head_dis*read_count*I4
        call MPI_FILE_SET_VIEW(fnum, read_head, MPI_REAL4, fview, "native", MPI_INFO_NULL, ierr)
        if (ierr /= MPI_SUCCESS) then
          if (st_mpi%rank == 0) then
            call write_err_stop("Set view "//err_mes//" in MPI program.")
          end if
        end if
      else
        call MPI_FILE_GET_POSITION(fnum, head_dis, ierr)
        if (ierr /= MPI_SUCCESS) then
          if (st_mpi%rank == 0) then
            call write_err_stop("Get position "//err_mes//" in MPI program.")
          end if
        end if
      end if
    end do

    deallocate(istat)

  end subroutine skip_mpi_file

  subroutine skip_mpi_file_int(bnum, err_mes, fmulti, fstep, finend, fnum, fetime)
  !*********************************************************************************************
  ! skip_mpi_file_int -- Skip mpi time interval list file
  !*********************************************************************************************
    ! -- module
    use constval_module, only: CHALEN
    ! -- inout
    integer(I4), intent(in) :: bnum
    character(*), intent(in) :: err_mes
    integer(I4), intent(inout) :: fnum
    real(SP), intent(in) :: fmulti, fstep, finend
    real(SP), intent(inout) :: fetime
    ! -- local
    integer(I4) :: ierr, time_flag, count_num, file_len
    character(CHALEN) :: intpath
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; time_flag = 0 ; count_num = 0

    do while (time_flag == 0)
      fetime = finend + fstep*count_num
      fetime = fetime*fmulti
      if (fetime <= st_init%rest_time .and. st_sim%res_type == 1) then
        call close_mpi_file(fnum)
        if (st_mpi%rank == 0) then
          read(unit=bnum,fmt='(a)',iostat=ierr) intpath
          if (ierr /= 0) then
            call write_err_read(bnum)
          end if
          file_len = len(intpath)
        end if
        call MPI_BCAST(intpath, file_len, MPI_CHARACTER, 0, st_mpi%comm, ierr)
        if (ierr /= MPI_SUCCESS) then
          if (st_mpi%rank == 0) then
            call write_err_stop("Broadcast file path "//err_mes//" in MPI program.")
          end if
        end if
        ! -- Open time interval mpi file (int_mpi)
          call open_int_mpi(trim(adjustl(intpath)), err_mes, fnum)
        count_num = count_num + 1
      else
        time_flag = 1
      end if
    end do

  end subroutine skip_mpi_file_int

  subroutine set_2dgrid_view(emess)
  !*********************************************************************************************
  ! set_2dgrid_view -- Set 2d grid file view
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: emess
    ! -- local
    integer(I4) :: ierr, tmptype
    integer(I4), allocatable :: xyblock(:), xytype(:)
    integer(KIND=MPI_ADDRESS_KIND) :: lb, extent
    integer(KIND=MPI_ADDRESS_KIND), allocatable :: xydis(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(xyblock(1), xydis(1), xytype(1))
    xyblock(1) = 1 ; xydis(1) = 0_MPI_ADDRESS_KIND ; xytype(1) = MPI_REAL4

    call MPI_TYPE_CREATE_STRUCT(1, xyblock, xydis, xytype, tmptype, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Create struct datatype "//emess//" in MPI program.")
      end if
    end if

    lb = 0_MPI_ADDRESS_KIND
    extent = int((st_grid%nx*st_grid%ny+1)*4, kind=MPI_ADDRESS_KIND)
    call MPI_TYPE_CREATE_RESIZED(tmptype, lb, extent, gview2d, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Resize struct datatype "//emess//" in MPI program.")
      end if
    end if

    call MPI_TYPE_COMMIT(gview2d, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Commit struct datatype "//emess//" in MPI program.")
      end if
    end if

    call MPI_TYPE_FREE(tmptype, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Free struct datatype "//emess//" in MPI program.")
      end if
    end if

    deallocate(xyblock, xydis, xytype)

  end subroutine set_2dgrid_view

  subroutine set_3dgrid_view(emess)
  !*********************************************************************************************
  ! set_3dgrid_view -- Set 3d grid file view
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: emess
    ! -- local
    integer(I4) :: ierr, tmptype
    integer(I4), allocatable :: xyzblock(:), xyztype(:)
    integer(KIND=MPI_ADDRESS_KIND) :: lb, extent
    integer(KIND=MPI_ADDRESS_KIND), allocatable :: xyzdis(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(xyzblock(1), xyzdis(1), xyztype(1))
    xyzblock(1) = 1 ; xyzdis(1) = 0_MPI_ADDRESS_KIND ; xyztype(1) = MPI_REAL4

    call MPI_TYPE_CREATE_STRUCT(1, xyzblock, xyzdis, xyztype, tmptype, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Create struct datatype "//emess//" in MPI program.")
      end if
    end if

    lb = 0_MPI_ADDRESS_KIND
    extent = int((st_grid%nxyz+1)*4, kind=MPI_ADDRESS_KIND)
    call MPI_TYPE_CREATE_RESIZED(tmptype, lb, extent, gview3d, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Resize struct datatype "//emess//" in MPI program.")
      end if
    end if

    call MPI_TYPE_COMMIT(gview3d, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Commit struct datatype "//emess//" in MPI program.")
      end if
    end if

    call MPI_TYPE_FREE(tmptype, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_stop("Free struct datatype "//emess//" in MPI program.")
      end if
    end if

    deallocate(xyzblock, xyzdis, xyztype)

  end subroutine set_3dgrid_view

  subroutine open_int_mpi(int_path, int_name, int_unit)
  !*********************************************************************************************
  ! open_int_mpi -- Open time interval mpi file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: int_path, int_name
    integer(I4), intent(out) :: int_unit
    ! -- local
    integer(I4) :: ierr, mpi_fh
    character(:), allocatable :: err_mes
    !-------------------------------------------------------------------------------------------
    call MPI_FILE_OPEN(st_mpi%comm, int_path, MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_fh, ierr)

    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        allocate(character(0) :: err_mes)
        err_mes = "Open input"//int_name//" time interval file."
        call write_err_stop(err_mes)
        deallocate(err_mes)
      end if
    end if

    int_unit = mpi_fh

  end subroutine open_int_mpi

  subroutine close_mpi_file(fileh)
  !*********************************************************************************************
  ! close_mpi_file -- Close mpi file
  !*********************************************************************************************
    ! -- module
    use utility_module, only: write_err_close
    ! -- inout
    integer(I4), intent(inout) :: fileh
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    call MPI_FILE_CLOSE(fileh, ierr)

    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        call write_err_close(fileh)
      end if
    end if

  end subroutine close_mpi_file

end module mpi_read
