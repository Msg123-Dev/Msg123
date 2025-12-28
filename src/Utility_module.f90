module utility_module
  ! -- modules
  use kind_module, only: I4

  implicit none
  private
  integer(I4), public :: log_fnum = 0, pro_totn = 1, my_rank = 0
  public :: get_file_stat, get_days
  public :: open_new_rtxt, open_new_rbin, open_new_wtxt, open_new_wbin
  public :: close_file
  public :: write_logf, write_success
  public :: write_err_read, write_err_write, write_err_close, write_err_stop
  public :: conv_unit, conv_i2s
  public :: iquick_sort

  ! -- local

  contains

  subroutine get_file_stat(file_name, file_num, is_file_opened)
  !***************************************************************************************
  ! get_file_stat -- Get file state
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: file_name
    integer(I4), intent(out) :: file_num
    logical, intent(out)  :: is_file_opened
    ! -- local

    !-------------------------------------------------------------------------------------
    inquire(file=file_name, opened=is_file_opened, number=file_num)

  end subroutine get_file_stat

  function get_days(nyear, nmonth) result(days)
  !***************************************************************************************
  ! get_days -- Get the number of days
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: nyear, nmonth
    ! -- local
    integer(I4) :: days
    integer(I4) :: mon(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
    !-------------------------------------------------------------------------------------
    days = mon(nmonth)
    if (nmonth == 2) then
      if (mod(nyear,400) == 0 .or. (mod(nyear,100) /= 0 .and. mod(nyear,4) == 0)) then
        days = 29
      end if
    end if

  end function get_days

  subroutine open_new_rtxt(stop_flag, write_flag, file_name, err_mes, file_num, file_ierr)
  !***************************************************************************************
  ! open_new_rtxt -- Open new read text file
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: stop_flag, write_flag
    character(*), intent(in) :: file_name, err_mes
    integer(I4), intent(out) :: file_num
    integer(I4), intent(out), optional :: file_ierr
    ! -- local
    integer(I4) :: ierr
    logical :: is_opened
    !-------------------------------------------------------------------------------------
    ! -- Get file state (file_stat)
      call get_file_stat(file_name, file_num, is_opened)

    ierr = 0
    open(newunit=file_num, file=file_name, form='formatted',&
         access='sequential', status='old', action='read', iostat=ierr)

    if (ierr == 0 .and. write_flag == 1) then
      call write_success("Open "//err_mes//" file", file_num)
    else if (stop_flag == 1) then
      call write_err_stop("Open "//err_mes//" file.")
    end if

    if (present(file_ierr)) then
      file_ierr = ierr
    end if

  end subroutine open_new_rtxt

  subroutine open_new_rbin(stop_flag, write_flag, file_name, err_mes, file_num, file_ierr)
  !***************************************************************************************
  ! open_new_rbin -- Open new read binary file
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: stop_flag, write_flag
    character(*), intent(in) :: file_name, err_mes
    integer(I4), intent(out) :: file_num
    integer(I4), intent(out), optional :: file_ierr
    ! -- local
    integer(I4) :: ierr
    logical :: is_opened
    !-------------------------------------------------------------------------------------
    ! -- Get file state (file_stat)
      call get_file_stat(file_name, file_num, is_opened)

    ierr = 0
    open(newunit=file_num, file=file_name, form='unformatted',&
         access='stream', status='old', action='read', iostat=ierr)

    if (ierr == 0 .and. write_flag == 1) then
      call write_success("Open "//err_mes//" file", file_num)
    else if (stop_flag == 1) then
      call write_err_stop("Open "//err_mes//" file.")
    end if

    if (present(file_ierr)) then
      file_ierr = ierr
    end if

  end subroutine open_new_rbin

  subroutine open_new_wtxt(file_name, err_mes, file_num)
  !***************************************************************************************
  ! open_new_wtxt -- Open new write text file
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: file_name, err_mes
    integer(I4), intent(out) :: file_num
    ! -- local
    integer(I4) :: ierr
    logical :: is_opened
    !-------------------------------------------------------------------------------------
    ! -- Get file state (file_stat)
      call get_file_stat(file_name, file_num, is_opened)

    ierr = 0
    if (.not. is_opened) then
      open(newunit=file_num, file=file_name, form='formatted',&
           access='sequential', status='replace', action='write', iostat=ierr)
    end if

    if (ierr == 0) then
      call write_success("Open "//err_mes//" file", file_num)
    else
      call write_err_stop("Open "//err_mes//" file.")
    end if

  end subroutine open_new_wtxt

  subroutine open_new_wbin(file_name, err_mes, file_num)
  !***************************************************************************************
  ! open_new_wbin -- Open new write binary file
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: file_name, err_mes
    integer(I4), intent(out) :: file_num
    ! -- local
    integer(I4) :: ierr
    logical :: is_opened
    !-------------------------------------------------------------------------------------
    ! -- Get file state (file_stat)
      call get_file_stat(file_name, file_num, is_opened)

    ierr = 0
    if (.not. is_opened) then
      open(newunit=file_num, file=file_name, form='unformatted',&
           access='stream', status='replace', action='write', iostat=ierr)
    end if

    if (ierr == 0) then
      call write_success("Open "//err_mes//" file", file_num)
    else
      call write_err_stop("Open "//err_mes//" file.")
    end if

  end subroutine open_new_wbin

  subroutine close_file(file_num)
  !***************************************************************************************
  ! close_file -- Close file
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    integer(I4) :: ierr
    logical :: open_check
    !-------------------------------------------------------------------------------------
    inquire(unit=file_num, opened=open_check)

    ierr = 0
    if (open_check) then
      close(unit=file_num, iostat=ierr)
    end if

    if (ierr /= 0) then
      call write_err_close(file_num)
    end if

  end subroutine close_file

  subroutine write_logf(err_mes)
  !***************************************************************************************
  ! write_logf -- Write log file
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: err_mes
    ! -- local

    !-------------------------------------------------------------------------------------
    write(log_fnum,'(a)') err_mes

  end subroutine write_logf

  subroutine write_success(suc_mes, file_num)
  !***************************************************************************************
  ! write_success -- Write success in log file
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: suc_mes
    integer(I4), intent(in) :: file_num
    ! -- local
    character(:), allocatable :: file_str
    !-------------------------------------------------------------------------------------
    file_str = conv_i2s(file_num)
    write(log_fnum,'(a)') "Succeed!! "//suc_mes//" as file number "//file_str//"."

    deallocate(file_str)

  end subroutine write_success

  subroutine write_err_read(file_num)
  !***************************************************************************************
  ! write_err_read -- Write error for read
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    character(:), allocatable :: file_str
    !-------------------------------------------------------------------------------------
    file_str = conv_i2s(file_num)
    call write_err_stop("Read file number "//file_str)

    deallocate(file_str)

  end subroutine write_err_read

  subroutine write_err_write(file_num)
  !***************************************************************************************
  ! write_err_write -- Write error for write
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    character(:), allocatable :: file_str
    !-------------------------------------------------------------------------------------
    file_str = conv_i2s(file_num)
    call write_err_stop("Write file number "//file_str)

    deallocate(file_str)

  end subroutine write_err_write

  subroutine write_err_close(file_num)
  !***************************************************************************************
  ! write_err_close -- Write error for close
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    character(:), allocatable :: file_str
    !-------------------------------------------------------------------------------------
    file_str = conv_i2s(file_num)
    call write_err_stop("Close file number "//file_str)

    deallocate(file_str)

  end subroutine write_err_close

  subroutine write_err_stop(err_mes)
  !***************************************************************************************
  ! write_err_stop -- Write error log file and stop
  !***************************************************************************************
    ! -- module
#ifdef MPI_MSG
    use mpi_initfin, only: abort_proc
#endif
    ! -- inout
    character(*), intent(in) :: err_mes
    ! -- local

    !-------------------------------------------------------------------------------------
    write(log_fnum,'(a)') "Error!! "//err_mes

#ifdef MPI_MSG
    ! -- Abort process (proc)
      call abort_proc(my_rank, log_fnum)
#endif
    stop

  end subroutine write_err_stop

  subroutine conv_unit(rank, time_char, mess_char, date, time_conv)
  !***************************************************************************************
  ! conv_unit -- Convert unit
  !***************************************************************************************
    ! -- modules
    use kind_module, only: SP
    use constval_module, only: SONE, MINSEC, HOURSEC, DAYSEC
    ! -- inout
    integer(I4), intent(in) :: rank
    character(*), intent(in) :: time_char, mess_char
    integer(I4), intent(in) :: date(:)
    real(SP), intent(out) :: time_conv
    ! -- local
    integer(I4) :: mday, year
    !-------------------------------------------------------------------------------------
    if (time_char == "SEC") then
        time_conv = SONE
    else if (time_char == "MIN") then
      time_conv = MINSEC
    else if (time_char == "HOU") then
      time_conv = HOURSEC
    else if (time_char == "DAY") then
      time_conv = DAYSEC
!    else if (time_char == "MON") then
!      if (date(1) == 0) then
!        mday = 30
!      else if (date(2) == 4 .or. date(2) == 6 .or. date(2) == 9 .or. date(2) == 11) then
!        mday = 30
!      else if (date(2) /= 2) then
!        mday = 31
!      else
!        mday = 28
!        if (mod(date(1),400) == 0 .or. (mod(date(1),100) /= 0 .and. mod(date(1),4) == 0)) then
!          mday = 29
!        end if
!      end if
!      time_conv = DAYSEC*mday
    else if (time_char == "YEA") then
      if (date(2) <= 2) then
        year = date(1)
      else
        year = date(1) + 1
      end if
      mday = 365
      if (year == 0) then
        mday = 365
      else if (mod(year,400) == 0 .or. (mod(year,100) /= 0 .and. mod(year,4) == 0)) then
        mday = 366
      end if
      time_conv = DAYSEC*mday
    else if (rank == 0) then
      call write_err_stop("Specified wrong time unit in "//mess_char//".")
    end if

  end subroutine conv_unit

  function conv_i2s(num) result(str)
  !***************************************************************************************
  ! conv_i2s -- Convert integer to string
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: num
    ! -- local
    integer(I4) :: num_digit
    character(:), allocatable :: str
    !-------------------------------------------------------------------------------------
    if (num /= 0) then
      if (num < 0) then
        num_digit = int(log10(dble(abs(num)))) + 2
      else
        num_digit = int(log10(dble(num))) + 1
      end if
    else
      num_digit = 1
    end if

    allocate(character(num_digit)::str)
    write(str,'(I0)') num

  end function conv_i2s

  recursive subroutine iquick_sort(in_x, first, last)
  !***************************************************************************************
  ! iquick_sort -- integer quick sort
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(inout) :: in_x(:)
    integer(I4), intent(in) :: first, last
    ! -- local
    integer(I4) :: i, j, k, t
    !-------------------------------------------------------------------------------------
    k = in_x((first+last)/2) ; i = first ; j = last
    do
      do while(in_x(i) < k)
        i = i + 1
      end do
      do while(k < in_x(j))
        j = j - 1
      end do
      if (i >= j) then
        exit
      end if
      t = in_x(i) ; in_x(i) = in_x(j) ; in_x(j) = t
      i = i + 1 ; j=j-1
    end do
    if (first < i-1) then
      call iquick_sort(in_x, first, i-1)
    end if
    if (j+1 < last) then
      call iquick_sort(in_x, j+1, last)
    end if

  end subroutine iquick_sort

end module utility_module
