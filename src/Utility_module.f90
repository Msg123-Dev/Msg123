module utility_module
  ! -- modules
  use kind_module, only: I4
  use types_module, only: mpi_set, gmap_set

  implicit none
  private
  integer(I4), public :: log_fnum = 0
  type(mpi_set), public :: st_mpi
  public :: get_file_stat, get_days
  public :: open_new_rtxt, open_new_rbin, open_new_wtxt, open_new_wbin
  public :: close_file
  public :: write_logf, write_success
  public :: write_err_read, write_err_write, write_err_close, write_err_stop
  public :: conv_unit, get_ilen, conv_i2s, conv_lower
  public :: iquick_sort, iquick_sort2
  public :: gmap_init, gmap_put, gmap_get, gmap_free

  ! -- local

  contains

  subroutine get_file_stat(file_name, file_num, is_file_opened)
  !*********************************************************************************************
  ! get_file_stat -- Get file state
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: file_name
    integer(I4), intent(out) :: file_num
    logical, intent(out)  :: is_file_opened
    ! -- local

    !-------------------------------------------------------------------------------------------
    inquire(file=file_name, opened=is_file_opened, number=file_num)

  end subroutine get_file_stat

  function get_days(nyear, nmonth) result(days)
  !*********************************************************************************************
  ! get_days -- Get the number of days
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: nyear, nmonth
    ! -- local
    integer(I4) :: days
    integer(I4) :: mon(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
    !-------------------------------------------------------------------------------------------
    days = mon(nmonth)
    if (nmonth == 2) then
      if (mod(nyear,400) == 0 .or. (mod(nyear,100) /= 0 .and. mod(nyear,4) == 0)) then
        days = 29
      end if
    end if

  end function get_days

  subroutine open_new_rtxt(stop_flag, write_flag, file_name, err_mes, file_num, file_ierr)
  !*********************************************************************************************
  ! open_new_rtxt -- Open new read text file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: stop_flag, write_flag
    character(*), intent(in) :: file_name, err_mes
    integer(I4), intent(out) :: file_num
    integer(I4), intent(out), optional :: file_ierr
    ! -- local
    integer(I4) :: ierr
    logical :: is_opened
    !-------------------------------------------------------------------------------------------
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
  !*********************************************************************************************
  ! open_new_rbin -- Open new read binary file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: stop_flag, write_flag
    character(*), intent(in) :: file_name, err_mes
    integer(I4), intent(out) :: file_num
    integer(I4), intent(out), optional :: file_ierr
    ! -- local
    integer(I4) :: ierr
    logical :: is_opened
    !-------------------------------------------------------------------------------------------
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
  !*********************************************************************************************
  ! open_new_wtxt -- Open new write text file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: file_name, err_mes
    integer(I4), intent(out) :: file_num
    ! -- local
    integer(I4) :: ierr
    logical :: is_opened
    !-------------------------------------------------------------------------------------------
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
  !*********************************************************************************************
  ! open_new_wbin -- Open new write binary file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: file_name, err_mes
    integer(I4), intent(out) :: file_num
    ! -- local
    integer(I4) :: ierr
    logical :: is_opened
    !-------------------------------------------------------------------------------------------
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
  !*********************************************************************************************
  ! close_file -- Close file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    integer(I4) :: ierr
    logical :: open_check
    !-------------------------------------------------------------------------------------------
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
  !*********************************************************************************************
  ! write_logf -- Write log file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: err_mes
    ! -- local

    !-------------------------------------------------------------------------------------------
    write(log_fnum,'(a)') err_mes

  end subroutine write_logf

  subroutine write_success(suc_mes, file_num)
  !*********************************************************************************************
  ! write_success -- Write success in log file
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(in) :: suc_mes
    integer(I4), intent(in) :: file_num
    ! -- local
    character(:), allocatable :: file_str
    !-------------------------------------------------------------------------------------------
    allocate(character(get_ilen(file_num)) :: file_str)
    call conv_i2s(file_num, file_str)
    write(log_fnum,'(a)') "Succeed!! "//suc_mes//" as file number "//file_str//"."

    deallocate(file_str)

  end subroutine write_success

  subroutine write_err_read(file_num)
  !*********************************************************************************************
  ! write_err_read -- Write error for read
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    character(:), allocatable :: file_str
    !-------------------------------------------------------------------------------------------
    allocate(character(get_ilen(file_num)) :: file_str)
    call conv_i2s(file_num, file_str)
    call write_err_stop("Read file number "//file_str)

    deallocate(file_str)

  end subroutine write_err_read

  subroutine write_err_write(file_num)
  !*********************************************************************************************
  ! write_err_write -- Write error for write
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    character(:), allocatable :: file_str
    !-------------------------------------------------------------------------------------------
    allocate(character(get_ilen(file_num)) :: file_str)
    call conv_i2s(file_num, file_str)
    call write_err_stop("Write file number "//file_str)

    deallocate(file_str)

  end subroutine write_err_write

  subroutine write_err_close(file_num)
  !*********************************************************************************************
  ! write_err_close -- Write error for close
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: file_num
    ! -- local
    character(:), allocatable :: file_str
    !-------------------------------------------------------------------------------------------
    allocate(character(get_ilen(file_num)) :: file_str)
    call conv_i2s(file_num, file_str)
    call write_err_stop("Close file number "//file_str)

    deallocate(file_str)

  end subroutine write_err_close

  subroutine write_err_stop(err_mes)
  !*********************************************************************************************
  ! write_err_stop -- Write error log file and stop
  !*********************************************************************************************
    ! -- module
#ifdef MPI_MSG
    use mpi_initfin, only: abort_proc
#endif
    ! -- inout
    character(*), intent(in) :: err_mes
    ! -- local

    !-------------------------------------------------------------------------------------------
    write(log_fnum,'(a)') "Error!! "//err_mes

#ifdef MPI_MSG
    ! -- Abort process (proc)
      call abort_proc(st_mpi%rank, log_fnum)
#endif
    error stop 1

  end subroutine write_err_stop

  subroutine conv_unit(rank, time_char, mess_char, date, time_conv)
  !*********************************************************************************************
  ! conv_unit -- Convert unit
  !*********************************************************************************************
    ! -- modules
    use kind_module, only: SP
    use constval_module, only: MINSEC, HOURSEC, DAYSEC, SONE
    ! -- inout
    integer(I4), intent(in) :: rank
    character(*), intent(in) :: time_char, mess_char
    integer(I4), intent(in) :: date(:)
    real(SP), intent(out) :: time_conv
    ! -- local
    integer(I4) :: mday, year
    !-------------------------------------------------------------------------------------------
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

  function get_ilen(num) result(num_digit)
  !*********************************************************************************************
  ! get_ilen -- Get the length of an integer
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: num
    ! -- local
    integer(I4) :: num_digit
    !-------------------------------------------------------------------------------------------
    if (num /= 0) then
      if (num < 0) then
        num_digit = int(log10(dble(abs(num)))) + 2
      else
        num_digit = int(log10(dble(num))) + 1
      end if
    else
      num_digit = 1
    end if

  end function get_ilen

  subroutine conv_i2s(num, str)
  !*********************************************************************************************
  ! conv_i2s -- Convert integer to string
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: num
    character(*), intent(out) :: str
    ! -- local
    !-------------------------------------------------------------------------------------------
    write(str,'(I0)') num

  end subroutine conv_i2s

  function conv_lower(in_char) result(low_char)
  !*********************************************************************************************
  ! conv_lower -- Convert a character string to lower case
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    character(*), intent(in) :: in_char
    ! -- local
    integer(I4) :: i, char_num
    character(len(in_char)) :: low_char
    !-------------------------------------------------------------------------------------------
    low_char(:) = in_char(:)
    do i = 1, len(in_char)
      char_num = iachar(low_char(i:i))
      if (char_num >= iachar("A") .and. char_num <= iachar("Z")) then
        low_char(i:i) = achar(char_num + 32)
      end if
    end do

  end function conv_lower

  recursive subroutine iquick_sort(in_x, first, last)
  !*********************************************************************************************
  ! iquick_sort -- integer quick sort
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(inout) :: in_x(:)
    integer(I4), intent(in) :: first, last
    ! -- local
    integer(I4) :: i, j, k, t
    !-------------------------------------------------------------------------------------------
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

  recursive subroutine iquick_sort2(key_x, val_x, first, last)
  !*********************************************************************************************
  ! iquick_sort2 -- integer quick sort of key_x, carrying val_x along with the same swaps
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(inout) :: key_x(:), val_x(:)
    integer(I4), intent(in) :: first, last
    ! -- local
    integer(I4) :: low_pos, high_pos, pivot, swap_tmp
    !-------------------------------------------------------------------------------------------
    pivot = key_x((first+last)/2) ; low_pos = first ; high_pos = last
    do
      do while(key_x(low_pos) < pivot)
        low_pos = low_pos + 1
      end do
      do while(pivot < key_x(high_pos))
        high_pos = high_pos - 1
      end do
      if (low_pos >= high_pos) then
        exit
      end if
      swap_tmp = key_x(low_pos) ; key_x(low_pos) = key_x(high_pos) ; key_x(high_pos) = swap_tmp
      swap_tmp = val_x(low_pos) ; val_x(low_pos) = val_x(high_pos) ; val_x(high_pos) = swap_tmp
      low_pos = low_pos + 1 ; high_pos = high_pos - 1
    end do
    if (first < low_pos-1) then
      call iquick_sort2(key_x, val_x, first, low_pos-1)
    end if
    if (high_pos+1 < last) then
      call iquick_sort2(key_x, val_x, high_pos+1, last)
    end if

  end subroutine iquick_sort2

  subroutine gmap_init(g_map, entry_num)
  !*********************************************************************************************
  ! gmap_init -- global to local map initialize
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    type(gmap_set), intent(inout) :: g_map
    integer(I4), intent(in) :: entry_num
    ! -- local
    integer(I4) :: slot
    !-------------------------------------------------------------------------------------------
    g_map%table_num = entry_num*2 + (entry_num/2) + 3

    if (mod(entry_num, 2) == 0) then
      g_map%table_num = g_map%table_num + 1
    end if

    if (allocated(g_map%table_key)) then
      deallocate(g_map%table_key, g_map%table_val)
    end if
    allocate(g_map%table_key(g_map%table_num), g_map%table_val(g_map%table_num))

    !$omp parallel do private(slot)
    do slot = 1, g_map%table_num
      g_map%table_key(slot) = 0 ; g_map%table_val(slot) = 0
    end do
    !$omp end parallel do

  end subroutine gmap_init

  subroutine gmap_put(g_map, map_key, map_val)
  !*********************************************************************************************
  ! gmap_put -- global to local map put
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    type(gmap_set), intent(inout) :: g_map
    integer(I4), intent(in) :: map_key, map_val
    ! -- local
    integer(I4) :: slot
    !-------------------------------------------------------------------------------------------
    slot = mod(map_key - 1, g_map%table_num) + 1

    do while (g_map%table_key(slot) /= 0 .and. g_map%table_key(slot) /= map_key)
      slot = mod(slot, g_map%table_num) + 1
    end do

    g_map%table_key(slot) = map_key ; g_map%table_val(slot) = map_val

  end subroutine gmap_put

  function gmap_get(g_map, map_key) result(map_val)
  !*********************************************************************************************
  ! gmap_get -- global to local map get
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    type(gmap_set), intent(in) :: g_map
    integer(I4), intent(in) :: map_key
    ! -- local
    integer(I4) :: map_val, slot
    !-------------------------------------------------------------------------------------------
    slot = mod(map_key - 1, g_map%table_num) + 1

    do while (g_map%table_key(slot) /= 0 .and. g_map%table_key(slot) /= map_key)
      slot = mod(slot, g_map%table_num) + 1
    end do

    if (g_map%table_key(slot) == map_key) then
      map_val = g_map%table_val(slot)
    else
      map_val = 0
    end if

  end function gmap_get

  subroutine gmap_free(g_map)
  !*********************************************************************************************
  ! gmap_free -- global to local map free
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    type(gmap_set), intent(inout) :: g_map
    ! -- local

    !-------------------------------------------------------------------------------------------
    if (allocated(g_map%table_key)) then
      deallocate(g_map%table_key, g_map%table_val)
    end if
    g_map%table_num = 0

  end subroutine gmap_free

end module utility_module
