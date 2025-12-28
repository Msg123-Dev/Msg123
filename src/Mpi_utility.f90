module mpi_utility
  ! -- modules
  use kind_module, only: I4, SP, DP
  use mpi_initfin, only: my_comm, abort_proc
  use utility_module, only: log_fnum, my_rank
  use mpi

  implicit none
  private
  public :: barrier_proc
  public :: mpisum_val, mpimax_val
  public :: bcast_val, bcast_char, bcast_file, bcast_extr_set
  public :: gather_val

  interface mpisum_val
    module procedure mpisum_i4_scalar
    module procedure mpisum_r4_scalar
    module procedure mpisum_r8_scalar
    module procedure mpisum_i4_array
    module procedure mpisum_r4_array
    module procedure mpisum_r8_array
  end interface

  interface mpimax_val
    module procedure mpimax_i4_scalar
    module procedure mpimax_r4_scalar
    module procedure mpimax_r8_scalar
    module procedure mpimax_i4_array
    module procedure mpimax_r4_array
    module procedure mpimax_r8_array
  end interface

  interface bcast_val
    module procedure bcast_i4_scalar
    module procedure bcast_r4_scalar
    module procedure bcast_r8_scalar
  end interface

  interface bcast_file
    module procedure bcast_file_path
    module procedure bcast_ftype_path
    module procedure bcast_path_unit
    module procedure bcast_ftype_path_unit
  end interface

  interface gather_val
    module procedure gather_i4_array
    module procedure gather_r4_array
    module procedure gather_r8_array
  end interface

  ! -- local

  contains

  subroutine barrier_proc()
  !***************************************************************************************
  ! barrier_proc -- Barrier process
  !***************************************************************************************
    ! -- module

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BARRIER(my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Barrier in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine barrier_proc

  subroutine mpisum_i4_scalar(loc_num, err_mes, sum_num)
  !***************************************************************************************
  ! mpisum_i4_scalar -- Sum integer value for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: sum_num
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_num, sum_num, 1, MPI_INTEGER, MPI_SUM, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpisum_i4_scalar

  subroutine mpisum_r4_scalar(loc_val, err_mes, sum_val)
  !***************************************************************************************
  ! mpisum_r4_scalar -- Sum real4 value for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: sum_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, sum_val, 1, MPI_REAL4, MPI_SUM, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpisum_r4_scalar

  subroutine mpisum_r8_scalar(loc_val, err_mes, sum_val)
  !***************************************************************************************
  ! mpisum_r8_scalar -- Sum real8 value for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: sum_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, sum_val, 1, MPI_REAL8, MPI_SUM, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpisum_r8_scalar

  subroutine mpisum_i4_array(loc_array, err_mes, sum_array)
  !***************************************************************************************
  ! mpisum_i4_array -- Sum integer array for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: sum_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, sum_array, a_len, MPI_INTEGER, MPI_SUM, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpisum_i4_array

  subroutine mpisum_r4_array(loc_array, err_mes, sum_array)
  !***************************************************************************************
  ! mpisum_r4_array -- Sum real4 array for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: sum_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, sum_array, a_len, MPI_REAL4, MPI_SUM, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpisum_r4_array

  subroutine mpisum_r8_array(loc_array, err_mes, sum_array)
  !***************************************************************************************
  ! mpisum_r8_array -- Sum real8 array for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: sum_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, sum_array, a_len, MPI_REAL8, MPI_SUM, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpisum_r8_array

  subroutine mpimax_i4_scalar(loc_num, err_mes, max_num)
  !***************************************************************************************
  ! mpimax_i4_scalar -- Max integer value for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: max_num
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_num, max_num, 1, MPI_INTEGER, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpimax_i4_scalar

  subroutine mpimax_r4_scalar(loc_val, err_mes, max_val)
  !***************************************************************************************
  ! mpimax_r4_scalar -- Max real4 value for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: max_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, max_val, 1, MPI_REAL4, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpimax_r4_scalar

  subroutine mpimax_r8_scalar(loc_val, err_mes, max_val)
  !***************************************************************************************
  ! mpimax_r8 -- Max real value for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: max_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, max_val, 1, MPI_REAL8, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpimax_r8_scalar

  subroutine mpimax_i4_array(loc_array, err_mes, max_array)
  !***************************************************************************************
  ! mpimax_i4_array -- Max integer array for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: max_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, max_array, a_len, MPI_INTEGER, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpimax_i4_array

  subroutine mpimax_r4_array(loc_array, err_mes, max_array)
  !***************************************************************************************
  ! mpimax_r4_array -- Max real4 array for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: max_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, max_array, a_len, MPI_REAL4, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpimax_r4_array

  subroutine mpimax_r8_array(loc_array, err_mes, max_array)
  !***************************************************************************************
  ! mpimax_r8_array -- Max real8 array for MPI
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: max_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, max_array, a_len, MPI_REAL8, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine mpimax_r8_array

  subroutine bcast_i4_scalar(iscalar, err_mes)
  !***************************************************************************************
  ! bcast_i4_scalar -- Bcast integer scalar value
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: iscalar
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(iscalar, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_i4_scalar

  subroutine bcast_r4_scalar(rscalar, err_mes)
  !***************************************************************************************
  ! bcast_r4_scalar -- Bcast real4 scalar value
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(inout) :: rscalar
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(rscalar, 1, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_r4_scalar

  subroutine bcast_r8_scalar(rscalar, err_mes)
  !***************************************************************************************
  ! bcast_r8_scalar -- Bcast real8 scalar value
  !***************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(inout) :: rscalar
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(rscalar, 1, MPI_REAL8, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_r8_scalar

  subroutine bcast_char(err_mes)
  !***************************************************************************************
  ! bcast_char -- Bcast character
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(inout) :: err_mes
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: rank_num, send_num
    integer(I4) :: str_len
    !-------------------------------------------------------------------------------------
    rank_num = 0 ; send_num = 0
    if (len_trim(adjustl(err_mes)) /= 0) then
      rank_num = my_rank
    end if

    ierr = 0
    call MPI_ALLREDUCE(rank_num, send_num, 1, MPI_INTEGER, MPI_MAX, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce rank number in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    str_len = len_trim(err_mes)
    call MPI_BCAST(err_mes, str_len, MPI_CHARACTER, send_num, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast message in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_char

  subroutine bcast_file_path(file_path, err_mes)
  !***************************************************************************************
  ! bcast_file_path -- Bcast file path
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(inout) :: file_path
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(file_path)
    call MPI_BCAST(file_path, str_len, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_file_path

  subroutine bcast_ftype_path(file_type, file_path, err_mes)
  !***************************************************************************************
  ! bcast_ftype_path -- Bcast file type and path
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: file_type
    character(*), intent(inout) :: file_path
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(file_path)
    call MPI_BCAST(file_type, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file type in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    call MPI_BCAST(file_path, str_len, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_ftype_path

  subroutine bcast_path_unit(file_path, file_unit, err_mes)
  !***************************************************************************************
  ! bcast_path_unit -- Bcast file path and convert unit
  !***************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(inout) :: file_path, file_unit
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len, uni_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(file_path) ; uni_len = len_trim(file_unit)
    call MPI_BCAST(file_path, str_len, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    call MPI_BCAST(file_unit, uni_len, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file unit in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_path_unit

  subroutine bcast_ftype_path_unit(file_type, file_path, file_unit, err_mes)
  !***************************************************************************************
  ! bcast_ftype_path_unit -- Bcast file type, path and convert unit
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: file_type
    character(*), intent(inout) :: file_path, file_unit
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len, uni_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(file_path) ; uni_len = len_trim(file_unit)

    call MPI_BCAST(file_type, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file type in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    call MPI_BCAST(file_path, str_len, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    call MPI_BCAST(file_unit, uni_len, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file unit in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_ftype_path_unit

  subroutine bcast_extr_set(extr_type, extr_step, extr_end, extr_path, err_mes)
  !***************************************************************************************
  ! bcast_extr_set -- Bcast extra setting
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: extr_type
    real(SP), intent(inout) :: extr_step, extr_end
    character(*), intent(inout) :: extr_path
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(extr_path)

    call MPI_BCAST(extr_type, 1, MPI_INTEGER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file type in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    call MPI_BCAST(extr_step, 1, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file step in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    call MPI_BCAST(extr_end, 1, MPI_REAL4, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file end time in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    call MPI_BCAST(extr_path, str_len, MPI_CHARACTER, 0, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

  end subroutine bcast_extr_set

  subroutine gather_i4_array(num_prot, in_num, loc_array, glo_array, err_mes)
  !***************************************************************************************
  ! gather_i4_array -- Gather integer array
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_prot, in_num
    integer(I4), intent(in) :: loc_array(:)
    integer(I4), intent(out) :: glo_array(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: i, sum_num, ierr
    integer(I4), allocatable :: rec_num(:), rec_count(:), rec_dis(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(rec_num(num_prot))
    !$omp parallel workshare
    rec_num(:) = 0
    !$omp end parallel workshare
    call MPI_ALLGATHER(in_num, 1, MPI_INTEGER, rec_num, 1, MPI_INTEGER, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    allocate(rec_count(0:num_prot-1), rec_dis(0:num_prot-1))
    !$omp parallel workshare
    rec_count(:) = 0 ; rec_dis(:) = 0
    !$omp end parallel workshare
    sum_num = 0
    do i = 1, num_prot
      rec_dis(i-1) = sum_num
      sum_num = sum_num + rec_num(i)
      rec_count(i-1) = rec_num(i)
    end do

    call MPI_ALLGATHERV(loc_array, in_num, MPI_INTEGER, glo_array, rec_count, rec_dis,&
                        MPI_INTEGER, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    deallocate(rec_num, rec_count, rec_dis)

  end subroutine gather_i4_array

  subroutine gather_r4_array(num_prot, in_num, loc_array, glo_array, err_mes)
  !***************************************************************************************
  ! gather_r4_array -- Gather real4 array
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_prot, in_num
    real(SP), intent(in) :: loc_array(:)
    real(SP), intent(out) :: glo_array(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: i, sum_num, ierr
    integer(I4), allocatable :: rec_num(:), rec_count(:), rec_dis(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(rec_num(num_prot))
    !$omp parallel workshare
    rec_num(:) = 0
    !$omp end parallel workshare
    call MPI_ALLGATHER(in_num, 1, MPI_INTEGER, rec_num, 1, MPI_INTEGER, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    allocate(rec_count(0:num_prot-1), rec_dis(0:num_prot-1))
    !$omp parallel workshare
    rec_count(:) = 0 ; rec_dis(:) = 0
    !$omp end parallel workshare
    sum_num = 0
    do i = 1, num_prot
      rec_dis(i-1) = sum_num
      sum_num = sum_num + rec_num(i)
      rec_count(i-1) = rec_num(i)
    end do

    call MPI_ALLGATHERV(loc_array, in_num, MPI_REAL8, glo_array, rec_count, rec_dis,&
                        MPI_REAL8, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    deallocate(rec_num, rec_count, rec_dis)

  end subroutine gather_r4_array

  subroutine gather_r8_array(num_prot, in_num, loc_array, glo_array, err_mes)
  !***************************************************************************************
  ! gather_r8_array -- Gather real8 array
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_prot, in_num
    real(DP), intent(in) :: loc_array(:)
    real(DP), intent(out) :: glo_array(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: i, sum_num, ierr
    integer(I4), allocatable :: rec_num(:), rec_count(:), rec_dis(:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(rec_num(num_prot))
    !$omp parallel workshare
    rec_num(:) = 0
    !$omp end parallel workshare
    call MPI_ALLGATHER(in_num, 1, MPI_INTEGER, rec_num, 1, MPI_INTEGER, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    allocate(rec_count(0:num_prot-1), rec_dis(0:num_prot-1))
    !$omp parallel workshare
    rec_count(:) = 0 ; rec_dis(:) = 0
    !$omp end parallel workshare
    sum_num = 0
    do i = 1, num_prot
      rec_dis(i-1) = sum_num
      sum_num = sum_num + rec_num(i)
      rec_count(i-1) = rec_num(i)
    end do

    call MPI_ALLGATHERV(loc_array, in_num, MPI_REAL8, glo_array, rec_count, rec_dis,&
                        MPI_REAL8, my_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (my_rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" array in MPI program."
      end if
      call abort_proc(my_rank, log_fnum)
    end if

    deallocate(rec_num, rec_count, rec_dis)

  end subroutine gather_r8_array

end module mpi_utility
