module mpi_utility
  ! -- modules
  use kind_module, only: I4, SP, DP
  use mpi_initfin, only: abort_proc
  use utility_module, only: log_fnum, st_mpi
  use mpi

  implicit none
  private
  public :: barrier_proc
  public :: mpisum_val, mpimax_val, mpimin_val, mpiexscan_val
  public :: bcast_val, bcast_char, bcast_file, bcast_extr_set
  public :: gather_val, scatterv_val
  public :: alltoall_val, alltoallv_val

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

  interface mpimin_val
    module procedure mpimin_i4_scalar
    module procedure mpimin_r4_scalar
    module procedure mpimin_r8_scalar
    module procedure mpimin_i4_array
    module procedure mpimin_r4_array
    module procedure mpimin_r8_array
  end interface

  interface mpiexscan_val
    module procedure mpiexscan_i4_array
    module procedure mpiexscan_r4_array
    module procedure mpiexscan_r8_array
  end interface

  interface bcast_val
    module procedure bcast_i4_scalar
    module procedure bcast_i4_array
    module procedure bcast_r4_scalar
    module procedure bcast_r4_array
    module procedure bcast_r8_scalar
    module procedure bcast_r8_array
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

  interface scatterv_val
    module procedure scatterv_i4_array
    module procedure scatterv_r4_array
  end interface

  interface alltoall_val
    module procedure alltoall_i4_array
    module procedure alltoall_r4_array
    module procedure alltoall_r8_array
  end interface

  interface alltoallv_val
    module procedure alltoallv_i4_array
    module procedure alltoallv_r4_array
    module procedure alltoallv_r8_array
  end interface

  ! -- local

  contains

  subroutine barrier_proc()
  !*********************************************************************************************
  ! barrier_proc -- Barrier process
  !*********************************************************************************************
    ! -- module

    ! -- inout

    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BARRIER(st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Barrier in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine barrier_proc

  subroutine mpisum_i4_scalar(loc_num, err_mes, sum_num)
  !*********************************************************************************************
  ! mpisum_i4_scalar -- Sum integer value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: sum_num
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_num, sum_num, 1, MPI_INTEGER, MPI_SUM, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpisum_i4_scalar

  subroutine mpisum_r4_scalar(loc_val, err_mes, sum_val)
  !*********************************************************************************************
  ! mpisum_r4_scalar -- Sum real4 value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: sum_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, sum_val, 1, MPI_REAL4, MPI_SUM, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpisum_r4_scalar

  subroutine mpisum_r8_scalar(loc_val, err_mes, sum_val)
  !*********************************************************************************************
  ! mpisum_r8_scalar -- Sum real8 value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: sum_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, sum_val, 1, MPI_REAL8, MPI_SUM, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpisum_r8_scalar

  subroutine mpisum_i4_array(loc_array, err_mes, sum_array)
  !*********************************************************************************************
  ! mpisum_i4_array -- Sum integer array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: sum_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, sum_array, a_len, MPI_INTEGER, MPI_SUM, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpisum_i4_array

  subroutine mpisum_r4_array(loc_array, err_mes, sum_array)
  !*********************************************************************************************
  ! mpisum_r4_array -- Sum real4 array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: sum_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, sum_array, a_len, MPI_REAL4, MPI_SUM, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpisum_r4_array

  subroutine mpisum_r8_array(loc_array, err_mes, sum_array)
  !*********************************************************************************************
  ! mpisum_r8_array -- Sum real8 array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: sum_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, sum_array, a_len, MPI_REAL8, MPI_SUM, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpisum_r8_array

  subroutine mpimax_i4_scalar(loc_num, err_mes, max_num)
  !*********************************************************************************************
  ! mpimax_i4_scalar -- Max integer value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: max_num
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_num, max_num, 1, MPI_INTEGER, MPI_MAX, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimax_i4_scalar

  subroutine mpimax_r4_scalar(loc_val, err_mes, max_val)
  !*********************************************************************************************
  ! mpimax_r4_scalar -- Max real4 value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: max_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, max_val, 1, MPI_REAL4, MPI_MAX, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimax_r4_scalar

  subroutine mpimax_r8_scalar(loc_val, err_mes, max_val)
  !*********************************************************************************************
  ! mpimax_r8 -- Max real value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: max_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, max_val, 1, MPI_REAL8, MPI_MAX, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimax_r8_scalar

  subroutine mpimax_i4_array(loc_array, err_mes, max_array)
  !*********************************************************************************************
  ! mpimax_i4_array -- Max integer array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: max_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, max_array, a_len, MPI_INTEGER, MPI_MAX, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimax_i4_array

  subroutine mpimax_r4_array(loc_array, err_mes, max_array)
  !*********************************************************************************************
  ! mpimax_r4_array -- Max real4 array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: max_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, max_array, a_len, MPI_REAL4, MPI_MAX, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimax_r4_array

  subroutine mpimax_r8_array(loc_array, err_mes, max_array)
  !*********************************************************************************************
  ! mpimax_r8_array -- Max real8 array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: max_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, max_array, a_len, MPI_REAL8, MPI_MAX, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimax_r8_array

  subroutine mpimin_i4_scalar(loc_num, err_mes, min_num)
  !*********************************************************************************************
  ! mpimin_i4_scalar -- Min integer value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_num
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: min_num
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_num, min_num, 1, MPI_INTEGER, MPI_MIN, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimin_i4_scalar

  subroutine mpimin_r4_scalar(loc_val, err_mes, min_val)
  !*********************************************************************************************
  ! mpimin_r4_scalar -- Min real4 value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: min_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, min_val, 1, MPI_REAL4, MPI_MIN, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimin_r4_scalar

  subroutine mpimin_r8_scalar(loc_val, err_mes, min_val)
  !*********************************************************************************************
  ! mpimin_r8_scalar -- Min real value for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_val
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: min_val
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLREDUCE(loc_val, min_val, 1, MPI_REAL8, MPI_MIN, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" value in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimin_r8_scalar

  subroutine mpimin_i4_array(loc_array, err_mes, min_array)
  !*********************************************************************************************
  ! mpimin_i4_array -- Min integer array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: min_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, min_array, a_len, MPI_INTEGER, MPI_MIN, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimin_i4_array

  subroutine mpimin_r4_array(loc_array, err_mes, min_array)
  !*********************************************************************************************
  ! mpimin_r4_array -- Min real4 array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: min_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, min_array, a_len, MPI_REAL4, MPI_MIN, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimin_r4_array

  subroutine mpimin_r8_array(loc_array, err_mes, min_array)
  !*********************************************************************************************
  ! mpimin_r8_array -- Min real8 array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: min_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_ALLREDUCE(loc_array, min_array, a_len, MPI_REAL8, MPI_MIN, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpimin_r8_array

  subroutine mpiexscan_i4_array(loc_array, err_mes, scan_array)
  !*********************************************************************************************
  ! mpiexscan_i4_array -- Exscan integer array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: scan_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_EXSCAN(loc_array, scan_array, a_len, MPI_INTEGER, MPI_SUM, st_mpi%comm, ierr)
    if (st_mpi%rank == 0) then
      scan_array(:) = 0
    end if
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Exscan "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpiexscan_i4_array

  subroutine mpiexscan_r4_array(loc_array, err_mes, scan_array)
  !*********************************************************************************************
  ! mpiexscan_r4_array -- Exscan real4 array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    integer(SP), intent(out) :: scan_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_EXSCAN(loc_array, scan_array, a_len, MPI_REAL4, MPI_SUM, st_mpi%comm, ierr)
    if (st_mpi%rank == 0) then
      scan_array(:) = 0
    end if
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Exscan "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpiexscan_r4_array

  subroutine mpiexscan_r8_array(loc_array, err_mes, scan_array)
  !*********************************************************************************************
  ! mpiexscan_r8_array -- Exscan real8 array for MPI
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: loc_array(:)
    character(*), intent(in) :: err_mes
    integer(DP), intent(out) :: scan_array(:)
    ! -- local
    integer(I4) :: a_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(loc_array(:))
    call MPI_EXSCAN(loc_array, scan_array, a_len, MPI_REAL8, MPI_SUM, st_mpi%comm, ierr)
    if (st_mpi%rank == 0) then
      scan_array(:) = 0
    end if
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Exscan "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine mpiexscan_r8_array

  subroutine bcast_i4_scalar(iscalar, err_mes)
  !*********************************************************************************************
  ! bcast_i4_scalar -- Bcast integer scalar value
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: iscalar
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(iscalar, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_i4_scalar

  subroutine bcast_i4_array(iarray, err_mes)
  !*********************************************************************************************
  ! bcast_i4_array -- Bcast integer array value
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: iarray(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr, a_len
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(iarray(:))
    call MPI_BCAST(iarray, a_len, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_i4_array

  subroutine bcast_r4_scalar(rscalar, err_mes)
  !*********************************************************************************************
  ! bcast_r4_scalar -- Bcast real4 scalar value
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(inout) :: rscalar
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(rscalar, 1, MPI_REAL4, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_r4_scalar

  subroutine bcast_r4_array(rarray, err_mes)
  !*********************************************************************************************
  ! bcast_r4_array -- Bcast real4 array value
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(inout) :: rarray(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr, a_len
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(rarray(:))
    call MPI_BCAST(rarray, a_len, MPI_REAL4, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_r4_array

  subroutine bcast_r8_scalar(rscalar, err_mes)
  !*********************************************************************************************
  ! bcast_r8_scalar -- Bcast real8 scalar value
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(inout) :: rscalar
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_BCAST(rscalar, 1, MPI_REAL8, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_r8_scalar

  subroutine bcast_r8_array(rarray, err_mes)
  !*********************************************************************************************
  ! bcast_r8_array -- Bcast real8 array value
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(inout) :: rarray(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: ierr, a_len
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; a_len = size(rarray(:))
    call MPI_BCAST(rarray, a_len, MPI_REAL8, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_r8_array

  subroutine bcast_char(err_mes)
  !*********************************************************************************************
  ! bcast_char -- Bcast character
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(inout) :: err_mes
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: rank_num, send_num
    integer(I4) :: str_len
    !-------------------------------------------------------------------------------------------
    rank_num = 0 ; send_num = 0
    if (len_trim(adjustl(err_mes)) /= 0) then
      rank_num = st_mpi%rank
    end if

    ierr = 0
    call MPI_ALLREDUCE(rank_num, send_num, 1, MPI_INTEGER, MPI_MAX, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allreduce rank number in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    str_len = len_trim(err_mes)
    call MPI_BCAST(str_len, 1, MPI_INTEGER, send_num, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path length in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    call MPI_BCAST(err_mes, str_len, MPI_CHARACTER, send_num, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast message in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_char

  subroutine bcast_file_path(file_path, err_mes)
  !*********************************************************************************************
  ! bcast_file_path -- Bcast file path
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(inout) :: file_path
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(file_path)
    call MPI_BCAST(str_len, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path length in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    call MPI_BCAST(file_path, str_len, MPI_CHARACTER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_file_path

  subroutine bcast_ftype_path(file_type, file_path, err_mes)
  !*********************************************************************************************
  ! bcast_ftype_path -- Bcast file type and path
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: file_type
    character(*), intent(inout) :: file_path
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(file_path)
    call MPI_BCAST(file_type, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file type in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    call MPI_BCAST(str_len, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path length in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    call MPI_BCAST(file_path, str_len, MPI_CHARACTER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_ftype_path

  subroutine bcast_path_unit(file_path, file_unit, err_mes)
  !*********************************************************************************************
  ! bcast_path_unit -- Bcast file path and convert unit
  !*********************************************************************************************
    ! -- module

    ! -- inout
    character(*), intent(inout) :: file_path, file_unit
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len, uni_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(file_path) ; uni_len = len_trim(file_unit)
    call MPI_BCAST(str_len, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path length in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    call MPI_BCAST(file_path, str_len, MPI_CHARACTER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    call MPI_BCAST(uni_len, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file unit length in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    call MPI_BCAST(file_unit, uni_len, MPI_CHARACTER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file unit in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_path_unit

  subroutine bcast_ftype_path_unit(file_type, file_path, file_unit, err_mes)
  !*********************************************************************************************
  ! bcast_ftype_path_unit -- Bcast file type, path and convert unit
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: file_type
    character(*), intent(inout) :: file_path, file_unit
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len, uni_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(file_path) ; uni_len = len_trim(file_unit)
    call MPI_BCAST(file_type, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file type in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    call MPI_BCAST(str_len, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path length in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    call MPI_BCAST(file_path, str_len, MPI_CHARACTER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    call MPI_BCAST(uni_len, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file unit length in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    call MPI_BCAST(file_unit, uni_len, MPI_CHARACTER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file unit in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_ftype_path_unit

  subroutine bcast_extr_set(extr_type, extr_step, extr_end, extr_path, err_mes)
  !*********************************************************************************************
  ! bcast_extr_set -- Bcast extra setting
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(inout) :: extr_type
    real(SP), intent(inout) :: extr_step, extr_end
    character(*), intent(inout) :: extr_path
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: str_len
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; str_len = len_trim(extr_path)

    call MPI_BCAST(extr_type, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file type in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    call MPI_BCAST(extr_step, 1, MPI_REAL4, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file step in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    call MPI_BCAST(extr_end, 1, MPI_REAL4, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file end time in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    call MPI_BCAST(str_len, 1, MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path length in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    call MPI_BCAST(extr_path, str_len, MPI_CHARACTER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Broadcast "//err_mes//" file path in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine bcast_extr_set

  subroutine gather_i4_array(num_prot, in_num, loc_array, glo_array, err_mes)
  !*********************************************************************************************
  ! gather_i4_array -- Gather integer array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_prot, in_num
    integer(I4), intent(in) :: loc_array(:)
    integer(I4), intent(out) :: glo_array(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: i, sum_num, ierr
    integer(I4), allocatable :: rec_num(:), rec_count(:), rec_dis(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(rec_num(num_prot))
    !$omp parallel do private(i)
    do i = 1, num_prot
      rec_num(i) = 0
    end do
    !$omp end parallel do
    call MPI_ALLGATHER(in_num, 1, MPI_INTEGER, rec_num, 1, MPI_INTEGER, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    allocate(rec_count(0:num_prot-1), rec_dis(0:num_prot-1))
    rec_count(0) = 0 ; rec_dis(0) = 0
    !$omp parallel do private(i)
    do i = 1, num_prot-1
      rec_count(i) = 0
      rec_dis(i) = 0
    end do
    !$omp end parallel do
    sum_num = 0
    do i = 1, num_prot
      rec_dis(i-1) = sum_num
      sum_num = sum_num + rec_num(i)
      rec_count(i-1) = rec_num(i)
    end do

    call MPI_ALLGATHERV(loc_array, in_num, MPI_INTEGER, glo_array, rec_count, rec_dis,&
                        MPI_INTEGER, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    deallocate(rec_num, rec_count, rec_dis)

  end subroutine gather_i4_array

  subroutine scatterv_i4_array(num_prot, out_num, glo_array, loc_array, err_mes)
  !*********************************************************************************************
  ! scatterv_i4_array -- Scatter an integer array held on rank 0 to each rank's contiguous
  !   range (the inverse of gather_i4_array). Only rank 0 needs glo_array to be valid.
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_prot, out_num
    integer(I4), intent(in) :: glo_array(:)
    integer(I4), intent(out) :: loc_array(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: i, sum_num, ierr
    integer(I4), allocatable :: sen_num(:), sen_count(:), sen_dis(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(sen_num(num_prot))
    !$omp parallel do private(i)
    do i = 1, num_prot
      sen_num(i) = 0
    end do
    !$omp end parallel do
    call MPI_ALLGATHER(out_num, 1, MPI_INTEGER, sen_num, 1, MPI_INTEGER, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" size in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    allocate(sen_count(0:num_prot-1), sen_dis(0:num_prot-1))
    sum_num = 0
    do i = 1, num_prot
      sen_dis(i-1) = sum_num
      sum_num = sum_num + sen_num(i)
      sen_count(i-1) = sen_num(i)
    end do

    call MPI_SCATTERV(glo_array, sen_count, sen_dis, MPI_INTEGER, loc_array, out_num,&
                      MPI_INTEGER, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Scatterv "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    deallocate(sen_num, sen_count, sen_dis)

  end subroutine scatterv_i4_array

  subroutine scatterv_r4_array(num_prot, out_num, glo_array, loc_array, err_mes)
  !*********************************************************************************************
  ! scatterv_r4_array -- Scatter a real4 array held on rank 0 to each rank's contiguous range
  !   (the real4 counterpart of scatterv_i4_array). Only rank 0 needs glo_array to be valid.
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_prot, out_num
    real(SP), intent(in) :: glo_array(:)
    real(SP), intent(out) :: loc_array(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: i, sum_num, ierr
    integer(I4), allocatable :: sen_num(:), sen_count(:), sen_dis(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(sen_num(num_prot))
    !$omp parallel do private(i)
    do i = 1, num_prot
      sen_num(i) = 0
    end do
    !$omp end parallel do
    call MPI_ALLGATHER(out_num, 1, MPI_INTEGER, sen_num, 1, MPI_INTEGER, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" size in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    allocate(sen_count(0:num_prot-1), sen_dis(0:num_prot-1))
    sum_num = 0
    do i = 1, num_prot
      sen_dis(i-1) = sum_num
      sum_num = sum_num + sen_num(i)
      sen_count(i-1) = sen_num(i)
    end do

    call MPI_SCATTERV(glo_array, sen_count, sen_dis, MPI_REAL4, loc_array, out_num,&
                      MPI_REAL4, 0, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Scatterv "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    deallocate(sen_num, sen_count, sen_dis)

  end subroutine scatterv_r4_array

  subroutine gather_r4_array(num_prot, in_num, loc_array, glo_array, err_mes)
  !*********************************************************************************************
  ! gather_r4_array -- Gather real4 array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_prot, in_num
    real(SP), intent(in) :: loc_array(:)
    real(SP), intent(out) :: glo_array(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: i, sum_num, ierr
    integer(I4), allocatable :: rec_num(:), rec_count(:), rec_dis(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(rec_num(num_prot))
    !$omp parallel do private(i)
    do i = 1, num_prot
      rec_num(i) = 0
    end do
    !$omp end parallel do
    call MPI_ALLGATHER(in_num, 1, MPI_INTEGER, rec_num, 1, MPI_INTEGER, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    allocate(rec_count(0:num_prot-1), rec_dis(0:num_prot-1))
    rec_count(0) = 0 ; rec_dis(0) = 0
    !$omp parallel do private(i)
    do i = 1, num_prot-1
      rec_count(i) = 0
      rec_dis(i) = 0
    end do
    !$omp end parallel do
    sum_num = 0
    do i = 1, num_prot
      rec_dis(i-1) = sum_num
      sum_num = sum_num + rec_num(i)
      rec_count(i-1) = rec_num(i)
    end do

    call MPI_ALLGATHERV(loc_array, in_num, MPI_REAL4, glo_array, rec_count, rec_dis,&
                        MPI_REAL4, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    deallocate(rec_num, rec_count, rec_dis)

  end subroutine gather_r4_array

  subroutine gather_r8_array(num_prot, in_num, loc_array, glo_array, err_mes)
  !*********************************************************************************************
  ! gather_r8_array -- Gather real8 array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_prot, in_num
    real(DP), intent(in) :: loc_array(:)
    real(DP), intent(out) :: glo_array(:)
    character(*), intent(in) :: err_mes
    ! -- local
    integer(I4) :: i, sum_num, ierr
    integer(I4), allocatable :: rec_num(:), rec_count(:), rec_dis(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(rec_num(num_prot))
    !$omp parallel do private(i)
    do i = 1, num_prot
      rec_num(i) = 0
    end do
    !$omp end parallel do
    call MPI_ALLGATHER(in_num, 1, MPI_INTEGER, rec_num, 1, MPI_INTEGER, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    allocate(rec_count(0:num_prot-1), rec_dis(0:num_prot-1))
    rec_count(0) = 0 ; rec_dis(0) = 0
    !$omp parallel do private(i)
    do i = 1, num_prot-1
      rec_count(i) = 0
      rec_dis(i) = 0
    end do
    !$omp end parallel do
    sum_num = 0
    do i = 1, num_prot
      rec_dis(i-1) = sum_num
      sum_num = sum_num + rec_num(i)
      rec_count(i-1) = rec_num(i)
    end do

    call MPI_ALLGATHERV(loc_array, in_num, MPI_REAL8, glo_array, rec_count, rec_dis,&
                        MPI_REAL8, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Allgather "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

    deallocate(rec_num, rec_count, rec_dis)

  end subroutine gather_r8_array

  subroutine alltoall_i4_array(send_array, err_mes, recv_array)
  !*********************************************************************************************
  ! alltoall_i4_array -- AlltoAll integer array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: send_array(:)
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: recv_array(:)
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLTOALL(send_array, 1, MPI_INTEGER, recv_array, 1, MPI_INTEGER, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Alltoall "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine alltoall_i4_array

  subroutine alltoall_r4_array(send_array, err_mes, recv_array)
  !*********************************************************************************************
  ! alltoall_r4_array -- AlltoAll real4 array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: send_array(:)
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: recv_array(:)
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLTOALL(send_array, 1, MPI_REAL, recv_array, 1, MPI_REAL, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Alltoall "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine alltoall_r4_array

  subroutine alltoall_r8_array(send_array, err_mes, recv_array)
  !*********************************************************************************************
  ! alltoall_r8_array -- AlltoAll real8 array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: send_array(:)
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: recv_array(:)
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_ALLTOALL(send_array, 1, MPI_REAL8, recv_array, 1, MPI_REAL8, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Alltoall "//err_mes//" in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if

  end subroutine alltoall_r8_array

  subroutine alltoallv_i4_array(send_buf, send_count, recv_count, err_mes, recv_buf)
  !*********************************************************************************************
  ! alltoallv_i4_array -- AlltoAll variable integer array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: send_buf(:)
    integer(I4), intent(in) :: send_count(:)
    integer(I4), intent(in) :: recv_count(:)
    character(*), intent(in) :: err_mes
    integer(I4), intent(out) :: recv_buf(:)
    ! -- local
    integer(I4) :: proc_id, send_pos, recv_pos, ierr
    integer(I4), allocatable :: send_disp(:), recv_disp(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(send_disp(0:st_mpi%totn-1), recv_disp(0:st_mpi%totn-1))
    send_pos = 0 ; recv_pos = 0
    do proc_id = 1, st_mpi%totn
      send_disp(proc_id-1) = send_pos ; send_pos = send_pos + send_count(proc_id)
      recv_disp(proc_id-1) = recv_pos ; recv_pos = recv_pos + recv_count(proc_id)
    end do
    call MPI_ALLTOALLV(send_buf, send_count, send_disp, MPI_INTEGER, recv_buf, recv_count,&
                       recv_disp, MPI_INTEGER, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Alltoallv "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    deallocate(send_disp, recv_disp)

  end subroutine alltoallv_i4_array

  subroutine alltoallv_r4_array(send_buf, send_count, recv_count, err_mes, recv_buf)
  !*********************************************************************************************
  ! alltoallv_r4_array -- AlltoAll variable real4 array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(SP), intent(in) :: send_buf(:)
    integer(I4), intent(in) :: send_count(:)
    integer(I4), intent(in) :: recv_count(:)
    character(*), intent(in) :: err_mes
    real(SP), intent(out) :: recv_buf(:)
    ! -- local
    integer(I4) :: proc_id, send_pos, recv_pos, ierr
    integer(I4), allocatable :: send_disp(:), recv_disp(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(send_disp(0:st_mpi%totn-1), recv_disp(0:st_mpi%totn-1))
    send_pos = 0 ; recv_pos = 0
    do proc_id = 1, st_mpi%totn
      send_disp(proc_id-1) = send_pos ; send_pos = send_pos + send_count(proc_id)
      recv_disp(proc_id-1) = recv_pos ; recv_pos = recv_pos + recv_count(proc_id)
    end do
    call MPI_ALLTOALLV(send_buf, send_count, send_disp, MPI_REAL, recv_buf, recv_count,&
                       recv_disp, MPI_REAL, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Alltoallv "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    deallocate(send_disp, recv_disp)

  end subroutine alltoallv_r4_array

  subroutine alltoallv_r8_array(send_buf, send_count, recv_count, err_mes, recv_buf)
  !*********************************************************************************************
  ! alltoallv_r8_array -- AlltoAll variable real8 array
  !*********************************************************************************************
    ! -- module

    ! -- inout
    real(DP), intent(in) :: send_buf(:)
    integer(I4), intent(in) :: send_count(:)
    integer(I4), intent(in) :: recv_count(:)
    character(*), intent(in) :: err_mes
    real(DP), intent(out) :: recv_buf(:)
    ! -- local
    integer(I4) :: proc_id, send_pos, recv_pos, ierr
    integer(I4), allocatable :: send_disp(:), recv_disp(:)
    !-------------------------------------------------------------------------------------------
    ierr = 0
    allocate(send_disp(0:st_mpi%totn-1), recv_disp(0:st_mpi%totn-1))
    send_pos = 0 ; recv_pos = 0
    do proc_id = 1, st_mpi%totn
      send_disp(proc_id-1) = send_pos ; send_pos = send_pos + send_count(proc_id)
      recv_disp(proc_id-1) = recv_pos ; recv_pos = recv_pos + recv_count(proc_id)
    end do
    call MPI_ALLTOALLV(send_buf, send_count, send_disp, MPI_REAL8, recv_buf, recv_count,&
                       recv_disp, MPI_REAL8, st_mpi%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (st_mpi%rank == 0) then
        write(log_fnum,'(a)') "Error!! Alltoallv "//err_mes//" array in MPI program."
      end if
      call abort_proc(st_mpi%rank, log_fnum)
    end if
    deallocate(send_disp, recv_disp)

  end subroutine alltoallv_r8_array

end module mpi_utility
