module mpi_initfin
  ! -- modules
  use kind_module, only: I4
  use mpi

  implicit none
  private
  public :: init_mpi, fin_mpi, abort_proc
  ! -- local
  integer(I4) :: my_comm

  contains

  subroutine init_mpi(num_log, mpi_env)
  !*********************************************************************************************
  ! init_mpi -- Initialize mpi
  !*********************************************************************************************
    ! -- module
    use types_module, only: mpi_set

    ! -- inout
    integer(I4), intent(in) :: num_log
    type(mpi_set), intent(inout) :: mpi_env
    ! -- local
    integer(I4) :: ierr, provided
    logical :: mpi_init_check
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_INITIALIZED(mpi_init_check, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (mpi_env%rank == 0) then
        write(num_log,'(a)') "Error!! Initialize MPI program."
      end if
      call abort_proc(mpi_env%rank, num_log)
    end if

    if (.not. mpi_init_check) then
      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, ierr)
      my_comm = MPI_COMM_WORLD
      mpi_env%comm = my_comm
      if (ierr /= MPI_SUCCESS) then
        if (mpi_env%rank == 0) then
          write(num_log,'(a)') "Error!! Start MPI program."
        end if
        call abort_proc(mpi_env%rank, num_log)
      end if
      if (provided < MPI_THREAD_MULTIPLE) then
        if (mpi_env%rank == 0) then
          write(num_log,'(a)') "Error!! MPI_THREAD_MULTIPLE not supported."
        end if
        call abort_proc(mpi_env%rank, num_log)
      end if
    end if

    call MPI_COMM_SIZE(my_comm, mpi_env%totn, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (mpi_env%rank == 0) then
        write(num_log,'(a)') "Error!! Get common size in MPI program."
      end if
      call abort_proc(mpi_env%rank, num_log)
    end if

    call MPI_COMM_RANK(my_comm, mpi_env%rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (mpi_env%rank == 0) then
        write(num_log,'(a)') "Error!! Get common rank in MPI program."
      end if
      call abort_proc(mpi_env%rank, num_log)
    end if

  end subroutine init_mpi

  subroutine fin_mpi(num_rank, num_log)
  !*********************************************************************************************
  ! fin_mpi -- Finalize mpi
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_rank, num_log
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------------
    ierr = 0
    call MPI_FINALIZE(ierr)
    if (ierr /= MPI_SUCCESS) then
      if (num_rank == 0) then
        write(num_log,'(a)') "Error!! Finalize MPI program."
      end if
    end if

  end subroutine fin_mpi

  subroutine abort_proc(num_rank, num_log)
  !*********************************************************************************************
  ! abort_proc -- Abort process
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_rank, num_log
    ! -- local
    integer(I4) :: ierr, errcode
    !-------------------------------------------------------------------------------------------
    ierr = 0 ; errcode = 1
    call MPI_ABORT(my_comm, errcode, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (num_rank == 0) then
        write(num_log,'(a)') "Error!! Abort MPI program."
      end if
    end if

    call MPI_FINALIZE(ierr)
    if (ierr /= MPI_SUCCESS) then
      if (num_rank == 0) then
        write(num_log,'(a)') "Error!! Finalize MPI program."
      end if
    end if

  end subroutine abort_proc

end module mpi_initfin
