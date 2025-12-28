module mpi_initfin
  ! -- modules
  use kind_module, only: I4
  use mpi

  implicit none
  private
  integer(I4), public :: my_comm
  public :: init_mpi, fin_mpi, abort_proc
  ! -- local

  contains

  subroutine init_mpi(num_log, num_prot, num_rank)
  !***************************************************************************************
  ! init_mpi -- Initialize mpi
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_log
    integer(I4), intent(out) :: num_prot, num_rank
    ! -- local
    integer(I4) :: ierr
    logical :: mpi_init_check
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_INITIALIZED(mpi_init_check, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (num_rank == 0) then
        write(num_log,'(a)') "Error!! Initialize MPI program."
      end if
      call abort_proc(num_rank, num_log)
    end if

    if (.not. mpi_init_check) then
      call MPI_INIT(ierr)
      my_comm = MPI_COMM_WORLD
      if (ierr /= MPI_SUCCESS) then
        if (num_rank == 0) then
          write(num_log,'(a)') "Error!! Start MPI program."
        end if
        call abort_proc(num_rank, num_log)
      end if
    end if

    call MPI_COMM_SIZE(my_comm, num_prot, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (num_rank == 0) then
        write(num_log,'(a)') "Error!! Get common size in MPI program."
      end if
      call abort_proc(num_rank, num_log)
    end if

    call MPI_COMM_RANK(my_comm, num_rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      if (num_rank == 0) then
        write(num_log,'(a)') "Error!! Get common rank in MPI program."
      end if
      call abort_proc(num_rank, num_log)
    end if

  end subroutine init_mpi

  subroutine fin_mpi(num_rank, num_log)
  !***************************************************************************************
  ! fin_mpi -- Finalize mpi
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_rank, num_log
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    call MPI_FINALIZE(ierr)
    if (ierr /= MPI_SUCCESS) then
      if (num_rank == 0) then
        write(num_log,'(a)') "Error!! Finalize MPI program."
      end if
    end if

  end subroutine fin_mpi

  subroutine abort_proc(num_rank, num_log)
  !***************************************************************************************
  ! abort_proc -- Abort process
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: num_rank, num_log
    ! -- local
    integer(I4) :: ierr, errcode
    !-------------------------------------------------------------------------------------
    ierr = 0
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
