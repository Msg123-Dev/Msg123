module write_module
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SNOVAL
  use utility_module, only: write_err_write
  use initial_module, only: st_grid
  use set_cell, only: get_cals_grid, get_calc_grid

  implicit none
  private
  public :: write_2dtxt, write_2dbin, write_3dtxt, write_3dbin
  public :: write_header_txt, write_header_bin

  ! -- local

  contains

  subroutine write_2dtxt(fnum, out_format, out_totn, calc_num, out_unit, out_val)
  !***************************************************************************************
  ! write_2dtxt -- Write 2D text file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, out_totn
    character(*) :: out_format
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, j, ierr
    integer(I4) :: xnum, ynum
    real(SP), allocatable :: array_out(:,:)
    !-------------------------------------------------------------------------------------
    allocate(array_out(st_grid%nx,st_grid%ny))
    !$omp parallel
    !$omp workshare
    array_out(:,:) = SNOVAL
    !$omp end workshare

    !$omp do private(i, xnum, ynum)
    do i = 1, out_totn
      call get_cals_grid(calc_num(i), xnum, ynum)
      array_out(xnum,ynum) = real(out_val(i)*out_unit, kind=SP)
    end do
    !$omp end do
    !$omp end parallel

    ierr = 0
    do j = 1, st_grid%ny
      write(unit=fnum,fmt=out_format,iostat=ierr) (array_out(i,j), i = 1, st_grid%nx)
      if (ierr /= 0) then
        call write_err_write(fnum)
      end if
    end do

    deallocate(array_out)

  end subroutine write_2dtxt

  subroutine write_2dbin(fnum, out_totn, calc_num, out_unit, out_val)
  !***************************************************************************************
  ! write_2dbin -- Write 2D binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, out_totn
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, j, ierr
    integer(I4) :: xnum, ynum
    real(SP), allocatable :: array_out(:,:)
    !-------------------------------------------------------------------------------------
    allocate(array_out(st_grid%nx,st_grid%ny))
    !$omp parallel
    !$omp workshare
    array_out(:,:) = SNOVAL
    !$omp end workshare

    !$omp do private(i, xnum, ynum)
    do i = 1, out_totn
      call get_cals_grid(calc_num(i), xnum, ynum)
      array_out(xnum,ynum) = real(out_val(i)*out_unit, kind=SP)
    end do
    !$omp end do
    !$omp end parallel

    ierr = 0
    write(unit=fnum,iostat=ierr) ((array_out(i,j), i = 1, st_grid%nx), j = 1, st_grid%ny)
    if (ierr /= 0) then
      call write_err_write(fnum)
    end if

    deallocate(array_out)

  end subroutine write_2dbin

  subroutine write_3dtxt(fnum, out_format, out_totn, calc_num, out_unit, out_val)
  !***************************************************************************************
  ! write_3dtxt -- Write 3D text file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, out_totn
    character(*) :: out_format
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, j, k, ierr
    integer(I4) :: xnum, ynum, znum
    real(SP), allocatable :: array_out(:,:,:)
    !-------------------------------------------------------------------------------------
    allocate(array_out(st_grid%nx,st_grid%ny,st_grid%nz))
    !$omp parallel
    !$omp workshare
    array_out(:,:,:) = SNOVAL
    !$omp end workshare

    !$omp do private(i, xnum, ynum, znum)
    do i = 1, out_totn
      call get_calc_grid(calc_num(i), xnum, ynum, znum)
      array_out(xnum,ynum,znum) = real(out_val(i)*out_unit, kind=SP)
    end do
    !$omp end do
    !$omp end parallel

    ierr = 0
    do k = 1, st_grid%nz
      do j = 1, st_grid%ny
        write(unit=fnum,fmt=out_format,iostat=ierr) (array_out(i,j,k), i = 1, st_grid%nx)
        if (ierr /= 0) then
          call write_err_write(fnum)
        end if
      end do
    end do

    deallocate(array_out)

  end subroutine write_3dtxt

  subroutine write_3dbin(fnum, out_totn, calc_num, out_unit, out_val)
  !***************************************************************************************
  ! write_3dbin -- Write 3D binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, out_totn
    integer(I4), intent(in) :: calc_num(:)
    real(SP), intent(in) :: out_unit
    real(DP), intent(in) :: out_val(:)
    ! -- local
    integer(I4) :: i, j, k, ierr
    integer(I4) :: xnum, ynum, znum
    real(SP), allocatable :: array_out(:,:,:)
    !-------------------------------------------------------------------------------------
    allocate(array_out(st_grid%nx,st_grid%ny,st_grid%nz))
    !$omp parallel
    !$omp workshare
    array_out(:,:,:) = SNOVAL
    !$omp end workshare

    !$omp do private(i, xnum, ynum, znum)
    do i = 1, out_totn
      call get_calc_grid(calc_num(i), xnum, ynum, znum)
      array_out(xnum,ynum,znum) = real(out_val(i)*out_unit, kind=SP)
    end do
    !$omp end do
    !$omp end parallel

    ierr = 0
    write(unit=fnum,iostat=ierr) (((array_out(i,j,k), i = 1, st_grid%nx), j = 1, st_grid%ny),&
                                  k = 1, st_grid%nz)
    if (ierr /= 0) then
      call write_err_write(fnum)
    end if

    deallocate(array_out)

  end subroutine write_3dbin

  subroutine write_header_txt(fnum, out_time)
  !***************************************************************************************
  ! write_header_txt -- Write header text file
  !***************************************************************************************
    ! -- modules
    use constval_module, only: OUTFORM
    ! -- inout
    integer(I4), intent(in) :: fnum
    real(SP), intent(in) :: out_time
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    write(unit=fnum,fmt="("//OUTFORM//")",iostat=ierr) out_time
    if (ierr /= 0) then
      call write_err_write(fnum)
    end if

  end subroutine write_header_txt

  subroutine write_header_bin(fnum, out_time)
  !***************************************************************************************
  ! write_header_txt -- Write header binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum
    real(SP), intent(in) :: out_time
    ! -- local
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    write(unit=fnum,iostat=ierr) out_time
    if (ierr /= 0) then
      call write_err_write(fnum)
    end if

  end subroutine write_header_bin

end module write_module
