module read_module
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: CHALEN
  use utility_module, only: open_new_rtxt, open_new_rbin, write_err_read, write_logf
  use initial_module, only: st_sim, st_grid, st_init, my_rank, in_type
#ifdef MPI_MSG
  use mpi_utility, only: bcast_val
  use utility_module, only: log_fnum
  use initial_module, only: pro_totn
#endif

  implicit none
  private
  public :: read_clasf, read_2dpointf, read_3dpointf, read_wpointf
  public :: read_2dtxt, read_2dbin, read_3dtxt, read_3dbin
  public :: skip_file, skip_file_int, flat_2dto2d, flat_2dto3d, flat_3dto3d
  public :: read_next, read_intn, read_2d_calcreg, read_3d_calcreg

  interface read_2dtxt
    module procedure read_2d_i4txt
    module procedure read_2d_r4txt
    module procedure read_2d_r8txt
  end interface

  interface read_2dbin
    module procedure read_2d_i4bin
    module procedure read_2d_r4bin
    module procedure read_2d_r8bin
  end interface

  interface read_3dtxt
    module procedure read_3d_i4txt
    module procedure read_3d_r4txt
    module procedure read_3d_r8txt
  end interface

  interface read_3dbin
    module procedure read_3d_i4bin
    module procedure read_3d_r4bin
    module procedure read_3d_r8bin
  end interface

  interface flat_2dto2d
    module procedure flat_2d_2di4
    module procedure flat_2d_2dr4
    module procedure flat_2d_2dr8
  end interface

  interface flat_2dto3d
    module procedure flat_2d_3di4
    module procedure flat_2d_3dr4
    module procedure flat_2d_3dr8
  end interface

  interface flat_3dto3d
    module procedure flat_3d_3di4
    module procedure flat_3d_3dr4
    module procedure flat_3d_3dr8
  end interface

  ! -- local

  contains

  subroutine read_clasf(fnum, targn, targ_name, targ_val)
  !***************************************************************************************
  ! read_clasf -- Read classification input file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, targn
    character(*), intent(out) :: targ_name(:)
    real(SP), intent(out) :: targ_val(:)
    ! -- local
    integer(I4) :: i
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do i = 1, targn
      read(unit=fnum,fmt=*,iostat=ierr) targ_name(i), targ_val(i)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine read_clasf

  subroutine read_2dpointf(fnum, targn, p_i, p_j, targ_val)
  !***************************************************************************************
  ! read_2dpointf -- Read 2d point input file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, targn
    integer(I4), intent(out) :: p_i(:), p_j(:)
    real(SP), intent(out) :: targ_val(:)
    ! -- local
    integer(I4) :: i
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do i = 1, targn
      read(unit=fnum,fmt=*,iostat=ierr) p_i(i), p_j(i), targ_val(i)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine read_2dpointf

  subroutine read_3dpointf(fnum, targn, p_i, p_j, p_k, targ_val)
  !***************************************************************************************
  ! read_3dpointf -- Read 3d point input file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, targn
    integer(I4), intent(out) :: p_i(:), p_j(:), p_k(:)
    real(SP), intent(out) :: targ_val(:)
    ! -- local
    integer(I4) :: i
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do i = 1, targn
      read(unit=fnum,fmt=*,iostat=ierr) p_i(i), p_j(i), p_k(i), targ_val(i)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine read_3dpointf

  subroutine read_wpointf(fnum, wnum, pid, pi, pj, pks, pke, pval)
  !***************************************************************************************
  ! read_wpointf -- Read 3d well point input file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, wnum
    integer(I4), intent(out) :: pid(:), pi(:), pj(:), pks(:), pke(:)
    real(SP), intent(out) :: pval(:)
    ! -- local
    integer(I4) :: i
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do i = 1, wnum
      read(unit=fnum,fmt=*,iostat=ierr) pid(i), pi(i), pj(i), pks(i), pke(i), pval(i)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine read_wpointf

  subroutine read_2d_i4txt(fnum, inx, iny, array_val)
  !***************************************************************************************
  ! read_2d_i4txt -- Read 2D integer text file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny
    integer(I4), intent(out) :: array_val(:,:)
    ! -- local
    integer(I4) :: i, j
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do j = 1, iny
      read(unit=fnum,fmt=*,iostat=ierr) (array_val(i,j), i = 1, inx)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine read_2d_i4txt

  subroutine read_2d_r4txt(fnum, inx, iny, array_val)
  !***************************************************************************************
  ! read_2d_r4txt -- Read 2D real4 text file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny
    real(SP), intent(out) :: array_val(:,:)
    ! -- local
    integer(I4) :: i, j
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do j = 1, iny
      read(unit=fnum,fmt=*,iostat=ierr) (array_val(i,j), i = 1, inx)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine read_2d_r4txt

  subroutine read_2d_r8txt(fnum, inx, iny, array_val)
  !***************************************************************************************
  ! read_2d_r8txt -- Read 2D real8 text file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny
    real(DP), intent(out) :: array_val(:,:)
    ! -- local
    integer(I4) :: i, j
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do j = 1, iny
      read(unit=fnum,fmt=*,iostat=ierr) (array_val(i,j), i = 1, inx)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine read_2d_r8txt

  subroutine read_2d_i4bin(fnum, inx, iny, no_val, array_val)
  !***************************************************************************************
  ! read_2d_i4bin -- Read 2D integer binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny
    integer(I4), intent(in) :: no_val
    integer(I4), intent(out) :: array_val(:,:)
    ! -- local
    integer(I4) :: i, j
    integer(I4) :: ierr
    integer(I4), allocatable :: array_in(:,:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(array_in(inx,iny))
    !$omp parallel workshare
    array_in(:,:) = no_val
    !$omp end parallel workshare

    read(unit=fnum,iostat=ierr) ((array_in(i,j), i = 1, inx), j = 1, iny)

    if (ierr /= 0) then
      call write_err_read(fnum)
    end if

    !$omp parallel workshare
    array_val(:,:) = array_in(:,:)
    !$omp end parallel workshare

    deallocate(array_in)

  end subroutine read_2d_i4bin

  subroutine read_2d_r4bin(fnum, inx, iny, no_val, array_val)
  !***************************************************************************************
  ! read_2d_r4bin -- Read 2D real4 binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny
    real(SP), intent(in) :: no_val
    real(SP), intent(out) :: array_val(:,:)
    ! -- local
    integer(I4) :: i, j
    integer(I4) :: ierr
    real(SP), allocatable :: array_in(:,:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(array_in(inx,iny))
    !$omp parallel workshare
    array_in(:,:) = no_val
    !$omp end parallel workshare

    read(unit=fnum,iostat=ierr) ((array_in(i,j), i = 1, inx), j = 1, iny)

    if (ierr /= 0) then
      call write_err_read(fnum)
    end if

    !$omp parallel workshare
    array_val(:,:) = array_in(:,:)
    !$omp end parallel workshare

    deallocate(array_in)

  end subroutine read_2d_r4bin

  subroutine read_2d_r8bin(fnum, inx, iny, no_val, array_val)
  !***************************************************************************************
  ! read_2d_r8bin -- Read 2D real8 binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny
    real(DP), intent(in) :: no_val
    real(DP), intent(out) :: array_val(:,:)
    ! -- local
    integer(I4) :: i, j
    integer(I4) :: ierr
    real(DP), allocatable :: array_in(:,:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(array_in(inx,iny))
    !$omp parallel workshare
    array_in(:,:) = no_val
    !$omp end parallel workshare

    read(unit=fnum,iostat=ierr) ((array_in(i,j), i = 1, inx), j = 1, iny)

    if (ierr /= 0) then
      call write_err_read(fnum)
    end if

    !$omp parallel workshare
    array_val(:,:) = array_in(:,:)
    !$omp end parallel workshare

    deallocate(array_in)

  end subroutine read_2d_r8bin

  subroutine read_3d_i4txt(fnum, inx, iny, inz, array_val)
  !***************************************************************************************
  ! read_3d_i4txt -- Read 3D integer text file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny, inz
    integer(I4), intent(out) :: array_val(:,:,:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do k = 1, inz
      do j = 1, iny
        read(unit=fnum,fmt=*,iostat=ierr) (array_val(i,j,k), i = 1, inx)
        if (ierr /= 0) then
          call write_err_read(fnum)
        end if
      end do
    end do

  end subroutine read_3d_i4txt

  subroutine read_3d_r4txt(fnum, inx, iny, inz, array_val)
  !***************************************************************************************
  ! read_3d_r4txt -- Read 3D real4 text file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny, inz
    real(SP), intent(out) :: array_val(:,:,:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do k = 1, inz
      do j = 1, iny
        read(unit=fnum,fmt=*,iostat=ierr) (array_val(i,j,k), i = 1, inx)
        if (ierr /= 0) then
          call write_err_read(fnum)
        end if
      end do
    end do

  end subroutine read_3d_r4txt

  subroutine read_3d_r8txt(fnum, inx, iny, inz, array_val)
  !***************************************************************************************
  ! read_3d_r8txt -- Read 3D real8 text file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny, inz
    real(DP), intent(out) :: array_val(:,:,:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do k = 1, inz
      do j = 1, iny
        read(unit=fnum,fmt=*,iostat=ierr) (array_val(i,j,k), i = 1, inx)
        if (ierr /= 0) then
          call write_err_read(fnum)
        end if
      end do
    end do

  end subroutine read_3d_r8txt

  subroutine read_3d_i4bin(fnum, inx, iny, inz, no_val, array_val)
  !***************************************************************************************
  ! read_3d_i4bin -- Read 3D integer binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny, inz
    integer(I4), intent(in) :: no_val
    integer(I4), intent(out) :: array_val(:,:,:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: ierr
    integer(I4), allocatable :: array_in(:,:,:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(array_in(inx,iny,inz))
    !$omp parallel workshare
    array_in(:,:,:) = no_val
    !$omp end parallel workshare

    read(unit=fnum,iostat=ierr) (((array_in(i,j,k), i = 1, inx), j = 1, iny), k = 1, inz)

    if (ierr /= 0) then
      call write_err_read(fnum)
    end if

    !$omp parallel workshare
    array_val(:,:,:) = array_in(:,:,:)
    !$omp end parallel workshare

    deallocate(array_in)

  end subroutine read_3d_i4bin

  subroutine read_3d_r4bin(fnum, inx, iny, inz, no_val, array_val)
  !***************************************************************************************
  ! read_3d_r4bin -- Read 3D real4 binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny, inz
    real(SP), intent(in) :: no_val
    real(SP), intent(out) :: array_val(:,:,:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: ierr
    real(SP), allocatable :: array_in(:,:,:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(array_in(inx,iny,inz))
    !$omp parallel workshare
    array_in(:,:,:) = no_val
    !$omp end parallel workshare

    read(unit=fnum,iostat=ierr) (((array_in(i,j,k), i = 1, inx), j = 1, iny), k = 1, inz)

    if (ierr /= 0) then
      call write_err_read(fnum)
    end if

    !$omp parallel workshare
    array_val(:,:,:) = array_in(:,:,:)
    !$omp end parallel workshare

    deallocate(array_in)

  end subroutine read_3d_r4bin

  subroutine read_3d_r8bin(fnum, inx, iny, inz, no_val, array_val)
  !***************************************************************************************
  ! read_3d_r8bin -- Read 3D real8 binary file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny, inz
    real(DP), intent(in) :: no_val
    real(DP), intent(out) :: array_val(:,:,:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: ierr
    real(DP), allocatable :: array_in(:,:,:)
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(array_in(inx,iny,inz))
    !$omp parallel workshare
    array_in(:,:,:) = no_val
    !$omp end parallel workshare

    read(unit=fnum,iostat=ierr) (((array_in(i,j,k), i = 1, inx), j = 1, iny), k = 1, inz)

    if (ierr /= 0) then
      call write_err_read(fnum)
    end if

    !$omp parallel workshare
    array_val(:,:,:) = array_in(:,:,:)
    !$omp end parallel workshare

    deallocate(array_in)

  end subroutine read_3d_r8bin

  subroutine skip_file(ftype, fnum, fname, fmulti, file_totn, fetime)
  !***************************************************************************************
  ! skip_file -- Skip input file
  !***************************************************************************************
    ! -- modules
    use utility_module, only: write_err_stop
    ! -- inout
    integer(I4), intent(in) :: ftype, fnum
    character(*), intent(in) :: fname
    real(SP), intent(in) :: fmulti
    integer(I4), intent(out) :: file_totn
    real(SP), intent(inout) :: fetime
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: time_flag
    integer(I4), allocatable :: type_txt(:), type_bin(:)
    character(:), allocatable :: err_mes
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(type_txt(2), type_bin(2))
    type_txt(:) = [in_type(3), in_type(5)] ; type_bin(:) = [in_type(4), in_type(6)]

    time_flag = 0
    do while (time_flag == 0)
      if (ftype == in_type(1) .or. ftype == in_type(2)) then
        read(unit=fnum,fmt=*,iostat=ierr) file_totn, fetime
        if (ierr /= 0) then
          call write_err_read(fnum)
        end if
        if (file_totn < 0) then
          err_mes = trim(adjustl(fname))//" number is wrong."
          call write_err_stop(err_mes)
          deallocate(err_mes)
        end if
      else if (any(ftype == type_txt(:))) then
        read(unit=fnum,fmt=*,iostat=ierr) fetime
        if (ierr /= 0) then
          call write_err_read(fnum)
        end if
      else if (any(ftype == type_bin(:))) then
        read(unit=fnum,iostat=ierr) fetime
        if (ierr /= 0) then
          call write_err_read(fnum)
        end if
      end if

      fetime = fetime*fmulti

      if (fetime <= st_init%rest_time .and. st_sim%res_type == 1) then
        if (ftype == in_type(1) .or. ftype == in_type(2)) then
          call skip_clas_point(fnum, file_totn)
        else if (ftype == in_type(3)) then
          call skip_2darray_txt(fnum, st_grid%nx, st_grid%ny)
        else if (ftype == in_type(5)) then
          call skip_3darray_txt(fnum, st_grid%nx, st_grid%ny, st_grid%nz)
        else if (ftype == in_type(4)) then
          call skip_2darray_bin(fnum, st_grid%nx, st_grid%ny)
        else if (ftype == in_type(6)) then
          call skip_3darray_bin(fnum, st_grid%nx, st_grid%ny, st_grid%nz)
        end if
      else
        time_flag = 1
      end if
    end do

    deallocate(type_txt, type_bin)

  end subroutine skip_file

  subroutine skip_file_int(ftype, bunit, fnum, fname, fmulti, fstep, finend, fetime)
  !***************************************************************************************
  ! skip_file_int -- Skip for time interval list file
  !***************************************************************************************
    ! -- modules
    use utility_module, only: close_file
    ! -- inout
    integer(I4), intent(in) :: ftype, bunit
    integer(I4), intent(inout) :: fnum
    character(*), intent(in) :: fname
    real(SP), intent(in) :: fmulti, fstep, finend
    real(SP), intent(inout) :: fetime
    ! -- local
    integer(I4) :: ierr
    integer(I4) :: time_flag, count_num
    integer(I4), allocatable :: type_txt(:), type_bin(:)
    character(CHALEN) :: intpath
    character(:), allocatable :: err_mes
    !-------------------------------------------------------------------------------------
    ierr = 0
    allocate(type_txt(2), type_bin(2))
    type_txt(:) = [in_type(3), in_type(5)] ; type_bin(:) = [in_type(4), in_type(6)]

    time_flag = 0 ; count_num = 0
    do while (time_flag == 0)
      fetime = finend + fstep*count_num
      fetime = fetime*fmulti
      if (fetime <= st_init%rest_time .and. st_sim%res_type == 1) then
        err_mes = "input"//trim(adjustl(fname))//" time interval"
        if (any(ftype == type_txt(:))) then
          call close_file(fnum)
          read(unit=bunit,fmt='(a)',iostat=ierr) intpath
          if (ierr /= 0) then
            call write_err_read(bunit)
          end if
          ! -- Open new read text file (new_rtxt)
            call open_new_rtxt(1, 0, trim(adjustl(intpath)), err_mes, fnum)
        else if (any(ftype == type_bin(:))) then
          call close_file(fnum)
          read(unit=bunit,fmt='(a)',iostat=ierr) intpath
          if (ierr /= 0) then
            call write_err_read(bunit)
          end if
          ! -- Open new read binary file (new_rbin)
            call open_new_rbin(1, 0, trim(adjustl(intpath)), err_mes, fnum)
        end if
        count_num = count_num + 1
      else
        time_flag = 1
        if (allocated(err_mes)) then
          deallocate(err_mes)
        end if
      end if
    end do

    deallocate(type_txt, type_bin)

  end subroutine skip_file_int

  subroutine skip_clas_point(fnum, targn)
  !***************************************************************************************
  ! skip_clas_point -- Skip classification and point input file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, targn
    ! -- local
    integer(I4) :: i
    integer(I4) :: ierr
    !-------------------------------------------------------------------------------------
    ierr = 0
    do i = 1, targn
      read(unit=fnum,fmt=*,iostat=ierr)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine skip_clas_point

  subroutine skip_2darray_txt(fnum, inx, iny)
  !***************************************************************************************
  ! skip_2darray_txt -- Skip 2d array text format file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny
    ! -- local
    integer(I4) :: i, j
    integer(I4) :: ierr
    real(DP) :: dummy
    !-------------------------------------------------------------------------------------
    ierr = 0
    do j = 1, iny
      read(unit=fnum,fmt=*,iostat=ierr) (dummy, i = 1, inx)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end do

  end subroutine skip_2darray_txt

  subroutine skip_3darray_txt(fnum, inx, iny, inz)
  !***************************************************************************************
  ! skip_3darray_txt -- Skip 3d array text format file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny, inz
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: ierr
    real(DP) :: dummy
    !-------------------------------------------------------------------------------------
    ierr = 0
    do k = 1, inz
      do j = 1, iny
        read(unit=fnum,fmt=*,iostat=ierr) (dummy, i = 1, inx)
        if (ierr /= 0) then
          call write_err_read(fnum)
        end if
      end do
    end do

  end subroutine skip_3darray_txt

  subroutine skip_2darray_bin(fnum, inx, iny)
  !***************************************************************************************
  ! skip_2darray_bin -- Skip 2d array binary format file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny
    ! -- local
    integer(I4) :: i, j
    integer(I4) :: ierr
    real(SP) :: dummy
    !-------------------------------------------------------------------------------------
    ierr = 0
    read(unit=fnum,iostat=ierr) ((dummy, i = 1, inx), j = 1, iny)
    if (ierr /= 0) then
      call write_err_read(fnum)
    end if

  end subroutine skip_2darray_bin

  subroutine skip_3darray_bin(fnum, inx, iny, inz)
  !***************************************************************************************
  ! skip_3darray_bin -- Skip 3d array binary format file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, inx, iny, inz
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: ierr
    real(SP) :: dummy
    !-------------------------------------------------------------------------------------
    ierr = 0
    read(unit=fnum,iostat=ierr) (((dummy, i = 1, inx), j = 1, iny), k = 1, inz)
    if (ierr /= 0) then
      call write_err_read(fnum)
    end if

  end subroutine skip_3darray_bin

  subroutine flat_2d_2di4(inx, iny, array_in, flat_out)
  !***************************************************************************************
  ! flat_2d_2di4 -- Flatten 2d array to 2d flat integer array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny
    integer(I4), intent(in) :: array_in(:,:)
    integer(I4), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, s_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, s_num)
    do j = 1, iny
      do i = 1, inx
        s_num = inx*(j-1) + i
        flat_out(s_num) = array_in(i,j)
      end do
    end do
    !$omp end parallel do

  end subroutine flat_2d_2di4

  subroutine flat_2d_2dr4(inx, iny, array_in, flat_out)
  !***************************************************************************************
  ! flat_2d_2dr4 -- Flatten 2d array to 2d flat real4 array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny
    real(SP), intent(in) :: array_in(:,:)
    real(SP), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, s_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, s_num)
    do j = 1, iny
      do i = 1, inx
        s_num = inx*(j-1) + i
        flat_out(s_num) = array_in(i,j)
      end do
    end do
    !$omp end parallel do

  end subroutine flat_2d_2dr4

  subroutine flat_2d_2dr8(inx, iny, array_in, flat_out)
  !***************************************************************************************
  ! flat_2d_2dr8 -- Flatten 2d array to 2d flat real8 array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny
    real(DP), intent(in) :: array_in(:,:)
    real(DP), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, s_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, s_num)
    do j = 1, iny
      do i = 1, inx
        s_num = inx*(j-1) + i
        flat_out(s_num) = array_in(i,j)
      end do
    end do
    !$omp end parallel do

  end subroutine flat_2d_2dr8

  subroutine flat_2d_3di4(inx, iny, inz, array_in, flat_out)
  !***************************************************************************************
  ! flat_2d_3di4 --  Flatten 2d array to 3d flat integer array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny, inz
    integer(I4), intent(in) :: array_in(:,:)
    integer(I4), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, k, c_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num)
    do k = 1, inz
      do j = 1, iny
        do i = 1, inx
          c_num = inx*(iny*(k-1) + (j-1)) + i
          flat_out(c_num) = array_in(i,j)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine flat_2d_3di4

  subroutine flat_2d_3dr4(inx, iny, inz, array_in, flat_out)
  !***************************************************************************************
  ! flat_2d_3dr4 --  Flatten 2d array to 3d flat real4 array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny, inz
    real(SP), intent(in) :: array_in(:,:)
    real(SP), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, k, c_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num)
    do k = 1, inz
      do j = 1, iny
        do i = 1, inx
          c_num = inx*(iny*(k-1) + (j-1)) + i
          flat_out(c_num) = array_in(i,j)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine flat_2d_3dr4

  subroutine flat_2d_3dr8(inx, iny, inz, array_in, flat_out)
  !***************************************************************************************
  ! flat_2d_3dr8 --  Flatten 2d array to 3d flat real8 array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny, inz
    real(DP), intent(in) :: array_in(:,:)
    real(DP), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, k, c_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num)
    do k = 1, inz
      do j = 1, iny
        do i = 1, inx
          c_num = inx*(iny*(k-1) + (j-1)) + i
          flat_out(c_num) = array_in(i,j)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine flat_2d_3dr8

  subroutine flat_3d_3di4(inx, iny, inz, array_in, flat_out)
  !***************************************************************************************
  ! flat_3d_3di4 -- Convert 3d array to 3d flat integer array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny, inz
    integer(I4), intent(in) :: array_in(:,:,:)
    integer(I4), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, k, c_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num)
    do k = 1, inz
      do j = 1, iny
        do i = 1, inx
          c_num = inx*(iny*(k-1) + (j-1)) + i
          flat_out(c_num) = array_in(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine flat_3d_3di4

  subroutine flat_3d_3dr4(inx, iny, inz, array_in, flat_out)
  !***************************************************************************************
  ! flat_3d_3dr4 -- Convert 3d array to 3d flat real4 array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny, inz
    real(SP), intent(in) :: array_in(:,:,:)
    real(SP), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, k, c_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num)
    do k = 1, inz
      do j = 1, iny
        do i = 1, inx
          c_num = inx*(iny*(k-1) + (j-1)) + i
          flat_out(c_num) = array_in(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine flat_3d_3dr4

  subroutine flat_3d_3dr8(inx, iny, inz, array_in, flat_out)
  !***************************************************************************************
  ! flat_3d_3dr8 -- Convert 3d array to 3d flat real8 array
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: inx, iny, inz
    real(DP), intent(in) :: array_in(:,:,:)
    real(DP), intent(out) :: flat_out(:)
    ! -- local
    integer(I4) :: i, j, k, c_num
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, c_num)
    do k = 1, inz
      do j = 1, iny
        do i = 1, inx
          c_num = inx*(iny*(k-1) + (j-1)) + i
          flat_out(c_num) = array_in(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine flat_3d_3dr8

  subroutine read_next(ftype, fnum, multi, mess, nx_totn, flag, ierr, etime)
  !***************************************************************************************
  ! read_next -- Read next time
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_read, only: read_mpi_head
#endif
    ! -- inout
    integer(I4), intent(in) :: ftype, fnum
    real(SP), intent(in) :: multi
    character(*), intent(in) :: mess
    integer(I4), intent(out) :: nx_totn, flag, ierr
    real(SP), intent(out) :: etime
    ! -- local
    integer(I4), allocatable :: type_txt(:), type_bin(:)
    real(SP) :: temp_etime
    character(:), allocatable :: err_mes
    !-------------------------------------------------------------------------------------
    allocate(type_txt(2), type_bin(2))
    type_txt(:) = [in_type(3), in_type(5)] ; type_bin(:) = [in_type(4), in_type(6)]

    ierr = 0 ; err_mes = "Read final step in "//mess//" file."
    if (ftype == in_type(1) .or. ftype == in_type(2)) then
      if (my_rank == 0) then
        read(unit=fnum,fmt=*,iostat=ierr) nx_totn, etime
        if (ierr /= 0) then
          call write_logf(err_mes)
          flag = 0 ; etime = st_sim%end_time
        else
          etime = etime*multi
        end if
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast scalar value (val)
          call bcast_val(nx_totn, mess//" total number")
      end if
#endif

    else if (any(ftype == type_txt(:))) then
      if (my_rank == 0) then
        read(unit=fnum,fmt=*,iostat=ierr) etime
        if (ierr /= 0) then
          call write_logf(err_mes)
          flag = 0 ; etime = st_sim%end_time
        else
          etime = etime*multi
        end if
      end if

    else if (any(ftype == type_bin(:))) then
#ifdef MPI_MSG
      ! -- Read mpi header (mpi_head)
        call read_mpi_head(fnum, ierr, temp_etime)
#else
      read(unit=fnum,iostat=ierr) temp_etime
#endif
      if (ierr /= 0) then
        if (my_rank == 0) then
          call write_logf(err_mes)
        end if
        flag = 0 ; etime = st_sim%end_time
      else
        etime = temp_etime*multi
      end if
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast scalar value (val)
        call bcast_val(flag, mess//" flag value")
      ! -- Bcast scalar value (val)
        call bcast_val(ierr, mess//" error value")
      ! -- Bcast scalar value (val)
        call bcast_val(etime, mess//" end time value")
    end if
#endif

    deallocate(type_txt, type_bin)

  end subroutine read_next

  subroutine read_intn(ftype, fnum, multi, fstep, mess, intfn, flag, ierr, etime)
  !***************************************************************************************
  ! read_intn -- Read next time for time interval file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
!    use mpi_utility, only: bcast_file_path
    use mpi_utility, only: bcast_file
    use mpi_read, only: open_mpi_read_file
#endif
    ! -- inout
    integer(I4), intent(in) :: ftype, fnum
    real(SP), intent(in) :: multi, fstep
    character(*), intent(in) :: mess
    integer(I4), intent(inout) :: intfn, flag, ierr
    real(SP), intent(inout) :: etime
    ! -- local
    integer(I4), allocatable :: type_txt(:), type_bin(:)
    character(CHALEN) :: nxi_path
    character(:), allocatable :: err_mes
    !-------------------------------------------------------------------------------------
    allocate(type_txt(2), type_bin(2))
    type_txt(:) = [in_type(3), in_type(5)] ; type_bin(:) = [in_type(4), in_type(6)]

    ierr = 0 ; err_mes = "Read final step in "//mess//" timeseries file."
    if (my_rank == 0) then
      read(unit=fnum,fmt='(a)',iostat=ierr) nxi_path
      if (ierr /= 0) then
        call write_logf(err_mes)
      end if
    end if
#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast scalar value (val)
        call bcast_val(ierr, mess//" error value")
      ! -- Bcast file (file)
        call bcast_file(nxi_path, mess//" timeseries")
    end if
#endif

    if (ierr /= 0) then
      flag = 0 ; etime = st_sim%end_time
      return
    end if

    if (any(ftype == type_txt(:))) then
      if (my_rank == 0) then
        ! -- Open new read text file (new_rtxt)
          call open_new_rtxt(0, 0, trim(adjustl(nxi_path)), err_mes, intfn, ierr)
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Bcast scalar value (val)
          call bcast_val(ierr, mess//" error value")
        ! -- Bcast scalar value (val)
          call bcast_val(intfn, mess//" file number")
      end if
#endif

    else if (any(ftype == type_bin(:))) then
#ifdef MPI_MSG
      ! -- Open mpi read file (mpi_read_file)
        call open_mpi_read_file(0, 0, trim(adjustl(nxi_path)), err_mes, intfn, ierr)
#else
      ! -- Open read new binary file (new_rbin)
        call open_new_rbin(0, 0, trim(adjustl(nxi_path)), err_mes, intfn, ierr)
#endif
    end if

    etime = etime + fstep*multi

    deallocate(type_txt, type_bin)

  end subroutine read_intn

  subroutine read_2d_calcreg(fnum, ftype, tg_flag, c_count)
  !***************************************************************************************
  ! read_2d_calcreg -- Read calcluation region file from 2d file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, ftype
    integer(I4), intent(out) :: tg_flag(:)
    integer(I4), intent(out) :: c_count
    ! -- local
    integer(I4) :: i, j, k, c_num, ierr
    integer(I4), allocatable :: array_read(:,:)
    !-------------------------------------------------------------------------------------
    allocate(array_read(st_grid%nx,st_grid%ny))
    !$omp parallel workshare
    array_read(:,:) = 0
    !$omp end parallel workshare

    ierr = 0
    if (ftype == in_type(3)) then
      do j = 1, st_grid%ny
        read(unit=fnum,fmt=*,iostat=ierr) (array_read(i,j), i = 1, st_grid%nx)
        if (ierr /= 0) then
          call write_err_read(fnum)
        end if
      end do
    else if (ftype == in_type(4)) then
      read(unit=fnum,iostat=ierr) ((array_read(i,j), i = 1, st_grid%nx), j = 1, st_grid%ny)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end if

    c_count = 0
    !$omp parallel do private(i, j, k, c_num) reduction(+:c_count)
    do k = 1, st_grid%nz
      do j = 1, st_grid%ny
        do i = 1, st_grid%nx
          c_num = st_grid%nx*st_grid%ny*(k-1) + st_grid%nx*(j-1) + i
          tg_flag(c_num) = array_read(i,j)
          if (tg_flag(c_num) > 0) then
            c_count = c_count + 1
          end if
        end do
      end do
    end do
    !$omp end parallel do

    deallocate(array_read)

  end subroutine read_2d_calcreg

  subroutine read_3d_calcreg(fnum, ftype, tg_flag, c_count)
  !***************************************************************************************
  ! read_3d_calcreg -- Read calcluation region file from 3d file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, ftype
    integer(I4), intent(out) :: tg_flag(:)
    integer(I4), intent(out) :: c_count
    ! -- local
    integer(I4) :: i, j, k, c_num, ierr
    integer(I4), allocatable :: array_read(:,:,:)
    !-------------------------------------------------------------------------------------
    allocate(array_read(st_grid%nx,st_grid%ny,st_grid%nz))
    !$omp parallel workshare
    array_read(:,:,:) = 0
    !$omp end parallel workshare

    ierr = 0
    if (ftype == in_type(5)) then
      do k = 1, st_grid%nz
        do j = 1, st_grid%ny
          read(unit=fnum,fmt=*,iostat=ierr) (array_read(i,j,k), i = 1, st_grid%nx)
          if (ierr /= 0) then
            call write_err_read(fnum)
          end if
        end do
      end do
    else if (ftype == in_type(6)) then
      read(unit=fnum,iostat=ierr) (((array_read(i,j,k), i = 1, st_grid%nx), j = 1, st_grid%ny), k = 1, st_grid%nz)
      if (ierr /= 0) then
        call write_err_read(fnum)
      end if
    end if

    c_count = 0
    !$omp parallel do private(i, j, k, c_num) reduction(+:c_count)
    do k = 1, st_grid%nz
      do j = 1, st_grid%ny
        do i = 1, st_grid%nx
          c_num = st_grid%nx*st_grid%ny*(k-1) + st_grid%ny*(j-1) + i
          tg_flag(c_num) = array_read(i,j,k)
          if (tg_flag(c_num) > 0) then
            c_count = c_count + 1
          end if
        end do
      end do
    end do
    !$omp end parallel do

    deallocate(array_read)

  end subroutine read_3d_calcreg

end module read_module
