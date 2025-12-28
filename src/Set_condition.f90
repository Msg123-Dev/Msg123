module set_condition
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SZERO, DZERO, DONE, SNOVAL
  use initial_module, only: my_rank, st_sim, st_grid, st_well, in_type
  use read_module, only: read_2dtxt, read_2dbin, read_3dtxt, read_3dbin, flat_2dto2d,&
                         flat_2dto3d, flat_3dto3d
  use set_cell, only: ncalc, ncals, loc2glo_ijk, loc2glo_ij, neib_ncalc
  use make_cell, only: cell_top, cell_bot
#ifdef MPI_MSG
  use mpi_utility, only: mpimax_val
  use utility_module, only: conv_i2s
  use mpi_read, only: read_mpi_file
  use mpi_set, only: scatter_xyval, scatter_xyzval
#endif

  implicit none
  private
  public :: set_clas2calc, set_point2seal, set_point2surf, set_bound2calc, set_mass2calc
  public :: set_2dfile2seal, set_3dfile2seal
  public :: set_2dfile2cals, set_2dfile2calc, set_3dfile2calc
  public :: set_point2well, set_2dwell, set_3dfile2well
  public :: set_well2index, set_well3d2index, set_wellprop
  public :: set_connect, set_srabyd, set_chabyd, set_wellconn

  interface set_2dfile2cals
    module procedure set_2di4_cals
    module procedure set_2dr4_cals
    module procedure set_2dr8_cals
  end interface

  interface set_2dfile2calc
    module procedure set_2di4_calc
    module procedure set_2dr4_calc
    module procedure set_2dr8_calc
  end interface

  interface set_3dfile2calc
    module procedure set_3di4_calc
    module procedure set_3dr4_calc
    module procedure set_3dr8_calc
  end interface

  integer(I4), public :: nseal, tconn_num
  integer(I4), allocatable, public :: off_row(:), off_index(:), conn_dir(:), sea_dir(:)
  real(DP), allocatable, public :: area_dis(:), sat_hydf(:), conn_dis(:)
  real(DP), allocatable, public :: sea_hydf(:), sea_dis(:)
  real(DP), allocatable, public :: surf_area(:), rech_area(:)
  real(DP), allocatable, public :: hydf_surf(:), abyd_surf(:), abyd_well(:), sea_abyd(:)
  integer(I4), allocatable, public :: well_index(:), well_conn(:)
  integer(I4), allocatable, public :: sea2cal(:), sea2sea(:)
  integer(I4), allocatable, public :: left_off(:), right_off(:)

  ! -- local
  integer(I4), allocatable :: well_nflag(:)
  real(DP), allocatable :: reg_dis(:,:)
  real(DP), allocatable :: surf_dis(:)

  contains

  subroutine set_clas2calc(tgn, tg_name, tg_val, calc_val, tg_flag, tg_num)
  !***************************************************************************************
  ! set_clas2calc -- Set calculation value from classification
  !***************************************************************************************
    ! -- modules
    use initial_module, only: st_clas
    use set_cell, only: clas_flag
    ! -- inout
    integer(I4), intent(in) :: tgn
    character(*), intent(in) :: tg_name(:)
    real(SP), intent(in) :: tg_val(:)
    real(SP), intent(inout) :: calc_val(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: out_num
    integer(I4), allocatable :: temp_flag(:)
    !-------------------------------------------------------------------------------------
    out_num = size(calc_val(:))
    allocate(temp_flag(out_num))
    !$omp parallel workshare
    temp_flag(:) = 0
    !$omp end parallel workshare

    do i = 1, tgn
      if (tg_val(i) /= SNOVAL) then
        do k = 1, st_clas%totn
          if (tg_name(i) == st_clas%name(k)) then
            !$omp parallel do private(j)
            do j = 1, out_num
              if (clas_flag(j,k) == 1) then
                calc_val(j) = tg_val(i)
                temp_flag(j) = 1
              end if
            end do
            !$omp end parallel do
          end if
        end do
      end if
    end do

    if (present(tg_flag)) then
      !$omp parallel workshare
      tg_flag(:) = temp_flag(:)
      !$omp end parallel workshare
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

    deallocate(temp_flag)

  end subroutine set_clas2calc

  subroutine set_point2seal(tgn, p_i, p_j, p_k, tg_val, cell_val, tg_flag, tg_num)
  !***************************************************************************************
  ! set_point2seal -- Set seal value from point
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: tgn
    integer(I4), intent(in) :: p_i(:), p_j(:), p_k(:)
    real(SP), intent(in) :: tg_val(:)
    real(SP), intent(out) :: cell_val(:)
    integer(I4), intent(out) :: tg_flag(:)
    integer(I4), intent(out) :: tg_num
    ! -- local
    integer(I4) :: i, c_num, ii, jj, kk, mpi_ncalc, cnum
    integer(I4) :: gridx, gridy, gridz
    integer(I4) :: pi, pj, pk, ist, ien, jst, jen, kst, ken, l_ijk
    !-------------------------------------------------------------------------------------
    !$omp parallel workshare
    tg_flag(:) = 0
    !$omp end parallel workshare
    mpi_ncalc = ncalc + neib_ncalc ; cnum = size(cell_val(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz
    !$omp parallel do private(i, ii, jj, kk, pi, pj, pk, ist, ien, jst, jen, kst, ken, c_num, l_ijk)
    do i = 1, tgn
      if (tg_val(i) /= SNOVAL) then
        pi = p_i(i) ; pj = p_j(i) ; pk = p_k(i)
        if (pi == -1 .and. pj == -1 .and. pk == -1) then !all cells
          ist = 1 ; ien = gridx ; jst = 1 ; jen = gridy ; kst = 1 ; ken = gridz
        else if (pi == -1 .and. pj == -1) then !i,j cells
          ist = 1 ; ien = gridx ; jst = 1 ; jen = gridy ; kst = pk ; ken = pk
        else if (pi == -1 .and. pk == -1) then !i,k cells
          ist = 1 ; ien = gridx ; jst = pj ; jen = pj ; kst = 1 ; ken = gridz
        else if (pj == -1 .and. pk == -1) then !j,k cells
          ist = pi ; ien = pi ; jst = 1 ; jen = gridy ; kst = 1 ; ken = gridz
        else if (pi == -1) then !only i cell
          ist = 1 ; ien = gridx ; jst = pj ; jen = pj ; kst = pk ; ken = pk
        else if (pj == -1) then !only j cell
          ist = pi ; ien = pi ; jst = 1 ; jen = gridy ; kst = pk ; ken = pk
        else if (pk == -1) then !only k cell
          ist = pi ; ien = pi ; jst = pj ; jen = pj ; kst = 1 ; ken = gridz
        else !others
          ist = pi ; ien = pi ; jst = pj ; jen = pj ; kst = pk ; ken = pk
        end if

        do ii = ist, ien
          do jj = jst, jen
            do kk = kst, ken
              c_num = gridx*(gridy*(kk-1) + (jj-1)) + ii
              l_ijk = findloc(loc2glo_ijk(:), value = c_num, dim = 1)
!              l_ijk = glo2loc_ijk(c_num)
              if (l_ijk > mpi_ncalc .and. l_ijk <= cnum) then
                cell_val(l_ijk-mpi_ncalc) = tg_val(i)
                tg_flag(l_ijk-mpi_ncalc) = 1
              end if
            end do
          end do
        end do
      end if
    end do
    !$omp end parallel do

    call count_flag(tg_flag, tg_num)

  end subroutine set_point2seal

  subroutine set_point2surf(tgn, p_i, p_j, tg_val, cell_val, tg_flag, tg_num)
  !***************************************************************************************
  ! set_point2surf -- Set surface value from point
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: tgn
    integer(I4), intent(in) :: p_i(:), p_j(:)
    real(SP), intent(in) :: tg_val(:)
    real(SP), intent(out) :: cell_val(:)
    integer(I4), intent(out) :: tg_flag(:)
    integer(I4), intent(out) :: tg_num
    ! -- local
    integer(I4) :: i, s_num, ii, jj
    integer(I4) :: gridx, gridy
    integer(I4) :: pi, pj, ist, ien, jst, jen, l_ij
    !-------------------------------------------------------------------------------------
    !$omp parallel workshare
    tg_flag(:) = 0
    !$omp end parallel workshare
    gridx = st_grid%nx ; gridy = st_grid%ny
    !$omp parallel do private(i, ii, jj, pi, pj, ist, ien, jst, jen, s_num, l_ij)
    do i = 1, tgn
      pi = p_i(i) ; pj = p_j(i)
      if (pi == -1 .and. pj == -1) then !all i,j cells
        ist = 1 ; ien = gridx ; jst = 1 ; jen = gridy
      else if (pi == -1) then !only i cell
        ist = 1 ; ien = gridx ; jst = pj ; jen = pj
      else if (pj == -1) then !only j cell
        ist = pi ; ien = pi ; jst = 1 ; jen = gridy
      else !others
        ist = pi ; ien = pi ; jst = pj ; jen = pj
      end if

      do ii = ist, ien
        do jj = jst, jen
          s_num = gridx*(jj-1) + ii
          l_ij = 0
          l_ij = findloc(loc2glo_ij(:), value = s_num, dim = 1)
          if (l_ij > 0 .and. l_ij <= ncals) then
            cell_val(l_ij) = tg_val(i)
            tg_flag(l_ij) = 1
          end if
        end do
      end do

    end do
    !$omp end parallel do

    call count_flag(tg_flag, tg_num)

  end subroutine set_point2surf

  subroutine set_bound2calc(bound_num, c_flag, cell_val, bound2calc, b_val)
  !***************************************************************************************
  ! set_bound2calc -- Set boundary value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: bound_num
    integer(I4), intent(in) :: c_flag(:)
    real(SP), intent(in) :: cell_val(:)
    integer(I4), intent(out) :: bound2calc(:)
    real(DP), intent(out) :: b_val(:)
    ! -- local
    integer(I4) :: i, bc
    !-------------------------------------------------------------------------------------
    bc = 0
    do i = 1, bound_num
      if (c_flag(i) == 1) then
        bc = bc + 1
        bound2calc(bc) = i
        b_val(bc) = cell_val(i)
      end if
    end do

  end subroutine set_bound2calc

  subroutine set_mass2calc(mass_num, c_flag, cell_val, mass2calc, mass_val)
  !***************************************************************************************
  ! set_mass2calc -- Set massbalance value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: mass_num
    integer(I4), intent(in) :: c_flag(:), cell_val(:)
    integer(I4), intent(out) :: mass2calc(:), mass_val(:)
    ! -- local
    integer(I4) :: i, bc
    !-------------------------------------------------------------------------------------
    bc = 0
    do i = 1, mass_num
      if (c_flag(i) == 1) then
        bc = bc + 1
        mass2calc(bc) = i
        mass_val(bc) = cell_val(i)
      end if
    end do

  end subroutine set_mass2calc

  subroutine set_2dfile2seal(fnum, ftype, int_ft, rsnum, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_2dfile2seal -- Set sea level value from 2d file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use set_cell, only: loc2glo_nos
#endif
    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft, rsnum
    real(SP), intent(in) :: no_val
    real(SP), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
#ifdef MPI_MSG
    integer(I4) :: i, j, k, s
    real(SP), allocatable :: array_seas(:), array_temp(:)
    character(:), allocatable :: str_fnum
#else
    integer(I4) :: dummy
#endif
    real(SP), allocatable :: array_flat(:)
    real(SP), allocatable :: array_read(:,:)
    !-------------------------------------------------------------------------------------
    cnum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    if (ftype == in_type(3) .or. int_ft == in_type(3)) then
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy), array_flat(gridxyz))
        !$omp parallel workshare
        array_read(:,:) = no_val ; array_flat(:) = no_val
        !$omp end parallel workshare
        call read_2dtxt(fnum, gridx, gridy, array_read)
        call flat_2dto3d(gridx, gridy, gridz, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(4) .or. int_ft == in_type(4)) then
#ifdef MPI_MSG
      allocate(array_seas(rsnum), array_temp(gridxyz), array_flat(gridxyz))
      !$omp parallel workshare
      array_seas(:) = no_val ; array_temp(:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_mpi_file(ftype, int_ft, fnum, rsnum, array_seas)
      !$omp parallel do private(i, j, k, s)
      do i = 1, rsnum
        s = 0 ; s = loc2glo_nos(i)
        if (s /= 0) then
          do j = 1, gridz
            k = gridx*gridy*(j-1) + s
            array_temp(k) = array_seas(i)
          end do
        end if
      end do
      !$omp end parallel do
      str_fnum = conv_i2s(fnum)
      ! -- MAX value for MPI (val)
        call mpimax_val(array_temp, "file number "//str_fnum, array_flat)
      deallocate(array_seas, array_temp)
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#else
      dummy = rsnum
      allocate(array_read(gridx,gridy), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_2dbin(fnum, gridx, gridy, no_val, array_read)
      call flat_2dto3d(gridx, gridy, gridz, array_read, array_flat)
      deallocate(array_read)
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_real4(cnum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_2dfile2seal

  subroutine set_3dfile2seal(fnum, ftype, int_ft, rcnum, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_3dfile2seal -- Set sea level from 3d file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use set_cell, only: loc2glo_noc
#endif
    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft, rcnum
    real(SP), intent(in) :: no_val
    real(SP), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
#ifdef MPI_MSG
    integer(I4) :: i, c
    real(SP), allocatable :: array_seac(:), array_temp(:)
    character(:), allocatable :: str_fnum
#else
    integer(I4) :: dummy
#endif
    real(SP), allocatable :: array_flat(:)
    real(SP), allocatable :: array_read(:,:,:)
    !-------------------------------------------------------------------------------------
    cnum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    if (ftype == in_type(5) .or. int_ft == in_type(5)) then
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy,gridz), array_flat(gridxyz))
        !$omp parallel workshare
        array_read(:,:,:) = no_val ; array_flat(:) = no_val
        !$omp end parallel workshare
        call read_3dtxt(fnum, gridx, gridy, gridz, array_read)
        call flat_3dto3d(gridx, gridy, gridz, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(6) .or. int_ft == in_type(6)) then
#ifdef MPI_MSG
      allocate(array_seac(rcnum), array_temp(gridxyz), array_flat(gridxyz))
      !$omp parallel workshare
      array_seac(:) = no_val ; array_temp(:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_mpi_file(ftype, int_ft, fnum, rcnum, array_seac)
      !$omp parallel do private(i, c)
      do i = 1, rcnum
        c = 0 ; c = loc2glo_noc(i)
        if (c /= 0) then
          array_temp(c) = array_seac(i)
        end if
      end do
      !$omp end parallel do
      str_fnum = conv_i2s(fnum)
      ! -- MAX value for MPI (val)
        call mpimax_val(array_temp, "file number "//str_fnum, array_flat)
      deallocate(array_seac, array_temp)
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#else
      dummy = rcnum
      allocate(array_read(gridx,gridy,gridz), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_3dbin(fnum, gridx, gridy, gridz, no_val, array_read)
      call flat_3dto3d(gridx, gridy, gridz, array_read, array_flat)
      deallocate(array_read)
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_real4(cnum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_3dfile2seal

  subroutine set_2di4_cals(fnum, ftype, int_ft, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_2di4_cals -- Set surface calculation value from 2d array integer file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_set, only: scatter_xyval
#endif
    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft, no_val
    integer(I4), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: snum
    integer(I4) :: gridx, gridy
    integer(I4), allocatable :: array_flat(:)
    integer(I4), allocatable :: array_read(:,:)
    !-------------------------------------------------------------------------------------
    snum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny
    if (ftype == in_type(3) .or. int_ft == in_type(3)) then
      allocate(array_flat(gridx*gridy))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy))
        !$omp parallel workshare
        array_read(:,:) = no_val
        !$omp end parallel workshare
        call read_2dtxt(fnum, gridx, gridy, array_read)
        !$omp parallel workshare
        array_flat(:) = reshape(array_read(:,:), shape(array_flat))
        !$omp end parallel workshare
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xy value (xyval)
        call scatter_xyval(snum, loc2glo_ij, array_flat, val_out)
#else
      call set_calci4(snum, loc2glo_ij, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(4) .or. int_ft == in_type(4)) then
#ifdef MPI_MSG
      allocate(array_flat(snum))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      call read_mpi_file(ftype, int_ft, fnum, snum, array_flat)
      !$omp parallel workshare
      val_out(:) = array_flat(:)
      !$omp end parallel workshare
#else
      allocate(array_read(gridx,gridy), array_flat(gridx*gridy))
      !$omp parallel workshare
      array_read(:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_2dbin(fnum, gridx, gridy, no_val, array_read)
      !$omp parallel workshare
      array_flat(:) = reshape(array_read(:,:), shape(array_flat))
      !$omp end parallel workshare
      deallocate(array_read)
      call set_calci4(snum, loc2glo_ij, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_int(snum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_2di4_cals

  subroutine set_2dr4_cals(fnum, ftype, int_ft, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_2dr4_cals -- Set surface calculation value from 2d real4 array file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_set, only: scatter_xyval
#endif
    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft
    real(SP), intent(in) :: no_val
    real(SP), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: snum
    integer(I4) :: gridx, gridy
    real(SP), allocatable :: array_flat(:)
    real(SP), allocatable :: array_read(:,:)
    !-------------------------------------------------------------------------------------
    snum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny
    if (ftype == in_type(3) .or. int_ft == in_type(3)) then
      allocate(array_flat(gridx*gridy))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy))
        !$omp parallel workshare
        array_read(:,:) = no_val
        !$omp end parallel workshare
        call read_2dtxt(fnum, gridx, gridy, array_read)
        call flat_2dto2d(gridx, gridy, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xy value (xyval)
        call scatter_xyval(snum, loc2glo_ij, array_flat, val_out)
#else
      call set_calcr4(snum, loc2glo_ij, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(4) .or. int_ft == in_type(4)) then
#ifdef MPI_MSG
      allocate(array_flat(snum))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      call read_mpi_file(ftype, int_ft, fnum, snum, array_flat)
      !$omp parallel workshare
      val_out(:) = array_flat(:)
      !$omp end parallel workshare
#else
      allocate(array_read(gridx,gridy), array_flat(gridx*gridy))
      !$omp parallel workshare
      array_read(:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_2dbin(fnum, gridx, gridy, no_val, array_read)
      call flat_2dto2d(gridx, gridy, array_read, array_flat)
      deallocate(array_read)
      call set_calcr4(snum, loc2glo_ij, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_real4(snum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_2dr4_cals

  subroutine set_2dr8_cals(fnum, ftype, int_ft, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_2dr8_cals -- Set surface calculation value from 2d real8 array file
  !***************************************************************************************
    ! -- modules
#ifdef MPI_MSG
    use mpi_set, only: scatter_xyval
#endif
    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft
    real(DP), intent(in) :: no_val
    real(DP), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: snum
    integer(I4) :: gridx, gridy
    real(DP), allocatable :: array_flat(:)
    real(DP), allocatable :: array_read(:,:)
    !-------------------------------------------------------------------------------------
    snum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny
    if (ftype == in_type(3) .or. int_ft == in_type(3)) then
      allocate(array_flat(gridx*gridy))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy))
        !$omp parallel workshare
        array_read(:,:) = no_val
        !$omp end parallel workshare
        call read_2dtxt(fnum, gridx, gridy, array_read)
        call flat_2dto2d(gridx, gridy, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter real xy value (rxyval)
        call scatter_xyval(snum, loc2glo_ij, array_flat, val_out)
#else
      call set_calcr8(snum, loc2glo_ij, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(4) .or. int_ft == in_type(4)) then
#ifdef MPI_MSG
      allocate(array_flat(snum))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      call read_mpi_file(ftype, int_ft, fnum, snum, array_flat)
      !$omp parallel workshare
      val_out(:) = array_flat(:)
      !$omp end parallel workshare
#else
      allocate(array_read(gridx,gridy), array_flat(gridx*gridy))
      !$omp parallel workshare
      array_read(:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_2dbin(fnum, gridx, gridy, no_val, array_read)
      call flat_2dto2d(gridx, gridy, array_read, array_flat)
      deallocate(array_read)
      call set_calcr8(snum, loc2glo_ij, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_real8(snum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_2dr8_cals

  subroutine set_2di4_calc(fnum, ftype, int_ft, snum, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_2di4_calc -- Set calculation value from 2d integer file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft, snum, no_val
    integer(I4), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: j, k, cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
    integer(I4), allocatable :: array_flat(:)
    integer(I4), allocatable :: array_read(:,:)
#ifdef MPI_MSG
    integer(I4) :: i, s
    integer(I4), allocatable :: array_surf(:), array_temp(:)
    character(:), allocatable :: str_fnum
#else
    integer(I4) :: dummy
#endif
    !-------------------------------------------------------------------------------------
    cnum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    if (ftype == in_type(3) .or. int_ft == in_type(3)) then
      allocate(array_flat(gridxyz))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy))
        !$omp parallel workshare
        array_read(:,:) = no_val
        !$omp end parallel workshare
        call read_2dtxt(fnum, gridx, gridy, array_read)
        !$omp parallel do private(j, k)
        do k = 1, gridz
          j = gridx*gridy*(k-1) + 1
          array_flat(j:j+gridx*gridy-1) = reshape(array_read(:,:), [gridx*gridy])
        end do
        !$omp end parallel do
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calci4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(4) .or. int_ft == in_type(4)) then
#ifdef MPI_MSG
      allocate(array_surf(snum), array_temp(gridxyz), array_flat(gridxyz))
      !$omp parallel workshare
      array_surf(:) = no_val ; array_temp(:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_mpi_file(ftype, int_ft, fnum, snum, array_surf)
      !$omp parallel do private(i, j, k, s)
      do i = 1, snum
        s = 0 ; s = loc2glo_ij(i)
        if (s /= 0) then
          do k = 1, gridz
            j = gridx*gridy*(k-1) + s
            array_temp(j) = array_surf(i)
          end do
        end if
      end do
      !$omp end parallel do
      str_fnum = conv_i2s(fnum)
      ! -- MAX value for MPI (val)
        call mpimax_val(array_temp, "file number "//str_fnum, array_flat)
      deallocate(array_surf, array_temp)
      call set_calci4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#else
      dummy = snum
      allocate(array_read(gridx,gridy), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_2dbin(fnum, gridx, gridy, no_val, array_read)
      !$omp parallel do private(j, k)
      do k = 1, gridz
        j = gridx*gridy*(k-1) + 1
        array_flat(j:j+gridx*gridy-1) = reshape(array_read(:,:), [gridx*gridy])
      end do
      !$omp end parallel do
      deallocate(array_read)
      call set_calci4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_int(cnum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_2di4_calc

  subroutine set_2dr4_calc(fnum, ftype, int_ft, snum, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_2dr4_calc -- Set calculation value from 2d real4 file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft, snum
    real(SP), intent(in) :: no_val
    real(SP), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
    real(SP), allocatable :: array_flat(:)
    real(SP), allocatable :: array_read(:,:)
#ifdef MPI_MSG
    integer(I4) :: i, j, k, s
    real(SP), allocatable :: array_surf(:), array_temp(:)
    character(:), allocatable :: str_fnum
#else
    integer(I4) :: dummy
#endif
    !-------------------------------------------------------------------------------------
    cnum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    if (ftype == in_type(3) .or. int_ft == in_type(3)) then
      allocate(array_flat(gridxyz))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy))
        !$omp parallel workshare
        array_read(:,:) = no_val
        !$omp end parallel workshare
        call read_2dtxt(fnum, gridx, gridy, array_read)
        call flat_2dto3d(gridx, gridy, gridz, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(4) .or. int_ft == in_type(4)) then
#ifdef MPI_MSG
      allocate(array_surf(snum), array_temp(gridxyz), array_flat(gridxyz))
      !$omp parallel workshare
      array_surf(:) = no_val ; array_temp(:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_mpi_file(ftype, int_ft, fnum, snum, array_surf)
      !$omp parallel do private(i, j, k, s)
      do i = 1, snum
        s = 0 ; s = loc2glo_ij(i)
        if (s /= 0) then
          do j = 1, gridz
            k = gridx*gridy*(j-1) + s
            array_temp(k) = array_surf(i)
          end do
        end if
      end do
      !$omp end parallel do
      str_fnum = conv_i2s(fnum)
      ! -- MAX value for MPI (val)
        call mpimax_val(array_temp, "file number "//str_fnum, array_flat)
      deallocate(array_surf, array_temp)
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#else
      dummy = snum
      allocate(array_read(gridx,gridy), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_2dbin(fnum, gridx, gridy, no_val, array_read)
      call flat_2dto3d(gridx, gridy, gridz, array_read, array_flat)
      deallocate(array_read)
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_real4(cnum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_2dr4_calc

  subroutine set_2dr8_calc(fnum, ftype, int_ft, snum, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_2dr8_calc -- Set calculation value from 2d real8 file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft, snum
    real(DP), intent(in) :: no_val
    real(DP), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
    real(DP), allocatable :: array_flat(:)
    real(DP), allocatable :: array_read(:,:)
#ifdef MPI_MSG
    integer(I4) :: i, j, k, s
    real(DP), allocatable :: array_surf(:), array_temp(:)
    character(:), allocatable :: str_fnum
#else
    integer(I4) :: dummy
#endif
    !-------------------------------------------------------------------------------------
    cnum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    if (ftype == in_type(3) .or. int_ft == in_type(3)) then
      allocate(array_flat(gridxyz))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy))
        !$omp parallel workshare
        array_read(:,:) = no_val
        !$omp end parallel workshare
        call read_2dtxt(fnum, gridx, gridy, array_read)
        call flat_2dto3d(gridx, gridy, gridz, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calcr8(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(4) .or. int_ft == in_type(4)) then
#ifdef MPI_MSG
      allocate(array_surf(snum), array_temp(gridxyz), array_flat(gridxyz))
      !$omp parallel workshare
      array_surf(:) = no_val ; array_temp(:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_mpi_file(ftype, int_ft, fnum, snum, array_surf)
      !$omp parallel do private(i, j, k, s)
      do i = 1, snum
        s = 0 ; s = loc2glo_ij(i)
        if (s /= 0) then
          do j = 1, gridz
            k = gridx*gridy*(j-1) + s
            array_temp(k) = array_surf(i)
          end do
        end if
      end do
      !$omp end parallel do
      str_fnum = conv_i2s(fnum)
      ! -- MAX value for MPI (val)
        call mpimax_val(array_temp, "file number "//str_fnum, array_flat)
      deallocate(array_surf, array_temp)
      call set_calcr8(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#else
      dummy = snum
      allocate(array_read(gridx,gridy), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_2dbin(fnum, gridx, gridy, no_val, array_read)
      call flat_2dto3d(gridx, gridy, gridz, array_read, array_flat)
      deallocate(array_read)
      call set_calcr8(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_real8(cnum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_2dr8_calc

  subroutine set_3di4_calc(fnum, ftype, int_ft, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_3di4_calc -- Set calculation value from 3d integer file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft, no_val
    integer(I4), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
    integer(I4), allocatable :: array_flat(:)
    integer(I4), allocatable :: array_read(:,:,:)
    !-------------------------------------------------------------------------------------
    cnum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    if (ftype == in_type(5) .or. int_ft == in_type(5)) then
      allocate(array_flat(gridxyz))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy,gridz))
        !$omp parallel workshare
        array_read(:,:,:) = no_val
        !$omp end parallel workshare
        call read_3dtxt(fnum, gridx, gridy, gridz, array_read)
        !$omp parallel workshare
        array_flat(:) = reshape(array_read(:,:,:), shape(array_flat))
        !$omp end parallel workshare
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calci4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(6) .or. int_ft == in_type(6)) then
#ifdef MPI_MSG
      call read_mpi_file(ftype, int_ft, fnum, cnum, val_out)
#else
      allocate(array_read(gridx,gridy,gridz), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_3dbin(fnum, gridx, gridy, gridz, no_val, array_read)
      !$omp parallel workshare
      array_flat(:) = reshape(array_read(:,:,:), shape(array_flat))
      !$omp end parallel workshare
      deallocate(array_read)
      call set_calci4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_int(cnum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_3di4_calc

  subroutine set_3dr4_calc(fnum, ftype, int_ft, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_3dr4_calc -- Set calculation value from 3d real4 file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft
    real(SP), intent(in) :: no_val
    real(SP), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
    real(SP), allocatable :: array_flat(:)
    real(SP), allocatable :: array_read(:,:,:)
    !-------------------------------------------------------------------------------------
    cnum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    if (ftype == in_type(5) .or. int_ft == in_type(5)) then
      allocate(array_flat(gridxyz))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy,gridz))
        !$omp parallel workshare
        array_read(:,:,:) = no_val
        !$omp end parallel workshare
        call read_3dtxt(fnum, gridx, gridy, gridz, array_read)
        call flat_3dto3d(gridx, gridy, gridz, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(6) .or. int_ft == in_type(6)) then
#ifdef MPI_MSG
      call read_mpi_file(ftype, int_ft, fnum, cnum, val_out)
#else
      allocate(array_read(gridx,gridy,gridz), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_3dbin(fnum, gridx, gridy, gridz, no_val, array_read)
      call flat_3dto3d(gridx, gridy, gridz, array_read, array_flat)
      deallocate(array_read)
      call set_calcr4(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_real4(cnum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_3dr4_calc

  subroutine set_3dr8_calc(fnum, ftype, int_ft, no_val, val_out, tg_flag, tg_num)
  !***************************************************************************************
  ! set_3dr8_calc -- Set calculation value from 3d real8 file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: fnum, ftype, int_ft
    real(DP), intent(in) :: no_val
    real(DP), intent(inout) :: val_out(:)
    integer(I4), intent(out), optional :: tg_flag(:)
    integer(I4), intent(out), optional :: tg_num
    ! -- local
    integer(I4) :: cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
    real(DP), allocatable :: array_flat(:)
    real(DP), allocatable :: array_read(:,:,:)
    !-------------------------------------------------------------------------------------
    cnum = size(val_out(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    if (ftype == in_type(5) .or. int_ft == in_type(5)) then
      allocate(array_flat(gridxyz))
      !$omp parallel workshare
      array_flat(:) = no_val
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy,gridz))
        !$omp parallel workshare
        array_read(:,:,:) = no_val
        !$omp end parallel workshare
        call read_3dtxt(fnum, gridx, gridy, gridz, array_read)
        call flat_3dto3d(gridx, gridy, gridz, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calcr8(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif

    else if (ftype == in_type(6) .or. int_ft == in_type(6)) then
#ifdef MPI_MSG
      call read_mpi_file(ftype, int_ft, fnum, cnum, val_out)
#else
      allocate(array_read(gridx,gridy,gridz), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:,:) = no_val ; array_flat(:) = no_val
      !$omp end parallel workshare
      call read_3dbin(fnum, gridx, gridy, gridz, no_val, array_read)
      call flat_3dto3d(gridx, gridy, gridz, array_read, array_flat)
      deallocate(array_read)
      call set_calcr8(cnum, loc2glo_ijk, no_val, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    if (present(tg_flag)) then
      call set_calc_flag_real8(cnum, no_val, val_out, tg_flag)
    end if

    if (present(tg_num)) then
      call count_flag(tg_flag, tg_num)
    end if

  end subroutine set_3dr8_calc

  subroutine set_point2well(mw_fnum, mw_totn)
  !***************************************************************************************
  ! set_point2well -- Set well from point file
  !***************************************************************************************
    ! -- modules
    use read_module, only: read_wpointf
#ifdef MPI_MSG
    use mpi_set, only: bcast_wellpoint
#endif
    ! -- inout
    integer(I4), intent(in) :: mw_fnum
    integer(I4), intent(out) :: mw_totn
    ! -- local
    integer(I4) :: i, c_num, s_num, rwn, locwn, cwn, oldwn, l_num
    integer(I4), allocatable :: rw_id(:), rw_i(:), rw_j(:), rw_ks(:), rw_ke(:)
    integer(I4), allocatable :: w_id(:), w_ij(:), w_ks(:), w_ke(:)
    real(SP), allocatable :: rw_val(:), w_val(:)
    !-------------------------------------------------------------------------------------
    rwn = st_well%totn
    allocate(rw_id(rwn), rw_i(rwn), rw_j(rwn), rw_ks(rwn), rw_ke(rwn))
    allocate(rw_val(rwn))
    !$omp parallel workshare
    rw_id(:) = 0 ; rw_i(:) = 0 ; rw_j(:) = 0 ; rw_ks(:) = 0 ; rw_ke(:) = 0
    rw_val(:) = SZERO
    !$omp end parallel workshare

    if (my_rank == 0) then
      call read_wpointf(mw_fnum, rwn, rw_id, rw_i, rw_j, rw_ks, rw_ke, rw_val)
    end if

#ifdef MPI_MSG
    call bcast_wellpoint(rwn, rw_id, rw_i, rw_j, rw_ks, rw_ke, rw_val)
#endif

    locwn = rwn
    allocate(w_id(locwn), w_ij(locwn), w_ks(locwn), w_ke(locwn))
    allocate(w_val(locwn))
    !$omp parallel workshare
    w_id(:) = 0 ; w_ij(:) = 0 ; w_ks(:) = 0 ; w_ke(:) = 0
    w_val(:) = SZERO
    !$omp end parallel workshare

    cwn = 0
    do i = 1, locwn
      c_num = st_grid%nx*(st_grid%ny*(rw_ks(i)-1) + (rw_j(i)-1)) + rw_i(i)
      l_num = 0 ; l_num = findloc(loc2glo_ijk(:), value = c_num, dim = 1)
      if (l_num /= 0 .and. l_num <= ncalc) then
        cwn = cwn + 1 ; s_num = st_grid%nx*(rw_j(i)-1) + rw_i(i)
        w_id(cwn) = rw_id(i) ; w_ij(cwn) = s_num
        w_ks(cwn) = rw_ks(i) ; w_ke(cwn) = rw_ke(i) ; w_val(cwn) = rw_val(i)
      end if
    end do

    deallocate(rw_id, rw_i, rw_j, rw_ks, rw_ke, rw_val)

    if (allocated(well_nflag)) then
      if (size(well_nflag) == 0) then
        oldwn = 0
      else
        oldwn = size(well_nflag)
      end if
      mw_totn = oldwn + cwn
      call set_multiwell(mw_totn, cwn, w_id, w_ij, w_ks, w_ke, w_val)
    else
      allocate(well_nflag(cwn), st_well%ij(cwn), st_well%ks(cwn), st_well%ke(cwn))
      allocate(st_well%value(cwn))
      !$omp parallel workshare
      well_nflag(:) = w_id(1:cwn) ; st_well%ij(:) = w_ij(1:cwn)
      st_well%ks(:) = w_ks(1:cwn) ; st_well%ke(:) = w_ke(1:cwn)
      st_well%value(:) = w_val(1:cwn)
      !$omp end parallel workshare
    end if
    mw_totn = size(well_nflag)
    deallocate(w_id, w_ij, w_ks, w_ke, w_val)

  end subroutine set_point2well

  subroutine set_2dwell(wfnum, wf_ftype, wf_int_ft, ws_ftype, we_ftype, mw_totn)
  !***************************************************************************************
  ! set_2dwell -- Set well from 2d file
  !***************************************************************************************
    ! -- modules
    use open_file, only: wells_fnum, welle_fnum
    ! -- inout
    integer(I4), intent(in) :: wfnum, wf_ftype, wf_int_ft, ws_ftype, we_ftype
    integer(I4), intent(inout) :: mw_totn
    ! -- local
    integer(I4) :: i, w, size_n
    integer(I4), save :: wswe_flag = 0
    integer(I4), allocatable :: w_ij(:), w_ks(:), w_ke(:), cals2well(:)
    real(SP), allocatable :: w_val(:)
    !-------------------------------------------------------------------------------------
    allocate(w_ij(ncals), w_ks(ncals), w_ke(ncals))
    allocate(w_val(ncals))
    !$omp parallel workshare
    w_ij(:) = 0 ; w_ks(:) = 0 ; w_ke(:) = 0
    w_val(:) = SZERO
    !$omp end parallel workshare

    call set_2dfile2cals(wfnum, wf_ftype, wf_int_ft, SZERO, w_val)

    if (wswe_flag == 0) then
      if (ws_ftype == in_type(3) .or. ws_ftype == in_type(4)) then
        call set_2dfile2cals(wells_fnum, ws_ftype, wf_int_ft, 0, w_ks)
      end if

      if (we_ftype == in_type(3) .or. we_ftype == in_type(4)) then
        call set_2dfile2cals(welle_fnum, we_ftype, wf_int_ft, 0, w_ke)
      end if

      if (my_rank == 0) then
        close(wells_fnum) ; close(welle_fnum)
      end if

      wswe_flag = 1
    end if

    allocate(cals2well(ncals))
    !$omp parallel workshare
    cals2well(:) = 0
    !$omp end parallel workshare

    call set_cals2well(w_ks, w_ke, w_val, mw_totn, cals2well, w_ij)

    if (allocated(well_nflag)) then
      size_n = size(well_nflag)
      if (mw_totn < size_n) then
        mw_totn = size_n
      end if
      call set_multiwell(mw_totn, ncals, cals2well, w_ij, w_ks, w_ke, w_val)
    else
      allocate(well_nflag(mw_totn), st_well%ij(mw_totn))
      allocate(st_well%ks(mw_totn), st_well%ke(mw_totn))
      allocate(st_well%value(mw_totn))
      !$omp parallel workshare
      well_nflag(:) = 0 ; st_well%ij(:) = 0 ; st_well%ks(:) = 0 ; st_well%ke(:) = 0
      st_well%value(:) = SZERO
      !$omp end parallel workshare

      !$omp parallel do private(i, w)
      do i = 1, ncals
        w = cals2well(i)
        if (w /= 0) then
          st_well%ij(w) = w_ij(i) ; st_well%ks(w) = w_ks(i) ; st_well%ke(w) = w_ke(i)
          st_well%value(w) = w_val(i) ; well_nflag(w) = w
        end if
      end do
      !$omp end parallel do
    end if

    deallocate(w_ij, w_ks, w_ke, w_val, cals2well)

  end subroutine set_2dwell

  subroutine set_cals2well(wks, wke, wval, wnum, s2w, wij)
  !***************************************************************************************
  ! set_cals2well -- Set relationship between cals and well
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wks(:), wke(:)
    real(SP), intent(in) :: wval(:)
    integer(I4), intent(inout) :: wnum
    integer(I4), intent(inout) :: s2w(:)
    integer(I4), intent(out) :: wij(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    do i = 1, ncals
      if (wval(i) /= SZERO .and. wks(i) /= 0 .and. wke(i) /= 0) then
        wij(i) = loc2glo_ij(i)
        if (s2w(i) == 0) then
          wnum = wnum + 1
          s2w(i) = wnum
        end if
      end if
    end do

  end subroutine set_cals2well

  subroutine set_3dfile2well(wfnum, ftype, int_ft, wnum, w2c)
  !***************************************************************************************
  ! set_3dfile2well -- Set well value from 3d file
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wfnum, ftype, int_ft
    integer(I4), intent(out) :: wnum
    integer(I4), intent(inout) :: w2c(:)
    ! -- local
    integer(I4) :: i, cnum
    integer(I4) :: gridx, gridy, gridz, gridxyz
    real(SP), allocatable :: array_flat(:), val_out(:), array_val(:)
    real(SP), allocatable :: array_read(:,:,:)
    !-------------------------------------------------------------------------------------
    cnum = size(w2c(:))
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz ; gridxyz = st_grid%nxyz
    allocate(val_out(cnum))
    !$omp parallel workshare
    val_out(:) = SZERO
    !$omp end parallel workshare

    if (ftype == in_type(5) .or. int_ft == in_type(5)) then
      allocate(array_flat(gridxyz))
      !$omp parallel workshare
      array_flat(:) = SZERO
      !$omp end parallel workshare
      if (my_rank == 0) then
        allocate(array_read(gridx,gridy,gridz))
        !$omp parallel workshare
        array_read(:,:,:) = SZERO
        !$omp end parallel workshare
        call read_3dtxt(wfnum, gridx, gridy, gridz, array_read)
        call flat_3dto3d(gridx, gridy, gridz, array_read, array_flat)
        deallocate(array_read)
      end if
#ifdef MPI_MSG
      ! -- Scatter xyz value (xyzval)
        call scatter_xyzval(cnum, loc2glo_ijk, array_flat, val_out)
#else
      call set_calcr4(cnum, loc2glo_ijk, SZERO, array_flat, val_out)
#endif

    else if (ftype == in_type(6) .or. int_ft == in_type(6)) then
#ifdef MPI_MSG
      call read_mpi_file(ftype, int_ft, wfnum, cnum, val_out)
#else
      allocate(array_read(gridx,gridy,gridz), array_flat(gridxyz))
      !$omp parallel workshare
      array_read(:,:,:) = SZERO ; array_flat(:) = SZERO
      !$omp end parallel workshare
      call read_3dbin(wfnum, gridx, gridy, gridz, SZERO, array_read)
      call flat_3dto3d(gridx, gridy, gridz, array_read, array_flat)
      deallocate(array_read)
      call set_calcr4(cnum, loc2glo_ijk, SZERO, array_flat, val_out)
#endif
    end if

    if (allocated(array_flat)) then
      deallocate(array_flat)
    end if

    allocate(array_val(cnum))
    !$omp parallel workshare
    array_val(:) = SZERO
    !$omp end parallel workshare

    wnum = 0
    do i = 1, cnum
      if (val_out(i) /= SZERO) then
        wnum = wnum + 1
        w2c(wnum) = i
        array_val(wnum) = val_out(i)
      end if
    end do

    allocate(st_well%value(wnum))
    !$omp parallel workshare
    st_well%value(:) = array_val(1:wnum)
    !$omp end parallel workshare

    deallocate(val_out, array_val)

  end subroutine set_3dfile2well

  subroutine set_well2index(wnum)
  !***************************************************************************************
  ! set_well2index -- Set well index
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wnum
    ! -- local
    integer(I4) :: i, k, index_count, c_num
    integer(I4), allocatable :: temp_wconn(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_wconn(st_grid%nz*wnum), well_index(0:wnum))
    !$omp parallel workshare
    temp_wconn(:) = 0 ; well_index(:) = 0
    !$omp end parallel workshare

    index_count = 0
    do i = 1, wnum
      do k = st_well%ks(i), st_well%ke(i)
        index_count = index_count + 1
        c_num = st_grid%nx*st_grid%ny*(k-1) + st_well%ij(i)
        temp_wconn(index_count) = findloc(loc2glo_ijk(:), value = c_num, dim = 1)
      end do
      well_index(i) = index_count
    end do

    allocate(well_conn(index_count))
    !$omp parallel workshare
    well_conn(:) = temp_wconn(1:index_count)
    !$omp end parallel workshare

    deallocate(temp_wconn)

  end subroutine set_well2index

  subroutine set_well3d2index(wnum, w2c)
  !***************************************************************************************
  ! set_well3d2index -- Set well index from 3d well
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wnum
    integer(I4), intent(in) :: w2c(:)
    ! -- local
    integer(I4) :: i, index_count
    !-------------------------------------------------------------------------------------
    allocate(well_conn(wnum), well_index(0:wnum))
    !$omp parallel workshare
    well_conn(:) = 0 ; well_index(:) = 0
    !$omp end parallel workshare

    index_count = 0
    do i = 1, wnum
      index_count = index_count + 1
      well_conn(index_count) = w2c(i)
      well_index(i) = index_count
    end do

  end subroutine set_well3d2index

  subroutine set_wellprop(wnum, wtop, wbott)
  !***************************************************************************************
  ! set_wellprop -- Set well property (top ,bottom)
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: wnum
    real(DP), intent(out) :: wtop(:), wbott(:)
    ! -- local
    integer(I4) :: i, bott_calc, surf_calc
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, surf_calc, bott_calc)
    do i = 1, wnum
      surf_calc = well_conn(well_index(i-1)+1)
      bott_calc = well_conn(well_index(i))
      wtop(i) = cell_top(surf_calc)
      wbott(i) = cell_bot(bott_calc)
    end do
    !$omp end parallel do

  end subroutine set_wellprop

  subroutine set_connect(cksx, cksy, cksz)
  !***************************************************************************************
  ! set_connect -- Set connectivity
  !***************************************************************************************
    ! -- modules
    use constval_module, only: FACE
    use initial_module, only: st_out_type
    use set_cell, only: calc2reg, seal_cnum, glo2loc_ijk, get_calc_grid
    use make_cell, only: dis2face, face_area, area_r
#ifdef MPI_MSG
    use initial_module, only: pro_totn
    use set_cell, only: neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                        recv_citem, send2recv
    use mpi_set, only: senrec_faceval, senrec_neibval
#endif
    ! -- inout
    real(SP), intent(in) :: cksx(:), cksy(:), cksz(:)
    ! -- local
    integer(I4) :: i, off_num, sflag, mpi_send, left_num, right_num, totn_clac
    integer(I4) :: i_num, j_num, k_num, g_num, n_num, w_num, e_num, s_num, u_num, d_num
    integer(I4) :: loc_n, loc_w, loc_e, loc_s, loc_u, loc_d
    integer(I4) :: gridx, gridy, gridz
    real(DP) :: dis1, dis2, dis12, hyd_face, tks, aks
    real(DP), allocatable :: reg_cksx(:), reg_cksy(:), reg_cksz(:)
    real(DP), allocatable :: reg_fare(:,:)
#ifdef MPI_MSG
    integer(I4) :: calc_num
    real(SP), allocatable :: mpi_cksx(:), mpi_cksy(:), mpi_cksz(:)
#endif
    !-------------------------------------------------------------------------------------
    totn_clac = ncalc + neib_ncalc
    gridx = st_grid%nx ; gridy = st_grid%ny ; gridz = st_grid%nz
    allocate(off_row(totn_clac*FACE), off_index(0:totn_clac), conn_dir(totn_clac*FACE))
    allocate(area_dis(totn_clac*FACE), sat_hydf(totn_clac*FACE))
    allocate(surf_area(ncals), rech_area(ncals))
    allocate(hydf_surf(ncals), abyd_surf(ncals))
    allocate(sea_abyd(seal_cnum*FACE), sea_hydf(seal_cnum*FACE))
    allocate(left_off(totn_clac*FACE), right_off(totn_clac*FACE))
    !$omp parallel workshare
    off_row(:) = 0 ; off_index(:) = 0 ; conn_dir(:) = 0
    area_dis(:) = DZERO ; sat_hydf(:) = DZERO
    surf_area(:) = DZERO ; rech_area(:) = DZERO
    hydf_surf(:) = DZERO ; abyd_surf(:) = DZERO
    sea_abyd(:) = DZERO ; sea_hydf(:) = DZERO
    left_off(:) = 0 ; right_off(:) = 0
    !$omp end parallel workshare

    if (st_out_type%velc > 0) then
      allocate(conn_dis(totn_clac*FACE))
      allocate(sea_dis(seal_cnum*FACE), sea_dir(seal_cnum*FACE))
      !$omp parallel workshare
      conn_dis(:) = DZERO ; sea_dis(:) = DZERO ; sea_dir(:) = 0
      !$omp end parallel workshare
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      allocate(reg_dis(totn_clac,FACE), reg_fare(totn_clac,FACE))
      !$omp parallel workshare
      reg_dis(:,:) = DZERO ; reg_fare(:,:) = DZERO
      !$omp end parallel workshare
      ! -- Send and Receive face value (faceval)
        call senrec_faceval(neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                            recv_citem, dis2face, reg_dis)
        call senrec_faceval(neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                            recv_citem, face_area, reg_fare)
      !$omp parallel workshare
      reg_dis(1:ncalc,:) = dis2face(1:ncalc,:)
      reg_fare(1:ncalc,:) = face_area(1:ncalc,:)
      !$omp end parallel workshare

      allocate(reg_cksx(totn_clac), mpi_cksx(totn_clac))
      !$omp parallel workshare
      reg_cksx(:) = DZERO ; mpi_cksx(:) = SZERO
      !$omp end parallel workshare
      ! -- Send and Receive neighbor value (neibval)
        call senrec_neibval(neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                            recv_citem, cksx, mpi_cksx)
      !$omp parallel workshare
      reg_cksx(1:ncalc) = real(cksx(1:ncalc), kind=DP)
      reg_cksx(ncalc+1:totn_clac) = real(mpi_cksx(ncalc+1:totn_clac), kind=DP)
      !$omp end parallel workshare
      deallocate(mpi_cksx)

      allocate(reg_cksy(totn_clac), mpi_cksy(totn_clac))
      !$omp parallel workshare
      reg_cksy(:) = DZERO ; mpi_cksy(:) = SZERO
      !$omp end parallel workshare
      ! -- Send and Receive neighbor value (neibval)
        call senrec_neibval(neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                            recv_citem, cksy, mpi_cksy)
      !$omp parallel workshare
      reg_cksy(1:ncalc) = real(cksy(1:ncalc), kind=DP)
      reg_cksy(ncalc+1:totn_clac) = real(mpi_cksy(ncalc+1:totn_clac), kind=DP)
      !$omp end parallel workshare
      deallocate(mpi_cksy)

      allocate(reg_cksz(totn_clac), mpi_cksz(totn_clac))
      !$omp parallel workshare
      reg_cksz(:) = DZERO ; mpi_cksz(:) = SZERO
      !$omp end parallel workshare
      ! -- Send and Receive neighbor value (neibval)
        call senrec_neibval(neib_mpi_totn, neib_num, send_cind, recv_cind, send_citem,&
                            recv_citem, cksz, mpi_cksz)
      !$omp parallel workshare
      reg_cksz(1:ncalc) = real(cksz(1:ncalc), kind=DP)
      reg_cksz(ncalc+1:totn_clac) = real(mpi_cksz(ncalc+1:totn_clac), kind=DP)
      !$omp end parallel workshare
      deallocate(mpi_cksz)
    else
      allocate(reg_dis(ncalc,FACE), reg_fare(ncalc,FACE))
      allocate(reg_cksx(ncalc), reg_cksy(ncalc), reg_cksz(ncalc))
      !$omp parallel workshare
      reg_dis(:,:) = dis2face(:,:) ; reg_fare(:,:) = face_area(:,:)
      reg_cksx(:) = cksx(:) ; reg_cksy(:) = cksy(:) ; reg_cksz(:) = cksz(:)
      !$omp end parallel workshare
      neib_ncalc = 0
    end if
#else
    allocate(reg_dis(ncalc,FACE), reg_fare(ncalc,FACE))
    allocate(reg_cksx(ncalc), reg_cksy(ncalc), reg_cksz(ncalc))
    !$omp parallel workshare
    reg_dis(:,:) = dis2face(:,:) ; reg_fare(:,:) = face_area(:,:)
    reg_cksx(:) = cksx(:) ; reg_cksy(:) = cksy(:) ; reg_cksz(:) = cksz(:)
    !$omp end parallel workshare
    neib_ncalc = 0
#endif

    allocate(surf_dis(ncals), sea2cal(seal_cnum*FACE), sea2sea(seal_cnum*FACE))
    !$omp parallel workshare
    surf_dis(:) = dis2face(1:ncals,1) ; sea2cal(:) = 0 ; sea2sea(:) = 0
    !$omp end parallel workshare
    deallocate(dis2face, face_area)

    nseal = 0 ; tconn_num = 0 ; mpi_send = 0 ; left_num = 0 ; right_num = 0
    do i = 1, totn_clac
      ! -- Get calculation number from grid number (calc_grid)
        call get_calc_grid(i, i_num, j_num, k_num)
      g_num = loc2glo_ijk(i)
      ! up direction
      if (k_num /= 1) then
        u_num = g_num-gridx*gridy ; loc_u = glo2loc_ijk(u_num)
        if (loc_u > 0) then
          off_num = loc_u ; tks = reg_cksz(i)
          if (off_num > totn_clac .and. i <= ncalc) then !sea grid
            sflag = 1 ; nseal = nseal + 1 ; sea2cal(nseal) = i
            sea2sea(nseal) = off_num - totn_clac
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(1, sflag, i, off_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            sea_abyd(nseal) = reg_fare(i,1)/dis12 ; sea_hydf(nseal) = hyd_face
            if (st_out_type%velc > 0) then
              sea_dis(nseal) = DONE/dis12 ; sea_dir(nseal) = 1
            end if
          else if (0 < off_num .and. off_num <= ncalc) then
            if (calc2reg(i) == calc2reg(off_num) .or. st_sim%reg_neib == 1) then
              sflag = 0 ; aks = reg_cksz(off_num) ; tconn_num = tconn_num + 1
              off_row(tconn_num) = off_num ; conn_dir(tconn_num) = 1
              left_num = left_num + 1 ; left_off(tconn_num) = left_num
              ! -- Set distance between adjacent cell (dis_adj)
                call set_dis_adj(1, sflag, i, off_num, dis1, dis2, dis12)
              ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
                call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
              area_dis(tconn_num) = reg_fare(i,1)/dis12 ; sat_hydf(tconn_num) = hyd_face
              if (st_out_type%velc > 0) then
                conn_dis(tconn_num) = DONE/dis12
              end if
            end if
          else if (i <= ncalc) then !for mpi
#ifdef MPI_MSG
            mpi_send = mpi_send + 1 ; calc_num = recv_citem(send2recv(mpi_send))
            tks = reg_cksz(i) ; aks = reg_cksz(calc_num) ; tconn_num = tconn_num + 1
            off_row(tconn_num) = calc_num ; sflag = 0 ; conn_dir(tconn_num) = 1
            left_num = left_num + 1 ; left_off(tconn_num) = left_num
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(1, sflag, i, calc_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            area_dis(tconn_num) = reg_fare(i,1)/dis12 ; sat_hydf(tconn_num) = hyd_face
            if (st_out_type%velc > 0) then
              conn_dis(tconn_num) = DONE/dis12
            end if
#endif
          end if
        end if
      else if (i <= ncals) then
        surf_area(i) = reg_fare(i,1)
        rech_area(i) = area_r(i)
        hydf_surf(i) = reg_cksz(i)
      end if
      ! north direction
      if (j_num /= 1) then
        n_num = g_num-gridx ; loc_n = glo2loc_ijk(n_num)
        if (loc_n > 0) then
          off_num = loc_n ; tks = reg_cksy(i)
          if (off_num > totn_clac .and. i <= ncalc) then !sea grid
            sflag = 1 ; nseal = nseal + 1 ; sea2cal(nseal) = i
            sea2sea(nseal) = off_num - totn_clac
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(2, sflag, i, off_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            sea_abyd(nseal) = reg_fare(i,2)/dis12 ; sea_hydf(nseal) = hyd_face
            if (st_out_type%velc > 0) then
              sea_dis(nseal) = DONE/dis12 ; sea_dir(nseal) = 2
            end if
          else if (0 < off_num .and. off_num <= ncalc) then
            if (calc2reg(i) == calc2reg(off_num) .or. st_sim%reg_neib == 1) then
              sflag = 0 ; aks = reg_cksy(off_num) ; tconn_num = tconn_num + 1
              off_row(tconn_num) = off_num ; conn_dir(tconn_num) = 2
              left_num = left_num + 1 ; left_off(tconn_num) = left_num
              ! -- Set distance between adjacent cell (dis_adj)
                call set_dis_adj(2, sflag, i, off_num, dis1, dis2, dis12)
              ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
                call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
              area_dis(tconn_num) = reg_fare(i,2)/dis12 ; sat_hydf(tconn_num) = hyd_face
              if (st_out_type%velc > 0) then
                conn_dis(tconn_num) = DONE/dis12
              end if
            end if
          else if (i <= ncalc) then !for mpi
#ifdef MPI_MSG
            mpi_send = mpi_send + 1 ; calc_num = recv_citem(send2recv(mpi_send))
            tks = reg_cksy(i) ; aks = reg_cksy(calc_num) ; tconn_num = tconn_num + 1
            off_row(tconn_num) = calc_num ; sflag = 0 ; conn_dir(tconn_num) = 2
            left_num = left_num + 1 ; left_off(tconn_num) = left_num
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(2, sflag, i, calc_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            area_dis(tconn_num) = reg_fare(i,2)/dis12 ; sat_hydf(tconn_num) = hyd_face
            if (st_out_type%velc > 0) then
              conn_dis(tconn_num) = DONE/dis12
            end if
#endif
          end if
        end if
      end if
      ! west direction
      if (i_num /= 1) then
        w_num = g_num-1 ; loc_w = glo2loc_ijk(w_num)
        if (loc_w > 0) then
          off_num = loc_w ; tks = reg_cksx(i)
          if (off_num > totn_clac .and. i <= ncalc) then !sea grid
            sflag = 1 ; nseal = nseal + 1 ; sea2cal(nseal) = i
            sea2sea(nseal) = off_num - totn_clac
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(3, sflag, i, off_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            sea_abyd(nseal) = reg_fare(i,3)/dis12 ; sea_hydf(nseal) = hyd_face
            if (st_out_type%velc > 0) then
              sea_dis(nseal) = DONE/dis12 ; sea_dir(nseal) = 3
            end if
          else if (0 < off_num .and. off_num <= ncalc) then
            if (calc2reg(i) == calc2reg(off_num) .or. st_sim%reg_neib == 1) then
              sflag = 0 ; aks = reg_cksx(off_num) ; tconn_num = tconn_num + 1
              off_row(tconn_num) = off_num ; conn_dir(tconn_num) = 3
              left_num = left_num + 1 ; left_off(tconn_num) = left_num
              ! -- Set distance between adjacent cell (dis_adj)
                call set_dis_adj(3, sflag, i, off_num, dis1, dis2, dis12)
              ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
                call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
              area_dis(tconn_num) = reg_fare(i,3)/dis12 ; sat_hydf(tconn_num) = hyd_face
              if (st_out_type%velc > 0) then
                conn_dis(tconn_num) = DONE/dis12
              end if
            end if
          else if (i <= ncalc) then !for mpi
#ifdef MPI_MSG
            mpi_send = mpi_send + 1 ; calc_num = recv_citem(send2recv(mpi_send))
            tks = reg_cksx(i) ; aks = reg_cksx(calc_num) ; tconn_num = tconn_num + 1
            off_row(tconn_num) = calc_num ; sflag = 0 ; conn_dir(tconn_num) = 3
            left_num = left_num + 1 ; left_off(tconn_num) = left_num
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(3, sflag, i, calc_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            area_dis(tconn_num) = reg_fare(i,3)/dis12 ; sat_hydf(tconn_num) = hyd_face
            if (st_out_type%velc > 0) then
              conn_dis(tconn_num) = DONE/dis12
            end if
#endif
          end if
        end if
      end if
      ! east direction
      if (i_num /= gridx) then
        e_num = g_num+1 ; loc_e = glo2loc_ijk(e_num)
        if (loc_e > 0) then
          off_num = loc_e ; tks = reg_cksx(i)
          if (off_num > totn_clac .and. i <= ncalc) then !sea grid
            sflag = 1 ; nseal = nseal + 1 ; sea2cal(nseal) = i
            sea2sea(nseal) = off_num - totn_clac
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(4, sflag, i, off_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            sea_abyd(nseal) = reg_fare(i,4)/dis12 ; sea_hydf(nseal) = hyd_face
            if (st_out_type%velc > 0) then
              sea_dis(nseal) = DONE/dis12 ; sea_dir(nseal) = 4
            end if
          else if (0 < off_num .and. off_num <= ncalc) then
            if (calc2reg(i) == calc2reg(off_num) .or. st_sim%reg_neib == 1) then
              sflag = 0 ; aks = reg_cksx(off_num) ; tconn_num = tconn_num + 1
              off_row(tconn_num) = off_num ; conn_dir(tconn_num) = 4
              right_num = right_num + 1 ; right_off(tconn_num) = right_num
              ! -- Set distance between adjacent cell (dis_adj)
                call set_dis_adj(4, sflag, i, off_num, dis1, dis2, dis12)
              ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
                call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
              area_dis(tconn_num) = reg_fare(i,4)/dis12 ; sat_hydf(tconn_num) = hyd_face
              if (st_out_type%velc > 0) then
                conn_dis(tconn_num) = DONE/dis12
              end if
            end if
          else if (i <= ncalc) then !for mpi
#ifdef MPI_MSG
            mpi_send = mpi_send + 1 ; calc_num = recv_citem(send2recv(mpi_send))
            tks = reg_cksx(i) ; aks = reg_cksx(calc_num) ; tconn_num = tconn_num + 1
            off_row(tconn_num) = calc_num ; sflag = 0 ; conn_dir(tconn_num) = 4
            right_num = right_num + 1 ; right_off(tconn_num) = right_num
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(4, sflag, i, calc_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            area_dis(tconn_num) = reg_fare(i,4)/dis12 ; sat_hydf(tconn_num) = hyd_face
            if (st_out_type%velc > 0) then
              conn_dis(tconn_num) = DONE/dis12
            end if
#endif
          end if
        end if
      end if
      ! south direction
      if (j_num /= gridy) then
        s_num = g_num+gridx ; loc_s = glo2loc_ijk(s_num)
        if (loc_s > 0) then
          off_num = loc_s ; tks = reg_cksy(i)
          if (off_num > totn_clac .and. i <= ncalc) then !sea grid
            sflag = 1 ; nseal = nseal + 1 ; sea2cal(nseal) = i
            sea2sea(nseal) = off_num - totn_clac
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(5, sflag, i, off_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            sea_abyd(nseal) = reg_fare(i,5)/dis12 ; sea_hydf(nseal) = hyd_face
            if (st_out_type%velc > 0) then
              sea_dis(nseal) = DONE/dis12 ; sea_dir(nseal) = 5
            end if
          else if (0 < off_num .and. off_num <= ncalc) then
            if (calc2reg(i) == calc2reg(off_num) .or. st_sim%reg_neib == 1) then
              sflag = 0 ; aks = reg_cksy(off_num) ; tconn_num = tconn_num + 1
              off_row(tconn_num) = off_num ; conn_dir(tconn_num) = 5
              right_num = right_num + 1 ; right_off(tconn_num) = right_num
              ! -- Set distance between adjacent cell (dis_adj)
                call set_dis_adj(5, sflag, i, off_num, dis1, dis2, dis12)
              ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
                call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
              area_dis(tconn_num) = reg_fare(i,5)/dis12 ; sat_hydf(tconn_num) = hyd_face
              if (st_out_type%velc > 0) then
                conn_dis(tconn_num) = DONE/dis12
              end if
            end if
          else if (i <= ncalc) then !for mpi
#ifdef MPI_MSG
            mpi_send = mpi_send + 1 ; calc_num = recv_citem(send2recv(mpi_send))
            tks = reg_cksy(i) ; aks = reg_cksy(calc_num) ; tconn_num = tconn_num + 1
            off_row(tconn_num) = calc_num ; sflag = 0 ; conn_dir(tconn_num) = 5
            right_num = right_num + 1 ; right_off(tconn_num) = right_num
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(5, sflag, i, calc_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            area_dis(tconn_num) = reg_fare(i,5)/dis12 ; sat_hydf(tconn_num) = hyd_face
            if (st_out_type%velc > 0) then
              conn_dis(tconn_num) = DONE/dis12
            end if
#endif
          end if
        end if
      end if
      ! down direction
      if (k_num /= gridz) then
        d_num = g_num+gridx*gridy ; loc_d = glo2loc_ijk(d_num)
        if (loc_d > 0) then
          off_num = loc_d ; tks = reg_cksz(i)
          if (off_num > totn_clac .and. i <= ncalc) then !sea grid
            sflag = 1 ; nseal = nseal + 1 ; sea2cal(nseal) = i
            sea2sea(nseal) = off_num - totn_clac
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(6, sflag, i, off_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            sea_abyd(nseal) = reg_fare(i,6)/dis12 ; sea_hydf(nseal) = hyd_face
            if (st_out_type%velc > 0) then
              sea_dis(nseal) = DONE/dis12 ; sea_dir(nseal) = 6
            end if
          else if (0 < off_num .and. off_num <= ncalc) then
            if (calc2reg(i) == calc2reg(off_num) .or. st_sim%reg_neib == 1) then
              sflag = 0 ; aks = reg_cksz(off_num) ; tconn_num = tconn_num + 1
              off_row(tconn_num) = off_num ; conn_dir(tconn_num) = 6
              right_num = right_num + 1 ; right_off(tconn_num) = right_num
              ! -- Set distance between adjacent cell (dis_adj)
                call set_dis_adj(6, sflag, i, off_num, dis1, dis2, dis12)
              ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
                call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
              area_dis(tconn_num) = reg_fare(i,6)/dis12 ; sat_hydf(tconn_num) = hyd_face
              if (st_out_type%velc > 0) then
                conn_dis(tconn_num) = DONE/dis12
              end if
            end if
          else if (i <= ncalc) then !for mpi
#ifdef MPI_MSG
            mpi_send = mpi_send + 1 ; calc_num = recv_citem(send2recv(mpi_send))
            tks = reg_cksz(i) ; aks = reg_cksz(calc_num) ; tconn_num = tconn_num + 1
            off_row(tconn_num) = calc_num ; sflag = 0 ; conn_dir(tconn_num) = 6
            right_num = right_num + 1 ; right_off(tconn_num) = right_num
            ! -- Set distance between adjacent cell (dis_adj)
              call set_dis_adj(6, sflag, i, calc_num, dis1, dis2, dis12)
            ! -- Set saturated hydradulic conductivity by harmonic mean (sat_hyd)
              call set_sat_hyd(sflag, tks, aks, dis1, dis2, hyd_face)
            area_dis(tconn_num) = reg_fare(i,6)/dis12 ; sat_hydf(tconn_num) = hyd_face
            if (st_out_type%velc > 0) then
              conn_dis(tconn_num) = DONE/dis12
            end if
#endif
          end if
        end if
      end if
      off_index(i) = tconn_num
    end do

    off_index(0) = 0

    deallocate(reg_cksx, reg_cksy, reg_cksz)
    deallocate(reg_dis, reg_fare)
    deallocate(area_r)
    deallocate(calc2reg)
    deallocate(glo2loc_ijk)

  end subroutine set_connect

  subroutine set_wellconn(wnum, wksx, wksy)
  !***************************************************************************************
  ! set_wellconn -- Set well connectivity
  !***************************************************************************************
    ! -- modules
    use constval_module, only: SHALF
    ! -- inout
    integer(I4), intent(in) :: wnum
    real(SP), intent(in) :: wksx(:), wksy(:)
    ! -- local
    integer(I4) :: i, j, k
    real(SP) :: cksxy
    !-------------------------------------------------------------------------------------
    allocate(abyd_well(well_index(wnum)))
    !$omp parallel
    !$omp workshare
    abyd_well(:) = DZERO
    !$omp end workshare

    !$omp do private(i, j, k, cksxy)
    do i = 1, wnum
      do k = well_index(i-1)+1, well_index(i)
        j = well_conn(k)
        cksxy = (wksx(j)*wksy(j))**SHALF
        abyd_well(k) = real(cksxy, kind=DP)*(cell_top(j)-cell_bot(j))
      end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine set_wellconn

  subroutine set_srabyd(surfn, bott_s, area_s, surf2surf, abyd_surftemp)
  !***************************************************************************************
  ! set_srabyd -- Set surface&recharge area and area by distance
  !***************************************************************************************
    ! -- modules
    use make_cell, only: cell_cent
    ! -- inout
    integer(I4), intent(in) :: surfn, surf2surf(:)
    real(DP), intent(in) :: bott_s(:), area_s(:)
    real(DP), intent(out) :: abyd_surftemp(:)
    ! -- local
    integer(I4) :: i, s
    real(DP) :: dis_s
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, s, dis_s)
    do i = 1, surfn
      s = surf2surf(i)
      dis_s = abs(bott_s(i) - cell_cent(s))
      abyd_surftemp(i) = abs(area_s(i))/dis_s
      surf_area(s) = surf_area(s) - area_s(i)
      rech_area(s) = rech_area(s) - area_s(i)
      if (surf_area(s) < DZERO) then
        surf_area(s) = DZERO
      end if
      if (rech_area(s) < DZERO) then
        rech_area(s) = DZERO
      end if
    end do
    !$omp end parallel do

  end subroutine set_srabyd

  subroutine set_chabyd()
  !***************************************************************************************
  ! set_chabyd -- Set charge area by distance
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, ncals
      abyd_surf(i) = surf_area(i)/surf_dis(i)
    end do
    !$omp end parallel do

  end subroutine set_chabyd

  subroutine set_dis_adj(dir, bf, targc, adjac, d1, d2, d12)
  !***************************************************************************************
  ! set_dis_adj -- Set distance between adjacent cell
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: dir, bf, targc, adjac
    real(DP), intent(out) :: d1, d2, d12
    ! -- local

    !-------------------------------------------------------------------------------------
    d1 = reg_dis(targc,dir)
    if (bf == 1) then
      d2 = d1
      d12 = d1
    else
      d2 = reg_dis(adjac,7-dir)
      d12 = d1 + d2
    end if

  end subroutine set_dis_adj

  subroutine set_sat_hyd(bf, hyd_c1, hyd_c2, d1, d2, hydf)
  !***************************************************************************************
  ! set_sat_hyd -- Set saturated hydradulic conductivity by harmonic mean
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: bf
    real(DP), intent(in) :: d1, d2, hyd_c1, hyd_c2
    real(DP), intent(out) :: hydf
    ! -- local
    real(DP) :: c1, c2, d12, cd, cd12
    !-------------------------------------------------------------------------------------
    if (bf == 1) then
      hydf = hyd_c1
    else
      d12 = d1 + d2
      c1 = hyd_c1*d12
      c2 = hyd_c2*d12
      cd12 = d1*c2+d2*c1
      cd = DONE/cd12
      hydf = c1*c2*cd
    end if

  end subroutine set_sat_hyd

  subroutine set_multiwell(max_num, new_num, mw_id, mwij, mwks, mwke, mwv)
  !***************************************************************************************
  ! set_multiwell -- Set existing multi well
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: max_num, new_num
    integer(I4), intent(in) :: mw_id(:), mwij(:), mwks(:), mwke(:)
    real(SP), intent(in) :: mwv(:)
    ! -- local
    integer(I4) :: i, nsize, n, m, plus, np
    integer(I4), allocatable :: temp_ij(:), temp_ks(:), temp_ke(:)
    integer(I4), allocatable :: temp_wflag(:)
    real(SP), allocatable :: temp_value(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_ij(max_num), temp_ks(max_num), temp_ke(max_num))
    allocate(temp_wflag(max_num), temp_value(max_num))
    !$omp parallel workshare
    temp_ij(:) = 0 ; temp_ks(:) = 0 ; temp_ke(:) = 0
    temp_wflag(:) = 0 ; temp_value(:) = SZERO
    !$omp end parallel workshare

    nsize = size(well_nflag)
    !$omp parallel do private(i)
    do i = 1, nsize
      temp_wflag(i) = well_nflag(i)
      temp_ij(i) = st_well%ij(i) ; temp_ks(i) = st_well%ks(i) ; temp_ke(i) = st_well%ke(i)
    end do
    !$omp end parallel do
    deallocate(well_nflag)
    deallocate(st_well%ij, st_well%ks, st_well%ke)

    plus = 0 ; np = nsize
    new: do m = 1, new_num
      if (mw_id(m) == 0) then
        cycle
      end if
      old: do n = 1, nsize
        if (temp_wflag(n) == mw_id(m)) then
          temp_value(n) = mwv(m)
          cycle new
        end if
      end do old
      plus = plus + 1
      np = nsize + plus
      temp_wflag(np) = mw_id(m) ; temp_ij(np) = mwij(m)
      temp_ks(np) = mwks(m) ; temp_ke(np) = mwke(m)
    end do new

    allocate(st_well%ij(np), st_well%ks(np), st_well%ke(np))
    allocate(st_well%value(np), well_nflag(np))
    !$omp parallel workshare
    st_well%ij(:) = temp_ij(1:np)
    st_well%ks(:) = temp_ks(1:np) ; st_well%ke(:) = temp_ke(1:np)
    st_well%value(:) = temp_value(1:np) ; well_nflag(:) = temp_wflag(1:np)
    !$omp end parallel workshare

    deallocate(temp_ij, temp_ks, temp_ke)
    deallocate(temp_value, temp_wflag)

  end subroutine set_multiwell

  subroutine set_calc_flag_int(val_num, nonval, inval, tg_flag)
  !***************************************************************************************
  ! set_calc_flag_int -- Set calculation target flag from integer
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: val_num
    integer(I4), intent(in) :: nonval
    integer(I4), intent(in) :: inval(:)
    integer(I4), intent(out) :: tg_flag(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, val_num
      if (inval(i) /= nonval) then
        tg_flag(i) = 1
      end if
    end do
    !$omp end parallel do

  end subroutine set_calc_flag_int

  subroutine set_calc_flag_real4(val_num, nonval, inval, tg_flag)
  !***************************************************************************************
  ! set_calc_flag_real4 -- Set calculation target flag from real4
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: val_num
    real(SP), intent(in) :: nonval
    real(SP), intent(in) :: inval(:)
    integer(I4), intent(out) :: tg_flag(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, val_num
      if (inval(i) /= nonval) then
        tg_flag(i) = 1
      end if
    end do
    !$omp end parallel do

  end subroutine set_calc_flag_real4

  subroutine set_calc_flag_real8(val_num, nonval, inval, tg_flag)
  !***************************************************************************************
  ! set_calc_flag_real8 -- Set calculation target flag from real8
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: val_num
    real(DP), intent(in) :: nonval
    real(DP), intent(in) :: inval(:)
    integer(I4), intent(out) :: tg_flag(:)
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, val_num
      if (inval(i) /= nonval) then
        tg_flag(i) = 1
      end if
    end do
    !$omp end parallel do

  end subroutine set_calc_flag_real8

  subroutine count_flag(tg_flag, count)
  !***************************************************************************************
  ! count_flag -- Count flag number
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: tg_flag(:)
    integer(I4), intent(out) :: count
    ! -- local

    !-------------------------------------------------------------------------------------
    count = 0
    !$omp parallel workshare
    count = sum(tg_flag)
    !$omp end parallel workshare

  end subroutine count_flag

  subroutine set_calci4(loc_num, cal2glo, noval, in_val, out_val)
  !***************************************************************************************
  ! set_calci4 -- Set calculation integer value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: cal2glo(:)
    integer(I4), intent(in) :: noval
    integer(I4), intent(in) :: in_val(:)
    integer(I4), intent(out) :: out_val(:)
    ! -- local
    integer(I4) :: i, c
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, c)
    do i = 1, loc_num
      c = cal2glo(i)
      if (in_val(c) /= noval) then
        out_val(i) = in_val(c)
      end if
    end do
    !$omp end parallel do

  end subroutine set_calci4

  subroutine set_calcr4(loc_num, cal2glo, noval, in_val, out_val)
  !***************************************************************************************
  ! set_calcr4 -- Set calculation real4 value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: cal2glo(:)
    real(SP), intent(in) :: noval
    real(SP), intent(in) :: in_val(:)
    real(SP), intent(out) :: out_val(:)
    ! -- local
    integer(I4) :: i, c
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, c)
    do i = 1, loc_num
      c = cal2glo(i)
      if (in_val(c) /= noval) then
        out_val(i) = in_val(c)
      end if
    end do
    !$omp end parallel do

  end subroutine set_calcr4

  subroutine set_calcr8(loc_num, cal2glo, noval, in_val, out_val)
  !***************************************************************************************
  ! set_calcr8 -- Set calculation real8 value
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: loc_num
    integer(I4), intent(in) :: cal2glo(:)
    real(DP), intent(in) :: noval
    real(DP), intent(in) :: in_val(:)
    real(DP), intent(out) :: out_val(:)
    ! -- local
    integer(I4) :: i, c
    !-------------------------------------------------------------------------------------
    !$omp parallel do private(i, c)
    do i = 1, loc_num
      c = cal2glo(i)
      if (in_val(c) /= noval) then
        out_val(i) = in_val(c)
      end if
    end do
    !$omp end parallel do

  end subroutine set_calcr8

end module set_condition
