module allocate_solution
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO
  use types_module, only: sol_set
  use set_cell, only: ncalc

  implicit none
  private
  public :: allocate_solvar
  integer(I4), public :: nreg_num
  integer(I4), allocatable, public :: dir_conn(:), dir_seal(:)
  integer(I4), allocatable, public :: left_offr(:), right_offr(:)

  type :: matrix_int
    integer(I4) :: unknow
    integer(I4) :: lunum
    integer(I4), allocatable :: offrow(:)
    integer(I4), allocatable :: offind(:)
  end type matrix_int
  type(matrix_int), allocatable, public :: crs_index(:)

  type :: matrix_real
    real(DP), allocatable :: dmat(:)
    real(DP), allocatable :: lumat(:)
    real(DP), allocatable :: x(:)
    real(DP), allocatable :: rhs(:)
  end type matrix_real
  type(matrix_real), allocatable, public :: array_var(:)

  type :: prolo_vec
    integer(I4), allocatable :: pindex(:)
    integer(I4), allocatable :: poffrow(:)
    real(DP), allocatable :: pval(:)
  end type prolo_vec
  type(prolo_vec), allocatable, public :: pro_var(:)

  type :: restr_vec
    integer(I4), allocatable :: rindex(:)
    integer(I4), allocatable :: roffrow(:)
    real(DP), allocatable :: rval(:)
  end type restr_vec
  type(restr_vec), allocatable, public :: res_var(:)

  ! -- local

  contains

  subroutine allocate_solvar(st_sol)
  !*********************************************************************************************
  ! allocate_solvar -- Allocate solution variable for time step
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local

    !-------------------------------------------------------------------------------------------
    ! -- Allocate timeupdate (timeup)
      call allocate_timeup(st_sol)

    ! -- Allocate matrix and vector (matvec)
      call allocate_matvec(st_sol)

  end subroutine allocate_solvar

  subroutine allocate_timeup(st_sol)
  !*********************************************************************************************
  ! allocate_timeup -- Allocate timeupdate
  !*********************************************************************************************
    ! -- modules
    use set_cell, only: ncals
    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    allocate(st_sol%head_old(ncalc), st_sol%srat_old(ncalc), st_sol%stor_old(ncalc))
    allocate(st_sol%surf_head(ncals), st_sol%surf_old(ncals), st_sol%surf_rati(ncals))
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      st_sol%head_old(i) = DZERO
      st_sol%srat_old(i) = DZERO ; st_sol%stor_old(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, ncals
      st_sol%surf_head(i) = DZERO
      st_sol%surf_old(i) = DZERO
      st_sol%surf_rati(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

  end subroutine allocate_timeup

  subroutine allocate_matvec(st_sol)
  !*********************************************************************************************
  ! allocate_matvec -- Allocate matrix and vector
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_ctrl, st_out_type
    use set_cell, only: neib_ncalc
    use set_condition, only: tconn_num, off_row, off_index, conn_dir, sea_dir, left_off, st_hydr
    use set_condition, only: st_bcnd
    use set_condition, only: right_off
    ! -- inout

    type(sol_set), intent(inout) :: st_sol
    ! -- local
    integer(I4) :: i, k
    !-------------------------------------------------------------------------------------------
    nreg_num = ncalc + neib_ncalc

    allocate(st_sol%head_new(nreg_num), st_sol%head_pre(nreg_num), st_sol%head_change(nreg_num))
    allocate(st_sol%srat_new(nreg_num), st_sol%rel_perm(nreg_num))
    allocate(st_sol%stor_new(nreg_num))
    allocate(st_hydr%abyd_conn(tconn_num), st_hydr%hydf_conn(tconn_num))
    allocate(st_hydr%abyd_seal(st_bcnd%seal_num), st_hydr%hydf_seal(st_bcnd%seal_num))
    allocate(st_bcnd%seal2calc(st_bcnd%seal_num), st_bcnd%seal2seal(st_bcnd%seal_num))
    allocate(dir_conn(tconn_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, nreg_num
      st_sol%head_new(i) = DZERO ; st_sol%head_pre(i) = DZERO ; st_sol%head_change(i) = DZERO
      st_sol%srat_new(i) = DZERO ; st_sol%rel_perm(i) = DZERO ; dir_conn(i) = 0
      st_sol%stor_new(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, tconn_num
      st_hydr%abyd_conn(i) = DZERO ; st_hydr%hydf_conn(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, st_bcnd%seal_num
      st_hydr%abyd_seal(i) = DZERO ; st_hydr%hydf_seal(i) = DZERO
      st_hydr%abyd_seal(i) = DZERO ; st_hydr%hydf_seal(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i, k)
    do i = 1, nreg_num
      do k = off_index(i-1)+1, off_index(i)
        st_hydr%abyd_conn(k) = st_hydr%area_dis(k) ; st_hydr%hydf_conn(k) = st_hydr%sat_hydf(k)
        dir_conn(k) = conn_dir(k)
      end do
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, st_bcnd%seal_num
      st_hydr%abyd_seal(i) = st_hydr%sea_abyd(i) ; st_hydr%hydf_seal(i) = st_hydr%sea_hydf(i)
      st_bcnd%seal2calc(i) = st_bcnd%sea2cal(i) ; st_bcnd%seal2seal(i) = st_bcnd%sea2sea(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(st_hydr%area_dis, st_hydr%sat_hydf, conn_dir)
    deallocate(st_hydr%sea_abyd, st_hydr%sea_hydf)
    deallocate(st_bcnd%sea2cal, st_bcnd%sea2sea)

    if (st_out_type%velc > 0) then
      allocate(st_hydr%inv_dis(tconn_num))
      allocate(st_hydr%dis_seal(st_bcnd%seal_num), dir_seal(st_bcnd%seal_num))
      !$omp parallel
      !$omp do private(i, k)
      do i = 1, ncalc
        do k = off_index(i-1)+1, off_index(i)
          st_hydr%inv_dis(k) = st_hydr%conn_dis(k)
        end do
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, st_bcnd%seal_num
        st_hydr%dis_seal(i) = st_hydr%sea_dis(i) ; dir_seal(i) = sea_dir(i)
      end do
      !$omp end do
      !$omp end parallel
      deallocate(st_hydr%conn_dis)
      deallocate(st_hydr%sea_dis, sea_dir)
    end if

    allocate(crs_index(st_ctrl%amg_nlevel), array_var(st_ctrl%amg_nlevel))
    allocate(crs_index(1)%offrow(tconn_num), crs_index(1)%offind(0:nreg_num))
    allocate(left_offr(tconn_num), right_offr(tconn_num))
    allocate(array_var(1)%dmat(nreg_num), array_var(1)%lumat(tconn_num))
    allocate(array_var(1)%rhs(nreg_num))
    !$omp parallel
    !$omp do private(i)
    do i = 1, tconn_num
      crs_index(1)%offrow(i) = 0
      crs_index(1)%offrow(i) = off_row(i)
      left_offr(i) = 0
      left_offr(i) = left_off(i)
      right_offr(i) = 0
      right_offr(i) = right_off(i)
      array_var(1)%lumat(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, nreg_num
      crs_index(1)%offind(i) = 0
      crs_index(1)%offind(i) = off_index(i)
      array_var(1)%dmat(i) = DZERO
      array_var(1)%rhs(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel
    crs_index(1)%unknow = ncalc ; crs_index(1)%lunum = tconn_num
    crs_index(1)%offind(0) = off_index(0)

    if (st_ctrl%precon_type == 1) then
      allocate(pro_var(st_ctrl%amg_nlevel), res_var(st_ctrl%amg_nlevel))
      allocate(array_var(1)%x(nreg_num))
      !$omp parallel do private(i)
      do i = 1, nreg_num
        array_var(1)%x(i) = DZERO
      end do
      !$omp end parallel do
    end if

    deallocate(off_row, off_index)

  end subroutine allocate_matvec

end module allocate_solution
