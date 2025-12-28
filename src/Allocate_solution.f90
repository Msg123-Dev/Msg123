module allocate_solution
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO
  use set_cell, only: ncalc

  implicit none
  private
  public :: allocate_solvar
  integer(I4), allocatable, public :: seal2calc(:), seal2seal(:)
  integer(I4), allocatable, public :: dir_conn(:), dir_seal(:)
  integer(I4), allocatable, public :: left_offr(:), right_offr(:)
  real(DP), allocatable, public :: head_old(:), srat_old(:)
  real(DP), allocatable, public :: head_new(:), head_pre(:), head_change(:)
  real(DP), allocatable, public :: srat_new(:)
  real(DP), allocatable, public :: surf_head(:), surf_old(:), surf_rati(:)
  real(DP), allocatable, public :: rel_perm(:)
  real(DP), allocatable, public :: abyd_conn(:), hydf_conn(:), inv_dis(:)
  real(DP), allocatable, public :: abyd_seal(:), hydf_seal(:), dis_seal(:)
  integer(I4), public :: nreg_num

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

  subroutine allocate_solvar()
  !***************************************************************************************
  ! allocate_solvar -- Allocate solution variable for time step
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    ! -- Allocate timeupdate (timeup)
      call allocate_timeup()

    ! -- Allocate matrix and vector (matvec)
      call allocate_matvec()

  end subroutine allocate_solvar

  subroutine allocate_timeup()
  !***************************************************************************************
  ! allocate_timeup -- Allocate timeupdate
  !***************************************************************************************
    ! -- modules
    use set_cell, only: ncals
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(head_old(ncalc), srat_old(ncalc))
    allocate(surf_head(ncals), surf_old(ncals), surf_rati(ncals))
    !$omp parallel workshare
    head_old(:) = DZERO ; srat_old(:) = DZERO
    surf_head(:) = DZERO ; surf_old(:) = DZERO ; surf_rati(:) = DZERO
    !$omp end parallel workshare

  end subroutine allocate_timeup

  subroutine allocate_matvec()
  !***************************************************************************************
  ! allocate_matvec -- Allocate matrix and vector
  !***************************************************************************************
    ! -- modules
    use initial_module, only: precon_type, st_out_type, amg_nlevel
    use set_cell, only : neib_ncalc
    use set_condition, only: tconn_num, nseal, off_row, off_index, area_dis, sat_hydf,&
                             conn_dis, conn_dir, sea_hydf, sea_abyd, sea_dis, sea_dir,&
                             sea2cal, sea2sea, left_off, right_off
    ! -- inout

    ! -- local
    integer(I4) :: i, k
    !-------------------------------------------------------------------------------------
    nreg_num = ncalc + neib_ncalc

    allocate(head_new(nreg_num), head_pre(nreg_num), head_change(nreg_num))
    allocate(srat_new(nreg_num), rel_perm(nreg_num))
    allocate(abyd_conn(tconn_num), hydf_conn(tconn_num))
    allocate(abyd_seal(nseal), hydf_seal(nseal))
    allocate(seal2calc(nseal), seal2seal(nseal))
    allocate(dir_conn(tconn_num))
    !$omp parallel
    !$omp workshare
    head_new(:) = DZERO ; head_pre(:) = DZERO ; head_change(:) = DZERO
    srat_new(:) = DZERO ; rel_perm(:) = DZERO
    abyd_conn(:) = DZERO ; hydf_conn(:) = DZERO
    abyd_seal(:) = DZERO ; hydf_seal(:) = DZERO
    seal2calc(:) = 0 ; seal2seal(:) = 0 ; dir_conn(:) = 0
    !$omp end workshare

    !$omp do private(i, k)
    do i = 1, nreg_num
      do k = off_index(i-1)+1, off_index(i)
        abyd_conn(k) = area_dis(k) ; hydf_conn(k) = sat_hydf(k)
        dir_conn(k) = conn_dir(k)
      end do
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, nseal
      abyd_seal(i) = sea_abyd(i) ; hydf_seal(i) = sea_hydf(i)
      seal2calc(i) = sea2cal(i) ; seal2seal(i) = sea2sea(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(area_dis, sat_hydf, conn_dir)
    deallocate(sea_abyd, sea_hydf)

    if (st_out_type%velc > 0) then
      allocate(inv_dis(tconn_num))
      allocate(dis_seal(nseal), dir_seal(nseal))
      !$omp parallel
      !$omp do private(i, k)
      do i = 1, ncalc
        do k = off_index(i-1)+1, off_index(i)
          inv_dis(k) = conn_dis(k)
        end do
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, nseal
        dis_seal(i) = sea_dis(i) ; dir_seal(i) = sea_dir(i)
      end do
      !$omp end do
      !$omp end parallel
      deallocate(conn_dis)
      deallocate(sea_dis, sea_dir)
    end if

    allocate(crs_index(amg_nlevel), array_var(amg_nlevel))
    allocate(crs_index(1)%offrow(tconn_num), crs_index(1)%offind(0:nreg_num))
    allocate(left_offr(tconn_num), right_offr(tconn_num))
    !$omp parallel workshare
    crs_index(1)%offrow(:) = 0 ; crs_index(1)%offind(:) = 0
    crs_index(1)%offrow(:) = off_row(1:tconn_num)
    crs_index(1)%offind(:) = off_index(:)
    left_offr(:) = 0 ; right_offr(:) = 0
    left_offr(:) = left_off(1:tconn_num) ; right_offr(:) = right_off(1:tconn_num)
    !$omp end parallel workshare
    crs_index(1)%unknow = ncalc ; crs_index(1)%lunum = tconn_num

    allocate(array_var(1)%dmat(nreg_num), array_var(1)%lumat(tconn_num))
    allocate(array_var(1)%rhs(nreg_num))
    !$omp parallel workshare
    array_var(1)%dmat(:) = DZERO ; array_var(1)%lumat(:) = DZERO
    array_var(1)%rhs(:) = DZERO
    !$omp end parallel workshare

    if (precon_type == 1) then
      allocate(pro_var(amg_nlevel), res_var(amg_nlevel))
      allocate(array_var(1)%x(nreg_num))
      !$omp parallel workshare
      array_var(1)%x(:) = DZERO
      !$omp end parallel workshare
    end if

    deallocate(off_row, off_index)

  end subroutine allocate_matvec

end module allocate_solution
