module calc_simulation
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO
  use allocate_solution, only: crs_index

  implicit none
  private
  public :: calc_l2norm2, calc_resi, calc_resl2norm

  ! -- local

  contains

  subroutine calc_l2norm2(l2level, inbr, l2norm2)
  !***************************************************************************************
  ! calc_l2norm2 -- Calculate l2 norm square
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: l2level
    real(DP), intent(in) :: inbr(:)
    real(DP), intent(out) :: l2norm2
    ! -- local
    integer(I4) :: l2_num
    !-------------------------------------------------------------------------------------
    l2_num = crs_index(l2level)%unknow

    !$omp parallel workshare
    l2norm2 = dot_product(inbr(1:l2_num), inbr(1:l2_num))
    !$omp end parallel workshare

  end subroutine calc_l2norm2

  subroutine calc_resi(rlevel, nres, estd, estlu, estx, curb, resb)
  !***************************************************************************************
  ! calc_resi -- Calculate residual
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: rlevel, nres
    real(DP), intent(in) :: estd(:), estlu(:), estx(:), curb(:)
    real(DP), intent(out) :: resb(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: off_sta, off_end
    real(DP), allocatable :: axsum(:)
    !-------------------------------------------------------------------------------------
    allocate(axsum(nres))
    !$omp parallel
    !$omp workshare
    axsum(:) = DZERO
    !$omp end workshare

    !$omp do private(i, j, k, off_sta, off_end)
    do i = 1, nres
      axsum(i) = estd(i)*estx(i)
      off_sta = crs_index(rlevel)%offind(i-1) + 1
      off_end = crs_index(rlevel)%offind(i)
      do k = off_sta, off_end
        j = crs_index(rlevel)%offrow(k)
        axsum(i) = axsum(i) + estlu(k)*estx(j)
      end do
      resb(i) = curb(i) - axsum(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(axsum)

  end subroutine calc_resi

  subroutine calc_resl2norm(rlevel, res2d, res2lu, res2x, res2b, resl2)
  !***************************************************************************************
  ! calc_resl2norm -- Calculate residual and l2 norm square
  !***************************************************************************************
    ! -- modules
    use constval_module, only: DTWO
    ! -- inout
    integer(I4), intent(in) :: rlevel
    real(DP), intent(in) :: res2d(:), res2lu(:), res2x(:), res2b(:)
    real(DP), intent(out) :: resl2
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: nres2, off_sta, off_end
    real(DP), allocatable :: res2sum(:), temp_l2(:)
    !-------------------------------------------------------------------------------------
    nres2 = crs_index(rlevel)%unknow
    allocate(res2sum(nres2), temp_l2(nres2))
    !$omp parallel
    !$omp workshare
    res2sum(:) = DZERO ; temp_l2(:) = DZERO
    !$omp end workshare

    !$omp do private(i, j, k, off_sta, off_end)
    do i = 1, nres2
      res2sum(i) = res2d(i)*res2x(i)
      off_sta = crs_index(rlevel)%offind(i-1) + 1
      off_end = crs_index(rlevel)%offind(i)
      do k = off_sta, off_end
        j = crs_index(rlevel)%offrow(k)
        res2sum(i) = res2sum(i) + res2lu(k)*res2x(j)
      end do
      temp_l2(i) = (res2b(i)-res2sum(i))**DTWO
    end do
    !$omp end do
    !$omp workshare
    resl2 = sum(temp_l2)
    !$omp end workshare
    !$omp end parallel

    deallocate(res2sum, temp_l2)

  end subroutine calc_resl2norm

end module calc_simulation
