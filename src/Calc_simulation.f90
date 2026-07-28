module calc_simulation
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO
  use allocate_solution, only: crs_index

  implicit none
  private
  public :: calc_l2norm2, calc_resi

  ! -- local

  contains

  subroutine calc_l2norm2(l2level, inbr, l2norm2)
  !*********************************************************************************************
  ! calc_l2norm2 -- Calculate l2 norm square
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: l2level
    real(DP), intent(in) :: inbr(:)
    real(DP), intent(out) :: l2norm2
    ! -- local
    integer(I4) :: i, l2_num
    !-------------------------------------------------------------------------------------------
    l2_num = crs_index(l2level)%unknow

    l2norm2 = DZERO
    !$omp parallel do private(i) reduction(+:l2norm2)
    do i = 1, l2_num
      l2norm2 = l2norm2 + inbr(i)*inbr(i)
    end do
    !$omp end parallel do

  end subroutine calc_l2norm2

  subroutine calc_resi(rlevel, nres, estd, estlu, estx, curb, resb)
  !*********************************************************************************************
  ! calc_resi -- Calculate residual
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: rlevel, nres
    real(DP), intent(in) :: estd(:), estlu(:), estx(:), curb(:)
    real(DP), intent(out) :: resb(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: off_sta, off_end
    real(DP) :: axsum
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, off_sta, off_end, axsum)
    do i = 1, nres
      axsum = estd(i)*estx(i)
      off_sta = crs_index(rlevel)%offind(i-1) + 1
      off_end = crs_index(rlevel)%offind(i)
      do k = off_sta, off_end
        j = crs_index(rlevel)%offrow(k)
        axsum = axsum + estlu(k)*estx(j)
      end do
      resb(i) = curb(i) - axsum
    end do
    !$omp end parallel do

  end subroutine calc_resi

end module calc_simulation
