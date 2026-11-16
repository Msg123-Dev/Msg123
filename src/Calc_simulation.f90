module calc_simulation
  ! -- modules
  use kind_module, only: I4, DP

  implicit none
  private
  public :: calc_l2norm2, calc_resi

  ! -- local

  contains

  subroutine calc_l2norm2(l2_num, inbr, l2norm2)
  !*********************************************************************************************
  ! calc_l2norm2 -- Calculate l2 norm square
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DZERO

    ! -- inout
    integer(I4), intent(in) :: l2_num
    real(DP), intent(in) :: inbr(:)
    real(DP), intent(out) :: l2norm2
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------

    l2norm2 = DZERO
    !$omp parallel do private(i) reduction(+:l2norm2)
    do i = 1, l2_num
      l2norm2 = l2norm2 + inbr(i)*inbr(i)
    end do
    !$omp end parallel do

  end subroutine calc_l2norm2

  subroutine calc_resi(nres, offind, offrow, estd, estlu, estx, curb, resb)
  !*********************************************************************************************
  ! calc_resi -- Calculate residual
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: nres
    integer(I4), intent(in) :: offind(0:), offrow(:)
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
      off_sta = offind(i-1) + 1
      off_end = offind(i)
      do k = off_sta, off_end
        j = offrow(k)
        axsum = axsum + estlu(k)*estx(j)
      end do
      resb(i) = curb(i) - axsum
    end do
    !$omp end parallel do

  end subroutine calc_resi

end module calc_simulation
