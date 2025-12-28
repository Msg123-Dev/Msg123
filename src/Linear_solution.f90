module linear_solution
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE, DINFI
  use initial_module, only: precon_type, maxinn_iter, nlevel, pro_totn
  use allocate_solution, only: array_var, crs_index, dir_conn
  use prep_calculation, only: conv_flag
  use check_simulation, only: check_insol
  use calc_simulation, only: calc_l2norm2, calc_resi
#ifdef MPI_MSG
  use mpi_utility, only: mpisum_val
  use mpi_solve, only: senrec_rvectv, precon_mpi_dilu, solve_mpi_ilu
#endif

  implicit none
  private
  public :: solve_linalg
  integer(I4), public :: in_iter

  ! -- local
  real(DP) :: bnorm, rnorm
  real(DP), allocatable :: resi(:)

  contains

  subroutine solve_linalg(init_norm, inx, last_norm)
  !***************************************************************************************
  ! solve_lnralg -- Solve linear algebra
  !***************************************************************************************
    ! -- module
    use allocate_solution, only: nreg_num
    use prep_calculation, only: form_switch
    ! -- inout
    real(DP), intent(in) :: init_norm
    real(DP), intent(inout) :: inx(:)
    real(DP), intent(out) :: last_norm
    ! -- local

    !-------------------------------------------------------------------------------------
    allocate(resi(nreg_num))
    !$omp parallel workshare
    resi(:) = DZERO
    !$omp end parallel workshare

    bnorm = init_norm
    if (form_switch == 0) then
    ! -- Solve Preconditioned Conjugate Gradient method (pcg)
      call solve_pcg(1, inx)
    else if (form_switch == 1) then
    ! -- Solve Preconditioned Bi-Conjugate Gradient Stabilized method (bicgs)
      call solve_bicgs(1, inx)
    end if

    last_norm = rnorm

    deallocate(resi)

  end subroutine solve_linalg

  subroutine solve_pcg(level, inx)
  !***************************************************************************************
  ! solve_pcg -- Solve Preconditioned Conjugate Gradient method
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: level
    real(DP), intent(inout) :: inx(:)
    ! -- local
    integer(I4) :: n, d_size, lu_size, reg_size
    real(DP) :: sk, sk0, sk2, alpha, beta
    real(DP), allocatable :: d(:), z(:), p(:), q(:)
    real(DP), allocatable :: level_d(:), level_lu(:), level_b(:)
#ifdef MPI_MSG
    real(DP) :: sum_sk, sum_rnorm
#endif
    !-------------------------------------------------------------------------------------
    d_size = crs_index(level)%unknow
    lu_size = crs_index(level)%lunum
    reg_size = size(inx)

    sk = DZERO ; sk0 = DINFI ; sk2 = DZERO ; alpha = DZERO ; beta = DZERO
    allocate(d(reg_size), z(reg_size), p(reg_size), q(d_size))
    allocate(level_d(reg_size), level_lu(lu_size), level_b(reg_size))
    !$omp parallel workshare
    d(:) = DZERO ; z(:) = DZERO ; p(:) = DZERO ; q(:) = DZERO
    d(:) = array_var(level)%dmat(:)
    level_d(:) = array_var(level)%dmat(:)
    level_lu(:) = array_var(level)%lumat(:)
    level_b(:) = array_var(level)%rhs(:)
    !$omp end parallel workshare

    ! -- Calculate residual (resi)
      call calc_resi(level, d_size, level_d, level_lu, inx, level_b, resi)

!    if (precon_type == 0) then
#ifdef MPI_MSG
      if (pro_totn /= 1 .and. level == 1) then
        ! -- Preconditon mpi incomplete lu diagonal (mpi_dilu)
          call precon_mpi_dilu(level_d, level_lu, d)
      else
!        ! -- Preconditon diagonal scaling (dscal)
!          call precon_dscal(d_size, level_d, d)
        ! -- Preconditon incomplete lu diagonal (dilu)
          call precon_dilu(level, d_size, level_d, level_lu, d)
      end if
#else
!      ! -- Preconditon diagonal scaling (dscal)
!        call precon_dscal(d_size, level_d, d)
      ! -- Preconditon incomplete lu diagonal (dilu)
        call precon_dilu(level, d_size, level_d, level_lu, d)
#endif
!    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(resi)
        call senrec_rvectv(d)
    end if
#endif

    pcg_inter: do in_iter = 1, maxinn_iter
      if (precon_type == 0) then
#ifdef MPI_MSG
        if (pro_totn /= 1 .and. level == 1) then
          ! -- Solve mpi ilu factorization (ilu)
            call solve_mpi_ilu(resi, d, level_lu, z)
        else
!          ! -- Solve preconditioned diagonal system (dsys)
!            call solve_dsys(d_size, resi, d, z)
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(level, d_size, resi, d, level_lu, z)
        end if
#else
!        ! -- Solve preconditioned diagonal system (dsys)
!          call solve_dsys(d_size, resi, d, z)
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(level, d_size, resi, d, level_lu, z)
#endif
      else if (precon_type == 1) then
        ! -- Loop cycle for amg (amg)
          call loop_amg(level, resi, d, z)
      end if

      sk = DZERO
      !$omp parallel workshare
      sk = dot_product(resi(1:d_size), z(1:d_size))
      !$omp end parallel workshare

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(sk, "sk l2-norm", sum_sk)
        sk = sum_sk
      end if
#endif

      if (in_iter == 1) then
        !$omp parallel workshare
        p(1:d_size) = z(1:d_size)
        !$omp end parallel workshare
      else
        beta = sk/sk0
        !$omp parallel do private(n)
        do n = 1, d_size
          p(n) = z(n) + beta*p(n)
        end do
        !$omp end parallel do
      end if

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(p)
      end if
#endif
      ! -- Calculate matrix-vector multiplication (matvec)
        call calc_matvec(level, d_size, p, level_d, level_lu, q)

      sk2 = DZERO
      !$omp parallel workshare
      sk2 = dot_product(p(1:d_size), q(1:d_size))
      !$omp end parallel workshare

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(sk2, "sk2 l2-norm", sum_sk)
        sk2 = sum_sk
      end if
#endif

      alpha = sk/sk2
      !$omp parallel do private(n)
      do n = 1, d_size
        inx(n) = inx(n) + alpha*p(n)
        resi(n) = resi(n) - alpha*q(n)
      end do
      !$omp end parallel do

      ! -- Calculate l2norm square (l2norm2)
        call calc_l2norm2(level, resi, rnorm)
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(rnorm, "residual l2-norm", sum_rnorm)
        rnorm = sum_rnorm
      end if
#endif
      ! -- Check inner solution (check_insol)
        call check_insol(bnorm, rnorm)

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(resi)
          call senrec_rvectv(inx)
      end if
#endif

      if (conv_flag == 1 .or. in_iter == maxinn_iter) then
        exit pcg_inter
      else
        sk0 = sk
      end if

    end do pcg_inter

    deallocate(d, z, p, q, level_d, level_lu, level_b)

  end subroutine solve_pcg

  subroutine solve_bicgs(level, inx)
  !***************************************************************************************
  ! solve_bicgs -- Solve Preconditioned Bi-Conjugate Gradient Stabilized method
  !***************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: level
    real(DP), intent(inout) :: inx(:)
    ! -- local
    integer(I4) :: n, d_size, lu_size, reg_size
    real(DP) :: sk, sk0, sk2, alpha, beta, bicgs_omega, ts, tt
    real(DP), allocatable :: d(:), z(:), p(:), v(:), rs(:), t(:)
    real(DP), allocatable :: level_d(:), level_lu(:), level_b(:)
#ifdef MPI_MSG
    real(DP) :: sum_sk, sum_rnorm
#endif
    !-------------------------------------------------------------------------------------
    d_size = crs_index(level)%unknow
    lu_size = crs_index(level)%lunum
    reg_size = size(inx)

    sk = DZERO ; sk0 = DINFI ; sk2 = DZERO ; alpha = DZERO ; beta = DZERO
    bicgs_omega = DONE ; ts = DZERO ; tt = DZERO
    allocate(d(reg_size), z(reg_size), p(reg_size), v(d_size), rs(d_size), t(d_size))
    allocate(level_d(reg_size), level_lu(lu_size), level_b(reg_size))
    !$omp parallel workshare
    z(:) = DZERO ; p(:) = DZERO ; v(:) = DZERO ; rs(:) = DZERO ; t(:) = DZERO
    d(:) = array_var(level)%dmat(:)
    level_d(:) = array_var(level)%dmat(:)
    level_lu(:) = array_var(level)%lumat(:)
    level_b(:) = array_var(level)%rhs(:)
    !$omp end parallel workshare

    ! -- Calculate residual (resi)
      call calc_resi(level, d_size, level_d, level_lu, inx, level_b, resi)

!    if (precon_type == 0) then
#ifdef MPI_MSG
      if (pro_totn /= 1 .and. level == 1) then
        ! -- Preconditon mpi incomplete lu diagonal (mpi_dilu)
          call precon_mpi_dilu(level_d, level_lu, d)
      else
!        ! -- Preconditon diagonal scaling (dscal)
!          call precon_dscal(d_size, level_d, d)
        ! -- Preconditon incomplete lu diagonal (dilu)
          call precon_dilu(level, d_size, level_d, level_lu, d)
      end if
#else
!      ! -- Preconditon diagonal scaling (dscal)
!        call precon_dscal(d_size, level_d, d)
      ! -- Preconditon incomplete lu diagonal (dilu)
        call precon_dilu(level, d_size, level_d, level_lu, d)
#endif
!    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(d)
    end if
#endif

    !$omp parallel workshare
    rs(1:d_size) = resi(1:d_size)
    !$omp end parallel workshare

    bicg_inter: do in_iter = 1, maxinn_iter
      sk = DZERO
      !$omp parallel workshare
      sk = dot_product(resi(1:d_size), rs(1:d_size))
      !$omp end parallel workshare
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(sk, "sk l2-norm", sum_sk)
        sk = sum_sk
      end if
#endif

      if (in_iter == 1) then
        !$omp parallel workshare
        p(1:d_size) = resi(1:d_size)
        !$omp end parallel workshare
      else
        beta = (sk/sk0)*(alpha/bicgs_omega)
        !$omp parallel do private(n)
        do n = 1, d_size
          p(n) = resi(n) + beta*(p(n)-bicgs_omega*v(n))
        end do
        !$omp end parallel do
      end if
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(p)
      end if
#endif

      if (precon_type == 0) then
#ifdef MPI_MSG
        if (pro_totn /= 1 .and. level == 1) then
          ! -- Solve mpi ilu factorization (ilu)
            call solve_mpi_ilu(p, d, level_lu, z)
        else
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(level, d_size, p, d, level_lu, z)
        end if
#else
!        ! -- Solve preconditioned diagonal system (dsys)
!          call solve_dsys(d_size, p, d, z)
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(level, d_size, p, d, level_lu, z)
#endif
      else if (precon_type == 1) then
        ! -- Loop cycle for amg (amg)
          call loop_amg(level, p, d, z)
      end if

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(z)
      end if
#endif

      ! -- Calculate matrix-vector multiplication (matvec)
        call calc_matvec(level, d_size, z, level_d, level_lu, v)

      sk2 = DZERO
      !$omp parallel workshare
      sk2 = dot_product(v(1:d_size), rs(1:d_size))
      !$omp end parallel workshare
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(sk2, "sk2 l2-norm", sum_sk)
        sk2 = sum_sk
      end if
#endif
      alpha = sk/sk2

      !$omp parallel do private(n)
      do n = 1, d_size
        inx(n) = inx(n) + alpha*z(n)
        resi(n) = resi(n) - alpha*v(n)
      end do
      !$omp end parallel do

      ! -- Calculate l2norm square (l2norm2)
        call calc_l2norm2(level, resi, rnorm)
#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(rnorm, "residual l2-norm", sum_rnorm)
        rnorm = sum_rnorm
      end if
#endif
      ! -- Check inner solution (check_insol)
        call check_insol(bnorm, rnorm)

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(resi)
          call senrec_rvectv(inx)
      end if
#endif

      if (conv_flag == 1) then
        exit bicg_inter
      else if (conv_flag /= 1) then
        if (precon_type == 0) then
#ifdef MPI_MSG
          if (pro_totn /= 1 .and. level == 1) then
            ! -- Solve mpi ilu factorization (ilu)
              call solve_mpi_ilu(resi, d, level_lu, z)
          else
            ! -- Solve ilu factorization (ilu)
              call solve_ilu(level, d_size, resi, d, level_lu, z)
          end if
#else
!        ! -- Solve preconditioned diagonal system (dsys)
!          call solve_dsys(d_size, resi, d, z)
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(level, d_size, resi, d, level_lu, z)
#endif
        else if (precon_type == 1) then
        ! -- Loop cycle for amg (amg)
          call loop_amg(level, resi, d, z)
        end if

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(z)
        end if
#endif

        ! -- Calculate matrix-vector multiplication (matvec)
          call calc_matvec(level, d_size, z, level_d, level_lu, t)

        ts = DZERO ; tt = DZERO
        !$omp parallel workshare
        ts = dot_product(t(1:d_size), resi(1:d_size))
        tt = dot_product(t(1:d_size), t(1:d_size))
        !$omp end parallel workshare

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Sum value for MPI (val)
            call mpisum_val(ts, "ts l2-norm", sum_sk)
          ts = sum_sk
          ! -- Sum value for MPI (val)
            call mpisum_val(tt, "tt l2-norm", sum_sk)
          tt = sum_sk
        end if
#endif

        bicgs_omega = ts/tt

        !$omp parallel do private(n)
        do n = 1, d_size
          inx(n) = inx(n) + bicgs_omega*z(n)
          resi(n) = resi(n) - bicgs_omega*t(n)
        end do
        !$omp end parallel do

        ! -- Calculate l2norm square (l2norm2)
          call calc_l2norm2(level, resi, rnorm)
#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Sum value for MPI (val)
            call mpisum_val(rnorm, "residual l2-norm", sum_rnorm)
          rnorm = sum_rnorm
        end if
#endif
        ! -- Check inner solution (check_insol)
          call check_insol(bnorm, rnorm)

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(resi)
            call senrec_rvectv(inx)
        end if
#endif

        if (conv_flag == 1 .or. in_iter == maxinn_iter) then
          exit bicg_inter
        else
          sk0 = sk
        end if
      end if

    end do bicg_inter

    deallocate(d, z, p, v, rs, t, level_d, level_lu, level_b)

  end subroutine solve_bicgs

  subroutine loop_amg(llevel, r, d, e)
  !***************************************************************************************
  ! loop_amg -- Loop cycle for amg
  !***************************************************************************************
    ! -- module
!    use initial_module, only: maxvcy_iter, max_sweep
    use initial_module, only: maxvcy_iter
    use allocate_solution, only: pro_var, res_var
    ! -- inout
    integer(I4), intent(in) :: llevel
    real(DP), intent(in) :: r(:), d(:)
    real(DP), intent(out) :: e(:)
    ! -- local
    integer(I4) :: i, j, k, d_size, lu_size, reg_size
    integer(I4) :: mgd_size, mgreg_size
    integer(I4) :: v_iter, vlevel, ncoa, nfin
!    integer(I4) :: sweep
    integer(I4) :: rst, ren, pst, pen
    real(DP), allocatable :: temp_d(:), temp_lu(:), temp_r(:), temp_x(:), temp_b(:)
    real(DP), allocatable :: trhs(:), tx(:)
    !-------------------------------------------------------------------------------------
    mgd_size = crs_index(llevel)%unknow ; mgreg_size = size(array_var(llevel)%x)
    !$omp parallel
    !$omp do private(i)
    do i = 1, mgd_size
      array_var(llevel)%x(i) = DZERO
      array_var(llevel)%dmat(i) = d(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, mgreg_size
      array_var(llevel)%rhs(i) = r(i)
    end do
    !$omp end do
    !$omp end parallel

    do v_iter = 1, maxvcy_iter  ! V-cycle loop
      do vlevel = llevel, nlevel-1

        d_size = crs_index(vlevel)%unknow
        lu_size = crs_index(vlevel)%lunum
        reg_size = size(array_var(vlevel)%x)

        allocate(temp_d(reg_size), temp_lu(lu_size), temp_x(reg_size), temp_b(reg_size))
        allocate(temp_r(reg_size))
        !$omp parallel workshare
        temp_d(:) = array_var(vlevel)%dmat(:) ; temp_lu(:) = array_var(vlevel)%lumat(:)
        temp_x(:) = array_var(vlevel)%x(:) ; temp_b(:) = array_var(vlevel)%rhs(:)
        temp_r(:) = DZERO
        !$omp end parallel workshare

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(temp_x)
            call senrec_rvectv(temp_b)
        end if
        if (pro_totn /= 1 .and. vlevel == 1) then
          ! -- Solve mpi ilu factorization (ilu)
            call solve_mpi_ilu(temp_b, temp_d, temp_lu, temp_x)
        else
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(vlevel, d_size, temp_b, temp_d, temp_lu, temp_x)
        end if
#else
!        ! -- Solve preconditioned diagonal system (dsys)
!          call solve_dsys(d_size, temp_b, temp_d, temp_x)
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(vlevel, d_size, temp_b, temp_d, temp_lu, temp_x)
!        do sweep = 1, max_sweep
!          call smooth_gs(vlevel, d_size, temp_d, temp_lu, temp_b, temp_x)
!        end do
#endif

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(temp_x)
        end if
#endif

        ! -- Calculate residual
          call calc_resi(vlevel, d_size, temp_d, temp_lu, temp_x, temp_b, temp_r)

        deallocate(temp_d, temp_lu, temp_x, temp_b)

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(temp_r)
        end if
#endif

        ncoa = crs_index(vlevel+1)%unknow
        allocate(trhs(ncoa))
        !$omp parallel
        !$omp workshare
        trhs(:) = DZERO
        !$omp end workshare
        !$omp do private(i, j, k, rst, ren)
        do i = 1, ncoa
          rst = res_var(vlevel+1)%rindex(i-1) + 1
          ren = res_var(vlevel+1)%rindex(i)
          do k = rst, ren
            j = res_var(vlevel+1)%roffrow(k)
            trhs(i) = trhs(i) + res_var(vlevel+1)%rval(k)*temp_r(j)
          end do
          array_var(vlevel+1)%rhs(i) = trhs(i)
        end do
        !$omp end do

        !$omp do private(i)
        do i = 1, ncoa
          array_var(vlevel+1)%x(i) = DZERO
        end do
        !$omp end do
        !$omp end parallel

        deallocate(temp_r, trhs)

      end do

      d_size = crs_index(nlevel)%unknow
      lu_size = crs_index(nlevel)%lunum
      reg_size = size(array_var(nlevel)%x)

      allocate(temp_d(d_size), temp_lu(lu_size), temp_x(reg_size), temp_b(reg_size))
      !$omp parallel workshare
      temp_d(:) = array_var(nlevel)%dmat(:) ; temp_lu(:) = array_var(nlevel)%lumat(:)
      temp_x(:) = array_var(nlevel)%x(:) ; temp_b(:) = array_var(nlevel)%rhs(:)
      !$omp end parallel workshare

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(temp_x)
          call senrec_rvectv(temp_b)
      end if
      if (pro_totn /= 1 .and. nlevel == 1) then
        ! -- Solve mpi ilu factorization (ilu)
          call solve_mpi_ilu(temp_b, temp_d, temp_lu, temp_x)
      else
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(nlevel, d_size, temp_b, temp_d, temp_lu, temp_x)
      end if
#else
!      ! -- Solve preconditioned diagonal system (dsys)
!        call solve_dsys(d_size, temp_b, temp_d, temp_x)
      ! -- Solve ilu factorization (ilu)
        call solve_ilu(nlevel, d_size, temp_b, temp_d, temp_lu, temp_x)
!      do sweep = 1, max_sweep
!        call smooth_gs(nlevel, d_size, temp_d, temp_lu, temp_b, temp_x)
!      end do
#endif

#ifdef MPI_MSG
      if (pro_totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(temp_x)
      end if
#endif
      !$omp parallel workshare
      array_var(nlevel)%x(:) = temp_x(:)
      !$omp end parallel workshare

      deallocate(temp_d, temp_lu, temp_x, temp_b)

      do vlevel = nlevel-1, llevel, -1
        nfin = crs_index(vlevel+1)%unknow
        allocate(tx(nfin))
        !$omp parallel
        !$omp workshare
        tx(:) = DZERO
        !$omp end workshare
        !$omp do private(i, j, k, pst, pen)
        do i = 1, nfin
          pst = pro_var(vlevel+1)%pindex(i-1) + 1
          pen = pro_var(vlevel+1)%pindex(i)
          do k = pst, pen
            j = pro_var(vlevel+1)%poffrow(k)
            tx(i) = tx(i) + pro_var(vlevel+1)%pval(k)*array_var(vlevel+1)%x(j)
          end do
          array_var(vlevel)%x(i) = array_var(vlevel)%x(i) + tx(i)
        end do
        !$omp end do
        !$omp end parallel

        deallocate(tx)

        d_size = crs_index(vlevel)%unknow
        lu_size = crs_index(vlevel)%lunum
        reg_size = size(array_var(vlevel)%x)

        allocate(temp_d(reg_size), temp_lu(lu_size), temp_x(reg_size), temp_b(reg_size))
        !$omp parallel workshare
        temp_d(:) = array_var(vlevel)%dmat(:) ; temp_lu(:) = array_var(vlevel)%lumat(:)
        temp_x(:) = array_var(vlevel)%x(:) ; temp_b(:) = array_var(vlevel)%rhs(:)
        !$omp end parallel workshare

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(temp_x)
            call senrec_rvectv(temp_b)
        end if
        if (pro_totn /= 1 .and. vlevel == 1) then
          ! -- Solve mpi ilu factorization (ilu)
            call solve_mpi_ilu(temp_b, temp_d, temp_lu, temp_x)
        else
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(vlevel, d_size, temp_b, temp_d, temp_lu, temp_x)
        end if
#else
!        ! -- Solve preconditioned diagonal system (dsys)
!          call solve_dsys(d_size, temp_b, temp_d, temp_x)
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(vlevel, d_size, temp_b, temp_d, temp_lu, temp_x)
!        do sweep = 1, max_sweep
!          call smooth_gs(vlevel, d_size, temp_d, temp_lu, temp_b, temp_x)
!        end do
#endif

        !$omp parallel workshare
        array_var(vlevel)%x(:) = temp_x(:) ; array_var(vlevel)%rhs(:) = temp_b(:)
        !$omp end parallel workshare

#ifdef MPI_MSG
        if (pro_totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(array_var(vlevel)%x)
            call senrec_rvectv(array_var(vlevel)%rhs)
        end if
#endif
        deallocate(temp_d, temp_lu, temp_x, temp_b)

      end do

    end do

    !$omp parallel do private(i)
    do i = 1, mgd_size
      e(i) = array_var(llevel)%x(i)
    end do
    !$omp end parallel do

  end subroutine loop_amg

  subroutine precon_dilu(plevel, npre, pre_ind, pre_inlu, pre_d)
  !***************************************************************************************
  ! precon_dilu -- Preconditon incomplete lu diagonal
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: plevel, npre
    real(DP), intent(in) :: pre_ind(:), pre_inlu(:)
    real(DP), intent(inout) :: pre_d(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: off_sta, off_end, off_sta2, off_end2
    integer(I4) :: offr, offr2
    real(DP) :: d_invk
    !-------------------------------------------------------------------------------------
    !$omp parallel workshare
    pre_d(1:npre) = pre_ind(1:npre)
    !$omp end parallel workshare
    do i = 1, npre
      off_sta = crs_index(plevel)%offind(i-1) + 1
      off_end = crs_index(plevel)%offind(i)
      do k = off_sta, off_end
        offr = crs_index(plevel)%offrow(k)
        if (3 >= dir_conn(k)) then
          d_invk = DONE/pre_d(offr)
          off_sta2 = crs_index(plevel)%offind(offr-1) + 1
          off_end2 = crs_index(plevel)%offind(offr)
          do j = off_sta2, off_end2
            offr2 = crs_index(plevel)%offrow(j)
            if (3 < dir_conn(j)) then
              if (offr2 == i) then
                pre_d(i) = pre_d(i) - pre_inlu(k)*d_invk*pre_inlu(j)
              end if
            end if
          end do
        end if
      end do
    end do

  end subroutine precon_dilu

!  subroutine precon_dscal(npre, pre_ind, pre_d)
!  !***************************************************************************************
!  ! precon_dscal -- Preconditon diagonal scaling
!  !***************************************************************************************
!    ! -- modules
!
!    ! -- inout
!    integer(I4), intent(in) :: npre
!    real(DP), intent(in) :: pre_ind(:)
!    real(DP), intent(out) :: pre_d(:)
!    ! -- local
!
!    !-------------------------------------------------------------------------------------
!    !$omp parallel workshare
!    pre_d(1:npre) = DONE/pre_ind(1:npre)
!    !$omp end parallel workshare
!
!  end subroutine precon_dscal

  subroutine solve_ilu(plevel, npre, inrhs, indmat, inlumat, outx)
  !***************************************************************************************
  ! solve_ilu -- Solve ilu factorization
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: plevel, npre
    real(DP), intent(in) :: inrhs(:), indmat(:), inlumat(:)
    real(DP), intent(out) :: outx(:)
    ! -- local
    integer(I4) :: i, k
    integer(I4) :: off_sta, off_end, offr
    real(DP), allocatable :: temp_outx(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_outx(npre))
    !$omp parallel workshare
    outx(1:npre) = inrhs(1:npre) ; temp_outx(:) = DZERO
    !$omp end parallel workshare

    ! Forward Substitution
!    !$omp parallel do private(i, k, s, off_sta, off_end, offr)
    do i = 1, npre
      temp_outx(i) = outx(i)
      off_sta = crs_index(plevel)%offind(i-1) + 1
      off_end = crs_index(plevel)%offind(i)
      do k = off_sta, off_end
        offr = crs_index(plevel)%offrow(k)
        if (3 >= dir_conn(k)) then
          temp_outx(i) = temp_outx(i) - inlumat(k)*outx(offr)
        end if
      end do
      outx(i) = temp_outx(i)/indmat(i)
    end do
!    !$omp end parallel do

    ! Backward Substitution
!    !$omp parallel do private(i, k, s, off_sta, off_end, offr)
    do i = npre, 1, -1
      temp_outx(i) = DZERO
      off_sta = crs_index(plevel)%offind(i-1) + 1
      off_end = crs_index(plevel)%offind(i)
      do k = off_sta, off_end
        offr = crs_index(plevel)%offrow(k)
        if (3 < dir_conn(k)) then
          temp_outx(i) = temp_outx(i) + inlumat(k)*outx(offr)
        end if
      end do
      outx(i) = outx(i) - temp_outx(i)/indmat(i)
    end do
!    !$omp end parallel do

  end subroutine solve_ilu

!  subroutine solve_dsys(nregpre, inrhs, indmat, outx)
!  !***************************************************************************************
!  ! solve_dsys -- Solve preconditioned diagonal system
!  !***************************************************************************************
!    ! -- modules
!
!    ! -- inout
!    integer(I4), intent(in) :: nregpre
!    real(DP), intent(in) :: inrhs(:), indmat(:)
!    real(DP), intent(out) :: outx(:)
!    ! -- local
!
!    !-------------------------------------------------------------------------------------
!    !$omp parallel workshare
!    outx(1:nregpre) = inrhs(1:nregpre)*indmat(1:nregpre)
!    !$omp end parallel workshare
!
!  end subroutine solve_dsys

!  subroutine smooth_gs(level, nsmo, ind, inlu, inb, inx)
!  !***************************************************************************************
!  ! smooth_gs -- Smooth by Gauss-seidel method
!  !***************************************************************************************
!    ! -- module
!
!    ! -- inout
!    integer(I4), intent(in) :: level, nsmo
!    real(DP), intent(in) :: ind(:), inlu(:), inb(:)
!    real(DP), intent(inout) :: inx(:)
!    ! -- local
!    integer(I4) :: i, n, k
!    integer(I4) :: ns, off_sta, off_end, offr
!    real(DP), allocatable :: temp_x(:), temp_d(:)
!    !-------------------------------------------------------------------------------------
!    if (level == nlevel) then
!      ns = maxinn_iter
!    else
!      ns = 1
!    end if
!    allocate(temp_x(nsmo), temp_d(nsmo))
!    !$omp parallel workshare
!    temp_x(:) = DZERO ; temp_d(:) = DONE/ind(1:nsmo)
!    !$omp end parallel workshare
!
!!    !$omp parallel do private(i, n, k, s, off_sta, off_end, offr)
!    do i = 1, ns
!      do n = 1, nsmo
!        off_sta = crs_index(level)%offind(n-1) + 1
!        off_end = crs_index(level)%offind(n)
!        do k = off_sta, off_end
!          offr = crs_index(level)%offrow(k)
!          temp_x(n) = temp_x(n) + inlu(k)*inx(offr)
!        end do
!        inx(n) = temp_d(n)*(inb(n)-temp_x(n))
!      end do
!    end do
!!    !$omp end parallel do
!
!  end subroutine smooth_gs

  subroutine calc_matvec(mvlevel, nmv, invec, indmat, inlumat, outvec)
  !***************************************************************************************
  ! calc_matvec -- Calculate matrix-vector multiplication
  !***************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: mvlevel, nmv
    real(DP), intent(in) :: invec(:), indmat(:), inlumat(:)
    real(DP), intent(out) :: outvec(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: off_sta, off_end
    real(DP), allocatable :: temp_vec(:)
    !-------------------------------------------------------------------------------------
    allocate(temp_vec(nmv))
    !$omp parallel
    !$omp workshare
    temp_vec(:) = DZERO
    !$omp end workshare
    !$omp do private(i, j, k, off_sta, off_end)
    do i = 1, nmv
      temp_vec(i) = indmat(i)*invec(i)
      off_sta = crs_index(mvlevel)%offind(i-1) + 1
      off_end = crs_index(mvlevel)%offind(i)
      do k = off_sta, off_end
        j = crs_index(mvlevel)%offrow(k)
        temp_vec(i) = temp_vec(i) + inlumat(k)*invec(j)
      end do
      outvec(i) = temp_vec(i)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine calc_matvec

end module linear_solution
