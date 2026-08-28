module linear_solution
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DONE, DINFI
  use types_module, only: kryl_set, amgt_set
  use utility_module, only: st_mpi, dilu_shift_num
  use initial_module, only: st_ctrl
  use prep_calculation, only: st_time
  use allocate_solution, only: nreg_num, dir_conn, crs_index, array_var
  use check_simulation, only: check_insol
  use calc_simulation, only: calc_l2norm2, calc_resi
#ifdef MPI_MSG
  use mpi_utility, only: mpisum_val
  use mpi_solve, only: senrec_rvectv, precon_mpi_dilu, solve_mpi_ilu
#endif

  implicit none
  private
  public :: allocate_krylov, allocate_amgalg, solve_linalg
  integer(I4), public :: in_iter

  ! -- local
  real(DP) :: bnorm, rnorm

  contains

  subroutine allocate_krylov(st_kryl)
  !*********************************************************************************************
  ! allocate_krylov -- Allocate for krylov work vectors
  !*********************************************************************************************
    ! -- module
    ! -- inout
    type(kryl_set), intent(inout) :: st_kryl
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    allocate(st_kryl%resi(nreg_num))
    allocate(st_kryl%d(nreg_num), st_kryl%z(nreg_num), st_kryl%p(nreg_num))
    allocate(st_kryl%q(nreg_num), st_kryl%v(nreg_num))
    allocate(st_kryl%rs(nreg_num), st_kryl%t(nreg_num))

    !$omp parallel do private(i)
    do i = 1, nreg_num
      st_kryl%resi(i) = DZERO
      st_kryl%d(i) = DZERO ; st_kryl%z(i) = DZERO ; st_kryl%p(i) = DZERO
      st_kryl%q(i) = DZERO ; st_kryl%v(i) = DZERO
      st_kryl%rs(i) = DZERO ; st_kryl%t(i) = DZERO
    end do
    !$omp end parallel do

  end subroutine allocate_krylov

  subroutine allocate_amgalg(st_amgt)
  !*********************************************************************************************
  ! allocate_amgalg -- Allocate for amg algebra
  !*********************************************************************************************
    ! -- module
    use set_condition, only: tconn_num
    ! -- inout
    type(amgt_set), intent(inout) :: st_amgt
    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    allocate(st_amgt%td(nreg_num), st_amgt%tx(nreg_num))
    allocate(st_amgt%tb(nreg_num), st_amgt%tr(nreg_num))
    allocate(st_amgt%trhs(nreg_num), st_amgt%tfx(nreg_num), st_amgt%tlu(tconn_num))
    allocate(st_amgt%save_rhs(nreg_num), st_amgt%dilu_d(nreg_num))

    !$omp parallel
    !$omp do private(i)
    do i = 1, nreg_num
      st_amgt%td(i) = DZERO ; st_amgt%tx(i) = DZERO
      st_amgt%tb(i) = DZERO ; st_amgt%tr(i) = DZERO
      st_amgt%trhs(i) = DZERO ; st_amgt%tfx(i) = DZERO
      st_amgt%save_rhs(i) = DZERO ; st_amgt%dilu_d(i) = DZERO
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, tconn_num
      st_amgt%tlu(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel

  end subroutine allocate_amgalg

  subroutine solve_linalg(init_norm, inx, st_kryl, st_amgt, last_norm)
  !*********************************************************************************************
  ! solve_lnralg -- Solve linear algebra
  !*********************************************************************************************
    ! -- module
    use allocate_solution, only: nreg_num
    ! -- inout
    real(DP), intent(in) :: init_norm
    real(DP), intent(inout) :: inx(:)
    type(kryl_set), intent(inout) :: st_kryl
    type(amgt_set), intent(inout) :: st_amgt
    real(DP), intent(out) :: last_norm
    ! -- local
    integer(I4) :: n
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(n)
    do n = 1, nreg_num
      st_kryl%resi(n) = DZERO
    end do
    !$omp end parallel do

    bnorm = init_norm
    if (st_time%form_switch == 0) then
    ! -- Solve Preconditioned Conjugate Gradient method (pcg)
      call solve_pcg(1, inx, st_kryl, st_amgt)
    else if (st_time%form_switch == 1) then
    ! -- Solve Preconditioned Bi-Conjugate Gradient Stabilized method (bicgs)
      call solve_bicgs(1, inx, st_kryl, st_amgt)
    end if

    last_norm = rnorm

  end subroutine solve_linalg

  subroutine solve_pcg(level, inx, st_kryl, st_amgt)
  !*********************************************************************************************
  ! solve_pcg -- Solve Preconditioned Conjugate Gradient method
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: level
    real(DP), intent(inout) :: inx(:)
    type(kryl_set), intent(inout) :: st_kryl
    type(amgt_set), intent(inout) :: st_amgt
    ! -- local
    integer(I4) :: n, d_size, lu_size, reg_size
    real(DP) :: sk, sk0, sk2, alpha, beta
#ifdef MPI_MSG
    real(DP) :: sum_sk, sum_rnorm
#endif
    !-------------------------------------------------------------------------------------------
    d_size = crs_index(level)%unknow
    lu_size = crs_index(level)%lunum
    reg_size = size(inx)

    sk = DZERO ; sk0 = DINFI ; sk2 = DZERO ; alpha = DZERO ; beta = DZERO
    !$omp parallel
    !$omp do private(n)
    do n = 1, reg_size
      st_kryl%z(n) = DZERO ; st_kryl%p(n) = DZERO
      st_kryl%d(n) = array_var(level)%dmat(n)
    end do
    !$omp end do
    !$omp do private(n)
    do n = 1, d_size
      st_kryl%q(n) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    ! -- Calculate residual (resi)
      call calc_resi(level, d_size, array_var(level)%dmat, array_var(level)%lumat, inx,&
                     array_var(level)%rhs, st_kryl%resi)

    if (st_ctrl%precon_type == 0) then
#ifdef MPI_MSG
      if (st_mpi%totn /= 1 .and. level == 1) then
        ! -- Preconditon mpi incomplete lu diagonal (mpi_dilu)
          call precon_mpi_dilu(array_var(level)%dmat, array_var(level)%lumat, st_kryl%d)
      else
        ! -- Preconditon incomplete lu diagonal (dilu)
          call precon_dilu(level, d_size, array_var(level)%dmat, array_var(level)%lumat, st_kryl%d)
      end if
#else
      ! -- Preconditon incomplete lu diagonal (dilu)
        call precon_dilu(level, d_size, array_var(level)%dmat, array_var(level)%lumat, st_kryl%d)
#endif
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_kryl%d)
      end if
#endif
    end if

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Send and Receive real vector value (rvectv)
        call senrec_rvectv(st_kryl%resi)
    end if
#endif

    pcg_inter: do in_iter = 1, st_ctrl%maxinn_iter
      if (st_ctrl%precon_type == 0) then
#ifdef MPI_MSG
        if (st_mpi%totn /= 1 .and. level == 1) then
          ! -- Solve mpi ilu factorization (ilu)
            call solve_mpi_ilu(st_kryl%resi, st_kryl%d, array_var(level)%lumat, st_kryl%z)
        else
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(level, d_size, st_kryl%resi, st_kryl%d, array_var(level)%lumat, st_kryl%z)
        end if
#else
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(level, d_size, st_kryl%resi, st_kryl%d, array_var(level)%lumat, st_kryl%z)
#endif
      else if (st_ctrl%precon_type == 1) then
        ! -- Loop cycle for amg (amg)
          call loop_amg(level, st_kryl%resi, st_kryl%z, st_amgt)
      end if

      sk = DZERO
      !$omp parallel do private(n) reduction(+:sk)
      do n = 1, d_size
        sk = sk + st_kryl%resi(n)*st_kryl%z(n)
      end do
      !$omp end parallel do

#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(sk, "sk l2-norm", sum_sk)
        sk = sum_sk
      end if
#endif

      if (in_iter == 1) then
        !$omp parallel do private(n)
        do n = 1, d_size
          st_kryl%p(n) = st_kryl%z(n)
        end do
        !$omp end parallel do
      else
        beta = sk/sk0
        !$omp parallel do private(n)
        do n = 1, d_size
          st_kryl%p(n) = st_kryl%z(n) + beta*st_kryl%p(n)
        end do
        !$omp end parallel do
      end if

#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_kryl%p)
      end if
#endif
      ! -- Calculate matrix-vector multiplication (matvec)
        call calc_matvec(level, d_size, st_kryl%p, array_var(level)%dmat, array_var(level)%lumat, st_kryl%q)

      sk2 = DZERO
      !$omp parallel do private(n) reduction(+:sk2)
      do n = 1, d_size
        sk2 = sk2 + st_kryl%p(n)*st_kryl%q(n)
      end do
      !$omp end parallel do

#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(sk2, "sk2 l2-norm", sum_sk)
        sk2 = sum_sk
      end if
#endif

      alpha = sk/sk2
      !$omp parallel do private(n)
      do n = 1, d_size
        inx(n) = inx(n) + alpha*st_kryl%p(n)
        st_kryl%resi(n) = st_kryl%resi(n) - alpha*st_kryl%q(n)
      end do
      !$omp end parallel do

      ! -- Calculate l2norm square (l2norm2)
        call calc_l2norm2(level, st_kryl%resi, rnorm)
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(rnorm, "residual l2-norm", sum_rnorm)
        rnorm = sum_rnorm
      end if
#endif
      ! -- Check inner solution (check_insol)
        call check_insol(bnorm, rnorm)

#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_kryl%resi)
          call senrec_rvectv(inx)
      end if
#endif

      if (st_time%conv_flag .or. in_iter == st_ctrl%maxinn_iter) then
        exit pcg_inter
      else
        sk0 = sk
      end if

    end do pcg_inter


  end subroutine solve_pcg

  subroutine solve_bicgs(level, inx, st_kryl, st_amgt)
  !*********************************************************************************************
  ! solve_bicgs -- Solve Preconditioned Bi-Conjugate Gradient Stabilized method
  !*********************************************************************************************
    ! -- module

    ! -- inout
    integer(I4), intent(in) :: level
    real(DP), intent(inout) :: inx(:)
    type(kryl_set), intent(inout) :: st_kryl
    type(amgt_set), intent(inout) :: st_amgt
    ! -- local
    integer(I4) :: n, d_size, lu_size, reg_size
    real(DP) :: sk, sk0, sk2, alpha, beta, bicgs_omega, ts, tt
#ifdef MPI_MSG
    real(DP) :: sum_sk, sum_rnorm
#endif
    !-------------------------------------------------------------------------------------------
    d_size = crs_index(level)%unknow
    lu_size = crs_index(level)%lunum
    reg_size = size(inx)

    sk = DZERO ; sk0 = DINFI ; sk2 = DZERO ; alpha = DZERO ; beta = DZERO
    bicgs_omega = DONE ; ts = DZERO ; tt = DZERO

    !$omp parallel
    !$omp do private(n)
    do n = 1, reg_size
      st_kryl%z(n) = DZERO ; st_kryl%p(n) = DZERO
      st_kryl%d(n) = array_var(level)%dmat(n)
    end do
    !$omp end do
    !$omp do private(n)
    do n = 1, d_size
      st_kryl%v(n) = DZERO ; st_kryl%rs(n) = DZERO ; st_kryl%t(n) = DZERO
    end do
    !$omp end do
    !$omp end parallel

    ! -- Calculate residual (resi)
      call calc_resi(level, d_size, array_var(level)%dmat, array_var(level)%lumat, inx,&
                     array_var(level)%rhs, st_kryl%resi)

    if (st_ctrl%precon_type == 0) then
#ifdef MPI_MSG
      if (st_mpi%totn /= 1 .and. level == 1) then
        ! -- Preconditon mpi incomplete lu diagonal (mpi_dilu)
          call precon_mpi_dilu(array_var(level)%dmat, array_var(level)%lumat, st_kryl%d)
      else
        ! -- Preconditon incomplete lu diagonal (dilu)
          call precon_dilu(level, d_size, array_var(level)%dmat, array_var(level)%lumat, st_kryl%d)
      end if
#else
      ! -- Preconditon incomplete lu diagonal (dilu)
        call precon_dilu(level, d_size, array_var(level)%dmat, array_var(level)%lumat, st_kryl%d)
#endif
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_kryl%d)
      end if
#endif
    end if

    !$omp parallel do private(n)
    do n = 1, d_size
      st_kryl%rs(n) = st_kryl%resi(n)
    end do
    !$omp end parallel do

    bicg_inter: do in_iter = 1, st_ctrl%maxinn_iter
      sk = DZERO
      !$omp parallel do private(n) reduction(+:sk)
      do n = 1, d_size
        sk = sk + st_kryl%resi(n)*st_kryl%rs(n)
      end do
      !$omp end parallel do
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(sk, "sk l2-norm", sum_sk)
        sk = sum_sk
      end if
#endif

      if (in_iter == 1) then
        !$omp parallel do private(n)
        do n = 1, d_size
          st_kryl%p(n) = st_kryl%resi(n)
        end do
        !$omp end parallel do
      else
        beta = (sk/sk0)*(alpha/bicgs_omega)
        !$omp parallel do private(n)
        do n = 1, d_size
          st_kryl%p(n) = st_kryl%resi(n) + beta*(st_kryl%p(n)-bicgs_omega*st_kryl%v(n))
        end do
        !$omp end parallel do
      end if
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_kryl%p)
      end if
#endif

      if (st_ctrl%precon_type == 0) then
#ifdef MPI_MSG
        if (st_mpi%totn /= 1 .and. level == 1) then
          ! -- Solve mpi ilu factorization (ilu)
            call solve_mpi_ilu(st_kryl%p, st_kryl%d, array_var(level)%lumat, st_kryl%z)
        else
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(level, d_size, st_kryl%p, st_kryl%d, array_var(level)%lumat, st_kryl%z)
        end if
#else
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(level, d_size, st_kryl%p, st_kryl%d, array_var(level)%lumat, st_kryl%z)
#endif
      else if (st_ctrl%precon_type == 1) then
        ! -- Loop cycle for amg (amg)
          call loop_amg(level, st_kryl%p, st_kryl%z, st_amgt)
      end if

#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_kryl%z)
      end if
#endif

      ! -- Calculate matrix-vector multiplication (matvec)
        call calc_matvec(level, d_size, st_kryl%z, array_var(level)%dmat, array_var(level)%lumat, st_kryl%v)

      sk2 = DZERO
      !$omp parallel do private(n) reduction(+:sk2)
      do n = 1, d_size
        sk2 = sk2 + st_kryl%v(n)*st_kryl%rs(n)
      end do
      !$omp end parallel do
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(sk2, "sk2 l2-norm", sum_sk)
        sk2 = sum_sk
      end if
#endif
      alpha = sk/sk2

      !$omp parallel do private(n)
      do n = 1, d_size
        inx(n) = inx(n) + alpha*st_kryl%z(n)
        st_kryl%resi(n) = st_kryl%resi(n) - alpha*st_kryl%v(n)
      end do
      !$omp end parallel do

      ! -- Calculate l2norm square (l2norm2)
        call calc_l2norm2(level, st_kryl%resi, rnorm)
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Sum value for MPI (val)
          call mpisum_val(rnorm, "residual l2-norm", sum_rnorm)
        rnorm = sum_rnorm
      end if
#endif
      ! -- Check inner solution (check_insol)
        call check_insol(bnorm, rnorm)

#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_kryl%resi)
          call senrec_rvectv(inx)
      end if
#endif

      if (st_time%conv_flag) then
        exit bicg_inter
      else
        if (st_ctrl%precon_type == 0) then
#ifdef MPI_MSG
          if (st_mpi%totn /= 1 .and. level == 1) then
            ! -- Solve mpi ilu factorization (ilu)
              call solve_mpi_ilu(st_kryl%resi, st_kryl%d, array_var(level)%lumat, st_kryl%z)
          else
            ! -- Solve ilu factorization (ilu)
              call solve_ilu(level, d_size, st_kryl%resi, st_kryl%d, array_var(level)%lumat, st_kryl%z)
          end if
#else
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(level, d_size, st_kryl%resi, st_kryl%d, array_var(level)%lumat, st_kryl%z)
#endif
        else if (st_ctrl%precon_type == 1) then
        ! -- Loop cycle for amg (amg)
          call loop_amg(level, st_kryl%resi, st_kryl%z, st_amgt)
        end if

#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(st_kryl%z)
        end if
#endif

        ! -- Calculate matrix-vector multiplication (matvec)
          call calc_matvec(level, d_size, st_kryl%z, array_var(level)%dmat, array_var(level)%lumat, st_kryl%t)

        ts = DZERO ; tt = DZERO
        !$omp parallel
        !$omp do private(n) reduction(+:ts)
        do n = 1, d_size
          ts = ts + st_kryl%t(n)*st_kryl%resi(n)
        end do
        !$omp end do
        !$omp do private(n) reduction(+:tt)
        do n = 1, d_size
          tt = tt + st_kryl%t(n)*st_kryl%t(n)
        end do
        !$omp end do
        !$omp end parallel

#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
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
          inx(n) = inx(n) + bicgs_omega*st_kryl%z(n)
          st_kryl%resi(n) = st_kryl%resi(n) - bicgs_omega*st_kryl%t(n)
        end do
        !$omp end parallel do

        ! -- Calculate l2norm square (l2norm2)
          call calc_l2norm2(level, st_kryl%resi, rnorm)
#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
          ! -- Sum value for MPI (val)
            call mpisum_val(rnorm, "residual l2-norm", sum_rnorm)
          rnorm = sum_rnorm
        end if
#endif
        ! -- Check inner solution (check_insol)
          call check_insol(bnorm, rnorm)

#ifdef MPI_MSG
        if (st_mpi%totn /= 1) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(st_kryl%resi)
            call senrec_rvectv(inx)
        end if
#endif

        if (st_time%conv_flag .or. in_iter == st_ctrl%maxinn_iter) then
          exit bicg_inter
        else
          sk0 = sk
        end if
      end if

    end do bicg_inter


  end subroutine solve_bicgs

  subroutine loop_amg(alevel, r, e, st_amgt)
  !*********************************************************************************************
  ! loop_amg -- Loop cycle for amg
  !*********************************************************************************************
    ! -- module
    use allocate_solution, only: pro_var, res_var
    ! -- inout
    integer(I4), intent(in) :: alevel
    real(DP), intent(in) :: r(:)
    real(DP), intent(out) :: e(:)
    type(amgt_set), intent(inout) :: st_amgt
    ! -- local
    integer(I4) :: i, j, k, d_size, lu_size, reg_size
    integer(I4) :: mgd_size, mgreg_size
    integer(I4) :: v_iter, vlevel, ncoa, nfin
    integer(I4) :: rst, ren, pst, pen
    !-------------------------------------------------------------------------------------------
    mgd_size = crs_index(alevel)%unknow ; mgreg_size = size(array_var(alevel)%x)
    !$omp parallel
    !$omp do private(i)
    do i = 1, mgd_size
      array_var(alevel)%x(i) = DZERO
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, mgreg_size
      st_amgt%save_rhs(i) = array_var(alevel)%rhs(i)
      array_var(alevel)%rhs(i) = r(i)
      st_amgt%dilu_d(i) = DZERO
    end do
    !$omp end do
    !$omp end parallel
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Preconditon mpi incomplete lu diagonal (mpi_dilu)
        call precon_mpi_dilu(array_var(alevel)%dmat, array_var(alevel)%lumat, st_amgt%dilu_d)
    else
      ! -- Preconditon incomplete lu diagonal (dilu)
        call precon_dilu(alevel, mgd_size, array_var(alevel)%dmat, array_var(alevel)%lumat, st_amgt%dilu_d)
    end if
    if (st_mpi%totn /= 1) then
      call senrec_rvectv(st_amgt%dilu_d)
    end if
#else
    ! -- Preconditon incomplete lu diagonal (dilu)
      call precon_dilu(alevel, mgd_size, array_var(alevel)%dmat, array_var(alevel)%lumat, st_amgt%dilu_d)
#endif

    do v_iter = 1, st_ctrl%maxvcy_iter  ! V-cycle loop
      do vlevel = alevel, st_ctrl%nlevel-1

        d_size = crs_index(vlevel)%unknow
        lu_size = crs_index(vlevel)%lunum
        reg_size = size(array_var(vlevel)%x)
        !$omp parallel
        !$omp do private(i)
        do i = 1, reg_size
          st_amgt%td(i) = array_var(vlevel)%dmat(i)
          st_amgt%tx(i) = array_var(vlevel)%x(i)
          st_amgt%tb(i) = array_var(vlevel)%rhs(i)
          st_amgt%tr(i) = DZERO
        end do
        !$omp end do
        !$omp do private(i)
        do i = 1, lu_size
          st_amgt%tlu(i) = array_var(vlevel)%lumat(i)
        end do
        !$omp end do
        if (vlevel == alevel) then
          !$omp do private(i)
          do i = 1, d_size
            st_amgt%td(i) = st_amgt%dilu_d(i)
          end do
          !$omp end do
        end if
        !$omp end parallel

#ifdef MPI_MSG
        if (st_mpi%totn /= 1 .and. vlevel == alevel) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(st_amgt%tx(1:reg_size))
            call senrec_rvectv(st_amgt%tb(1:reg_size))
        end if
        if (st_mpi%totn /= 1 .and. vlevel == alevel) then
          ! -- Solve mpi ilu factorization (ilu)
            call solve_mpi_ilu(st_amgt%tb(1:reg_size), st_amgt%td(1:reg_size),&
                               st_amgt%tlu(1:lu_size),&
                               st_amgt%tx(1:reg_size))
        else
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(vlevel, d_size, st_amgt%tb(1:reg_size), st_amgt%td(1:reg_size),&
                           st_amgt%tlu(1:lu_size), st_amgt%tx(1:reg_size))
        end if
#else
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(vlevel, d_size, st_amgt%tb(1:reg_size), st_amgt%td(1:reg_size),&
                         st_amgt%tlu(1:lu_size), st_amgt%tx(1:reg_size))
#endif

#ifdef MPI_MSG
        if (st_mpi%totn /= 1 .and. vlevel == alevel) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(st_amgt%tx(1:reg_size))
        end if
#endif

        ! -- Calculate residual
          call calc_resi(vlevel, d_size, st_amgt%td(1:reg_size), st_amgt%tlu(1:lu_size),&
                         st_amgt%tx(1:reg_size), st_amgt%tb(1:reg_size), st_amgt%tr(1:reg_size))

#ifdef MPI_MSG
        if (st_mpi%totn /= 1 .and. vlevel == alevel) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(st_amgt%tr(1:reg_size))
        end if
#endif

        ncoa = crs_index(vlevel+1)%unknow
        !$omp parallel
        !$omp do private(i)
        do i = 1, ncoa
          st_amgt%trhs(i) = DZERO
        end do
        !$omp end do
        !$omp do private(i, j, k, rst, ren)
        do i = 1, ncoa
          rst = res_var(vlevel+1)%rindex(i-1) + 1
          ren = res_var(vlevel+1)%rindex(i)
          do k = rst, ren
            j = res_var(vlevel+1)%roffrow(k)
            st_amgt%trhs(i) = st_amgt%trhs(i) + res_var(vlevel+1)%rval(k)*st_amgt%tr(j)
          end do
          array_var(vlevel+1)%rhs(i) = st_amgt%trhs(i)
        end do
        !$omp end do

        !$omp do private(i)
        do i = 1, ncoa
          array_var(vlevel+1)%x(i) = DZERO
        end do
        !$omp end do
        !$omp end parallel
      end do

      d_size = crs_index(st_ctrl%nlevel)%unknow
      lu_size = crs_index(st_ctrl%nlevel)%lunum
      reg_size = size(array_var(st_ctrl%nlevel)%x)

      !$omp parallel
      !$omp do private(i)
      do i = 1, d_size
        st_amgt%td(i) = array_var(st_ctrl%nlevel)%dmat(i)
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, lu_size
        st_amgt%tlu(i) = array_var(st_ctrl%nlevel)%lumat(i)
      end do
      !$omp end do
      !$omp do private(i)
      do i = 1, reg_size
        st_amgt%tx(i) = array_var(st_ctrl%nlevel)%x(i)
        st_amgt%tb(i) = array_var(st_ctrl%nlevel)%rhs(i)
      end do
      !$omp end do
      if (st_ctrl%nlevel == alevel) then
        !$omp do private(i)
        do i = 1, d_size
          st_amgt%td(i) = st_amgt%dilu_d(i)
        end do
        !$omp end do
      end if
      !$omp end parallel

#ifdef MPI_MSG
      if (st_mpi%totn /= 1 .and. st_ctrl%nlevel == alevel) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_amgt%tx(1:reg_size))
          call senrec_rvectv(st_amgt%tb(1:reg_size))
      end if
      if (st_mpi%totn /= 1 .and. st_ctrl%nlevel == alevel) then
        ! -- Solve mpi ilu factorization (ilu)
          call solve_mpi_ilu(st_amgt%tb(1:reg_size), st_amgt%td(1:d_size),&
                             st_amgt%tlu(1:lu_size),&
                             st_amgt%tx(1:reg_size))
      else
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(st_ctrl%nlevel, d_size, st_amgt%tb(1:reg_size), st_amgt%td(1:d_size),&
                         st_amgt%tlu(1:lu_size), st_amgt%tx(1:reg_size))
      end if
#else
        ! -- Solve ilu factorization (ilu)
        call solve_ilu(st_ctrl%nlevel, d_size, st_amgt%tb(1:reg_size), st_amgt%td(1:d_size),&
                       st_amgt%tlu(1:lu_size), st_amgt%tx(1:reg_size))
#endif

#ifdef MPI_MSG
      if (st_mpi%totn /= 1 .and. st_ctrl%nlevel == alevel) then
        ! -- Send and Receive real vector value (rvectv)
          call senrec_rvectv(st_amgt%tx(1:reg_size))
      end if
#endif
      array_var(st_ctrl%nlevel)%x(:) = st_amgt%tx(1:reg_size)

      do vlevel = st_ctrl%nlevel-1, alevel, -1
        nfin = crs_index(vlevel+1)%unknow
        !$omp parallel
        !$omp do private(i)
        do i = 1, nfin
          st_amgt%tfx(i) = DZERO
        end do
        !$omp end do
        !$omp do private(i, j, k, pst, pen)
        do i = 1, nfin
          pst = pro_var(vlevel+1)%pindex(i-1) + 1
          pen = pro_var(vlevel+1)%pindex(i)
          do k = pst, pen
            j = pro_var(vlevel+1)%poffrow(k)
            st_amgt%tfx(i) = st_amgt%tfx(i) + pro_var(vlevel+1)%pval(k)*array_var(vlevel+1)%x(j)
          end do
          array_var(vlevel)%x(i) = array_var(vlevel)%x(i) + st_amgt%tfx(i)
        end do
        !$omp end do
        !$omp end parallel

        d_size = crs_index(vlevel)%unknow
        lu_size = crs_index(vlevel)%lunum
        reg_size = size(array_var(vlevel)%x)

        !$omp parallel
        !$omp do private(i)
        do i = 1, reg_size
          st_amgt%td(i) = array_var(vlevel)%dmat(i)
          st_amgt%tx(i) = array_var(vlevel)%x(i)
          st_amgt%tb(i) = array_var(vlevel)%rhs(i)
        end do
        !$omp end do
        !$omp do private(i)
        do i = 1, lu_size
          st_amgt%tlu(i) = array_var(vlevel)%lumat(i)
        end do
        !$omp end do
        if (vlevel == alevel) then
          !$omp do private(i)
          do i = 1, d_size
            st_amgt%td(i) = st_amgt%dilu_d(i)
          end do
          !$omp end do
        end if
        !$omp end parallel

#ifdef MPI_MSG
        if (st_mpi%totn /= 1 .and. vlevel == alevel) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(st_amgt%tx(1:reg_size))
            call senrec_rvectv(st_amgt%tb(1:reg_size))
        end if
        if (st_mpi%totn /= 1 .and. vlevel == alevel) then
          ! -- Solve mpi ilu factorization (ilu)
            call solve_mpi_ilu(st_amgt%tb(1:reg_size), st_amgt%td(1:reg_size),&
                               st_amgt%tlu(1:lu_size),&
                               st_amgt%tx(1:reg_size))
        else
          ! -- Solve ilu factorization (ilu)
            call solve_ilu(vlevel, d_size, st_amgt%tb(1:reg_size), st_amgt%td(1:reg_size),&
                           st_amgt%tlu(1:lu_size), st_amgt%tx(1:reg_size))
        end if
#else
        ! -- Solve ilu factorization (ilu)
          call solve_ilu(vlevel, d_size, st_amgt%tb(1:reg_size), st_amgt%td(1:reg_size),&
                         st_amgt%tlu(1:lu_size), st_amgt%tx(1:reg_size))
#endif

        array_var(vlevel)%x(:) = st_amgt%tx(1:reg_size)
        array_var(vlevel)%rhs(:) = st_amgt%tb(1:reg_size)

#ifdef MPI_MSG
        if (st_mpi%totn /= 1 .and. vlevel == alevel) then
          ! -- Send and Receive real vector value (rvectv)
            call senrec_rvectv(array_var(vlevel)%x)
            call senrec_rvectv(array_var(vlevel)%rhs)
        end if
#endif

      end do
    end do

    !$omp parallel
    !$omp do private(i)
    do i = 1, mgd_size
      e(i) = array_var(alevel)%x(i)
    end do
    !$omp end do
    !$omp do private(i)
    do i = 1, mgreg_size
      array_var(alevel)%rhs(i) = st_amgt%save_rhs(i)
    end do
    !$omp end do
    !$omp end parallel


  end subroutine loop_amg

  subroutine precon_dilu(plevel, npre, pre_ind, pre_inlu, pre_d)
  !*********************************************************************************************
  ! precon_dilu -- Preconditon incomplete lu diagonal
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: plevel, npre
    real(DP), intent(in) :: pre_ind(:), pre_inlu(:)
    real(DP), intent(inout) :: pre_d(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: off_sta, off_end, off_sta2, off_end2
    integer(I4) :: offr, offr2
    real(DP) :: d_invk, d_flor
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, npre
      pre_d(i) = pre_ind(i)
    end do
    !$omp end parallel do
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
      if (st_ctrl%dilu_shift > DZERO) then
        d_flor = st_ctrl%dilu_shift*abs(pre_ind(i))
        if (abs(pre_d(i)) < d_flor) then
          pre_d(i) = sign(d_flor, pre_ind(i))
          dilu_shift_num = dilu_shift_num + 1
        end if
      end if
    end do

  end subroutine precon_dilu

  subroutine solve_ilu(plevel, npre, inrhs, indmat, inlumat, outx)
  !*********************************************************************************************
  ! solve_ilu -- Solve ilu factorization
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: plevel, npre
    real(DP), intent(in) :: inrhs(:), indmat(:), inlumat(:)
    real(DP), intent(out) :: outx(:)
    ! -- local
    integer(I4) :: i, k
    integer(I4) :: off_sta, off_end, offr
    real(DP) :: temp_outx
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i)
    do i = 1, npre
      outx(i) = inrhs(i)
    end do
    !$omp end parallel do

    ! Forward Substitution
!    !$omp parallel do private(i, k, off_sta, off_end, offr, temp_outx)
    do i = 1, npre
      temp_outx = outx(i)
      off_sta = crs_index(plevel)%offind(i-1) + 1
      off_end = crs_index(plevel)%offind(i)
      do k = off_sta, off_end
        offr = crs_index(plevel)%offrow(k)
        if (3 >= dir_conn(k)) then
          temp_outx = temp_outx - inlumat(k)*outx(offr)
        end if
      end do
      outx(i) = temp_outx/indmat(i)
    end do
!    !$omp end parallel do

    ! Backward Substitution
!    !$omp parallel do private(i, k, off_sta, off_end, offr, temp_outx)
    do i = npre, 1, -1
      temp_outx = DZERO
      off_sta = crs_index(plevel)%offind(i-1) + 1
      off_end = crs_index(plevel)%offind(i)
      do k = off_sta, off_end
        offr = crs_index(plevel)%offrow(k)
        if (3 < dir_conn(k)) then
          temp_outx = temp_outx + inlumat(k)*outx(offr)
        end if
      end do
      outx(i) = outx(i) - temp_outx/indmat(i)
    end do
!    !$omp end parallel do

  end subroutine solve_ilu

  subroutine calc_matvec(mvlevel, nmv, invec, indmat, inlumat, outvec)
  !*********************************************************************************************
  ! calc_matvec -- Calculate matrix-vector multiplication
  !*********************************************************************************************
    ! -- modules

    ! -- inout
    integer(I4), intent(in) :: mvlevel, nmv
    real(DP), intent(in) :: invec(:), indmat(:), inlumat(:)
    real(DP), intent(out) :: outvec(:)
    ! -- local
    integer(I4) :: i, j, k
    integer(I4) :: off_sta, off_end
    real(DP) :: temp_vec
    !-------------------------------------------------------------------------------------------
    !$omp parallel do private(i, j, k, off_sta, off_end, temp_vec)
    do i = 1, nmv
      temp_vec = indmat(i)*invec(i)
      off_sta = crs_index(mvlevel)%offind(i-1) + 1
      off_end = crs_index(mvlevel)%offind(i)
      do k = off_sta, off_end
        j = crs_index(mvlevel)%offrow(k)
        temp_vec = temp_vec + inlumat(k)*invec(j)
      end do
      outvec(i) = temp_vec
    end do
    !$omp end parallel do

  end subroutine calc_matvec

end module linear_solution
