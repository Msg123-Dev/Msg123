module make_cell
  ! -- modules
  use kind_module, only: I4, DP
  use constval_module, only: DZERO, DHALF, FACE
  use read_input, only: glob_x, glob_y, glob_z
  use set_cell, only: ncals, ncalc, get_calc_grid

  implicit none
  private
  public :: make_cell_info
  real(DP), allocatable, public :: dis2face(:,:), face_area(:,:), area_r(:)
  real(DP), allocatable, public :: cell_vol(:), surf_elev(:)
  real(DP), allocatable, public :: cell_top(:), cell_cent(:), cell_bot(:)

  ! -- local
  real(DP), allocatable :: fp_xw(:), fp_yw(:), fp_xs(:), fp_ys(:), fp_xe(:), fp_ye(:)
  real(DP), allocatable :: fp_xn(:), fp_yn(:), fp_zw(:), fp_zs(:), fp_ze(:), fp_zn(:)
  real(DP), allocatable :: fp_zt(:), fp_zb(:), cp_x(:), cp_y(:), cp_z(:)

  contains

  subroutine make_cell_info()
  !***************************************************************************************
  ! make_cell_info -- Make cell information
  !***************************************************************************************
    ! -- modules
    use initial_module, only: my_rank, st_grid
#ifdef MPI_MSG
    use initial_module, only: pro_totn
    use mpi_set, only: bcast_glob_xyzv
#endif
    ! -- inout

    ! -- local

    !-------------------------------------------------------------------------------------
    if (my_rank /= 0) then
      allocate(glob_x(st_grid%nx+1,st_grid%ny+1), glob_y(st_grid%nx+1,st_grid%ny+1))
      allocate(glob_z(st_grid%nx+1,st_grid%ny+1,st_grid%nz+1))
      !$omp parallel workshare
      glob_x(:,:) = DZERO ; glob_y(:,:) = DZERO ; glob_z(:,:,:) = DZERO
      !$omp end parallel workshare
    end if

#ifdef MPI_MSG
    if (pro_totn /= 1) then
      ! -- Bcast global xyz value (glob_xyzv)
        call bcast_glob_xyzv()
    end if
#endif

    ! -- Make center point on each face and cell (cent_point)
      call make_cent_point()

    ! -- Make distance between point and face center (dis_point2face)
      call make_dis_point2face()

    ! -- Make face area (face_area)
      call make_face_area()

    ! -- Make cell volume (cell_vol)
      call make_cell_vol()

    ! -- Make cell surface elevation (cell_surfelev)
      call make_cell_surfelev()

    ! -- Make recharge area (rech_area)
      call make_rech_area()

      deallocate(glob_x, glob_y, glob_z)

  end subroutine make_cell_info

  subroutine make_cent_point()
  !***************************************************************************************
  ! make_cent_point -- Make center point information
  !***************************************************************************************
    ! -- modules
    use constval_module, only: DQUA
    ! -- inout

    ! -- local
    integer(I4) :: i, xn, yn, zn
    !-------------------------------------------------------------------------------------
    allocate(fp_xw(ncalc), fp_yw(ncalc), fp_xs(ncalc), fp_ys(ncalc))
    allocate(fp_xe(ncalc), fp_ye(ncalc), fp_xn(ncalc), fp_yn(ncalc))
    allocate(fp_zw(ncalc), fp_zs(ncalc), fp_ze(ncalc), fp_zn(ncalc))
    allocate(fp_zt(ncalc), fp_zb(ncalc))
    allocate(cp_x(ncalc), cp_y(ncalc), cp_z(ncalc))
    !$omp parallel
    !$omp workshare
    fp_xw(:) = DZERO ; fp_yw(:) = DZERO ; fp_xs(:) = DZERO ; fp_ys(:) = DZERO
    fp_xe(:) = DZERO ; fp_ye(:) = DZERO ; fp_xn(:) = DZERO ; fp_yn(:) = DZERO
    fp_zw(:) = DZERO ; fp_zs(:) = DZERO ; fp_ze(:) = DZERO ; fp_zn(:) = DZERO
    fp_zt(:) = DZERO ; fp_zb(:) = DZERO
    cp_x(:) = DZERO ; cp_y(:) = DZERO ; cp_z(:) = DZERO
    !$omp end workshare

    ! make center point for face and cell
    !$omp do private(i, xn, yn, zn)
    do i = 1, ncalc
      call get_calc_grid(i, xn, yn, zn)
      fp_xw(i) = DHALF*(glob_x(xn,yn)+glob_x(xn,yn+1))
      fp_yw(i) = DHALF*(glob_y(xn,yn)+glob_y(xn,yn+1))
      fp_xe(i) = DHALF*(glob_x(xn+1,yn)+glob_x(xn+1,yn+1))
      fp_ye(i) = DHALF*(glob_y(xn+1,yn)+glob_y(xn+1,yn+1))
      fp_xn(i) = DHALF*(glob_x(xn,yn)+glob_x(xn+1,yn))
      fp_yn(i) = DHALF*(glob_y(xn,yn)+glob_y(xn+1,yn))
      fp_xs(i) = DHALF*(glob_x(xn,yn+1)+glob_x(xn+1,yn+1))
      fp_ys(i) = DHALF*(glob_y(xn,yn+1)+glob_y(xn+1,yn+1))
      cp_x(i) = DQUA*(glob_x(xn,yn)+glob_x(xn,yn+1)+glob_x(xn+1,yn)+glob_x(xn+1,yn+1))
      cp_y(i) = DQUA*(glob_y(xn,yn)+glob_y(xn,yn+1)+glob_y(xn+1,yn)+glob_y(xn+1,yn+1))
      fp_zw(i) = DQUA*(glob_z(xn,yn,zn)+glob_z(xn,yn,zn+1)+glob_z(xn,yn+1,zn)+glob_z(xn,yn+1,zn+1))
      fp_ze(i) = DQUA*(glob_z(xn+1,yn,zn)+glob_z(xn+1,yn+1,zn)+glob_z(xn+1,yn,zn+1)+glob_z(xn+1,yn+1,zn+1))
      fp_zn(i) = DQUA*(glob_z(xn,yn,zn)+glob_z(xn,yn,zn+1)+glob_z(xn+1,yn,zn)+glob_z(xn+1,yn,zn+1))
      fp_zs(i) = DQUA*(glob_z(xn,yn+1,zn)+glob_z(xn+1,yn+1,zn)+glob_z(xn,yn+1,zn+1)+glob_z(xn+1,yn+1,zn+1))
      fp_zt(i) = DQUA*(glob_z(xn,yn,zn)+glob_z(xn,yn+1,zn)+glob_z(xn+1,yn,zn)+glob_z(xn+1,yn+1,zn))
      fp_zb(i) = DQUA*(glob_z(xn,yn,zn+1)+glob_z(xn,yn+1,zn+1)+glob_z(xn+1,yn,zn+1)+glob_z(xn+1,yn+1,zn+1))
      cp_z(i) = DHALF*(fp_zt(i)+fp_zb(i))
    end do
    !$omp end do
    !$omp end parallel

  end subroutine make_cent_point

  subroutine make_dis_point2face()
  !***************************************************************************************
  ! make_dis_point2face -- Make distance between point and face center
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i
    real(DP) :: dis2, dis3, dis4, dis5
    !-------------------------------------------------------------------------------------
    allocate(dis2face(ncalc,FACE))
    !$omp parallel
    !$omp workshare
    dis2face(:,:) = DZERO
    !$omp end workshare

    ! make distance between point and face center
    !$omp do private(i, dis2, dis3, dis4, dis5)
    do i = 1, ncalc
      dis2face(i,1) = fp_zt(i) - cp_z(i)
      dis2face(i,6) = cp_z(i) - fp_zb(i)
      dis3 = (cp_x(i)-fp_xw(i))**2 + (cp_y(i)-fp_yw(i))**2 + (cp_z(i)-fp_zw(i))**2
      dis4 = (cp_x(i)-fp_xe(i))**2 + (cp_y(i)-fp_ye(i))**2 + (cp_z(i)-fp_ze(i))**2
      dis2 = (cp_x(i)-fp_xn(i))**2 + (cp_y(i)-fp_yn(i))**2 + (cp_z(i)-fp_zn(i))**2
      dis5 = (cp_x(i)-fp_xs(i))**2 + (cp_y(i)-fp_ys(i))**2 + (cp_z(i)-fp_zs(i))**2
      dis2face(i,3) = sqrt(dis3) ; dis2face(i,4) = sqrt(dis4)
      dis2face(i,2) = sqrt(dis2) ; dis2face(i,5) = sqrt(dis5)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(fp_xw, fp_yw, fp_xs, fp_ys, fp_xe, fp_ye, fp_xn, fp_yn)
    deallocate(fp_zw, fp_ze, fp_zs, fp_zn)

  end subroutine make_dis_point2face

  subroutine make_face_area()
  !***************************************************************************************
  ! meke_face_area -- Make face area
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, xn, yn, zn
    real(DP) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, z5, z6, z7, z8
    !-------------------------------------------------------------------------------------
    allocate(face_area(ncalc,FACE))
    !$omp parallel
    !$omp workshare
    face_area(:,:) = DZERO
    !$omp end workshare

    ! make face area
    !$omp do private(i, xn, yn, zn, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, z5, z6, z7, z8)
    do i = 1, ncalc
      call get_calc_grid(i, xn, yn, zn)
      x1 = glob_x(xn,yn) ; x2 = glob_x(xn+1,yn) ; x3 = glob_x(xn+1,yn+1)
      x4 = glob_x(xn,yn+1) ; y1 = glob_y(xn,yn) ; y2 = glob_y(xn+1,yn)
      y3 = glob_y(xn+1,yn+1) ; y4 = glob_y(xn,yn+1) ; z1 = glob_z(xn,yn,zn)
      z2 = glob_z(xn+1,yn,zn) ; z3 = glob_z(xn+1,yn+1,zn) ; z4 = glob_z(xn,yn+1,zn)
      z5 = glob_z(xn,yn,zn+1) ; z6 = glob_z(xn+1,yn,zn+1)
      z7 = glob_z(xn+1,yn+1,zn+1) ; z8 = glob_z(xn,yn+1,zn+1)

      face_area(i,1) = fa_2tri(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4)
      face_area(i,6) = fa_2tri(x1, x2, x3, x4, y1, y2, y3, y4, z5, z6, z7, z8)
      face_area(i,3) = fa_2tri(x1, x1, x4, x4, y1, y1, y4, y4, z1, z5, z8, z4)
      face_area(i,4) = fa_2tri(x2, x2, x3, x3, y2, y2, y3, y3, z2, z6, z7, z3)
      face_area(i,2) = fa_2tri(x1, x2, x2, x1, y1, y2, y2, y1, z1, z2, z6, z5)
      face_area(i,5) = fa_2tri(x4, x3, x3, x4, y4, y3, y3, y4, z4, z3, z7, z8)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine make_face_area

!  function fa_diag(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4) result(quadarea)
!  !***************************************************************************************
!  ! fa_diag -- Calculate face area using cross product with diagonal vector
!  !***************************************************************************************
!    ! -- modules
!
!    ! -- inout
!    real(DP), intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
!    ! -- local
!    real(DP) :: quadarea
!    real(DP) :: vec3d1(3), vec3d2(3)
!    !-------------------------------------------------------------------------------------
!    vec3d1(1) = x3 - x1 ; vec3d1(2) = y3 - y1 ; vec3d1(3) = z3 - z1
!    vec3d2(1) = x4 - x2 ; vec3d2(2) = y4 - y2 ; vec3d2(3) = z4 - z2
!    quadarea = (vec3d1(2)*vec3d2(3) - vec3d1(3)*vec3d2(2))**2&
!             + (vec3d1(3)*vec3d2(1) - vec3d1(1)*vec3d2(3))**2&
!             + (vec3d1(1)*vec3d2(2) - vec3d1(2)*vec3d2(1))**2
!
!    quadarea = sqrt(quadarea)*DHALF
!
!  end function fa_diag

  function fa_2tri(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4) result(triarea)
  !***************************************************************************************
  ! fa_2tri -- Calculate face area using cross product with two triangle
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
    ! -- local
    real(DP) :: triarea, triarea1, triarea2
    real(DP) :: vec3d1(3), vec3d2(3)
    !-------------------------------------------------------------------------------------
    vec3d1(1) = x1 - x2 ; vec3d1(2) = y1 - y2 ; vec3d1(3) = z1 - z2
    vec3d2(1) = x3 - x2 ; vec3d2(2) = y3 - y2 ; vec3d2(3) = z3 - z2
    triarea1 = (vec3d1(2)*vec3d2(3) - vec3d1(3)*vec3d2(2))**2&
             + (vec3d1(3)*vec3d2(1) - vec3d1(1)*vec3d2(3))**2&
             + (vec3d1(1)*vec3d2(2) - vec3d1(2)*vec3d2(1))**2

    vec3d1(1) = x1 - x4 ; vec3d1(2) = y1 - y4 ; vec3d1(3) = z1 - z4
    vec3d2(1) = x3 - x4 ; vec3d2(2) = y3 - y4 ; vec3d2(3) = z3 - z4
    triarea2 = (vec3d1(2)*vec3d2(3) - vec3d1(3)*vec3d2(2))**2&
             + (vec3d1(3)*vec3d2(1) - vec3d1(1)*vec3d2(3))**2&
             + (vec3d1(1)*vec3d2(2) - vec3d1(2)*vec3d2(1))**2

    triarea = (sqrt(triarea1)+sqrt(triarea2))*DHALF

  end function fa_2tri

  subroutine make_cell_vol()
  !***************************************************************************************
  ! meke_cell_vol -- Make cell volume
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i, xn, yn, zn
    real(DP) :: cx, cy, cz
    real(DP) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, z5, z6, z7, z8
    !-------------------------------------------------------------------------------------
    allocate(cell_vol(ncalc))
    !$omp parallel
    !$omp workshare
    cell_vol(:) = DZERO
    !$omp end workshare

    ! make cell volume
    !$omp do private(i, xn, yn, zn, cx, cy, cz, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, z5, z6, z7, z8)
    do i = 1, ncalc
      call get_calc_grid(i, xn, yn, zn)
      cx = cp_x(i) ; cy = cp_y(i) ; cz = cp_z(i)
      x1 = glob_x(xn,yn) ; x2 = glob_x(xn+1,yn) ; x3 = glob_x(xn+1,yn+1)
      x4 = glob_x(xn,yn+1) ; y1 = glob_y(xn,yn) ; y2 = glob_y(xn+1,yn)
      y3 = glob_y(xn+1,yn+1) ; y4 = glob_y(xn,yn+1) ; z1 = glob_z(xn,yn,zn)
      z2 = glob_z(xn+1,yn,zn) ; z3 = glob_z(xn+1,yn+1,zn) ; z4 = glob_z(xn,yn+1,zn)
      z5 = glob_z(xn,yn,zn+1) ; z6 = glob_z(xn+1,yn,zn+1)
      z7 = glob_z(xn+1,yn+1,zn+1) ; z8 = glob_z(xn,yn+1,zn+1)

      cell_vol(i) = vol(cx, x1, x4, x4, cy, y1, y4, y4, cz, z1, z4, z8) + &
                    vol(cx, x4, x1, x1, cy, y4, y1, y1, cz, z8, z5, z1) + &
                    vol(cx, x2, x3, x3, cy, y2, y3, y3, cz, z2, z3, z7) + &
                    vol(cx, x3, x2, x2, cy, y3, y2, y2, cz, z7, z6, z2) + &
                    vol(cx, x1, x1, x2, cy, y1, y1, y2, cz, z1, z5, z6) + &
                    vol(cx, x2, x2, x1, cy, y2, y2, y1, cz, z6, z2, z1) + &
                    vol(cx, x4, x4, x3, cy, y4, y4, y3, cz, z4, z8, z7) + &
                    vol(cx, x3, x3, x4, cy, y3, y3, y4, cz, z7, z3, z4) + &
                    vol(cx, x1, x4, x3, cy, y1, y4, y3, cz, z1, z4, z3) + &
                    vol(cx, x3, x2, x1, cy, y3, y2, y1, cz, z3, z2, z1) + &
                    vol(cx, x1, x4, x3, cy, y1, y4, y3, cz, z5, z8, z7) + &
                    vol(cx, x3, x2, x1, cy, y3, y2, y1, cz, z7, z6, z5)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(cp_x, cp_y, cp_z)

  end subroutine make_cell_vol

  function vol(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4) result(trivol)
  !***************************************************************************************
  ! vol -- Calculate cell vol using cross product with a focus on cell center
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
    ! -- local
    real(DP) :: trivol
    real(DP) :: vec3d1(3), vec3d2(3), vec3d3(3), cropro12(3)
    !-------------------------------------------------------------------------------------
    vec3d1(1) = x2 - x3 ; vec3d1(2) = y2 - y3 ; vec3d1(3) = z2 - z3
    vec3d2(1) = x4 - x3 ; vec3d2(2) = y4 - y3 ; vec3d2(3) = z4 - z3
    vec3d3(1) = x1 - x3 ; vec3d3(2) = y1 - y3 ; vec3d3(3) = z1 - z3
    cropro12(1) = vec3d1(2)*vec3d2(3) - vec3d1(3)*vec3d2(2)
    cropro12(2) = vec3d1(3)*vec3d2(1) - vec3d1(1)*vec3d2(3)
    cropro12(3) = vec3d1(1)*vec3d2(2) - vec3d1(2)*vec3d2(1)

    trivol = abs(dot_product(cropro12(1:3), vec3d3(1:3)))/(6.00_DP)

  end function vol

  subroutine make_cell_surfelev()
  !***************************************************************************************
  ! make_cell_surfelev -- Make cell surface elevation
  !***************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------
    allocate(cell_top(ncalc), cell_cent(ncalc), cell_bot(ncalc))
    allocate(surf_elev(ncals))
    !$omp parallel
    !$omp workshare
    cell_top(:) = DZERO ; cell_cent(:) = DZERO ; cell_bot(:) = DZERO
    surf_elev(:) = DZERO
    !$omp end workshare

    ! make cell top and bottom elevation
    !$omp do private(i)
    do i = 1, ncals
      surf_elev(i) = fp_zt(i)
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, ncalc
      cell_top(i) = fp_zt(i) ; cell_bot(i) = fp_zb(i)
      cell_cent(i) = DHALF*(fp_zt(i)+fp_zb(i))
    end do
    !$omp end do
    !$omp end parallel

    deallocate(fp_zt, fp_zb)

  end subroutine make_cell_surfelev

  subroutine make_rech_area()
  !***************************************************************************************
  ! make_rech_area -- Make recharge area
  !***************************************************************************************
    ! -- modules
    use set_cell, only: get_cals_grid
    ! -- inout

    ! -- local
    integer(I4) :: i, xn, yn
    real(DP) :: x1, x2, x3, x4, y1, y2, y3, y4
    !-------------------------------------------------------------------------------------
    allocate(area_r(ncals))
    !$omp parallel
    !$omp workshare
    area_r(:) = DZERO
    !$omp end workshare

    ! make recharge area
    !$omp do private(i, xn, yn, x1, x2, x3, x4, y1, y2, y3, y4)
    do i = 1, ncals
      call get_cals_grid(i, xn, yn)
      x1 = glob_x(xn,yn) ; x2 = glob_x(xn+1,yn) ; x3 = glob_x(xn+1,yn+1)
      x4 = glob_x(xn,yn+1) ; y1 = glob_y(xn,yn) ; y2 = glob_y(xn+1,yn)
      y3 = glob_y(xn+1,yn+1) ; y4 = glob_y(xn,yn+1)

      area_r(i) = rarea(x1, x2, x3, x4, y1, y2, y3, y4)
    end do
    !$omp end do
    !$omp end parallel

  end subroutine make_rech_area

  function rarea(x1, x2, x3, x4, y1, y2, y3, y4) result(r_area)
  !***************************************************************************************
  ! rarea -- Calculate recharge area using Heron's formula
  !***************************************************************************************
    ! -- modules

    ! -- inout
    real(DP), intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4
    ! -- local
    real(DP) :: r_area, side1, side2, side3, side4, dialog, s1, s2
    real(DP) :: tri(2)
    !-------------------------------------------------------------------------------------
    side1 = sqrt((x1-x2)**2 + (y1-y2)**2) ; side2 = sqrt((x2-x3)**2 + (y2-y3)**2)
    side3 = sqrt((x3-x4)**2 + (y3-y4)**2) ; side4 = sqrt((x4-x1)**2 + (y4-y1)**2)
    dialog = sqrt((x1-x3)**2 + (y1-y3)**2)
    s1 = (side1+side2+dialog)*DHALF ; s2 = (side3+side4+dialog)*DHALF

    tri(1) = sqrt(s1*(s1-side1)*(s1-side2)*(s1-dialog))
    tri(2) = sqrt(s2*(s2-side3)*(s2-side4)*(s2-dialog))
    r_area = tri(1) + tri(2)

  end function rarea

end module make_cell
