module set_boundary
  ! -- modules
  use kind_module, only: I4, SP, DP
  use constval_module, only: SNOVAL, DNOVAL
  use types_module, only: rlbc_set, bfview_set, bound_fview
  use utility_module, only: st_mpi, write_logf, write_err_stop, get_ilen, conv_i2s
  use initial_module, only: in_type, st_in_type, st_rivf_type, st_lakf_type, st_in_path,&
                            st_in_unit
  use set_cell, only: ncals
  use set_condition, only: st_hydr, st_bcnd
  use assign_boundary, only: assign_surfbv, assign_rilav, st_forc
  use calc_boundary, only: conv_rech2calc, calc_wlbd, calc_blld, calc_lsurf
#ifdef MPI_MSG
  use mpi_utility, only: mpisum_val, bcast_file
  use mpi_set, only: cals_r4view, cals_r4hview
#endif

  implicit none
  private
  public :: set_bound

  type(rlbc_set), public :: st_rive, st_lake
#ifdef MPI_MSG
  type(bound_fview), public :: bfview
  type(bfview_set), public :: rfview, lfview
#endif

  ! -- local
  integer(I4) :: sum_riwln, sum_ribln, sum_riwdn, sum_riden, sum_riwin, sum_rilen
  integer(I4) :: sum_riarn, sum_lawln, sum_labln, sum_lawdn, sum_laarn

  contains

  subroutine set_bound()
  !*********************************************************************************************
  ! set_bound -- Set boundary
  !*********************************************************************************************
    ! -- modules
    use constval_module, only: DZERO
    use open_file, only: open_in_rivef, open_in_lakef
    use set_cell, only: ncalc
    use set_condition, only: set_connect, set_srabyd, set_chabyd, set_wellconn
    use calc_boundary, only: calc_reprev, count_rivecalc, count_lakecalc, calc_rivea
#ifdef MPI_MSG
   use mpi_set, only: bcast_bound_ftype, bcast_solval
#endif
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4) :: sum_rechn, sum_precn, sum_evapn
    integer(I4) :: rfv_wl, rfv_wd, rfv_bl, rfv_de, rfv_wi, rfv_le
    integer(I4) :: lfv_wl, lfv_wd, lfv_bl, lfv_ar
    character(:), allocatable :: num_str, err_mes
    !-------------------------------------------------------------------------------------------
#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Bcast boundary file type (bound_ftype)
        call bcast_bound_ftype()
    end if
#endif

    ! -- Set sea level information (seal_info)
      call set_seal_info()

    ! -- Set recharge information (rech_info)
      call set_rech_info()

    ! -- Set well information (well_info)
      call set_well_info()

    ! -- Set precipitation information (prec_info)
      call set_prec_info()

    ! -- Set evapotranspiration information (evap_info)
      call set_evap_info()

#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(st_bcnd%rech_num, "recharge", sum_rechn)
      call mpisum_val(st_bcnd%prec_num, "precipitation", sum_precn)
      call mpisum_val(st_bcnd%evap_num, "evapotranspiration", sum_evapn)
#else
    sum_rechn = st_bcnd%rech_num ; sum_precn = st_bcnd%prec_num ; sum_evapn = st_bcnd%evap_num
#endif

    if (sum_rechn == 0 .and. sum_precn /= 0 .and. sum_evapn /= 0) then
      ! -- Calculate recharge from precipitation and evapotranspiration (reprev)
        call calc_reprev(st_bcnd%rech_num)
        call conv_rech2calc(st_bcnd%rech_num)
      if (st_mpi%rank == 0) then
        call write_logf("Recharge is calculated from precipitation and evapotranspiration.")
        allocate(character(get_ilen(st_bcnd%rech_num)) :: num_str)
        allocate(character(0) :: err_mes)
        call conv_i2s(st_bcnd%rech_num, num_str)
        err_mes = "Set "//num_str//" recharge rate."
        call write_logf(err_mes)
        deallocate(num_str, err_mes)
      end if
    else if (st_mpi%rank == 0) then
      if (sum_rechn == 0 .and. sum_precn == 0 .and. sum_evapn /= 0) then
        call write_logf("Caution!! Specified only evapotranspiration in input file.")
      else if (sum_rechn == 0 .and. sum_precn /= 0 .and. sum_evapn == 0) then
        call write_logf("Caution!! Specified only precipitation in input file.")
      end if
    end if

    if (st_in_type%rive == in_type(0)) then
#ifdef MPI_MSG
      ! -- Read input river file (inrivef)
        call open_in_rivef(st_in_path%rive, cals_r4view, cals_r4hview, rfv_wl, rfv_wd, rfv_bl,&
                           rfv_de, rfv_wi, rfv_le)
      rfview%wl = rfv_wl ; rfview%wd = rfv_wd ; rfview%bl = rfv_bl ; rfview%de = rfv_de
      rfview%wi = rfv_wi ; rfview%le = rfv_le
#else
      ! -- Read input river file (inrivef)
        call open_in_rivef(st_in_path%rive, 0, 0, rfv_wl, rfv_wd, rfv_bl, rfv_de, rfv_wi,&
                           rfv_le)
#endif
    end if

    ! -- Set river water level information (riwl_info)
      call set_riwl_info()

    ! -- Set river bottom level information (ribl_info)
      call set_ribl_info()

#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(st_rive%num%wl, "river water level", sum_riwln)
      call mpisum_val(st_rive%num%bl, "river bottom level", sum_ribln)
#else
    sum_riwln = st_rive%num%wl ; sum_ribln = st_rive%num%bl
#endif

    if (sum_riwln == 0 .or. sum_ribln == 0) then
      ! -- Set river water depth information (riwd_info)
        call set_riwd_info()
#ifdef MPI_MSG
      ! -- Sum value for MPI (val)
        call mpisum_val(st_rive%num%wd, "river water depth", sum_riwdn)
#else
      sum_riwdn = st_rive%num%wd
#endif
    end if

    if (sum_ribln == 0) then
      ! -- Set river depth information (ride_info)
        call set_ride_info()

#ifdef MPI_MSG
      ! -- Sum value for MPI (val)
        call mpisum_val(st_rive%num%de, "river depth", sum_riden)
#else
      sum_riden = st_rive%num%de
#endif
      ! -- Set river bottom level (rive_bott)
        call set_rive_bott()

#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(st_rive%num%bl, "river bottom level", sum_ribln)
#else
    sum_ribln = st_rive%num%bl
#endif
    end if

    ! -- Set river water level (rive_wlevel)
      call set_rive_wlevel()

    ! -- Set river width information (riwi_info)
      call set_riwi_info()

    ! -- Set river length information (rile_info)
      call set_rile_info()

#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(st_rive%num%wi, "river width", sum_riwin)
      call mpisum_val(st_rive%num%le, "river length", sum_rilen)
#else
    sum_riwin = st_rive%num%wi ; sum_rilen = st_rive%num%le
#endif

    if (sum_riwin > 0 .and. sum_rilen > 0) then
      allocate(st_rive%cflag%ar(ncals))
      allocate(st_rive%calc%ar(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_rive%cflag%ar(i) = 0
        st_rive%calc%ar(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Calculate river area (rivea)
        call calc_rivea(st_rive%cflag%wi, st_rive%cflag%le, st_rive%calc%wi, st_rive%calc%le,&
                        st_rive%cflag%ar, st_rive%calc%ar, st_rive%num%ar)
    end if

#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(st_rive%num%ar, "river area", sum_riarn)
#else
    sum_riarn = st_rive%num%ar
#endif

    if (sum_riwln /= 0 .and. sum_ribln /= 0 .and. sum_riarn /= 0) then
      ! -- Count river calculation (rivecalc)
        call count_rivecalc(st_rive%cflag%wl, st_rive%cflag%bl, st_rive%cflag%ar,&
                            st_rive%calc%wl, st_rive%calc%bl, st_rive%calc%ar,&
                            st_bcnd%rive_num)
    end if

    st_rive%num%wl = 0

    if (st_in_type%lake == in_type(0)) then
#ifdef MPI_MSG
      ! -- Open input lake file (in_lakef)
        call open_in_lakef(st_in_path%lake, cals_r4view, cals_r4hview, lfv_wl, lfv_wd, lfv_bl,&
                           lfv_ar)
      lfview%wl = lfv_wl ; lfview%wd = lfv_wd ; lfview%bl = lfv_bl ; lfview%ar = lfv_ar
#else
      ! -- Open input lake file (in_lakef)
        call open_in_lakef(st_in_path%lake, 0, 0, lfv_wl, lfv_wd, lfv_bl, lfv_ar)
#endif
    end if

    ! -- Set lake water level information (lawl_info)
      call set_lawl_info()

    ! -- Set lake bottom level information (labl_info)
      call set_labl_info()

#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(st_lake%num%wl, "lake water level", sum_lawln)
      call mpisum_val(st_lake%num%bl, "lake bottom level", sum_labln)
#else
    sum_lawln = st_lake%num%wl ; sum_labln = st_lake%num%bl
#endif

    if (sum_lawln == 0 .or. sum_labln == 0) then
      ! -- Set lake water depth information (lawd_info)
        call set_lawd_info()
#ifdef MPI_MSG
      ! -- Sum value for MPI (val)
        call mpisum_val(st_lake%num%wd, "lake water depth", sum_lawdn)
#else
      sum_lawdn = st_lake%num%wd
#endif
    end if

    ! -- Set lake water or bottom level (lake_wblevel)
      call set_lake_wblevel()

    ! -- Set lake area information (laar_info)
      call set_laar_info()

#ifdef MPI_MSG
    ! -- Sum value for MPI (val)
      call mpisum_val(st_lake%num%ar, "lake area", sum_laarn)
#else
    sum_laarn = st_lake%num%ar
#endif

    if (sum_lawln /= 0 .and. sum_labln /= 0 .and. sum_laarn /= 0) then
      ! -- Count lake calculation cell (lakecalc)
        call count_lakecalc(st_lake%cflag%wl, st_lake%cflag%bl, st_lake%cflag%ar,&
                            st_lake%calc%wl, st_lake%calc%bl, st_lake%calc%ar,&
                            st_bcnd%lake_num)
    end if

    st_lake%num%wl = 0

    ! -- Set connectivity (connect)
      call set_connect(st_hydr%read_hydx, st_hydr%read_hydy, st_hydr%read_hydz)

    if (st_bcnd%well_num /= 0) then
      ! -- Set well connectivity (wellconn)
        call set_wellconn(st_bcnd%well_num, st_hydr%read_hydx, st_hydr%read_hydy)
    end if

    if (st_bcnd%rive_num /= 0) then
      allocate(st_forc%abyd_rive(st_bcnd%rive_num))
      !$omp parallel do private(i)
      do i = 1, st_bcnd%rive_num
        st_forc%abyd_rive(i) = DZERO
      end do
      !$omp end parallel do
      ! -- Set surface&recharge area and area by distance (srabyd)
        call set_srabyd(st_bcnd%rive_num, st_forc%rive_bott, st_forc%rive_area,&
                        st_bcnd%rive2cals, st_forc%abyd_rive)
    end if

    if (st_bcnd%lake_num /= 0) then
      allocate(st_forc%abyd_lake(st_bcnd%lake_num))
      !$omp parallel do private(i)
      do i = 1, st_bcnd%lake_num
        st_forc%abyd_lake(i) = DZERO
      end do
      !$omp end parallel do
      ! -- Set surface&recharge area and area by distance (srabyd)
        call set_srabyd(st_bcnd%lake_num, st_forc%lake_bott, st_forc%lake_area,&
                        st_bcnd%lake2cals, st_forc%abyd_lake)
    end if

    ! -- Set charge area by distance (chabyd)
      call set_chabyd()

    allocate(st_forc%read_head(ncalc))
    !$omp parallel
    !$omp do private(i)
    do i = 1, ncalc
      st_forc%read_head(:) = DZERO
    end do
    !$omp end do

    !$omp do private(i)
    do i = 1, ncalc
      st_forc%read_head(i) = st_hydr%read_init(i)
    end do
    !$omp end do
    !$omp end parallel

    deallocate(st_hydr%read_init)

#ifdef MPI_MSG
    if (st_mpi%totn /= 1) then
      ! -- Bcast solution value (solval)
        call bcast_solval()
    end if
#endif

  end subroutine set_bound

  subroutine set_seal_info()
  !*********************************************************************************************
  ! set_seal_info -- Set sea level information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_seal
    use open_file, only: open_in_sealf
    use assign_boundary, only: assign_sealv
#ifdef MPI_MSG
    use open_file, only: st_intse
    use mpi_set, only: surf_r4view, surf_r4hview, cell_r4view, cell_r4hview
#endif
    ! -- inout

    ! -- local
    integer(I4), allocatable :: all_seal_type(:)
    logical, allocatable :: all_seal_mask(:)
    !-------------------------------------------------------------------------------------------
    allocate(all_seal_type(7), all_seal_mask(7))
    all_seal_type(:) = [in_type(1:7)]
    all_seal_mask(:) = (st_in_type%seal == all_seal_type(:))

    if (any(all_seal_mask)) then
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Bcast file (file)
          call bcast_file(st_in_path%seal, st_in_unit%seal, "sea level")
      end if

      if (st_in_type%seal == in_type(4)) then
        if (len_trim(adjustl(st_in_unit%seal)) == 0) then
            bfview%seal = surf_r4view ; st_intse%type = 0
        else
            bfview%seal = surf_r4hview ; st_intse%type = st_in_type%seal
        end if
      else if (st_in_type%seal == in_type(6)) then
        if (len_trim(adjustl(st_in_unit%seal)) == 0) then
            bfview%seal = cell_r4view
        else
            bfview%seal = cell_r4hview
        end if
      else if (st_in_type%seal == in_type(7)) then
        if (st_intse%type == in_type(4)) then
          bfview%seal = surf_r4view
        else if (st_intse%type == in_type(6)) then
          bfview%seal = cell_r4view
        end if
      end if

      ! -- Open input sea level file (in_sealf)
        call open_in_sealf(st_in_type%seal, st_in_path%seal, st_in_unit%seal, bfview%seal)
#else
      ! -- Open input sea level file (in_sealf)
        call open_in_sealf(st_in_type%seal, st_in_path%seal, st_in_unit%seal)
#endif

    else
      st_seal%totn = 0
      if (st_mpi%rank == 0) then
        call write_logf("Set closed boundary problem.")
      end if
    end if

    ! -- Assign sea level value (sealv)
      call assign_sealv(st_in_type%seal)

    deallocate(all_seal_type, all_seal_mask)

  end subroutine set_seal_info

  subroutine set_rech_info()
  !*********************************************************************************************
  ! set_rech_info -- Set recharge information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_rech
    use open_file, only: open_in_rechf, st_intre
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: all_rech_type(:)
    logical, allocatable :: all_rech_mask(:)
    !-------------------------------------------------------------------------------------------
    allocate(all_rech_type(4), all_rech_mask(4))
    all_rech_type(:) = [in_type(1), in_type(3:4), in_type(7)]
    all_rech_mask(:) = (st_in_type%rech == all_rech_type(:))

    if (any(all_rech_mask)) then
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Bcast file (file)
          call bcast_file(st_in_path%rech, st_in_unit%rech, "recharge")
      end if

      if (st_in_type%rech == in_type(4)) then
        if (len_trim(adjustl(st_in_unit%rech)) == 0) then
          bfview%rech = cals_r4view ; st_intre%type = 0
        else
          bfview%rech = cals_r4hview ; st_intre%type = st_in_type%rech
        end if
      else if (st_in_type%rech == in_type(7)) then
        bfview%rech = cals_r4view
      end if

      ! -- Open input recharge file (in_rechf)
        call open_in_rechf(st_in_type%rech, st_in_path%rech, st_in_unit%rech, bfview%rech)
#else
      ! -- Open input recharge file (in_rechf)
        call open_in_rechf(st_in_type%rech, st_in_path%rech, st_in_unit%rech)
#endif
    end if

    if (st_rech%totn > 0) then
      allocate(st_bcnd%rech_cflag(ncals))
      allocate(st_forc%read_rech(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_bcnd%rech_cflag(i) = 0
        st_forc%read_rech(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign recharge value
        call assign_surfbv(st_in_type%rech, st_intre%type, st_rech, st_bcnd%rech_num,&
                           st_bcnd%rech_cflag,&
                           st_forc%read_rech)

      call conv_rech2calc(st_bcnd%rech_num)
    end if

    deallocate(all_rech_type, all_rech_mask)

  end subroutine set_rech_info

  subroutine set_well_info()
  !*********************************************************************************************
  ! set_well_info -- Set well information
  !*********************************************************************************************
    ! -- modules
    use open_file, only: open_in_wellf, open_in_wlayf
    use assign_boundary, only: assign_wellv
#ifdef MPI_MSG
    use open_file, only: st_intwe
    use mpi_set, only: cals_i4view, calc_r4view, calc_r4hview
#endif
    ! -- inout

    ! -- local
    integer(I4), allocatable :: all_well_type(:)
    logical, allocatable :: all_well_mask(:)
    !-------------------------------------------------------------------------------------------
    allocate(all_well_type(6), all_well_mask(6))
    all_well_type(:) = [in_type(2:7)]
    all_well_mask(:) = (st_in_type%well == all_well_type(:))

    if (any(all_well_mask)) then
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Bcast file (file)
          call bcast_file(st_in_path%well, st_in_unit%well, "well")
      end if

      if (st_in_type%well == in_type(4)) then
        if (len_trim(adjustl(st_in_unit%well)) == 0) then
          bfview%well = cals_r4view ; st_intwe%type = 0
        else
          bfview%well = cals_r4hview ; st_intwe%type = st_in_type%well
        end if
      else if (st_in_type%well == in_type(6)) then
        if (len_trim(adjustl(st_in_unit%well)) == 0) then
          bfview%well = calc_r4view ; st_intwe%type = 0
        else
          bfview%well = calc_r4hview ; st_intwe%type = st_in_type%well
        end if
      else if (st_in_type%well == in_type(7)) then
        bfview%well = calc_r4view
      end if

      ! -- Open input well file (in_wellf)
        call open_in_wellf(st_in_type%well, st_in_path%well, st_in_unit%well, bfview%well)
#else
      ! -- Open input well file (in_wellf)
        call open_in_wellf(st_in_type%well, st_in_path%well, st_in_unit%well)
#endif

      if (st_in_type%well == in_type(3) .or. st_in_type%well == in_type(4)) then
#ifdef MPI_MSG
        ! -- Open input well layer file (in_wlayf)
          call open_in_wlayf(st_in_type%weks, st_in_type%weke, st_in_path%weks,&
                             st_in_path%weke, cals_i4view)
#else
        ! -- Open input well layer file (in_wlayf)
          call open_in_wlayf(st_in_type%weks, st_in_type%weke, st_in_path%weks, st_in_path%weke)
#endif
      end if
    end if

    st_bcnd%well_num = 0
    ! -- Assign well value (wellv)
      call assign_wellv(st_in_type%well, st_in_type%weks, st_in_type%weke, st_bcnd%well_num)

    deallocate(all_well_type, all_well_mask)

  end subroutine set_well_info

  subroutine set_prec_info()
  !*********************************************************************************************
  ! set_prec_info -- Set precipitation information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_prec
    use open_file, only: open_in_precf, st_intpr
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: all_prec_type(:)
    logical, allocatable :: all_prec_mask(:)
    !-------------------------------------------------------------------------------------------
    allocate(all_prec_type(4), all_prec_mask(4))
    all_prec_type(:) = [in_type(1), in_type(3:4), in_type(7)]
    all_prec_mask(:) = (st_in_type%prec == all_prec_type(:))

    if (any(all_prec_mask)) then
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Bcast file (file)
          call bcast_file(st_in_path%prec, st_in_unit%prec, "precipitation")
      end if

      if (st_in_type%prec == in_type(4)) then
        if (len_trim(adjustl(st_in_unit%prec)) == 0) then
          bfview%prec = cals_r4view ; st_intpr%type = 0
        else
          bfview%prec = cals_r4hview ; st_intpr%type = st_in_type%prec
        end if
      else if (st_in_type%prec == in_type(7)) then
        bfview%prec = cals_r4view
      end if

      ! -- Open input precipitation file (in_precf)
        call open_in_precf(st_in_type%prec, st_in_path%prec, st_in_unit%prec, bfview%prec)
#else
      ! -- Open input precipitation file (in_precf)
        call open_in_precf(st_in_type%prec, st_in_path%prec, st_in_unit%prec)
#endif
    end if

    if (st_prec%totn > 0) then
      allocate(st_bcnd%prec_cflag(ncals))
      allocate(st_forc%read_prec(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_bcnd%prec_cflag(i) = 0
        st_forc%read_prec(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign precipitation value
        call assign_surfbv(st_in_type%prec, st_intpr%type, st_prec, st_bcnd%prec_num,&
                           st_bcnd%prec_cflag,&
                           st_forc%read_prec)
      deallocate(st_bcnd%prec_cflag)
    end if

    deallocate(all_prec_type, all_prec_mask)

  end subroutine set_prec_info

  subroutine set_evap_info()
  !*********************************************************************************************
  ! set_evap_info -- Set evapotranspiration information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_evap
    use open_file, only: open_in_evapf, st_intev
    ! -- inout

    ! -- local
    integer(I4) :: i
    integer(I4), allocatable :: all_evap_type(:)
    logical, allocatable :: all_evap_mask(:)
    !-------------------------------------------------------------------------------------------
    allocate(all_evap_type(4), all_evap_mask(4))
    all_evap_type(:) = [in_type(1), in_type(3:4), in_type(7)]
    all_evap_mask(:) = (st_in_type%evap == all_evap_type(:))

    if (any(all_evap_mask)) then
#ifdef MPI_MSG
      if (st_mpi%totn /= 1) then
        ! -- Bcast file (file)
          call bcast_file(st_in_path%evap, st_in_unit%evap, "evapotranspiration")
      end if

      if (st_in_type%evap == in_type(4)) then
        if (len_trim(adjustl(st_in_unit%evap)) == 0) then
          bfview%evap = cals_r4view ; st_intev%type = 0
        else
          bfview%evap = cals_r4hview ; st_intev%type = st_in_type%evap
        end if
      else if (st_in_type%evap == in_type(7)) then
        bfview%evap = cals_r4view
      end if

      ! -- Open input evapotranspiration file (in_evapf)
        call open_in_evapf(st_in_type%evap, st_in_path%evap, st_in_unit%evap, bfview%evap)
#else
      ! -- Open input evapotranspiration file (in_evapf)
        call open_in_evapf(st_in_type%evap, st_in_path%evap, st_in_unit%evap)
#endif
    end if

    if (st_evap%totn > 0) then
      allocate(st_bcnd%evap_cflag(ncals))
      allocate(st_forc%read_evap(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_bcnd%evap_cflag(i) = 0
        st_forc%read_evap(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign evapotranspiration value
        call assign_surfbv(st_in_type%evap, st_intev%type, st_evap, st_bcnd%evap_num,&
                           st_bcnd%evap_cflag,&
                           st_forc%read_evap)
      deallocate(st_bcnd%evap_cflag)
    end if

    deallocate(all_evap_type, all_evap_mask)

  end subroutine set_evap_info

  subroutine set_riwl_info()
  !*********************************************************************************************
  ! set_riwl_info -- Set river water level information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_riwl
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    allocate(st_rive%cflag%wl(ncals))
    allocate(st_rive%calc%wl(ncals))
    !$omp parallel do private(i)
    do i = 1, ncals
      st_rive%cflag%wl(i) = 0
      st_rive%calc%wl(i) = SNOVAL
    end do
    !$omp end parallel do
    ! -- Assign river water level value
      call assign_rilav(st_rivf_type%wlev, 0, st_riwl, st_rive%num%wl, st_rive%cflag%wl,&
                        st_rive%calc%wl)

  end subroutine set_riwl_info

  subroutine set_ribl_info()
  !*********************************************************************************************
  ! set_ribl_info -- Set river bottom level information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_ribl
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    allocate(st_rive%cflag%bl(ncals))
    allocate(st_rive%calc%bl(ncals))
    !$omp parallel do private(i)
    do i = 1, ncals
      st_rive%cflag%bl(i) = 0
      st_rive%calc%bl(i) = SNOVAL
    end do
    !$omp end parallel do
    ! -- Assign river bottom level value
      call assign_rilav(st_rivf_type%blev, 0, st_ribl, st_rive%num%bl, st_rive%cflag%bl,&
                        st_rive%calc%bl)

  end subroutine set_ribl_info

  subroutine set_riwd_info()
  !*********************************************************************************************
  ! set_riwd_info -- Set river water depth information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_riwd
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    if (st_riwd%totn > 0) then
      allocate(st_rive%cflag%wd(ncals))
      allocate(st_rive%calc%wd(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_rive%cflag%wd(i) = 0
        st_rive%calc%wd(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign river water depth value
        call assign_rilav(st_rivf_type%wdep, 0, st_riwd, st_rive%num%wd, st_rive%cflag%wd,&
                          st_rive%calc%wd)
    end if

  end subroutine set_riwd_info

  subroutine set_ride_info()
  !*********************************************************************************************
  ! set_ride_info -- Set river depth information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_ride
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    if (st_ride%totn > 0) then
      allocate(st_rive%cflag%de(ncals))
      allocate(st_rive%calc%de(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_rive%cflag%de(i) = 0
        st_rive%calc%de(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign river depth value
        call assign_rilav(st_rivf_type%dept, 0, st_ride, st_rive%num%de, st_rive%cflag%de,&
                          st_rive%calc%de)
    end if

  end subroutine set_ride_info

  subroutine set_riwi_info()
  !*********************************************************************************************
  ! set_riwi_info -- Set river width information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_riwi
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    if (st_riwi%totn > 0) then
      allocate(st_rive%cflag%wi(ncals))
      allocate(st_rive%calc%wi(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_rive%cflag%wi(i) = 0
        st_rive%calc%wi(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign river width value
        call assign_rilav(st_rivf_type%widt, 0, st_riwi, st_rive%num%wi, st_rive%cflag%wi,&
                          st_rive%calc%wi)
    end if

  end subroutine set_riwi_info

  subroutine set_rile_info()
  !*********************************************************************************************
  ! set_rile_info -- Set river length information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_rile
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    if (st_rile%totn > 0) then
      allocate(st_rive%cflag%le(ncals))
      allocate(st_rive%calc%le(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_rive%cflag%le(i) = 0
        st_rive%calc%le(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign river length value
        call assign_rilav(st_rivf_type%leng, 0, st_rile, st_rive%num%le, st_rive%cflag%le,&
                          st_rive%calc%le)
    end if

  end subroutine set_rile_info

  subroutine set_lawl_info()
  !*********************************************************************************************
  ! set_lawl_info -- Set lake water level information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_lawl
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    allocate(st_lake%cflag%wl(ncals))
    allocate(st_lake%calc%wl(ncals))
    !$omp parallel do private(i)
    do i = 1, ncals
      st_lake%cflag%wl(i) = 0
      st_lake%calc%wl(i) = SNOVAL
    end do
    !$omp end parallel do
    ! -- Assign lake water level value
      call assign_rilav(st_lakf_type%wlev, 0, st_lawl, st_lake%num%wl, st_lake%cflag%wl,&
                        st_lake%calc%wl)

  end subroutine set_lawl_info

  subroutine set_labl_info()
  !*********************************************************************************************
  ! set_labl_info -- Set lake bottom level information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_labl
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    allocate(st_lake%cflag%bl(ncals))
    allocate(st_lake%calc%bl(ncals))
    !$omp parallel do private(i)
    do i = 1, ncals
      st_lake%cflag%bl(i) = 0
      st_lake%calc%bl(i) = SNOVAL
    end do
    !$omp end parallel do
    ! -- Assign lake bottom level value
      call assign_rilav(st_lakf_type%blev, 0, st_labl, st_lake%num%bl, st_lake%cflag%bl,&
                        st_lake%calc%bl)

  end subroutine set_labl_info

  subroutine set_lawd_info()
  !*********************************************************************************************
  ! set_lawd_info -- Set lake water depth information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_lawd
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    if (st_lawd%totn > 0) then
      allocate(st_lake%cflag%wd(ncals))
      allocate(st_lake%calc%wd(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_lake%cflag%wd(i) = 0
        st_lake%calc%wd(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign lake water depth value
        call assign_rilav(st_lakf_type%wdep, 0, st_lawd, st_lake%num%wd, st_lake%cflag%wd,&
                          st_lake%calc%wd)
    end if

  end subroutine set_lawd_info

  subroutine set_laar_info()
  !*********************************************************************************************
  ! set_laar_info -- Set lake area information
  !*********************************************************************************************
    ! -- modules
    use initial_module, only: st_laar
    ! -- inout

    ! -- local
    integer(I4) :: i
    !-------------------------------------------------------------------------------------------
    if (st_laar%totn > 0) then
      allocate(st_lake%cflag%ar(ncals))
      allocate(st_lake%calc%ar(ncals))
      !$omp parallel do private(i)
      do i = 1, ncals
        st_lake%cflag%ar(i) = 0
        st_lake%calc%ar(i) = SNOVAL
      end do
      !$omp end parallel do
      ! -- Assign lake area value
        call assign_rilav(st_lakf_type%area, 1, st_laar, st_lake%num%ar, st_lake%cflag%ar,&
                          st_lake%calc%ar)
    end if

  end subroutine set_laar_info

  subroutine set_rive_bott()
  !*********************************************************************************************
  ! set_rive_bott -- Set river bottom
  !*********************************************************************************************
    ! -- modules
    use calc_boundary, only: calc_blsl
    ! -- inout

    ! -- local
    character(:), allocatable :: err_mes
    !-------------------------------------------------------------------------------------------
    if (sum_riden /= 0) then
      ! -- Calculate bottom level from surface level (blsl)
        call calc_blsl(st_rive%cflag%de, st_rive%calc%de, st_rive%cflag%bl, st_rive%calc%bl,&
                       st_rive%num%bl)
      deallocate(st_rive%cflag%de, st_rive%calc%de)
      if (st_mpi%rank == 0) then
        allocate(character(0) :: err_mes)
        err_mes = "River bottom level is calculated from surface elevation and river depth."
        call write_logf(err_mes)
      end if
    else if (sum_riwln /= 0 .and. sum_riwdn /= 0) then
      ! -- Calculate bottom level from water level and water depth (blld)
        call calc_blld(st_rive%cflag%wl, st_rive%calc%wl, st_rive%cflag%wd, st_rive%calc%wd,&
                       st_rive%cflag%bl, st_rive%calc%bl, st_rive%num%bl)
      deallocate(st_rive%cflag%wd, st_rive%calc%wd)
      if (st_mpi%rank == 0) then
        allocate(character(0) :: err_mes)
        err_mes = "River bottom level is calculated from water level and water depth."
        call write_logf(err_mes)
      end if
    else if (sum_riwln /= 0 .and. sum_riwdn == 0) then
      ! -- Calculate level from surface (lsurf)
        call calc_lsurf(st_rive%cflag%wl, st_rive%cflag%bl, st_rive%calc%bl, st_rive%num%bl)
      if (st_mpi%rank == 0) then
        allocate(character(0) :: err_mes)
        err_mes = "River bottom level is setted to surface elevation."
        call write_logf(err_mes)
      end if
    else if (sum_riwln == 0 .and. sum_riwdn /= 0 .and. st_mpi%rank == 0) then
      allocate(character(0) :: err_mes)
      err_mes = "Only specified river water depth."
      call write_err_stop(err_mes)
    end if

    if (allocated(err_mes)) then
      deallocate(err_mes)
    end if

  end subroutine set_rive_bott

  subroutine set_rive_wlevel()
  !*********************************************************************************************
  ! set_rive_wlevel -- Set river water level
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    character(:), allocatable :: num_str, err_mes
    !-------------------------------------------------------------------------------------------
    if (sum_riwln == 0 .and. sum_ribln /= 0 .and. sum_riwdn /= 0) then
      ! -- Calculate water level from bottom level and water depth (wlbd)
        call calc_wlbd(st_rive%cflag%bl, st_rive%calc%bl, st_rive%cflag%wd, st_rive%calc%wd,&
                       st_rive%cflag%wl, st_rive%calc%wl, st_rive%num%wl)
      deallocate(st_rive%cflag%wd, st_rive%calc%wd)
#ifdef MPI_MSG
      ! -- Sum value for MPI (val)
        call mpisum_val(st_rive%num%wl, "river water level", sum_riwln)
#else
      sum_riwln = st_rive%num%wl
#endif
      if (st_mpi%rank == 0) then
        call write_logf("River water level is calculated from bottom level.")
        allocate(character(get_ilen(sum_riwln)) :: num_str)
        allocate(character(0) :: err_mes)
        call conv_i2s(sum_riwln, num_str)
        err_mes = "Set "//num_str//" river water level."
        call write_logf(err_mes)
      end if
    else if (sum_riwln == 0 .and. sum_ribln /= 0 .and. sum_riden /= 0) then
      if (st_mpi%rank == 0) then
        allocate(character(0) :: err_mes)
        err_mes = "Not calculated river water level from river bottom level and river depth."
        call write_err_stop(err_mes)
      end if
    else if (sum_riwln == 0 .and. sum_ribln /= 0 .and. sum_riden == 0) then
      if (st_mpi%rank == 0) then
        allocate(character(0) :: err_mes)
        err_mes = "Not calculated river water level from only river bottom level."
        call write_err_stop(err_mes)
      end if
    end if

    if (allocated(err_mes)) then
      deallocate(err_mes)
    end if

    if (allocated(num_str)) then
      deallocate(num_str)
    end if

  end subroutine set_rive_wlevel

  subroutine set_lake_wblevel()
  !*********************************************************************************************
  ! set_lake_wblevel -- Set lake water or bottom level
  !*********************************************************************************************
    ! -- modules

    ! -- inout

    ! -- local
    character(:), allocatable :: num_str, err_mes
    !-------------------------------------------------------------------------------------------
    if (sum_lawln == 0 .and. sum_labln /= 0 .and. sum_lawdn /= 0) then
      ! -- Calculate water level from bottom level and water depth (wlbd)
        call calc_wlbd(st_lake%cflag%bl, st_lake%calc%bl, st_lake%cflag%wd, st_lake%calc%wd,&
                       st_lake%cflag%wl, st_lake%calc%wl, st_lake%num%wl)
      deallocate(st_lake%cflag%wd, st_lake%calc%wd)
      if (st_mpi%rank == 0) then
        call write_logf("Lake water level is calculated from bottom level.")
        allocate(character(get_ilen(sum_lawln)) :: num_str)
        allocate(character(0) :: err_mes)
        call conv_i2s(sum_lawln, num_str)
        err_mes = "Set "//num_str//" lake water level."
        call write_logf(err_mes)
      end if
    else if (sum_lawln /= 0 .and. sum_labln == 0 .and. sum_lawdn /= 0) then
      ! -- Calculate bottom level from water level and water depth (blld)
        call calc_blld(st_lake%cflag%wl, st_lake%calc%wl, st_lake%cflag%wd, st_lake%calc%wd,&
                       st_lake%cflag%bl, st_lake%calc%bl, st_lake%num%bl)
      deallocate(st_lake%cflag%wd, st_lake%calc%wd)
      if (st_mpi%rank == 0) then
        call write_logf("Lake bottom level is calculated from water level and water depth.")
        allocate(character(get_ilen(sum_labln)) :: num_str)
        allocate(character(0) :: err_mes)
        call conv_i2s(sum_labln, num_str)
        err_mes = "Set "//num_str//" lake bottom level."
        call write_logf(err_mes)
      end if
    else if (sum_lawln == 0 .and. sum_labln == 0 .and. sum_lawdn /= 0) then
      if (st_mpi%rank == 0) then
        allocate(character(0) :: err_mes)
        err_mes = "Only specified lake water depth."
        call write_err_stop(err_mes)
      end if
    else if (sum_lawln == 0 .and. sum_labln /= 0 .and. sum_lawdn == 0) then
      ! -- Set level from surface (levsurf)
        call calc_lsurf(st_lake%cflag%bl, st_lake%cflag%wl, st_lake%calc%wl, st_lake%num%wl)
      if (st_mpi%rank == 0) then
        call write_logf("Lake water level is setted to surface elevation.")
        allocate(character(get_ilen(sum_lawln)) :: num_str)
        allocate(character(0) :: err_mes)
        call conv_i2s(sum_lawln, num_str)
        err_mes = "Set "//num_str//" lake water level."
        call write_logf(err_mes)
      end if
    else if (sum_lawln /= 0 .and. sum_labln == 0 .and. sum_lawdn == 0) then
      if (st_mpi%rank == 0) then
        allocate(character(0) :: err_mes)
        err_mes = "Not calculated river bottom level."
        call write_err_stop(err_mes)
      end if
    end if

    if (allocated(err_mes)) then
      deallocate(err_mes)
    end if

    if (allocated(num_str)) then
      deallocate(num_str)
    end if

  end subroutine set_lake_wblevel

end module set_boundary
