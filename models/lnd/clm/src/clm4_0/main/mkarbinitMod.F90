module mkarbinitMod
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: mkarbinitMod
!
! !DESCRIPTION:
!
!
!---------------------------------------------------------------------------

! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_varctl   , only : iulog, use_cn, use_cndv, &
                              use_vancouver, use_mexicocity
    use shr_sys_mod  , only : shr_sys_flush
    use spmdMod      , only : masterproc

    implicit none

    SAVE
    private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

    public mkarbinit   ! Make arbitrary initial conditions

!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkarbinit
!
! !INTERFACE:
  subroutine mkarbinit()
!
! !DESCRIPTION:
! Initializes the following time varying variables:
! water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
! snow       : snowdp, snl, dz, z, zi
! temperature: t_soisno, t_veg, t_grnd
!
! !USES:
    use shr_const_mod, only : SHR_CONST_TKFRZ
    use clmtype
    use clm_varpar   , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use clm_varcon   , only : bdsno, istice, istwet, istsoil, isturb, &
                              denice, denh2o, spval, sb, icol_road_perv, &
                              icol_road_imperv, icol_roof, icol_sunwall, &
                              icol_shadewall
    use clm_varcon   , only : istcrop
    use clm_varcon   , only : istice_mec, h2osno_max
    use clm_varctl   , only : iulog
    use spmdMod      , only : masterproc
    use decompMod    , only : get_proc_bounds
    use SNICARMod    , only : snw_rds_min
    use HydrologyTracer , only : pwtrc, Rstnd, TracerCheckEqual, pdbg, mdbg, cdbg

!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 3/07/08 Keith Oleson: initialize h2osoi_vol for all soil layers to 0.3
! 3/18/08 David Lawrence, initialize deep layers
! 03/28/08 Mark Flanner, initialize snow aerosols and grain size
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)        ! column index associated with each pft
    integer , pointer :: ctype(:)          ! column type
    integer , pointer :: clandunit(:)      ! landunit index associated with each column
    integer , pointer :: ltype(:)          ! landunit type
    logical , pointer :: lakpoi(:)         ! true => landunit is a lake point
    integer , pointer :: plandunit(:)      ! landunit index associated with each pft
    logical , pointer :: urbpoi(:)         ! true => landunit is an urban point
    logical , pointer :: ifspecial(:)      ! true => landunit is not vegetated
    real(r8), pointer :: dz(:,:)           ! layer thickness depth (m)
    real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(r8), pointer :: bsw2(:,:)         ! Clapp and Hornberger "b" for CN code
    real(r8), pointer :: psisat(:,:)       ! soil water potential at saturation for CN code (MPa)
    real(r8), pointer :: vwcsat(:,:)       ! volumetric water content at saturation for CN code (m3/m3)
    real(r8), pointer :: zi(:,:)           ! interface level below a "z" level (m)
    real(r8), pointer :: wa(:)             ! water in the unconfined aquifer (mm)
    real(r8), pointer :: wt(:)             ! total water storage (unsaturated soil water + groundwater) (mm)
    real(r8), pointer :: zwt(:)            ! water table depth (m)
    real(r8), pointer :: qflx_snow_melt(:)  ! snow melt (net)

    real(r8), pointer :: wtr_h2osoi_ice(:,:,:)   ! tracer ice lens (kg/m2)
    real(r8), pointer :: wtr_h2osoi_liq(:,:,:)   ! tracer liquid water (kg/m2)
    real(r8), pointer :: wtr_wa(:,:)       ! tracer water in the unconfined aquifer (mm)
    real(r8), pointer :: wtr_wt(:,:)       ! tracer total water storage (unsaturated soil water + groundwater) (mm)
    real(r8), pointer :: wtr_qflx_snow_melt(:,:) ! tracer snow melt (net)
!
! local pointers to implicit out arguments
!
    integer , pointer :: snl(:)             ! number of snow layers
    real(r8), pointer :: t_soisno(:,:)      ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: t_lake(:,:)        ! lake temperature (Kelvin)  (1:nlevlak)
    real(r8), pointer :: t_grnd(:)          ! ground temperature (Kelvin)
    real(r8), pointer :: t_veg(:)           ! vegetation temperature (Kelvin)
    real(r8), pointer :: t_ref2m(:)         ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_ref2m_u(:)       ! Urban 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: t_ref2m_r(:)       ! Rural 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: h2osoi_vol(:,:)    ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: h2ocan_col(:)      ! canopy water (mm H2O) (column-level)
    real(r8), pointer :: h2ocan_pft(:)      ! canopy water (mm H2O) (pft-level)
    real(r8), pointer :: h2osno(:)          ! snow water (mm H2O)
    real(r8), pointer :: snowdp(:)          ! snow height (m)
    real(r8), pointer :: qflx_irrig(:)      ! irrigation flux (mm H2O/s)
    real(r8), pointer :: eflx_lwrad_out(:)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: soilpsi(:,:)       ! soil water potential in each soil layer (MPa)
    real(r8), pointer :: snw_rds(:,:)       ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(r8), pointer :: snw_rds_top(:)     ! snow grain size, top (col) [microns]
    real(r8), pointer :: sno_liq_top(:)     ! liquid water fraction (mass) in top snow layer (col) [frc]
    real(r8), pointer :: mss_bcpho(:,:)     ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcphi(:,:)     ! mass of hydrophillic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bctot(:,:)     ! total mass of BC (pho+phi) (col,lyr) [kg]
    real(r8), pointer :: mss_bc_col(:)      ! total mass of BC in snow column (col) [kg]
    real(r8), pointer :: mss_bc_top(:)      ! total mass of BC in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_bcphi(:,:) ! mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_bcpho(:,:) ! mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_ocpho(:,:)     ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)     ! mass of hydrophillic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_octot(:,:)     ! total mass of OC (pho+phi) (col,lyr) [kg]
    real(r8), pointer :: mss_oc_col(:)      ! total mass of OC in snow column (col) [kg]
    real(r8), pointer :: mss_oc_top(:)      ! total mass of OC in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_ocphi(:,:) ! mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_ocpho(:,:) ! mass concentration of OC species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_dst1(:,:)      ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)      ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)      ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)      ! mass of dust species 4 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dsttot(:,:)    ! total mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_col(:)     ! total mass of dust in snow column (col) [kg]
    real(r8), pointer :: mss_dst_top(:)     ! total mass of dust in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_dst1(:,:)  ! mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst2(:,:)  ! mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst3(:,:)  ! mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst4(:,:)  ! mass concentration of dust species 4 (col,lyr) [kg/kg]
    real(r8), pointer :: irrig_rate(:)         ! current irrigation rate [mm/s]
    integer,  pointer :: n_irrig_steps_left(:) ! number of time steps for which we still need to irrigate today (if 0, ignore irrig_rate)

    real(r8), pointer :: wtr_h2osoi_vol(:,:,:)    ! tracer volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: wtr_h2ocan_col(:,:)      ! tracer canopy water (mm H2O) (column-level)
    real(r8), pointer :: wtr_h2ocan_pft(:,:)      ! tracer canopy water (mm H2O) (pft-level)
    real(r8), pointer :: wtr_h2osno(:,:)          ! tracer snow water (mm H2O)
    real(r8), pointer :: wtr_qflx_irrig(:,:)      ! tracer irrigation flux (mm H2O/s)
    real(r8), pointer :: RLeafWaterSunPsn(:,:)     ! tracer ratio of sunlit leaf water, scaled by photosynthesis
    real(r8), pointer :: RLeafWaterShaPsn(:,:)     ! tracer ratio of shaded leaf water, scaled by photosynthesis
    real(r8), pointer :: RCelluloseSun(:,:)     ! tracer ratio of sunlit cellulose ()
    real(r8), pointer :: RCelluloseSha(:,:)     ! tracer ratio of sunlit cellulose ()
    real(r8), pointer :: RCelluloseSunPsn(:,:)     ! tracer ratio of sunlit cellulose times photosynthesis ()
    real(r8), pointer :: RCelluloseShaPsn(:,:)     ! tracer ratio of sunlit cellulose times photosynthesis ()
    real(r8), pointer :: RCanopyVapor_col(:,:)     ! tracer ratio of canopy vapor () (column-level)
    real(r8), pointer :: RCanopyVapor_pft(:,:)     ! tracer ratio of canopy vapor () (pft-level)
    real(r8), pointer :: RCanopyWater_col(:,:)     ! tracer ratio of canopy water () (column-level)
    real(r8), pointer :: RCanopyWater_pft(:,:)     ! tracer ratio of canopy water () (pft-level)
    real(r8), pointer :: RLeafWaterSunStdy(:,:)! tracer ratio of sunlit leaf water (steady state) () (pft-level)
    real(r8), pointer :: RLeafWaterShaStdy(:,:)! tracer ratio of shaded leaf water (steady state) () (pft-level)
    real(r8), pointer :: RLeafWaterSunDiag_col(:,:)! tracer ratio of sunlit leaf water (incl Peclet) () (column-level)
    real(r8), pointer :: RLeafWaterSunDiag_pft(:,:)! tracer ratio of sunlit leaf water (incl Peclet) () (pft-level)
    real(r8), pointer :: RLeafWaterShaDiag_col(:,:)! tracer ratio of shaded leaf water (incl Peclet) () (column-level)
    real(r8), pointer :: RLeafWaterShaDiag_pft(:,:)! tracer ratio of shaded leaf water (incl Peclet) () (pft-level)
    real(r8), pointer :: RLeafWaterSun_col(:,:)    ! tracer ratio of sunlit leaf water () (column-level)
    real(r8), pointer :: RLeafWaterSun_pft(:,:)    ! tracer ratio of sunlit leaf water () (pft-level)
    real(r8), pointer :: RLeafWaterSha_col(:,:)    ! tracer ratio of shaded leaf water () (column-level)
    real(r8), pointer :: RLeafWaterSha_pft(:,:)    ! tracer ratio of shaded leaf water () (pft-level)
    real(r8), pointer :: RXylem_pft(:,:)           ! tracer ratio of xylem water () (pft-level)
    real(r8), pointer :: RXylem_col(:,:)           ! tracer ratio of xylem water () (column-level)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer :: j,l,c,p      ! indices
    integer :: m            ! tracer index
    integer :: nlevs        ! number of levels
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8):: vwc,psi      ! for calculating soilpsi

    logical :: ldbg=.false.
!-----------------------------------------------------------------------

    if ( masterproc )then
        write(iulog,*) 'Setting initial data to non-spun up values'
    end if

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => lun%itype
    lakpoi     => lun%lakpoi
    ifspecial  => lun%ifspecial
    urbpoi     => lun%urbpoi

    ! Assign local pointers to derived subtypes components (column-level)

    ctype            => col%itype
    clandunit        => col%landunit
    snl              => cps%snl
    dz               => cps%dz
    watsat           => cps%watsat
    bsw2             => cps%bsw2
    vwcsat           => cps%vwcsat
    psisat           => cps%psisat
    soilpsi          => cps%soilpsi
    h2osoi_ice       => cws%h2osoi_ice
    h2osoi_liq       => cws%h2osoi_liq
    h2osoi_vol       => cws%h2osoi_vol
    h2ocan_col       => pws_a%h2ocan
    qflx_irrig       => cwf%qflx_irrig
    qflx_snow_melt   => cwf%qflx_snow_melt
    snowdp           => cps%snowdp
    h2osno           => cws%h2osno
    t_soisno         => ces%t_soisno
    t_lake           => ces%t_lake
    t_grnd           => ces%t_grnd
    zi               => cps%zi
    wa               => cws%wa
    wt               => cws%wt
    zwt              => cws%zwt
    snw_rds          => cps%snw_rds
    snw_rds_top      => cps%snw_rds_top
    sno_liq_top      => cps%sno_liq_top
    mss_bcpho        => cps%mss_bcpho
    mss_bcphi        => cps%mss_bcphi
    mss_bctot        => cps%mss_bctot
    mss_bc_col       => cps%mss_bc_col
    mss_bc_top       => cps%mss_bc_top
    mss_cnc_bcphi    => cps%mss_cnc_bcphi
    mss_cnc_bcpho    => cps%mss_cnc_bcpho
    mss_ocpho        => cps%mss_ocpho
    mss_ocphi        => cps%mss_ocphi
    mss_octot        => cps%mss_octot
    mss_oc_col       => cps%mss_oc_col
    mss_oc_top       => cps%mss_oc_top
    mss_cnc_ocphi    => cps%mss_cnc_ocphi
    mss_cnc_ocpho    => cps%mss_cnc_ocpho
    mss_dst1         => cps%mss_dst1
    mss_dst2         => cps%mss_dst2
    mss_dst3         => cps%mss_dst3
    mss_dst4         => cps%mss_dst4
    mss_dsttot       => cps%mss_dsttot
    mss_dst_col      => cps%mss_dst_col
    mss_dst_top      => cps%mss_dst_top
    mss_cnc_dst1     => cps%mss_cnc_dst1
    mss_cnc_dst2     => cps%mss_cnc_dst2
    mss_cnc_dst3     => cps%mss_cnc_dst3
    mss_cnc_dst4     => cps%mss_cnc_dst4
    n_irrig_steps_left => cps%n_irrig_steps_left
    irrig_rate       => cps%irrig_rate

    wtr_h2osoi_ice       => cws%wtr_h2osoi_ice
    wtr_h2osoi_liq       => cws%wtr_h2osoi_liq
    wtr_h2osoi_vol       => cws%wtr_h2osoi_vol
    wtr_h2ocan_col       => pws_a%wtr_h2ocan
    wtr_qflx_irrig       => cwf%wtr_qflx_irrig
    wtr_qflx_snow_melt   => cwf%wtr_qflx_snow_melt
    wtr_h2osno           => cws%wtr_h2osno
    wtr_wa               => cws%wtr_wa
    wtr_wt               => cws%wtr_wt

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn        => pft%column
    h2ocan_pft     => pws%h2ocan
    t_veg          => pes%t_veg
    t_ref2m        => pes%t_ref2m
    t_ref2m_u      => pes%t_ref2m_u
    t_ref2m_r      => pes%t_ref2m_r
    plandunit      => pft%landunit
    eflx_lwrad_out => pef%eflx_lwrad_out  

    wtr_h2ocan_pft     => pws%wtr_h2ocan
    RCanopyVapor_col  => pws_a%RCanopyVapor
    RCanopyWater_col  => pws_a%RCanopyWater
    RLeafWaterSunDiag_col => pws_a%RLeafWaterSunDiag
    RLeafWaterShaDiag_col => pws_a%RLeafWaterShaDiag
    RLeafWaterSun_col => pws_a%RLeafWaterSun
    RLeafWaterSha_col => pws_a%RLeafWaterSha
    RCanopyVapor_pft  => pws%RCanopyVapor
    RCanopyWater_pft  => pws%RCanopyWater
    RLeafWaterSunPsn  => pws%RLeafWaterSunPsn
    RLeafWaterShaPsn  => pws%RLeafWaterShaPsn
    RLeafWaterSunStdy => pws%RLeafWaterSunStdy
    RLeafWaterShaStdy => pws%RLeafWaterShaStdy
    RLeafWaterSunDiag_pft => pws%RLeafWaterSunDiag
    RLeafWaterShaDiag_pft => pws%RLeafWaterShaDiag
    RLeafWaterSun_pft => pws%RLeafWaterSun
    RLeafWaterSha_pft => pws%RLeafWaterSha
    RCelluloseSun     => pws%RCelluloseSun
    RCelluloseSha     => pws%RCelluloseSha
    RCelluloseSunPsn  => pws%RCelluloseSunPsn
    RCelluloseShaPsn  => pws%RCelluloseShaPsn
    RXylem_pft        => pws%RXylem
    RXylem_col        => pws_a%RXylem

    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

!    if(endp>pdbg .and. begp<pdbg) ldbg=.true.

    ! NOTE: h2ocan, h2osno, and snowdp has valid values everywhere
    ! canopy water (pft level)

    do p = begp, endp
       h2ocan_pft(p) = 0._r8
       do m = 1, pwtrc
          wtr_h2ocan_pft(p,m) = Rstnd(m)*h2ocan_pft(p)
          RCanopyVapor_pft(p,m) = Rstnd(m)
          RCanopyWater_pft(p,m) = Rstnd(m)
          RLeafWaterSunPsn(p,m) = 0._r8
          RLeafWaterShaPsn(p,m) = 0._r8
          RLeafWaterSunStdy(p,m) = Rstnd(m)
          RLeafWaterShaStdy(p,m) = Rstnd(m)
          RLeafWaterSunDiag_pft(p,m) = Rstnd(m)
          RLeafWaterShaDiag_pft(p,m) = Rstnd(m)
          RLeafWaterSun_pft(p,m) = Rstnd(m)
          RLeafWaterSha_pft(p,m) = Rstnd(m)
          RXylem_pft(p,m) = Rstnd(m)
          RCelluloseSun(p,m) = Rstnd(m)
          RCelluloseSha(p,m) = Rstnd(m)
          RCelluloseSunPsn(p,m) = 0._r8
          RCelluloseShaPsn(p,m) = 0._r8
       end do
       
       ! added for canopy water mass balance under dynamic pft weights
       !pps%tlai(p) = 0._r8
       !pps%tsai(p) = 0._r8
       !pps%elai(p) = 0._r8
       !pps%esai(p) = 0._r8
       !pps%htop(p) = 0._r8
       !pps%hbot(p) = 0._r8
       !pps%frac_veg_nosno_alb(p) = 0._r8
    end do

    do c = begc,endc

       ! canopy water (column level)

       h2ocan_col(c) = 0._r8
       qflx_snow_melt(c) = 0._r8
       do m = 1, pwtrc
          wtr_h2ocan_col(c,m) = Rstnd(m)*h2ocan_col(c)
          wtr_qflx_snow_melt(c,m) = Rstnd(m)*qflx_snow_melt(c)
          RCanopyVapor_col(c,m) = Rstnd(m)
          RCanopyWater_col(c,m) = Rstnd(m)
          RLeafWaterSunDiag_col(c,m) = Rstnd(m)
          RLeafWaterShaDiag_col(c,m) = Rstnd(m)
          RLeafWaterSun_col(c,m) = Rstnd(m)
          RLeafWaterSha_col(c,m) = Rstnd(m)
          RXylem_col(c,m) = Rstnd(m)
       end do

       ! snow water

       l = clandunit(c)

       ! Note: Glacier_mec columns are initialized with half the maximum snow cover.
       ! This gives more realistic values of qflx_glcice sooner in the simulation
       !  for columns with net ablation, at the cost of delaying ice formation
       !  in columns with net accumulation.
       if (ltype(l)==istice) then
          h2osno(c) = h2osno_max
          do m = 1, pwtrc
             wtr_h2osno(c,m) = Rstnd(m)*h2osno(c)
          end do
       elseif (ltype(l)==istice_mec) then
          h2osno(c) = 0.5_r8 * h2osno_max   ! 50 cm if h2osno_max = 1 m
          do m = 1, pwtrc
             wtr_h2osno(c,m) = Rstnd(m)*h2osno(c)
          end do
       else
          h2osno(c) = 0._r8
          do m = 1, pwtrc
             wtr_h2osno(c,m) = Rstnd(m)*h2osno(c)
          end do
       endif

       ! snow depth

       snowdp(c)  = h2osno(c) / bdsno

       ! Initialize Irrigation to zero
       if (ltype(l)==istsoil) then
          n_irrig_steps_left(c) = 0
          irrig_rate(c)         = 0.0_r8
       end if

    end do

    ! Set snow layer number, depth and thickiness

    call snowdp2lev(begc, endc)

    ! Set snow/soil temperature, note:
    ! t_soisno only has valid values over non-lake
    ! t_lake   only has valid values over lake
    ! t_grnd has valid values over all land
    ! t_veg  has valid values over all land

    ! NOTE: THESE MEMORY COPIES ARE INEFFICIENT -- SINCE nlev LOOP IS NESTED FIRST!!!!
    do c = begc,endc

       t_soisno(c,-nlevsno+1:nlevgrnd) = spval
       t_lake(c,1:nlevlak) = spval

       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake
          t_soisno(c,-nlevsno+1:0) = spval
          if (snl(c) < 0) then    !snow layer temperatures
             do j = snl(c)+1, 0
                t_soisno(c,j) = 250._r8
             enddo
          endif
          if (ltype(l)==istice .or. ltype(l)==istice_mec) then
             do j = 1, nlevgrnd
                t_soisno(c,j) = 250._r8

             end do
          else if (ltype(l) == istwet) then
             do j = 1, nlevgrnd
                t_soisno(c,j) = 277._r8
             end do
          else if (ltype(l) == isturb) then
             if (use_vancouver) then
                if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
                   ! Set road top layer to initial air temperature and interpolate other
                   ! layers down to 20C in bottom layer
                   do j = 1, nlevurb
                      t_soisno(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevurb-1)) 
                   end do
                   ! Set wall and roof layers to initial air temperature
                else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
                   do j = 1, nlevurb
                      t_soisno(c,j) = 297.56
                   end do
                else
                   do j = 1, nlevurb
                      t_soisno(c,j) = 283._r8
                   end do
                end if
             else if (use_mexicocity) then
                if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
                   ! Set road top layer to initial air temperature and interpolate other
                   ! layers down to 22C in bottom layer
                   do j = 1, nlevurb
                      t_soisno(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevurb-1)) 
                   end do
                else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
                   ! Set wall and roof layers to initial air temperature
                   do j = 1, nlevurb
                      t_soisno(c,j) = 289.46
                   end do
                else
                   do j = 1, nlevurb
                      t_soisno(c,j) = 283._r8
                   end do
                end if
             else
                if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
                   do j = 1, nlevurb
                      t_soisno(c,j) = 274._r8
                   end do
                else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
                   ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
                   ! shock from large heating/air conditioning flux
                   do j = 1, nlevurb
                      t_soisno(c,j) = 292._r8
                   end do
                end if
             end if
          else
             do j = 1, nlevgrnd
                t_soisno(c,j) = 274._r8
             end do
          endif
          t_grnd(c) = t_soisno(c,snl(c)+1)
       else                     !lake
          t_lake(c,1:nlevlak) = 277._r8
          t_grnd(c) = t_lake(c,1)
       endif

    end do

    do p = begp, endp
       c = pcolumn(p)
       l = plandunit(p)

       ! Initialize Irrigation to zero
       if (ltype(l)==istsoil) then
          qflx_irrig(c)      = 0.0_r8
          do m = 1, pwtrc
             wtr_qflx_irrig(c,m) = Rstnd(m)*qflx_irrig(c)
          end do
       end if

       if (use_vancouver) then
          t_veg(p) = 297.56
          t_ref2m(p) = 297.56
          if (urbpoi(l)) then
             t_ref2m_u(p) = 297.56
          else
             t_ref2m_u(p) = spval
          end if
          if (ifspecial(l)) then
             t_ref2m_r(p) = spval
          else
             t_ref2m_r(p) = 297.56
          end if
       else if (use_mexicocity) then
          t_veg(p) = 289.46
          t_ref2m(p) = 289.46
          if (urbpoi(l)) then
             t_ref2m_u(p) = 289.46
          else
             t_ref2m_u(p) = spval
          end if
          if (ifspecial(l)) then
             t_ref2m_r(p) = spval
          else
             t_ref2m_r(p) = 289.46
          end if
       else
          t_veg(p) = 283._r8
          t_ref2m(p) = 283._r8
          if (urbpoi(l)) then
             t_ref2m_u(p) = 283._r8
          else
             t_ref2m_u(p) = spval
          end if
          if (ifspecial(l)) then
             t_ref2m_r(p) = spval
          else
             t_ref2m_r(p) = 283._r8
          end if
       end if
       eflx_lwrad_out(p) = sb * (t_grnd(c))**4
    end do

    ! Set snow/soil ice and liquid mass

    ! volumetric water is set first and liquid content and ice lens are obtained
    ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
    ! and urban pervious road (other urban columns have zero soil water)

    h2osoi_vol(begc:endc,         1:) = spval
    h2osoi_liq(begc:endc,-nlevsno+1:) = spval
    h2osoi_ice(begc:endc,-nlevsno+1:) = spval

    wa(begc:endc)  = 5000._r8
    wt(begc:endc)  = 5000._r8
    zwt(begc:endc) = 0._r8

    do m=1, pwtrc
       wtr_h2osoi_vol(begc:endc,         1:nlevsoi,m) = Rstnd(m)*h2osoi_vol(begc:endc,         1:)
       wtr_h2osoi_liq(begc:endc,-nlevsno+1:nlevsoi,m) = Rstnd(m)*h2osoi_liq(begc:endc,-nlevsno+1:)
       wtr_h2osoi_ice(begc:endc,-nlevsno+1:nlevsoi,m) = Rstnd(m)*h2osoi_ice(begc:endc,-nlevsno+1:)
       wtr_wa(begc:endc,m)  = Rstnd(m)*5000._r8
       wtr_wt(begc:endc,m)  = Rstnd(m)*5000._r8
    end do

    do c = begc,endc
       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake
          if (ltype(l) == isturb) then
             if (ctype(c) == icol_road_perv) then
                wa(c)  = 4800._r8
                wt(c)  = wa(c)
                do m = 1, pwtrc
                   wtr_wa(c,m)  = Rstnd(m)*4800._r8
                   wtr_wt(c,m)  = wtr_wa(c,m)
                end do
                zwt(c) = (25._r8 + zi(c,nlevsoi)) - wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
             else
                wa(c)  = spval
                wt(c)  = spval
                do m = 1, pwtrc
                   wtr_wa(c,m)  = Rstnd(m)*spval
                   wtr_wt(c,m)  = Rstnd(m)*spval
                end do
                zwt(c) = spval
             end if
          else
             wa(c)  = 4800._r8
             wt(c)  = wa(c)
             do m = 1, pwtrc
                wtr_wa(c,m)  = Rstnd(m)*4800._r8
                wtr_wt(c,m)  = wtr_wa(c,m)
             end do
             zwt(c) = (25._r8 + zi(c,nlevsoi)) - wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
          end if
       end if
    end do

    do c = begc,endc
       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake

          ! volumetric water
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevsoi) then
                   h2osoi_vol(c,j) = 0.0_r8
                   do m = 1, pwtrc
                     wtr_h2osoi_vol(c,j,m) = Rstnd(m)*h2osoi_vol(c,j)
                   end do
                else
                   h2osoi_vol(c,j) = 0.3_r8
                   do m = 1, pwtrc
                     wtr_h2osoi_vol(c,j,m) = Rstnd(m)*h2osoi_vol(c,j)
                   end do
                endif
             end do
          else if (ltype(l) == isturb) then 
             nlevs = nlevurb
             do j = 1, nlevs
                if (ctype(c) == icol_road_perv .and. j <= nlevsoi) then
                   h2osoi_vol(c,j) = 0.3_r8
                   do m = 1, pwtrc
                     wtr_h2osoi_vol(c,j,m) = Rstnd(m)*h2osoi_vol(c,j)
                   end do
                else
                   h2osoi_vol(c,j) = 0.0_r8
                   do m = 1, pwtrc
                     wtr_h2osoi_vol(c,j,m) = Rstnd(m)*h2osoi_vol(c,j)
                   end do
                end if
             end do
          else if (ltype(l) == istwet) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevsoi) then
                   h2osoi_vol(c,j) = 0.0_r8
                   do m = 1, pwtrc
                     wtr_h2osoi_vol(c,j,m) = Rstnd(m)*h2osoi_vol(c,j)
                   end do
                else
                   h2osoi_vol(c,j) = 1.0_r8
                   do m = 1, pwtrc
                     wtr_h2osoi_vol(c,j,m) = Rstnd(m)*h2osoi_vol(c,j)
                   end do
                endif
             end do
          else if (ltype(l) == istice .or. ltype(l) == istice_mec) then
             nlevs = nlevgrnd 
             do j = 1, nlevs
                h2osoi_vol(c,j) = 1.0_r8
                do m = 1, pwtrc
                   wtr_h2osoi_vol(c,j,m) = Rstnd(m)*h2osoi_vol(c,j)
                end do
             end do
          endif
          do j = 1, nlevs
             h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
             do m = 1, pwtrc
                wtr_h2osoi_vol(c,j,m) = Rstnd(m)*h2osoi_vol(c,j)
             end do
        
             ! soil layers
             if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                h2osoi_ice(c,j)  = dz(c,j)*denice*h2osoi_vol(c,j)
                h2osoi_liq(c,j) = 0._r8
                do m = 1, pwtrc
                  wtr_h2osoi_ice(c,j,m) = Rstnd(m)*h2osoi_ice(c,j)
                  wtr_h2osoi_liq(c,j,m) = Rstnd(m)*h2osoi_liq(c,j)
                end do
             else
                h2osoi_ice(c,j) = 0._r8
                h2osoi_liq(c,j) = dz(c,j)*denh2o*h2osoi_vol(c,j)
                do m = 1, pwtrc
                  wtr_h2osoi_ice(c,j,m) = Rstnd(m)*h2osoi_ice(c,j)
                  wtr_h2osoi_liq(c,j,m) = Rstnd(m)*h2osoi_liq(c,j)
                end do
             endif
          end do

          if (use_cn) then
             ! soil water potential (added 10/21/03, PET)
             ! required for CN code
             if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
                nlevs = nlevgrnd
                do j = 1, nlevs
                   if (h2osoi_liq(c,j) > 0._r8) then
                      vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
                      psi = psisat(c,j) * (vwc/vwcsat(c,j))**bsw2(c,j)
                      soilpsi(c,j) = max(psi, -15.0_r8)
                      soilpsi(c,j) = min(soilpsi(c,j),0.0_r8)
                   end if
                end do
             end if
          end if
       end if
    end do

    ! Set snow

    do j = -nlevsno+1, 0
       do c = begc,endc
          l = clandunit(c)
          if (.not. lakpoi(l)) then  !not lake
             if (j > snl(c)) then
                h2osoi_ice(c,j) = dz(c,j)*250._r8
                h2osoi_liq(c,j) = 0._r8
                do m = 1, pwtrc
                  wtr_h2osoi_ice(c,j,m) = Rstnd(m)*h2osoi_ice(c,j)
                  wtr_h2osoi_liq(c,j,m) = Rstnd(m)*h2osoi_liq(c,j)
                end do
             end if
          end if
       end do
    end do


    ! initialize SNICAR fields:
    do c = begc,endc
       mss_bctot(c,:) = 0._r8
       mss_bcpho(c,:) = 0._r8
       mss_bcphi(c,:) = 0._r8
       mss_cnc_bcphi(c,:)=0._r8
       mss_cnc_bcpho(c,:)=0._r8

       mss_octot(c,:) = 0._r8
       mss_ocpho(c,:) = 0._r8
       mss_ocphi(c,:) = 0._r8
       mss_cnc_ocphi(c,:)=0._r8
       mss_cnc_ocpho(c,:)=0._r8
       
       mss_dst1(c,:) = 0._r8
       mss_dst2(c,:) = 0._r8
       mss_dst3(c,:) = 0._r8
       mss_dst4(c,:) = 0._r8
       mss_dsttot(c,:) = 0._r8
       mss_cnc_dst1(c,:)=0._r8
       mss_cnc_dst2(c,:)=0._r8
       mss_cnc_dst3(c,:)=0._r8
       mss_cnc_dst4(c,:)=0._r8
       
       if (snl(c) < 0) then
          snw_rds(c,snl(c)+1:0)        = snw_rds_min
          snw_rds(c,-nlevsno+1:snl(c)) = 0._r8
          snw_rds_top(c)               = snw_rds_min
          sno_liq_top(c) = h2osoi_liq(c,snl(c)+1) / (h2osoi_liq(c,snl(c)+1)+h2osoi_ice(c,snl(c)+1))
       elseif (h2osno(c) > 0._r8) then
          snw_rds(c,0)             = snw_rds_min
          snw_rds(c,-nlevsno+1:-1) = 0._r8
          snw_rds_top(c)           = spval
          sno_liq_top(c)           = spval
       else
          snw_rds(c,:)   = 0._r8
          snw_rds_top(c) = spval
          sno_liq_top(c) = spval
       endif
    enddo


  end subroutine mkarbinit

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: snowdp2lev
!
! !INTERFACE:
  subroutine snowdp2lev(lbc, ubc)
!
! !DESCRIPTION:
! Create snow layers and interfaces given snow depth.
! Note that cps%zi(0) is set in routine iniTimeConst.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use clm_varpar  , only : nlevsno
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: lbc, ubc                    ! column bounds
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: clandunit(:)  ! landunit index associated with each column
  real(r8), pointer :: snowdp(:)     ! snow height (m)
  logical , pointer :: lakpoi(:)     ! true => landunit is a lake point
!
! local pointers to implicit out arguments
!
  integer , pointer :: snl(:)        ! number of snow layers
  real(r8), pointer :: z(:,:)        ! layer depth  (m) over snow only
  real(r8), pointer :: dz(:,:)       ! layer thickness depth (m) over snow only
  real(r8), pointer :: zi(:,:)       ! interface depth (m) over snow only
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: c,l,j      !indices
!-----------------------------------------------------------------------
 
  ! Assign local pointers to derived subtypes components (landunit-level)

  lakpoi => lun%lakpoi

  ! Assign local pointers to derived type members (column-level)

  clandunit => col%landunit
  snowdp    => cps%snowdp
  snl       => cps%snl
  zi        => cps%zi
  dz        => cps%dz
  z         => cps%z

  ! Initialize snow levels and interfaces (lake and non-lake points)

  do c = lbc, ubc
     dz(c,-nlevsno+1: 0) = 1.e36_r8
     z (c,-nlevsno+1: 0) = 1.e36_r8
     zi(c,-nlevsno  :-1) = 1.e36_r8
  end do

  ! Determine snow levels and interfaces for non-lake points

  do c = lbc,ubc
     l = clandunit(c)
     if (.not. lakpoi(l)) then
        if (snowdp(c) < 0.01_r8) then
           snl(c) = 0
           dz(c,-nlevsno+1:0) = 0._r8
           z (c,-nlevsno+1:0) = 0._r8
           zi(c,-nlevsno+0:0) = 0._r8
        else
           if ((snowdp(c) >= 0.01_r8) .and. (snowdp(c) <= 0.03_r8)) then
              snl(c) = -1
              dz(c,0)  = snowdp(c)
           else if ((snowdp(c) > 0.03_r8) .and. (snowdp(c) <= 0.04_r8)) then
              snl(c) = -2
              dz(c,-1) = snowdp(c)/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.04_r8) .and. (snowdp(c) <= 0.07_r8)) then
              snl(c) = -2
              dz(c,-1) = 0.02_r8
              dz(c, 0) = snowdp(c) - dz(c,-1)
           else if ((snowdp(c) > 0.07_r8) .and. (snowdp(c) <= 0.12_r8)) then
              snl(c) = -3
              dz(c,-2) = 0.02_r8
              dz(c,-1) = (snowdp(c) - 0.02_r8)/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.12_r8) .and. (snowdp(c) <= 0.18_r8)) then
              snl(c) = -3
              dz(c,-2) = 0.02_r8
              dz(c,-1) = 0.05_r8
              dz(c, 0) = snowdp(c) - dz(c,-2) - dz(c,-1)
           else if ((snowdp(c) > 0.18_r8) .and. (snowdp(c) <= 0.29_r8)) then
              snl(c) = -4
              dz(c,-3) = 0.02_r8
              dz(c,-2) = 0.05_r8
              dz(c,-1) = (snowdp(c) - dz(c,-3) - dz(c,-2))/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.29_r8) .and. (snowdp(c) <= 0.41_r8)) then
              snl(c) = -4
              dz(c,-3) = 0.02_r8
              dz(c,-2) = 0.05_r8
              dz(c,-1) = 0.11_r8
              dz(c, 0) = snowdp(c) - dz(c,-3) - dz(c,-2) - dz(c,-1)
           else if ((snowdp(c) > 0.41_r8) .and. (snowdp(c) <= 0.64_r8)) then
              snl(c) = -5
              dz(c,-4) = 0.02_r8
              dz(c,-3) = 0.05_r8
              dz(c,-2) = 0.11_r8
              dz(c,-1) = (snowdp(c) - dz(c,-4) - dz(c,-3) - dz(c,-2))/2._r8
              dz(c, 0) = dz(c,-1)
           else if (snowdp(c) > 0.64_r8) then
              snl(c) = -5
              dz(c,-4) = 0.02_r8
              dz(c,-3) = 0.05_r8
              dz(c,-2) = 0.11_r8
              dz(c,-1) = 0.23_r8
              dz(c, 0)=snowdp(c)-dz(c,-4)-dz(c,-3)-dz(c,-2)-dz(c,-1)
           endif
        end if
     end if
  end do

  ! The following loop is currently not vectorized

  do c = lbc,ubc
     l = clandunit(c)
     if (.not. lakpoi(l)) then
        do j = 0, snl(c)+1, -1
           z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
           zi(c,j-1) = zi(c,j) - dz(c,j)
        end do
     end if
  end do

  ! Determine snow levels and interfaces for lake points

  do c = lbc,ubc
     l = clandunit(c)
     if (lakpoi(l)) then
        snl(c) = 0
        dz(c,-nlevsno+1:0) = 0._r8
        z (c,-nlevsno+1:0) = 0._r8
        zi(c,-nlevsno+0:0) = 0._r8
     end if
  end do

  end subroutine snowdp2lev

!-----------------------------------------------------------------------

end module mkarbinitMod
