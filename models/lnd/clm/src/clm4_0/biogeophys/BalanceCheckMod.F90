
module BalanceCheckMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: BalanceCheckMod
!
! !DESCRIPTION:
! Water and energy balance check.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
  use clm_varctl,   only: iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BeginWaterBalance  ! Initialize water balance check
  public :: BalanceCheck       ! Water and energy balance check
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BeginWaterBalance
!
! !INTERFACE:
  subroutine BeginWaterBalance(lbc, ubc, lbp, ubp, &
             num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
             num_hydrologyc, filter_hydrologyc)
!
! !DESCRIPTION:
! Initialize column-level water balance at beginning of time step
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nlevgrnd, nlevsoi
    use subgridAveMod, only : p2c
    use clm_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, &
                              icol_road_imperv
    use HydrologyTracer, only : pwtrc, Rstnd, mdbg, pdbg, cdbg
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column-index bounds
    integer, intent(in) :: lbp, ubp                    ! pft-index bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in) :: num_lakec                   ! number of column non-lake points in column filter
    integer, intent(in) :: filter_lakec(ubc-lbc+1)     ! column filter for non-lake points
    integer , intent(in)  :: num_hydrologyc               ! number of column soil points in column filter
    integer , intent(in)  :: filter_hydrologyc(ubc-lbc+1) ! column filter for soil points
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in variables
!
    real(r8), pointer :: h2osno(:)             ! snow water (mm H2O)
    real(r8), pointer :: h2osoi_ice(:,:)       ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)       ! liquid water (kg/m2)
    real(r8), pointer :: h2ocan_pft(:)         ! canopy water (mm H2O) (pft-level) 
    real(r8), pointer :: wa(:)                 ! water in the unconfined aquifer (mm)
    integer , pointer :: ctype(:)              ! column type 
    real(r8), pointer :: zwt(:)                ! water table depth (m)
    real(r8), pointer :: zi(:,:)               ! interface level below a "z" level (m)

    real(r8), pointer :: wtr_h2osno(:,:)             ! tracer snow water (mm H2O)
    real(r8), pointer :: wtr_h2osoi_ice(:,:,:)       ! tracer ice lens (kg/m2)
    real(r8), pointer :: wtr_h2osoi_liq(:,:,:)       ! tracer liquid water (kg/m2)
    real(r8), pointer :: wtr_h2ocan_pft(:,:)         ! tracer canopy water (mm H2O) (pft-level) 
!
! local pointers to original implicit out variables
!
    real(r8), pointer :: h2ocan_col(:)         ! canopy water (mm H2O) (column level)
    real(r8), pointer :: begwb(:)              ! water mass begining of the time step

    real(r8), pointer :: wtr_h2ocan_col(:,:)         ! tracer canopy water (mm H2O) (column level)
    real(r8), pointer :: wtr_begwb(:,:)              ! tracer water mass begining of the time step
    real(r8), pointer :: wtr_wa(:,:)           !tracer water in the unconfined aquifer  !JN
!
! !OTHER LOCAL VARIABLES:
!
    integer :: c, p, f, j, fc, m            ! indices
    !DEBUG
    logical :: ldbg=.false.
!    if(lbc<=cdbg .and. ubc>=cdbg) ldbg=.true.
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    h2osno             => cws%h2osno
    h2osoi_ice         => cws%h2osoi_ice
    h2osoi_liq         => cws%h2osoi_liq
    begwb              => cwbal%begwb
    h2ocan_col         => pws_a%h2ocan
    wa                 => cws%wa
    ctype              => col%itype
    zwt                => cws%zwt
    zi                 => cps%zi
    wtr_h2osno             => cws%wtr_h2osno
    wtr_h2osoi_ice         => cws%wtr_h2osoi_ice
    wtr_h2osoi_liq         => cws%wtr_h2osoi_liq
    wtr_begwb              => cwbal%wtr_begwb
    wtr_wa             => cws%wtr_wa !JN  
    wtr_h2ocan_col         => pws_a%wtr_h2ocan

    ! Assign local pointers to derived type members (pft-level)

    h2ocan_pft         => pws%h2ocan
    wtr_h2ocan_pft         => pws%wtr_h2ocan

    ! Determine beginning water balance for time step
    ! pft-level canopy water averaged to column
    call p2c(num_nolakec, filter_nolakec, h2ocan_pft, h2ocan_col)
    call p2c(pwtrc, num_nolakec, filter_nolakec, wtr_h2ocan_pft, wtr_h2ocan_col)

    do f = 1, num_hydrologyc
       c = filter_hydrologyc(f)
       if(zwt(c) <= zi(c,nlevsoi)) then
          wa(c) = 5000._r8
          do m = 1, pwtrc
            !!wtr_wa(c,m) = Rstnd(m)*wa(c) !set tracer aquifer to match bulk water. -JN
            wtr_wa(c,m) = cws%RGNIP(c,m)*wa(c) !set tracer aquifer to match bulk water. -JN
          end do
       end if
    end do

    do f = 1, num_nolakec
       c = filter_nolakec(f)
       if (ctype(c) == icol_roof .or. ctype(c) == icol_sunwall &
          .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_road_imperv) then
         begwb(c) = h2ocan_col(c) + h2osno(c)
         do m = 1, pwtrc
            wtr_begwb(c,m) = wtr_h2ocan_col(c,m) + wtr_h2osno(c,m)
         end do
!DEBUG
if(c.eq.cdbg .and. ldbg) then
write(6,*) 'begwb: roof/sunwall/shadewall/road imperv'
write(6,*) 'begwb:',h2ocan_col(cdbg),wtr_h2ocan_col(cdbg,mdbg)
write(6,*) 'begwb:',h2osno(cdbg),wtr_h2osno(cdbg,mdbg)
endif
       else
         begwb(c) = h2ocan_col(c) + h2osno(c) + wa(c)
         do m = 1, pwtrc
            wtr_begwb(c,m) = wtr_h2ocan_col(c,m) + wtr_h2osno(c,m) + wtr_wa(c,m) !JN
         end do
!DEBUG
if(c.eq.cdbg .and. ldbg) then
write(6,*) 'begwb: NOT roof/sunwall/shadewall/road imperv'
write(6,*) 'begwb:',h2ocan_col(cdbg),wtr_h2ocan_col(cdbg,mdbg)
write(6,*) 'begwb:',h2osno(cdbg),wtr_h2osno(cdbg,mdbg)
write(6,*) 'begwb:',wa(cdbg),wtr_wa(cdbg,mdbg)
endif
       end if
    end do
    do j = 1, nlevgrnd
      do f = 1, num_nolakec
         c = filter_nolakec(f)
         begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
         do m = 1, pwtrc
            wtr_begwb(c,m) = wtr_begwb(c,m) + wtr_h2osoi_ice(c,j,m) + wtr_h2osoi_liq(c,j,m)
         end do
!DEBUG
if(c.eq.cdbg .and. ldbg) then
write(6,*) 'begwb:'
write(6,*) j,h2osoi_ice(cdbg,j),wtr_h2osoi_ice(cdbg,j,mdbg)
write(6,*) j,h2osoi_liq(cdbg,j),wtr_h2osoi_liq(cdbg,j,mdbg)
endif
      end do
    end do

    do f = 1, num_lakec
       c = filter_lakec(f)
       begwb(c) = h2osno(c)
       do m = 1, pwtrc
          wtr_begwb(c,m) = wtr_h2osno(c,m)
       end do
    end do

  end subroutine BeginWaterBalance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BalanceCheck
!
! !INTERFACE:
  subroutine BalanceCheck(lbp, ubp, lbc, ubc, lbl, ubl, lbg, ubg)
!
! !DESCRIPTION:
! This subroutine accumulates the numerical truncation errors of the water
! and energy balance calculation. It is helpful to see the performance of
! the process of integration.
!
! The error for energy balance:
!
! error = abs(Net radiation - change of internal energy - Sensible heat
!             - Latent heat)
!
! The error for water balance:
!
! error = abs(precipitation - change of water storage - evaporation - runoff)
!
! !USES:
    use clmtype
    use clm_atmlnd   , only : clm_a2l
    use subgridAveMod
    use clm_time_manager , only : get_step_size, get_nstep
    use clm_varcon   , only : isturb, icol_roof, icol_sunwall, icol_shadewall, &
                              spval, icol_road_perv, icol_road_imperv, istice_mec, &
                              istdlak, istslak, istwet, istcrop, istsoil
    use clm_varctl   , only : glc_dyntopo
    use HydrologyTracer, only : pwtrc, cdbg, mdbg, pdbg, Rstnd, get_wratio, h2otiny
!!DEBUG uses
    use decompMod, only : get_proc_clumps
    use filterMod, only : filter
    use clm_varpar,only : nlevgrnd, nlevsno, nlevsoi
!
! !ARGUMENTS:
    implicit none
    integer :: lbp, ubp ! pft-index bounds
    integer :: lbc, ubc ! column-index bounds
    integer :: lbl, ubl ! landunit-index bounds
    integer :: lbg, ubg ! grid-index bounds
!
! !CALLED FROM:
! subroutine clm_driver
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 10 November 2000: Mariana Vertenstein
! Migrated to new data structures by Mariana Vertenstein and
! Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    logical , pointer :: do_capsnow(:)         ! true => do snow capping
    real(r8), pointer :: qflx_floodc(:)     ! total runoff due to flooding
    real(r8), pointer :: qflx_snow_melt(:)  ! snow melt (net)
    real(r8), pointer :: qflx_rain_grnd_col(:) ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd_col(:) ! snow on ground after interception (mm H2O/s) [+]
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    integer , pointer :: plandunit(:)       ! pft's landunit index
    integer , pointer :: cgridcell(:)       ! column's gridcell index
    integer , pointer :: clandunit(:)       ! column's landunit index
    integer , pointer :: ltype(:)           ! landunit type 
    integer , pointer :: ctype(:)           ! column type 
    real(r8), pointer :: pwtgcell(:)        ! pft's weight relative to corresponding gridcell
    real(r8), pointer :: cwtgcell(:)        ! column's weight relative to corresponding gridcell
    real(r8), pointer :: forc_rain(:)       ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)       ! snow rate [mm/s]
    real(r8), pointer :: forc_lwrad(:)      ! downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: endwb(:)           ! water mass end of the time step
    real(r8), pointer :: begwb(:)           ! water mass begining of the time step
    real(r8), pointer :: fsa(:)             ! solar radiation absorbed (total) (W/m**2)
    real(r8), pointer :: fsr(:)             ! solar radiation reflected (W/m**2)
    real(r8), pointer :: eflx_lwrad_out(:)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: sabv(:)            ! solar radiation absorbed by vegetation (W/m**2)
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: eflx_sh_tot(:)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_totg(:)    ! total sensible heat flux at grid level (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_dynbal(:)     ! energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot(:)     ! total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: eflx_soil_grnd(:)  ! soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: qflx_evap_tot(:)   ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: qflx_irrig(:)      ! irrigation flux (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)       ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)      ! qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_drain(:)      ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_runoff(:)     ! total runoff (mm H2O /s)
    real(r8), pointer :: qflx_runoffg(:)    ! total runoff at gridcell level inc land cover change flux (mm H2O /s)
    real(r8), pointer :: qflx_liq_dynbal(:) ! liq runoff due to dynamic land cover change (mm H2O /s)
    real(r8), pointer :: qflx_snwcp_ice(:)  ! excess snowfall due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: qflx_glcice(:)     ! flux of new glacier ice (mm H2O /s) [+ if ice grows]
    real(r8), pointer :: qflx_glcice_frz(:) ! ice growth (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_iceg(:) ! excess snowfall due to snow cap inc land cover change flux (mm H20/s)
    real(r8), pointer :: qflx_ice_dynbal(:) ! ice runoff due to dynamic land cover change (mm H2O /s)
    real(r8), pointer :: forc_solad(:,:)    ! direct beam radiation (vis=forc_sols , nir=forc_soll )
    real(r8), pointer :: forc_solai(:,:)    ! diffuse radiation     (vis=forc_solsd, nir=forc_solld)
    real(r8), pointer :: eflx_traffic_pft(:)    ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_wasteheat_pft(:)  ! sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: canyon_hwr(:)      ! ratio of building height to street width
    real(r8), pointer :: eflx_heat_from_ac_pft(:) !sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: h2osno(:)             ! snow water (mm H2O)
    real(r8), pointer :: h2osno_old(:)         ! snow water (mm H2O) at previous time step
    real(r8), pointer :: qflx_dew_snow(:)      ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow(:)      ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_top_soil(:)      ! net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_dew_grnd(:)      ! ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_grnd(:)     ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_prec_grnd(:)     ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_snwcp_liq(:)     ! excess liquid water due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: qflx_sl_top_soil(:)   ! liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
    integer , pointer :: snl(:)                ! number of snow layers

    real(r8), pointer :: wtr_qflx_floodc(:,:)  ! tracer total runoff due to flooding
    real(r8), pointer :: wtr_qflx_snow_melt(:,:)  !tracer snow melt (net)
    real(r8), pointer :: forc_wtr_rain(:,:) ! tracer rain rate [mm/s]
    real(r8), pointer :: forc_wtr_snow(:,:) ! tracer snow rate [mm/s]
    real(r8), pointer :: wtr_endwb(:,:)     ! tracer water mass end of the time step
    real(r8), pointer :: wtr_begwb(:,:)     ! tracer water mass begining of the time step
    real(r8), pointer :: wtr_qflx_evap_tot(:,:)   ! tracer qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: wtr_qflx_irrig(:,:)      ! tracer irrigation flux (mm H2O /s)
    real(r8), pointer :: wtr_qflx_surf(:,:)       ! tracer surface runoff (mm H2O /s)
    real(r8), pointer :: wtr_qflx_qrgwl(:,:)      ! tracer qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: wtr_qflx_drain(:,:)      ! tracer sub-surface runoff (mm H2O /s)
    real(r8), pointer :: wtr_qflx_glcice(:,:)     ! tracer flux of new glacier ice (mm H2O /s) [+ if ice grows]
    real(r8), pointer :: wtr_qflx_glcice_frz(:,:) ! tracer ice growth (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_ice(:,:)  ! tracer excess snowfall due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: wtr_qflx_runoff(:,:)     ! tracer total runoff inc land cover change flux (mm H2O /s)
    real(r8), pointer :: wtr_qflx_runoffg(:,:)    ! tracer total runoff at gridcell level inc land cover change flux (mm H2O /s)
    real(r8), pointer :: wtr_qflx_snwcp_iceg(:,:) ! tracer excess snowfall due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: wtr_qflx_ice_dynbal(:,:) ! tracer ice runoff due to dynamic land cover change (mm H2O /s)
    real(r8), pointer :: wtr_qflx_liq_dynbal(:,:) ! tracer liq runoff due to dynamic land cover change (mm H2O /s)
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: errh2o(:)          ! water conservation error (mm H2O)
    real(r8), pointer :: errsol(:)          ! solar radiation conservation error (W/m**2)
    real(r8), pointer :: errlon(:)          ! longwave radiation conservation error (W/m**2)
    real(r8), pointer :: errseb(:)          ! surface energy conservation error (W/m**2)
    real(r8), pointer :: netrad(:)          ! net radiation (positive downward) (W/m**2)
    real(r8), pointer :: errsoi_col(:)      ! column-level soil/lake energy conservation error (W/m**2)
    real(r8), pointer :: snow_sources(:)    ! snow sources (mm H2O /s)
    real(r8), pointer :: snow_sinks(:)      ! snow sinks (mm H2O /s)
    real(r8), pointer :: errh2osno(:)       ! error in h2osno (kg m-2)

    real(r8), pointer :: wtr_errh2o(:,:)    ! tracer water conservation error (mm H2O)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer  :: p,c,l,g                     ! indices
    integer  :: m                           ! tracer index
    real(r8) :: dtime                       ! land model time step (sec)
    integer  :: nstep                       ! time step number
    logical  :: found                       ! flag in search loop
    integer  :: indexp,indexc,indexl,indexg ! index of first found in search loop
    integer  :: index                       !
    real(r8) :: forc_rain_col(lbc:ubc)      ! column level rain rate [mm/s]
    real(r8) :: forc_snow_col(lbc:ubc)      ! column level snow rate [mm/s]
    real(r8) :: forc_wtr_rain_col(lbc:ubc,pwtrc)  ! tracer column level rain rate [mm/s]
    real(r8) :: forc_wtr_snow_col(lbc:ubc,pwtrc)  ! tracer column level snow rate [mm/s]

    real(r8) :: totice(2), totliq(2)  !for debugging water balance
    real(r8) :: maxwaterr             !more tracer water balance debugging
    integer  :: indwaterr             ! column index of each max wat err
    integer :: j,fc,nclumps,nc,idbg
    real(r8), pointer :: h2osoi_liq(:,:)
    real(r8), pointer :: wtr_h2osoi_liq(:,:,:)
    real(r8), pointer :: h2osoi_ice(:,:)
    real(r8), pointer :: wtr_h2osoi_ice(:,:,:)

    !DEBUG
    logical :: ldbg=.false.
!    if(lbc<=cdbg .and. ubc>=cdbg) ldbg=.true.
    h2osoi_liq => cws%h2osoi_liq
    wtr_h2osoi_liq => cws%wtr_h2osoi_liq
    h2osoi_ice => cws%h2osoi_ice
    wtr_h2osoi_ice => cws%wtr_h2osoi_ice

!-----------------------------------------------------------------------

    ! Assign local pointers to derived type scalar members (gridcell-level)

    do_capsnow          => cps%do_capsnow
    qflx_floodc         => cwf%qflx_floodc
    qflx_snow_melt      => cwf%qflx_snow_melt
    qflx_rain_grnd_col  => pwf_a%qflx_rain_grnd
    qflx_snow_grnd_col  => pwf_a%qflx_snow_grnd
    clandunit           => col%landunit
    forc_rain           => clm_a2l%forc_rain
    forc_snow           => clm_a2l%forc_snow
    forc_lwrad          => clm_a2l%forc_lwrad
    forc_solad          => clm_a2l%forc_solad
    forc_solai          => clm_a2l%forc_solai

    wtr_qflx_floodc         => cwf%wtr_qflx_floodc
    wtr_qflx_snow_melt      => cwf%wtr_qflx_snow_melt
    forc_wtr_rain     => clm_a2l%forc_wtr_rain
    forc_wtr_snow     => clm_a2l%forc_wtr_snow

    ! Assign local pointers to derived type scalar members (landunit-level)

    ltype             => lun%itype
    canyon_hwr        => lun%canyon_hwr

    ! Assign local pointers to derived type scalar members (column-level)

    ctype             => col%itype
    cgridcell         => col%gridcell
    cwtgcell          => col%wtgcell
    endwb             => cwbal%endwb
    begwb             => cwbal%begwb
    qflx_irrig        => cwf%qflx_irrig
    qflx_surf         => cwf%qflx_surf
    qflx_qrgwl        => cwf%qflx_qrgwl
    qflx_drain        => cwf%qflx_drain
    qflx_runoff       => cwf%qflx_runoff
    qflx_snwcp_ice    => pwf_a%qflx_snwcp_ice
    qflx_evap_tot     => pwf_a%qflx_evap_tot
    qflx_glcice       => cwf%qflx_glcice
    qflx_glcice_frz   => cwf%qflx_glcice_frz
    errh2o            => cwbal%errh2o
    errsoi_col        => cebal%errsoi
    h2osno             => cws%h2osno
    h2osno_old         => cws%h2osno_old
    qflx_dew_snow      => pwf_a%qflx_dew_snow
    qflx_sub_snow      => pwf_a%qflx_sub_snow
    qflx_top_soil      => cwf%qflx_top_soil
    qflx_evap_grnd     => pwf_a%qflx_evap_grnd
    qflx_dew_grnd      => pwf_a%qflx_dew_grnd
    qflx_prec_grnd     => pwf_a%qflx_prec_grnd
    qflx_snwcp_liq     => pwf_a%qflx_snwcp_liq
    qflx_sl_top_soil   => cwf%qflx_sl_top_soil
    snow_sources       => cws%snow_sources
    snow_sinks         => cws%snow_sinks
    errh2osno          => cws%errh2osno
    snl                => cps%snl

    wtr_endwb         => cwbal%wtr_endwb
    wtr_begwb         => cwbal%wtr_begwb
    wtr_qflx_irrig     => cwf%wtr_qflx_irrig
    wtr_qflx_surf     => cwf%wtr_qflx_surf
    wtr_qflx_qrgwl    => cwf%wtr_qflx_qrgwl
    wtr_qflx_drain    => cwf%wtr_qflx_drain
    wtr_qflx_snwcp_ice=> pwf_a%wtr_qflx_snwcp_ice
    wtr_qflx_evap_tot => pwf_a%wtr_qflx_evap_tot
    wtr_qflx_glcice   => cwf%wtr_qflx_glcice
    wtr_qflx_glcice_frz=> cwf%wtr_qflx_glcice_frz
    wtr_errh2o        => cwbal%wtr_errh2o
    wtr_qflx_runoff    => cwf%wtr_qflx_runoff

    ! Assign local pointers to derived type scalar members (pft-level)

    pgridcell         => pft%gridcell
    plandunit         => pft%landunit
    pwtgcell          => pft%wtgcell
    fsa               => pef%fsa
    fsr               => pef%fsr
    eflx_lwrad_out    => pef%eflx_lwrad_out
    eflx_lwrad_net    => pef%eflx_lwrad_net
    sabv              => pef%sabv
    sabg              => pef%sabg
    eflx_sh_tot       => pef%eflx_sh_tot
    eflx_lh_tot       => pef%eflx_lh_tot
    eflx_soil_grnd    => pef%eflx_soil_grnd
    errsol            => pebal%errsol
    errseb            => pebal%errseb
    errlon            => pebal%errlon
    netrad            => pef%netrad
    eflx_wasteheat_pft => pef%eflx_wasteheat_pft
    eflx_heat_from_ac_pft => pef%eflx_heat_from_ac_pft
    eflx_traffic_pft  => pef%eflx_traffic_pft

    ! Assign local pointers to derived type scalar members (gridcell-level)

    qflx_runoffg       => gwf%qflx_runoffg
    qflx_liq_dynbal    => gwf%qflx_liq_dynbal
    qflx_snwcp_iceg    => gwf%qflx_snwcp_iceg
    qflx_ice_dynbal    => gwf%qflx_ice_dynbal
    eflx_sh_totg       => gef%eflx_sh_totg
    eflx_dynbal        => gef%eflx_dynbal

    wtr_qflx_runoffg   => gwf%wtr_qflx_runoffg
    wtr_qflx_snwcp_iceg => gwf%wtr_qflx_snwcp_iceg
    wtr_qflx_liq_dynbal => gwf%wtr_qflx_liq_dynbal
    wtr_qflx_ice_dynbal => gwf%wtr_qflx_ice_dynbal

    ! Get step size and time step

    nstep = get_nstep()
    dtime = get_step_size()

    ! Determine column level incoming snow and rain
    ! Assume no incident precipitation on urban wall columns (as in Hydrology1Mod.F90).

    do c = lbc,ubc
       g = cgridcell(c)
       if (ctype(c) == icol_sunwall .or.  ctype(c) == icol_shadewall) then
          forc_rain_col(c) = 0.
          forc_snow_col(c) = 0.
       else
          forc_rain_col(c) = forc_rain(g)
          forc_snow_col(c) = forc_snow(g)
       end if
    end do

    ! Water balance check

    do c = lbc, ubc
       g = cgridcell(c)
       l = clandunit(c)

       ! Note: Some glacier_mec cols may have zero weight
       if (cwtgcell(c) > 0._r8 .or. ltype(l)==istice_mec)then
          errh2o(c) = endwb(c) - begwb(c) &
               - (forc_rain_col(c) + forc_snow_col(c) + qflx_irrig(c) + qflx_floodc(c) &
               - qflx_evap_tot(c) - qflx_surf(c) &
               - qflx_qrgwl(c) - qflx_drain(c) - qflx_snwcp_ice(c)) * dtime

          ! Suppose glc_dyntopo = T:   
          ! (1) We have qflx_snwcp_ice = 0, and excess snow has been incorporated in qflx_glcice.  
          !     This flux must be included here to complete the water balance.
          ! (2) Meltwater from ice is allowed to run off and is included in qflx_qrgwl,
          !     but the water content of the ice column has not changed (at least for now) because
          !     an equivalent ice mass has been "borrowed" from the base of the column.  That
          !     meltwater is included in qflx_glcice.
          !
          ! Note that qflx_glcice is only valid over ice_mec landunits; elsewhere it is spval

          if (glc_dyntopo .and. ltype(l)==istice_mec) then
             errh2o(c) = errh2o(c) + qflx_glcice(c)*dtime
          end if

       else

          errh2o(c) = 0.0_r8

       end if

    end do

    found = .false.
    do c = lbc, ubc
       if (abs(errh2o(c)) > 1e-7_r8) then
          found = .true.
          indexc = c
       end if
    end do
    if ( found ) then
       write(iulog,*)'WARNING:  water balance error ',&
            ' nstep = ',nstep,' indexc= ',indexc,' errh2o= ',errh2o(indexc),' landunit type= ',ltype(clandunit(indexc))
       if ((ctype(indexc) .eq. icol_roof .or. ctype(indexc) .eq. icol_road_imperv .or. &
            ctype(indexc) .eq. icol_road_perv) .and. abs(errh2o(indexc)) > 1.e-1 .and. (nstep > 2) ) then
          write(iulog,*)'clm urban model is stopping - error is greater than 1.e-1'
          write(iulog,*)'nstep = ',nstep,' indexc= ',indexc,' errh2o= ',errh2o(indexc)
          write(iulog,*)'ctype(indexc): ',ctype(indexc)
          write(iulog,*)'forc_rain    = ',forc_rain_col(indexc)
          write(iulog,*)'forc_snow    = ',forc_snow_col(indexc)
          write(iulog,*)'endwb        = ',endwb(indexc)
          write(iulog,*)'begwb        = ',begwb(indexc)
          write(iulog,*)'qflx_evap_tot= ',qflx_evap_tot(indexc)
          write(iulog,*)'qflx_irrig   = ',qflx_irrig(indexc)
          write(iulog,*)'qflx_surf    = ',qflx_surf(indexc)
          write(iulog,*)'qflx_qrgwl   = ',qflx_qrgwl(indexc)
          write(iulog,*)'qflx_drain   = ',qflx_drain(indexc)
          write(iulog,*)'qflx_flood   = ',qflx_floodc(indexc)
          write(iulog,*)'qflx_snwcp_ice   = ',qflx_snwcp_ice(indexc)
          write(iulog,*)'clm model is stopping'
          call endrun()
       else if (abs(errh2o(indexc)) > .10_r8 .and. (nstep > 2) ) then
          write(iulog,*)'clm model is stopping - error is greater than .10'
          write(iulog,*)'nstep = ',nstep,' indexc= ',indexc,' errh2o= ',errh2o(indexc)
          write(iulog,*)'ctype(indexc): ',ctype(indexc)
          write(iulog,*)'forc_rain    = ',forc_rain_col(indexc)
          write(iulog,*)'forc_snow    = ',forc_snow_col(indexc)
          write(iulog,*)'endwb        = ',endwb(indexc)
          write(iulog,*)'begwb        = ',begwb(indexc)
          write(iulog,*)'qflx_evap_tot= ',qflx_evap_tot(indexc)
          write(iulog,*)'qflx_irrig   = ',qflx_irrig(indexc)
          write(iulog,*)'qflx_surf    = ',qflx_surf(indexc)
          write(iulog,*)'qflx_qrgwl   = ',qflx_qrgwl(indexc)
          write(iulog,*)'qflx_drain   = ',qflx_drain(indexc)
          write(iulog,*)'qflx_flood   = ',qflx_floodc(indexc)
          write(iulog,*)'qflx_snwcp_ice   = ',qflx_snwcp_ice(indexc)
          write(iulog,*)'clm model is stopping'
          call endrun()
       end if
    end if

    ! Snow balance check
    do c = lbc, ubc
       g = cgridcell(c)
       l = clandunit(c)
       ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
       ! any given time step but only if there is at least one snow layer.  h2osno 
       ! also includes snow that is part of the soil column (an initial snow layer is 
       ! only created if h2osno > 10mm).

       ! --------------------------------------------------------------------- !
       ! SPM - brought in qflx_snow_melt to get snow 
       ! balance working after the flooding modifications were in place.
       ! This new check is based on a perfrostsims branch of S. Swenson.
       ! --------------------------------------------------------------------- !

       if (snl(c) .lt. 0) then
          snow_sources(c) = qflx_prec_grnd(c) + qflx_dew_snow(c) + qflx_dew_grnd(c)
          snow_sinks(c)   = qflx_sub_snow(c) + qflx_evap_grnd(c) + qflx_snow_melt(c) &
                          + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) + qflx_sl_top_soil(c)

          if (ltype(l) == istdlak) then 
             if ( do_capsnow(c) ) then
                snow_sources(c) = qflx_snow_grnd_col(c) &
                     + qflx_dew_snow(c) + qflx_dew_grnd(c)  
                
                snow_sinks(c)   = qflx_sub_snow(c) + qflx_evap_grnd(c)  &
                     + (qflx_snwcp_ice(c) + qflx_snwcp_liq(c) - qflx_prec_grnd(c))  &
                     + qflx_snow_melt(c)  + qflx_sl_top_soil(c)
             else
                snow_sources(c) = qflx_snow_grnd_col(c) &
                     + qflx_rain_grnd_col(c) &
                     + qflx_dew_snow(c) + qflx_dew_grnd(c) 
                   
                snow_sinks(c)   = qflx_sub_snow(c) + qflx_evap_grnd(c)  &
                     + qflx_snow_melt(c)  + qflx_sl_top_soil(c)
             endif
          endif

          if (ltype(l) == istsoil .or. ltype(l) == istcrop .or. ltype(l) == istwet ) then
              if ( do_capsnow(c) ) then
                 snow_sources(c) = qflx_dew_snow(c) + qflx_dew_grnd(c)  &
                     + qflx_prec_grnd(c)

                 snow_sinks(c)   = qflx_sub_snow(c) + qflx_evap_grnd(c) &
                     + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) &
                     + qflx_snow_melt(c) + qflx_sl_top_soil(c)
              else
                 snow_sources(c) = qflx_snow_grnd_col(c)  &
                     + qflx_rain_grnd_col(c) &
                     +  qflx_dew_snow(c) + qflx_dew_grnd(c)

                 snow_sinks(c)   = qflx_sub_snow(c) + qflx_evap_grnd(c) &
                     + qflx_snow_melt(c) + qflx_sl_top_soil(c)
             endif
          endif

          if (ltype(l) == istice_mec .and. glc_dyntopo) then
             snow_sinks(c) = snow_sinks(c) + qflx_glcice_frz(c)
          end if

          errh2osno(c) = (h2osno(c) - h2osno_old(c)) - (snow_sources(c) - snow_sinks(c)) * dtime
       else
          snow_sources(c) = 0._r8
          snow_sinks(c) = 0._r8
          errh2osno(c) = 0._r8
       end if
    end do

    found = .false.
    do c = lbc, ubc 
       if (cwtgcell(c) > 0._r8 .and. abs(errh2osno(c)) > 1.0e-7_r8) then
          found = .true.
          indexc = c 
       end if
    end do

    if ( found ) then
       write(iulog,*)'WARNING:  snow balance error ',&
            ' nstep = ',nstep,' indexc= ',indexc,' errh2osno= ',errh2osno(indexc)
       if (abs(errh2osno(indexc)) > 0.1_r8 .and. (nstep > 2) ) then
          write(iulog,*)'clm model is stopping - error is greater than .10'
          write(iulog,*)'nstep = ',nstep,' indexc= ',indexc,' errh2osno= ',errh2osno(indexc)
          write(iulog,*)'ltype: ', ltype(clandunit(indexc))
          write(iulog,*)'ctype(indexc): ',ctype(indexc)
          write(iulog,*)'snl: ',snl(indexc)
          write(iulog,*)'h2osno: ',h2osno(indexc)
          write(iulog,*)'h2osno_old: ',h2osno_old(indexc)
          write(iulog,*)'snow_sources: ', snow_sources(indexc)
          write(iulog,*)'snow_sinks: ', snow_sinks(indexc)
          write(iulog,*)'qflx_prec_grnd: ',qflx_prec_grnd(indexc)*dtime
          write(iulog,*)'qflx_sub_snow: ',qflx_sub_snow(indexc)*dtime
          write(iulog,*)'qflx_evap_grnd: ',qflx_evap_grnd(indexc)*dtime
          write(iulog,*)'qflx_top_soil: ',qflx_top_soil(indexc)*dtime
          write(iulog,*)'qflx_dew_snow: ',qflx_dew_snow(indexc)*dtime
          write(iulog,*)'qflx_dew_grnd: ',qflx_dew_grnd(indexc)*dtime
          write(iulog,*)'qflx_snwcp_ice: ',qflx_snwcp_ice(indexc)*dtime
          write(iulog,*)'qflx_snow_melt: ',qflx_snow_melt(indexc)*dtime
          write(iulog,*)'qflx_snwcp_liq: ',qflx_snwcp_liq(indexc)*dtime
          write(iulog,*)'qflx_sl_top_soil: ',qflx_sl_top_soil(indexc)*dtime
          write(iulog,*)'qflx_glcice_frz: ',qflx_glcice_frz(indexc)*dtime
          write(iulog,*)'clm model is stopping'
          call endrun()
       end if
    end if

    ! Tracer water balance check

    do m = 1, pwtrc

       ! Determine column level incoming tracer snow and rain
       ! Assume no incident precipitation on urban wall columns (as in Hydrology1Mod.F90).
       do c = lbc, ubc
         g = cgridcell(c)
         if (ctype(c) == icol_sunwall .or.  ctype(c) == icol_shadewall) then
            forc_wtr_rain_col(c,m) = 0._r8
            forc_wtr_snow_col(c,m) = 0._r8
         else
            forc_wtr_rain_col(c,m) = forc_wtr_rain(g,m)
            forc_wtr_snow_col(c,m) = forc_wtr_snow(g,m)
         end if
       end do

       do c = lbc, ubc
         g = cgridcell(c)
         l = clandunit(c)

         if (cwtgcell(c) > 0._r8 .or. ltype(l)==istice_mec)then
!!          errh2o(c) = endwb(c) - begwb(c) &
!!               - (forc_rain_col(c) + forc_snow_col(c) + qflx_irrig(c) + qflx_floodc(c) &
!!               - qflx_evap_tot(c) - qflx_surf(c) &
!!               - qflx_qrgwl(c) - qflx_drain(c) - qflx_snwcp_ice(c)) * dtime

           wtr_errh2o(c,m) = wtr_endwb(c,m) - wtr_begwb(c,m) &
              - (forc_wtr_rain_col(c,m) + forc_wtr_snow_col(c,m) + wtr_qflx_irrig(c,m) + wtr_qflx_floodc(c,m) &
              - wtr_qflx_evap_tot(c,m) - wtr_qflx_surf(c,m) &
              - wtr_qflx_qrgwl(c,m) - wtr_qflx_drain(c,m) - wtr_qflx_snwcp_ice(c,m)) * dtime

           if (glc_dyntopo .and. ltype(l)==istice_mec) then
              wtr_errh2o(c,m) = wtr_errh2o(c,m) + wtr_qflx_glcice(c,m)*dtime
           end if

         else

            wtr_errh2o(c,m) = 0.0_r8

         end if

       end do

        !for water balance debugging output only
        totice(:)=0._r8
        totliq(:)=0._r8
        nclumps = get_proc_clumps()
        do nc=1,nclumps
           do j = 1, nlevgrnd
              do fc = 1, filter(nc)%num_nolakec
                 c = filter(nc)%nolakec(fc)
                 totliq(1)=totliq(1)+h2osoi_liq(c,j)
                 totice(1)=totice(1)+h2osoi_ice(c,j)
                 totliq(2)=totliq(2)+wtr_h2osoi_liq(c,j,1)
                 totice(2)=totice(2)+wtr_h2osoi_ice(c,j,1)
              end do
           end do
        end do

      found = .false.
      do c = lbc, ubc
         if (cwtgcell(c) > 0._r8 .and. abs(wtr_errh2o(c,m)) > 1.e-2_r8*Rstnd(m)) then ! 1e-2 for production; 1e-5 for debugging
            found = .true.
            indexc = c
            exit
         end if
      end do

      if ( found .and. (nstep > 2) ) then
!! Suppress this output for large runs
!!         write(iulog,*)'WARNING: tracer water balance error',&
!!                        ' nstep=',nstep,' tracer index=',m,' wtr_errh2o=',wtr_errh2o(indexc,m),&
!!                        ' landunit type= ',ltype(clandunit(indexc)),' c=',index,' p=',col%pfti(indexc),col%pftf(indexc)

        if(cwtgcell(c) > 0._r8 .and. abs(wtr_errh2o(c,m)) > 10._r8*Rstnd(m)) then !model is REALLY screwing up - kill it

          write(iulog,*) 'error - BalanceCheck: Tracer water balance'
          write(iulog,*) 'snl=            ',snl(c)
          write(iulog,*) 'water tracer index and time step= ',m,nstep !JN
          if(snl(c)<0) then
            do j = snl(c)+1,0
              write(iulog,*) j,'snoliq=',cws%h2osoi_liq(c,j),cws%wtr_h2osoi_liq(c,j,m)
              write(iulog,*) j,'snoice=',cws%h2osoi_ice(c,j),cws%wtr_h2osoi_ice(c,j,m)
            enddo
          endif
          write(iulog,*) 'totliq=         ',totliq(1),totliq(2)
          write(iulog,*) 'totice=         ',totice(1),totice(2)
          write(iulog,*) '---------------------------------------------------------'
          write(iulog,*) 'errh2o=         ',errh2o(c),wtr_errh2o(c,m)
          write(iulog,*) 'endwb=          ',endwb(c),wtr_endwb(c,m)
          write(iulog,*) 'begwb=          ',begwb(c),wtr_begwb(c,m)
          write(iulog,*) 'forc_rain=      ',forc_rain_col(c),forc_wtr_rain_col(c,m)
          write(iulog,*) 'forc_snow=      ',forc_snow_col(c),forc_wtr_snow_col(c,m)
          write(iulog,*) 'qflx_irrig=     ',qflx_irrig(c),wtr_qflx_irrig(c,m)
          write(iulog,*) 'qflx_floodc=    ',qflx_floodc(c),wtr_qflx_floodc(c,m)
          write(iulog,*) 'qflx_evap_tot=  ',qflx_evap_tot(c),wtr_qflx_evap_tot(c,m)
          write(iulog,*) 'qflx_surf=      ',qflx_surf(c),wtr_qflx_surf(c,m)
          write(iulog,*) 'qflx_qrgwl=     ',qflx_qrgwl(c),wtr_qflx_qrgwl(c,m)
          write(iulog,*) 'qflx_drain=     ',qflx_drain(c),wtr_qflx_drain(c,m)
          write(iulog,*) 'qflx_snwcp_ice= ',qflx_snwcp_ice(c),wtr_qflx_snwcp_ice(c,m)
          if (glc_dyntopo) then
             write(iulog,*) 'qflx_glcice=    ',qflx_glcice(c),wtr_qflx_glcice(c,m)
          end if
 
          c = indexc
          g = cgridcell(c)
          write(iulog,*) 'ENDWB        ', endwb(c)              , wtr_endwb(c,m) 
          write(iulog,*) 'BEGWB        ', begwb(c)              , wtr_begwb(c,m) 
          write(iulog,*) 'QFLX_EVAP_TOT', qflx_evap_tot(c)*dtime, wtr_qflx_evap_tot(c,m)*dtime
          write(iulog,*) 'QFLX_SURF    ', qflx_surf(c)    *dtime, wtr_qflx_surf(c,m)    *dtime
          write(iulog,*) 'QFLX_QRGWL   ', qflx_qrgwl(c)   *dtime, wtr_qflx_qrgwl(c,m)   *dtime
          write(iulog,*) 'QFLX_DRAIN   ', qflx_drain(c)   *dtime, wtr_qflx_drain(c,m)   *dtime
          write(iulog,*) 'FORC_RAIN    ', forc_rain_col(c)*dtime, forc_wtr_rain_col(c,m)*dtime
          write(iulog,*) 'FORC_SNOW    ', forc_snow_col(c)*dtime, forc_wtr_snow_col(c,m)*dtime
 
         !Prevent iCLM4 from killing itself. -JN
          !write(iulog,*)'clm model is stopping'
          !call endrun()
        end if
      end if
    end do

    ! Energy balance checks

    do p = lbp, ubp
       l = plandunit(p)
       ! Note: Some glacier_mec pfts may have zero weight
       if (pwtgcell(p)>0._r8 .or. ltype(l)==istice_mec) then
          g = pgridcell(p)

          ! Solar radiation energy balance
          ! Do not do this check for an urban pft since it will not balance on a per-column
          ! level because of interactions between columns and since a separate check is done
          ! in the urban radiation module
          if (ltype(l) /= isturb) then
             errsol(p) = fsa(p) + fsr(p) &
                  - (forc_solad(g,1) + forc_solad(g,2) + forc_solai(g,1) + forc_solai(g,2))
          else
             errsol(p) = spval
          end if
          
          ! Longwave radiation energy balance
          ! Do not do this check for an urban pft since it will not balance on a per-column
          ! level because of interactions between columns and since a separate check is done
          ! in the urban radiation module
          if (ltype(l) /= isturb) then
             errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(g)
          else
             errlon(p) = spval
          end if
          
          ! Surface energy balance
          ! Changed to using (eflx_lwrad_net) here instead of (forc_lwrad - eflx_lwrad_out) because
          ! there are longwave interactions between urban columns (and therefore pfts). 
          ! For surfaces other than urban, (eflx_lwrad_net) equals (forc_lwrad - eflx_lwrad_out),
          ! and a separate check is done above for these terms.
          
          if (ltype(l) /= isturb) then
             errseb(p) = sabv(p) + sabg(p) + forc_lwrad(g) - eflx_lwrad_out(p) &
                         - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p)
          else
             errseb(p) = sabv(p) + sabg(p) &
                         - eflx_lwrad_net(p) &
                         - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p) &
                         + eflx_wasteheat_pft(p) + eflx_heat_from_ac_pft(p) + eflx_traffic_pft(p)
          end if
          netrad(p) = fsa(p) - eflx_lwrad_net(p)
       end if
    end do

    ! Solar radiation energy balance check

    found = .false.
    do p = lbp, ubp
       l = plandunit(p)
       if (pwtgcell(p)>0._r8 .or. ltype(l)==istice_mec) then
          if ( (errsol(p) /= spval) .and. (abs(errsol(p)) > .10_r8) ) then
             found = .true.
             indexp = p
             indexg = pgridcell(p)
          end if
       end if
    end do
    if ( found  .and. (nstep > 2) ) then
       write(iulog,100)'BalanceCheck: solar radiation balance error', nstep, indexp, errsol(indexp)
       write(iulog,*)'fsa          = ',fsa(indexp)
       write(iulog,*)'fsr          = ',fsr(indexp)
       write(iulog,*)'forc_solad(1)= ',forc_solad(indexg,1)
       write(iulog,*)'forc_solad(2)= ',forc_solad(indexg,2)
       write(iulog,*)'forc_solai(1)= ',forc_solai(indexg,1)
       write(iulog,*)'forc_solai(2)= ',forc_solai(indexg,2)
       write(iulog,*)'forc_tot     = ',forc_solad(indexg,1)+forc_solad(indexg,2)&
                                  +forc_solai(indexg,1)+forc_solai(indexg,2)
       write(iulog,*)'clm model is stopping'
       call endrun()
    end if

    ! Longwave radiation energy balance check

    found = .false.
    do p = lbp, ubp
       l = plandunit(p)
       if (pwtgcell(p)>0._r8 .or. ltype(l)==istice_mec) then
          if ( (errlon(p) /= spval) .and. (abs(errlon(p)) > .10_r8) ) then
             found = .true.
             indexp = p
          end if
       end if
    end do
    if ( found  .and. (nstep > 2) ) then
       write(iulog,100)'BalanceCheck: longwave enery balance error',nstep,indexp,errlon(indexp)
       write(iulog,*)'clm model is stopping'
       call endrun()
    end if

    ! Surface energy balance check

    found = .false.
    do p = lbp, ubp
       l = plandunit(p)
       if (pwtgcell(p)>0._r8 .or. ltype(l)==istice_mec) then
          if (abs(errseb(p)) > .10_r8 ) then
             found = .true.
             indexp = p
          end if
       end if
    end do
    if ( found  .and. (nstep > 2) ) then
       write(iulog,100)'BalanceCheck: surface flux energy balance error',nstep,indexp,errseb(indexp)
       write(iulog,*)' sabv           = ',sabv(indexp)
       write(iulog,*)' sabg           = ',sabg(indexp)
       write(iulog,*)' eflx_lwrad_net = ',eflx_lwrad_net(indexp)
       write(iulog,*)' eflx_sh_tot    = ',eflx_sh_tot(indexp)
       write(iulog,*)' eflx_lh_tot    = ',eflx_lh_tot(indexp)
       write(iulog,*)' eflx_soil_grnd = ',eflx_soil_grnd(indexp)
       write(iulog,*)'clm model is stopping'
       call endrun()
    end if

    ! Soil energy balance check

    found = .false.
    do c = lbc, ubc
       if (abs(errsoi_col(c)) > 1.0e-7_r8 ) then
          found = .true.
          indexc = c
       end if
    end do
    if ( found ) then
       write(iulog,100)'BalanceCheck: soil balance error',nstep,indexc,errsoi_col(indexc)
       if (abs(errsoi_col(indexc)) > .10_r8 .and. (nstep > 2) ) then
          write(iulog,*)'clm model is stopping'
          call endrun()
       end if
    end if

    ! Update SH and RUNOFF for dynamic land cover change energy and water fluxes
    call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                &
              qflx_runoff(lbc:ubc), qflx_runoffg(lbg:ubg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do m = 1, pwtrc
       call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                              &
                 wtr_qflx_runoff(lbc:ubc,m), wtr_qflx_runoffg(lbg:ubg,m),   &
                 c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    enddo

    do g = lbg, ubg
       qflx_runoffg(g) = qflx_runoffg(g) - qflx_liq_dynbal(g)
!! uncomment once the dynlandMod.F90 accounts for tracer water in the land use
!! changes - TW (9 Sept 2014)
       do m = 1, pwtrc
          wtr_qflx_runoffg(g,m) = wtr_qflx_runoffg(g,m) - wtr_qflx_liq_dynbal(g,m)
       end do
    enddo

    call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                      &
              qflx_snwcp_ice(lbc:ubc), qflx_snwcp_iceg(lbg:ubg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do m = 1, pwtrc
       call c2g( lbc, ubc, lbl, ubl, lbg, ubg, &
                 wtr_qflx_snwcp_ice(lbc:ubc,m), wtr_qflx_snwcp_iceg(lbg:ubg,m), &
                 c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    enddo

    do g = lbg, ubg
       qflx_snwcp_iceg(g) = qflx_snwcp_iceg(g) - qflx_ice_dynbal(g)
!! uncomment once the dynlandMod.F90 accounts for tracer water in the land use
!! changes - TW (9 Sept 2014)
       do m = 1, pwtrc
          wtr_qflx_snwcp_iceg(g,m) = wtr_qflx_snwcp_iceg(g,m) - wtr_qflx_ice_dynbal(g,m)
       end do
    enddo

    call p2g( lbp, ubp, lbc, ubc, lbl, ubl, lbg, ubg,      &
              eflx_sh_tot(lbp:ubp), eflx_sh_totg(lbg:ubg), &
              p2c_scale_type='unity',c2l_scale_type='urbanf',l2g_scale_type='unity')
    do g = lbg, ubg
       eflx_sh_totg(g) =  eflx_sh_totg(g) - eflx_dynbal(g)
    enddo

100 format (1x,a,' nstep =',i10,' point =',i6,' imbalance =',f12.6,' W/m2')
200 format (1x,a,' nstep =',i10,' point =',i6,' imbalance =',f12.6,' mm')

  end subroutine BalanceCheck

end module BalanceCheckMod
