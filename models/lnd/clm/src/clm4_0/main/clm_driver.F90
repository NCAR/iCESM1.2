module clm_driver

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_driver
!
! !DESCRIPTION:
! This module provides the main CLM driver physics calling sequence.  Most
! computations occurs over ``clumps'' of gridcells (and associated subgrid
! scale entities) assigned to each MPI process. Computation is further
! parallelized by looping over clumps on each process using shared memory OpenMP.
!
! The main CLM driver physics calling sequence for clm_driver1 is as follows:
! \begin{verbatim}
!
!     + interpMonthlyVeg      interpolate monthly vegetation data        [! CN or ! CNDV]
!       + readMonthlyVegetation read vegetation data for two months      [! CN or ! CNDV]
!
! ==== Begin Loop over clumps ====
!  -> dynland_hwcontent   Get initial heat, water content
!     + pftdyn_interp                                                    [pftdyn]
!     + dynland_hwcontent   Get new heat, water content                  [pftdyn]
! ==== End Loop over clumps  ====
!
! ==== Begin Loop over clumps ====
!  -> clm_driverInit      save of variables from previous time step
!  -> Hydrology1          canopy interception and precip on ground
!     -> FracWet          fraction of wet vegetated surface and dry elai
!  -> SurfaceRadiation    surface solar radiation
!  -> UrbanRadiation      surface solar and longwave radiation for Urban landunits
!  -> Biogeophysics1      leaf temperature and surface fluxes
!  -> BareGroundFluxes    surface fluxes for bare soil or snow-covered
!                         vegetation patches
!  -> UrbanFluxes         surface fluxes for urban landunits
!     -> MoninObukIni     first-guess Monin-Obukhov length and wind speed
!     -> FrictionVelocity friction velocity and potential temperature and
!                         humidity profiles
!  -> CanopyFluxes        leaf temperature and surface fluxes for vegetated
!                         patches
!     -> QSat             saturated vapor pressure, specific humidity, &
!                         derivatives at leaf surface
!     -> MoninObukIni     first-guess Monin-Obukhov length and wind speed
!     -> FrictionVelocity friction velocity and potential temperature and
!                         humidity profiles
!     -> Stomata          stomatal resistance and photosynthesis for
!                         sunlit leaves
!     -> Stomata          stomatal resistance and photosynthesis for
!                         shaded leaves
!     -> QSat             recalculation of saturated vapor pressure,
!                         specific humidity, & derivatives at leaf surface
!   + DustEmission        Dust mobilization
!   + DustDryDep          Dust dry deposition
!  -> Biogeophysics_Lake  lake temperature and surface fluxes
!   + VOCEmission         compute VOC emission                          [VOC]
!  -> Biogeophysics2      soil/snow & ground temp and update surface fluxes
!  -> pft2col             Average from PFT level to column level
!  -> Hydrology2          surface and soil hydrology
!  -> Hydrology_Lake      lake hydrology
!  -> SnowAge_grain       update snow effective grain size for snow radiative transfer
!   + CNEcosystemDyn      Carbon Nitrogen model ecosystem dynamics:     [CN]
!                         vegetation phenology and soil carbon  
!   + EcosystemDyn        "static" ecosystem dynamics:                  [! CN ]
!                         vegetation phenology and soil carbon  
!  -> BalanceCheck        check for errors in energy and water balances
!  -> SurfaceAlbedo       albedos for next time step
!  -> UrbanAlbedo         Urban landunit albedos for next time step
!  ====  End Loop over clumps  ====
!
! Second phase of the clm main driver, for handling history and restart file output.
!
!  -> write_diagnostic    output diagnostic if appropriate
!  -> updateAccFlds       update accumulated fields
!  -> hist_update_hbuf    accumulate history fields for time interval
!  -> htapes_wrapup       write history tapes if appropriate
!  -> restFile_write      write restart file if appropriate
! \end{verbatim}
!
! Optional subroutines are denoted by an plus (+) with the associated
! CPP token or variable in brackets at the end of the line.
!
! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clmtype
  use clm_varctl          , only : wrtdia, fpftdyn, iulog, create_glacier_mec_landunit, &
                                   use_cn, use_cndv, use_exit_spinup
  use spmdMod             , only : masterproc,mpicom
  use decompMod           , only : get_proc_clumps, get_clump_bounds, get_proc_bounds
  use filterMod           , only : filter, setFilters
  use CNDVMod             , only : dv, histCNDV
  use pftdynMod           , only : pftwt_interp
  use pftdynMod           , only : pftdyn_interp, pftdyn_wbal_init, pftdyn_wbal
  use pftdynMod           , only : pftdyn_cnbal
  use dynlandMod          , only : dynland_hwcontent
  use clm_varcon          , only : zlnd, isturb
  use clm_time_manager    , only : get_step_size,get_curr_date,get_ref_date,get_nstep,get_curr_calday,is_perpetual
  use histFileMod         , only : hist_update_hbuf, hist_htapes_wrapup, htapes_check
  use restFileMod         , only : restFile_write, restFile_filename
  use accFldsMod          , only : updateAccFlds
  use clm_driverInitMod   , only : clm_driverInit
  use BalanceCheckMod     , only : BeginWaterBalance, BalanceCheck
  use SurfaceRadiationMod , only : SurfaceRadiation
  use Hydrology1Mod       , only : Hydrology1
  use Hydrology2Mod       , only : Hydrology2
  use HydrologyLakeMod    , only : HydrologyLake
  use Biogeophysics1Mod   , only : Biogeophysics1
  use BareGroundFluxesMod , only : BareGroundFluxes
  use CanopyFluxesMod     , only : CanopyFluxes
  use Biogeophysics2Mod   , only : Biogeophysics2
  use BiogeophysicsLakeMod, only : BiogeophysicsLake
  use SurfaceAlbedoMod    , only : SurfaceAlbedo
  use pft2colMod          , only : pft2col
  use CNSetValueMod       , only : CNZeroFluxes_dwt
  use CNEcosystemDynMod   , only : CNEcosystemDyn
  use CNAnnualUpdateMod   , only : CNAnnualUpdate
  use CNBalanceCheckMod   , only : BeginCBalance, BeginNBalance, &
                                   CBalanceCheck, NBalanceCheck
  use ndepStreamMod       , only : ndep_interp
  use STATICEcosysDynMod  , only : EcosystemDyn
  use DUSTMod             , only : DustDryDep, DustEmission
  use VOCEmissionMod      , only : VOCEmission
  use seq_drydep_mod      , only : n_drydep, drydep_method, DD_XLND
  use STATICEcosysDynMod  , only : interpMonthlyVeg
  use DryDepVelocity      , only : depvel_compute
  use abortutils          , only : endrun
  use UrbanMod            , only : UrbanAlbedo, UrbanRadiation, UrbanFluxes 
  use SNICARMod           , only : SnowAge_grain
  use clm_atmlnd          , only : clm_map2gcell
  use clm_glclnd          , only : create_clm_s2x
  use perf_mod
  use HydrologyTracer     , only : TracerCheckEqual, HydrologyTracerRescale, pwtrc, &
                                   cdbg, pdbg, jdbg, mdbg, pdbgi, pdbgf, Rstnd, ludbg
  use clm_atmlnd          , only : clm_a2l
  use clm_varpar          , only : nlevsoi

!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_drv                 ! clm physics,history, restart writes
  !!public :: clm_driver1              ! Phase one of the clm driver (clm physics)
  !!public :: clm_driver2              ! Phase two of the clm driver (history, restart writes updates etc.)
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: write_diagnostic       ! Write diagnostic information to log file
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: clm_drv
!
! !INTERFACE:
subroutine clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
!
! !DESCRIPTION:
!
! First phase of the clm driver calling the clm physics. An outline of
! the calling tree is given in the description of this module.
!
! !USES:

! !ARGUMENTS:
  implicit none
  logical,         intent(in) :: doalb       ! true if time for surface albedo calc
  real(r8),        intent(in) :: nextsw_cday ! calendar day for nstep+1
  real(r8),        intent(in) :: declinp1    ! declination angle for next time step
  real(r8),        intent(in) :: declin      ! declination angle for current time step
  logical,         intent(in) :: rstwr       ! true => write restart file this step
  logical,         intent(in) :: nlend       ! true => end of run on this step
  character(len=*),intent(in) :: rdate       ! restart file time stamp for name
!
! !REVISION HISTORY:
! 2002.10.01  Mariana Vertenstein latest update to new data structures
! 11/26/03, Peter Thornton: Added new call for SurfaceRadiationSunShade when
!  cpp directive SUNSHA is set, for sunlit/shaded canopy radiation.
! 4/25/05, Peter Thornton: Made the sun/shade routine the default, no longer
!  need to have SUNSHA defined.  
! Oct/05 & Jul/07 Sam Levis: Starting dates of CNDV and crop model work
! 2/29/08, Dave Lawrence: Revised snow cover fraction according to Niu and Yang, 2007
! 3/6/09, Peter Thornton: Added declin as new argument, for daylength control on Vcmax
! 2008.11.12  B. Kauffman: morph routine casa() in casa_ecosytemDyn(), so casa
!    is more similar to CN & DGVM
! 2/25/2012 M. Vertenstein: Removed CASA references 
!
!EOP
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: clandunit(:) ! landunit index associated with each column
  integer , pointer :: itypelun(:)  ! landunit type
!
! !OTHER LOCAL VARIABLES:
  integer  :: nstep                    ! time step number
  real(r8) :: dtime                    ! land model time step (sec)
  real(r8) :: t1, t2, t3               ! temporary for mass balance checks
  integer  :: nc, fc, c, fp, p, l, g   ! indices
  integer  :: m                        ! tracer index
  integer  :: nclumps                  ! number of clumps on this processor
  integer  :: begg, endg               ! clump beginning and ending gridcell indices
  integer  :: begl, endl               ! clump beginning and ending landunit indices
  integer  :: begc, endc               ! clump beginning and ending column indices
  integer  :: begp, endp               ! clump beginning and ending pft indices
  integer  :: begg_proc, endg_proc     ! proc beginning and ending gridcell indices
  integer  :: begl_proc, endl_proc     ! proc beginning and ending landunit indices
  integer  :: begc_proc, endc_proc     ! proc beginning and ending column indices
  integer  :: begp_proc, endp_proc     ! proc beginning and ending pft indices
  type(column_type), pointer :: cptr   ! pointer to column derived subtype
  integer  :: yrp1                     ! year (0, ...) for nstep+1
  integer  :: monp1                    ! month (1, ..., 12) for nstep+1
  integer  :: dayp1                    ! day of month (1, ..., 31) for nstep+1
  integer  :: secp1                    ! seconds into current date for nstep+1
  integer  :: yr                       ! year (0, ...)
  integer  :: mon                      ! month (1, ..., 12)
  integer  :: day                      ! day of month (1, ..., 31)
  integer  :: sec                      ! seconds of the day
  integer  :: ncdate                   ! current date
  integer  :: nbdate                   ! base date (reference date)
  integer  :: kyr                      ! thousand years, equals 2 at end of first year
  character(len=256) :: filer          ! restart file name
  integer :: ier                       ! error code
!!
  logical :: ldbg=.false.    ! debugging output to ccsm.log.XXXXXX-XXXXXX ?
  integer  :: idbg         ! debugging index
  integer, pointer :: snl(:)        ! number of snow layers (j=snl(c)+1 is top)
  real(r8) :: tmp1,tmp2             ! temporary debugging guys
  real(r8), pointer :: dbg_ptr(:,:)
!!
!-----------------------------------------------------------------------

  ! Assign local pointers to derived subtypes components (landunit-level)

  itypelun            => lun%itype

  ! Assign local pointers to derived subtypes components (column-level)

  clandunit           => col%landunit

  ! Set pointers into derived type

  cptr => col

!DEBUG
snl => cps%snl

  if (use_cn) then
     ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
     if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
        call t_startf('interpMonthlyVeg')
        call interpMonthlyVeg()
        call t_stopf('interpMonthlyVeg')
     endif
  else
     ! Determine weights for time interpolation of monthly vegetation data.
     ! This also determines whether it is time to read new monthly vegetation and
     ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
     ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
     ! weights obtained here are used in subroutine ecosystemdyn to obtain time
     ! interpolated values.
     if (doalb .or. ( n_drydep > 0 .and. drydep_method == DD_XLND )) then
        call t_startf('interpMonthlyVeg')
        call interpMonthlyVeg()
        call t_stopf('interpMonthlyVeg')
     end if
  end if

  ! ============================================================================
  ! Loop over clumps
  ! ============================================================================

  nclumps = get_proc_clumps()

  !$OMP PARALLEL DO PRIVATE (nc,g,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps

     ! ============================================================================
     ! Determine clump boundaries
     ! ============================================================================

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
!debugging?
!if(pdbg>=begp .and. pdbg<=endp) ldbg=.true.
!do c=begc,endc
!if(cwf%qflx_surf(c)>0._r8 .and. cwf%qflx_surf(c).eq.cwf%qflx_surf(c)) write(6,*) 'QSURF:',c,cwf%qflx_surf(c)-(cwf%wtr_qflx_surf(c,3))/Rstnd(3)
!enddo

!DEBUG
#if 0
if(pdbg>=begp .and. pdbg<=endp) then
write(6,*) 'DEBUGGING THINGS FOR PDBG=',pdbg,pft%wtgcell(pdbg)
write(6,*) pft%gridcell(pdbg),pft%landunit(pdbg),pft%column(pdbg)
write(6,*) lun%itype(pft%landunit(pdbg)),col%itype(pft%column(pdbg)),pft%itype(pdbg),lun%canyon_hwr(pft%landunit(pdbg))
endif
#endif

if(ldbg) then
write(6,*) begp,'ClmDrv00: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) begp,'ClmDrv00: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) begp,'ClmDrv00: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) begp,'ClmDrv00: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) begp,'ClmDrv00: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) begp,'ClmDrv00: RCanopyVapor(1)=',pws%RCanopyVapor(1,mdbg)/Rstnd(mdbg)
write(6,*) begp,'ClmDrv00: RLeafWaterSun(1)=',pws%RLeafWaterSun(1,mdbg)/Rstnd(mdbg)
write(6,*) begp,'ClmDrv00: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
endif

idbg=85 !use idbg here first to avoid recompiling HydrologyTracer, which takes much longer
!if(idbg>=begc .and. idbg<=endc) write(6,*) 'c=',idbg,'p=',col%pfti(idbg),col%pftf(idbg)
if(ldbg) then
write(6,*) 'clmdrv: nstep',get_nstep(),'nc=',nc,'begg=',begg,'endg=',endg,'begc=',begc,'endc=',endc,'begp=',begp,'endp=',endp
write(6,*) 'cdbg=',cdbg,'pdbg=',col%pfti(cdbg),pdbg,col%pftf(cdbg),'c(pdbg)=',pft%column(pdbg)
elseif(pdbg>=begp .and. pdbg<=endp) then
!write(6,*) 'clmdrv: nstep',get_nstep(),'cdbg=',cdbg,'pdbg=',col%pfti(cdbg),pdbg,col%pftf(cdbg)
endif
     ! ============================================================================
     ! change pft weights and compute associated heat & water fluxes
     ! ============================================================================

     ! initialize heat and water content and dynamic balance fields to zero
     do g = begg,endg
        gwf%qflx_liq_dynbal(g) = 0._r8
        gws%gc_liq2(g)         = 0._r8
        gws%gc_liq1(g)         = 0._r8
        gwf%qflx_ice_dynbal(g) = 0._r8
        gws%gc_ice2(g)         = 0._r8 
        gws%gc_ice1(g)         = 0._r8
        gef%eflx_dynbal(g)     = 0._r8
        ges%gc_heat2(g)        = 0._r8
        ges%gc_heat1(g)        = 0._r8
        do m = 1, pwtrc
           gws%wtr_gc_liq2(g,m)         = 0._r8
           gws%wtr_gc_liq1(g,m)         = 0._r8
           gws%wtr_gc_ice2(g,m)         = 0._r8 
           gws%wtr_gc_ice1(g,m)         = 0._r8
           gwf%wtr_qflx_liq_dynbal(g,m) = 0._r8
           gwf%wtr_qflx_ice_dynbal(g,m) = 0._r8
        end do
     enddo

     !--- get initial heat,water content ---
      call dynland_hwcontent( begg, endg, gws%gc_liq1(begg:endg), &
                              gws%gc_ice1(begg:endg), ges%gc_heat1(begg:endg), gws%wtr_gc_liq1(begg:endg,:), gws%wtr_gc_ice1(begg:endg,:) )
   end do
   !$OMP END PARALLEL DO

   if (.not. use_cndv) then
      if (fpftdyn /= ' ') then
         call pftdyn_interp  ! change the pft weights
      
         !$OMP PARALLEL DO PRIVATE (nc,g,begg,endg,begl,endl,begc,endc,begp,endp)
         do nc = 1,nclumps
            call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
            
            !--- get new heat,water content: (new-old)/dt = flux into lnd model ---
            call dynland_hwcontent( begg, endg, gws%gc_liq2(begg:endg), &
                 gws%gc_ice2(begg:endg), ges%gc_heat2(begg:endg) , gws%wtr_gc_liq2(begg:endg,:), gws%wtr_gc_ice2(begg:endg,:) )
            dtime = get_step_size()
            do g = begg,endg
               gwf%qflx_liq_dynbal(g) = (gws%gc_liq2 (g) - gws%gc_liq1 (g))/dtime
               gwf%qflx_ice_dynbal(g) = (gws%gc_ice2 (g) - gws%gc_ice1 (g))/dtime
               gef%eflx_dynbal    (g) = (ges%gc_heat2(g) - ges%gc_heat1(g))/dtime
               do m = 1,pwtrc
                  gwf%wtr_qflx_liq_dynbal(g,m) = (gws%wtr_gc_liq2 (g,m) - gws%wtr_gc_liq1 (g,m))/dtime
                  gwf%wtr_qflx_ice_dynbal(g,m) = (gws%wtr_gc_ice2 (g,m) - gws%wtr_gc_ice1 (g,m))/dtime
               end do
            enddo
         end do
         !$OMP END PARALLEL DO
      end if
   end if
      
   !$OMP PARALLEL DO PRIVATE (nc,g,begg,endg,begl,endl,begc,endc,begp,endp)
   do nc = 1,nclumps
     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! Initialize the mass balance checks: water, carbon, and nitrogen
     ! ============================================================================

     call t_startf('begwbal')
     call BeginWaterBalance(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec, &
          filter(nc)%num_hydrologyc, filter(nc)%hydrologyc)
     call t_stopf('begwbal')

     if (use_cn) then
        call t_startf('begcnbal')
        call BeginCBalance(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
        call BeginNBalance(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
        call t_stopf('begcnbal')
     end if

  end do
  !$OMP END PARALLEL DO

  ! ============================================================================
  ! Initialize h2ocan_loss to zero
  ! ============================================================================

  call t_startf('pftdynwts')

  !$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps
     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
     call pftdyn_wbal_init( begc, endc )

     if (use_cndv) then
        ! NOTE: Currently CNDV and fpftdyn /= ' ' are incompatible
        call CNZeroFluxes_dwt( begc, endc, begp, endp )
        call pftwt_interp( begp, endp )
        call pftdyn_wbal( begg, endg, begc, endc, begp, endp )
        call pftdyn_cnbal( begc, endc, begp, endp )
        call setFilters(nc)
     else
        ! ============================================================================
        ! Update weights and reset filters if dynamic land use
        ! This needs to be done outside the clumps loop, but after BeginWaterBalance()
        ! The call to CNZeroFluxes_dwt() is needed regardless of fpftdyn
        ! ============================================================================
        if (use_cn) then
           call CNZeroFluxes_dwt( begc, endc, begp, endp )
        end if
        if (fpftdyn /= ' ') then
           if (use_cn) then
              call pftdyn_cnbal( begc, endc, begp, endp )
           end if
           call setFilters(nc)
        end if
     end if

  end do
  !$OMP END PARALLEL DO


  if (use_cn) then
     ! ============================================================================
     ! Update dynamic N deposition field, on albedo timestep
     ! currently being done outside clumps loop, but no reason why it couldn't be
     ! re-written to go inside.
     ! ============================================================================
     ! PET: switching CN timestep
     call ndep_interp()
  end if
  call t_stopf('pftdynwts')

  !$OMP PARALLEL DO PRIVATE (nc,l,c,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps

     ! ============================================================================
     ! Determine clump boundaries
     ! ============================================================================

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! Initialize variables from previous time step and
     ! Determine canopy interception and precipitation onto ground surface.
     ! Determine the fraction of foliage covered by water and the fraction
     ! of foliage that is dry and transpiring. Initialize snow layer if the
     ! snow accumulation exceeds 10 mm.
     ! ============================================================================
     
     ! initialize intracellular CO2 (Pa) parameters each timestep for use in VOCEmission
     pps%cisun(begp:endp) = -999._r8
     pps%cisha(begp:endp) = -999._r8

     ! initialize declination for current timestep
     do c = begc,endc
        cps%decl(c) = declin
     end do


if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv01: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv01: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_dew_snow=',pwf_a%qflx_dew_snow(cdbg),pwf_a%wtr_qflx_dew_snow(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv01: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)
enddo
!if(cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg) < 0._r8) write(6,*) 'here1 clm_driver h2osoi_liq < 0'
endif

     call t_startf('drvinit')
     call clm_driverInit(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec)
     call t_stopf('drvinit')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv02: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv02: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv02: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
endif

     ! ============================================================================
     ! Hydrology1
     ! ============================================================================

     call t_startf('hydro1')
     call Hydrology1(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('hydro1')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv03: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv03: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv03: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
endif

     ! ============================================================================
     ! Surface Radiation
     ! ============================================================================

     call t_startf('surfrad')

     ! Surface Radiation for non-urban columns

     call SurfaceRadiation(begp, endp, &
                           filter(nc)%num_nourbanp, filter(nc)%nourbanp)

     ! Surface Radiation for urban columns

     call UrbanRadiation(nc, begl, endl, begc, endc, begp, endp, &
                         filter(nc)%num_nourbanl, filter(nc)%nourbanl, &
                         filter(nc)%num_urbanl, filter(nc)%urbanl, &
                         filter(nc)%num_urbanc, filter(nc)%urbanc, &
                         filter(nc)%num_urbanp, filter(nc)%urbanp)

     call t_stopf('surfrad')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv04: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: h2osno=',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv04: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: q_ref2m=',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv04: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
!if(cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg) < 0._r8) write(6,*) 'here4 clm_driver h2osoi_liq < 0'
if( (cws%h2osno(cdbg)==0._r8 .and. cws%wtr_h2osno(cdbg,mdbg)/=0._r8) .or. &
    (cws%h2osno(cdbg)/=0._r8 .and. cws%wtr_h2osno(cdbg,mdbg)==0._r8)) write(6,*) 'here4 clm_driver sno wat zero mismatch'
do idbg=1,nlevsoi
if( (cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) .or. &
    (cws%h2osoi_liq(cdbg,idbg)/=0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)==0._r8) ) write(6,*) 'ClmDrv04:',idbg,'soi wat zero mismatch'
if( (cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) .or. &
    (cws%h2osoi_ice(cdbg,idbg)/=0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)==0._r8) ) write(6,*) 'ClmDrv04:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! Determine leaf temperature and surface fluxes based on ground
     ! temperature from previous time step.
     ! ============================================================================

     call t_startf('bgp1')
     call Biogeophysics1(begg, endg, begc, endc, begp, endp, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp1')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv05: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv05: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv05: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
!if(cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg) < 0._r8) write(6,*) 'here5 clm_driver h2osoi_liq < 0'
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv05:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv05:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! Determine bare soil or snow-covered vegetation surface temperature and fluxes
     ! Calculate Ground fluxes (frac_veg_nosno is either 1 or 0)
     ! ============================================================================

     call t_startf('bgflux')

     ! BareGroundFluxes for all pfts except lakes and urban landunits

     call BareGroundFluxes(begp, endp, &
                           filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp)
     call t_stopf('bgflux')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv06: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv06: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv06: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv06:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv06:',idbg,'soi ice zero mismatch'
enddo
endif

     ! Fluxes for all Urban landunits

     call t_startf('uflux')
     call UrbanFluxes(nc, begp, endp, begl, endl, begc, endc, &
                      filter(nc)%num_nourbanl, filter(nc)%nourbanl, &
                      filter(nc)%num_urbanl, filter(nc)%urbanl, &
                      filter(nc)%num_urbanc, filter(nc)%urbanc, &
                      filter(nc)%num_urbanp, filter(nc)%urbanp)
     call t_stopf('uflux')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv07: snl=',cps%snl(cdbg)
write(6,*) 'ClmDrv07: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv07: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv07: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv07:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv07:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! Determine non snow-covered vegetation surface temperature and fluxes
     ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
     ! and leaf water change by evapotranspiration
     ! ============================================================================

     call t_startf('canflux')
     call CanopyFluxes(begg, endg, begc, endc, begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('canflux')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv08: snl=',cps%snl(cdbg)
write(6,*) 'ClmDrv08: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv08: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv08: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv08:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv08:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! Determine lake temperature and surface fluxes
     ! ============================================================================

     call t_startf('bgplake')
     call BiogeophysicsLake(begc, endc, begp, endp, &
                            filter(nc)%num_lakec, filter(nc)%lakec, &
                            filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('bgplake')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv09: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv09: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv09: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv09:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv09:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! DUST and VOC emissions
     ! ============================================================================

     call t_startf('bgc')

     ! Dust mobilization (C. Zender's modified codes)
     call DustEmission(begp, endp, begc, endc, begl, endl, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)

     ! Dust dry deposition (C. Zender's modified codes)
     call DustDryDep(begp, endp)

     ! VOC emission (A. Guenther's MEGAN (2006) model)
     call VOCEmission(begp, endp, &
                      filter(nc)%num_soilp, filter(nc)%soilp)

     call t_stopf('bgc')

     ! ============================================================================
     ! Determine soil/snow temperatures including ground temperature and
     ! update surface fluxes for new ground temperature.
     ! ============================================================================

     call t_startf('bgp2')
     call Biogeophysics2(begl, endl, begc, endc, begp, endp, &
                         filter(nc)%num_urbanl,  filter(nc)%urbanl, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp2')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv10: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv10: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv10: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv10:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv10:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! Perform averaging from PFT level to column level
     ! ============================================================================

     call t_startf('pft2col')
     call pft2col(begc, endc, filter(nc)%num_nolakec, filter(nc)%nolakec)
     call t_stopf('pft2col')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv11: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv11: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv11: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv11:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv11:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! Vertical (column) soil and surface hydrology
     ! ============================================================================

     call t_startf('hydro2')
     call Hydrology2(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
                     filter(nc)%num_urbanc, filter(nc)%urbanc, &
                     filter(nc)%num_snowc, filter(nc)%snowc, &
                     filter(nc)%num_nosnowc, filter(nc)%nosnowc)
     call t_stopf('hydro2')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv12: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv12: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv12: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv12:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv12:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! Lake hydrology
     ! ============================================================================

     call t_startf('hylake')
     call HydrologyLake(begp, endp, &
                        filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('hylake')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv13: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv13: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_dew_grnd=',pwf%qflx_dew_grnd(pdbg),pwf%wtr_qflx_dew_grnd(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_dew_grnd=',pwf_a%qflx_dew_grnd(cdbg),pwf_a%wtr_qflx_dew_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv13: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv13:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv13:',idbg,'soi ice zero mismatch'
enddo
endif

     ! ============================================================================
     ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
     ! ============================================================================

     do c = begc,endc
        l = clandunit(c)
        if (itypelun(l) == isturb) then
           ! Urban landunit use Bonan 1996 (LSM Technical Note)
           cps%frac_sno(c) = min( cps%snowdp(c)/0.05_r8, 1._r8)
        else
           ! snow cover fraction in Niu et al. 2007
           cps%frac_sno(c) = 0.0_r8
           if(cps%snowdp(c) .gt. 0.0_r8)  then
             cps%frac_sno(c) = tanh(cps%snowdp(c)/(2.5_r8*zlnd* &
               (min(800._r8,cws%h2osno(c)/cps%snowdp(c))/100._r8)**1._r8) )
           endif
        end if
     end do

     ! ============================================================================
     ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack 
     ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation of 
     ! wet-snow metamorphism in respect of liquid-water content, Ann. Glaciol.
     ! ============================================================================
     call SnowAge_grain(begc, endc, &
          filter(nc)%num_snowc, filter(nc)%snowc, &
          filter(nc)%num_nosnowc, filter(nc)%nosnowc)

     ! ============================================================================
     ! Ecosystem dynamics: Uses CN, CNDV, or static parameterizations
     ! ============================================================================
     call t_startf('ecosysdyn')

     if (use_cn) then
        ! fully prognostic canopy structure and C-N biogeochemistry
        ! - CNDV defined: prognostic biogeography; else prescribed
        ! - crop model:   crop algorithms called from within CNEcosystemDyn
        call CNEcosystemDyn(begc,endc,begp,endp,filter(nc)%num_soilc,&
             filter(nc)%soilc, filter(nc)%num_soilp, &
             filter(nc)%soilp, filter(nc)%num_pcropp, &
             filter(nc)%pcropp, doalb)
        call CNAnnualUpdate(begc,endc,begp,endp,filter(nc)%num_soilc,&
             filter(nc)%soilc, filter(nc)%num_soilp, &
             filter(nc)%soilp)
     else
        ! Prescribed biogeography,
        ! prescribed canopy structure, some prognostic carbon fluxes
        call EcosystemDyn(begp, endp, &
                          filter(nc)%num_nolakep, filter(nc)%nolakep, &
                          doalb)
     end if
     call t_stopf('ecosysdyn')

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv14: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv14: h2ocan(c) =',pws_a%h2ocan(cdbg),pws_a%wtr_h2ocan(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_evap_soi=',pwf%qflx_evap_soi(pdbg),pwf%wtr_qflx_evap_soi(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv14: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
if( (pwf%qflx_evap_soi(idbg)>0._r8 .and. pwf%wtr_qflx_evap_soi(idbg,mdbg)<0._r8) .or. &
    (pwf%qflx_evap_soi(idbg)<0._r8 .and. pwf%wtr_qflx_evap_soi(idbg,mdbg)>0._r8)) write(6,*) 'ClmDrv14:',idbg,'qsoil sign mismatch'
if( (pwf%qflx_evap_veg(idbg)>0._r8 .and. pwf%wtr_qflx_evap_veg(idbg,mdbg)<0._r8) .or. &
    (pwf%qflx_evap_veg(idbg)<0._r8 .and. pwf%wtr_qflx_evap_veg(idbg,mdbg)>0._r8)) write(6,*) 'ClmDrv14:',idbg,'qvege sign mismatch'
if( (pwf%qflx_tran_veg(idbg)>0._r8 .and. pwf%wtr_qflx_tran_veg(idbg,mdbg)<0._r8) .or. &
    (pwf%qflx_tran_veg(idbg)<0._r8 .and. pwf%wtr_qflx_tran_veg(idbg,mdbg)>0._r8)) write(6,*) 'ClmDrv14:',idbg,'qvegt sign mismatch'
if( (pwf%qflx_evap_soi(idbg)==0._r8 .and. pwf%wtr_qflx_evap_soi(idbg,mdbg)/=0._r8) .or. &
    (pwf%qflx_evap_soi(idbg)/=0._r8 .and. pwf%wtr_qflx_evap_soi(idbg,mdbg)==0._r8)) write(6,*) 'ClmDrv14:',idbg,'qsoil zero mismatch'
if( (pwf%qflx_evap_veg(idbg)==0._r8 .and. pwf%wtr_qflx_evap_veg(idbg,mdbg)/=0._r8) .or. &
    (pwf%qflx_evap_veg(idbg)/=0._r8 .and. pwf%wtr_qflx_evap_veg(idbg,mdbg)==0._r8)) write(6,*) 'ClmDrv14:',idbg,'qvege zero mismatch'
if( (pwf%qflx_tran_veg(idbg)==0._r8 .and. pwf%wtr_qflx_tran_veg(idbg,mdbg)/=0._r8) .or. &
    (pwf%qflx_tran_veg(idbg)/=0._r8 .and. pwf%wtr_qflx_tran_veg(idbg,mdbg)==0._r8)) write(6,*) 'ClmDrv14:',idbg,'qvegt zero mismatch'
enddo
do idbg=1,nlevsoi
if(cws%h2osoi_liq(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_liq(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv14:',idbg,'soi wat zero mismatch'
if(cws%h2osoi_ice(cdbg,idbg)==0._r8 .and. cws%wtr_h2osoi_ice(cdbg,idbg,mdbg)/=0._r8) write(6,*) 'ClmDrv14:',idbg,'soi ice zero mismatch'
enddo
endif

     ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
     call depvel_compute(begp,endp)

     ! ============================================================================
     ! Correct the tracer water for numerical imprecision
     ! ============================================================================

     call HydrologyTracerRescale(begc, endc, begp, endp, &
                filter(nc)%num_lakec, filter(nc)%nolakec)

if(ldbg) then
!if(endp.gt.pdbg .and. ldbg) then
write(6,*) 'ClmDrv15: snl=',cps%snl(cdbg)
write(6,*) 'ClmDrv15: snotop_liq=',cws%h2osoi_liq(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_liq(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: snotop_ice=',cws%h2osoi_ice(cdbg,snl(cdbg)+1),cws%wtr_h2osoi_ice(cdbg,snl(cdbg)+1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: h2osoi_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: h2osoi_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: h2osno    =',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: sabv,t_veg=',pef%sabv(pdbg),pes%t_veg(pdbg)
write(6,*) 'ClmDrv15: q_ref2m   =',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_evap_veg=',pwf%qflx_evap_veg(pdbg),pwf%wtr_qflx_evap_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_tran_veg=',pwf%qflx_tran_veg(pdbg),pwf%wtr_qflx_tran_veg(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_sub_snow=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_dew_snow=',pwf%qflx_dew_snow(pdbg),pwf%wtr_qflx_dew_snow(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_infl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_top_soil=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_surf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_drain=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: qflx_evap_grnd=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: RXylem(p)=',pws%RXylem(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: RCanopyVapor(p)=',pws%RCanopyVapor(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: RLeafWaterSun(p)=',pws%RLeafWaterSun(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: RLeafWaterSha(p)=',pws%RLeafWaterSha(pdbg,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: RXylem(1)=',pws%RXylem(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: RLeafWaterSha(1)=',pws%RLeafWaterSha(1,mdbg)/Rstnd(mdbg)
write(6,*) 'ClmDrv15: PFTGUY:'
do idbg=pdbgi,pdbgf
write(6,*) idbg,'wgt=',pft%wtcol(idbg)
write(6,*) idbg,pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg)
write(6,*) idbg,pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg)
enddo
write(6,*) '*** nstep=',get_nstep(),'differences ***'
write(6,*) 'w_liq=',cws%h2osoi_liq(cdbg,jdbg),cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg),cws%h2osoi_liq(cdbg,jdbg)-(cws%wtr_h2osoi_liq(cdbg,jdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'w_ice=',cws%h2osoi_ice(cdbg,jdbg),cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg),cws%h2osoi_ice(cdbg,jdbg)-(cws%wtr_h2osoi_ice(cdbg,jdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'w_sno=',cws%h2osno(cdbg),cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg),cws%h2osno(cdbg)-(cws%wtr_h2osno(cdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'w_can=',pws%h2ocan(pdbg),pws%wtr_h2ocan(pdbg,mdbg)/Rstnd(mdbg),pws%h2ocan(pdbg)-(pws%wtr_h2ocan(pdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'q_ref=',pes%q_ref2m(pdbg),pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg),pes%q_ref2m(pdbg)-(pes%wtr_q_ref2m(pdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qsbsn=',pwf%qflx_sub_snow(pdbg),pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg),pwf%qflx_sub_snow(pdbg)-(pwf%wtr_qflx_sub_snow(pdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qinfl=',cwf%qflx_infl(cdbg),cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg),cwf%qflx_infl(cdbg)-(cwf%wtr_qflx_infl(cdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qtops=',cwf%qflx_top_soil(cdbg),cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg),cwf%qflx_top_soil(cdbg)-(cwf%wtr_qflx_top_soil(cdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qsurf=',cwf%qflx_surf(cdbg),cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg),cwf%qflx_surf(cdbg)-(cwf%wtr_qflx_surf(cdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qdrai=',cwf%qflx_drain(cdbg),cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg),cwf%qflx_drain(cdbg)-(cwf%wtr_qflx_drain(cdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qevgr=',pwf_a%qflx_evap_grnd(cdbg),pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg),pwf_a%qflx_evap_grnd(cdbg)-(pwf_a%wtr_qflx_evap_grnd(cdbg,mdbg)/Rstnd(mdbg))
write(6,*) 'PFT differences:'
do idbg=pdbgi,pdbgf
write(6,*) '*** p=',idbg,', wtcol=',pft%wtcol(idbg),'***'
write(6,*) 'h2oca=',pws%h2ocan(idbg),pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg),pws%h2ocan(idbg)-(pws%wtr_h2ocan(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qsoil=',pwf%qflx_evap_soi(idbg),pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_evap_soi(idbg)-(pwf%wtr_qflx_evap_soi(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qvege=',pwf%qflx_evap_veg(idbg),pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_evap_veg(idbg)-(pwf%wtr_qflx_evap_veg(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qvegt=',pwf%qflx_tran_veg(idbg),pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_tran_veg(idbg)-(pwf%wtr_qflx_tran_veg(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qcane=',pwf%qflx_evap_can(idbg),pwf%wtr_qflx_evap_can(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_evap_can(idbg)-(pwf%wtr_qflx_evap_can(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qtote=',pwf%qflx_evap_tot(idbg),pwf%wtr_qflx_evap_tot(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_evap_tot(idbg)-(pwf%wtr_qflx_evap_tot(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qevgr=',pwf%qflx_evap_grnd(idbg),pwf%wtr_qflx_evap_grnd(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_evap_grnd(idbg)-(pwf%wtr_qflx_evap_grnd(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qdwgr=',pwf%qflx_dew_grnd(idbg),pwf%wtr_qflx_dew_grnd(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_dew_grnd(idbg)-(pwf%wtr_qflx_dew_grnd(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qsbsn=',pwf%qflx_sub_snow(idbg),pwf%wtr_qflx_sub_snow(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_sub_snow(idbg)-(pwf%wtr_qflx_sub_snow(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qdwsn=',pwf%qflx_dew_snow(idbg),pwf%wtr_qflx_dew_snow(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_dew_snow(idbg)-(pwf%wtr_qflx_dew_snow(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qintr=',pwf%qflx_prec_intr(idbg),pwf%wtr_qflx_prec_intr(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_prec_intr(idbg)-(pwf%wtr_qflx_prec_intr(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qgrnd=',pwf%qflx_prec_grnd(idbg),pwf%wtr_qflx_prec_grnd(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_prec_grnd(idbg)-(pwf%wtr_qflx_prec_grnd(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qrain=',pwf%qflx_rain_grnd(idbg),pwf%wtr_qflx_rain_grnd(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_rain_grnd(idbg)-(pwf%wtr_qflx_rain_grnd(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qsnow=',pwf%qflx_snow_grnd(idbg),pwf%wtr_qflx_snow_grnd(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_snow_grnd(idbg)-(pwf%wtr_qflx_snow_grnd(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qscpi=',pwf%qflx_snwcp_ice(idbg),pwf%wtr_qflx_snwcp_ice(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_snwcp_ice(idbg)-(pwf%wtr_qflx_snwcp_ice(idbg,mdbg)/Rstnd(mdbg))
write(6,*) 'qscpl=',pwf%qflx_snwcp_liq(idbg),pwf%wtr_qflx_snwcp_liq(idbg,mdbg)/Rstnd(mdbg),pwf%qflx_snwcp_liq(idbg)-(pwf%wtr_qflx_snwcp_liq(idbg,mdbg)/Rstnd(mdbg))
enddo
endif

     ! ============================================================================
     ! Check the energy and water balance, also carbon and nitrogen balance
     ! ============================================================================

     call t_startf('balchk')
     call BalanceCheck(begp, endp, begc, endc, begl, endl, begg, endg)
     call t_stopf('balchk')

     if (use_exit_spinup) then
        ! skip calls to C and N balance checking during EXIT_SPINUP
        ! because the system is (intentionally) not conserving mass
        ! on the first EXIT_SPINUP doalb timestep     
     else
        if (use_cn) then
           nstep = get_nstep()
           if (nstep > 2) then
              call t_startf('cnbalchk')
              call CBalanceCheck(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
              call NBalanceCheck(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
              call t_stopf('cnbalchk')
           end if
        end if
     end if

!DEBUG
#if 0
if(begl<=ludbg .and. ludbg<=endl) then
write(6,*) 'QVEGT(p) for ldbg=ludbg pfts:',lun%pfti(ludbg),lun%pftf(ludbg)
do p=lun%pfti(ludbg),lun%pftf(ludbg)
  write(6,*) p,pwf%qflx_tran_veg(p),pwf%wtr_qflx_tran_veg(p,1)
enddo
write(6,*) 'loc(pwf%qflx_tran_veg(lun%pfti(ludbg)))      =',loc(pwf%qflx_tran_veg(lun%pfti(ludbg)))
write(6,*) 'loc(pwf%wtr_qflx_tran_veg(lun%pfti(ludbg),1))=',loc(pwf%wtr_qflx_tran_veg(lun%pfti(ludbg),1))
endif
#endif
#if 0
do p=begp,endp
  if((pwf%wtr_qflx_tran_veg(p,1).ne.pwf%wtr_qflx_tran_veg(p,1)) .and. &
     (pwf%qflx_tran_veg(p).eq.pwf%qflx_tran_veg(p)) .and. &
     (pwf%qflx_tran_veg(p) > 0._r8)) write(6,*) 'QVEGT nan problem: p',p,pwf%qflx_tran_veg(p),pwf%wtr_qflx_tran_veg(p,1)
enddo
do c=begc,endc
  if((pwf_a%wtr_qflx_tran_veg(c,1).ne.pwf_a%wtr_qflx_tran_veg(c,1)) .and. &
     (pwf_a%qflx_tran_veg(c).eq.pwf_a%qflx_tran_veg(c)) .and. &
     (pwf_a%qflx_tran_veg(c) > 0._r8)) write(6,*) 'QVEGT nan problem: c',c,pwf_a%qflx_tran_veg(c),pwf_a%wtr_qflx_tran_veg(c,1)
  if(cws%h2osno(c)*1.5 < cws%wtr_h2osno(c,2)) write(6,*) 'H2OSNO enrich problem: c',c,cws%h2osno(c),cws%wtr_h2osno(c,2)
enddo
#endif

     ! ============================================================================
     ! Determine albedos for next time step
     ! ============================================================================

     if (doalb) then
        call t_startf('surfalb')

        ! Albedos for non-urban columns

        call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, &
                           filter(nc)%num_nourbanc, filter(nc)%nourbanc, &
                           filter(nc)%num_nourbanp, filter(nc)%nourbanp, &
                           nextsw_cday, declinp1)

        call t_stopf('surfalb')

        ! Albedos for urban columns

        call t_startf('urbsurfalb')

        if (filter(nc)%num_urbanl > 0) then
           call UrbanAlbedo(nc, begl, endl, begc, endc, begp, endp,   &
                            filter(nc)%num_urbanl, filter(nc)%urbanl, &
                            filter(nc)%num_urbanc, filter(nc)%urbanc, &
                            filter(nc)%num_urbanp, filter(nc)%urbanp)
        end if

        call t_stopf('urbsurfalb')

     end if

  end do
  !$OMP END PARALLEL DO

  ! ============================================================================
  ! Determine gridcell averaged properties to send to atm (l2as and l2af derived types)
  ! ============================================================================

  call t_startf('clm_map2gcell')
  call clm_map2gcell( )
  call t_stopf('clm_map2gcell')

  ! ============================================================================
  ! Determine fields to send to glc
  ! ============================================================================
  
  if (create_glacier_mec_landunit) then
     call t_startf('create_s2x')
     call create_clm_s2x(init=.false.)
     call t_stopf('create_s2x')
  end if
  

  ! ============================================================================
  ! Write global average diagnostics to standard output
  ! ============================================================================

  nstep = get_nstep()
  if (wrtdia) call mpi_barrier(mpicom,ier)
  call t_startf('wrtdiag')
  call write_diagnostic(wrtdia, nstep)
  call t_stopf('wrtdiag')

  ! ============================================================================
  ! Update accumulators
  ! ============================================================================

  call t_startf('accum')
  call updateAccFlds()
  call t_stopf('accum')

  ! ============================================================================
  ! Update history buffer
  ! ============================================================================

  call t_startf('hbuf')
  call hist_update_hbuf()
  call t_stopf('hbuf')

  ! ============================================================================
  ! Call dv (dynamic vegetation) at last time step of year
  ! NOTE: monp1, dayp1, and secp1 correspond to nstep+1
  ! ============================================================================

  if (use_cndv) then
     call t_startf('d2dgvm')
     dtime = get_step_size()
     call get_curr_date(yrp1, monp1, dayp1, secp1, offset=int(dtime))
     if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
        
        ! Get date info.  kyr is used in lpj().  At end of first year, kyr = 2.
        call get_curr_date(yr, mon, day, sec)
        ncdate = yr*10000 + mon*100 + day
        call get_ref_date(yr, mon, day, sec)
        nbdate = yr*10000 + mon*100 + day
        kyr = ncdate/10000 - nbdate/10000 + 1
        
        if (masterproc) write(iulog,*) 'End of year. CNDV called now: ncdate=', &
             ncdate,' nbdate=',nbdate,' kyr=',kyr,' nstep=', nstep
        
        nclumps = get_proc_clumps()
        
        !$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
        do nc = 1,nclumps
           call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
           call dv(begg, endg, begp, endp,                  &
                filter(nc)%num_natvegp, filter(nc)%natvegp, kyr)
        end do
        !$OMP END PARALLEL DO
     end if
     call t_stopf('d2dgvm')
  end if

  ! ============================================================================
  ! Create history and write history tapes if appropriate
  ! ============================================================================

  call t_startf('clm_drv_io')

  call t_startf('clm_drv_io_htapes')
  call hist_htapes_wrapup( rstwr, nlend )
  call t_stopf('clm_drv_io_htapes')

  ! ============================================================================
  ! Write to CNDV history buffer if appropriate
  ! ============================================================================

  if (use_cndv) then
     if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
        call t_startf('clm_drv_io_hdgvm')
        call histCNDV()
        if (masterproc) write(iulog,*) 'Annual CNDV calculations are complete'
        call t_stopf('clm_drv_io_hdgvm')
     end if
  end if

  ! ============================================================================
  ! Write restart/initial files if appropriate
  ! ============================================================================

  if (rstwr) then
     call t_startf('clm_drv_io_wrest')
     filer = restFile_filename(rdate=rdate)
     call restFile_write( filer, nlend, rdate=rdate )
     call t_stopf('clm_drv_io_wrest')
  end if

  call t_stopf('clm_drv_io')

end subroutine clm_drv

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diagnostic
!
! !INTERFACE:
subroutine write_diagnostic (wrtdia, nstep)
!
! !DESCRIPTION:
! Write diagnostic surface temperature output each timestep.  Written to
! be fast but not bit-for-bit because order of summations can change each
! timestep.
!
! !USES:
  use clm_atmlnd , only : clm_l2a
  use decompMod  , only : get_proc_bounds, get_proc_global
  use spmdMod    , only : masterproc, npes, MPI_REAL8, MPI_ANY_SOURCE, &
                          MPI_STATUS_SIZE, mpicom, MPI_SUM
  use shr_sys_mod, only : shr_sys_flush
  use abortutils , only : endrun
!
! !ARGUMENTS:
  implicit none
  logical, intent(in) :: wrtdia     !true => write diagnostic
  integer, intent(in) :: nstep      !model time step
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: p                       ! loop index
  integer :: begp, endp              ! per-proc beginning and ending pft indices
  integer :: begc, endc              ! per-proc beginning and ending column indices
  integer :: begl, endl              ! per-proc beginning and ending landunit indices
  integer :: begg, endg              ! per-proc gridcell ending gridcell indices
  integer :: numg                    ! total number of gridcells across all processors
  integer :: numl                    ! total number of landunits across all processors
  integer :: numc                    ! total number of columns across all processors
  integer :: nump                    ! total number of pfts across all processors
  integer :: ier                     ! error status
  real(r8):: psum                    ! partial sum of ts
  real(r8):: tsum                    ! sum of ts
  real(r8):: tsxyav                  ! average ts for diagnostic output
  integer :: status(MPI_STATUS_SIZE) ! mpi status
  logical,parameter :: old_sendrecv = .false.  ! Flag if should use old send/receive method rather than MPI reduce
!------------------------------------------------------------------------

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)

  if (wrtdia) then

     call t_barrierf('sync_write_diag', mpicom)
     psum = sum(clm_l2a%t_rad(begg:endg))
     if (old_sendrecv) then
        if (masterproc) then
           tsum = psum
           do p = 1, npes-1
              call mpi_recv(psum, 1, MPI_REAL8, p, 999, mpicom, status, ier)
              if (ier/=0) then
                 write(iulog,*) 'write_diagnostic: Error in mpi_recv()'
                 call endrun
              end if
              tsum = tsum + psum
           end do
        else
           call mpi_send(psum, 1, MPI_REAL8, 0, 999, mpicom, ier)
           if (ier/=0) then
              write(iulog,*) 'write_diagnostic: Error in mpi_send()'
              call endrun
           end if
        end if
     else
        call mpi_reduce(psum, tsum, 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
        if (ier/=0) then
           write(iulog,*) 'write_diagnostic: Error in mpi_reduce()'
           call endrun
        end if
     endif
     if (masterproc) then
        tsxyav = tsum / numg
        write(iulog,1000) nstep, tsxyav
        call shr_sys_flush(iulog)
     end if

  else

     if (masterproc) then
        write(iulog,*)'clm2: completed timestep ',nstep
        call shr_sys_flush(iulog)
     end if

  endif

1000 format (1x,'nstep = ',i10,'   TS = ',f21.15)

end subroutine write_diagnostic

end module clm_driver
