
module BareGroundFluxesMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: BareGroundFluxesMod
!
! !DESCRIPTION:
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BareGroundFluxes
!
! !INTERFACE:
  subroutine BareGroundFluxes(lbp, ubp, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
    use clmtype
    use clm_atmlnd         , only : clm_a2l
    use clm_varpar         , only : nlevgrnd, nlevsoi
    use clm_varcon         , only : cpair, vkc, grav, denice, denh2o, istsoil, watmin
    use clm_varcon         , only : istcrop, tfrz
    use clm_varctl         , only : use_c13
    use shr_const_mod      , only : SHR_CONST_RGAS
    use clm_time_manager   , only : get_step_size, get_nstep
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use QSatMod            , only : QSat
    use HydrologyTracer    , only : pwtrc, ixbase, get_wratio, voltiny, h2otiny, &
                                    cdbg, pdbg, ispec, TracerCheckEqual, mdbg, rsfcmin
    use HydrologyIsotope   , only : hydro_isotope, wiso_alpliq, wiso_alpice, &
                                    wiso_alpkin_mb96, wiso_alpkin_m78, difrm, isphdo, isph218o, &
                                    tracer_constant_ratio
    use decompMod          , only : get_proc_bounds

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                     ! pft bounds
    integer, intent(in) :: num_nolakep                  ! number of pft non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)    ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! This routine originally had a long list of parameters, and also a reference to
! the entire clm derived type.  For consistency, only the derived type reference
! is passed (now pointing to the current column and pft), and the other original
! parameters are initialized locally. Using t_grnd instead of tg (tg eliminated
! as redundant).
! 1/23/02, PET: Added pft reference as parameter. All outputs will be written
! to the pft data structures, and averaged to the column level outside of
! this routine.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)        ! pft's column index
    integer , pointer :: pgridcell(:)      ! pft's gridcell index
    integer , pointer :: plandunit(:)      ! pft's landunit index
    integer , pointer :: ltype(:)          ! landunit type
    integer , pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
    real(r8), pointer :: t_grnd(:)         ! ground surface temperature [K]
    real(r8), pointer :: thm(:)            ! intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(r8), pointer :: qg(:)             ! specific humidity at ground surface [kg/kg]
    real(r8), pointer :: wtr_qg(:,:)       ! tracer specific humidity at ground surface [kg/kg]
    real(r8), pointer :: thv(:)            ! virtual potential temperature (kelvin)
    real(r8), pointer :: dqgdT(:)          ! temperature derivative of "qg"
    real(r8), pointer :: htvp(:)           ! latent heat of evaporation (/sublimation) [J/kg]
    real(r8), pointer :: beta(:)           ! coefficient of conective velocity [-]
    real(r8), pointer :: zii(:)            ! convective boundary height [m]
    real(r8), pointer :: forc_u(:)         ! atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v(:)         ! atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_th(:)        ! atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_q(:)         ! atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_wtr_q(:,:)   ! tracer atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_rho(:)       ! density (kg/m**3)
    real(r8), pointer :: forc_pbot(:)      ! atmospheric pressure (Pa)
    real(r8), pointer :: forc_hgt_u_pft(:) ! observational height of wind at pft level [m]
    real(r8), pointer :: psnsun(:)         ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: psnsha(:)         ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(r8), pointer :: z0mg_col(:)       ! roughness length, momentum [m]
    real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(r8), pointer :: wtr_h2osoi_ice(:,:,:)   ! tracer ice lens (kg/m2)
    real(r8), pointer :: wtr_h2osoi_liq(:,:,:)   ! tracer liquid water (kg/m2)
    real(r8), pointer :: dz(:,:)           ! layer depth (m)
    real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: frac_sno(:)       ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: soilbeta(:)       ! soil wetness relative to field capacity
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: z0hg_col(:)       ! roughness length, sensible heat [m]
    real(r8), pointer :: z0qg_col(:)       ! roughness length, latent heat [m]
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: dlrad(:)         ! downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad(:)         ! upward longwave radiation above the canopy [W/m2]
    real(r8), pointer :: cgrnds(:)        ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl(:)        ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8), pointer :: cgrnd(:)         ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: taux(:)          ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy(:)          ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8), pointer :: eflx_sh_grnd(:)  ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot(:)   ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:) ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:) ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: wtr_qflx_evap_soi(:,:) ! tracer soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_tot(:,:) ! tracer qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: t_ref2m(:)       ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: q_ref2m(:)       ! 2 m height surface specific humidity (kg/kg)
    real(r8), pointer :: wtr_q_ref2m(:,:) ! tracer 2 m height surface specific humidity (kg/kg)
    real(r8), pointer :: t_ref2m_r(:)     ! Rural 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: rh_ref2m_r(:)    ! Rural 2 m height surface relative humidity (%)
    real(r8), pointer :: rh_ref2m(:)      ! 2 m height surface relative humidity (%)
    real(r8), pointer :: t_veg(:)         ! vegetation temperature (Kelvin)
    real(r8), pointer :: btran(:)         ! transpiration wetness factor (0 to 1)
    real(r8), pointer :: rsfc(:)          ! surface boundary layer resistance [s/m]
    real(r8), pointer :: rlitter(:)       ! litter layer resistance [s/m]
    real(r8), pointer :: rlb(:)           ! leaf boundary layer resistance [s/m]
    real(r8), pointer :: raw1(:)          ! canopy-GCM moisture resistance (history file save) [s/m]
    real(r8), pointer :: raw2(:)          ! surface-canopy moisture resistance (history file save) [s/m]
    real(r8), pointer :: rssun(:)         ! sunlit stomatal resistance (s/m)
    real(r8), pointer :: rssha(:)         ! shaded stomatal resistance (s/m)
    real(r8), pointer :: ram1(:)          ! aerodynamical resistance (s/m)
    real(r8), pointer :: fpsn(:)          ! photosynthesis (umol CO2 /m**2 /s)
    real(r8), pointer :: rootr(:,:)       ! effective fraction of roots in each soil layer
    real(r8), pointer :: rresis(:,:)      ! root resistance by layer (0-1)  (nlevgrnd)	
    real(r8), pointer :: fsoieq(:)        ! soil water tracer equilibration fraction
    real(r8), pointer :: wtcol(:)         ! weight of pft relative to column
    integer , pointer :: snl(:)           ! number of snow layers (snl(c)+1 is surface)
    real(r8), pointer :: RXylem(:,:)      ! tracer ratio of xylem water (root-weighted soil)

!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer, parameter  :: niters = 3  ! maximum number of iterations for surface temperature
    integer  :: p,c,g,f,j,l            ! indices
    integer  :: m                      ! tracer index
    integer  :: filterp(ubp-lbp+1)     ! pft filter for vegetated pfts
    integer  :: fn                     ! number of values in local pft filter
    integer  :: fp                     ! lake filter pft index
    integer  :: iter                   ! iteration index
    real(r8) :: zldis(lbp:ubp)         ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(lbp:ubp)        ! displacement height [m]
    real(r8) :: zeta                   ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc                     ! convective velocity [m/s]
    real(r8) :: dth(lbp:ubp)           ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                   ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(lbp:ubp)           ! diff of humidity between ref. height and surface
    real(r8) :: wtr_dqh(lbp:ubp,pwtrc) ! tracer diff of humidity between ref. height and surface
    real(r8) :: obu(lbp:ubp)           ! Monin-Obukhov length (m)
    real(r8) :: ur(lbp:ubp)            ! wind speed at reference height [m/s]
    real(r8) :: um(lbp:ubp)            ! wind speed including the stablity effect [m/s]
    real(r8) :: temp1(lbp:ubp)         ! relation for potential temperature profile
    real(r8) :: temp12m(lbp:ubp)       ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(lbp:ubp)         ! relation for specific humidity profile
    real(r8) :: temp22m(lbp:ubp)       ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(lbp:ubp)         ! friction velocity [m/s]
    real(r8) :: tstar                  ! temperature scaling parameter
    real(r8) :: qstar                  ! moisture scaling parameter
    real(r8) :: thvstar                ! virtual potential temperature scaling parameter
    real(r8) :: cf                     ! heat transfer coefficient from leaves [-]
    real(r8) :: ram                    ! aerodynamical resistance [s/m]
    real(r8) :: rah                    ! thermal resistance [s/m]
    real(r8) :: raw                    ! moisture resistance [s/m]
    real(r8) :: raih                   ! temporary variable [kg/m2/s]
    real(r8) :: raiw                   ! temporary variable [kg/m2/s]
    real(r8) :: raww                   ! portion of moisture resistance which does not fractionate [s/m]
    real(r8) :: wtr_raiw(lbp:ubp,pwtrc)! tracer conductance [kg/m2/s]
    real(r8) :: fm(lbp:ubp)            ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: z0mg_pft(lbp:ubp)
    real(r8) :: z0hg_pft(lbp:ubp)
    real(r8) :: z0qg_pft(lbp:ubp)
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m 
    real(r8) :: www                    ! surface soil wetness [-]
    real(r8) :: dtime                  ! land model time step (sec)
    real(r8) :: RDew(pwtrc)            ! water tracer ratio for dew
    real(r8) :: RSoilWater(pwtrc)      ! tracer ratio of vapour in soil
    real(r8) :: RAtmVapor(pwtrc)       ! Ratio of atmospheric vapor 
    real(r8) :: h2onew                 ! soil water with evap removed
    real(r8) :: h2otot                 ! soil water remaining after evap
    real(r8) :: alpeq                  ! equilibrium fractionation for evaporation
    real(r8) :: alpk                   ! kinetic fractionation for evaporation
    real(r8) :: wtr_qflx_exchange      ! isotopic exchange flux
    real(r8) :: h2otop,h2osat,h2oeql    ! input to Mathieu and Bariac (1996) kinetic fractionation
    real(r8) :: zwtmm,zimm(6)        ! h2otop: total 'superficial' water content (first 5 layers)
                                        ! h2osat: total saturated water content (first 5 layers)
                                        ! h2oeql: total water content in eq w atm
    real(r8) :: vol_eq                   ! equilibrium volumetric water content
    real(r8) :: tempi                       ! temp variable for calculating vol_eq
    real(r8) :: temp0                       ! temp variable for calculating vol_eq
    real(r8) :: voleq1                      ! temp variable for calculating vol_eq


    logical  :: lfail                  ! fail falg for TracerCheckEqual
! for get_proc_bounds:
    integer :: begp, endp              ! per-proc beginning and ending pft indices
    integer :: lbc, ubc
    integer :: begl, endl              ! per-proc beginning and ending landunit indices
    integer :: begg, endg              ! per-proc gridcell ending gridcell indices

    real(r8) :: RSoilLiq_tmp(lbp:ubp,1:nlevsoi,pwtrc) ! temporary tracer soil water ratios
    real(r8) :: sumx,summ,allr
    real(r8) :: zbnd           ! boundary layer height (m)
    real(r8) :: resrat(lbp:ubp)         ! resistance ratio
    real(r8) :: Ren(lbp:ubp)            ! Reynolds number for surface
    real(r8) :: dd                      ! diffusivity of water vapor in air
    real(r8) :: kk(lbp:ubp,pwtrc), alpknew(lbp:ubp,pwtrc)  ! different alpk to use... a la Merlivat and Jouzel (1979)
    real(r8) :: tmp1,tmp2               ! for calculating kinetic effect
    logical :: ldbg=.false.
!if(pdbg>=lbp .and. pdbg<ubp) ldbg=.true.

!------------------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    forc_u     => clm_a2l%forc_u
    forc_v     => clm_a2l%forc_v

    ! Assign local pointers to derived type members (landunit-level)

    ltype      => lun%itype

    ! Assign local pointers to derived type members (column-level)

    forc_th    => ces%forc_th
    forc_t     => ces%forc_t
    forc_pbot  => cps%forc_pbot
    forc_rho   => cps%forc_rho
    forc_q     => cws%forc_q
    forc_wtr_q => cws%forc_wtr_q
    pcolumn    => pft%column
    pgridcell  => pft%gridcell
    frac_veg_nosno => pps%frac_veg_nosno
    dlrad  => pef%dlrad
    ulrad  => pef%ulrad
    t_grnd => ces%t_grnd
    qg     => cws%qg
    wtr_qg     => cws%wtr_qg
    z0mg_col => cps%z0mg
    z0hg_col => cps%z0hg
    z0qg_col => cps%z0qg
    thv    => ces%thv
    beta   => cps%beta
    zii    => cps%zii
    ram1   => pps%ram1
    cgrnds => pef%cgrnds
    cgrndl => pef%cgrndl
    cgrnd  => pef%cgrnd
    dqgdT  => cws%dqgdT
    htvp   => cps%htvp
    watsat         => cps%watsat
    h2osoi_ice     => cws%h2osoi_ice
    dz             => cps%dz
    h2osoi_liq     => cws%h2osoi_liq
    wtr_h2osoi_ice     => cws%wtr_h2osoi_ice
    wtr_h2osoi_liq     => cws%wtr_h2osoi_liq
    frac_sno       => cps%frac_sno
    soilbeta       => cws%soilbeta
    wtcol      => pft%wtcol
    snl     => cps%snl

    ! Assign local pointers to derived type members (pft-level)

    taux => pmf%taux
    tauy => pmf%tauy
    eflx_sh_grnd => pef%eflx_sh_grnd
    eflx_sh_tot => pef%eflx_sh_tot
    qflx_evap_soi => pwf%qflx_evap_soi
    qflx_evap_tot => pwf%qflx_evap_tot
    wtr_qflx_evap_soi => pwf%wtr_qflx_evap_soi
    wtr_qflx_evap_tot => pwf%wtr_qflx_evap_tot
    t_ref2m => pes%t_ref2m
    q_ref2m => pes%q_ref2m
    wtr_q_ref2m => pes%wtr_q_ref2m
    t_ref2m_r => pes%t_ref2m_r
    rh_ref2m_r => pes%rh_ref2m_r
    plandunit => pft%landunit
    rh_ref2m => pes%rh_ref2m
    t_veg => pes%t_veg
    thm => pes%thm
    btran => pps%btran
    rsfc  => pps%rsfc
    rlitter  => pps%rlitter
    rlb  => pps%rlb
    raw1  => pps%raw1
    raw2  => pps%raw2
    rssun => pps%rssun
    rssha => pps%rssha
    rootr => pps%rootr
    rresis => pps%rresis
    psnsun => pcf%psnsun
    psnsha => pcf%psnsha
    fpsn => pcf%fpsn
    forc_hgt_u_pft => pps%forc_hgt_u_pft
    fsoieq   => pps%fsoieq
    RXylem => pws%RXylem

    call get_proc_bounds(begg, endg, begl, endl, lbc, ubc, begp, endp)

    ! Get step size

    dtime = get_step_size()

    ! Filter pfts where frac_veg_nosno is zero

    fn = 0
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       if (frac_veg_nosno(p) == 0) then
          fn = fn + 1
          filterp(fn) = p
       end if
    end do

    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step

    do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)

       ! Initialization variables

       displa(p) = 0._r8
       dlrad(p)  = 0._r8
       ulrad(p)  = 0._r8

       ur(p) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       dth(p) = thm(p)-t_grnd(c)
       dqh(p) = forc_q(c) - qg(c)
       do m = 1, pwtrc
          wtr_dqh(p,m) = forc_wtr_q(c,m) - wtr_qg(c,m)
       end do
       dthv = dth(p)*(1._r8+0.61_r8*forc_q(c))+0.61_r8*forc_th(c)*dqh(p)
       zldis(p) = forc_hgt_u_pft(p)

       ! Copy column roughness to local pft-level arrays

       z0mg_pft(p) = z0mg_col(c)
       z0hg_pft(p) = z0hg_col(c)
       z0qg_pft(p) = z0qg_col(c)

       ! Initialize Monin-Obukhov length and wind speed

       call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg_pft(p), um(p), obu(p))

    end do

    ! Perform stability iteration
    ! Determine friction velocity, and potential temperature and humidity
    ! profiles of the surface boundary layer

    do iter = 1, niters

       call FrictionVelocity(lbp, ubp, fn, filterp, &
                             displa, z0mg_pft, z0hg_pft, z0qg_pft, &
                             obu, iter, ur, um, ustar, &
                             temp1, temp2, temp12m, temp22m, fm)

       do f = 1, fn
          p = filterp(f)
          c = pcolumn(p)
          g = pgridcell(p)

          tstar = temp1(p)*dth(p)
          qstar = temp2(p)*dqh(p)
          z0hg_pft(p) = z0mg_pft(p)/exp(0.13_r8 * (ustar(p)*z0mg_pft(p)/1.5e-5_r8)**0.45_r8)
          z0qg_pft(p) = z0hg_pft(p)
          thvstar = tstar*(1._r8+0.61_r8*forc_q(c)) + 0.61_r8*forc_th(c)*qstar
          zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

          if (zeta >= 0._r8) then                   !stable
             zeta = min(2._r8,max(zeta,0.01_r8))
             um(p) = max(ur(p),0.1_r8)
          else                                      !unstable
             zeta = max(-100._r8,min(zeta,-0.01_r8))
             wc = beta(c)*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333_r8
             um(p) = sqrt(ur(p)*ur(p) + wc*wc)
          end if
          obu(p) = zldis(p)/zeta
       end do

    end do ! end stability iteration

     do j = 1, nlevgrnd
       do f = 1, fn
          p = filterp(f)
          rootr(p,j) = 0._r8
          rresis(p,j) = 0._r8
        end do
     end do

    do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)
       l = plandunit(p)

       ! Determine aerodynamic resistances

       ram     = 1._r8/(ustar(p)*ustar(p)/um(p))
       rah     = 1._r8/(temp1(p)*ustar(p))
       raw     = 1._r8/(temp2(p)*ustar(p))
       raih    = forc_rho(c)*cpair/rah

       ! Soil evaporation resistance
       www     = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
       www     = min(max(www,0.0_r8),1._r8)

       ! Surface boundary layer resistance (TW 26 Feb 2015)
       ! (separate raw into raww (no frac) and rsfc (frac), for isotopes only)
       rsfc(p)= rsfcmin
       raww= raw-rsfc(p)
       if(soilbeta(c)>rsfcmin) then
          ! Option 1: like ISOLSM
          !rsfc(p) = 150._r8*frac_sno(c) + (1._r8-frac_sno(c))*raw*(1._r8-soilbeta(c))/soilbeta(c)
          !if(rsfc(p)<raw) then
          !   raww= raw-rsfc(p)
          !   rsfc(p)= raw-raww
          !else
          !   rsfc(p)= rsfcmin
          !   raww= raw-rsfc(p)
          !end if
          ! Option 2: soilbeta-weighted
          raww = soilbeta(c)*raw
          rsfc(p) = (1._r8-soilbeta(c))*raw
       end if

       ! Save resistances for history files
       raw2(p) = raw

       ! TW - 6 March 2015
       ! Calculate a kinetic effect based on Merlivat and Jouzel, 1979 [Option 1: from Brutsaert's original model]
#if 0
       zbnd=10._r8     ! boundary layer height... don't really have this
       dd=0.244_r8     ! diffusivity of water vapor in air (cm2/s) [Merlivat 1978]
       Ren(p)=ustar(p)*z0mg_pft(p)/1.5e-5_r8
       if(Ren(p)>=1._r8) then  !use rough surface equation
          tmp1 = 7.3_r8*(Ren(p)**0.25_r8)*((1.5e-5_r8/dd)**0.5_r8)
          resrat(p) = (1._r8/vkc)*log(zbnd/z0mg_pft(p)) - (5._r8/tmp1)
          do m=1,pwtrc
             tmp1 = (1._r8/difrm(m))**0.5_r8 -1._r8
             tmp2 = tmp1+1._r8 +resrat(p)
             kk(p,m) = tmp1/tmp2
          end do
       else                    !use smooth surface equation
          tmp1 = 13.6_r8*((1.5e-5_r8/dd)**0.67_r8)
          resrat(p) = (1._r8/vkc)*log(Ren(p)*zbnd/(30._r8*z0mg_pft(p))) /tmp1
          do m=1,pwtrc
             tmp1 = (1._r8/difrm(m))**0.67_r8 -1._r8
             tmp2 = tmp1+1._r8 +resrat(p)
             kk(p,m) = tmp1/tmp2
          end do
       end if
#endif
       ! Calculate a kinetic effect based on Merlivat and Jouzel, 1979 [Option 2: using fit from MJ79 Fig 2]
       Ren(p)=ustar(p)*z0mg_pft(p)/1.5e-5_r8
       if(Ren(p)>=1._r8) then  !use rough surface equation
         do m=1,pwtrc
           if(ispec(m)==ixbase)   kk(p,m)=0._r8                       !16O
           if(ispec(m)==isph218o) kk(p,m)=8.82_r8*ustar(p) + 0.472_r8 !18O
         end do
       else                    !use smooth surface equation
         do m=1,pwtrc
           if(ispec(m)==ixbase)   kk(p,m)=0._r8   !16O
           if(ispec(m)==isph218o) kk(p,m)=6._r8   !18O
         end do
       end if
       kk(p,isphdo)=kk(p,isph218o)*0.88_r8 ! D

       if(tracer_constant_ratio) kk(p,:)=0._r8

       !changed by K.Sakaguchi. Soilbeta is used for evaporation
       if (dqh(p) .gt. 0._r8) then   !dew  (beta is not applied, just like rsoil used to be)
          raiw    = forc_rho(c)/(raw)
          do m = 1, pwtrc
             ! M78 or MB96 kinetic factor
             !if(cws%h2osno(c)>0._r8) then
                alpk = wiso_alpkin_m78(ispec(m),.67_r8,.true.)
             !else
             !   alpk = wiso_alpkin_mb96(ispec(m),cws%h2osoi_vol(c,1),watmin,watsat(c,1),.true.)
             !end if
             !wtr_raiw(p,m) = forc_rho(c)/(raww+alpk*rsfc(p))
             ! MJ79 kinetic factor
             !alpk = 1000._r8/(1000._r8-kk(p,m))
             wtr_raiw(p,m) = forc_rho(c)/(alpk*raw)
          end do
       else
       ! Lee and Pielke 1992 beta is applied
          raiw    = soilbeta(c)*forc_rho(c)/(raw)
          do m = 1, pwtrc
             ! M78 or MB96 kinetic factor
             !if(cws%h2osno(c)>0._r8) then
                alpk = wiso_alpkin_m78(ispec(m),.67_r8,.true.)
             !else
             !   alpk = wiso_alpkin_mb96(ispec(m),cws%h2osoi_vol(c,1),watmin,watsat(c,1),.true.)
             !end if
             !wtr_raiw(p,m) = soilbeta(c)*forc_rho(c)/(raww+alpk*rsfc(p))
             ! MJ79 kinetic factor
             !alpk = 1000._r8/(1000._r8-kk(p,m))
             wtr_raiw(p,m) = soilbeta(c)*forc_rho(c)/(alpk*raw)
          end do
       end if

       ram1(p) = ram  !pass value to global variable

       ! Output to pft-level data structures
       ! Derivative of fluxes with respect to ground temperature

       cgrnds(p) = raih
       cgrndl(p) = raiw*dqgdT(c)
       cgrnd(p)  = cgrnds(p) + htvp(c)*cgrndl(p)

       ! Surface fluxes of momentum, sensible and latent heat
       ! using ground temperatures from previous time step

       taux(p)          = -forc_rho(c)*forc_u(g)/ram
       tauy(p)          = -forc_rho(c)*forc_v(g)/ram
       eflx_sh_grnd(p)  = -raih*dth(p)
       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       qflx_evap_soi(p) = -raiw*dqh(p)
       qflx_evap_tot(p) = qflx_evap_soi(p)

       ! Compute the soil water flushing time scale and factor
       ! Must be evaporating liquid.
       ! (modify for kinetic effects?)

       if (qflx_evap_soi(p) > 0._r8 .and. t_grnd(c) > tfrz) then
         h2onew = max(wtcol(p)*h2osoi_liq(c,1) - dtime*qflx_evap_soi(p) ,0.0_r8)
         fsoieq(p)   = raiw*forc_q(c)/max(h2onew,h2otiny)
         fsoieq(p)   = exp(-dtime*fsoieq(p))
       else
         fsoieq(p)   = 1._r8
       end if
       fsoieq(p)   = 1._r8      ! nothing happens

       ! Tracer hydrology
       do m = 1, pwtrc

         ! Do a first guess, as the simplest default case
         ! Grab isotope ratio of first snow or soil layer with enough water to count
         RAtmVapor(m) = get_wratio(forc_wtr_q(c,m),forc_wtr_q(c,ixbase),voltiny)
         if(RAtmVapor(m)==0._r8 .and. qflx_evap_soi(p).ne.0._r8) RAtmVapor(m)=get_wratio(wtr_dqh(p,m),wtr_dqh(p,ixbase),voltiny)
         RSoilWater(m) = RAtmVapor(m)   ! initialize
         searchforwater: do j=snl(c)+1, nlevsoi
            if(wtr_h2osoi_liq(c,j,ixbase)+wtr_h2osoi_ice(c,j,ixbase) >= watmin) then
               ! not exactly the ratio of the water, but it should do - TW
               RSoilWater(m) = get_wratio(wtr_h2osoi_liq(c,j,m) + &
                               wtr_h2osoi_ice(c,j,m),wtr_h2osoi_liq(c,j,ixbase) + &
                               wtr_h2osoi_ice(c,j,ixbase),h2otiny)
               exit searchforwater
            endif
         end do searchforwater

         ! Craig-Gordon...
         !  regular evaporation (eq and kin fractionation); start off constant ratio
         wtr_qflx_evap_soi(p,m) = RSoilWater(m)*qflx_evap_soi(p)
         if     (qflx_evap_soi(p) > 0._r8 .and. t_grnd(c) >  tfrz) then
            alpeq = wiso_alpliq(ispec(m),t_grnd(c),.true.)
            wtr_qflx_evap_soi(p,m) = wtr_raiw(p,m)*( (fsoieq(p)*RSoilWater(m)*qg(c)/alpeq) + &
                                     ((1._r8-fsoieq(p))*forc_wtr_q(c,m)) - forc_wtr_q(c,m))
         !  'sublimation' of soil liquid (no fractionation)
         elseif (qflx_evap_soi(p) > 0._r8 .and. t_grnd(c) <= tfrz) then
            alpeq = 1._r8
            wtr_qflx_evap_soi(p,m) = wtr_raiw(p,ixbase)*( (fsoieq(p)*RSoilWater(m)*qg(c)/alpeq) + &
                                     ((1._r8-fsoieq(p))*forc_wtr_q(c,m)) - forc_wtr_q(c,m))
         !  dew/frost (eq and kin fractionation)
         elseif (qflx_evap_soi(p) < 0._r8) then
            alpeq = wiso_alpliq(ispec(m),t_grnd(c),.true.)
            wtr_qflx_evap_soi(p,m) = wtr_raiw(p,m)*( (fsoieq(p)*RSoilWater(m)*qg(c)/alpeq) + &
                                     ((1._r8-fsoieq(p))*forc_wtr_q(c,m)) - forc_wtr_q(c,m))

         end if

         ! Take into account partial equilibration, as is calulated in
         ! HydroTracersCanopy (just like canopy water done below)
         ! Notice this is pft dependent because of fsoieq. Sum to column later, 
         ! where it will be used to actually update the soil water.

         wtr_qflx_evap_tot(p,m) = wtr_qflx_evap_soi(p,m)
       end do   ! end 'tracer hydrology'

       ! 2 m height air temperature

       t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

       ! 2 m height specific humidity

       q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))
       do m = 1, pwtrc
         wtr_q_ref2m(p,m) = forc_wtr_q(c,m) + temp2(p)*wtr_dqh(p,m)*(1./temp22m(p) - 1./temp2(p))
       end do

       ! 2 m height relative humidity
                                                                                
       call QSat(t_ref2m(p), forc_pbot(c), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

       rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)

       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
         rh_ref2m_r(p) = rh_ref2m(p)
         t_ref2m_r(p) = t_ref2m(p)
       end if

       ! Variables needed by history tape

       t_veg(p) = forc_t(c)
       btran(p) = 0._r8
       cf = forc_pbot(c)/(SHR_CONST_RGAS*0.001_r8*thm(p))*1.e06_r8
       rssun(p) = 1._r8/1.e15_r8 * cf
       rssha(p) = 1._r8/1.e15_r8 * cf
       raw1(p) = 1._r8/1.e15_r8 * cf
       rlitter(p) = 1._r8/1.e15_r8 * cf
       rlb(p) = 1._r8/1.e15_r8 * cf

       ! Add the following to avoid NaN

       psnsun(p) = 0._r8
       psnsha(p) = 0._r8
       fpsn(p) = 0._r8
       pps%lncsun(p) = 0._r8
       pps%lncsha(p) = 0._r8
       pps%vcmxsun(p) = 0._r8
       pps%vcmxsha(p) = 0._r8
       ! adding code for isotopes, 8/17/05, PET
       pps%cisun(p) = 0._r8
       pps%cisha(p) = 0._r8
       if (use_c13) then
          pps%alphapsnsun(p) = 0._r8
          pps%alphapsnsha(p) = 0._r8
          pepv%rc13_canair(p) = 0._r8
          pepv%rc13_psnsun(p) = 0._r8
          pepv%rc13_psnsha(p) = 0._r8
          pc13f%psnsun(p) = 0._r8
          pc13f%psnsha(p) = 0._r8
          pc13f%fpsn(p) = 0._r8
       end if

       ! Add RXylem for bare patches too - need for averaging
       ! rootr should = 0 for all layers and bare pfts, but just make sure...
       do m = 1, pwtrc
         allr = 0._r8
         sumx = 0._r8
         summ = 0._r8
         do j = 1, nlevsoi
           RSoilLiq_tmp(p,j,m) = get_wratio(wtr_h2osoi_liq(c,j,m),wtr_h2osoi_liq(c,j,ixbase),h2otiny)
           sumx = sumx + rootr(p,j)*wtr_h2osoi_liq(c,j,ixbase)
           summ = summ + rootr(p,j)*wtr_h2osoi_liq(c,j,m)
           allr = allr + rootr(p,j)
         end do
         if (allr < 0.01) then
           RXylem(p,m) = RSoilLiq_tmp(p,1,m)
         else
           summ = summ/allr
           sumx = sumx/allr
           RXylem(p,m) = get_wratio(summ, sumx, h2otiny)
         end if
         pws%RLeafWaterSunPsn(p,m) =RXylem(p,m)*pcf%psnsun(p)
         pws%RLeafWaterShaPsn(p,m) =RXylem(p,m)*pcf%psnsha(p)
         pws%RLeafWaterSunStdy(p,m)=RXylem(p,m)
         pws%RLeafWaterShaStdy(p,m)=RXylem(p,m)
         pws%RLeafWaterSunDiag(p,m)=RXylem(p,m)
         pws%RLeafWaterShaDiag(p,m)=RXylem(p,m)
         pws%RCelluloseSun(p,m)    =RXylem(p,m)
         pws%RCelluloseSha(p,m)    =RXylem(p,m)
         pws%RCelluloseSunPsn(p,m) =RXylem(p,m)*pcf%psnsun(p)
         pws%RCelluloseShaPsn(p,m) =RXylem(p,m)*pcf%psnsha(p)
       end do

    end do

  end subroutine BareGroundFluxes

end module BareGroundFluxesMod
