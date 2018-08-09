module HydrologyTracer

!-----------------------------------------------------------------------
!BOP
! !MODULE: HydrologyTracer
!
! (Notice name not "WaterTracer" which is used in CAM)
!
! !DESCRIPTION:
!   Handles all water tracer hydrology.
!
!   This implimentation works from the ground up to interface with the water 
!   isotopes in CAM and for working with the carbon isotopes like in isolsm.
!   Tracer code is generic to allow "tagging" experiment, and addition
!   of equeous chemistry. This code relies on an external module to 
!   handle isotope fractionation.
!
!   For now, tarcers are H2OTR, H218O, HDO. The H2OTR is an
!   exact replica of the normal water used in the model, and is used to
!   check that the tracer scheme is working exactly as expected.
!   Adding H217O and HTO is a trivial extension if this is implimented
!   using the indexing and registration code from CAM.
!
!   This module operates at a hierler level than the "type" of water
!   (for instance, isotopes), and as such will interface with the 
!   more primitive code when needed.
!
! !REVISION HISTORY
!
!   Created David Noone <dcn@colorado.edu> - Thu Jun 16 18:45:47 MDT 2005
!   Updated David Noone <dcn@colorado.edu> - Fri Jul 13 15:40:29 MDT 2007
!   Updated Tony Wong <anthony.e.wong@colorado.edu - Wed Aug 14 10:23:38 MDT 2013
!       (for CLM4 use)
!
!EOP
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
!-----------------------------------------------------------------------
  use abortutils, only: endrun
!-----------------------------------------------------------------------
  implicit none
  private
  save
!-----------------------------------------------------------------------
!
! Module variables
!
  logical, public :: hydro_tracer = .true.  ! (namelist) flag to use this module
  logical, parameter, public :: tracer_forc = .false.   ! is this CLM standalone, but with isotopic forcing?
!
! Public interfaces
!
  public :: HydrologyTracerInit             ! Initializes module (sets tracers)
  public :: HydrologyTracerRegister         ! Sets tracer specific array values
  public :: HydrologyTracerRatios           ! Assigns ratios of pools
  public :: HydrologyTracerCanopy           ! Computes ratios for canopy fluxs
  public :: HydrologyTracerSoilWater        ! Updates tracer soilwater profile
  public :: HydrologyTracerRescale          ! Rescales state variables to match prognostic
  public :: HydrologyTracerCheck            ! Check all tracer arrays
  public :: TracerCheckPFT                  ! Check pft fluxes mostly
  public :: TracerCheckColumn               ! Check column state and fluxes
  public :: TracerCheckGrid                 ! Check grid level input/output
  public :: TracerCheckPFTDelta             ! Check pft fluxes mostly delta values
  public :: TracerCheckColumnDelta          ! Check column state and fluxes delta values
  public :: get_wratio                      ! computes tracer ratio
  public :: TracerCheckEqual                ! Check if 2 1D arrays are equal
!
! Module variables
!
  integer, parameter, public :: pwtrc = 3   ! number of tracers

! Indices of tracers in tracer array
!
  integer , parameter, public :: ixbase = 1         ! index of tracer same as prognostic
  integer , parameter, public :: ixtest = ixbase    ! index of tracer used in test
!  integer , parameter, public :: ixtest = 2		! isotope and NOFRAC
!  integer , parameter, public :: ixtest = 3		! isotope and NOFRAC

!
! Useful trivaially small values (for demoninator in ratio calulations)
! (Balance between smaller for more accurate and bigger for more robust ratios)
!
  real(r8), parameter, public :: h2otiny = 1.e-10_r8   ! mm (kg/m2)
  real(r8), parameter, public :: voltiny = 1.e-10_r8   ! m3/m3 (kg/kg)
  real(r8), parameter, public :: qflxtiny = 1.e-12_r8  ! mm/s (kg/m2/s)
!
! Array checking variables
!
  logical , parameter :: ldorescale = .TRUE.   ! flag to rescale for numerical precision
  logical , parameter :: lchecknone = .false.   ! flag bypass checking (faster when confident)
  logical , parameter :: lcheckstop = .true.    ! flag to terminate run on failed check
  logical , parameter :: lcheckverb = .false.   ! flag verbose checking (diagnostics)
  logical , parameter :: ldeltadiag = .true.    ! flag to output delta values for all variables
!
! Tolerances for error checking. Smaller is stricter.
!
  real(r8), parameter :: etols = 1.0e-6     ! check error tolerance on state (mm, kg/m2)
  real(r8), parameter :: etolf = 1.0e-9     ! check error tolerance on fluxes (kg/m2/sec)
  real(r8), parameter :: etolr = 1.0e-6     ! check error tolerance on ratio (kg/kg,m3/m3)
  real(r8), parameter :: errtol = 1.e-6     ! allowed rescaling correction (mm, kg/m2)

!
! Physics control
!
  real(r8), parameter :: taulef = 7200.    ! leaf water adjustment time scale
  real(r8), parameter :: taucan = 7200.    ! canopy vapor adjustment time scale
  real(r8), parameter :: taumin = 1800.     ! minimum adjustment time scale

!  real(r8), parameter :: flim = 0.33333     ! limit fraction of soil flux of available
  real(r8), parameter :: flim = 1._r8     ! for debugging
!
! Indicies useful for code development
!
  integer, parameter, public :: gdbg = 34      ! gridcell index
  integer, parameter, public :: ludbg = 50      ! landunit index
  integer, parameter, public :: cdbg = 1      ! column index
  integer, parameter, public :: pdbg = 2      ! pft index
  integer, parameter, public :: pdbgi= 1     ! pft starting index
  integer, parameter, public :: pdbgf= 17     ! pft ending index
!  integer, parameter, public :: cdbg = 777      ! column index
!  integer, parameter, public :: pdbg = 3611     ! pft index
!  integer, parameter, public :: pdbgi= 3609     ! pft starting index
!  integer, parameter, public :: pdbgf= 3625     ! pft ending index
  integer, parameter, public :: mdbg = 2        ! tracer index
  integer, parameter, public :: jdbg = -1        ! soil level index

  real(r8), parameter, public :: rsfcmin = 1.e-3_r8  ! minimum surface boundary layer resistance [s/m]
!
! Tracer names to append to variable names on output tape
!
  character(len=8), public :: wtrcnam(pwtrc)
!
! Standard ratio for numerics (need not be natural abundance)
!  Notice index for ixbase must be 1.0.
!
  real(r8), public :: Rstnd(pwtrc)
!
! Natural abundance, say, Rsmow (only used to normalize inputs if needed)
!
  real(r8), public :: Rntrl(pwtrc)
!
! Species identifier - what "type" of water. 
! Use of this means you have another module to work out the differences
! between species. Isotopes for instance.
!
  integer , public :: ispec(pwtrc)
!
!=======================================================================

CONTAINS

!=======================================================================
  subroutine HydrologyTracerInit()
!-----------------------------------------------------------------------
!  Initializes module by registering the tracer names and a standard
!  ratio.
!-----------------------------------------------------------------------
    use clm_varcon, only: spval
    use HydrologyIsotope
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
    write(6,*) 'HydrologyTracerInit'
    if (.not. hydro_tracer) call endrun('HydroTracer not set')

    Rstnd(:) = spval

    if (hydro_isotope) then
      if(pwtrc==1) then
        call HydrologyTracerRegister(1,'H2OTR', 1._r8, wiso_get_rnat(ispwater),ispwater) ! BASE 
      elseif(pwtrc==2) then
        call HydrologyTracerRegister(1,'H2OTR', 1._r8, wiso_get_rnat(ispwater),ispwater) ! BASE 
        call HydrologyTracerRegister(2,'HDO'  , 1._r8, wiso_get_rnat(isphdo)  ,isphdo  )
      elseif(pwtrc==3) then
!        call HydrologyTracerRegister(1,'H2OTR', 1._r8, wiso_get_rnat(ispwater),ispwater) ! BASE 
!        call HydrologyTracerRegister(2,'HDO'  , 1._r8, wiso_get_rnat(isphdo)  ,isphdo  )
!        call HydrologyTracerRegister(3,'H218O', 1._r8, wiso_get_rnat(isph218o),isph218o)
        call HydrologyTracerRegister(1,'H2OTR', wiso_get_rstd(ispwater), wiso_get_rnat(ispwater),ispwater) ! BASE 
        call HydrologyTracerRegister(2,'HDO'  , wiso_get_rstd(isphdo)  , wiso_get_rnat(isphdo)  ,isphdo  )
        call HydrologyTracerRegister(3,'H218O', wiso_get_rstd(isph218o), wiso_get_rnat(isph218o),isph218o)
      end if

    else
      if(pwtrc==1) then
        call HydrologyTracerRegister(1,'H2OTR', 1._r8, 1._r8, 0) ! BASE 
      elseif(pwtrc==2) then
        call HydrologyTracerRegister(1,'H2OTR', 1._r8, 1._r8, 0) ! BASE 
        call HydrologyTracerRegister(2,'HTR01', 1._r8, 1._r8, 0)
      elseif(pwtrc==3) then
        call HydrologyTracerRegister(1,'H2OTR', 1._r8, 1._r8, 0) ! BASE 
        call HydrologyTracerRegister(2,'HTR01', 1._r8, 1._r8, 0)
        call HydrologyTracerRegister(3,'HTR02', 1._r8, 1._r8, 0)
      end if
    end if

    if (count(Rstnd == spval) /= 0) then
      call endrun('(HydrologyTracerInit) Failed to register all tracers.')
    end if

    ! Check base tracer index
    write(6,*) 'base tracer name: '//wtrcnam(ixbase)
    if (Rstnd(ixbase) /= 1.0_r8) then
      call endrun('(HydrologyTracerIni) Base ratio is not unity.')
    endif

    write(6,*) 'HydrologyTracerInit: done.'
    return
  end subroutine HydrologyTracerInit


!=======================================================================
  subroutine HydrologyTracerRegister(m,vname,rs,rn,isp)
!-----------------------------------------------------------------------
! Assigns name and ratio to module arrays
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
    integer, intent(in) :: m	              ! tracer index
    character(len=*), intent(in) :: vname     ! tracer name
    real(r8), intent(in) :: rs                ! tracer standard ratio
    real(r8), intent(in) :: rn                ! tracer natural abundance ratio
    integer , intent(in) :: isp               ! tracer species index
!-----------------------------------------------------------------------
    if (m <     1) call endrun('Given tracer index < 1')
    if (m > pwtrc) call endrun('Given tracer index > pwtrc')

    wtrcnam(m) = vname
    Rstnd(m)   = rs
    Rntrl(m)   = rn
    ispec(m)   = isp

    write(6,1) m,wtrcnam(m),ispec(m),Rstnd(m),rntrl(m)
 1  format('  --> HydrologyTracerRegister: ',i3,x,a8,i3, &
                           '  ( Rstd=',EN14.4,' Rnat=',EN14.4,')')
    return
  end subroutine HydrologyTracerRegister


!=======================================================================
  subroutine HydrologyTracerRatios(fn, filterp, lbp, ubp, &
                     RCanopyWater, RSnowWater, RXylemWater, &
                     RSoilLiq, RSoilIce, RSoilVol )
!-----------------------------------------------------------------------
!
! Computes the tracer/prognostic ratio of various water pools.
! If your are looking for  the isotopic composition of canopy air and 
! leaf water, look in the HydrologyIsotope. 
! Ratios are defined at PFT level since they are needed for fluxes. 
! This is despite the fact that the pools are defied at column level
! (except for the canopy liquid).
!
!-----------------------------------------------------------------------
   use clmtype
   use clm_varpar,       only: nlevsoi
   use clm_varcon,       only: tfrz
   use HydrologyIsotope, only: wiso_alpliq, wiso_alpice
!-----------------------------------------------------------------------
   implicit none
!------------------------- Input Arguments -----------------------------
   integer , intent(in)  :: fn                   ! number of points 
   integer , intent(in)  :: filterp(:)            ! PFT point indicies
   integer, intent(in)   :: lbp, ubp	              ! PFT array bounds
!------------------------ Output Arguments -----------------------------
   real(r8), intent(out) :: RCanopyWater(lbp:ubp,pwtrc)  ! canopy water
   real(r8), intent(out) :: RSnowWater(lbp:ubp,pwtrc)    ! snow water
   real(r8), intent(out) :: RXylemWater(lbp:ubp,pwtrc)   ! water takeb by roots
   real(r8), intent(out) :: RSoilLiq(lbp:ubp,nlevsoi,pwtrc)    ! soil liquid
   real(r8), intent(out) :: RSoilIce(lbp:ubp,nlevsoi,pwtrc)    ! soil ice
   real(r8), intent(out) :: RSoilVol(lbp:ubp,nlevsoi,pwtrc)    ! net soil content
!-----------------------------------------------------------------------
!
! Local pointers to CLM type variable
!
   integer , pointer :: pcolumn(:)            !pft's column index
   real(r8), pointer :: rootr_pft(:,:)        !effective fraction of roots in each soil layer

   real(r8), pointer :: t_veg(:)              !Vegetation temperature (K)

   real(r8), pointer :: h2ocan(:)             !canopy water (mm H2O)
   real(r8), pointer :: h2osno(:)             !snow water (mm H2O)
   real(r8), pointer :: h2osoi_liq(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: h2osoi_ice(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: h2osoi_vol(:,:)       !volumetric soil water [m3/m3]  (nlevsoi)

   real(r8), pointer :: wtr_h2ocan(:,:)       !tracer canopy tracer water (mm H2O)
   real(r8), pointer :: wtr_h2osno(:,:)       !tracer snow water (mm H2O)
   real(r8), pointer :: wtr_h2osoi_liq(:,:,:) !tracer liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: wtr_h2osoi_liq_save(:,:,:) !tracer liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: wtr_h2osoi_ice(:,:,:) !tracer ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: wtr_h2osoi_vol(:,:,:) !tracer volumetric soil water  [m3/m3]  (nlevsoi)
!
   real(r8), pointer :: RCanopyVapor (:,:) ! vapour in canopy air space
   real(r8), pointer :: RXylem (:,:) ! xylem water, saved for hist files
!
   real(r8) allr		! sum of root weight
   real(r8) sumx,summ  	        ! sum of soil water for tracer and prognostic
   real(r8) alpeq               ! isotope fractionation factor
!
   integer f,p,c,j,m 		! array indicies
!
    logical :: ldbg=.false.
!    if(lbp<=pdbg .and. pdbg<=ubp) ldbg=.true.
!-----------------------------------------------------------------------

! Assign local pointers to CLM type variables

   pcolumn           => pft%column
   rootr_pft         => pps%rootr

   t_veg             => pes%t_veg

   h2ocan            => pws%h2ocan
   h2osno            => cws%h2osno
   h2osoi_liq        => cws%h2osoi_liq
   h2osoi_ice        => cws%h2osoi_ice
   h2osoi_vol        => cws%h2osoi_vol
!
   wtr_h2ocan        => pws%wtr_h2ocan
   wtr_h2osno        => cws%wtr_h2osno
   wtr_h2osoi_liq    => cws%wtr_h2osoi_liq
   wtr_h2osoi_liq_save    => cws%wtr_h2osoi_liq_save
   wtr_h2osoi_ice    => cws%wtr_h2osoi_ice
   wtr_h2osoi_vol    => cws%wtr_h2osoi_vol

   RCanopyVapor      => pws%RCanopyVapor
   RXylem            => pws%RXylem

!dir$ concurrent
!cdir nodep
   do m = 1, pwtrc
     do f  = 1, fn
         p = filterp(f)
         c = pcolumn(p)

         ! SNOW
         RSnowWater(p,m)   = get_wratio(wtr_h2osno(c,m),wtr_h2osno(c,ixbase),h2otiny)

         ! CANOPY INTERCEPTED WATER
         RCanopyWater(p,m) = get_wratio(wtr_h2ocan(p,m),wtr_h2ocan(p,ixbase),h2otiny)

         ! When no water on canopy, assigne value in equilibrium with
         ! canopy. This gives a useful "initial" value in case of dew
         ! onto dry canopy with fractional equilibration.
         if (abs(wtr_h2ocan(p,ixbase)) < h2otiny) then
           if (t_veg(p) > tfrz) then
!            alpeq = get_alpliq(ispec(m),t_veg(p),.true.)  ! look-up
             alpeq = wiso_alpliq(ispec(m),t_veg(p),.true.) ! compute
           else
!            alpeq = get_alpice(ispec(m),t_veg(p),.true.)  ! look-up
             alpeq = wiso_alpice(ispec(m),t_veg(p),.true.) ! compute
           endif
           RCanopyWater(p,m) = alpeq*RCanopyVapor(p,m)
         end if

         ! SOIL, AND XYLEM
         allr = 0._r8
         sumx = 0._r8
         summ = 0._r8
         do j = 1, nlevsoi
           RSoilLiq(p,j,m) = get_wratio(wtr_h2osoi_liq(c,j,m),wtr_h2osoi_liq(c,j,ixbase),h2otiny)
           RSoilIce(p,j,m) = get_wratio(wtr_h2osoi_ice(c,j,m),wtr_h2osoi_ice(c,j,ixbase),h2otiny)
           RSoilVol(p,j,m) = get_wratio(wtr_h2osoi_vol(c,j,m),wtr_h2osoi_vol(c,j,ixbase),voltiny)
           sumx = sumx + rootr_pft(p,j)*wtr_h2osoi_liq(c,j,ixbase)
           summ = summ + rootr_pft(p,j)*wtr_h2osoi_liq(c,j,m)
           allr = allr + rootr_pft(p,j)
           wtr_h2osoi_liq_save(c,j,m) = wtr_h2osoi_liq(c,j,m)
         end do

         if (allr < 0.01) then
           RXylemWater(p,m) = RSoilLiq(p,1,m)
         else
           summ = summ/allr
           sumx = sumx/allr
           RXylemWater(p,m) = get_wratio(summ, sumx, h2otiny)
         end if
         RXylem(p,m)=RXylemWater(p,m)
     end do
   end do

   return
  end subroutine HydrologyTracerRatios


!=======================================================================
  subroutine HydrologyTracerCanopy(fn, filterp, lbp, ubp, lbg, ubg, &
                       qaf, qsatl   , wtaq0   , wtgq0 , &
                       wtlsunq0, wtlshaq0, wtlcanq0, rb    , forc_pbot , forc_t , &
                       RXylemWater, RCanopyWater, RSoilVapor , forc_rho , wtl,wta,wtal)
!-----------------------------------------------------------------------
!
! Computes the tracer ratio in canopy air, and an effective ratio
! for vapour above leaves (leaf water plus canopy liquid). 
! Notice that this is really just a placeholder for the more useful
! calculation that would be done for isotopes.
!
! The canopy vapour and leave water assume steady state, much like the
! Craig and Gordon (1965) assumptions. This the values predicted here
! can in fact form the first part of the isotope calulations.
!
!-----------------------------------------------------------------------
   use clmtype
   use clm_varcon      , only: tfrz
   use clm_time_manager, only: get_step_size, get_nstep
   use HydrologyIsotope, only: wiso_alpliq, wiso_alpice, wiso_alpkin_m78
!-----------------------------------------------------------------------
   implicit none
!------------------------- Input Arguments -----------------------------
   integer , intent(in)  :: fn			  ! number of filtered points 
   integer , intent(in)  :: filterp(:)            ! filterd PFT point indicies
   integer , intent(in)  :: lbp, ubp              ! lower and upper pft bounds
   integer , intent(in)  :: lbg, ubg              ! lower and upper grid bounds
   real(r8), intent(in)  :: qaf(lbp:ubp)          ! canopy specific humidity [kg/kg]
   real(r8), intent(in)  :: qsatl(lbp:ubp)        ! leaf specific humidity [kg/kg]
   real(r8), intent(in)  :: rb(lbp:ubp)           ! leaf boundary layer resistance [s/m]
!   real(r8), intent(in)  :: wtaq0(lbp:ubp)        ! normalized latent heat conductance for air [-]
!   real(r8), intent(in)  :: wtgq0(lbp:ubp)        ! normalized conductance for ground [-]
!   real(r8), intent(in)  :: wtlsunq0(lbp:ubp)     ! normalized conoductance for sunlit leaves
!   real(r8), intent(in)  :: wtlshaq0(lbp:ubp)     ! normalized conductance for shaded leaves
!   real(r8), intent(in)  :: wtlcanq0(lbp:ubp)     ! normalized conductance for wet canopy
   real(r8), intent(in)  :: wtaq0(lbp:ubp,pwtrc)        ! normalized latent heat conductance for air [-]
   real(r8), intent(in)  :: wtgq0(lbp:ubp,pwtrc)        ! normalized conductance for ground [-] CRT
!  real(r8), intent(in)  :: wtgq0(pwtrc)        ! normalized conductance for ground [-]
   real(r8), intent(in)  :: wtlsunq0(lbp:ubp,pwtrc)     ! normalized conoductance for sunlit leaves
   real(r8), intent(in)  :: wtlshaq0(lbp:ubp,pwtrc)     ! normalized conductance for shaded leaves
   real(r8), intent(in)  :: wtlcanq0(lbp:ubp,pwtrc)     ! normalized conductance for wet canopy
   real(r8), intent(in)  :: forc_pbot(lbg:ubg)    ! atmospheric pressure (Pa)
   real(r8), intent(in)  :: forc_t(lbg:ubg)       ! atmospheric temperature (K)
   real(r8), intent(in)  :: forc_rho(lbg:ubg)     ! atmospheric density
   real(r8), intent(in)  :: RXylemWater  (lbp:ubp,pwtrc) ! water at evaporating surface
   real(r8), intent(in)  :: wtl(lbp:ubp) ! wtlq0(:) from CanopyFluxes?
   real(r8), intent(in)  :: wta(lbp:ubp) ! wtaq0(:) from CanFlx?
   real(r8), intent(in)  :: wtal(lbp:ubp)! wtalq(:) from CanFlx
!---------------------- Input/Output Arguments -------------------------
   real(r8), intent(inout) :: RCanopyWater(lbp:ubp,pwtrc)  ! intercepted liquid
   real(r8), intent(inout) :: RSoilVapor(lbp:ubp,pwtrc)    ! Vapour in top soil
!------------------------- Local Variables -----------------------------
!
! Local pointers to clm type variables: (in) input
!
   integer , pointer :: pcolumn(:)         ! pft's column index
   integer , pointer :: pgridcell(:)       ! pft's gridcell index
   real(r8), pointer :: wtcol(:)           ! pft fractional contribution to each column

   real(r8), pointer :: t_veg(:)           ! vegetation temperature (Kelvin)
   real(r8), pointer :: rssun(:)           ! sunlit stomatal resistance (s/m)
   real(r8), pointer :: rssha(:)           ! shaded stomatal resistance (s/m)
   real(r8), pointer :: fcaneq(:)          ! canopy water equilibration factor
   real(r8), pointer :: fsoieq(:)          ! soil water equilibration factor
   real(r8), pointer :: psnsun(:)          ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha(:)          ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: laisun(:)          ! sunlit leaf area
   real(r8), pointer :: laisha(:)          ! shaded leaf area

   real(r8), pointer :: qg(:)              ! specific humidity at ground surface [kg/kg]
   real(r8), pointer :: wtr_qg(:,:)        ! specific humidity at ground surface [kg/kg]
   real(r8), pointer :: forc_q(:)          ! atmospheric specific humidity (kg/kg)
   real(r8), pointer :: forc_wtr_q(:,:)    ! atmospheric specific humidity (kg/kg)
   real(r8), pointer :: qflx_tran_veg(:)   ! transpiration flux (mm H2O/s)
!
! Local pointers to clm type variables: (inout) output
!
   real(r8), pointer :: leaf_mr(:)             ! leaf maintenance respiration
   real(r8), pointer :: RCelluloseSun(:,:) ! tracer ratio of sunlit cellulose (incl Peclet)
   real(r8), pointer :: RCelluloseSha(:,:) ! tracer ratio of shaded cellulose (incl Peclet)
   real(r8), pointer :: RCelluloseSunPsn(:,:) ! tracer ratio of sunlit cellulose times photosynthesis (incl Peclet)
   real(r8), pointer :: RCelluloseShaPsn(:,:) ! tracer ratio of shaded cellulose times photosynthesis (incl Peclet)
   real(r8), pointer :: RLeafWaterSun(:,:) ! water at evaporating surface
   real(r8), pointer :: RLeafWaterSha(:,:) ! water at evaporating surface
   real(r8), pointer :: RLeafWaterSunPsn(:,:) ! tracer ratio of sunlit leaf water at evaporating surface times photosynthesis (incl Peclet)
   real(r8), pointer :: RLeafWaterShaPsn(:,:) ! tracer ratio of shaded leaf water at evaporating surface times photosynthesis (incl Peclet)
   real(r8), pointer :: RLeafWaterSunStdy(:,:) ! steady state leaf water at evaporating surface (no Peclet) 
   real(r8), pointer :: RLeafWaterShaStdy(:,:) ! steady state leaf water at evaporating surface (no Peclet)
   real(r8), pointer :: RLeafWaterSunDiag(:,:) ! water at evaporating surface (incl Peclet)
   real(r8), pointer :: RLeafWaterShaDiag(:,:) ! water at evaporating surface (incl Peclet)
   real(r8), pointer :: RCanopyVapor (:,:) ! vapour in canopy air space
   real(r8), pointer :: RCanopyWaterHist (:,:) ! tracer canopy liquid water ratio
!
! Local variables
!
   integer c,f,g,p,m                    ! indicies

   real(r8) dtime		        ! time step (seconds)

   real(r8) RCanopyVaporSteady          ! steady state (C-G) vapour in canopy air space

   real(r8) alpeq	                ! equilibrium fractionation for water
   real(r8) cgxysun	                ! Craig-Gordon weight for sunlit xylem water
   real(r8) cgxysha	                ! Craig-Gordon weight for shaded xylem water
   real(r8) cgcnvap	                ! Craig-Gordon weight for canopy vapor 
   real(r8) fex_Lsun	                ! implict factor for sun leaves
   real(r8) fex_Lsha	                ! implict factor for shaded leaves
   real(r8) fex_Vcan	                ! implict factor for canopy vapor

   real(r8) qafx                        ! diagnosed specific humidity in canopy (local version)
   real(r8) qlsun                       ! specific humidity depression at sun leaf surface "qs(1-h)"
   real(r8) qlsha                       ! specific humidity depression at shad leaf surface "qs(1-h)"
   real(r8) wtr_dqlsun                  ! tracr s.hum depression at sun leaf surface "qs(1-h)"
   real(r8) wtr_dqlsha                  ! tracer s. hum depression at shaded leaf surface "qs(1-h)"
   real(r8) alpk_lf, alpk_st            ! leaf boundary layer and leaf stomata kinetic fractionation factors, respectively

   real(r8) tau_vcan            ! residence time for canopy airspace
   real(r8) tau_lsun            ! residence time for sunlit leaf water
   real(r8) tau_lsha            ! residence time for shaded leaf water

   real(r8) numer,denom         ! numerator and denominator for canopy ratio
   real(r8) :: cf                       ! conversion factor for var_taulef
   real(r8) :: mass_l=10._r8            ! mass of leaf water per m2 of LAI

   real(r8) :: ll=0.015_r8         ! effective path length [m] (PECLET)
   real(r8) :: cc=5.55e4_r8        ! molar concentration of water [mol/m3] (PECLET)
   real(r8) :: pecsun,pecsha            ! Peclet number [] (PECLET)
   real(r8) :: des,dls                  ! evaporation site and leaf isotopic discriminations [] (PECLET)
   real(r8),dimension(pwtrc) :: &       ! diffusivity of tracers in water (H218O from Wang, 1954) [m2/s] (PECLET)
                    dd=(/2.30e-9_r8,2.66e-9_r8,2.66e-9_r8/)    ! (H216O from Easteal et al, 1984)
   real(r8) :: psnsun_frac(lbp:ubp),psnsha_frac(lbp:ubp)  ! sunlit and shaded fractions of photosynthesis
   real(r8) :: assun(lbp:ubp),assha(lbp:ubp)                ! assimilation (lai*psn-resp) (CELLULOSE)
   real(r8) :: RLeaf(lbp:ubp,pwtrc),DLeaf(lbp:ubp,pwtrc),DXylem(lbp:ubp,pwtrc) ! (CELLULOSE)
   real(r8),dimension(pwtrc) :: &       ! kinetic effect for cellulose heterotrophic metabolism
                    epsh=(/0._r8 , -171._r8 , 27._r8/)          ! Yakir and DeNiro
   real(r8),dimension(pwtrc) :: &       ! kinetic effect for cellulose autotrophic metabolism
                    epsa=(/0._r8 ,  158._r8 , 27._r8/)          ! Yakir and DeNiro, and Roden
   real(r8),dimension(pwtrc) :: &       ! fraction of 
                    fiso=(/0.5_r8 ,  0.36_r8 , 0.42_r8/)          ! Roden
   real(r8) :: fcell=0.5_r8             ! cellulose fraction of PSN actually used

   logical :: lvar_taulef=.false.        ! use variable tau_lef?
   integer :: idbg
   logical :: ldbg=.false.
!   if(ubp>=pdbg .and. lbp<=pdbg) ldbg=.true.
!-----------------------------------------------------------------------

   ! Assign local pointers of clmtype variables

   pcolumn        => pft%column
   pgridcell      => pft%gridcell
   wtcol          => pft%wtcol

   t_veg          => pes%t_veg
   rssun          => pps%rssun
   rssha          => pps%rssha
   fcaneq         => pps%fcaneq
   fsoieq         => pps%fsoieq
   psnsun         => pcf%psnsun
   psnsha         => pcf%psnsha
   laisun         => pps%laisun
   laisha         => pps%laisha
   qflx_tran_veg  => pwf%qflx_tran_veg

   qg             => cws%qg

   forc_wtr_q     => cws%forc_wtr_q  !clm4

   wtr_qg         => cws%wtr_qg

   leaf_mr      => pcf%leaf_mr
   RCelluloseSun    => pws%RCelluloseSun
   RCelluloseSha    => pws%RCelluloseSha
   RCelluloseSunPsn => pws%RCelluloseSunPsn
   RCelluloseShaPsn => pws%RCelluloseShaPsn
   RLeafWaterSunDiag=> pws%RLeafWaterSunDiag
   RLeafWaterShaDiag=> pws%RLeafWaterShaDiag
   RLeafWaterSunStdy=> pws%RLeafWaterSunStdy
   RLeafWaterShaStdy=> pws%RLeafWaterShaStdy
   RLeafWaterSun    => pws%RLeafWaterSun
   RLeafWaterSha    => pws%RLeafWaterSha
   RLeafWaterSunPsn => pws%RLeafWaterSunPsn
   RLeafWaterShaPsn => pws%RLeafWaterShaPsn
   RCanopyVapor     => pws%RCanopyVapor
   RCanopyWaterHist => pws%RCanopyWater

   ! Obtain model time step size

   dtime = get_step_size()

   ! Initialize arrays for "total" water
!dir$ concurrent
!cdir nodep
   do f = 1, fn
      p = filterp(f)
      c = pcolumn(p)
      g = pgridcell(p)

      ! Reestablish the canopy specific humidity using tracer variables
      ! (Must match wtr_qafx(ixbase), helps maintain numerical precision)
      qafx = (wtlcanq0(p,ixbase)+wtlsunq0(p,ixbase)+wtlshaq0(p,ixbase))*qsatl(p) + &
                  wtgq0(p,ixbase)*wtr_qg(c,ixbase) + wtaq0(p,ixbase)*forc_wtr_q(c,ixbase) !CRT
!                 wtgq0(p,ixbase)*wtr_qg(c,ixbase) + wtaq0(p,ixbase)*forc_wtr_q(g,ixbase) !BEN
!                 wtgq0(ixbase)*wtr_qg(c,ixbase) + wtaq0(p,ixbase)*forc_wtr_q(c,ixbase)

      ! Compute sturation vapour pressure at leaf surface boundary 
      ! (Conductance weighted average of canopy and leaf interior values)
      qlsun = (qafx*rssun(p) + qsatl(p)*rb(p))/(rssun(p) + rb(p))
      qlsha = (qafx*rssha(p) + qsatl(p)*rb(p))/(rssha(p) + rb(p))

      ! Set leaf water  residence time scale and implict adjustment factor
      tau_lsun = taulef
      tau_lsha = taulef

      if (lvar_taulef) then
         cf = forc_pbot(g) / (8.314 * forc_t(g))   ! convert to mol/m2/s
         tau_lsun = mass_l*rssun(p)*forc_pbot(g)/(cf*qlsun)
         tau_lsun = max(tau_lsun, taumin)
         tau_lsha = mass_l*rssha(p)*forc_pbot(g)/(cf*qlsha)
         tau_lsha = max(tau_lsha, taumin)
      end if

      fex_lsun = exp(-dtime/tau_lsun)
      fex_lsha = exp(-dtime/tau_lsha)

      ! Set canopy residence time scale and implict adjustment factor
      tau_vcan = taucan
      fex_vcan = exp(-dtime/tau_vcan)

      ! Initialize species dependent variables
      do m = 1, pwtrc

        if(t_veg(p) > tfrz) then
           alpeq = wiso_alpliq(ispec(m),t_veg(p),.true.)      ! compute
        else
           alpeq = wiso_alpice(ispec(m),t_veg(p),.true.)      ! compute
        end if

        ! Assign a leaf surface effective specific humidity using kinetic
        ! fractionation through stomate and leaf boundary layer.

        alpk_lf = wiso_alpkin_m78(ispec(m),.67_r8,.true.)
        alpk_st = wiso_alpkin_m78(ispec(m),1.0_r8,.true.)

        wtr_dqlsun = (qsatl(p) - qlsun)*alpk_st + &
                     (qlsun    - qafx )*alpk_lf
        wtr_dqlsha = (qsatl(p) - qlsha)*alpk_st + &
                     (qlsha    - qafx )*alpk_lf

        ! Set up weighting factors for solving the steady state leaf model 

        cgxysun = alpeq*wtr_dqlsun/qsatl(p)
        cgxysha = alpeq*wtr_dqlsha/qsatl(p)
        cgcnvap = alpeq*qafx      /qsatl(p)

        ! Steady state canopy vapor: Solve the three simultaneous equations... (YES)
        !    1) 5 way balance: 2x transp, ground+canopy evap, and atmosphere
        !    2) Graig-Gordon type mass balance for sun and shaded leaf water
        !    3) Partial equilibration of canopy intercepted liquid water
        !    4) Partial equilibration of upper layer soil water

        numer = ( wtlsunq0(p,m)*cgxysun  *RXylemWater (p,m)  + &
                  wtlshaq0(p,m)*cgxysha  *RXylemWater (p,m)  + &
                  wtlcanq0(p,m)*fcaneq(p)*RCanopyWater(p,m) )*qsatl(p)/alpeq + &
                  wtgq0(p,m)   *fsoieq(p)*RSoilVapor  (p,m)  *qg(c) + & !CRT
!                 wtgq0(m)   *fsoieq(p)*RSoilVapor  (p,m)  *qg(c) + &
!                 wtaq0(p,m)*forc_wtr_q(g,m) !BEN
                  wtaq0(p,m)*forc_wtr_q(c,m)

        denom = qafx    *(1.- wtlsunq0(p,m) - wtlshaq0(p,m)) - &
                qsatl(p)*(1.-fcaneq(p))*wtlcanq0(p,m) - &
                qg(c)   *(1.-fsoieq(p))*wtgq0(p,m) !CRT
!               qg(c)   *(1.-fsoieq(p))*wtgq0(m)

        RCanopyVaporSteady = get_wratio(numer,denom,voltiny)

        ! Steady state leaf water: "Craig and Gordon"
        RLeafWaterSunStdy(p,m) = cgxysun*RXylemWater(p,m) + &
                                 cgcnvap*RCanopyVaporSteady

        RLeafWaterShaStdy(p,m) = cgxysha*RXylemWater(p,m) + &
                                 cgcnvap*RCanopyVaporSteady

        ! Default to steady state. Update next if you want...
        RCanopyVapor(p,m) = RCanopyVaporSteady
        RLeafWaterSun(p,m) = RLeafWaterSunStdy(p,m)
        RLeafWaterSha(p,m) = RLeafWaterShaStdy(p,m)

        ! Apply implicit correction for non-steady state (prognostic) model
        ! Removed, replace with diagnostic leaf calculation (Peclet) (TW 26 Feb 2015)
#if 0
        RCanopyVapor(p,m)  = (        fex_vcan)*RCanopyVapor(p,m) + &
                             (1._r8 - fex_vcan)*RCanopyVaporSteady
        RLeafWaterSun(p,m) = (        fex_lsun)*RLeafWaterSun(p,m) + &
                             (1._r8 - fex_lsun)*RLeafWaterSunStdy(p,m)
        RLeafWaterSha(p,m) = (        fex_lsha)*RLeafWaterSha(p,m) + &
                             (1._r8 - fex_lsha)*RLeafWaterShaStdy(p,m)
#endif

        ! Peclet effect (diagnostic calculation for leaf water ratios)
        RLeafWaterSunDiag(p,m) = RLeafWaterSun(p,m)
        RLeafWaterShaDiag(p,m) = RLeafWaterSha(p,m)
        if(psnsun(p)> 0._r8) then 
           psnsun_frac(p)=psnsun(p)/(psnsun(p)+psnsha(p))
        else
           psnsun_frac(p)=0.5_r8
        end if
        if(laisun(p)> 0._r8) then 
           pecsun = psnsun_frac(p)*qflx_tran_veg(p)*ll / (laisun(p)*18e-3_r8*cc*dd(ispec(m)))  ! day
        else
           pecsun = psnsun_frac(p)*qflx_tran_veg(p)*ll / (18e-3_r8*cc*dd(ispec(m))) ! night
        end if
        if(laisha(p)> 0._r8) then 
           pecsha = (1._r8-psnsun_frac(p))*qflx_tran_veg(p)*ll / (laisha(p)*18e-3_r8*cc*dd(ispec(m)))  ! day
        else
           pecsha = (1._r8-psnsun_frac(p))*qflx_tran_veg(p)*ll / (18e-3_r8*cc*dd(ispec(m)))            ! night
        end if
        if(pecsun   /=0._r8) then    !sun    !EK 25 Jun 2015 (to "Cellulose model")
           if ( RXylemWater(p,m) == 0.0_r8 ) then
              RLeafWaterSunDiag(p,m) = 0.0_r8
           else
              des = RLeafWaterSun(p,m)/RXylemWater(p,m) -1._r8     ! isotopic discrimination w/o Peclet effect
              dls = des*(1._r8-exp(-pecsun))/pecsun                ! leaf water enrichment w/ Peclet effect
              RLeafWaterSunDiag(p,m) = (dls+1._r8)*RXylemWater(p,m)
           end if
        end if
        if(pecsha   /=0._r8) then    !sha
           if ( RXylemWater(p,m) == 0.0_r8 ) then
              RLeafWaterSunDiag(p,m) = 0.0_r8
           else
              des = RLeafWaterSun(p,m)/RXylemWater(p,m) -1._r8     ! isotopic discrimination w/o Peclet effect
              dls = des*(1._r8-exp(-pecsha))/pecsha                ! leaf water enrichment w/ Peclet effect
              RLeafWaterSunDiag(p,m) = (dls+1._r8)*RXylemWater(p,m)
           end if
        end if


        ! Cellulose model (Roden et al, 1999)
        !!assun(p) = laisun(p)*(psnsun(p) - leaf_mr(p))/(laisun(p)+laisha(p))
        !!assha(p) = laisha(p)*(psnsha(p) - leaf_mr(p))/(laisun(p)+laisha(p))
        assun(p) = psnsun(p)
        assha(p) = psnsha(p)
        DXylem(p,m)= 1000._r8*( (RXylemWater(p,m)/Rstnd(m)) -1._r8)             ! convert to delta
        ! ...sunlit leaves
        DLeaf(p,m) = 1000._r8*( (RLeafWaterSunDiag(p,m)/Rstnd(m)) -1._r8)       ! convert to delta
        RCelluloseSun(p,m) = fiso(ispec(m))*(DXylem(p,m)+epsh(ispec(m))) + &
                            (1._r8-fiso(ispec(m)))*(DLeaf(p,m)+epsa(ispec(m)))  ! this is actually a delta!
        RCelluloseSun(p,m) = (RCelluloseSun(p,m)/1000._r8 + 1._r8)*Rstnd(m)     ! convert to ratio
        ! ...shaded leaves
        DLeaf(p,m) = 1000._r8*( (RLeafWaterShaDiag(p,m)/Rstnd(m)) -1._r8)       ! convert to delta
        RCelluloseSha(p,m) = fiso(ispec(m))*(DXylem(p,m)+epsh(ispec(m))) + &
                            (1._r8-fiso(ispec(m)))*(DLeaf(p,m)+epsa(ispec(m)))  ! this is actually a delta!
        RCelluloseSha(p,m) = (RCelluloseSha(p,m)/1000._r8 + 1._r8)*Rstnd(m)     ! convert to ratio
        ! Photosynthesis-scaled leaf water isotope ratios, for the history tape
        RLeafWaterSunPsn(p,m) = RLeafWaterSunDiag(p,m)*psnsun(p)
        RLeafWaterShaPsn(p,m) = RLeafWaterShaDiag(p,m)*psnsha(p)
        RCelluloseSunPsn(p,m) = RCelluloseSun(p,m)*psnsun(p)
        RCelluloseShaPsn(p,m) = RCelluloseSha(p,m)*psnsha(p)

        ! Use the Peclet-included leaf water ratios?
#if 0
        RLeafWaterSun(p,m) = RLeafWaterSunDiag(p,m)
        RLeafWaterSha(p,m) = RLeafWaterShaDiag(p,m)
#endif

        ! Apply partial equilibration of canopy intercepted liquid with vapour
        ! There is no net exchange of mass, but the molecules involved change.

        RCanopyWater(p,m) = (        fcaneq(p))*RCanopyWater(p,m) + &
                            (1._r8 - fcaneq(p))*RCanopyVapor(p,m)*alpeq 

        ! Output for history files canopy water, in partial equilibrium with
        ! canopy vapor/canopy water
        !RCanopyWaterHist(p,m) = RCanopyVapor(p,m)*alpeq
        RCanopyWaterHist(p,m) = RCanopyWater(p,m)

        ! Apply partial eqilibration of soil water with canopy vapour
        ! As above, not net change in mass, just the molecules.

        RSoilVapor(p,m)   = (        fsoieq(p))*RSoilVapor(p,m) + &
                            (1._r8 - fsoieq(p))*RCanopyVapor(p,m)

      end do
   end do

   return
  end subroutine HydrologyTracerCanopy


!=======================================================================
  subroutine HydrologyTracerSoilWater(lbc,ubc,num_hydrologyc,filter_hydrologyc, &
                 h2osoi_old )
!-----------------------------------------------------------------------
!
! Backs out the soil water budget for tracer water. 
! Just does the soil water transport. Drainage is done separately.
! (This is because it is easier to test/code, not for science reasons)
!
! Solves for qflx given:
!
!  [dq/dt] = qflx(j+1) - qflx(j) + qsrc
!
!
! under the boundary condition that qflx(1) = qflx_infl
! (which is negative when evaporation dominates)
!
! Since this is an explict calulation it is conditionally stable, and is
! most conviniently stabilized by substepping. This however, means the 
! Xylem water ratio will only be aproximate as calulated in Biogeophys1.
!
!-----------------------------------------------------------------------
    use clmtype
    use clm_varpar       , only : nlevsoi
    use clm_time_manager , only : get_step_size
    use HydrologyIsotope , only : hydro_isotope, wiso_alpliq
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    integer , intent(in) :: lbc, ubc                ! column bounds
    integer , intent(in) :: num_hydrologyc               ! number of column soil points in column filter
    integer , intent(in) :: filter_hydrologyc(ubc-lbc+1) ! column filter for soil points
    real(r8), intent(in) :: h2osoi_old(lbc:ubc,nlevsoi)   ! initial soil water [kg/m2]
!-----------------------------------------------------------------------
!
! Implict input variables stored in clmtype
!
    integer , pointer :: pfti(:)                    !beginning pft index for each column
    integer , pointer :: pftf(:)                    !ending pft index for each column
    integer , pointer :: npfts(:)                   !column's number of pfts 

    integer , pointer :: cgridcell(:)		    !grid indicies for each column
    integer , pointer :: snl(:)                     !number of snow layers
    real(r8), pointer :: rootr_col(:,:)             !effective fraction of roots in each soil layer
    real(r8), pointer :: wtcol(:)                   !weight relative to column for each pft

    real(r8), pointer :: RGNIP(:,:)                 !tracer water ratio climatology
    real(r8), pointer :: t_grnd(:)                  !ground temperature (Kelvin)

    real(r8), pointer :: h2osoi_liq(:,:)            !liquid water (kg/m2)
    real(r8), pointer :: wtr_h2osoi_liq(:,:,:)      !tracer liquid water (kg/m2)
    real(r8), pointer :: wtr_h2osoi_liq_save(:,:,:)      !tracer liquid water (kg/m2)

    real(r8), pointer :: qflx_drain(:)              !sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_infl(:)               !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_tran_veg_col(:)       !vegetation transpiration (mm H2O/s) (+ = to atm)


    real(r8), pointer :: wtr_qflx_drain(:,:)        !sub-surface runoff (mm H2O /s)
    real(r8), pointer :: wtr_qflx_infl(:,:)         !infiltration (mm H2O /s)
    real(r8), pointer :: wtr_qflx_tran_veg_pft(:,:) !(pft) vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_tran_veg_col(:,:) !(col) vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_tot_col(:,:) !Col average Total evaporation (passed to atmosphere)
    real(r8), pointer :: wtr_qflx_evap_soi_col(:,:) !Col average Evaporation from soil

    real(r8), pointer :: RCanopyVapor_pft(:,:)      ! vapour in canopy air space PFT level
    real(r8), pointer :: RCanopyVapor_col(:,:)      ! vapour in canopy air space column level
    real(r8), pointer :: fsoieq_pft(:)              ! soil water tracer equilibration fraction
    real(r8), pointer :: fsoieq_col(:)              ! soil water tracer equilibration fraction
!
! Local variables
!
    integer fc, c, g, j, m, p                   ! indices
    integer  :: isub                         ! substep counter
    integer  :: nsub                         ! number of substeps
    real(r8) :: dtsub                        ! duration of substeps
    real(r8) :: dtime                        ! land model time step (sec)
    real(r8) :: delh2o(lbc:ubc,nlevsoi)      ! change in soil liquid [kg/m2]
    real(r8) :: wtr_delh2o(lbc:ubc,nlevsoi)  ! change in soil liquid [kg/m2]
    real(r8) :: h2onew(lbc:ubc,nlevsoi)      ! substepped h2osoi_liq [kg/m2]
    real(r8) :: qflx_soi(lbc:ubc,nlevsoi+1)  ! interface water fluxes [kg/m2/sec]
    real(r8) :: RFlx(lbc:ubc,nlevsoi+1,pwtrc)! tracer ratio of interface fluxes
    real(r8) :: RSoil(lbc:ubc,nlevsoi,pwtrc) ! tracer ratio of soil liquid
    real(r8) :: wratio(pwtrc)                ! water tracer ratio for drainage
    real(r8) :: Rnew(pwtrc)                  ! new isotope ratio
    real(r8) :: dflux			     ! flux divergence
    real(r8) :: maxdt
    logical :: ldbg=.false.
!    if(lbc<=cdbg .and. ubc>=cdbg) ldbg=.true.
!-----------------------------------------------------------------------

    ! Assign local pointers to CLM type

    npfts                 => grc%npfts
    pfti                  => col%pfti
    pftf                  => col%pftf
    wtcol                 => pft%wtcol
    cgridcell             => col%gridcell
    rootr_col             => cps%rootr_column
    snl                   => cps%snl

    t_grnd                => ces%t_grnd
    RGNIP                 => cws%RGNIP

    h2osoi_liq            => cws%h2osoi_liq
    wtr_h2osoi_liq        => cws%wtr_h2osoi_liq
    wtr_h2osoi_liq_save        => cws%wtr_h2osoi_liq_save

    qflx_drain            => cwf%qflx_drain
    qflx_infl             => cwf%qflx_infl
    qflx_tran_veg_col     => pwf_a%qflx_tran_veg

    wtr_qflx_tran_veg_pft => pwf%wtr_qflx_tran_veg

    wtr_qflx_drain        => cwf%wtr_qflx_drain
    wtr_qflx_infl         => cwf%wtr_qflx_infl
    wtr_qflx_tran_veg_col => pwf_a%wtr_qflx_tran_veg
    wtr_qflx_evap_soi_col => pwf_a%wtr_qflx_evap_soi
    wtr_qflx_evap_tot_col => pwf_a%wtr_qflx_evap_tot

    RCanopyVapor_pft      => pws%RCanopyVapor
    fsoieq_pft            => pps%fsoieq
    RCanopyVapor_col      => pws_a%RCanopyVapor
    fsoieq_col            => pps_a%fsoieq

    ! Get time step

    dtime = get_step_size()

    ! Assign a timestep for substepping.
    dtsub = dtime
!    dtsub = dtime/4.		! match with flim
!    dtsub = 600.           ! 10 min
!    dtsub = 180.            ! 3 min
!    dtsub = 60.            ! 1 min
!    dtsub = 30.            ! 30 sec
    dtsub = 6.            ! 6 sec
    if (mod(dtime,dtsub) /= 0) call endrun('HydrologyTracer: non integer substeps')
    nsub = nint(dtime/dtsub)

    ! Back out interface fluxes from (corrected) tendency 
    ! Note sign convention: < 0 is up.

!dir$ concurrent
!cdir nodep
    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      qflx_soi(c,1) = 0._r8
      delh2o(c,1) = (h2osoi_liq(c,1)-h2osoi_old(c,1))/dtime
      qflx_soi(c,2) = qflx_infl(c) - qflx_tran_veg_col(c)*rootr_col(c,1) - delh2o(c,1)
    end do
    do j = 2, nlevsoi
!dir$ concurrent
!cdir nodep
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         delh2o(c,j) = (h2osoi_liq(c,j)-h2osoi_old(c,j))/dtime
         qflx_soi(c,j+1) = qflx_soi(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j) - delh2o(c,j)
      end do
    end do

    ! zero the drainage flux array so we can get it by summation
    ! set up an new "total" water change for substeps
!dir$ concurrent
!cdir nodep
    do fc = 1, num_hydrologyc
      c = filter_hydrologyc(fc)
      do m = 1,pwtrc
         wtr_qflx_drain(c,m) = 0._r8
      end do
      do j = 1, nlevsoi
        h2onew(c,j) = h2osoi_old(c,j)
      end do
    end do

    ! Iterate over substeps

    do isub = 1, nsub

      ! Assign a partial updated water content

      do j = 1, nlevsoi
        do fc = 1, num_hydrologyc
           c = filter_hydrologyc(fc)
           h2onew(c,j) = h2onew(c,j) + dtsub*delh2o(c,j)
        end do
      end do

    ! Set soil ratios. Notice use ixbase for ratios, since incremental.
!dir$ concurrent
!cdir nodep
      do m = 1, pwtrc
        do j = 1, nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             if (j.eq.1) then
               RSoil(c,j,m) = get_wratio(wtr_h2osoi_liq(c,j,m)+dtsub*wtr_qflx_infl(c,m), &
                                         wtr_h2osoi_liq(c,j,ixbase)+dtsub*wtr_qflx_infl(c,ixbase),h2otiny)
             else
               RSoil(c,j,m) = get_wratio(wtr_h2osoi_liq(c,j,m),wtr_h2osoi_liq(c,j,ixbase),h2otiny)
             endif
          end do
        end do
      end do

    ! assign tracer flux ratios (special case for below-model recharge)
    ! (These could be done as interface mean values for second order accuracy)
      maxdt = dtime
!dir$ concurrent
!cdir nodep
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         do m = 1, pwtrc
           RFlx(c,1      ,m) = 0._r8                ! at top
         end do
         if (qflx_soi(c,nlevsoi+1) > 0.) then       ! (down)
            do m = 1, pwtrc
              RFlx(c,nlevsoi+1,m) = qflx_soi(c,nlevsoi+1)*RSoil(c,nlevsoi,m)
              if (RFlx(c,nlevsoi+1,m) > flim*wtr_h2osoi_liq(c,nlevsoi,m)/dtsub ) then
!!                 write(*,*) 'LIMITING RFLUX at base:',RFlx(c,nlevsoi+1,m),  &
!!                                     flim*wtr_h2osoi_liq(c,nlevsoi,m)/dtsub
                  if ( RFlx(c,nlevsoi+1,m) > 0.0_r8 ) &    !EK 25 Jun 2015
                  maxdt = min(maxdt,flim*wtr_h2osoi_liq(c,nlevsoi,m)/RFlx(c,nlevsoi+1,m))
                  RFlx(c,nlevsoi+1,m) = flim*wtr_h2osoi_liq(c,nlevsoi,m)/dtsub
              end if
            end do
         else                                       ! (up)
            do m = 1, pwtrc
              RFlx(c,nlevsoi+1,m) = qflx_soi(c,nlevsoi+1)*RGNIP(c,m)
            end do
         end if
      end do
!
      do j = 2, nlevsoi
!dir$ concurrent
!cdir nodep
        do fc = 1, num_hydrologyc
           c = filter_hydrologyc(fc)
           if (qflx_soi(c,j) > 0.) then             ! (down)
              do m = 1, pwtrc
                RFlx(c,j,m) = qflx_soi(c,j)*RSoil(c,j-1,m)
                if (RFlx(c,j,m) > flim*wtr_h2osoi_liq(c,j-1,m)/dtsub) then
!!                  write(*,*) 'LIMIT DOWN RFLUX at i,j,m=',c,j,m, &
!!                          RFlx(c,j,m),flim*wtr_h2osoi_liq(c,j-1,m)/dtsub
                  if ( RFlx(c,j,m) > 0.0_r8 ) &    !EK 25 Jun 2015
                  maxdt = min(maxdt,flim*wtr_h2osoi_liq(c,j-1,m)/RFlx(c,j,m))
                  RFlx(c,j,m) = flim*wtr_h2osoi_liq(c,j-1,m)/dtsub
                end if
              end do
           else                                     ! (upward)
              do m = 1, pwtrc
                RFlx(c,j,m) = qflx_soi(c,j)*RSoil(c,j  ,m)
                if (RFlx(c,j,m) < -flim*wtr_h2osoi_liq(c,j,m)/dtsub ) then
!!                  write(*,*) 'LIMIT UP RFLUX at i,j,m=',c,j,m, &
!!                        RFlx(c,j,m),flim*wtr_h2osoi_liq(c,j,m)/dtsub
                  if ( RFlx(c,j,m) /= 0.0_r8 ) &    !EK 25 Jun 2015
                  maxdt = min(maxdt,-flim*wtr_h2osoi_liq(c,j,m)/RFlx(c,j,m))
                  RFlx(c,j,m) = -flim*wtr_h2osoi_liq(c,j,m)/dtsub
                end if
              end do
           end if
        end do
      end do

      if (maxdt < dtsub) then
        write(*,*) 'TracerSoilWater: dtsub too big. Maxdt =',maxdt
      end if

      ! Update the profile in the top layer
      ! Solve forward problem for all tracers. 
      ! (Should correct save xylem water if substepping, as per drainage)
!dir$ concurrent
!cdir nodep
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        do m = 1, pwtrc
          dflux = wtr_qflx_infl(c,m) + (RFlx(c,1,m) - RFlx(c,2,m))
          wtr_delh2o(c,1) = dflux - RSoil(c,1,m)*qflx_tran_veg_col(c)*rootr_col(c,1)
          wtr_h2osoi_liq(c,1,m) = wtr_h2osoi_liq(c,1,m) + dtsub*wtr_delh2o(c,1)
          Rnew(m) = get_wratio(wtr_h2osoi_liq(c,1,m),wtr_h2osoi_liq(c,1,ixbase),h2otiny)
        end do
      end do

      ! Continue updating for the simpler soil (no infiltration or equilibration)
      do j = 2, nlevsoi
!dir$ concurrent
!cdir nodep
        do fc = 1, num_hydrologyc
           c = filter_hydrologyc(fc)
           do m = 1, pwtrc
              dflux = RFlx(c,j,m) - RFlx(c,j+1,m) 
              wtr_delh2o(c,j) = dflux - RSoil(c,j,m)*qflx_tran_veg_col(c)*rootr_col(c,j)
              wtr_h2osoi_liq(c,j,m) = wtr_h2osoi_liq(c,j,m) + dtsub*wtr_delh2o(c,j)
           end do
        end do
      end do

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         do m =1,pwtrc
           RSoil(c,nlevsoi,m) = get_wratio(wtr_h2osoi_liq(c,nlevsoi,m),wtr_h2osoi_liq(c,nlevsoi,ixbase),h2otiny)
         end do

         if (qflx_soi(c,nlevsoi+1) > 0.) then       ! (down)
            do m = 1, pwtrc
              RFlx(c,nlevsoi+1,m) = qflx_soi(c,nlevsoi+1)*RSoil(c,nlevsoi,m)
            end do
         else                       ! (up)
            do m = 1, pwtrc
              RFlx(c,nlevsoi+1,m) = qflx_soi(c,nlevsoi+1)*RGNIP(c,m)
            end do
         end if
      end do

!dir$ concurrent
!cdir nodep
      do m = 1, pwtrc
        do fc = 1, num_hydrologyc
           c = filter_hydrologyc(fc)
           wtr_qflx_drain(c,m) = wtr_qflx_drain(c,m) + dtsub*RFlx(c,nlevsoi+1,m)
        end do
      end do

    end do      ! isub, substeps

! Finalize the drainage as the flux deficit at the base.
! (Notice CLM does not conserve mass here, so it is better to just
!  apply a scaling so wtr_qflx_drain remains consistent with qflx_drain)
!
!dir$ concurrent
!cdir nodep
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       do m = 1, pwtrc
         wtr_qflx_drain(c,m) = wtr_qflx_drain(c,m)/dtime
       end do

       ! A first guess on tracer drainage, modified in SoilHydrology > Drainage
       ! after the fixes for negative tracer water.
       do m = 1, pwtrc
         wratio(m) = get_wratio(wtr_h2osoi_liq(c,nlevsoi,m),wtr_h2osoi_liq(c,nlevsoi,ixbase),h2otiny)
         wtr_qflx_drain(c,m) = wratio(m)*qflx_drain(c)
       end do
    end do

    return
  end subroutine HydrologyTracerSoilWater

 
!=======================================================================
  subroutine HydrologyTracerRescale(begc,endc,begp,endp, &
                    num_nolakec,filter_nolakec)
!-----------------------------------------------------------------------
!
! Rescales tracer quantities such that the ratio is preserved and the
! tracer H2O exactly matches the prognostic water. This avoids any
! possible problems arising from numerical imprecision, and slow
! accumulation of small numerical error.
!
! Only need to scale the state variables. These are those which appear
! on the restart file.
!
! It is assumed that this routine in called within the clump scan within
! driver.
!
!-----------------------------------------------------------------------
    use clmtype
    use clm_varpar, only : nlevsoi
    use clm_varcon, only : denh2o, denice
    use subgridAveMod
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    integer, intent(in) :: begc, endc		! column array bounds
    integer, intent(in) :: begp, endp		! pft array bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(endc-begc+1)   ! column filter for non-lake points

!-----------------------------------------------------------------------
   
! Local pointers to clmtype variables: column state
   
    integer , pointer :: pcolumn(:)            !pft's column index
    real(r8), pointer :: wtcol(:)              !weight of pft relative to column
    real(r8), pointer :: dz(:,:)               !layer depth (m)

    real(r8), pointer :: h2ocan(:)             !canopy water (mm H2O)
    real(r8), pointer :: h2ocan_col(:)         !canopy water (mm H2O)
    real(r8), pointer :: h2osno(:)             !snow water (mm H2O)
    real(r8), pointer :: h2osoi_liq(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: h2osoi_ice(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: h2osoi_vol(:,:)       !volumetric soil water [m3/m3]  (nlevsoi)
!
    real(r8), pointer :: wtr_h2ocan(:,:)       !tracer canopy tracer water (mm H2O)
    real(r8), pointer :: wtr_h2ocan_col(:,:)   !tracer canopy tracer water (mm H2O)
    real(r8), pointer :: wtr_h2osno(:,:)       !tracer snow water (mm H2O)
    real(r8), pointer :: wtr_h2osoi_liq(:,:,:) !tracer liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: wtr_h2osoi_ice(:,:,:) !tracer ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: wtr_h2osoi_vol(:,:,:) !tracer volumetric soil water  [m3/m3]  (nlevsoi)
!
    real(r8), pointer :: ptrp(:)         ! pointer to input pft array
    real(r8), pointer :: ptrc(:)         ! pointer to output column array

! Local variables

    integer c, fc, p, j, m

    real(r8) :: errcol(begc:endc)		!error on columns
    real(r8) :: errpft(begp:endp)		!error on pfts
    real(r8) :: ratcol(begc:endc,pwtrc)		!ratio of error on columns
    real(r8) :: ratpft(begp:endp,pwtrc)		!ration on error for pfts

!-----------------------------------------------------------------------

    ! Return if no rescale requested

    if (.not. ldorescale) return

    ! Assign local pointers

    pcolumn           => pft%column
    wtcol             => pft%wtcol
    dz                => cps%dz

    h2ocan            => pws%h2ocan
    h2ocan_col        => pws_a%h2ocan
    h2osno            => cws%h2osno
    h2osoi_liq        => cws%h2osoi_liq
    h2osoi_ice        => cws%h2osoi_ice
    h2osoi_vol        => cws%h2osoi_vol

    wtr_h2ocan        => pws%wtr_h2ocan
    wtr_h2ocan_col    => pws_a%wtr_h2ocan
    wtr_h2osno        => cws%wtr_h2osno
    wtr_h2osoi_liq    => cws%wtr_h2osoi_liq
    wtr_h2osoi_ice    => cws%wtr_h2osoi_ice
    wtr_h2osoi_vol    => cws%wtr_h2osoi_vol


    ! H2OSNO (column)

    do c = begc, endc
       errcol(c) =  wtr_h2osno(c,ixbase) - h2osno(c)
!       if (abs(errcol(c)) > errtol) then
!         write(*,*) 'HydrologyTracerRescale: H2OSNO:', c, col%pfti(c),col%pftf(c), errcol(c)
!         write(*,*) wtr_h2osno(c,ixbase), h2osno(c), errcol(c)
!         call endrun(' error too big to ignore')
!       endif
       do m = 1, pwtrc
         ratcol(c,m) = get_wratio(wtr_h2osno(c,m),wtr_h2osno(c,ixbase),h2otiny)
       end do
       ! Must do in 2 loops, otherwise altering the ixbase
       do m = 1, pwtrc
          wtr_h2osno(c,m) = wtr_h2osno(c,m) - ratcol(c,m)*errcol(c)
       end do
    end do

    ! H2OSOI_LIQ (column)

    do j = 1, nlevsoi
      do c = begc, endc
         errcol(c) =  wtr_h2osoi_liq(c,j,ixbase) - h2osoi_liq(c,j)
!         if (abs(errcol(c)) > errtol) then
!           write(*,*) 'HydrologyTracerRescale: H2OSOI_LIQ:', c,col%pfti(c),col%pftf(c), j, errcol(c)
!           write(*,*) wtr_h2osoi_liq(c,j,ixbase), h2osoi_liq(c,j), errcol(c)
!           call endrun(' error too big to ignore')
!         endif
         do m = 1, pwtrc
           ratcol(c,m) = get_wratio(wtr_h2osoi_liq(c,j,m), wtr_h2osoi_liq(c,j,ixbase),h2otiny)
         end do
         ! Must do in 2 loops, otherwise altering the ixbase
         do m = 1, pwtrc
            wtr_h2osoi_liq(c,j,m) = wtr_h2osoi_liq(c,j,m) - ratcol(c,m)*errcol(c)
         end do
      end do
    end do

    ! H2OSOI_ICE (column)

    do j = 1, nlevsoi
      do c = begc, endc
         errcol(c) =  wtr_h2osoi_ice(c,j,ixbase) - h2osoi_ice(c,j)
!         if (abs(errcol(c)) > errtol) then
!           write(*,*) 'HydrologyTracerRescale: H2OSOI_ICE:', c,col%pfti(c),col%pftf(c), j, errcol(c)
!           write(*,*) wtr_h2osoi_ice(c,j,ixbase), h2osoi_ice(c,j), errcol(c)
!           call endrun(' error too big to ignore')
!         endif
         do m = 1, pwtrc
           ratcol(c,m) = get_wratio(wtr_h2osoi_ice(c,j,m),wtr_h2osoi_ice(c,j,ixbase),h2otiny)
         end do
         do m = 1, pwtrc
            wtr_h2osoi_ice(c,j,m) = wtr_h2osoi_ice(c,j,m) - ratcol(c,m)*errcol(c)
         end do
      end do
    end do

    ! Redo the volumetric soil water where not lake (else should be "spval")

    do m = 1, pwtrc
      do j = 1, nlevsoi
        do fc = 1, num_nolakec
           c = filter_nolakec(fc)
           wtr_h2osoi_vol(c,j,m) = wtr_h2osoi_liq(c,j,m)/(dz(c,j)*denh2o) + &
                                   wtr_h2osoi_ice(c,j,m)/(dz(c,j)*denice)
        end do
      end do
    end do

#if 0
    ! H2OCAN (pft)

    do p = begp, endp
       errpft(p) =  wtr_h2ocan(p,ixbase) - h2ocan(p)
!       if (abs(errpft(p)) > errtol) then
!         write(*,*) 'HydrologyTracerRescale: H2OCAN:', p, pft%column(p), col%pfti(pft%column(p)), col%pftf(pft%column(p)), errpft(p)
!         write(*,*) wtr_h2ocan(p,ixbase), h2ocan(p)
!         call endrun(' error too big to ignore')
!       endif
       do m = 1, pwtrc
         ratpft(p,m) = get_wratio(wtr_h2ocan(p,m),wtr_h2ocan(p,ixbase),h2otiny)
       end do
       do m = 1, pwtrc
          wtr_h2ocan(p,m) = wtr_h2ocan(p,m) - ratpft(p,m)*errpft(p)
       end do
    end do

    ! H2OCAN (col)
    ! just do this instead of separate calculation below
    ! since have corrected h2ocan_pft from above

    call p2c(num_nolakec, filter_nolakec, h2ocan, h2ocan_col)
    call p2c(pwtrc, num_nolakec, filter_nolakec, wtr_h2ocan, wtr_h2ocan_col)
#if 0
    do c = begc, endc
       errcol(c) =  wtr_h2ocan_col(c,ixbase) - h2ocan_col(c)
!       if (abs(errcol(c)) > errtol) then
!         write(*,*) 'HydrologyTracerRescale: H2OCAN_A:', c,col%pfti(c),col%pftf(c), errcol(c)
!         write(*,*) wtr_h2ocan_col(c,ixbase), h2ocan_col(c)
!         call endrun(' error too big to ignore')
!       endif
       do m = 1, pwtrc
         ratcol(p,m) = get_wratio(wtr_h2ocan_col(c,m),wtr_h2ocan_col(c,ixbase),h2otiny)
       end do
       do m = 1, pwtrc
          wtr_h2ocan_col(c,m) = wtr_h2ocan_col(c,m) - ratcol(c,m)*errcol(c)
       end do
    end do
#endif
!!   ! Reassign colum average H2OCAN (just like p2c calls in pft2colMod)
!!
!!    do m = 1, pwtrc
!!      ptrp => clm3%g%l%c%p%pws%wtr_h2ocan(:,m)
!!      ptrc => clm3%g%l%c%cws%pws_a%wtr_h2ocan(:,m)
!!      call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)
!!    end do

#endif
    
    return
  end subroutine HydrologyTracerRescale


!=======================================================================
  subroutine HydrologyTracerCheck(subr,lbg, ubg, lbc, ubc, lbp, ubp)
!-----------------------------------------------------------------------
!
! Checks all tracer variables for consistency with prognostic water.
!
! These check access the global clm type variables for an entire
! "clump".   As such, it is expected that this routine be called from 
! "driver" or immediate decendents.
!
!-----------------------------------------------------------------------
    use perf_mod, only : t_startf, t_stopf
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*),intent(in) :: subr 	! name of calling routine
    integer, intent(in)  :: lbp, ubp    ! clump beginning and ending pft indices
    integer, intent(in)  :: lbc, ubc    ! clump beginning and ending column indices
    integer, intent(in)  :: lbg, ubg    ! clump beginning and ending gridcell indices
!------------------------- Local Variables -----------------------------
    logical lfailany			  ! flag for any failed check
!-----------------------------------------------------------------------
    if (lchecknone) return
    call t_startf('HydTrcCheck')
    if (lcheckverb) write(6,*) 'TracerCheck: (',trim(subr),')'

    ! Reset error reporting flag as inocent until proven guilty
    lfailany = .false.

    ! Do all the cheks
    call t_startf('WtrTrcCheckPFT')
    call TracerCheckPFT   (subr,lbp,ubp,lfailany)		! fluxes mostly
    call t_stopf('WtrTrcCheckPFT')

    call t_startf('WtrTrcCheckCol')
    call TracerCheckColumn(subr,lbc,ubc,lfailany)		! state and fluxes
    call t_stopf('WtrTrcCheckCol')

    call t_startf('WtrTrcCheckGrid')
    call TracerCheckGrid  (subr,lbg,ubg,lfailany)		! input/output
    call t_stopf('WtrTrcCheckGrid')


    ! Same but check delta values  (for HDO and H218O)

    call t_startf('WtrTrcCheckPFTDel')
    call TracerCheckPFTDelta(subr,lbc,ubc,lfailany)		! PFT isotope check
    call t_stopf('WtrTrcCheckPFTDel')

    call t_startf('WtrTrcCheckColDel')
    call TracerCheckColumnDelta(subr,lbc,ubc,lfailany)		! Column isotope check
    call t_stopf('WtrTrcCheckColDel')

    ! Bail out as needed
    if (lfailany .and. lcheckstop) &
       call endrun('TracerCheck: Terminating on failed tracer check.')

    if (lcheckverb) write(6,*) 'TracerCheck: (',trim(subr),') done.'

    call t_stopf('HydTrcCheck')
    return
  end subroutine HydrologyTracerCheck


!=======================================================================
  subroutine TracerCheckGrid(subr,lbg,ubg,lfailany)
!-----------------------------------------------------------------------
! Check forcing variables are exactly the same as the base state
!-----------------------------------------------------------------------
    use clmtype
!    use clm_atmlnd,       only: clm_a2l
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*),intent(in) :: subr 	! name of calling routine
    integer, intent(in)  :: lbg, ubg    ! clump beginning and ending gridcell indices
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout) ::  lfailany		! did not pass check
!------------------------- Local Variables -----------------------------
!
! Local pointers to clmtype variables: ATM2LND
!
    real(r8), pointer :: forc_q(:)          !atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_rain(:)       ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)       ! snow rate [mm/s]
!
    real(r8), pointer :: forc_wtr_q(:,:)    ! tracer atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_wtr_rain(:,:) ! tracer rain rate [mm/s]
    real(r8), pointer :: forc_wtr_snow(:,:) ! tracer snow rate [mm/s]
!
! Local pointers to clmtype variables: LND2ATM
!
    real(r8), pointer :: h2osno(:)               !snow water (mm H2O)
    real(r8), pointer :: qflx_evap_tot(:)        !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg

    real(r8), pointer :: wtr_h2osno(:,:)         !tracer snow water (mm H2O)
    real(r8), pointer :: wtr_qflx_evap_tot(:,:)  !tracer qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
!
! Flag if any of the checks fail
!
    logical lfail
!
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members
    forc_q            => cws%forc_q  !clm4
    forc_rain         => cws%forc_rain
    forc_snow         => cws%forc_snow

    qflx_evap_tot     => pwf_a%qflx_evap_tot
    h2osno            => cws%h2osno

    forc_wtr_q        => cws%forc_wtr_q  !clm4
    forc_wtr_rain     => cws%forc_wtr_rain
    forc_wtr_snow     => cws%forc_wtr_snow

    wtr_qflx_evap_tot => pwf_a%wtr_qflx_evap_tot
    wtr_h2osno        => cws%wtr_h2osno


    ! Set local checking flag
    lfail = .false.

    ! Do the checks
    call Check1d('FORC_Q'    ,lbg,ubg,etols,lfail,forc_q   ,forc_wtr_q   (:,ixtest))
    call Check1d('FORC_RAIN' ,lbg,ubg,etolf,lfail,forc_rain,forc_wtr_rain(:,ixtest))
    call Check1d('FORC_SNOW' ,lbg,ubg,etolf,lfail,forc_snow,forc_wtr_snow(:,ixtest))
    call Check1d('Q_EVAP_TOT',lbg,ubg,etolf,lfail,qflx_evap_tot,wtr_qflx_evap_tot(:,ixtest))
    call Check1d('H2OSNO'    ,lbg,ubg,etols,lfail,h2osno   ,wtr_h2osno   (:,ixtest))

    ! Report
    if (lfail) then
      write(6,*) 'TracerCheckGrid: ('//trim(subr)//')  ***** FAILED *****'
!!      if (lcheckstop) call endrun('ABORTED')
    else
      if (lcheckverb) write(6,*) 'TracerCheckGrid: ('//trim(subr)//')  All OK'
    endif

    ! update universal flag
    lfailany = lfailany .or. lfail
    return
  end subroutine TracerCheckGrid


!=======================================================================
  subroutine TracerCheckColumn(subr,lbc,ubc,lfailany)
!-----------------------------------------------------------------------
! Check column level tracer state is exactly the same as the base state
!-----------------------------------------------------------------------
    use clmtype
    use clm_varpar, only: nlevsoi, nlevsno
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*),intent(in) :: subr 	! name of calling routine
    integer, intent(in) :: lbc, ubc		! begining/end indicies
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout) ::  lfailany		! did not pass check
!------------------------- Local Variables -----------------------------
!
! Local pointers to clmtype variables: column state
!   
    real(r8), pointer :: h2osno(:)             !snow water (mm H2O)
    real(r8), pointer :: h2osoi_liq(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: h2osoi_ice(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: h2osoi_vol(:,:)       !volumetric soil water  [m3/m3]  (nlevsoi)
    real(r8), pointer :: qg(:)                 !ground specific humidity [kg/kg]
    real(r8), pointer :: snowice(:)            !average snow ice lens
    real(r8), pointer :: snowliq(:)            !average snow liquid water
!
    real(r8), pointer :: wtr_h2osno(:,:)       !tracer snow water (mm H2O)
    real(r8), pointer :: wtr_h2osoi_liq(:,:,:) !tracer liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: wtr_h2osoi_ice(:,:,:) !tracer ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: wtr_h2osoi_vol(:,:,:) !tracer volumetric soil water  [m3/m3]  (nlevsoi)
    real(r8), pointer :: wtr_qg(:,:)           !tracer ground specific humidity [kg/kg]
    real(r8), pointer :: wtr_snowice(:,:)      !tracer average snow ice lens
!
! Local pointers to clmtype variables: column flux
!   
    real(r8), pointer :: qflx_infl(:)          !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)          !surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_drain(:)         !sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_top_soil(:)      !net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_snomelt(:)       !snow melt (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)         !qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qmelt(:)              !snow melt [mm/s]
!
    real(r8), pointer :: wtr_snowliq(:,:)      !tracer average snow liquid water
    real(r8), pointer :: wtr_qflx_infl(:,:)    !tracer infiltration (mm H2O /s)
    real(r8), pointer :: wtr_qflx_surf(:,:)    !tracer surface runoff (mm H2O /s)
    real(r8), pointer :: wtr_qflx_drain(:,:)   !tracer sub-surface runoff (mm H2O /s)
    real(r8), pointer :: wtr_qflx_top_soil(:,:)!tracer net water input into soil from top (mm/s)
    real(r8), pointer :: wtr_qflx_snomelt(:,:) !tracer snow melt (mm H2O /s)
    real(r8), pointer :: wtr_qflx_qrgwl(:,:)   !tracer qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: wtr_qmelt(:,:)        !tracer snow melt [mm/s]
!
! PFT level quantities averages up to column
!
    real(r8), pointer :: h2ocan_a(:)           !canopy water (mm H2O)
    real(r8), pointer :: qflx_prec_intr(:) !interception of precipitation [mm/s]
    real(r8), pointer :: qflx_prec_grnd(:) !water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_rain_grnd(:) !rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd(:) !snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_liq(:) !excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_ice(:) !excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_veg(:)  !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:)  !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can(:)  !evaporation from leaves and stems
    real(r8), pointer :: qflx_evap_soi(:)  !soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:)  !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_grnd(:) !ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)  !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow(:)  !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow(:)  !sublimation rate from snow pack (mm H2O /s) [+]
!
    real(r8), pointer :: wtr_h2ocan_a(:,:)       !tracer canopy water (mm H2O)
    real(r8), pointer :: wtr_qflx_prec_intr(:,:) !tracer interception of precipitation [mm/s]
    real(r8), pointer :: wtr_qflx_prec_grnd(:,:) !tracer water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: wtr_qflx_rain_grnd(:,:) !tracer rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snow_grnd(:,:) !tracer snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_liq(:,:)   !tracer excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_ice(:,:)   !tracer excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_evap_veg(:,:)  !tracer vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_tran_veg(:,:)  !tracer vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_can(:,:)  !tracer evaporation from leaves and stems
    real(r8), pointer :: wtr_qflx_evap_soi(:,:)  !tracer soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_tot(:,:)  !tracer qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: wtr_qflx_evap_grnd(:,:) !tracer ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_dew_grnd(:,:)  !tracer ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_dew_snow(:,:)  !tracer surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_sub_snow(:,:)  !tracer sublimation rate from snow pack (mm H2O /s) [+]

! Column level water balance variables
    real(r8), pointer :: endwb(:)           ! water mass end of the time step
    real(r8), pointer :: begwb(:)           ! water mass begining of the time step
    real(r8), pointer :: wtr_endwb(:,:)     ! water mass end of the time step
    real(r8), pointer :: wtr_begwb(:,:)     ! water mass begining of the time step

! Dummies for snow and soil layers
    real(r8) h2osoi_soi(lbc:ubc,nlevsoi),wtr_h2osoi_soi(lbc:ubc,nlevsoi,pwtrc)
!
    integer c,j,m
!
    logical lfail				! local failure flag
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members
    h2osno            => cws%h2osno
    h2osoi_liq        => cws%h2osoi_liq
    h2osoi_ice        => cws%h2osoi_ice
    h2osoi_vol        => cws%h2osoi_vol
    qg                => cws%qg
    snowice           => cws%snowice
    snowliq           => cws%snowliq

    qflx_infl         => cwf%qflx_infl
    qflx_surf         => cwf%qflx_surf
    qflx_drain        => cwf%qflx_drain
    qflx_top_soil     => cwf%qflx_top_soil
    qflx_snomelt      => cwf%qflx_snomelt
    qflx_qrgwl        => cwf%qflx_qrgwl
    qmelt             => cwf%qmelt
!
    wtr_h2osno        => cws%wtr_h2osno
    wtr_h2osoi_liq    => cws%wtr_h2osoi_liq
    wtr_h2osoi_ice    => cws%wtr_h2osoi_ice
    wtr_h2osoi_vol    => cws%wtr_h2osoi_vol
    wtr_qg            => cws%wtr_qg
    wtr_snowice       => cws%wtr_snowice
    wtr_snowliq       => cws%wtr_snowliq

    wtr_qflx_infl     => cwf%wtr_qflx_infl
    wtr_qflx_surf     => cwf%wtr_qflx_surf
    wtr_qflx_drain    => cwf%wtr_qflx_drain
    wtr_qflx_top_soil => cwf%wtr_qflx_top_soil
    wtr_qflx_snomelt  => cwf%wtr_qflx_snomelt
    wtr_qflx_qrgwl    => cwf%wtr_qflx_qrgwl
    wtr_qmelt         => cwf%wtr_qmelt

    ! Quantities averaged from PFT level

    h2ocan_a           => pws_a%h2ocan
    qflx_prec_intr     => pwf_a%qflx_prec_intr
    qflx_prec_grnd     => pwf_a%qflx_prec_grnd
    qflx_rain_grnd     => pwf_a%qflx_rain_grnd
    qflx_snow_grnd     => pwf_a%qflx_snow_grnd
    qflx_snwcp_liq     => pwf_a%qflx_snwcp_liq
    qflx_snwcp_ice     => pwf_a%qflx_snwcp_ice
    qflx_tran_veg      => pwf_a%qflx_tran_veg
    qflx_evap_veg      => pwf_a%qflx_evap_veg
    qflx_evap_soi      => pwf_a%qflx_evap_soi
    qflx_evap_can      => pwf_a%qflx_evap_can
    qflx_evap_tot      => pwf_a%qflx_evap_tot
    qflx_evap_grnd     => pwf_a%qflx_evap_grnd
    qflx_dew_grnd      => pwf_a%qflx_dew_grnd
    qflx_dew_snow      => pwf_a%qflx_dew_snow
    qflx_sub_snow      => pwf_a%qflx_sub_snow

    wtr_h2ocan_a       => pws_a%wtr_h2ocan
    wtr_qflx_prec_intr => pwf_a%wtr_qflx_prec_intr
    wtr_qflx_prec_grnd => pwf_a%wtr_qflx_prec_grnd
    wtr_qflx_rain_grnd => pwf_a%wtr_qflx_rain_grnd
    wtr_qflx_snow_grnd => pwf_a%wtr_qflx_snow_grnd
    wtr_qflx_snwcp_liq => pwf_a%wtr_qflx_snwcp_liq
    wtr_qflx_snwcp_ice => pwf_a%wtr_qflx_snwcp_ice
    wtr_qflx_tran_veg  => pwf_a%wtr_qflx_tran_veg
    wtr_qflx_evap_veg  => pwf_a%wtr_qflx_evap_veg
    wtr_qflx_evap_can  => pwf_a%wtr_qflx_evap_can
    wtr_qflx_evap_tot  => pwf_a%wtr_qflx_evap_tot
    wtr_qflx_evap_soi  => pwf_a%wtr_qflx_evap_soi
    wtr_qflx_evap_grnd => pwf_a%wtr_qflx_evap_grnd
    wtr_qflx_dew_grnd  => pwf_a%wtr_qflx_dew_grnd
    wtr_qflx_dew_snow  => pwf_a%wtr_qflx_dew_snow
    wtr_qflx_sub_snow  => pwf_a%wtr_qflx_sub_snow

    ! Water balance variables

    endwb             => cwbal%endwb
    begwb             => cwbal%begwb
    wtr_endwb         => cwbal%wtr_endwb
    wtr_begwb         => cwbal%wtr_begwb


    ! Set local checking flag
    lfail = .false.

    ! Check the state variables

    call Check2d('H2OSNO_LIQ',lbc,ubc,1,nlevsno, &
                                       etols,lfail,h2osoi_liq(:,:)   , wtr_h2osoi_liq(:,:,ixtest))
    call Check2d('H2OSNO_ICE',lbc,ubc,1,nlevsno, &
                                       etols,lfail,h2osoi_ice(:,:)   , wtr_h2osoi_ice(:,:,ixtest))

    ! (pack into local arrays to enable passing via arguments, rather
    ! than expicit.
    do j = 1, nlevsoi
      do c = lbc, ubc
        h2osoi_soi(c,j) = h2osoi_liq(c,j)
        do m = 1, pwtrc
          wtr_h2osoi_soi(c,j,m) = wtr_h2osoi_liq(c,j,m)
        end do
      end do
    end do
    call Check2d('H2OSOI_LIQ',lbc,ubc,1,nlevsoi, &
                                  etols,lfail,h2osoi_soi(:,:), wtr_h2osoi_soi(:,:,ixtest))
    do j = 1, nlevsoi
      do c = lbc, ubc
        h2osoi_soi(c,j) = h2osoi_ice(c,j)
        do m = 1, pwtrc
          wtr_h2osoi_soi(c,j,m) = wtr_h2osoi_ice(c,j,m)
        end do
      end do
    end do
    call Check2d('H2OSOI_ICE',lbc,ubc,1,nlevsoi, &
                                  etols,lfail,h2osoi_soi(:,:), wtr_h2osoi_soi(:,:,ixtest))
    call Check2d('H2OSOI_VOL',lbc,ubc,1,nlevsoi, &
                                       etols,lfail,h2osoi_vol   , wtr_h2osoi_vol(:,:,ixtest))

    call Check1d('QG    '   ,lbc,ubc,etolr,lfail,qg           , wtr_qg           (:,ixtest))
    call Check1d('H2OSNO'   ,lbc,ubc,etols,lfail,h2osno       , wtr_h2osno       (:,ixtest))
    call Check1d('SNOWICE'  ,lbc,ubc,etols,lfail,snowice      , wtr_snowice      (:,ixtest))
    call Check1d('SNOWLIQ'  ,lbc,ubc,etols,lfail,snowliq      , wtr_snowliq      (:,ixtest))

    ! Check the flux variabes
    call Check1d('Q_INFL'   ,lbc,ubc,etolf,lfail,qflx_infl    , wtr_qflx_infl    (:,ixtest))
    call Check1d('Q_SURF'   ,lbc,ubc,etolf,lfail,qflx_surf    , wtr_qflx_surf    (:,ixtest))
    call Check1d('Q_DRAIN'  ,lbc,ubc,etolf,lfail,qflx_drain   , wtr_qflx_drain   (:,ixtest))
    call Check1d('Q_TOPSOIL',lbc,ubc,etolf,lfail,qflx_top_soil, wtr_qflx_top_soil(:,ixtest))
    call Check1d('Q_SNOMELT',lbc,ubc,etolf,lfail,qflx_snomelt , wtr_qflx_snomelt (:,ixtest))
    call Check1d('Q_QRGWL'  ,lbc,ubc,etolf,lfail,qflx_qrgwl   , wtr_qflx_qrgwl   (:,ixtest))
    call Check1d('QMELT'    ,lbc,ubc,etolf,lfail,qmelt        , wtr_qmelt        (:,ixtest))


    ! Check PFT level fluxes averaged to column

    call Check1d('H2OCAN_A' ,lbc,ubc,etols,lfail,h2ocan_a      , wtr_h2ocan_a      (:,ixtest))

    call Check1d('QPINTR_A' ,lbc,ubc,etolf,lfail,qflx_prec_intr, wtr_qflx_prec_intr(:,ixtest))
    call Check1d('QPGRND_A' ,lbc,ubc,etolf,lfail,qflx_prec_grnd, wtr_qflx_prec_grnd(:,ixtest))
    call Check1d('QRINTR_A' ,lbc,ubc,etolf,lfail,qflx_rain_grnd, wtr_qflx_rain_grnd(:,ixtest))
    call Check1d('QRSNOW_A' ,lbc,ubc,etolf,lfail,qflx_snow_grnd, wtr_qflx_snow_grnd(:,ixtest))
    call Check1d('QSNCAPLIQ_A' ,lbc,ubc,etolf,lfail,qflx_snwcp_liq, wtr_qflx_snwcp_liq(:,ixtest))
    call Check1d('QSNCAPICE_A' ,lbc,ubc,etolf,lfail,qflx_snwcp_ice, wtr_qflx_snwcp_ice(:,ixtest))
    call Check1d('QTRVEG_A' ,lbc,ubc,etolf,lfail,qflx_tran_veg , wtr_qflx_tran_veg (:,ixtest))
    call Check1d('QEVVEG_A' ,lbc,ubc,etolf,lfail,qflx_evap_veg , wtr_qflx_evap_veg (:,ixtest))
    call Check1d('QEVCAN_A' ,lbc,ubc,etolf,lfail,qflx_evap_can , wtr_qflx_evap_can (:,ixtest))
    call Check1d('QEVSOI_A' ,lbc,ubc,etolf,lfail,qflx_evap_soi , wtr_qflx_evap_soi (:,ixtest))
    call Check1d('QEVTOT_A' ,lbc,ubc,etolf,lfail,qflx_evap_tot , wtr_qflx_evap_tot (:,ixtest))
    call Check1d('QEVGND_A' ,lbc,ubc,etolf,lfail,qflx_evap_grnd, wtr_qflx_evap_grnd(:,ixtest))
    call Check1d('QDWGND_A' ,lbc,ubc,etolf,lfail,qflx_dew_grnd , wtr_qflx_dew_grnd (:,ixtest))
    call Check1d('QSBSNO_A' ,lbc,ubc,etolf,lfail,qflx_sub_snow , wtr_qflx_sub_snow (:,ixtest))
    call Check1d('QDWSNO_A' ,lbc,ubc,etolf,lfail,qflx_dew_snow , wtr_qflx_dew_snow (:,ixtest))

    ! Check water balance variables

    call Check1d('BEGWB'    ,lbc,ubc,etols,lfail,begwb         , wtr_begwb(:,ixtest))
    call Check1d('ENDWB'    ,lbc,ubc,etols,lfail,endwb         , wtr_endwb(:,ixtest))

    ! Report
    if (lfail) then
      write(6,*) 'TracerCheckColumn: ('//trim(subr)//')  ***** FAILED *****'
!!      if (lcheckstop) call endrun('ABORTED')
    else
      if (lcheckverb) write(6,*) 'TracerCheckColumn: ('//trim(subr)//')  All OK'
    endif

    ! update universal flag
    lfailany = lfailany .or. lfail
    return
  end subroutine TracerCheckColumn


!=======================================================================
  subroutine TracerCheckPFT(subr,lbp,ubp,lfailany)
!-----------------------------------------------------------------------
! Check PFT level tracer state is exactly the same as the base state
!-----------------------------------------------------------------------
    use clmtype
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*),intent(in) :: subr 	! name of calling routine
    integer, intent(in) :: lbp, ubp		! begining/end indicies
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout) ::  lfailany		! did not pass check
!------------------------- Local Variables -----------------------------
!
! Local pointers to clmtype variables: PFT state
!   
    real(r8), pointer :: h2ocan(:)         !canopy water (mm H2O)
    real(r8), pointer :: wtr_h2ocan(:,:)   !tracer canopy tracer water (mm H2O
!
! Local pointers to clmtype variables: PFT flux
!
    real(r8), pointer :: qflx_prec_intr(:) !interception of precipitation [mm/s]
    real(r8), pointer :: qflx_prec_grnd(:) !water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_rain_grnd(:) !rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_liq(:) !excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_ice(:) !excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snow_grnd(:) !snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_evap_veg(:)  !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:)  !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can(:)  !evaporation from leaves and stems
    real(r8), pointer :: qflx_evap_soi(:)  !soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:)  !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_grnd(:) !ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)  !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow(:)  !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow(:)  !sublimation rate from snow pack (mm H2O /s) [+]
!
    real(r8), pointer :: wtr_qflx_prec_intr(:,:) !tracer interception of precipitation [mm/s]
    real(r8), pointer :: wtr_qflx_prec_grnd(:,:) !tracer water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: wtr_qflx_rain_grnd(:,:) !tracer rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snow_grnd(:,:) !tracer snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_liq(:,:) !tracer excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_ice(:,:) !tracer excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_evap_veg(:,:)  !tracer vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_tran_veg(:,:)  !tracer vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_can(:,:)  !tracer evaporation from leaves and stems
    real(r8), pointer :: wtr_qflx_evap_soi(:,:)  !tracer soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_tot(:,:)  !tracer qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: wtr_qflx_evap_grnd(:,:) !tracer ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_dew_grnd(:,:)  !tracer ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_dew_snow(:,:)  !tracer surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_sub_snow(:,:)  !tracer sublimation rate from snow pack (mm H2O /s) [+]
!
    logical lfail				! local failure flag
!-----------------------------------------------------------------------

    ! Assign local pointers to derive type arrays

    h2ocan             => pws%h2ocan
    qflx_prec_intr     => pwf%qflx_prec_intr
    qflx_prec_grnd     => pwf%qflx_prec_grnd
    qflx_rain_grnd     => pwf%qflx_rain_grnd
    qflx_snow_grnd     => pwf%qflx_snow_grnd
    qflx_snwcp_liq     => pwf%qflx_snwcp_liq
    qflx_snwcp_ice     => pwf%qflx_snwcp_ice
    qflx_tran_veg      => pwf%qflx_tran_veg
    qflx_evap_veg      => pwf%qflx_evap_veg
    qflx_evap_soi      => pwf%qflx_evap_soi
    qflx_evap_can      => pwf%qflx_evap_can
    qflx_evap_tot      => pwf%qflx_evap_tot
    qflx_evap_grnd     => pwf%qflx_evap_grnd
    qflx_dew_grnd      => pwf%qflx_dew_grnd
    qflx_dew_snow      => pwf%qflx_dew_snow
    qflx_sub_snow      => pwf%qflx_sub_snow

    wtr_h2ocan         => pws%wtr_h2ocan
    wtr_qflx_prec_intr => pwf%wtr_qflx_prec_intr
    wtr_qflx_prec_grnd => pwf%wtr_qflx_prec_grnd
    wtr_qflx_rain_grnd => pwf%wtr_qflx_rain_grnd
    wtr_qflx_snow_grnd => pwf%wtr_qflx_snow_grnd
    wtr_qflx_snwcp_liq => pwf%wtr_qflx_snwcp_liq
    wtr_qflx_snwcp_ice => pwf%wtr_qflx_snwcp_ice
    wtr_qflx_tran_veg  => pwf%wtr_qflx_tran_veg
    wtr_qflx_evap_veg  => pwf%wtr_qflx_evap_veg
    wtr_qflx_evap_can  => pwf%wtr_qflx_evap_can
    wtr_qflx_evap_tot  => pwf%wtr_qflx_evap_tot
    wtr_qflx_evap_soi  => pwf%wtr_qflx_evap_soi
    wtr_qflx_evap_grnd => pwf%wtr_qflx_evap_grnd
    wtr_qflx_dew_grnd  => pwf%wtr_qflx_dew_grnd
    wtr_qflx_dew_snow  => pwf%wtr_qflx_dew_snow
    wtr_qflx_sub_snow  => pwf%wtr_qflx_sub_snow

    ! Set local checking flag
    lfail = .false.

    ! Check flux variables
    call Check1d('QPINTR',lbp,ubp,etolf,lfail,qflx_prec_intr, wtr_qflx_prec_intr(:,ixtest))
    call Check1d('QPGRND',lbp,ubp,etolf,lfail,qflx_prec_grnd, wtr_qflx_prec_grnd(:,ixtest))
    call Check1d('QRINTR',lbp,ubp,etolf,lfail,qflx_rain_grnd, wtr_qflx_rain_grnd(:,ixtest))
    call Check1d('QRSNOW',lbp,ubp,etolf,lfail,qflx_snow_grnd, wtr_qflx_snow_grnd(:,ixtest))
    call Check1d('QSNCAPLIQ',lbp,ubp,etolf,lfail,qflx_snwcp_liq, wtr_qflx_snwcp_liq(:,ixtest))
    call Check1d('QSNCAPICE',lbp,ubp,etolf,lfail,qflx_snwcp_ice, wtr_qflx_snwcp_ice(:,ixtest))
    call Check1d('QTRVEG',lbp,ubp,etolf,lfail,qflx_tran_veg , wtr_qflx_tran_veg (:,ixtest))
    call Check1d('QEVVEG',lbp,ubp,etolf,lfail,qflx_evap_veg , wtr_qflx_evap_veg (:,ixtest))
    call Check1d('QEVCAN',lbp,ubp,etolf,lfail,qflx_evap_can , wtr_qflx_evap_can (:,ixtest))
    call Check1d('QEVSOI',lbp,ubp,etolf,lfail,qflx_evap_soi , wtr_qflx_evap_soi (:,ixtest))
    call Check1d('QEVTOT',lbp,ubp,etolf,lfail,qflx_evap_tot , wtr_qflx_evap_tot (:,ixtest))
    call Check1d('QEVGND',lbp,ubp,etolf,lfail,qflx_evap_grnd, wtr_qflx_evap_grnd(:,ixtest))
    call Check1d('QDWGND',lbp,ubp,etolf,lfail,qflx_dew_grnd , wtr_qflx_dew_grnd (:,ixtest))
    call Check1d('QSBSNO',lbp,ubp,etolf,lfail,qflx_sub_snow , wtr_qflx_sub_snow (:,ixtest))
    call Check1d('QDWSNO',lbp,ubp,etolf,lfail,qflx_dew_snow , wtr_qflx_dew_snow (:,ixtest))

    ! Check state variables
    call Check1d('H2OCAN',lbp,ubp,etols,lfail,h2ocan        , wtr_h2ocan(:,ixtest))

    ! Report
    if (lfail) then
      write(6,*) 'TracerCheckPFT: ('//trim(subr)//')  ***** FAILED *****'
!!      if (lcheckstop) call endrun('ABORTED')
    else
      if (lcheckverb) write(6,*) 'TracerCheckPFT: ('//trim(subr)//')  All OK'
    endif

    ! update universal flag
    lfailany = lfailany .or. lfail

    return
  end subroutine TracerCheckPFT


!=======================================================================
  subroutine TracerCheckColumnDelta(subr,lbc,ubc,lfailany)
!-----------------------------------------------------------------------
! Checks delta values of tracers fall with reasonable bounds.
!-----------------------------------------------------------------------
    use clmtype
    use clm_varpar, only: nlevsoi, nlevsno
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*),intent(in) :: subr 	! name of calling routine
    integer, intent(in) :: lbc, ubc		! begining/end indicies
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout) ::  lfailany		! did not pass check
!------------------------- Local Variables -----------------------------
!
! Local pointers to clmtype variables: column state
!   
    real(r8), pointer :: h2osno(:)             !snow water (mm H2O)
    real(r8), pointer :: h2osoi_liq(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: h2osoi_ice(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: h2osoi_vol(:,:)       !volumetric soil water  [m3/m3]  (nlevsoi)
    real(r8), pointer :: qg(:)                 !ground specific humidity [kg/kg]
    real(r8), pointer :: snowice(:)            !average snow ice lens
    real(r8), pointer :: snowliq(:)            !average snow liquid water
!
    real(r8), pointer :: wtr_h2osno(:,:)       !tracer snow water (mm H2O)
    real(r8), pointer :: wtr_h2osoi_liq(:,:,:) !tracer liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: wtr_h2osoi_ice(:,:,:) !tracer ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
    real(r8), pointer :: wtr_h2osoi_vol(:,:,:) !tracer volumetric soil water  [m3/m3]  (nlevsoi)
    real(r8), pointer :: wtr_qg(:,:)           !tracer ground specific humidity [kg/kg]
    real(r8), pointer :: wtr_snowice(:,:)      !tracer average snow ice lens
!
! Local pointers to clmtype variables: column flux
!   
    real(r8), pointer :: qflx_infl(:)          !infiltration (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)          !surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_drain(:)         !sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_top_soil(:)      !net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_snomelt(:)       !snow melt (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)         !qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qmelt(:)              !snow melt [mm/s]
!
    real(r8), pointer :: wtr_snowliq(:,:)      !tracer average snow liquid water
    real(r8), pointer :: wtr_qflx_infl(:,:)    !tracer infiltration (mm H2O /s)
    real(r8), pointer :: wtr_qflx_surf(:,:)    !tracer surface runoff (mm H2O /s)
    real(r8), pointer :: wtr_qflx_drain(:,:)   !tracer sub-surface runoff (mm H2O /s)
    real(r8), pointer :: wtr_qflx_top_soil(:,:)!tracer net water input into soil from top (mm/s)
    real(r8), pointer :: wtr_qflx_snomelt(:,:) !tracer snow melt (mm H2O /s)
    real(r8), pointer :: wtr_qflx_qrgwl(:,:)   !tracer qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: wtr_qmelt(:,:)        !tracer snow melt [mm/s]
!
! PFT level quantities averages up to column
!
    real(r8), pointer :: h2ocan_a(:)           !canopy water (mm H2O)
    real(r8), pointer :: qflx_prec_intr(:) !interception of precipitation [mm/s]
    real(r8), pointer :: qflx_prec_grnd(:) !water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_rain_grnd(:) !rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd(:) !snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_liq(:) !excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_ice(:) !excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_veg(:)  !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:)  !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can(:)  !evaporation from leaves and stems
    real(r8), pointer :: qflx_evap_soi(:)  !soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:)  !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_grnd(:) !ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)  !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow(:)  !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow(:)  !sublimation rate from snow pack (mm H2O /s) [+]
!
    real(r8), pointer :: wtr_h2ocan_a(:,:)       !tracer canopy water (mm H2O)
    real(r8), pointer :: wtr_qflx_prec_intr(:,:) !tracer interception of precipitation [mm/s]
    real(r8), pointer :: wtr_qflx_prec_grnd(:,:) !tracer water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: wtr_qflx_rain_grnd(:,:) !tracer rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snow_grnd(:,:) !tracer snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_liq(:,:)   !tracer excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_ice(:,:)   !tracer excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_evap_veg(:,:)  !tracer vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_tran_veg(:,:)  !tracer vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_can(:,:)  !tracer evaporation from leaves and stems
    real(r8), pointer :: wtr_qflx_evap_soi(:,:)  !tracer soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_tot(:,:)  !tracer qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: wtr_qflx_evap_grnd(:,:) !tracer ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_dew_grnd(:,:)  !tracer ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_dew_snow(:,:)  !tracer surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_sub_snow(:,:)  !tracer sublimation rate from snow pack (mm H2O /s) [+]

    real(r8), pointer :: RleafWaterSun(:,:) ! water at evaporating surface
    real(r8), pointer :: RleafWaterSha(:,:) ! water at evaporating surface
    real(r8), pointer :: RCanopyVapor (:,:) ! vapour in canopy air space 

! Column level water balance variables
    real(r8), pointer :: endwb(:)           ! water mass end of the time step
    real(r8), pointer :: begwb(:)           ! water mass begining of the time step
    real(r8), pointer :: wtr_endwb(:,:)     ! water mass end of the time step
    real(r8), pointer :: wtr_begwb(:,:)     ! water mass begining of the time step

! Dummies for snow and soil layers
    real(r8) h2osoi_soi(lbc:ubc,nlevsoi),wtr_h2osoi_soi(lbc:ubc,nlevsoi,pwtrc)
!
    integer c,j,m


    real(r8) one1d(lbc:ubc) 		    ! and array of ones
!
    logical lfail				! local failure flag
!-----------------------------------------------------------------------

    write(*,*) 'CheckColumnDelta: ('//trim(subr)//')'

    ! Assign local pointers to derived type members
    h2osno            => cws%h2osno
    h2osoi_liq        => cws%h2osoi_liq
    h2osoi_ice        => cws%h2osoi_ice
    h2osoi_vol        => cws%h2osoi_vol
    qg                => cws%qg
    snowice           => cws%snowice
    snowliq           => cws%snowliq

    qflx_infl         => cwf%qflx_infl
    qflx_surf         => cwf%qflx_surf
    qflx_drain        => cwf%qflx_drain
    qflx_top_soil     => cwf%qflx_top_soil
    qflx_snomelt      => cwf%qflx_snomelt
    qflx_qrgwl        => cwf%qflx_qrgwl
    qmelt             => cwf%qmelt
!
    wtr_h2osno        => cws%wtr_h2osno
    wtr_h2osoi_liq    => cws%wtr_h2osoi_liq
    wtr_h2osoi_ice    => cws%wtr_h2osoi_ice
    wtr_h2osoi_vol    => cws%wtr_h2osoi_vol
    wtr_qg            => cws%wtr_qg
    wtr_snowice       => cws%wtr_snowice
    wtr_snowliq       => cws%wtr_snowliq

    wtr_qflx_infl     => cwf%wtr_qflx_infl
    wtr_qflx_surf     => cwf%wtr_qflx_surf
    wtr_qflx_drain    => cwf%wtr_qflx_drain
    wtr_qflx_top_soil => cwf%wtr_qflx_top_soil
    wtr_qflx_snomelt  => cwf%wtr_qflx_snomelt
    wtr_qflx_qrgwl    => cwf%wtr_qflx_qrgwl
    wtr_qmelt         => cwf%wtr_qmelt

    ! Quantities averaged from PFT level

    h2ocan_a           => pws_a%h2ocan
    qflx_prec_intr     => pwf_a%qflx_prec_intr
    qflx_prec_grnd     => pwf_a%qflx_prec_grnd
    qflx_rain_grnd     => pwf_a%qflx_rain_grnd
    qflx_snow_grnd     => pwf_a%qflx_snow_grnd
    qflx_snwcp_liq     => pwf_a%qflx_snwcp_liq
    qflx_snwcp_ice     => pwf_a%qflx_snwcp_ice
    qflx_tran_veg      => pwf_a%qflx_tran_veg
    qflx_evap_veg      => pwf_a%qflx_evap_veg
    qflx_evap_soi      => pwf_a%qflx_evap_soi
    qflx_evap_can      => pwf_a%qflx_evap_can
    qflx_evap_tot      => pwf_a%qflx_evap_tot
    qflx_evap_grnd     => pwf_a%qflx_evap_grnd
    qflx_dew_grnd      => pwf_a%qflx_dew_grnd
    qflx_dew_snow      => pwf_a%qflx_dew_snow
    qflx_sub_snow      => pwf_a%qflx_sub_snow

    wtr_h2ocan_a       => pws_a%wtr_h2ocan
    wtr_qflx_prec_intr => pwf_a%wtr_qflx_prec_intr
    wtr_qflx_prec_grnd => pwf_a%wtr_qflx_prec_grnd
    wtr_qflx_rain_grnd => pwf_a%wtr_qflx_rain_grnd
    wtr_qflx_snow_grnd => pwf_a%wtr_qflx_snow_grnd
    wtr_qflx_snwcp_liq => pwf_a%wtr_qflx_snwcp_liq
    wtr_qflx_snwcp_ice => pwf_a%wtr_qflx_snwcp_ice
    wtr_qflx_tran_veg  => pwf_a%wtr_qflx_tran_veg
    wtr_qflx_evap_veg  => pwf_a%wtr_qflx_evap_veg
    wtr_qflx_evap_can  => pwf_a%wtr_qflx_evap_can
    wtr_qflx_evap_tot  => pwf_a%wtr_qflx_evap_tot
    wtr_qflx_evap_soi  => pwf_a%wtr_qflx_evap_soi
    wtr_qflx_evap_grnd => pwf_a%wtr_qflx_evap_grnd
    wtr_qflx_dew_grnd  => pwf_a%wtr_qflx_dew_grnd
    wtr_qflx_dew_snow  => pwf_a%wtr_qflx_dew_snow
    wtr_qflx_sub_snow  => pwf_a%wtr_qflx_sub_snow

    ! Water balance variables

    endwb             => cwbal%endwb
    begwb             => cwbal%begwb
    wtr_endwb         => cwbal%wtr_endwb
    wtr_begwb         => cwbal%wtr_begwb

    ! Tracer only variables

    RCanopyVapor      => pws_a%RCanopyVapor
    RLeafWaterSun     => pws_a%RleafWaterSun
    RLeafWaterSha     => pws_a%RleafWaterSha

    ! set up ones
    one1d(:)   = 1._r8


    ! Set local checking flag
    lfail = .false.

    ! Check the state variables
    call Check2ddelta('H2OSNO_LIQ',cdbg,lbc,ubc,1,nlevsno, &
                                       etols,lfail,h2osoi_liq   , wtr_h2osoi_liq(:,:,:))
    call Check2ddelta('H2OSNO_ICE',cdbg,lbc,ubc,1,nlevsno, &
                                       etols,lfail,h2osoi_ice   , wtr_h2osoi_ice(:,:,:))

    ! Pack into local arrays to allow passing via arguments (not implicit)
    do j = 1, nlevsoi
      do c = lbc, ubc
        h2osoi_soi(c,j) = h2osoi_liq(c,j)
        do m = 1, pwtrc
          wtr_h2osoi_soi(c,j,m) = wtr_h2osoi_liq(c,j,m)
        end do
      end do
    end do
    call Check2ddelta('H2OSOI_LIQ',cdbg,lbc,ubc,1,nlevsoi, &
                                       etols,lfail,h2osoi_soi   , wtr_h2osoi_soi(:,:,:))
    do j = 1, nlevsoi
      do c = lbc, ubc
        h2osoi_soi(c,j) = h2osoi_ice(c,j)
        do m = 1, pwtrc
          wtr_h2osoi_soi(c,j,m) = wtr_h2osoi_ice(c,j,m)
        end do
      end do
    end do
    call Check2ddelta('H2OSOI_ICE',cdbg,lbc,ubc,1,nlevsoi, &
                                       etols,lfail,h2osoi_soi   , wtr_h2osoi_soi(:,:,:))

    call Check2ddelta('H2OSOI_VOL',cdbg,lbc,ubc,1,         nlevsoi, &
                                       etols,lfail,h2osoi_vol   , wtr_h2osoi_vol(:,:,:))

    call Check1ddelta('QG    '   ,cdbg,lbc,ubc,etolr,lfail,qg           , wtr_qg             (:,:),.false.)
    call Check1ddelta('H2OSNO'   ,cdbg,lbc,ubc,etols,lfail,h2osno       , wtr_h2osno         (:,:),.false.)
    call Check1ddelta('SNOWICE'  ,cdbg,lbc,ubc,etols,lfail,snowice      , wtr_snowice        (:,:),.false.)
    call Check1ddelta('SNOWLIQ'  ,cdbg,lbc,ubc,etols,lfail,snowliq      , wtr_snowliq        (:,:),.false.)


    ! Check the flux variabes
    call Check1ddelta('Q_INFL'   ,cdbg,lbc,ubc,etolf,lfail,qflx_infl    , wtr_qflx_infl      (:,:),.true.)
    call Check1ddelta('Q_SURF'   ,cdbg,lbc,ubc,etolf,lfail,qflx_surf    , wtr_qflx_surf      (:,:),.false.)
    call Check1ddelta('Q_DRAIN'  ,cdbg,lbc,ubc,etolf,lfail,qflx_drain   , wtr_qflx_drain     (:,:),.false.)
    call Check1ddelta('Q_TOPSOIL',cdbg,lbc,ubc,etolf,lfail,qflx_top_soil, wtr_qflx_top_soil  (:,:),.false.)
    call Check1ddelta('Q_SNOMELT',cdbg,lbc,ubc,etolf,lfail,qflx_snomelt , wtr_qflx_snomelt   (:,:),.false.)
    call Check1ddelta('Q_QRGWL'  ,cdbg,lbc,ubc,etolf,lfail,qflx_qrgwl   , wtr_qflx_qrgwl     (:,:),.false.)
    call Check1ddelta('QMELT'    ,cdbg,lbc,ubc,etolf,lfail,qmelt        , wtr_qmelt          (:,:),.false.)


    ! Check PFT level fluxes averaged to column

    call Check1ddelta('H2OCAN_A' ,cdbg,lbc,ubc,etols,lfail,h2ocan_a      , wtr_h2ocan_a      (:,:),.false.)

    ! Check tracer only state variables

    call Check1ddelta('RCanVap_a' ,cdbg,lbc,ubc,etols,lfail,one1d        , RCanopyVapor      (:,:),.false.)
    call Check1ddelta('RLeafSun_a',cdbg,lbc,ubc,etols,lfail,one1d        , RLeafWaterSun     (:,:),.false.)
    call Check1ddelta('RLeafSha_a',cdbg,lbc,ubc,etols,lfail,one1d        , RLeafWaterSha     (:,:),.false.)

    ! PFT level fluxes

    call Check1ddelta('QPINTR_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_prec_intr, wtr_qflx_prec_intr(:,:),.false.)
    call Check1ddelta('QPGRND_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_prec_grnd, wtr_qflx_prec_grnd(:,:),.false.)
    call Check1ddelta('QRINTR_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_rain_grnd, wtr_qflx_rain_grnd(:,:),.false.)
    call Check1ddelta('QRSNOW_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_snow_grnd, wtr_qflx_snow_grnd(:,:),.false.)
    call Check1ddelta('QSNCAPLIQ_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_snwcp_liq, wtr_qflx_snwcp_liq(:,:),.false.)
    call Check1ddelta('QSNCAPICE_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_snwcp_ice, wtr_qflx_snwcp_ice(:,:),.false.)
    call Check1ddelta('QTRVEG_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_tran_veg , wtr_qflx_tran_veg (:,:),.true.)
    call Check1ddelta('QEVVEG_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_evap_veg , wtr_qflx_evap_veg (:,:),.false.)
    call Check1ddelta('QEVCAN_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_evap_can , wtr_qflx_evap_can (:,:),.false.)
    call Check1ddelta('QEVSOI_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_evap_soi , wtr_qflx_evap_soi (:,:),.false.)
    call Check1ddelta('QEVTOT_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_evap_tot , wtr_qflx_evap_tot (:,:),.true.)
    call Check1ddelta('QEVGND_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_evap_grnd, wtr_qflx_evap_grnd(:,:),.true.)
    call Check1ddelta('QDWGND_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_dew_grnd , wtr_qflx_dew_grnd (:,:),.false.)
    call Check1ddelta('QSBSNO_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_sub_snow , wtr_qflx_sub_snow (:,:),.false.)
    call Check1ddelta('QDWSNO_A' ,cdbg,lbc,ubc,etolf,lfail,qflx_dew_snow , wtr_qflx_dew_snow (:,:),.false.)

    ! Check water balance variables

    call Check1ddelta('BEGWB'    ,cdbg,lbc,ubc,etols,lfail,begwb         , wtr_begwb(:,:),.false.)
    call Check1ddelta('ENDWB'    ,cdbg,lbc,ubc,etols,lfail,endwb         , wtr_endwb(:,:),.false.)

    ! Report
    if (lfail) then
      write(6,*) 'TracerCheckColumnDelta: ('//trim(subr)//')  ***** FAILED *****'
      if (lcheckstop) call endrun('ABORTED')
    else
!!      if (lcheckverb) write(6,*) 'TracerCheckColumnDelta: ('//trim(subr)//')  All OK'
      write(6,*) 'TracerCheckColumnDelta: ('//trim(subr)//')  All OK'
    endif

    ! update universal flag
    lfailany = lfailany .or. lfail
    return
  end subroutine TracerCheckColumnDelta


!=======================================================================
  subroutine TracerCheckPFTDelta(subr,lbp,ubp,lfailany)
!-----------------------------------------------------------------------
! Check PFT level tracer state is exactly the same as the base state
!-----------------------------------------------------------------------
    use clmtype
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*),intent(in) :: subr 	! name of calling routine
    integer, intent(in) :: lbp, ubp		! begining/end indicies
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout) ::  lfailany		! did not pass check
!------------------------- Local Variables -----------------------------
!
! Local pointers to clmtype variables: PFT state
!   
    real(r8), pointer :: h2ocan(:)         !canopy water (mm H2O)
    real(r8), pointer :: wtr_h2ocan(:,:)   !tracer canopy tracer water (mm H2O
!
! Local pointers to clmtype variables: PFT flux
!
    real(r8), pointer :: qflx_prec_intr(:) !interception of precipitation [mm/s]
    real(r8), pointer :: qflx_prec_grnd(:) !water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_rain_grnd(:) !rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd(:) !snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_liq(:) !excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_ice(:) !excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_veg(:)  !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:)  !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can(:)  !evaporation from leaves and stems
    real(r8), pointer :: qflx_evap_soi(:)  !soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:)  !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_grnd(:) !ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)  !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow(:)  !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow(:)  !sublimation rate from snow pack (mm H2O /s) [+]
!
    real(r8), pointer :: wtr_qflx_prec_intr(:,:) !tracer interception of precipitation [mm/s]
    real(r8), pointer :: wtr_qflx_prec_grnd(:,:) !tracer water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: wtr_qflx_rain_grnd(:,:) !tracer rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snow_grnd(:,:) !tracer snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_liq(:,:) !tracer excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_snwcp_ice(:,:) !tracer excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_evap_veg(:,:)  !tracer vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_tran_veg(:,:)  !tracer vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_can(:,:)  !tracer evaporation from leaves and stems
    real(r8), pointer :: wtr_qflx_evap_soi(:,:)  !tracer soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtr_qflx_evap_tot(:,:)  !tracer qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: wtr_qflx_evap_grnd(:,:) !tracer ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: wtr_qflx_dew_grnd(:,:)  !tracer ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_dew_snow(:,:)  !tracer surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: wtr_qflx_sub_snow(:,:)  !tracer sublimation rate from snow pack (mm H2O /s) [+]
!
    real(r8), pointer :: RleafWaterSun(:,:)      ! water at evaporating surface
    real(r8), pointer :: RleafWaterSha(:,:)      ! water at evaporating surface
    real(r8), pointer :: RCanopyVapor (:,:)      ! vapour in canopy air space 
!
    real(r8) one1d(lbp:ubp) 		         ! and array of ones
!
    logical lfail				 ! local failure flag
!-----------------------------------------------------------------------

    write(*,*) 'CheckPFTDelta: ('//trim(subr)//')'

    ! Assign local pointers to derive type arrays

    h2ocan             => pws%h2ocan
    qflx_prec_intr     => pwf%qflx_prec_intr
    qflx_prec_grnd     => pwf%qflx_prec_grnd
    qflx_rain_grnd     => pwf%qflx_rain_grnd
    qflx_snow_grnd     => pwf%qflx_snow_grnd
    qflx_snwcp_liq     => pwf%qflx_snwcp_liq
    qflx_snwcp_ice     => pwf%qflx_snwcp_ice
    qflx_tran_veg      => pwf%qflx_tran_veg
    qflx_evap_veg      => pwf%qflx_evap_veg
    qflx_evap_soi      => pwf%qflx_evap_soi
    qflx_evap_can      => pwf%qflx_evap_can
    qflx_evap_tot      => pwf%qflx_evap_tot
    qflx_evap_grnd     => pwf%qflx_evap_grnd
    qflx_dew_grnd      => pwf%qflx_dew_grnd
    qflx_dew_snow      => pwf%qflx_dew_snow
    qflx_sub_snow      => pwf%qflx_sub_snow

    wtr_h2ocan         => pws%wtr_h2ocan
    wtr_qflx_prec_intr => pwf%wtr_qflx_prec_intr
    wtr_qflx_prec_grnd => pwf%wtr_qflx_prec_grnd
    wtr_qflx_rain_grnd => pwf%wtr_qflx_rain_grnd
    wtr_qflx_snow_grnd => pwf%wtr_qflx_snow_grnd
    wtr_qflx_snwcp_liq => pwf%wtr_qflx_snwcp_liq
    wtr_qflx_snwcp_ice => pwf%wtr_qflx_snwcp_ice
    wtr_qflx_tran_veg  => pwf%wtr_qflx_tran_veg
    wtr_qflx_evap_veg  => pwf%wtr_qflx_evap_veg
    wtr_qflx_evap_can  => pwf%wtr_qflx_evap_can
    wtr_qflx_evap_tot  => pwf%wtr_qflx_evap_tot
    wtr_qflx_evap_soi  => pwf%wtr_qflx_evap_soi
    wtr_qflx_evap_grnd => pwf%wtr_qflx_evap_grnd
    wtr_qflx_dew_grnd  => pwf%wtr_qflx_dew_grnd
    wtr_qflx_dew_snow  => pwf%wtr_qflx_dew_snow
    wtr_qflx_sub_snow  => pwf%wtr_qflx_sub_snow

    ! Tracer only variables

    RCanopyVapor       => pws%RCanopyVapor
    RLeafWaterSun      => pws%RleafWaterSun
    RLeafWaterSha      => pws%RleafWaterSha

    ! set up ones
    one1d(:)   = 1._r8


    ! Set local checking flag
    lfail = .false.

    ! Check state variables
    call Check1ddelta('H2OCAN',pdbg,lbp,ubp,etols,lfail,h2ocan        , wtr_h2ocan        (:,:) ,.false.)

    ! Check tracer only state variables
    call Check1ddelta('RCanVap' ,pdbg,lbp,ubp,etols,lfail,one1d       , RCanopyVapor      (:,:) ,.false.)
    call Check1ddelta('RLeafSun',pdbg,lbp,ubp,etols,lfail,one1d       , RLeafWaterSun     (:,:) ,.false.)
    call Check1ddelta('RLeafSha',pdbg,lbp,ubp,etols,lfail,one1d       , RLeafWaterSha     (:,:) ,.false.)

    ! Check flux variables
    call Check1ddelta('QPINTR',pdbg,lbp,ubp,etolf,lfail,qflx_prec_intr, wtr_qflx_prec_intr(:,:),.false.)
    call Check1ddelta('QPGRND',pdbg,lbp,ubp,etolf,lfail,qflx_prec_grnd, wtr_qflx_prec_grnd(:,:),.false.)
    call Check1ddelta('QRINTR',pdbg,lbp,ubp,etolf,lfail,qflx_rain_grnd, wtr_qflx_rain_grnd(:,:),.false.)
    call Check1ddelta('QRSNOW',pdbg,lbp,ubp,etolf,lfail,qflx_snow_grnd, wtr_qflx_snow_grnd(:,:),.false.)
    call Check1ddelta('QSNCAPLIQ',pdbg,lbp,ubp,etolf,lfail,qflx_snwcp_liq, wtr_qflx_snwcp_liq(:,:),.false.)
    call Check1ddelta('QSNCAPICE',pdbg,lbp,ubp,etolf,lfail,qflx_snwcp_ice, wtr_qflx_snwcp_ice(:,:),.false.)
    call Check1ddelta('QTRVEG',pdbg,lbp,ubp,etolf,lfail,qflx_tran_veg , wtr_qflx_tran_veg (:,:),.true.)
    call Check1ddelta('QEVVEG',pdbg,lbp,ubp,etolf,lfail,qflx_evap_veg , wtr_qflx_evap_veg (:,:),.true.)
    call Check1ddelta('QEVCAN',pdbg,lbp,ubp,etolf,lfail,qflx_evap_can , wtr_qflx_evap_can (:,:),.true.)
    call Check1ddelta('QEVSOI',pdbg,lbp,ubp,etolf,lfail,qflx_evap_soi , wtr_qflx_evap_soi (:,:),.true.)
    call Check1ddelta('QEVTOT',pdbg,lbp,ubp,etolf,lfail,qflx_evap_tot , wtr_qflx_evap_tot (:,:),.true.)
    call Check1ddelta('QEVGND',pdbg,lbp,ubp,etolf,lfail,qflx_evap_grnd, wtr_qflx_evap_grnd(:,:),.true.)
    call Check1ddelta('QDWGND',pdbg,lbp,ubp,etolf,lfail,qflx_dew_grnd , wtr_qflx_dew_grnd (:,:),.false.)
    call Check1ddelta('QSBSNO',pdbg,lbp,ubp,etolf,lfail,qflx_sub_snow , wtr_qflx_sub_snow (:,:),.false.)
    call Check1ddelta('QDWSNO',pdbg,lbp,ubp,etolf,lfail,qflx_dew_snow , wtr_qflx_dew_snow (:,:),.false.)

    ! Report
    if (lfail) then
      write(6,*) 'TracerCheckPFTDelta: ('//trim(subr)//')  ***** FAILED *****'
      if (lcheckstop) call endrun('ABORTED')
    else
!!      if (lcheckverb) write(6,*) 'TracerCheckPFTDelta: ('//trim(subr)//')  All OK'
      write(6,*) 'TracerCheckPFTDelta: ('//trim(subr)//')  All OK'
    endif

    ! update universal flag
    lfailany = lfailany .or. lfail
    return
  end subroutine TracerCheckPFTDelta


!=======================================================================
  subroutine Check1d(vname,begi,endi,etol,lfail,qprog,qtrac)
!-----------------------------------------------------------------------
! Performs a check on a one-dimensional array.
! Notice special check if both are NaN, which is possible on nstep=1.
!-----------------------------------------------------------------------
    use clm_varcon, only: spval
    use nanMod
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*), intent(in) :: vname	! name of variable
    integer,  intent(in)         :: begi,endi   ! lower/upper bounds
    real(r8), intent(in)         :: etol        ! tolerance on error
    real(r8), intent(in)         :: qprog(:)    ! prgnostic values
    real(r8), intent(in)         :: qtrac(:)    ! tracer values
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout)       :: lfail		! flag for failure
!------------------------- Local Variables -----------------------------
    logical lfailthis
    real(r8) error		! error (difference) 
    integer i                   ! array index
!-----------------------------------------------------------------------
    lfailthis = .false.
    if (lcheckverb) write(*,*) 'Checking: '//trim(vname)
    do i = begi,endi
       if (isnan(qprog(i)) .and. isnan(qtrac(i))) then
         error = 0.0			! Both are NaN
       else if (isnan(qprog(i))) then
         write(*,*) 'Check1d: qprog is NaN, qtrac is not'
         error = 10*etol		! ensure fail
       else if (isnan(qtrac(i))) then
         write(*,*) 'Check1d: qtrac is NaN, qprog is not'
         error = 10*etol		! ensure fail
       else if (qprog(i)==spval .and. qtrac(i)== spval) then
         error = 0.0
       else
         error = abs(qprog(i) - qtrac(i)/Rstnd(ixtest))
       end if
       if (error > etol) then
          lfailthis = .true.
          write(6,*) trim(vname)//' (Check failed): error =',i,error
          write(6,*) '   qprog,qtrac:',qprog(i),qtrac(i)/Rstnd(ixtest)
          go to 1
!          call endrun(' Check1d aborted!')
       end if
    end do
 1  continue
    if (lfailthis) write(6,*) 'TRACER CHECK FAILED: '//trim(vname)
    lfail = lfail .or. lfailthis
    return
  end subroutine Check1d

!=======================================================================
  subroutine Check2d(vname,begi,endi,begj,endj,etol,lfail,qprog,qtrac)
!-----------------------------------------------------------------------
! Performs a check on a two-dimensional array. 
! Notice special check if both are NaN, which is possible on nstep=1.
!-----------------------------------------------------------------------
    use clm_varcon, only: spval
    use nanmod
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*), intent(in) :: vname       ! name of variable
    integer,  intent(in)         :: begi,endi   ! dim 1 lower/upper bounds
    integer,  intent(in)         :: begj,endj   ! dim 2 lower/upper bounds
    real(r8), intent(in)         :: etol        ! tolerance on error
    real(r8), intent(in)         :: qprog(:,:)  ! prgnostic values
    real(r8), intent(in)         :: qtrac(:,:)  ! tracer values
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout) :: lfail           ! flag for failure
!------------------------- Local Variables -----------------------------
    logical lfailthis
    real(r8) error              ! error (difference) 
    integer i,j                 ! array indices
!-----------------------------------------------------------------------
    lfailthis = .false.
    if (lcheckverb) write(*,*) 'Checking: '//trim(vname)
    do j = begj,endj
      do i = begi,endi
         if (isnan(qprog(i,j)) .and. isnan(qtrac(i,j))) then
           error = 0.0			! Both are NaN
         else if (isnan(qprog(i,j))) then
           write(*,*) 'Check2d: qprog is NaN, qtrac is not'
           error = 10*etol              ! ensure fail 
         else if (isnan(qtrac(i,j))) then
           write(*,*) 'Check2d: qtrac is NaN, qprog is not'
           error = 10*etol              ! ensure fail
         else if (qprog(i,j)==spval .and. qtrac(i,j)== spval) then
           error = 0.
         else
           error = abs(qprog(i,j) - qtrac(i,j)/Rstnd(ixtest))
         end if
         if (error > etol) then
            write(6,*) trim(vname)//' (Check failed): error =',i,j,error
            write(6,*) ' qprog,qtrac:',qprog(i,j),qtrac(i,j)/Rstnd(ixtest)
            lfailthis = .true.
            go to 1
!!            call endrun(' Check2d aborted!')
         end if
      end do
    end do
 1  continue
    if (lfailthis) write(6,*) 'TRACER CHECK FAILED: '//trim(vname)
    lfail = lfail .or. lfailthis
    return
  end subroutine Check2d

!=======================================================================
  subroutine Check1ddelta(vname,idbg,begi,endi,etol,lfail,qprog,qtrac,ldiag)
!-----------------------------------------------------------------------
! Performs a check on a one-dimensional array.
! Notice special check if both are NaN, which is possible on nstep=1.
!-----------------------------------------------------------------------
    use clm_varcon, only: spval
    use nanMod
    use HydrologyIsotope, only : hydro_isotope, wiso_delta, pwspc
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
! Reasonable upper and lower bounds on delta values (for checking only)
  real(r8), parameter, dimension(pwspc) :: &
       delmin  = (/ -0.01    ,-500.      ,-80.     ,-40.     ,-1000.  /), &
       delmax  = (/ +0.01    ,+500.      ,+80.     ,+40.     ,+9000.  /)
!!       delmin  = (/ -0.01    ,-2000.      ,-2000.     ,-2000.     ,-1000.  /), &
!!       delmax  = (/ +0.01    ,+2000.      ,+2000.     ,+2000.     ,+9000.  /)
!------------------------- Input Arguments -----------------------------
    character(len=*), intent(in) :: vname	! name of variable
    integer,  intent(in)         :: begi,endi   ! lower/upper bounds
    integer , intent(in)	 :: idbg        ! index to write out
    real(r8), intent(in)         :: etol        ! tolerance on error
    real(r8), intent(in)         :: qprog(:)    ! prgnostic values
    real(r8), intent(in)         :: qtrac(:,:)  ! tracer values
    logical , intent(in)         :: ldiag       ! flag for only diagnostic ,no check
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout)       :: lfail		! flag for failure
!------------------------- Local Variables -----------------------------
    logical lfailthis
    real(r8) del(pwtrc)		! error (difference) 
    integer i,m                 ! array index
!-----------------------------------------------------------------------
    lfailthis = .false.
    if (.not. hydro_isotope) return
    if (lcheckverb) write(*,*) 'Checking: '//trim(vname)
    do i = begi,endi
 
       do m = 1, pwtrc
         if (isnan(qprog(i)) .and. isnan(qtrac(i,m))) then	! both are NaN, so OK
         else if (isnan(qprog(i)) .or. isnan(qtrac(i,m))) then	! only 1 is nan... problem
              write(6,*) 'Check1ddelta: NaN at point:',i
         else
           del(m) = wiso_delta(ispec(m),qtrac(i,m),qtrac(i,ixbase))
!!           del(m) = 1000.*(get_wratio(qtrac(i,m),qtrac(i,ixbase),qflxtiny)/Rstnd(ispec(m)) - 1.)
           if (abs(qprog(i)) .lt. etol .and. del(m) <-990.) then
                ! just nothing there
           else
             if (.not. ldiag) then
               if ((del(m) < delmin(ispec(m)) .or. del(m) > delmax(ispec(m))) &
                        .and. abs(qprog(i)) > 1000.*etol) then
                  if (.not. lfailthis) then	! just write first instance
                  write(6,*) trim(vname)//' (Delta check failed): error =', &
                      i,m,del(m),qprog(i),qtrac(i,ixbase),qtrac(i,m),etol
                  end if
                  lfailthis = .true.
                  go to 1
               end if
             end if
           end if
         end if
 1     continue
       end do

       if (ldeltadiag .and. i.eq.idbg) then
           if (.not. isnan(qprog(i))) then
               write(6,11) vname,qtrac(i,ixbase),del(2),del(3),del(2)-8.*del(3)
           endif
       endif

    end do 	
 11 format(a16,'   ','   q=',EN14.4,' dD=',f10.3,' d18O=',f10.3,' d=',f10.3)
    if (lfailthis) write(6,*) 'DELTA CHECK FAILED: '//trim(vname)
    lfail = lfail .or. lfailthis
    return
  end subroutine Check1ddelta

!=======================================================================
  subroutine Check2ddelta(vname,idbg,begi,endi,begj,endj,etol,lfail,qprog,qtrac)
!-----------------------------------------------------------------------
! Performs a check on a one-dimensional array.
! Notice special check if both are NaN, which is possible on nstep=1.
!-----------------------------------------------------------------------
    use clm_varcon, only: spval
    use nanMod
    use HydrologyIsotope, only : hydro_isotope, wiso_delta, pwspc
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
! Reasonable upper and lower bounds on delta values (for checking only)
  real(r8), parameter, dimension(pwspc) :: &
       delmin  = (/ -5.    ,-500.      ,-80.     ,-40.     ,-1000.  /), &
       delmax  = (/ +5.    ,+500.      ,+80.     ,+40.     ,+9000.  /)
!!       delmin  = (/ -0.01    ,-2000.      ,-2000.     ,-2000.     ,-1000.  /), &
!!       delmax  = (/ +0.01    ,+2000.      ,+2000.     ,+2000.     ,+9000.  /)
!------------------------- Input Arguments -----------------------------
    character(len=*), intent(in) :: vname	! name of variable
    integer,  intent(in)         :: begi,endi   ! lower/upper bounds
    integer,  intent(in)         :: begj,endj   ! dim 2 lower/upper bounds
    integer , intent(in)	 :: idbg        ! index to write out
    real(r8), intent(in)         :: etol        ! tolerance on error
    real(r8), intent(in)         :: qprog(:,:)  ! prognostic values
    real(r8), intent(in)         :: qtrac(:,:,:)! tracer values
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout)       :: lfail		! flag for failure
!------------------------- Local Variables -----------------------------
    logical lfailthis
    real(r8) del(pwtrc,begj:endj)		! error (difference) 
    integer i,j,m                 ! array index
!-----------------------------------------------------------------------
    lfailthis = .false.
    if (.not. hydro_isotope) return
    if (lcheckverb) write(*,*) 'Checking: '//trim(vname)
    do i = begi,endi
 
       do j = begj, endj
         do m = 1, pwtrc
           if (isnan(qprog(i,j)) .and. isnan(qtrac(i,j,m))) then	! both are NaN, so OK
           else if (isnan(qprog(i,j)) .or. isnan(qtrac(i,j,m))) then	! both are NaN, so OK
                write(6,*) 'Check2ddelta: NaN at point:',i
           else
             del(m,j) = wiso_delta(ispec(m),qtrac(i,j,m),qtrac(i,j,ixbase))
             if (abs(qprog(i,j)) .lt. etol .and. abs(del(m,j) + 1000.) < 1. ) then
                  ! just nothing there
             else
               if (del(m,j) < delmin(ispec(m)) .or. del(m,j) > delmax(ispec(m))) then
                  if (.not. lfailthis) then	! just write first instance
                    write(6,*) trim(vname)//' (Delta check failed): error =', &
                              i,j,m,del(m,j),qprog(i,j),qtrac(i,j,ixbase),qtrac(i,j,m)
                  end if
                  lfailthis = .true.
                  go to 1
               end if
             endif
           end if
 1       continue
         end do
       end do


       
       if (ldeltadiag .and. i.eq.idbg) then
         do j = begj,endj
           if (.not. isnan(qprog(i,j))) then
               write(6,11) vname,j,qtrac(i,j,ixbase),del(2,j),del(3,j),del(2,j)-8.*del(3,j)
           else
               write(6,12) vname,j
           endif
         end do
       endif

    end do 	
 11 format(a16,i3,'   q=',EN14.4,' dD=',f10.3,' d18O=',f10.3,' d=',f10.3)
 12 format(a16,i3,'   q= NaN')
    if (lfailthis) write(6,*) 'DELTA CHECK FAILED: '//trim(vname)
    lfail = lfail .or. lfailthis
    return
  end subroutine Check2ddelta

!=======================================================================
  function get_wratio(qtrac, qprog, qtiny)
!-----------------------------------------------------------------------
! Computes tracer ratio using "tiny" value
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
    real(r8) :: get_wratio
!---------------------------- Arguments --------------------------------
    real(r8), intent(in) :: qprog
    real(r8), intent(in) :: qtrac
    real(r8), intent(in) :: qtiny
!-----------------------------------------------------------------------

    ! Set tracer ratio to trivial for small mass
    ! THIS is the THE trick to getting high numericla precision

    if (abs(qprog) < qtiny) then
      get_wratio = 0._r8
      return
    end if

    if (qprog > 0._r8) then
!!       get_wratio = qtrac / (qprog+qtiny)	! introduces cumulative error
       get_wratio = qtrac / max(qprog,qtiny)
    else
!!       get_wratio = qtrac / (qprog-qtiny)	! introduces cumulative error
       get_wratio = qtrac / min(qprog,-qtiny)
    end if

    return
  end function get_wratio

!=======================================================================
  subroutine TracerCheckEqual(vname,sname,begi,endi,etol,lfail,qprog,qtrac)
!-----------------------------------------------------------------------
! Performs a check on 2 one-dimensional arrays, seeing if they're equal.
! Intention is to write out an error file with the name of the
! subroutine which dropped the ball.
! Notice special check if both are NaN, which is possible on nstep=1.
!-----------------------------------------------------------------------
    use clm_varcon, only: spval
    use nanMod
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    character(len=*), intent(in) :: vname       ! name of variable
    character(len=*), intent(in) :: sname       ! name of calling subroutine
    integer,  intent(in)         :: begi,endi   ! lower/upper bounds
    real(r8), intent(in)         :: etol        ! tolerance on error
    real(r8), intent(in)         :: qprog(:)    ! prgnostic values
    real(r8), intent(in)         :: qtrac(:)    ! tracer values
!---------------------- Input/Output Arguments -------------------------
    logical, intent(inout)       :: lfail       ! flag for failure
!------------------------- Local Variables -----------------------------
    logical lfailthis
    real(r8) error              ! error (difference) 
    integer i                   ! array index
!-----------------------------------------------------------------------
    lfailthis = .false.
    if (lcheckverb) write(*,*) 'Checking: '//trim(vname)
    do i = begi,endi
       if (isnan(qprog(i)) .and. isnan(qtrac(i))) then
         error = 0.0            ! Both are NaN
       else if (isnan(qprog(i))) then
         write(*,*) 'TracerCheckEqual: qprog is NaN, qtrac is not'
         error = 10*etol        ! ensure fail
       else if (isnan(qtrac(i))) then
         write(*,*) 'TracerCheckEqual: qtrac is NaN, qprog is not'
         error = 10*etol        ! ensure fail
       else if (qprog(i)==spval .and. qtrac(i)== spval) then
         error = 0.0
       else
         error = abs(qprog(i) - qtrac(i)/Rstnd(ixtest))
       end if
       if (error > etol) then
          lfailthis = .true.
          write(6,*) trim(vname)//' (Check failed): error =',i,error
          write(6,*) '   qprog,qtrac:',qprog(i),qtrac(i)/Rstnd(ixtest)
          go to 1
!          call endrun(' TracerCheckEqual aborted!')
       end if
    end do
 1  continue
    if (lfailthis) write(6,*) 'TRACERCHECKEQUAL FAILED: '//trim(vname)//' in '//trim(sname)
    lfail = lfail .or. lfailthis
    return

  end subroutine TracerCheckEqual

!=======================================================================
end module HydrologyTracer
