module HydrologyIsotope

!-----------------------------------------------------------------------
!BOP
! !MODULE: HydrologyIsotope
!
!  Notice not called WaterIsotope which is used in CAM.
!  Functions have SAME NAMES, and SAME INTERFACE!
!  
!
! !DESCRIPTION:
!   Handles all water isotope physics.
!
!   Interfaces with HydrologyTracer, for special case with isotope
!   fractionation. This is primitive and CAN NOT depend on
!   HydrologyTracer. If needed, can be done externally.
!
!   With the cpp definition NOFRAC set, all fractionations are set
!   to 1.0_r8 such that we have a simple unit test for the code.
!
! !REVISION HISTORY
!
!   Created David Noone <dcn@colorado.edu> - Thu Jun 16 18:45:47 MDT 2005
!   Updated David Noone <dcn@colorado.edu> - Fri Jul 13 15:40:56 MDT 2007
!
!EOP
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
!-----------------------------------------------------------------------
  use abortutils, only: endrun
!-----------------------------------------------------------------------
#undef NOFRAC
!-----------------------------------------------------------------------
  implicit none
  private
  save
!-----------------------------------------------------------------------

!
! Interfaces
!
  public :: HydrologyIsotopeInit    ! initializes this module

  public :: get_alpliq              ! look-up value for alpliq
  public :: get_alpice              ! look-up value for alpice

  public :: wiso_alpliq                 ! liquid/vapor quilibrium fractionation
  public :: wiso_alpice                 ! ice/vapor equilibrium fractionation
  public :: wiso_alpsol                 ! ice/liquid equilibrium fractionation
  public :: wiso_alpkin_m78             ! molecular/kinetic fractionation (Merlivat, 1978)
  public :: wiso_alpkin_mb96            ! molecular/kinetic fractionation (Mathieu and Bariac, 1996)
  public :: wiso_alpkin_mb96_2          ! molecular/kinetic fractionation (Mathieu and Bariac, 1996)
  public :: wiso_decay                  ! tendency from radioactive decay (HTO)
  public :: wiso_delta                  ! compute the delta value

  public :: wiso_get_rnat               ! retrive internal Rnat
  public :: wiso_get_rstd               ! retrive internal Rstd
  public :: wiso_get_rao2               ! retrive internal Rao2
  public :: wiso_get_fisub              ! retrive internal fisub
  public :: wiso_get_mwisp              ! retrive internal molecular weight
  public :: wiso_get_epsmw              ! retrive internal molecular weight ratio

!
! Module variables
!
  logical, public :: hydro_isotope = .true.     ! (namelist) flag to use this module
  logical, public :: tracer_constant_ratio = .false.  ! (namelist?)flag to do a constant ratio test for tracers
                        !! if tracer_constant_ratio = .true., you need to modify
                        !! rstd and rnat below to match (keep m=1 tracer to match
                        !! bulk water), and difrm for kinetic fractionation
                        !! this is only used to tell the fractionations to all
                        !! be 1
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  integer, parameter, public :: pwspc = 5       ! number of water species
!
! Species indicies
!
  integer, parameter, public :: ispundef = 0    ! Undefined (not this module)
  integer, parameter, public :: ispwater = 1    ! total water (same as H2O)
  integer, parameter, public :: isph2o   = 1    ! H2O
  integer, parameter, public :: isphdo   = 2    ! HDO
  integer, parameter, public :: isph218o = 3    ! H218O
  integer, parameter, public :: isph217o = 4    ! H217O
  integer, parameter, public :: isphto   = 5    ! HTO

!
! Exponents for kinetic fractionation (Riley et al 2002)
!
  real(r8), parameter, public :: ff_km = 1.0     ! kinetic exponent for molecular diffusion
  real(r8), parameter, public :: ff_kb = 0.67    ! kinetic exponent for boundary layer

!
! Private isotopic constants
!
!!constant fraction to preserve in debugging code:
!!  real(r8), dimension(pwspc), parameter, public :: &
!!      confrac = (/ 1._r8, 0.5_r8, 0.5_r8 , 1._r8, 1._r8 /)
!!
  character(len=8), dimension(pwspc), parameter :: & ! species names
      spnam  = (/'H216O'   ,'HD16O'   ,'H218O'   ,'H217O'   ,'HTO  '/)
!
! Physical constants for isotopic molecules (not atoms)
!
  integer , dimension(pwspc), parameter :: &  ! posible isotopic substitutions
       fisub = (/ 1._r8    ,2._r8     ,1._r8     ,1._r8     ,2._r8  /)

  real(r8), dimension(pwspc), parameter :: &  ! molecular weights
       mwisp = (/ 18.      ,19.       ,20.       ,19.       ,20.    /)  !original
!       mwisp = (/ 18.      ,18.       ,18.       ,19.       ,20.    /)  !constant frac testing

  real(r8), dimension(pwspc), parameter :: &  ! mol. weight ratio
!       epsmw = (/ 1.    ,19./18.   ,20./18.   ,19./18.   ,20./18 /)     !original
       epsmw = (/ 1.    ,1.   ,1.   ,19./18.   ,20./18 /)     !const frac testing

  real(r8), dimension(pwspc), parameter, public :: &  ! diffusivity ratio (note D/H, not HDO/H2O)
       difrm = (/ 1._r8    ,0.9757_r8    ,0.9727_r8 , 0._r8 , 0._r8/)       ! Merlivat 1978 (direct from paper) (DEFAULT)
!       difrm = (/ 1.    ,0.9836504 ,0.9686999 ,0.9836504 ,0.9686999  /)    ! kinetic theory
!       difrm = (/ 1.    ,1. ,1. ,0.9836504 ,0.9686999  /)    ! const frac test (above is orig)
!       difrm = (/ 1._r8    ,0.9836504_r8 ,0.9686999_r8 ,0.9686999_r8 ,0.9836504_r8  /)    ! this with expk
!       difrm = (/ 1._r8    ,0.9839    ,0.9691    ,0.9873    ,0.96796    /)    ! Cappa etal 2003

! Isotopic ratios in natural abundance (SMOW)
  real(r8), dimension(pwspc), parameter :: &  ! SMOW isotope ratios
!       rnat  = (/ 1.    ,155.76e-6,2005.20e-6  ,402.00e-6 ,77.88e-06/) !orig
       rnat  = (/ 1.    ,1. ,1.  ,402.00e-6 ,77.88e-06/)  !equal numerics
!       rnat  = (/ 1.    ,0.7 , 0.5  ,402.00e-6 ,77.88e-06/)  ! constant fraction test

! Prescribed isotopic ratios (largely arbitrary and tunable)
  real(r8), dimension(pwspc), parameter :: &  ! model standard isotope ratio
!       rstd  = (/ 1.    ,155.76e-6,2005.20e-6  ,402.00e-6 , 77.88e-06/)   ! natural abundance
       rstd  = (/ 1._r8    ,1.0_r8    ,1.0_r8    ,1.0_r8    ,1.0_r8     /)  ! equal numerics
!       rstd  = (/ 1._r8    ,0.7_r8    ,0.5_r8    ,1.0_r8    ,1.0_r8     /)  ! constant fraction test

! Isotope enrighment of atmospheric oxygen (Bender 1999, triple-isotope)
  real(r8), dimension(pwspc), parameter :: &  ! mean ocean surface enrichent 
       bao2  = (/ 1.    ,1.0       ,0.97704   ,0.988222  ,1.0        /)

! Coefficients for fractionation: ln(alpha) = A/T2 + B/T + C
! From Majoube, 1971a:
!  real(r8), parameter, dimension(pwspc) :: &  ! liquid/vapour
!       alpal_real = (/0.        , 24.844e+3, 1.137e+3 , 1.137e+3 , 24.844e+3 /), & !orig
!       alpbl_real = (/0.        ,-76.248   ,-0.4156   ,-0.4156   ,-76.248    /), &
!       alpcl_real = (/0.        , 52.612e-3,-2.0667e-3,-2.0667e-3, 52.612e-3 /), &
!       alpal = (/0.        , 0. , 0. , 1.137e+3 , 24.844e+3 /), & !const frac test
!       alpbl = (/0.        , 0. , 0. ,-0.4156   ,-76.248    /), &
!       alpcl = (/0.        , 0. , 0. ,-2.0667e-3, 52.612e-3 /)

!From Horita and Wesolowski, 1994:
  real(r8), parameter, dimension(pwspc) :: &  ! liquid/vapour
      alpal_real = (/ 0._r8,  1158.8e-12_r8   , 0.35041e+6_r8, 0.35041e+6_r8, 1158.8e-12_r8 /), &
      alpbl_real = (/ 0._r8, -1620.1e-9_r8    , -1.6664e+3_r8, -1.6664e+3_r8, -1620.1e-9_r8 /), &
      alpcl_real = (/ 0._r8,   794.84e-6_r8   ,     6.7123_r8, 6.7123_r8, 794.84e-6_r8 /), &
      alpdl_real = (/ 0._r8,  -161.04e-3_r8   ,  -7.685e-3_r8, -7.685e-3_r8, -161.04e-3_r8   /), &
      alpel_real = (/ 0._r8,     2.9992e+6_r8 ,         0._r8, 0._r8, 2.9992e+6_r8 /), &
      alpal = (/ 0._r8, 0._r8, 0._r8, 0.35041e+6_r8, 1158.8e-12_r8 /), &
      alpbl = (/ 0._r8, 0._r8, 0._r8, -1.6664e+3_r8, -1620.1e-9_r8 /), &
      alpcl = (/ 0._r8, 0._r8, 0._r8, 6.7123_r8, 794.84e-6_r8      /), &
      alpdl = (/ 0._r8, 0._r8, 0._r8, -7.685e-3_r8, -161.04e-3_r8  /), &
      alpel = (/ 0._r8, 0._r8, 0._r8, 0._r8, 2.9992e+6_r8          /)  

! From ISOCAM3:
!  real(r8), parameter, dimension(pwspc) :: &  ! ice/vapour
!       alpai_real = (/0.        , 16288.   , 0.       , 0.       , 16288.    /), & !orig
!       alpbi_real = (/0.        , 0.       , 11.839   , 11.839   , 0.        /), &
!       alpci_real = (/0.        ,-9.34e-2  ,-28.224e-3,-28.224e-3,-9.34e-2   /), &
!       alpai = (/0.        , 0.  , 0.   , 0.       , 16288.    /), & !const frac test
!       alpbi = (/0.        , 0.  , 0.   , 11.839   , 0.        /), &
!       alpci = (/0.        , 0.  , 0.   ,-28.224e-3,-9.34e-2   /)

!From Merlivat & Nief, 1967 for HDO, and Majoube, 1971b for H218O:
  real(r8), parameter, dimension(pwspc) :: &  ! ice/vapour
       alpai_real = (/0.        , 16289.   , 0.       , 0.       , 16289.    /), & !orig
       alpbi_real = (/0.        , 0.       , 11.839   , 11.839   , 0.        /), &
       alpci_real = (/0.        ,-9.45e-2  ,-28.224e-3,-28.224e-3,-9.45e-2   /), &
       alpai = (/0.        , 0.  , 0.   , 0.       , 16288.    /), & !const frac test
       alpbi = (/0.        , 0.  , 0.   , 11.839   , 0.        /), &
       alpci = (/0.        , 0.  , 0.   ,-28.224e-3,-9.34e-2   /)

  real(r8), parameter, dimension(pwspc) :: &  ! solid/liquid
       alpas_real = (/1.0       ,1.02      ,1.003     ,1.003     ,1.02       /), & !orig
       alpas = (/1.0       ,1.0      ,1.0     ,1.003     ,1.02       /)  !const frac test

! Modifier for non-standard species EQUILIBRIUM 
  real(r8), parameter, dimension(pwspc) :: &  ! specices fractionation modifier
       expa  = (/ 1.    ,1.0       ,1.0      ,0.525      ,2.0        /)
!       expa  = (/ 1._r8    ,1.0       ,1.0      ,0.52441    ,2.0        /)

! Other special constants
! 
  real(r8), parameter :: hlhto = 3.88e+8 ! halflife of HTO (12.43 years) [s]
!
! Internal storage arrays for fractionation factors
!
  integer , parameter :: ptres = 10000	! look-up table resolution
  real(r8), parameter :: tmin = 180.
  real(r8), parameter :: tmax = 320.
  real(r8), parameter :: ptdel = ptres/(tmax-tmin)

  real(r8), public :: alpkm(pwspc)
  real(r8), public :: alpkb(pwspc)
  real(r8)         :: alpliq(ptres,pwspc)
  real(r8)         :: alpice(ptres,pwspc)

!=======================================================================
CONTAINS

!=======================================================================
  subroutine HydrologyIsotopeInit()
!-----------------------------------------------------------------------
! Initializes module
!-----------------------------------------------------------------------
    use perf_mod
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
    integer isp, i
    real(r8) dt, tk
!-----------------------------------------------------------------------
    if (.not. hydro_isotope) then
       call endrun('HydrologyIsotopeInit: hydroisotope not set .true.')
    end if

#ifdef NOFRAC
    write(6,*) '***************************************************'
    write(6,*) ' '
    write(6,*) '          NO FRACTIONATION OPTION ENABLED '
    write(6,*) ' '
    write(6,*) '***************************************************'
#endif

    ! Initialize kinetic fractionation factors
    do isp = 1, pwspc
      if(tracer_constant_ratio) then
        alpkm(isp) = 1._r8
        alpkb(isp) = 1._r8
      else
        alpkm(isp) = wiso_alpkin_m78(isp,ff_km,.false.)
        alpkb(isp) = wiso_alpkin_m78(isp,ff_kb,.false.)
      end if
    end do

    ! Initialize equilibrium fractionation coefficients look-up table
    dt = (tmax - tmin)/dble(ptres)
    do isp = 1, pwspc
      do i = 1, ptres
        tk = tmin + (i-1)*dt
        if(tracer_constant_ratio) then
          alpliq(i,isp) = 1._r8
          alpice(i,isp) = 1._r8
        else
          alpliq(i,isp) = wiso_alpliq(isp,tk,.true.)
          alpice(i,isp) = wiso_alpice(isp,tk,.true.)
        end if
      end do
    end do

    write(6,*) 'HydrologyIsotopeInit: Using water isotopes for tracers.'
    return
  end subroutine HydrologyIsotopeInit

!=======================================================================
  function get_alpliq(isp,tk)
!-----------------------------------------------------------------------
! Gets fractionation factor from look-up table
!-----------------------------------------------------------------------
    integer , intent(in) :: isp
    real(r8), intent(in) :: tk		! tempersture
    real(r8) get_alpliq
    integer i
!-----------------------------------------------------------------------
    i = ptdel*(tk-tmin)
    i = max(i,1)
    i = min(i,ptres)
    get_alpliq = alpliq(i,isp)
    return
  end function get_alpliq

!=======================================================================
  function get_alpice(isp,tk)
!-----------------------------------------------------------------------
! Gets fractionation factor from look-up table
!-----------------------------------------------------------------------
    integer , intent(in) :: isp
    real(r8), intent(in) :: tk		! tempersture
    real(r8) get_alpice
    integer i
!-----------------------------------------------------------------------
    i = ptdel*(tk-tmin)
    i = max(i,1)
    i = min(i,ptres)
    get_alpice = alpice(i,isp)
    return
  end function get_alpice

!=======================================================================
  function wiso_alpliq(isp,tk,lfrac)
!-----------------------------------------------------------------------
! Purpose: return liquid/vapour fractionation from equation
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 10:59:13 MDT 2003
!-----------------------------------------------------------------------
    use perf_mod
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp  ! species indes
    real(r8), intent(in)        :: tk   ! temperature (k)
    logical , intent(in), optional :: lfrac ! flag for actually fractionating (true=fractionate)
    real(r8) :: wiso_alpliq             ! return fractionation
!-----------------------------------------------------------------------
    call t_startf('alpliq')
!
    ! If lfrac is not present OR it is present and it is false, frac fact = 1
    if(.not. present(lfrac) .or. (present(lfrac) .and. .not. lfrac)) then
      if (isp <= isph2o) then
        wiso_alpliq = 1._r8
      else
        !Majoube, 1971:
!        wiso_alpliq = exp(alpal(isp)/(tk**2) + alpbl(isp)/tk + alpcl(isp))
        !Horita and Wesolowski, 1994:
        if(isp == isphdo) then !HDO has different formulation:
          wiso_alpliq = exp(alpal(isp)*tk**3 + alpbl(isp)*tk**2 + alpcl(isp)*tk + alpdl(isp) + alpel(isp)/tk**3)
        else
          wiso_alpliq = exp(alpal(isp)/tk**3 + alpbl(isp)/tk**2 + alpcl(isp)/tk + alpdl(isp))
        end if
        wiso_alpliq = wiso_alpliq**expa(isp)
#ifdef NOFRAC
        wiso_alpliq = 1._r8
#endif
      end if
    ! only other possibility is lfrac is present and true, in which fractionate away!
    else
      if (isp <= isph2o) then
        wiso_alpliq = 1._r8
      else
        !Majoube, 1971:
!        wiso_alpliq = exp(alpal_real(isp)/(tk**2) + alpbl_real(isp)/tk + alpcl_real(isp))
        !Horita and Wesolowski, 1994:
        if(isp == isphdo) then !HDO has different formulation:
          wiso_alpliq = exp(alpal_real(isp)*tk**3 + alpbl_real(isp)*tk**2 + alpcl_real(isp)*tk + alpdl_real(isp) + alpel_real(isp)/tk**3)
        else
          wiso_alpliq = exp(alpal_real(isp)/tk**3 + alpbl_real(isp)/tk**2 + alpcl_real(isp)/tk + alpdl_real(isp))
        end if
        wiso_alpliq = wiso_alpliq**expa(isp)
#ifdef NOFRAC
        wiso_alpliq = 1._r8
#endif
      end if
    end if
    
!DEBUG
!wiso_alpliq=1._r8

    if(tracer_constant_ratio) wiso_alpliq=1._r8

    call t_stopf('alpliq')
    return
  end function wiso_alpliq

!=======================================================================
  function wiso_alpice(isp,tk,lfrac)
!-----------------------------------------------------------------------
! Purpose: return ice/vapour fractionation from equation
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
    use perf_mod
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp  ! species indes
    real(r8), intent(in)        :: tk   ! temperature (k)
    logical , intent(in), optional :: lfrac ! flag for actually fractionating (true=fractionate)
    real(r8) :: wiso_alpice             ! return fractionation
!-----------------------------------------------------------------------
    call t_startf('alpice')

    ! If lfrac is not present OR it is present and it is false, frac fact = 1
    if(.not. present(lfrac) .or. (present(lfrac) .and. .not. lfrac)) then
      if (isp <= isph2o) then
        wiso_alpice = 1._r8
      else
        wiso_alpice = exp(alpai(isp)/(tk**2) + alpbi(isp)/tk + alpci(isp))
        wiso_alpice = wiso_alpice**expa(isp)
#ifdef NOFRAC
        wiso_alpice = 1._r8
#endif
      endif
    else ! only other case possible is lfrac present and true, so fractionate!
      if (isp <= isph2o) then
        wiso_alpice = 1._r8
      else
        wiso_alpice = exp(alpai_real(isp)/(tk**2) + alpbi_real(isp)/tk + alpci_real(isp))
        wiso_alpice = wiso_alpice**expa(isp)
#ifdef NOFRAC
        wiso_alpice = 1._r8
#endif
      end if
    end if

!DEBUG
!wiso_alpice=1._r8

    if(tracer_constant_ratio) wiso_alpice=1._r8

    call t_stopf('alpice')
    return
  end function wiso_alpice

!=======================================================================
  function wiso_alpsol(isp,lfrac)
!-----------------------------------------------------------------------
! Purpose: return ice/liquid fractionation from look-up tables
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp  ! species indes
!    real(r8), intent(in)        :: tk   ! temperature (k)
    logical , intent(in), optional :: lfrac ! flag for actually fractionating (true=fractionate)
    real(r8) :: wiso_alpsol               ! return fractionation
!-----------------------------------------------------------------------
    ! If lfrac is not present OR it is present and it is false, frac fact = 1
    if(.not. present(lfrac) .or. (present(lfrac) .and. .not. lfrac)) then
      if (isp <= isph2o) then
        wiso_alpsol = 1._r8
      else
        wiso_alpsol = alpas(isp)
        wiso_alpsol = wiso_alpsol**expa(isp)
#ifdef NOFRAC
        wiso_alpsol = 1._r8
#endif
      end if
    else !fractionate away! (lfrac must be present and true)
      if (isp <= isph2o) then
        wiso_alpsol = 1._r8
      else
        wiso_alpsol = alpas_real(isp)
        wiso_alpsol = wiso_alpsol**expa(isp)
#ifdef NOFRAC
        wiso_alpsol = 1._r8
#endif
      end if
    end if

!DEBUG
!wiso_alpsol=1._r8

    if(tracer_constant_ratio) wiso_alpsol=1._r8

    return
  end function wiso_alpsol

!=======================================================================
  function wiso_alpkin_m78(isp,ff,lfrac)
!-----------------------------------------------------------------------
! Calculate a kinetic fractionation based on the ratio of diffusivities
! raised to some power. this basically exists so we can use comiso.h in
! lsm without buggering up things by having params.h in lsm code.
! use value ff = 1.0 for pure diffusion, while for mixed
! molecular/turbulent boundary layers, ff = 2/3 (0.6667, etc).
!   Calculates based on Merlivat, 1978
!-----------------------------------------------------------------------
    use perf_mod
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
    integer , intent(in) :: isp		! species index (isp(m))
    real(r8), intent(in) ::  ff		! exponent (typically 1 or 2/3)
    logical , intent(in), optional :: lfrac ! flag for actually fractionating (true=fractionate)
    real(r8) :: wiso_alpkin_m78             ! return value
!-----------------------------------------------------------------------
    call t_startf('alpkin_m78')

!HERE2
    ! If lfrac is not present OR it is present and it is false, frac = 1
    if(.not. present(lfrac) .or. (present(lfrac) .and. .not. lfrac)) then
      wiso_alpkin_m78 = 1._r8
    else ! only other case possible is lfrac present and true, so fractionate!
      if (isp == isph2o) then
        wiso_alpkin_m78 = 1._r8

      else
        wiso_alpkin_m78 = (1/difrm(isp))**ff

#ifdef NOFRAC           /* diagnostic mode no fractionation */
        wiso_alpkin_m78 = 1.0_r8
#endif
      endif
    endif

!DEBUG
!wiso_alpkin_m78=1._r8

    if(tracer_constant_ratio) wiso_alpkin_m78=1._r8

    call t_stopf('alpkin_m78')
    return
  end function wiso_alpkin_m78

!=======================================================================
  function wiso_alpkin_mb96(isp,watsurf,watmin,watsat,lfrac)
!-----------------------------------------------------------------------
! Calculate a kinetic fractionation based on Mathieu and Bariac, 1996:
!   alp = (D/Di)^nk
!    nk = ((W1-Wr)*0.5+(Wsat-W1)*1.0) / (Wsat-Wr)
!   W1 = volumetric top soil layer water
!   Wr = residual (minimum) water content (watmin)
!   Wsat = saturated water content (cps%watsat(c,1))
!-----------------------------------------------------------------------
    use perf_mod
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
    integer , intent(in) :: isp     ! species index (isp(m))
    real(r8), intent(in) :: watsurf ! volumetric surface water
    real(r8), intent(in) :: watmin  ! volumetric minimum water
    real(r8), intent(in) :: watsat  ! volumetric saturated water
    logical , intent(in), optional :: lfrac ! flag for actually fractionating (true=fractionate)
    real(r8) :: wiso_alpkin_mb96    ! return value

    real(r8) :: nk                  ! exponent

!-----------------------------------------------------------------------
    call t_startf('alpkin_mb96')
    
    ! If lfrac is not present OR it is present and it is false, frac = 1
    if(.not. present(lfrac) .or. (present(lfrac) .and. .not. lfrac)) then
      wiso_alpkin_mb96 = 1._r8
    else ! only other case possible is lfrac present and true, so fractionate!
      if (isp == isph2o) then
        wiso_alpkin_mb96 = 1._r8

      else
        nk = ((watsurf-watmin)*0.5 + (watsat-watsurf)*1.0 ) / (watsat-watmin)
        wiso_alpkin_mb96 = (1/difrm(isp))**nk

      endif
    endif

!DEBUG
!wiso_alpkin_mb96 = 1._r8

    if(tracer_constant_ratio) wiso_alpkin_mb96=1._r8
!DEBUG
!write(6,*) 'MB96',isp,wiso_alpkin_mb96

    call t_stopf('alpkin_mb96')
    return
  end function wiso_alpkin_mb96

!=======================================================================
  function wiso_alpkin_mb96_2(isp,r1,n1,r2,n2,r3,n3)
!-----------------------------------------------------------------------
! Calculate a kinetic fractionation based on Mathieu and Bariac, 1996:
!   alp_k = r_total,iso / r_total = sum( ri*(D/D_iso)^ni ) / sum( ri )
! Input:  water species, isp; resistances, r1,2,3; exponents, n1,2,3,
! representing turbulent (0.5) to diffusive (1.0) balance (.67 = laminar
! boundary layer)
! Output: wiso_alpkin_mb96_2, kinetic fractionation factor
!           (heavy isotopes, should be > 1)
!-----------------------------------------------------------------------
    use perf_mod
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
    integer , intent(in) :: isp         ! species index (isp(m))
    real(r8), intent(in) :: r1,n1       ! need at least one resistance and exponent
    real(r8), intent(in), optional :: r2,n2 ! additional resistance and exponent
    real(r8), intent(in), optional :: r3,n3 ! additional resistance and exponent
    real(r8) :: wiso_alpkin_mb96_2      ! return value
    real(r8) :: num,den                 ! make it easier 
!-----------------------------------------------------------------------
    call t_startf('alpkin_mb96_2')
    if (isp == isph2o) then
      wiso_alpkin_mb96_2 = 1._r8

    else
      num = 1._r8
      den = 1._r8
      wiso_alpkin_mb96_2 = num/den

    endif

!DEBUG
!wiso_alpkin_mb96_2 = 1._r8

    if(tracer_constant_ratio) wiso_alpkin_mb96_2=1._r8

    call t_stopf('alpkin_mb96_2')
    return
  end function wiso_alpkin_mb96_2

!=======================================================================
  subroutine wiso_decay(isp,dtime,q,dqdcy)
!-----------------------------------------------------------------------
! Impliments radioactive decay (for tritium, etc)
!-----------------------------------------------------------------------
  integer , intent(in)  :: isp		! species index (MUST BE isphto, for now)
  real(r8), intent(in)  :: q		! mass of stuff
  real(r8), intent(in)  :: dtime	! time interval of decay
  real(r8), intent(out) :: dqdcy	! change in mass due to decay
!-----------------------------------------------------------------------
    if (isp /= isphto) call endrun('(wiso_decay) isp /= isphto: TRITIUM ONLY')
    dqdcy = q * (0.5**(dtime/hlhto) - 1.) / dtime
    return
  end subroutine wiso_decay

!=======================================================================
  function wiso_delta(isp,qiso,qtot)
!-----------------------------------------------------------------------
! Purpose: Compute isotopic delta value from masses
! Author David Noone <dcn@caltech.edu> - Tue Jul  1 08:32:45 MDT 2003
!-----------------------------------------------------------------------
    integer, intent(in)  :: isp         ! species index
    real(r8),intent(in)  :: qiso        ! isotopic mass
    real(r8),intent(in)  :: qtot        ! isotopic mass
    real(r8) :: wiso_delta              ! return value
!-----------------------------------------------------------------------
    wiso_delta = 1000.*(wiso_ratio(isp,qiso,qtot)/Rstd(isp) - 1.)
    return
  end function wiso_delta

!=======================================================================
  function wiso_ratio(isp,qiso,qtot)
!-----------------------------------------------------------------------
! Purpose: Compute isotopic ratio from masses, with numerical checks
! Author David Noone <dcn@caltech.edu> - Tue Jul  1 08:32:45 MDT 2003
!-----------------------------------------------------------------------
    integer, intent(in)  :: isp         ! species index
    real(r8),intent(in)  :: qiso        ! isotopic mass
    real(r8),intent(in)  :: qtot        ! isotopic mass
    real(r8) :: wiso_ratio              ! return value
!-----------------------------------------------------------------------
    real(r8) :: qtiny = 1.e-16
!-----------------------------------------------------------------------
    if (abs(qtot) < qtiny) then
       wiso_ratio = 0.0
       return
    endif

    if (qtot > 0) then
      wiso_ratio = qiso/max(qtot, qtiny)
    else
      wiso_ratio = qiso/min(qtot,-qtiny)
    end if
!!    wiso_ratio = espmw(isp)*wiso_ratio/fisum(isp)      ! correct!
    return
  end function wiso_ratio


!=======================================================================
  function wiso_get_rnat(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal rnat variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:14 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_rnat             ! return isotope ratio
!-----------------------------------------------------------------------
    wiso_get_rnat = rnat(isp)
    return
  end function wiso_get_rnat

!=======================================================================
  function wiso_get_rstd(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal Rstd variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:14 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_rstd             ! return isotope ratio
!-----------------------------------------------------------------------
    wiso_get_rstd = rstd(isp)
    return
  end function wiso_get_rstd

!=======================================================================
  function wiso_get_rao2(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal Rao2 variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:04 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_rao2             ! return isotope ratio
!-----------------------------------------------------------------------
    wiso_get_rao2 = bao2(isp)*rstd(isp)
    return
  end function wiso_get_rao2

!=======================================================================
  function wiso_get_fisub(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal fisub variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:28:52 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp         ! species index
    real(r8) :: wiso_get_fisub           ! return number of substitutions
!-----------------------------------------------------------------------
    wiso_get_fisub = fisub(isp)
    return
  end function wiso_get_fisub

!=======================================================================
  function wiso_get_mwisp(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal mwwsp variable, based on species index
! Author: David Noone <dcn@colorado.edu> - Tue Feb  3 14:49:08 MST 2004
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp         ! species index
    real(r8) :: wiso_get_mwisp           ! return molecular weight
!-----------------------------------------------------------------------
    wiso_get_mwisp = mwisp(isp)
    return
  end function wiso_get_mwisp

!=======================================================================
  function wiso_get_epsmw(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal epsmw variable, based on species index
! Author: David Noone <dcn@colorado.edu> - Tue Feb  3 14:49:08 MST 2004
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp         ! species index
    real(r8) :: wiso_get_epsmw           ! return molecular weight
!-----------------------------------------------------------------------
    wiso_get_epsmw = epsmw(isp)
    return
  end function wiso_get_epsmw


!=======================================================================
end module HydrologyIsotope
