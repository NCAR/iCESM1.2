module water_tracer_vars
!-----------------------------------------------------------------------
!This module contains all of the generic variables and arrays needed for
!water tracers and water isotopes.  The only real reason this module
!exists is to prevent dependency errors during compliation that occur if 
!these variables are included in water_tracers.F90
!
!Written by:  Jesse Nusbaumer <nusbaume@colorado.edu> - July, 2013
!
!Code originally from authors of water_tracers.F90
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst
  use water_isotopes, only: pwtspec
  use water_types,    only: pwtype

  implicit none

  private
  save

!------------------- Module Variable Declarations -----------------------
integer, parameter, public    :: WTRC_MAX_CNST  = 700   ! Maximum number of water tracers allowed
                                                        !NOTE:  Currently set-up for 100 different water
                                                        !species and/or tags (7*100 = 700). - JN
integer, parameter, public    :: WTRC_WSET_STD  = 1     ! water set index for "regular" water


! Namelist variables
  logical, public    :: trace_water               = .false.     ! set true to activate [off]
  logical, public    :: wisotope                  = .false.     ! activate water isotopes [off]

  logical, public    :: wtrc_lh2oadj              = .false.     ! adjust tracer H20 to Q
!  logical, public    :: wtrc_lnomfix              = .true.      ! do not apply usual mass fixer (eul core)
  logical, public    :: wtrc_lzmlin               = .true.      ! linear interpolation for zm midpoints (else log)
!  logical, public    :: wtrc_cldw_adv             = .false.     ! true for advected, false for non-advected
  logical, public    :: wtrc_warn_only            = .true.      ! true for message only, no endrun
  logical, public    :: wtrc_add_cvprecip         = .false.     ! true to add QRAINC and QSNOWC, if not done by microphysics
  logical, public    :: wtrc_add_stprecip         = .false.     ! true to add QRAINS and QSNOWS, if not done by microphysics
  logical, public    :: wtrc_alpha_kinetic        = .false.     ! include kinetic effects in fractionation
  logical, public    :: wtrc_check_total_h2o      = .false.     ! check total mass conservation
  logical, public    :: wtrc_check_show_types     = .false.     ! check total mass conservation
  logical, public    :: wtrc_detrain_in_macrop    = .false.     ! do detrainment of isotopes in macrop (approx only)

!
  integer, public    :: wtrc_niter                = 10          ! number of iterations to use when applying process rates
  integer, public    :: wtrc_citer                = 20          ! number of iterations in dicm (10 < nitr < 1000) 
  real(r8), public   :: wtrc_qchkmin              = 1.e-14_r8   ! minimum relative difference to trigger check failure
  real(r8), public   :: wtrc_qmin                 = 1.e-18_r8   ! smalles support mixing ratio (matches qsmall in MG microphysics)
  real(r8), public   :: wtrc_fixed_alpha(pwtspec) = 1._r8       ! default standard fractionation factor, used when wisotope is false
  real(r8), public   :: wtrc_fixed_rstd(pwtspec)  = 1._r8       ! default standard isotope ratio, used when wisotope is false

  character(len=32), public      :: water_tracer_model = "none"

!
! The following tables define the water isotopes tracers that will be created. A maximum
! of WTRC_MAX_CNST tracers can be defined.
!
! Tracer names, no more than 5 characters for history files
  character(len=8), dimension(WTRC_MAX_CNST), public  :: wtrc_names  = ""  !Made public in order to generate output
  character(len=8), dimension(WTRC_MAX_CNST), public  :: wtrc_species_names  = ""
  character(len=8), dimension(WTRC_MAX_CNST), public  :: wtrc_type_names     = ""
  character(len=8), dimension(WTRC_MAX_CNST), public  :: wtrc_srfvap_names   = ""
  character(len=8), dimension(WTRC_MAX_CNST), public  :: wtrc_srfpcp_names   = ""
  character(len=8), dimension(WTRC_MAX_CNST), public  :: wtrc_tag_names      = ""

  character(len=8), dimension(WTRC_MAX_CNST), public  :: wtrc_out_names      = ""

!
! Constituent names, for the bulk (non-isotopic) water
  character(len=8), dimension(pwtype), parameter, public :: & ! constituent names 
      wtrc_bulk_names =   (/ 'Q       ', 'CLDLIQ  ', 'CLDICE  ', 'QRAINS  ', 'QSNOWS  ', 'QRAINC  ', 'QSNOWC  ' /)
!
! Constituent indices for the bulk (non-isotopic) water
  integer, public   :: wtrc_bulk_indices(pwtype)


! These are derived off of the fields provided in the namelist files.
!
! The variables specify a certain set of water tracers.
  integer, public :: wtrc_ncnst                         ! Number of water tracers
  integer, public :: wtrc_indices(WTRC_MAX_CNST)        ! Tracer constituent indices
  integer, public :: wtrc_species(WTRC_MAX_CNST)        ! Tracer species enumerations
  integer, public :: wtrc_types(WTRC_MAX_CNST)          ! Tracer water type enumerations
  logical, public :: wtrc_is_tag(WTRC_MAX_CNST)         ! Tracer is a tag rather than an isotopolgue?

!
! Index arrays for all specified water tracers
!
! Sort by Water Type
  integer, public :: wtrc_ntype(pwtype)                  ! number of each water type
  integer, public :: wtrc_iatype(WTRC_MAX_CNST, pwtype)  ! index arrays for water type

! Sort by Species
  integer, public :: wtrc_nspec(pwtspec)                 ! number of each species
  integer, public :: wtrc_iaspec(WTRC_MAX_CNST, pwtspec) ! index arrays for each species

! Sort by Tracer water sets (groups of all water types for the same species)
  integer, public :: wtrc_nwset                                    ! number of each water set
  integer, public :: wtrc_iawset(pwtype, WTRC_MAX_CNST / pwtype)   ! index arrays for the water sets  
  integer, public :: wtrc_srfpcp_indices(pwtype, WTRC_MAX_CNST / pwtype)  ! pbuf index arrays for surface precipitation tracers

! Surface fields for coupling to surface models.
  integer, public :: wtrc_nsrfvap                       ! number of surface water vapor tracers 
  integer, public :: wtrc_iasrfvap(WTRC_MAX_CNST)       ! index arrays for surface vapor tracers 
  integer, public :: wtrc_nsrfpcp                       ! number of surface precipitation tracers
  integer, public :: wtrc_iasrfpcp(WTRC_MAX_CNST)       ! index arrays for surface precipitation tracers

!
! Configuration pointers/indices
!
! NOTE: These should probably be moved to constituents.
  integer, public :: iwater(pcnst)     ! flag for water type (see water_types)
  integer, public :: iwspec(pcnst)     ! flag for water (isotope) species (see water_isotopes)
  logical, public :: iwistag(pcnst)    ! flag for tagged water

end module water_tracer_vars
