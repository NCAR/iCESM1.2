module water_tracers
!-----------------------------------------------------------------------
!
! Provide core functionality for water tracers.
!
! This module works in tandem with "waterisotope" which specifically
! teats the istopic fractionation.
!
! All interface routine are identified by wtrc_*, etc.
!
! Indexing ASSUMES (as everywhere else in CAM), normal water vapour is
! m=1. Cloud liquid and ice are probably m=2 and m=3... but not
! assumed, as they are correctly registered.
!
! DEFAULT CONFIGURATION is for 3 additional tracers with each of 3 phases, to
! parallel the base CAM water prognosis. This default is invoked by setting the
! namelist variable "wisotope". If not set, don't know what to do, so 
! complain and crash.
!
! Note total vapout (Q) is registered as wet, even though it is treated
! as dry. This means the PD coupling is diferent from the vertical
! diffusion. Do get around this we make all water species wet, but have
! spacial cased for vapour, as in main cam code.
!
!
! Code here based on chemistry.F90 and cldcond.F90.
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 15:22:48 MDT 2003
!
! Modified (slightly) for CAM5 by Jesse Nusbaumer <nusbaume@colorado.edu> - 4-3-11
!
! 2/4/2012 - bardeenc : Significant refactoring to accomplish several goals:
!   - broke out water_types
!   - configuration is driven by configure and build_namelist, rather than
!     being hard coded
!   - added support for rain and snow
!   - made an explicit distinction between tagging and isotopes
!       to allow for better checking
!   - added explicit definition of surface fields
!   - made more granular, but extensible list of tracers to ease
!     looping over subsets of the tracers
!
!   For precipitation, still need to decide if we need to track
!   convective and large scale precipitation separately. Also need
!   deal with precipitation totals at the surface and getting those
!   coupled to the surface models.
!
!   Will the interface for the process rates from stratiform also
!   work for convection and any other place where phase changes or
!   transport is happening?
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver
  use constituents,   only: pcnst
  use abortutils,     only: endrun
  use cam_logfile,    only: iulog
  use water_isotopes, only: pwtspec, ispundef, isph2o, isphdo, isph218o
  use water_types,    only: pwtype, iwtundef, iwtvap, iwtliq, iwtice, iwtstrain, iwtstsnow, iwtcvrain, iwtcvsnow
  use water_tracer_vars

  implicit none
  
  private
  save

!------------------------ Module Interfaces -----------------------------
!
! Public interfaces
!
! CAM physics interface
  public :: wtrc_implements_cnst  ! checks tracer implementation
  public :: wtrc_init             ! initialize module parameters
  public :: wtrc_init_cnst        ! sets data if not on IC file
  public :: wtrc_readnl           ! read the water_trace namelist
  public :: wtrc_register         ! register constituents

! Queries
  public :: wtrc_get_icnst        ! gents cnst index (into wtrc tables) from tracer name
  public :: wtrc_is_ice           ! logical function for m = ice
  public :: wtrc_is_liq           ! logical function for m = liquid
  public :: wtrc_is_cvrain        ! logical function for m = convective rain
  public :: wtrc_is_cvsnow        ! logical function for m = convective snow
  public :: wtrc_is_strain        ! logical function for m = stratiform rain
  public :: wtrc_is_stsnow        ! logical function for m = stratiform snow
  public :: wtrc_is_tagged        ! logical function for m is tagged
  public :: wtrc_is_vap           ! logical function for m = vapour
  public :: wtrc_is_wtrc          ! logical function for m = water tracer

! Water tracer processing
  public :: wtrc_add_rate         ! add process rate to rates matrix
  public :: wtrc_add_rates        ! add process rate to rates matrix
  public :: wtrc_apply_rates      ! calculates isotope tendencies
  public :: wtrc_sediment         ! calculates tendencies due to condensate sedimentation
  public :: wtrc_get_alpha        ! returns the fractionation factor
  public :: wtrc_get_rstd         ! returns the standard isotopic ratios of water vapor
  public :: wtrc_clear_precip     ! sediment cloud condensate to surface as precipitation
  public :: wtrc_collect_precip   ! sediment cloud condensate to surface as precipitation
  public :: wtrc_diagnose_bulk_precip ! diagnose surface precipitation from bulk fields
  public :: wtrc_diagnose_precip  ! diagnose surface precipitation
  public :: wtrc_output_precip    ! output precip fields to history file
  public :: wtrc_init_rates       ! initialize process rate matrix
  public :: wtrc_ratio            ! calulates ratio to precision
  public :: wtrc_ratio_all        ! calulates ratio for all tracers
  public :: wtrc_liqvap_equil     ! equilibrates liquid and vapour <-Might need to be in csm_share
  public :: wtrc_vap_distil       ! distils vapour for some increment
  public :: wtrc_dicm             ! core fractionation routine (iterative)
  public :: wtrc_precip_evap      ! Calculates (convective) precipitation evaporation tendencies
  public :: wtrc_efac             ! equilibrium implicit factor
  public :: wtrc_q1q2_pjr         ! calculates water tracer tendency from ZM deep scheme
  public :: wtrc_mass_fixer       ! resets the Standard (H2O) tracer to Q, and adjusts others accordingly
  public :: wtrc_equil_time       ! calculates the fraction of rain that has experienced equilibration
!  public :: wtrc_rescale          ! scaler routine

  ! Diagnostics
  public :: wtrc_setup_diag       ! write tracer configuration 
  public :: wtrc_check_h2o        ! check for water mass conservation 
  public :: wtrc_adjust_h2o       ! adjust water values to avoid negative and/or small numbers generated by the numerics 

  public :: wtrc_check            ! checks tracer with prognostic
  public :: wtrc_chkdelta         ! checks delta values
  public :: wtrc_init_qpert       ! initialize boundary layer perturbation
  public :: wtrc_qchk1            ! compare 1d tracer with prognostic
  public :: wtrc_qchk2            ! compare 2d tracer with prognostic
  public :: wtrc_qchk3            ! compare all 2d with prognostsic


!------------------- Module Variable Declarations -----------------------
!Moved to water_tracer_vars.F90
!-----------------------------------------------------------------------
contains

!=======================================================================
subroutine wtrc_readnl(nlfile)
!-----------------------------------------------------------------------
!
! Purpose: Read in setting from the namelist file that control the operation
!          of this module.
!
! Method:
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------

    ! Read water tracer namelist group.

    use abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use spmd_utils,      only: masterproc

    ! args

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars

    integer :: unitn, ierr

    ! read namelist for water tracers
    namelist /water_tracer_nl/ &
      water_tracer_model, trace_water, wisotope, wtrc_lh2oadj, wtrc_warn_only, &
      wtrc_qchkmin, wtrc_add_cvprecip, wtrc_add_stprecip, wtrc_names, wtrc_species_names, &
      wtrc_type_names, wtrc_srfvap_names, wtrc_srfpcp_names, &
      wtrc_tag_names, wtrc_alpha_kinetic, wtrc_check_total_h2o, &
      wtrc_detrain_in_macrop, wtrc_niter, wtrc_citer, wtrc_qmin, wtrc_check_show_types, &
      wtrc_fixed_alpha, wtrc_fixed_rstd, wtrc_lzmlin


     if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'water_tracer_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, water_tracer_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('water_tracer_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
     end if

#ifdef SPMD
    call mpibcast(trace_water,            1,                    mpilog,  0, mpicom)
    call mpibcast(wisotope,               1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_lh2oadj,           1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_warn_only,         1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_qchkmin,           1,                    mpir8,   0, mpicom)
    call mpibcast(wtrc_qmin,              1,                    mpir8,   0, mpicom)
    call mpibcast(wtrc_fixed_alpha,       pwtspec,              mpir8,   0, mpicom)
    call mpibcast(wtrc_fixed_rstd,        pwtspec,              mpir8,   0, mpicom)
    call mpibcast(wtrc_add_cvprecip,      1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_add_stprecip,      1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_alpha_kinetic,     1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_check_total_h2o,   1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_check_show_types,  1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_detrain_in_macrop, 1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_lzmlin,            1,                    mpilog,  0, mpicom)
    call mpibcast(wtrc_niter,             1,                    mpiint,  0, mpicom)
    call mpibcast(wtrc_citer,             1,                    mpiint,  0, mpicom)
    call mpibcast(water_tracer_model,     len(water_tracer_model),                  mpichar, 0, mpicom)
    call mpibcast(wtrc_names,             len(wtrc_names(1))*WTRC_MAX_CNST,         mpichar, 0, mpicom)
    call mpibcast(wtrc_species_names,     len(wtrc_species_names(1))*WTRC_MAX_CNST, mpichar, 0, mpicom)
    call mpibcast(wtrc_type_names,        len(wtrc_type_names(1))*WTRC_MAX_CNST,    mpichar, 0, mpicom)
    call mpibcast(wtrc_srfvap_names,      len(wtrc_srfvap_names(1))*WTRC_MAX_CNST,  mpichar, 0, mpicom)
    call mpibcast(wtrc_srfpcp_names,      len(wtrc_srfpcp_names(1))*WTRC_MAX_CNST,  mpichar, 0, mpicom)
    call mpibcast(wtrc_tag_names,         len(wtrc_tag_names(1))*WTRC_MAX_CNST,     mpichar, 0, mpicom)
#endif

end subroutine wtrc_readnl

!=======================================================================
  subroutine wtrc_init
!-----------------------------------------------------------------------
!
! Purpose: initialize water_tracer parameterizations and indexing
!          (declare additional history field)
!
! Method:
!   Initialize modules related to water tracers (water_isotopes and
!   water_types) and add output fields for water tracer state variables.
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 18:01:52 MDT 2003
! 2/5/2012 - bardeenc : Moved indexing to register and added definition
!            of output fields.
!
!-----------------------------------------------------------------------
  use cam_history,    only: addfld, add_default, phys_decomp
  use water_isotopes, only: wiso_init, wiso_get_ispec
  use water_types,    only: wtype_init, wtype_get_itype
  use constituents,   only: cnst_name

  integer            :: icnst
  integer            :: iwset
  character(len=128) :: longname
!-----------------------------------------------------------------------

    if (water_tracer_model /= "none") then

      ! Initialize the water type module
      call wtype_init

      ! If necessary, initialize isotope module
      if (wisotope) then
        call wiso_init
      end if

      write(iulog,*) 'WTRC_INIT: Initializing water tracers.'

      ! Add diagnostic fields for the isotope constituents, and make them
      ! default fields.
      do icnst = 1, wtrc_ncnst
        longname = trim(wtrc_species_names(icnst)) // ' isotope mmr for ' // trim(wtrc_type_names(icnst))

        call addfld(trim(wtrc_names(icnst)), 'kg/kg ', pver, 'A', longname, phys_decomp)
        call add_default(trim(wtrc_names(icnst)), 1, ' ')
      end do

      ! Currently, rain and snow are not prognostic variables, so add then
      ! to the state structure if necessary.
      if (wtrc_add_cvprecip) then
        longname = ' Pseudo-prognostic bulk mmr for ' // trim(wtrc_type_names(iwtcvrain)) // ' - temporary'
        call addfld(trim(wtrc_bulk_names(iwtcvrain)), 'kg/kg ', pver, 'A', longname, phys_decomp)
        call add_default(trim(wtrc_bulk_names(iwtcvrain)), 1, ' ')

        longname = ' Pseudo-prognostic bulk mmr for ' // trim(wtrc_type_names(iwtcvsnow)) // ' - temporary'
        call addfld(trim(wtrc_bulk_names(iwtcvsnow)), 'kg/kg ', pver, 'A', longname, phys_decomp)
        call add_default(trim(wtrc_bulk_names(iwtcvsnow)), 1, ' ')
      end if

      if (wtrc_add_stprecip) then
        longname = ' Pseudo-prognostic bulk mmr for ' // trim(wtrc_type_names(iwtstrain)) // ' - temporary'
        call addfld(trim(wtrc_bulk_names(iwtstrain)), 'kg/kg ', pver, 'A', longname, phys_decomp)
        call add_default(trim(wtrc_bulk_names(iwtstrain)), 1, ' ')

        longname = ' Pseudo-prognostic bulk mmr for ' // trim(wtrc_type_names(iwtstsnow)) // ' - temporary'
        call addfld(trim(wtrc_bulk_names(iwtstsnow)), 'kg/kg ', pver, 'A', longname, phys_decomp)
        call add_default(trim(wtrc_bulk_names(iwtstsnow)), 1, ' ')
      end if

      ! Add output fields for the surface precipitation
      do iwset = 1, wtrc_nwset
        longname = 'Convective rain rate for ' // trim(cnst_name(wtrc_iawset(iwtcvrain, iwset)))
        call addfld('PRECRC_' // trim(cnst_name(wtrc_iawset(iwtcvrain, iwset))), 'm/s', 1, 'A', longname, phys_decomp)
        call add_default('PRECRC_' // trim(cnst_name(wtrc_iawset(iwtcvrain, iwset))), 1, ' ')

        longname = 'Convective snow rate (water equivalent) for ' // trim(cnst_name(wtrc_iawset(iwtcvsnow, iwset)))
        call addfld('PRECSC_' // trim(cnst_name(wtrc_iawset(iwtcvsnow, iwset))), 'm/s', 1, 'A', longname, phys_decomp)
        call add_default('PRECSC_' // trim(cnst_name(wtrc_iawset(iwtcvsnow, iwset))), 1, ' ')

        longname = 'Large-scale (stable) rain rate for ' // trim(cnst_name(wtrc_iawset(iwtstrain, iwset)))
        call addfld('PRECRL_' // trim(cnst_name(wtrc_iawset(iwtstrain, iwset))), 'm/s', 1, 'A', longname, phys_decomp)
        call add_default('PRECRL_' // trim(cnst_name(wtrc_iawset(iwtstrain, iwset))), 1, ' ')

        longname = 'Large-scale (stable) snow rate (water equivalent) for ' // trim(cnst_name(wtrc_iawset(iwtstsnow, iwset)))
        call addfld('PRECSL_' // trim(cnst_name(wtrc_iawset(iwtstsnow, iwset))), 'm/s', 1, 'A', longname, phys_decomp)
        call add_default('PRECSL_' // trim(cnst_name(wtrc_iawset(iwtstsnow, iwset))), 1, ' ')
      end do

      write(iulog,*) 'WTRC_INIT: done.'
      write(iulog,*) ''
    end if

    return
end subroutine wtrc_init

!=======================================================================
  subroutine wtrc_register
!-----------------------------------------------------------------------
!
! Purpose: register advected water tracer constituents
!
! Method:
!  Calls CAM constituent registration routines based on 
!  water tracer species and phase indexing. Create a number of tracer
!  lists, and do some consistency checking of the tracers variables
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 15:31:56 MDT 2003
! Modified for CAM5 by Jesse Nusbaumer <nusbaume@colorado.edu> - 3/30/2011
! 2/3/2012 - bardeenc : Refactored, extended for rain and snow.
!
!-----------------------------------------------------------------------
    use constituents,   only: cnst_get_ind, cnst_add
    use water_isotopes, only: wiso_get_ispec
    use water_types,    only: wtype_get_itype
    use physconst,      only: cpair, cpwv, mwh2o, mwdry
    use constituents,   only: cnst_name
    use physics_buffer, only: pbuf_add_field, dtype_r8

    real(r8)             :: mw  ! molecular weight for the constituent
    real(r8)             :: cp  ! specific heat for the constituent

    integer              :: itmpndx, icnst, itag, iwset, pidx

    character(len=128)     :: longname

   !For wtrc_out_names -JN:
    character(len=8) :: wtstrng = "" !variable name
    integer          :: wtstsz       !name length

!-----------------------------------------------------------------------
!
! Initialize all tracers as non-water, with unknown species
!
    iwater(:)  = iwtundef
    iwspec(:)  = ispundef
    iwistag(:) = .false.

  ! Set the species of the total water as H2O, but DONT set them 
  ! as water tracers, as they are prognostic (not "tracers")
  !
  ! Set names of variable tendencies and declare them as history variables
    if (water_tracer_model /= "none") then

      write(iulog,*) ''
      write(iulog,*) 'WTRC_REGISTER: starting ...'

      ! Count up the number of water tracers
      icnst = 1
      do while(wtrc_names(icnst) /= "")
        icnst = icnst + 1
      end do
      wtrc_ncnst  = icnst-1
      wtrc_nwset  = wtrc_ncnst / pwtype

      ! The number of constituents should be a perfect multiple of the number
      ! of water tracers.
      if (mod(wtrc_ncnst, pwtype) /= 0) then
        call endrun("WTRC_INIT: Error - The number of water traces should be an even multiple of the number of water types.")
      end if

      if (wtrc_ncnst <= 0) then
        call endrun("WTRC_INIT: Error - Water tracers have been enabled, but no tracers have been specified.")
      end if

      if (wtrc_ncnst > WTRC_MAX_CNST) then
        call endrun("WTRC_INIT: Error - Too many water tracers have been enabled, increase WTRC_MAX_CNST and recompile.")
      end if

      ! Look up all of the species and types
      do icnst = 1, wtrc_ncnst
        wtrc_species(icnst) = wiso_get_ispec(trim(wtrc_species_names(icnst)))
        if (wtrc_species(icnst) == ispundef) then
          call endrun('WTRC_INIT: Error - Tracer species id not found for ' // trim(wtrc_species_names(icnst)))
        end if

        wtrc_types(icnst)   = wtype_get_itype(trim(wtrc_type_names(icnst)))
        if (wtrc_types(icnst) == iwtundef) then
          call endrun('WTRC_INIT: Error - Tracer type id not found for ' // trim(wtrc_type_names(icnst)))
        end if
      end do

      ! Currently, rain and snow are not prognostic variables, so add then
      ! to the state structure if necessary.
      if (wtrc_add_cvprecip) then
        longname = ' Pseudo-prognostic bulk mmr for ' // trim(wtrc_type_names(iwtcvrain)) // ' - temporary'
        call cnst_add(trim(wtrc_bulk_names(iwtcvrain)), mwdry, cpair, 0._r8, itmpndx, longname, readiv=.false.)

        longname = ' Pseudo-prognostic bulk mmr for ' // trim(wtrc_type_names(iwtcvsnow)) // ' - temporary'
        call cnst_add(trim(wtrc_bulk_names(iwtcvsnow)), mwdry, cpair, 0._r8, itmpndx, longname, readiv=.false.)
      end if

      if (wtrc_add_stprecip) then
        longname = ' Pseudo-prognostic bulk mmr for ' // trim(wtrc_type_names(iwtstrain)) // ' - temporary'
        call cnst_add(trim(wtrc_bulk_names(iwtstrain)), mwdry, cpair, 0._r8, itmpndx, longname, readiv=.false.)

        longname = ' Pseudo-prognostic bulk mmr for ' // trim(wtrc_type_names(iwtstsnow)) // ' - temporary'
        call cnst_add(trim(wtrc_bulk_names(iwtstsnow)), mwdry, cpair, 0._r8, itmpndx, longname, readiv=.false.)
      end if

      ! Make a list of the bulk (parent) constituents for the water types.
      ! It might be better if iwspec and iwater were added as variables in constituent.F90
      ! and the calls that create these constituents updated to set the variables correctly.
      do icnst = 1, pwtype
        call cnst_get_ind(wtrc_bulk_names(icnst), wtrc_bulk_indices(icnst))
        iwspec(wtrc_bulk_indices(icnst)) = isph2o
      end do

      ! Make a list of the bulk (parent) constituents for the water types.
      ! It might be better if iwspec and iwater were added as variables in
      ! constituent.F90
      ! and the calls that create these constituents updated to set the
      ! variables correctly.
      do icnst = 1, pwtype
        call cnst_get_ind(wtrc_bulk_names(icnst), wtrc_bulk_indices(icnst))
        iwspec(wtrc_bulk_indices(icnst)) = isph2o
      end do

      ! Is this constituent a tag or an isotopologue?
      wtrc_is_tag(:) = .false.
      do icnst = 1,wtrc_ncnst
        do itag = 1,wtrc_ncnst
          if (wtrc_tag_names(itag) == wtrc_names(icnst)) then
            wtrc_is_tag(icnst) = .true.
            exit
          end if
        end do
      end do

      ! Add constituents for each of the isotope tracers.
      wtrc_nspec(:) = 0
      wtrc_ntype(:) = 0

      do icnst = 1, wtrc_ncnst

        ! The weight of the condensate is too large for the standard diffusion
        ! routines, so treat them as just part of the air parcel and adjust
        ! later
        ! in the microphysics.
        if (wtrc_types(icnst) == iwtvap) then
          mw = mwh2o
          cp = cpwv
        else
          mw = mwdry
          cp = cpair
        end if

        longname = trim(wtrc_species_names(icnst)) // ' isotope mmr for ' // trim(wtrc_type_names(icnst))

        if (wtrc_is_tag(icnst)) then
          longname = trim(longname) // ' (tagged)'
        end if

        ! Register the water isotope tracer.          
        call wtrc_cnst_add(wtrc_names(icnst), wtrc_types(icnst), wtrc_species(icnst), &
          wtrc_is_tag(icnst), mw, cp, 0._r8, wtrc_indices(icnst), longname, readiv=.false.)

        ! Make some indices and counts to make it easier to loop over the isotopes.
        wtrc_ntype(wtrc_types(icnst)) = wtrc_ntype(wtrc_types(icnst)) + 1
        wtrc_iatype(wtrc_ntype(wtrc_types(icnst)), wtrc_types(icnst)) = wtrc_indices(icnst)

        wtrc_nspec(wtrc_species(icnst)) = wtrc_nspec(wtrc_species(icnst)) + 1
!        wtrc_iaspec(wtrc_nspec(wtrc_types(icnst)), wtrc_species(icnst)) = wtrc_indices(icnst)
        wtrc_iaspec(wtrc_nspec(wtrc_species(icnst)), wtrc_species(icnst)) = wtrc_indices(icnst)

        ! The constituents are assumed to be ordered by wset.
        wtrc_iawset(wtrc_types(icnst), ((icnst - 1) / pwtype) + 1) = wtrc_indices(icnst)

        ! Make sure they are all the same species.
        if (mod((icnst - 1), pwtype) /= 0) then
          if (wtrc_species(icnst) /= wtrc_species(((icnst - 1) / pwtype) * pwtype + 1)) then
            call endrun("WTRC_INIT: Error - Water tracers must be ordered by wsets of the same isotopic species.")
          end if
        end if
      end do

      ! Determine which specifies are coupled to the surface.      
      icnst = 1
      do while(wtrc_srfvap_names(icnst) /= "")
        wtrc_iasrfvap(icnst) = wtrc_get_icnst(wtrc_srfvap_names(icnst))
        if (wtrc_iasrfvap(icnst) == -1) then
          call endrun('WTRC_INIT: Error - Surface vapor tracer not found for ' // trim(wtrc_srfvap_names(icnst)))
        end if
        icnst = icnst + 1
      end do
      wtrc_nsrfvap = icnst-1

      ! Set output names:
      do icnst=1, wtrc_nsrfvap
        wtstrng = wtrc_srfvap_names(icnst)          !copy string
        wtstsz  = scan(wtstrng,'V',BACK=.true.)     !determine size of string (to cut off trailing "V").
        wtrc_out_names(icnst) = wtstrng(1:wtstsz-1) !Remove trailing V.
      end do

      icnst = 1
      do while(wtrc_srfpcp_names(icnst) /= "")
        wtrc_iasrfpcp(icnst) = wtrc_get_icnst(wtrc_srfpcp_names(icnst))
        if (wtrc_iasrfpcp(icnst) == -1) then
          call endrun('WTRC_INIT: Error - Surface precipitation tracer not found for ' // trim(wtrc_srfpcp_names(icnst)))
        end if
        icnst = icnst + 1
      end do
      wtrc_nsrfpcp = icnst-1

      ! Add 2D surface precipitation variables in the physics buffer for each of the 
      ! precipitation types and wsets.
      wtrc_srfpcp_indices(:,:) = -1

      do iwset = 1, wtrc_nwset
        call pbuf_add_field('P_' // trim(cnst_name(wtrc_iawset(iwtcvrain, iwset))), 'physpkg', dtype_r8, (/pcols/), pidx)
        wtrc_srfpcp_indices(iwtcvrain,iwset) = pidx

        call pbuf_add_field('P_' // trim(cnst_name(wtrc_iawset(iwtcvsnow, iwset))), 'physpkg', dtype_r8, (/pcols/), pidx)
        wtrc_srfpcp_indices(iwtcvsnow,iwset) = pidx

        call pbuf_add_field('P_' // trim(cnst_name(wtrc_iawset(iwtstrain, iwset))), 'physpkg', dtype_r8, (/pcols/), pidx)
        wtrc_srfpcp_indices(iwtstrain,iwset) = pidx

        call pbuf_add_field('P_' // trim(cnst_name(wtrc_iawset(iwtstsnow, iwset))), 'physpkg', dtype_r8, (/pcols/), pidx)
        wtrc_srfpcp_indices(iwtstsnow,iwset) = pidx
      end do

      ! Once the registration is done, report what we actually have just to make sure
      call wtrc_setup_diag

      write(iulog,*) 'WTRC_REGISTER: done.'
      write(iulog,*) ''
    end if


    return
end subroutine wtrc_register

!=======================================================================
  function wtrc_is_wtrc(m)
!-----------------------------------------------------------------------
! Returns true if tracer is water tracer
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_wtrc
!-----------------------------------------------------------------------
    wtrc_is_wtrc = .false.
    if (iwater(m) /= iwtundef) wtrc_is_wtrc = .true.
  return
  end function wtrc_is_wtrc

!=======================================================================
  function wtrc_is_vap(m)
!-----------------------------------------------------------------------
! Returns true if tracer is vapour
!-----------------------------------------------------------------------
  integer, intent(in) :: m    ! constituent index
  logical wtrc_is_vap
!-----------------------------------------------------------------------
    wtrc_is_vap = .false.
    if (iwater(m) == iwtvap) wtrc_is_vap = .true.
  return
  end function wtrc_is_vap

!=======================================================================
  function wtrc_is_liq(m)
!-----------------------------------------------------------------------
! Returns true if tracer is cloud liquid
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_liq
!-----------------------------------------------------------------------
    wtrc_is_liq = .false.
    if (iwater(m) == iwtliq) wtrc_is_liq = .true.
  return
  end function wtrc_is_liq

!=======================================================================
  function wtrc_is_ice(m)
!-----------------------------------------------------------------------
! Returns true if tracer is cloud ice
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_ice
!-----------------------------------------------------------------------
    wtrc_is_ice = .false.
    if (iwater(m) == iwtice) wtrc_is_ice = .true.
  return
  end function wtrc_is_ice

!=======================================================================
  function wtrc_is_cvrain(m)
!-----------------------------------------------------------------------
! Returns true if tracer is convective rain
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_cvrain
!-----------------------------------------------------------------------
    wtrc_is_cvrain = .false.
    if (iwater(m) == iwtcvrain) wtrc_is_cvrain = .true.
  return
  end function wtrc_is_cvrain

!=======================================================================
  function wtrc_is_strain(m)
!-----------------------------------------------------------------------
! Returns true if tracer is stratiform rain
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_strain
!-----------------------------------------------------------------------
    wtrc_is_strain = .false.
    if (iwater(m) == iwtstrain) wtrc_is_strain = .true.
  return
  end function wtrc_is_strain

!=======================================================================
  function wtrc_is_cvsnow(m)
!-----------------------------------------------------------------------
! Returns true if tracer is convective snow
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_cvsnow
!-----------------------------------------------------------------------
    wtrc_is_cvsnow = .false.
    if (iwater(m) == iwtcvsnow) wtrc_is_cvsnow = .true.
  return
  end function wtrc_is_cvsnow

!=======================================================================
  function wtrc_is_stsnow(m)
!-----------------------------------------------------------------------
! Returns true if tracer is stratiform snow
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_stsnow
!-----------------------------------------------------------------------
    wtrc_is_stsnow = .false.
    if (iwater(m) == iwtstsnow) wtrc_is_stsnow = .true.
  return
  end function wtrc_is_stsnow

!=======================================================================
  function wtrc_is_tagged(m)
!-----------------------------------------------------------------------
! Returns true if tracer is tagged water
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_tagged
!-----------------------------------------------------------------------
    wtrc_is_tagged = iwistag(m)
  return
  end function wtrc_is_tagged

!=======================================================================
  function wtrc_get_icnst(name)
!-----------------------------------------------------------------------
! Purpose: Retrieve wtrc constituent index, based on tracer name
! Author: Chuck Bardeen
!-----------------------------------------------------------------------
    character(len=*),  intent(in)  :: name              ! tracer name
    integer                        :: wtrc_get_icnst    ! return constituent index
!-----------------------------------------------------------------------
    do wtrc_get_icnst = 1, wtrc_ncnst
      if (name == wtrc_names(wtrc_get_icnst)) then
        return
      end if
    end do
    
    wtrc_get_icnst = -1
    
    return
  end function wtrc_get_icnst
  

!=======================================================================
  subroutine wtrc_set_ptendlq(ptend)
!-----------------------------------------------------------------------
! Purpose: Retrieve wtrc constituent index, based on tracer name
! Author: Chuck Bardeen
!-----------------------------------------------------------------------
    use physics_types,  only: physics_ptend

    type(physics_ptend), intent(inout) :: ptend    ! State tendencies

    integer                            :: icnst    ! constituent index
!-----------------------------------------------------------------------
    do icnst = 1, wtrc_ncnst
      ptend%lq(wtrc_indices(icnst)) = .true.
    end do
    
    return
  end subroutine wtrc_set_ptendlq
  

!=======================================================================
  subroutine wtrc_cnst_add(name, iwt, isp, is_tag, mwc, cpc, qminc, ind, &
                           longname, readiv, mixtype)
!-----------------------------------------------------------------------
! Purpose: provide a wrapper for cnst_add with added index condifuration
!          for more details registration of water tracers
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 21:02:25 MDT 2003
!-----------------------------------------------------------------------
  use constituents, only: cnst_add 
!---------------------------- Arguments --------------------------------
    character(len=*), intent(in) :: &
       name      ! constituent name for variable name in history file(8 char max)
    character(len=*), intent(in), optional :: &
       longname  ! long_name attribute in netcdf output (128 char max) [name]
    logical,          intent(in), optional :: &
       readiv    ! true => read initial values from initial file (default: true)
    character(len=*),         intent(in), optional :: &
       mixtype    ! mixing ratio type (dry, wet)

!    integer, intent(in)    :: type   ! flag indicating advected or nonadvected
    integer, intent(in)    :: iwt    ! water type indicator
    integer, intent(in)    :: isp    ! water species indicator
    logical, intent(in)    :: is_tag ! tagged tracer indicator
    real(r8),intent(in)    :: mwc    ! const. molecular weight (kg/kmol)
    real(r8),intent(in)    :: cpc    ! const. spcfic heat  const press (J/kg/K)
    real(r8),intent(in)    :: qminc  ! minimum  mass mixing ratio (kg/kg)
!                                        normally 0., except water 1.E-12, for
!                                        radiation.

    integer, intent(out)   :: ind    ! global constituent index (in q array)

!-----------------------------------------------------------------------
!
! Pass arguments on to normal code
!
    call cnst_add(name, mwc, cpc, qminc, ind, longname, readiv, mixtype)
!
! Knowing the tracer index assign water type and species
!
    iwater(ind)  = iwt
    iwspec(ind)  = isp
    iwistag(ind) = is_tag
!
    return
  end subroutine wtrc_cnst_add


!=======================================================================
  function wtrc_implements_cnst(name)
!-----------------------------------------------------------------------
!
! Purpose: return true if specified constituent is implemented by this package
! Notice wtrc_names should be the same as hard coded calls to cnst_resister.
!  
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 16:10:29 MDT 2003
! 2/3/2012 - bardeenc : Refactored, put in stratiform interfaces and extended for rain and snow.
!
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------
     character(len=*), intent(in) :: name   ! constituent name
     logical :: wtrc_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: icnst
!-----------------------------------------------------------------------
    wtrc_implements_cnst = .false.
    
    do icnst = 1, wtrc_ncnst
      if (name == wtrc_names(icnst)) then
        wtrc_implements_cnst = .true.
        return
      end if
    end do
    
    return
  end function wtrc_implements_cnst


!=======================================================================
  subroutine wtrc_init_cnst(name, q, gcid, qq, ql, qi)
!-----------------------------------------------------------------------
!
! Initializes water tracers if not read from initial conditions file.
! Assign as some standard mass fraction of the prognostic waters
! (which  can assumed are set correctly). If using water isotope
! Set the standard ration, else set zero
!
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 18:24:57 MDT 2003
! 2/3/2012 - bardeenc : Changed to a way compatible with CAM5's initialization
! interface; however, it requires some special coding in initdat.F90.
!
!-----------------------------------------------------------------------
    use constituents, only: pcnst, cnst_get_ind
!---------------------------- Arguments --------------------------------
    character(len=*),intent(in)  :: name               ! tracer name
    real(r8), intent(out)        :: q(:,:)             ! mass mixing ratio (gcol, lev)
    integer, intent(in)          :: gcid(:)            ! global column id
    real(r8), intent(in)         :: qq(:,:)            ! Q mass mixing ratio (gcol, lev)
    real(r8), intent(in)         :: ql(:,:)            ! CLDLIQ mass mixing ratio (gcol, lev)
    real(r8), intent(in)         :: qi(:,:)            ! CLDICE mass mixing ratio (gcol, lev)
 !------------------------- Local Variables -----------------------------
    integer ixwtrc              ! index of water tracer
    integer ixwprg              ! intext of water prognostic
    real(r8) rat                ! an isotope ratio
!-----------------------------------------------------------------------
!

    ! Retrieve the tracer index, and work out index of equivilent prognostic
    call cnst_get_ind(name, ixwtrc)      ! this SHOULD be m in calling routine
    if (.not. wtrc_is_wtrc(ixwtrc)) then
      call endrun( 'WTRC_INIT_CNST: non water tracer detected.')
    else
   
      ! Assign tracer to be total, scaled by some standard ratio 
      if (trace_water) then
      
        !standard ratio or isotope run:     
        rat = wtrc_get_rstd(iwspec(ixwtrc))
        !----------------------
        if(wtrc_is_tagged(ixwtrc)) then !not standard (H2O) water
          rat = 0._r8
        end if
        !----------------------       
        if (iwater(ixwtrc) == iwtvap) then
          q(:,:) = rat * qq(:,:)
        else if (iwater(ixwtrc) == iwtliq) then
          q(:,:) = rat * ql(:,:)
        else if (iwater(ixwtrc) == iwtice) then
          q(:,:) = rat * qi(:,:)
        else
          q(:,:) = 0._r8
        end if
      else
        q(:,:) = 0._r8
      end if
    end if
    
    return
  end subroutine wtrc_init_cnst


!=======================================================================
  subroutine wtrc_add_rate(process_rates, icol, iz, isrctype, idsttype, rtype, rate, do_reverse)
!-----------------------------------------------------------------------
!
! Purpose: Add process rate to matrix
!
! Method:
!   Add a process rate for the production of idsttype from isrctype
!   to a matrix of process rates. This routine will add the reverse
!   rate for the loss term for isrctype to idsttype.
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use water_types,    only: pwtype
  
  real(r8), intent(inout)         :: process_rates(pcols,pver,pwtype,pwtype,pwtype) ! process rates (kg/kg/sec)
  integer, intent(in)             :: icol                                           ! column index
  integer, intent(in)             :: iz                                             ! vertical index
  integer, intent(in)             :: isrctype                                       ! source water type index
  integer, intent(in)             :: idsttype                                       ! destination water type index
  integer, intent(in)             :: rtype                                          ! water type used to calculate ratio - JN
  real(r8), intent(in)            :: rate                                           ! process rate  (kg/kg/s) for src -> dst
  logical, intent(in), optional   :: do_reverse                                     ! apply reverse rates too (default=true) 
    
  logical                         :: ldo_reverse
!-----------------------------------------------------------------------
    
    ! Determine whether to apply the reverse rates.
    if (present(do_reverse)) then
      ldo_reverse = do_reverse
    else
      ldo_reverse = .true.
    end if
    
    process_rates(icol, iz, idsttype, isrctype, rtype) = process_rates(icol, iz, idsttype, isrctype, rtype) + rate

    if ((isrctype /= idsttype) .and. (ldo_reverse)) then
      process_rates(icol, iz, isrctype, idsttype, rtype) = process_rates(icol, iz, isrctype, idsttype, rtype) - rate
    end if
    
    return
  end subroutine wtrc_add_rate


!=======================================================================
  subroutine wtrc_add_rates(process_rates, ncol, isrctype, idsttype, rtype, rate, do_reverse)
!-----------------------------------------------------------------------
!
! Purpose: Add process rates to matrix
!
! Method:
!   Add process rates for the production of idsttype from isrctype
!   to a matrix of process rates. This routine will also add the
!   rate for the loss term for  for isrctype to idsttype.
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use water_types,    only: pwtype
  
  real(r8), intent(inout)         :: process_rates(pcols,pver,pwtype,pwtype,pwtype)  ! process rates (kg/kg/sec)
  integer, intent(in)             :: ncol                                     ! number of columns
  integer, intent(in)             :: isrctype                                 ! source water type index
  integer, intent(in)             :: idsttype                                 ! destination water type index
  integer, intent(in)             :: rtype                                    ! water type used to calculate Ratio - JN
  real(r8), intent(in)            :: rate(pcols,pver)                         ! process rate  (kg/kg/s) for src -> dst
  logical, intent(in), optional   :: do_reverse                               ! apply reverse rates too (default=true) 

  integer                         :: k
  integer                         :: icol
!-----------------------------------------------------------------------
    
    do k = 1, pver
      do icol = 1, ncol
        call wtrc_add_rate(process_rates, icol, k, isrctype, idsttype, rtype, rate(icol,k), do_reverse)
      end do
    end do
    
    return
  end subroutine wtrc_add_rates


!=======================================================================
  subroutine wtrc_init_rates(process_rates)
!-----------------------------------------------------------------------
!
! Purpose: Initialize process rate matrix
!
! Method:
!   Initializes the matrix of process rates.
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use water_types,    only: pwtype
  
  real(r8), intent(out)   :: process_rates(pcols,pver,pwtype,pwtype,pwtype)   ! Process rates (kg/kg/sec)
  
!-----------------------------------------------------------------------
    
    process_rates(:, :, :, :, :) = 0._r8
    
    return
  end subroutine wtrc_init_rates


!=======================================================================
  subroutine wtrc_apply_rates(pstate, ptend_sum, pbuf, dtime, micro, pre_rates, sed_rates, post_rates, do_stprecip, &
                              liqcldf, icecldf, fc, fi, prelat, postlat, frzro, meltso)
!-----------------------------------------------------------------------
!
! Purpose: Calculate water tracer tendencies from process rates
!
! Method:
!   Given the microphysics process rates (pre-sedimentation, sedimentation
!   and/or post-sedimentation), calculate the effect of these rater on the
!   mass mixing ratio of water tracers.
!
!   An iterative technique is used with pre, sed and post rates applied
!   during each iteration.
!
!   NOTE: The physics state state supplied should be the same one
!   used for the bulk water properties calculation that generated the
!   rates, and the tendency should only include the results of the
!   tendencies from those processes.
!
! Author: Chuck Bardeen
!
! NOTE:  I have commented out the adjust_H2O calls in this routine because
! they seemed to produce errors.  However, the testing was done with niter
! set to 1.  Thus if looping is turned on, and the values become substantially
! worse, then the adjust_H2O subroutines might need to be reinstated. - JN
! NOTE:  Adding addtional iterations appeared to make no difference.  Still,
! the adjust_H2O routine calls will remain just in case, at least until the code is 
! cleaned-up. - JN
!
! NOTE: The subroutine here makes the assumption that the tendency rates are a
! linear function of the state, e.g.: dq/dt = f(q_i) = f(Rq) = R*f(q).
! However, in the actual microphysics routine, the relationship is often-times
! non-linear, and thus this subroutine may not be able to simulate those cloud
! physical tendencies with full accuracy.  If it is found that the stratiform
! isotope values are not as they should be, this could be the culprit. - JN
!
! NOTE:  The order of operations has been changed such that each phase change happens
! sequentially, instead of all at once.  Although this is not as physically
! realistic, it does eliminate an error in the isotopic stratiform precipitation
! ,and it allows one to use substantially less code.  However, if its found that
! accurate isotope values simply cannot be produced with this subroutine, then the
! modified order of operations may be to blame. - JN
!
!-----------------------------------------------------------------------
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init, &
                            physics_ptend_sum, physics_ptend_reset
  use water_isotopes, only: pwtspec
  use water_types,    only: pwtype
  use constituents,   only: cnst_name, pcnst
  use physconst,      only: cpair,gravit,rhoh2o
  use physics_buffer, only: physics_buffer_desc

  type(physics_state), intent(in)    :: pstate                                      ! State of the atmosphere
  type(physics_ptend), intent(inout) :: ptend_sum                                   ! State tendencies
  type(physics_buffer_desc), pointer :: pbuf(:)                                     ! physics buffer
  real(r8), intent(in)               :: dtime                                       ! length of timestep (s)
  logical, intent(in)                :: micro                                       ! Microphysics?
  real(r8), intent(in), optional     :: pre_rates(pcols,pver,pwtype,pwtype,pwtype)  ! Pre sedimentation process rates (kg/kg/sec)
  real(r8), intent(in), optional     :: sed_rates(pcols,pver,pwtype)                ! Sedimentation rates (kg/kg/sec)
  real(r8), intent(in), optional     :: post_rates(pcols,pver,pwtype,pwtype,pwtype) ! Post sedimentation process rates (kg/kg/sec)

  logical, intent(in), optional      :: do_stprecip
 
  !Needed for wtrc_sediment:
  real(r8), intent(in), optional     :: liqcldf(pcols,pver) !liquid cloud fraction
  real(r8), intent(in), optional     :: icecldf(pcols,pver) !ice cloud fraction
  real(r8), intent(in), optional     :: fc(pcols,pver)      !initial fall velocity of cloud liquid
  real(r8), intent(in), optional     :: fi(pcols,pver)      !initial fall velocity of cloud ice       

  !latent heating terms:
  real(r8), intent(in), optional     :: prelat(pcols,pver)  !latent heating due to pre-rates
  real(r8), intent(in), optional     :: postlat(pcols,pver) !latent heating due to post-rates
  
  !precipitation phase changes:
  real(r8), intent(in), optional     :: frzro(pcols,pver)   !fraction of rain that freezes
  real(r8), intent(in), optional     :: meltso(pcols,pver)  !fraction of snow that melts

  integer            :: iwset
  integer            :: itype
  integer            :: icol
  integer            :: icnst
  integer            :: iter
  integer            :: i
  integer            :: k
  integer            :: m
  integer            :: mbase
  integer            :: ncol
  integer            :: msrc
  integer            :: isrctype
  integer            :: idsttype
  integer            :: rtype 
  integer            :: mdst
  integer            :: ispec
  logical            :: isOk
  logical            :: ldo_stprecip
  real(r8)           :: R
  real(r8)           :: Rs                       !used for sedimenting qloc
  real(r8)           :: alpha                    !fractionation factor
  real(r8)           :: fr                       !used for Rayleigh fractionation
  real(r8)           :: qloc(pcols,pver,pcnst)
  real(r8)           :: qloc0(pcols,pver,pcnst)
  real(r8)           :: tloc(pcols,pver)

  real(r8)           :: rmass(pcols,wtrc_nwset)  !Vertically-integrated rain
  real(r8)           :: smass(pcols,wtrc_nwset)  !Vertically-integrated snow
  real(r8)           :: rmass0(pcols,wtrc_nwset) !copy of rain
  real(r8)           :: smass0(pcols,wtrc_nwset) !copy of snow 

  real(r8)           :: diff(pcols,pver,pwtype)  !Numerical precision fixer
  real(r8)           :: pdiff                    !precipitation numerical fixer
  real(r8)           :: dliqiso                  !change in liquid due to isotopic equilibration
  real(r8)           :: qtmp                     !temporary variables needed for isotopic equilibration
  real(r8)           :: itmp
  real(r8)           :: rtmp
  real(r8)           :: irtmp

  real(r8)           :: dz                       !layer thickness in meters

  !For partial rain equilibration:
  real(r8)           :: difrm(pwtspec) !should be read from water_isotopes.F90
  real(r8)           :: radius         !raindrop radius in meters
  real(r8)           :: fequil         !fraction of rain that is equilibrated.


!-----------------------------------------------------------------------
!

  if (trace_water) then
  
    if (present(do_stprecip)) then
      ldo_stprecip = do_stprecip
    else
      ldo_stprecip = .false.
    end if
    
    ncol = pstate%ncol

    !Copy initial state:
    qloc(:ncol,:,:)  = pstate%q(:ncol,:,:)
    qloc0(:ncol,:,:) = qloc(:ncol,:,:)
    tloc(:ncol,:)    = pstate%t(:ncol,:)

   !Clear out the stratiform precipitation.
    if (ldo_stprecip) then
      call wtrc_clear_precip(pstate, pbuf, iwtstrain)
      call wtrc_clear_precip(pstate, pbuf, iwtstsnow)
     !Initialize precipitation integrals:
      rmass(:,:)  = 0._r8 
      smass(:,:)  = 0._r8
      rmass0(:,:) = 0._r8
      smass0(:,:) = 0._r8
    end if
    
    ! Check to see if the total water mass in the bulk fields matches the
    ! water traces. For the h216o model, this indicates that bulk fields are
    ! out of sync with the isotopes, indicating that some water process was
    ! missed or perhaps handled incorrectly.
!    isOk = wtrc_check_h2o("WTRC_APPLY_RATES(PRE-" // trim(ptend_sum%name) // ")", pstate, qloc, dtime)

    ! Iterate over each of the isotopologues and calculate tendencies
    ! based upon the sedimentation rates from the bulk prognostic species.

      !***************************     
      !Pre-sedimentation processes
      !***************************
      if(present(pre_rates)) then
        !Do this top down, like is done in MG microphysics.
        do iter=1, wtrc_niter
          do k = 1, pver
            do i = 1, ncol
              !----------------------------------------------------
              !calculate average temperature during cloud processes
              !----------------------------------------------------
              tloc(i,k) = (tloc(i,k) + (tloc(i,k)+prelat(i,k)/cpair*dtime))/2._r8
              !-----------------------------
              do isrctype=1,pwtype !loop through moisture sources
                do idsttype=1,pwtype !loop through water types that are modified
                  rtype = isrctype  !the moisture source is what determines the value of R
                  if(pre_rates(i,k,idsttype,isrctype,rtype) .gt. 0._r8) then !make sure destination type is increasing
                    do iwset=1,wtrc_nwset !loop over water tracers/isotopes

                      msrc = wtrc_iawset(isrctype,iwset) !source water index
                      mbase= wtrc_iawset(isrctype,1)     !H2O source water index
                      mdst = wtrc_iawset(idsttype,iwset) !destination water index

                      if(isrctype .gt. iwtice) then !Is precipitation the source?
                       !calculate ratio:
                        if(isrctype .eq. iwtstrain) then !rain
                          R = wtrc_ratio(iwspec(msrc),rmass0(i,iwset),rmass0(i,1)) 
                        else !snow
                          R = wtrc_ratio(iwspec(msrc),smass0(i,iwset),smass0(i,1))
                        end if 
                      else
                        R = wtrc_ratio(iwspec(msrc),qloc0(i,k,msrc),qloc0(i,k,mbase)) !calculate ratio for single level
                      end if

                      if(wisotope .and. (iwset .ne. 1) .and. (isrctype .eq. iwtvap) .and. &
                         ((idsttype .eq. iwtice) .or. (idsttype .eq. iwtstsnow))) then !only do Rayleigh distillation for water isotopes.

                        ispec = iwspec(mdst) !determine isotope species of ice or snow (and thus vapoor)
                        alpha = wtrc_get_alpha(qloc0(i,k,wtrc_iawset(iwtvap,WTRC_WSET_STD)), tloc(i,k), ispec, isrctype, idsttype, &
                                            .true., pstate%pmid(i,k))

                        !calculate change in vapor mass:
                        fr = qloc(i,k,mbase)/qloc0(i,k,mbase) 
                        !Constrain ratio fraction:
                        fr = min(1._r8,max(0._r8,fr))

                        qloc(i,k,msrc) = qloc0(i,k,msrc)*(fr**alpha)
                        if(idsttype .eq. iwtstsnow) then !snow?
                          smass(i,iwset) = smass(i,iwset) + (qloc0(i,k,msrc)-qloc(i,k,msrc))*pstate%pdel(i,k) 
                        else !ice?
                          qloc(i,k,mdst) = qloc(i,k,mdst) + (qloc0(i,k,msrc)-qloc(i,k,msrc))
                        end if

                      else if(((isrctype .eq. iwtliq) .or. (isrctype .eq. iwtice)) .and. (isrctype .eq. idsttype)) then 
                        !Apply Bergeron process to ice:
                        !-----------------------------
                        if(isrctype .eq. iwtliq) then !liquid type used for Bergeron to ice
                          mdst = wtrc_iawset(iwtice,iwset)
                        else !ice type used for Bergeron to snow
                          msrc = wtrc_iawset(iwtliq,iwset)          !change source term to liquid
                          mbase = wtrc_iawset(iwtliq,WTRC_WSET_STD) !needed for ratio calculation
                          mdst = wtrc_iawset(iwtstsnow,iwset)
                         !Re-calculate source ratio:
                          R = wtrc_ratio(iwspec(msrc),qloc0(i,k,msrc),qloc0(i,k,mbase)) !calculate ratio for single level
                        end if
                        if(wisotope .and. (iwset .gt. 1)) then !Are you doing water isotopes?
                          !Fractionation variables:
                          ispec = iwspec(wtrc_iawset(iwtice,iwset)) 
                   
                         !cloud ice and snow use the same fractionation factors, so just use ice:
                          alpha = wtrc_get_alpha(qloc0(i,k,wtrc_iawset(iwtvap,WTRC_WSET_STD)), tloc(i,k), ispec, iwtvap, iwtice, &
                                            .true., pstate%pmid(i,k))

                          !cloud liquid to cloud vapor:
                          qloc(i,k,wtrc_iawset(iwtvap,iwset)) = qloc(i,k,wtrc_iawset(iwtvap,iwset))+&
                               R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter
                          qloc(i,k,msrc) = qloc(i,k,msrc)-R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter

                          !save variables
                           qloc0(i,k,wtrc_iawset(iwtvap,iwset)) = qloc(i,k,wtrc_iawset(iwtvap,iwset))

                          !Assume Rayleigh distillation during vapor to ice formation:

                          !calculate change in vapor mass:
                          fr = (qloc0(i,k,wtrc_iawset(iwtvap,1))/(qloc0(i,k,wtrc_iawset(iwtvap,1))+&
                              pre_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter))
                          !Constrain ratio fraction:
                          fr = min(1._r8,max(0._r8,fr))

                          qloc(i,k,wtrc_iawset(iwtvap,iwset)) = qloc0(i,k,wtrc_iawset(iwtvap,iwset))*(fr**alpha) !vapor
                          if(mdst .eq. wtrc_iawset(iwtstsnow,iwset)) then !snow?
                            smass(i,iwset) = smass(i,iwset) + (qloc0(i,k,wtrc_iawset(iwtvap,iwset))-&
                                qloc(i,k,wtrc_iawset(iwtvap,iwset)))*pstate%pdel(i,k)
                          else !ice
                            qloc(i,k,mdst) = qloc(i,k,mdst) + (qloc0(i,k,wtrc_iawset(iwtvap,iwset))-&
                                qloc(i,k,wtrc_iawset(iwtvap,iwset)))
                          end if            
                        else
                          !liquid directly to ice:
                          if(mdst .eq. wtrc_iawset(iwtstsnow,iwset)) then !snow?
                            smass(i,iwset) = smass(i,iwset) + (R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime*pstate%pdel(i,k))/wtrc_niter
                          else !ice
                            qloc(i,k,mdst) = qloc(i,k,mdst)+R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter
                          end if
                          qloc(i,k,msrc) = qloc(i,k,msrc)-R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter
                        end if !wisotope
                        !------------------------------
                     ! else if((iwset .gt. 1) .and. (isrctype .eq. iwtstrain) .and. (idsttype .eq. iwtvap)) then !rain re-evaporation?
                     !   if((rmass(i,iwset) .gt. 0._r8) .and. ((icecldf(i,k) + liqcldf(i,k)) .le. 1.e-4_r8)) then !rain exists but no clouds (i.e. below cloud base

                          !determine isotope species:
                     !     ispec = iwspec(wtrc_iawset(iwtvap,iwset))

                          !calculate layer thickness in meters:
                     !     dz = pstate%zi(i,k) - pstate%zi(i,k+1)

                          !apply Stewart/equilibration:
                     !     call stewart_isoevap(ispec, R, rmass0(i,1), rmass(i,1), rmass0(i,iwset), qloc0(i,k,wtrc_iawset(iwtvap,1)), &
                     !                          qloc(i,k,wtrc_iawset(iwtvap,1)), qloc0(i,k,mdst),tloc(i,k), pstate%pmid(i,k), &
                     !                          pstate%pdel(i,k), dz, dtime, rmass(i,iwset), qloc(i,k,mdst))
                     !   else
                          !apply tendencies normally:
                     !     qloc(i,k,mdst) = qloc(i,k,mdst)+R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter
                     !     rmass(i,iwset) = rmass(i,iwset)-(R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime*pstate%pdel(i,k))/wtrc_niter
                     !   end if
                      else
                        alpha = 1._r8
                        !Add tendency to destination:
                        if(idsttype .eq. iwtstrain) then !rain?
                          rmass(i,iwset) = rmass(i,iwset) + (R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime*pstate%pdel(i,k))/wtrc_niter
                        else if(idsttype .eq. iwtstsnow) then !snow?
                          smass(i,iwset) = smass(i,iwset) + (R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime*pstate%pdel(i,k))/wtrc_niter
                        else !not precip?
                          qloc(i,k,mdst) = qloc(i,k,mdst)+alpha*R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter
                        end if 
                        !Subtract tendency from source:
                        if(isrctype .ne. idsttype) then
                          if(isrctype .eq. iwtstrain) then !rain?
                            rmass(i,iwset) = rmass(i,iwset) - (R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime*pstate%pdel(i,k))/wtrc_niter
                          else if(isrctype .eq. iwtstsnow) then !snow?
                            smass(i,iwset) = smass(i,iwset) - (R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime*pstate%pdel(i,k))/wtrc_niter
                          else !not precip?
                            qloc(i,k,msrc)  = qloc(i,k,msrc)-alpha*R*pre_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter
                          end if
                        end if !subtraction
                      end if !ice/snow distillation
                    end do !iwset
!                    call wtrc_check_tracer_mass(pbuf,pstate,qloc,isrctype,idsttype,micro,iter,dtime,1e-10_r8,.false.) !check mass
                    qloc0(i,k,:) = qloc(i,k,:) !update state
                    rmass0(i,:)  = rmass(i,:)  !update precip
                    smass0(i,:)  = smass(i,:)
                  end if
                end do !idsttype
                !---------------------------------
                !equilibrate vapor and cloud water
                !---------------------------------
                do iwset = 1,wtrc_nwset
                  if(wisotope) then !are you using water isotopes?
                    if(iwset .ne. 1) then !not H2O tracer?
                      !set fractionation factor
                      !NOTE:  Assumes RH = 100%
                       alpha = wtrc_get_alpha(qloc0(i,k,wtrc_iawset(iwtvap,WTRC_WSET_STD)), tloc(i,k), iwspec(wtrc_iawset(iwtvap,iwset)), &
                                              iwtvap, iwtliq, .false.,1._r8,.false.)
                      !qloc:
                       call wtrc_liqvap_equil(alpha, 1._r8, qloc(i,k,wtrc_iawset(iwtvap,1)), qloc(i,k,wtrc_iawset(iwtliq,1)),&
                                             qloc(i,k,wtrc_iawset(iwtvap,iwset)), qloc(i,k,wtrc_iawset(iwtliq,iwset)), dliqiso)

                      !ECHAM5 fraction (debugging):
                      ! call wtrc_liqvap_equil(alpha, 0.95_r8, qloc(i,k,wtrc_iawset(iwtvap,1)), qloc(i,k,wtrc_iawset(iwtliq,1)),&
                      !                       qloc(i,k,wtrc_iawset(iwtvap,iwset)), qloc(i,k,wtrc_iawset(iwtliq,iwset)), dliqiso)

                    end if
                  end if
                end do
                qloc0(i,k,:) = qloc(i,k,:) !update state
                !---------------------------------
              end do !isrctype
              !--------------------------------------------
              !Calculate precipitation melting and freezing 
              !--------------------------------------------
              if(ldo_stprecip) then !only do in microphysics
                !NOTE:  melting and freezing in-theory should never occur at the same time, but occasionally it does anyways. -JN
                do iwset = 1,wtrc_nwset
                  !melting:
                  R = wtrc_ratio(iwspec(wtrc_iawset(iwtstrain,iwset)),smass0(i,iwset),smass0(i,1))
                  rmass(i,iwset) = rmass(i,iwset) + R*meltso(i,k)*dtime
                  smass(i,iwset) = smass(i,iwset) - R*meltso(i,k)*dtime
                end do
                rmass0(i,:)  = rmass(i,:)  !update precip
                smass0(i,:)  = smass(i,:)
                do iwset = 1,wtrc_nwset 
                  !freezing:
                  R = wtrc_ratio(iwspec(wtrc_iawset(iwtstrain,iwset)),rmass0(i,iwset),rmass0(i,1))
                  smass(i,iwset) = smass(i,iwset) + R*frzro(i,k)*dtime
                  rmass(i,iwset) = rmass(i,iwset) - R*frzro(i,k)*dtime
                end do
                rmass0(i,:)  = rmass(i,:)  !update precip
                smass0(i,:)  = smass(i,:)
              end if
              !---------------------------------
              !Equilibrate Rain below Clouds/LCL
              !---------------------------------
              if(ldo_stprecip) then !only do precipitation 
                do iwset = 1,wtrc_nwset
                  if(wisotope .and. (iwset .ne. 1)) then
                     
                    !NOTE:  0.0001 (1.e-4) appears to be the cutoff in MG1, such that
                    !any cloud fraction value equal to or below that number is assumed to be clear sky. - JN

                    if((rmass(i,iwset) .gt. 0._r8) .and. ((icecldf(i,k) + liqcldf(i,k)) .le. 1.e-4_r8)) then !rain exists but no clouds (i.e. below cloud base)

                      !determine isotope species:
                      ispec = iwspec(wtrc_iawset(iwtvap,iwset))

                      !calculate layer thickness in meters:
                      dz = pstate%zi(i,k) - pstate%zi(i,k+1)

                      !calculate equilibrium fractionation factor:
                      alpha = wtrc_get_alpha(qloc0(i,k,wtrc_iatype(1,iwtvap)), tloc(i,k), &
                                             iwspec(wtrc_iatype(iwset,iwtvap)), &
                                             iwtvap, iwtliq, .false.,1._r8,.false.)

                      !set molecular diffusion ratios:
                      difrm = (/ 1._r8, 1._r8, 0.9757_r8, 0.9727_r8 /)         ! Merlivat 1978 (direct from paper)

                      if (rmass(i,1) .gt. 0._r8) then
                        !calculate raindrop radius based off rain rate:
                        !NOTE: mean diameter in mm = 4/(4.1*R^-0.21), with R in mm/hr, assuming
                        !a Marshall-Palmer distribution. -JN
                        radius = 2._r8/(4.1_r8*((rmass(i,1)/gravit/dtime)*60._r8*60._r8)**-0.21_r8)

                        !convert radius to meters:
                        radius = radius/1000._r8

                        !calculate equilibration fraction:
                        call wtrc_equil_time(ispec,tloc(i,k),pstate%pmid(i,k),radius,dz,alpha,difrm(ispec),fequil)

                        !Constrain possible values:
                        fequil = min(1._r8,max(0._r8,fequil))
                      else
                        fequil = 1._r8
                      end if

                      !make copies of water mass quantities (to avoid changes from wtrc_liqvap_equil call):
                      qtmp = qloc(i,k,wtrc_iatype(1,iwtvap)) * pstate%pdel(i,k) 
                      itmp = qloc(i,k,wtrc_iatype(iwset,iwtvap)) * pstate%pdel(i,k)
                      rtmp = rmass(i,1) 
                      irtmp= rmass(i,iwset)

                      call wtrc_liqvap_equil(alpha, fequil, qtmp, rtmp,&
                                             itmp, irtmp, dliqiso) 

                      !vapor:
                      qloc(i,k,wtrc_iatype(iwset,iwtvap)) = qloc(i,k,wtrc_iatype(iwset,iwtvap)) - &
                                                            (dliqiso / pstate%pdel(i,k))
                      !rain:
                      rmass(i,iwset) = rmass(i,iwset) + dliqiso
               
                    end if !cloud base
                    
                  end if
                end do
                qloc0(i,k,:) = qloc(i,k,:) !update state
                rmass0(i,:)  = rmass(i,:)  !update precip
                smass0(i,:)  = smass(i,:)
              end if
              !-------------------------
              !correct any precip errors
              !-------------------------
              !NOTE:  May not be needed. - JN
              if(wisotope .and. ldo_stprecip) then
               !rain:
                pdiff = rmass(i,2)-rmass(i,1) !2=H216O,1=H2O
                do iwset=2,wtrc_nwset         !don't loop over H2O
                  R = wtrc_ratio(iwspec(wtrc_iawset(iwtstrain,iwset)),rmass0(i,iwset),rmass0(i,2)) !use H216O as denominator
                  rmass(i,iwset) = rmass(i,iwset) - R*pdiff                !adjust rain mass
                  mdst = wtrc_iawset(iwtvap,iwset)                         !assume extra mass goes to vapor (largest atmospheric reservoir)
                  qloc(i,k,mdst) = qloc(i,k,mdst)+R*pdiff/pstate%pdel(i,k) !add mass to vapor
                end do
               !snow:
                pdiff = smass(i,2)-smass(i,1) !2=H216O,1=H2O
                do iwset=2,wtrc_nwset         !don't loop over H2O
                  R = wtrc_ratio(iwspec(wtrc_iawset(iwtstrain,iwset)),smass0(i,iwset),smass0(i,2)) !use H216O as denominator
                  smass(i,iwset) = smass(i,iwset) - R*pdiff                  !adjust rain mass
                  mdst = wtrc_iawset(iwtvap,iwset)                           !assume extra mass goes to vapor (largest atmospheric reservoir)
                  qloc(i,k,mdst) = qloc(i,k,mdst)+R*pdiff/pstate%pdel(i,k)   !add mass to vapor
                end do
                qloc0(i,k,:) = qloc(i,k,:) !update state
                rmass0(i,:)  = rmass(i,:)  !update precip
                smass0(i,:)  = smass(i,:)
              end if
              !------------------
              !update temperature
              !------------------
              !NOTE ON TEMPERATURE CHANGE: 
              !Arguably, the most physically realistic method of temperature integration would
              !be to loop over all processes each time the temperature is changed, as this would better
              !reflect the fact that in reality, all of these processes are occuring at once.  However,
              !doing so appears difficult to accomplish given the current setup, and it also doesn't match the
              !actual MG scheme (where the pre-rates are not impacted by the post-rates, which wouldn't
              !be true given a full loop). Instead, each section (excluding sedimentation at the moment)
              !will be looped independently, updating the final temperature by the amount the latent heating
              !is updated in the actual cloud physics routines. This allows the water tracer/
              !isotope routines to be algorithmically as close as possible to the actual model physics - JN

              !NOTE: Temperature currently calculated as an average value over the entire parameterization time period. - JN
!              tloc(i,k) = tloc(i,k) + prelat(i,k)/cpair*dtime/wtrc_niter
              !NOTE: This temperature call simply updates the temperature to its value after pre-rate effects.
               tloc(i,k) = pstate%t(i,k) + prelat(i,k)/cpair*dtime/wtrc_niter        
              !------------------  
            end do!atmospheric columns (i)
          end do   !atmospheric levels (k)
        end do !iter
      end if !pre-rates

     !*************
     !Sedimentation
     !*************
     if (present(sed_rates)) then
       call wtrc_sediment(wtrc_niter,ncol,pstate%lchnk,dtime,fc,fi,liqcldf,icecldf,pstate%pdel,pbuf,tloc,qloc)
     end if
  
     qloc0(:ncol,:,:) = qloc(:ncol,:,:) !reset control quantity

     !****************************
     !Post-Sedimentation Processes
     !**************************** 
      if (present(post_rates)) then
        !Do this top down, like is done in MG microphysics.
        do iter=1,wtrc_niter !temperature iterator
          do k = 1, pver
            do i = 1, ncol
              !----------------------------------------------------
              !calculate average temperature during cloud processes
              !----------------------------------------------------
              tloc(i,k) = (tloc(i,k) + (tloc(i,k)+postlat(i,k)/cpair*dtime))/2._r8
              !----------------------------------------------------
              do isrctype=1,pwtype !loop through moisture sources
                do idsttype=1,pwtype !loop through water types that are modified
                  rtype = isrctype  !the moisture source is what determines the value of R
                  if(post_rates(i,k,idsttype,isrctype,rtype) .gt. 0._r8) then !make sure destination type is increasing
                    do iwset=1,wtrc_nwset !loop over water tracers/isotopes

                      msrc = wtrc_iawset(isrctype,iwset) !source water index
                      mbase= wtrc_iawset(isrctype,1)     !H2O source water index
                      mdst = wtrc_iawset(idsttype,iwset) !destination water index

                      if(isrctype .gt. iwtice) then !Is precipitation the source?
                        !--------------------
                        !use column integrals
                        !--------------------
                        !NOTE:  This section is currently never called. - JN
                        if(k .gt. 1) then !don't do for top level
                          if(sum(qloc0(i,1:k,mbase)*pstate%pdel(i,1:k)) .gt. 0._r8) then  !prevent NaNs
                            R = sum(qloc0(i,1:k,msrc)*pstate%pdel(i,1:k))/sum(qloc0(i,1:k,mbase)*pstate%pdel(i,1:k)) !calculate integral ratio
                          else
                            R = wtrc_ratio(iwspec(msrc), qloc0(i,k,msrc),qloc0(i,k,mbase)) !just do for single level
                          end if
                        else
                          R = wtrc_ratio(iwspec(msrc), qloc0(i,k,msrc),qloc0(i,k,mbase)) !calculate for single level
                        end if
                        !----------------------      
                      else
                        R = wtrc_ratio(iwspec(msrc),qloc0(i,k,msrc),qloc0(i,k,mbase)) !calculate ratio for single level
                      end if

                      if(wisotope .and. (iwset .ne. 1) .and. (isrctype .eq. iwtvap) .and. &
                         ((idsttype .eq. iwtice) .or. (idsttype .eq. iwtstsnow))) then !only do Rayleigh distillation on water isotopes.

                        ispec = iwspec(mdst) !determine isotope species of ice or snow (and thus vapoor)
                        alpha = wtrc_get_alpha(qloc0(i,k,wtrc_iawset(iwtvap,WTRC_WSET_STD)), tloc(i,k), ispec, isrctype, idsttype, &
                                            .true., pstate%pmid(i,k))

                        !calculate change in vapor mass:
                        fr = qloc(i,k,mbase)/qloc0(i,k,mbase)
                        !Constrain ratio fraction:
                        fr = min(1._r8,max(0._r8,fr))

                        qloc(i,k,msrc) = qloc0(i,k,msrc)*(fr**alpha)
                        qloc(i,k,mdst) = qloc(i,k,mdst) + (qloc0(i,k,msrc)-qloc(i,k,msrc))
                      else
                        alpha = 1._r8
                        !Add tendency to destination:
                        qloc(i,k,mdst)  = qloc(i,k,mdst) + alpha*R*post_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter
                        !Subtract tendency from source:
                        if(isrctype .ne. idsttype) then
                          qloc(i,k,msrc)  = qloc(i,k,msrc) - alpha*R*post_rates(i,k,idsttype,isrctype,rtype)*dtime/wtrc_niter
                        end if
                      end if !ice/snow distillation

                    end do !iwset
!                    call wtrc_check_tracer_mass(pbuf,pstate,qloc,isrctype,idsttype,micro,iter,dtime,1e-10_r8,.true.) !check mass
                    qloc0(i,k,:) = qloc(i,k,:) !update state
                  end if
                end do
                !---------------------------------
                !equilibrate vapor and cloud water
                !---------------------------------
                do iwset = 1,wtrc_nwset
                  if(wisotope) then !are you using water isotopes?
                    if(iwset .ne. 1) then !not H2O tracer?
                      !set fractionation factor
                      !NOTE:  Assumes RH = 100%
                      alpha = wtrc_get_alpha(qloc0(i,k,wtrc_iawset(iwtvap,WTRC_WSET_STD)), tloc(i,k), iwspec(wtrc_iawset(iwtvap,iwset)), &
                                             iwtvap, iwtliq, .false.,1._r8,.false.)
                      !qloc:
                       call wtrc_liqvap_equil(alpha, 1._r8, qloc(i,k,wtrc_iawset(iwtvap,1)), qloc(i,k,wtrc_iawset(iwtliq,1)),&
                                             qloc(i,k,wtrc_iawset(iwtvap,iwset)), qloc(i,k,wtrc_iawset(iwtliq,iwset)), dliqiso)
                    end if
                  end if
                end do
                qloc0(i,k,:) = qloc(i,k,:) !update state
                !--------------------------------
              end do !isrctype
              !------------------
              !update temperature
              !------------------
              !NOTE: Temperature currently calculated as an average value over the entire parameterization time period. - JN
              !tloc(i,k) = tloc(i,k) + postlat(i,k)/cpair*dtime/wtrc_niter 
              !NOTE: This temperature call simply updates the temperature to its value after pre-rate effects.
              tloc(i,k) = pstate%t(i,k) + (prelat(i,k)+postlat(i,k))/cpair*dtime/wtrc_niter
              !------------------
              !-----------------------------------------
              !partially equilibrate vapor and cloud ice (to represent slow equilibration of ice crystals in atmosphere)
              !-----------------------------------------
!              do iwset = 1,wtrc_nwset
!                if(wisotope) then !are you using water isotopes?
!                  if(iwset .ne. 1) then !not H2O tracer?
                    !set fractionation factor
                    !NOTE:  Assumes RH = 100%
!                    alpha = wtrc_get_alpha(qloc0(i,k,wtrc_iawset(iwtvap,WTRC_WSET_STD)), tloc(i,k), iwspec(wtrc_iawset(iwtvap,iwset)), &
!                                           iwtvap, iwtice, .false.,1._r8,.false.)
!                    !qloc:
!                     call wtrc_liqvap_equil(alpha, 0.02_r8, qloc(i,k,wtrc_iawset(iwtvap,1)), qloc(i,k,wtrc_iawset(iwtice,1)),&
!                                           qloc(i,k,wtrc_iawset(iwtvap,iwset)), qloc(i,k,wtrc_iawset(iwtice,iwset)), dliqiso)
!                  end if
!                end if
!              end do
              !--------------------------------
            end do!atmospheric columns (i)
          end do   !atmospheric levels (k)
        end do !iter loop
      end if !post-rates

      ! Use the tendencies to update the state for the next iteration, and also
      ! round any small numbers to avoid negative numbers and numerical problems.
!      do iwset = 1, wtrc_nwset

        ! Eliminate negative mixing ratios.
        ! NOTE: Use 0 as a threshold rather than wtrc_qmin during the substepping.
!        call wtrc_adjust_h2o("WTRC_APPLY_RATES(" // trim(ptend_sum%name) // ")", iwset, ncol, qloc, qmin=0._r8) 
!      end do
  
      do idsttype = 1, pwtype
        qloc(:ncol, :, wtrc_bulk_indices(idsttype)) = qloc(:ncol, :, wtrc_bulk_indices(idsttype)) + &
          ptend_sum%q(:ncol, :, wtrc_bulk_indices(idsttype)) * dtime 
          
!        call wtrc_adjust_h2o_bulk("WTRC_APPLY_RATES(" // trim(ptend_sum%name) // ")", ncol, qloc) 
      end do

    ! Now that we have the final state, apply q limits like those used by the
    ! MG microphysics.
    !
    ! NOTE: There is some inconsistency here, since other parameterizations and
    ! physics_update may use different limits. It would be nice if all
    ! parameterizations used cnst_qmin.
!    do iwset = 1, wtrc_nwset
!      call wtrc_adjust_h2o("WTRC_APPLY_RATES(" // trim(ptend_sum%name) // ")", iwset, ncol, qloc) 
!    end do
      
    ! Sedimentation of cloud liquid and cloud ice is stratiform precipitation.
    ! out of the bottom layer is precipitation.
    !
    ! TBD - Can't do this until the sedimentation rates are provided for motion to
    ! and from layer rather than change in the layer.
    ! NOTE:  No longer needed (handled in wtrc_sediment) - JN.
  !  if (present(sed_rates)) then
  !    if (ldo_stprecip) then
  !      call wtrc_collect_precip(pstate, pbuf, iwtliq, iwtstrain, sed_rates, dtime)
  !      call wtrc_collect_precip(pstate, pbuf, iwtice, iwtstsnow, sed_rates, dtime)
  !    end if
  !  end if

    ! Diagnose stratiform precipitation, since CAM currently does diagnostic rather
    ! than prognostic precipitation. If it was prognostic, then sedimentation of
    ! rain and snow would also generate surface precipitation.
    if (ldo_stprecip) then
     !original version:
!      call wtrc_diagnose_precip(pstate, qloc, pbuf, iwtstrain, dtime)
!      call wtrc_diagnose_precip(pstate, qloc, pbuf, iwtstsnow, dtime)

     !modified version:
      call wtrc_diagnose_precip(pstate, rmass, pbuf, iwtstrain, dtime)
      call wtrc_diagnose_precip(pstate, smass, pbuf, iwtstsnow, dtime)

      call wtrc_diagnose_bulk_precip(pstate, ptend_sum, iwtstrain, dtime)
      call wtrc_diagnose_bulk_precip(pstate, ptend_sum, iwtstsnow, dtime)
    end if

    !Calculate net tendencies:
    diff(:,:,:) = 0._r8 !set initial difference to zero!
    do i=1,ncol
      do k=1,pver
        do icnst = 1, wtrc_ncnst
         
          !Allow for tendency to change state:
!          ptend_sum%lq(wtrc_indices(icnst)) = .true. 

          !Calculate state:
          ptend_sum%q(i,k,wtrc_indices(icnst)) = (qloc(i,k,wtrc_indices(icnst)) - pstate%q(i,k,wtrc_indices(icnst))) / dtime
 
          if(icnst .le. pwtype) then !H2O only
            diff(i,k,icnst) = ptend_sum%q(i,k,wtrc_indices(icnst))-(ptend_sum%q(i,k,wtrc_bulk_indices(icnst)))
          end if

        end do
      end do
    end do

!---------------------------------
!Numerical corrections:
!---------------------------------

  !check diff value to make sure it isn't too large:
  if(sum(abs(ptend_sum%q(:ncol,:,wtrc_indices(1:pwtype)))) .ne. 0._r8) then
    if(sum(abs(diff(:ncol,:pver,:)))/sum(abs(ptend_sum%q(:ncol,:,wtrc_indices(1:pwtype)))) .gt. 1.e-10_r8) then
      write(iulog,*) 'WARNING: BIG ERROR, diff value:',sum(abs(diff(:ncol,:,:)))/sum(abs(ptend_sum%q(:ncol,:,wtrc_indices(1:pwtype)))),&
      sum(abs(diff(:ncol,:,:))),sum(abs(diff(:ncol,:,1:3))),sum(qloc(:ncol,:pver,wtrc_iatype(1,:))-qloc(:ncol,:pver,wtrc_bulk_indices(:))),&
      sum(abs(ptend_sum%q(:ncol,:,wtrc_indices(1:pwtype))))
    end if
  end if

!Apply correction and divide by time:
  do i=1,ncol
    do k=1,pver
      do icnst = 1,pwtype  !only loop over bulk_water!
        do m = 1,wtrc_nwset    !loop over water species

          if(m .eq. 1) then !is it the standard water tracer?
            qtmp = ptend_sum%q(i,k,wtrc_iatype(m,icnst))
          end if  

          !calculate water tracer ratio:
          R = wtrc_ratio(iwspec(wtrc_iatype(m,icnst)),ptend_sum%q(i,k,wtrc_iatype(m,icnst)),qtmp)

          ptend_sum%q(i,k,wtrc_iatype(m,icnst)) = ptend_sum%q(i,k,wtrc_iatype(m,icnst))-R*diff(i,k,icnst) 

        end do
      end do
    end do
  end do
 
!Apply correction again (I don't know why this is needed, but for some reason it is, as apparently the fix above produces
!another numerical error for cloud liquid - JN)
  do i=1,ncol
    do k=1,pver
      do icnst = 1,pwtype  !only loop over bulk_water!
        do m = 1,wtrc_nwset    !loop over water species

          if(m .eq. 1) then !Only do for H2O
            diff(i,k,icnst) = ptend_sum%q(i,k,wtrc_iatype(m,icnst))-ptend_sum%q(i,k,wtrc_bulk_indices(icnst))
            qtmp = ptend_sum%q(i,k,wtrc_iatype(m,icnst))
          end if

          !calculate water tracer ratio:
          R = wtrc_ratio(iwspec(wtrc_iatype(m,icnst)), ptend_sum%q(i,k,wtrc_iatype(m,icnst)),qtmp)

          ptend_sum%q(i,k,wtrc_iatype(m,icnst)) = ptend_sum%q(i,k,wtrc_iatype(m,icnst))-R*diff(i,k,icnst)

        end do
      end do
    end do
  end do

  !verify again that the error isn't hiding an actual error:
  if(sum(abs(ptend_sum%q(:ncol,:,wtrc_indices(1:pwtype)))) .ne. 0._r8) then
    if(sum(abs(diff(:ncol,:pver,:)))/sum(abs(ptend_sum%q(:ncol,:,wtrc_indices(1:pwtype)))) .gt. 1.e-20_r8) then
      write(iulog,*) 'WARNING: BIG ERROR AGAIN, diff value:',sum(abs(diff(:ncol,:,:)))/sum(abs(ptend_sum%q(:ncol,:,wtrc_indices(1:pwtype)))),&
      sum(abs(diff(:ncol,:,:))),sum(abs(diff(:ncol,:,1:3))),sum(qloc(:ncol,:pver,wtrc_iatype(1,:))-qloc(:ncol,:pver,wtrc_bulk_indices(:)))
    end if
  end if

!------------------------------
!End Numerical corrections
!------------------------------

    ! Optionally check to see if the total water is conserved by all
    ! non-tagged water wsets.
    !
    ! NOTE: For the h216o model, if the pre-test indicated an error in
    ! total mass, but the post-test indicates the same error, then the
    ! processes being applied were mass conserving and the problem is
    ! with the initial state.
   ! isOk = wtrc_check_h2o("WTRC_APPLY_RATES(POST-" // trim(ptend_sum%name) // ")", pstate, qloc, dtime, ptend=ptend_sum)

end if !tracer_water
  
  return
end subroutine wtrc_apply_rates

!====================================================================
subroutine wtrc_check_tracer_mass(pbuf,pstate,qloc,stype,dtype,micro,iter,dtime,level,prec)
!
!Purpose:  To determine if any mass was lost or gained during the application of a tendency.
!
!NOTE:  Currently only uses the last water tracer (wtrc_nwset).  Will be fixed later... - JN
!
!Written by:  Jesse Nusbaumer <nusbaume@colorado.edu> - April, 2013
!
!**************
!Use statements:
!**************

use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
use physconst,      only: gravit, rhoh2o
use constituents,   only: cnst_name

!*****************
!Delcare variables:
!*****************

!Input variables:
type(physics_buffer_desc), pointer :: pbuf(:)  !physics buffer     
type(physics_state), intent(in)                      :: pstate          ! State of the atmosphere

real(r8),intent(in)   :: qloc(pcols,pver,pcnst) !state of the atmosphere after tendency change.

integer, intent(in)   :: stype !source water type 
integer, intent(in)   :: dtype !destination water type
integer, intent(in)   :: iter  !iteration loop number
real(r8), intent(in)  :: dtime !model time step
real(r8), intent(in)  :: level !mass change cutoff
logical, intent(in)   :: micro !micro or macro physics?
logical, intent(in)   :: prec  !call precip values from pbuff?

!Local variables:
integer               :: sedrainidx !physics buffer index for tracer rain
integer               :: sedsnowidx !physics buffer index for tracer snow
integer               :: sdrnh2oidx !physics buffer index for H2O rain
integer               :: sdswh2oidx !physics buffer index for h2o snow

real(r8)              :: qtot,qorg,qprc,qtoth,qorgh,qprch !Mass checkers
real(r8)              :: sedrain(pcols) !tracer sed. rain
real(r8)              :: sedsnow(pcols) !tracer sed. snow
real(r8)              :: sdrnh2o(pcols) !H2O sed. rain
real(r8)              :: sdswh2o(pcols) !H2O sed. snow

real(r8), pointer, dimension(:)    :: sedrainpt !tracer sed. rain pointer
real(r8), pointer, dimension(:)    :: sedsnowpt !tracer sed. snow pointer
real(r8), pointer, dimension(:)    :: sdrnh2opt !H2O sed. rain pointer
real(r8), pointer, dimension(:)    :: sdswh2opt !H2O sed. snow pointer

integer :: i  !loop control variable

!*****************************************
!read in precipitation from physics buffer
!*****************************************

 if(prec) then !grab precipitation in buffer

   !tracer water:
    call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstrain, wtrc_nwset), sedrainpt)
    call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstsnow, wtrc_nwset), sedsnowpt)
    sedrain = sedrainpt * gravit * rhoh2o * dtime !convert to kg/m/s2
    sedsnow = sedsnowpt * gravit * rhoh2o * dtime
   !H2O water:
    call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstrain, 1), sdrnh2opt)
    call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstsnow, 1), sdswh2opt)
    sdrnh2o = sdrnh2opt * gravit * rhoh2o * dtime
    sdswh2o = sdswh2opt * gravit * rhoh2o * dtime

 else !set everything to zero
   sedrain(:) = 0._r8
   sedsnow(:) = 0._r8
   sdrnh2o(:) = 0._r8
   sdswh2o(:) = 0._r8
 end if 

!*****************
!Loop over columns
!*****************

do i=1,pstate%ncol

!**************************
!calculate column integrals
!**************************

  !Isotope/water tracer:
   qtot = 0._r8 !initalize variables
   qorg = 0._r8
   qprc = 0._r8
   qtot = sum(qloc(i,:,wtrc_iawset(iwtvap,wtrc_nwset))*pstate%pdel(i,:))+&
          sum(qloc(i,:,wtrc_iawset(iwtliq,wtrc_nwset))*pstate%pdel(i,:))+& !sum of HDO 
          sum(qloc(i,:,wtrc_iawset(iwtice,wtrc_nwset))*pstate%pdel(i,:))
   qorg = sum(pstate%q(i,:,wtrc_iawset(iwtvap,wtrc_nwset))*pstate%pdel(i,:))+&
          sum(pstate%q(i,:,wtrc_iawset(iwtliq,wtrc_nwset))*pstate%pdel(i,:))+& !sum of 
          sum(pstate%q(i,:,wtrc_iawset(iwtice,wtrc_nwset))*pstate%pdel(i,:))
   qprc = sedrain(i)+sedsnow(i)+sum(qloc(i,:,wtrc_iawset(iwtstrain,wtrc_nwset))*pstate%pdel(i,:))+&
          sum(qloc(i,:,wtrc_iawset(iwtstsnow,wtrc_nwset))* pstate%pdel(i,:))
  !H2O (standard) tracer:
   qtoth = 0._r8 !initalize variables
   qorgh = 0._r8
   qprch = 0._r8
   qtoth = sum(qloc(i,:,wtrc_iawset(iwtvap,1))*pstate%pdel(i,:))+&
           sum(qloc(i,:,wtrc_iawset(iwtliq,1))*pstate%pdel(i,:))+& !sum of HDO 
           sum(qloc(i,:,wtrc_iawset(iwtice,1))*pstate%pdel(i,:))
   qorgh = sum(pstate%q(i,:,wtrc_iawset(iwtvap,1))*pstate%pdel(i,:))+&
           sum(pstate%q(i,:,wtrc_iawset(iwtliq,1))*pstate%pdel(i,:))+& !sum of
           sum(pstate%q(i,:,wtrc_iawset(iwtice,1))*pstate%pdel(i,:))
   qprch = sdrnh2o(i)+sdswh2o(i)+sum(qloc(i,:,wtrc_iawset(iwtstrain,1))*pstate%pdel(i,:))+&
           sum(qloc(i,:,wtrc_iawset(iwtstsnow,1))*pstate%pdel(i,:))

!******************
!Perform mass check
!******************

   !mass check:
   if(qorg .ne. 0._r8) then
     if(abs(qorg-qtot-qprc)/qorg .gt. level) write(iulog,*) 'conservation error here for HDO!',stype,dtype,micro,iter,(qorg-qtot-qprc)/qorg,qorg-qtot,&
                                                                    qorgh-qtoth,qprc,qprch,qorg-qorgh,qprc-qprch,&
                                                                    sum(qloc(i,:,wtrc_iawset(iwtstrain,2))*&
                                                                    pstate%pdel(i,:))-sum(qloc(i,:,wtrc_iawset(iwtstrain,1))*&
                                                                    pstate%pdel(i,:)),&
                                                                    sum(qloc(i,:,wtrc_iawset(iwtstsnow,2))* pstate%pdel(i,:))-&
                                                                    sum(qloc(i,:,wtrc_iawset(iwtstsnow,1))* pstate%pdel(i,:)),&
                                                                    sum(qloc(i,:,wtrc_iawset(iwtvap,2))*pstate%pdel(i,:))-&
                                                                    sum(qloc(i,:,wtrc_iawset(iwtvap,1))*pstate%pdel(i,:)),&
                                                                    sum(qloc(i,:,wtrc_iawset(iwtice,2))*pstate%pdel(i,:))-&
                                                                    sum(qloc(i,:,wtrc_iawset(iwtice,1))*pstate%pdel(i,:)),&
                                                                    sedsnow(i)-sdswh2o(i),sedrain(i)-sdrnh2o(i),(qprc-sedsnow(i)-&
                                                                    sedrain(i))-(qprch-sdswh2o(i)-sdrnh2o(i)),prec,i
   end if
   if(qorgh .ne. 0._r8) then
     if(abs(qorgh-qtoth-qprch)/qorgh .gt. level) write(iulog,*) 'conservation error here for H2O too!',iter,(qorgh-qtoth-qprch)/qorgh,prec
   end if

!***************
!End column loop
!***************

 end do

end subroutine wtrc_check_tracer_mass

!=============================================================================================
subroutine wtrc_sediment(wtrc_niter,ncol,lchnk,dtime,wtfc,wtfi,liqcldf,icecldf,pdel,pbuf,tloc,qloc)
!--------------------------------------
!
!Purpose:  To calculate the sedimentation rates of water tracers and its
!impact on different water tracer quantities (vapor, condensate, precipitation).
!
!NOTE:  This code is needed because the tendency at a particular grid point
!due to sedimentation is caused by the flux out of the grid box minus the flux
!coming in from the grid box above.  Thus both the condensate ratio of the grid
!box and the grid box above are needed, along with the sediment fall velocities.
!These quantities are difficult to extract from the output tendencies, so
!instead the sedimentation rates for water tracers are calculated directly, via
!this subroutine. - JN
!
!Written by:  Jesse Nusbaumer <nusbaume@colorado.edu> - January 2013
!
!--------------------------------------

!***************************
!use statements/header files
!***************************

use physconst,      only: gravit, cpair, latvap, latice
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
use constituents,   only: cnst_name

!*****************
!Declare variables:
!*****************

!Input arguments:
integer, intent(in)  :: wtrc_niter          !water tracer temperature iteration  
integer, intent(in)  :: ncol                !number of columns in chunk
integer, intent(in)  :: lchnk               !chunk number
real(r8), intent(in) :: dtime               !model time step
real(r8), intent(in) :: wtfc(pcols,pver)    !initial fall velocity of cloud liquid
real(r8), intent(in) :: wtfi(pcols,pver)    !initial fall velocity of cloud ice
real(r8), intent(in) :: liqcldf(pcols,pver) !liquid cloud fraction (unitless)
real(r8), intent(in) :: icecldf(pcols,pver) !ice cloud fraction    (unitless)
real(r8), intent(in) :: pdel(pcols,pver)    !change in pressure per vertical level

!Input/Output arguments:
type(physics_buffer_desc), pointer :: pbuf(:) !physics buffer

real(r8), pointer, dimension(:)    :: srfpcp      ! Surface precipitation (m/s)

real(r8), intent(inout) :: tloc(pcols,pver)       !air temperature
real(r8), intent(inout) :: qloc(pcols,pver,pcnst) !current state of water tracers (kg/kg)

!Local variables:
integer  :: i,m,n,k                  !loop control variables
integer  :: nstep                    !CFL-stable sub-step
integer  :: srfpcidx                 !physics buffer index
real(r8) :: rgvm                     !fastest velocity (for CFL calculation)
real(r8) :: dum                      !dummy variable used for liquid cloud fraction ratios (unitless)
real(r8) :: dum1                     !dummy variable used for ice cloud fraction ratios (unitless)
real(r8) :: faltndqce                !change in liquid due to sedimentation (no-evap)
real(r8) :: faltndqie                !change in ice due to sedimentation (no-evap)
real(r8) :: faltndc                  !change in liquid due to sedimentation
real(r8) :: faltndi                  !change in ice due to sedimentation
real(r8) :: fc(pver)                 !cloud liquid fall velocity
real(r8) :: fi(pver)                 !cloud ice fall velocity  
real(r8) :: lcldm(pcols,pver)        !liquid cloud fraction
real(r8) :: icldm(pcols,pver)        !ice cloud fraction
real(r8) :: faloutc(pver,wtrc_nwset) !falling liquid amount
real(r8) :: falouti(pver,wtrc_nwset) !falling ice amount
real(r8) :: precr(pcols,wtrc_nwset)  !rain amount (m/s)
real(r8) :: preci(pcols,wtrc_nwset)  !snow amount (m/s)

!water isotopes:
real(r8) :: alpha                    !fractionation factor
real(r8) :: dliqiso                  !change in liquid due to isotopic equilibration
real(r8) :: tloc0(pcols,pver)        !copy of air temperature 

real(r8), parameter :: mincld = 0.0001_r8 !minimum cloud fraction (unitless)

!************************
!Save initial temperature
!************************

tloc0(:,:) = tloc(:,:) !copy initial temperature

!***********************
!calculate cloud amounts
!***********************

do i=1,ncol 
  do k=1,pver
    lcldm(i,k)=max(liqcldf(i,k),mincld)
    icldm(i,k)=max(icecldf(i,k),mincld)
  end do
end do

!*******************
!Initalize variables
!*******************

precr(:,:) = 0._r8
preci(:,:) = 0._r8

!******************
!Loop over columns:
!******************

do i=1,ncol   !loop over columns
 
  nstep = 1 !initalize variables

!********************************
!set initial fall velocity values
!********************************    

  fc(:) = wtfc(i,:)
  fi(:) = wtfi(i,:)

!******************
!calculate sub-step
!******************

  do k=1,pver
    rgvm = max(fi(k),fc(k))
    nstep = max(int(rgvm*dtime/pdel(i,k)+1._r8),nstep)
  end do

!******************
!Loop over substeps
!******************

  do n = 1,nstep  !loop over sub-time step to ensure stability  

!******************************
!Calculate sedimentation fluxes
!******************************
    do k = 1,pver
      do m=1,wtrc_nwset 
        falouti(k,m) = fi(k)*qloc(i,k,wtrc_iatype(m,iwtice))
        faloutc(k,m) = fc(k)*qloc(i,k,wtrc_iatype(m,iwtliq))
      end do
    end do

!******************************
!Calculate changes at model top
!******************************

    k = 1
    do m=1,wtrc_nwset
      faltndi = falouti(k,m)/pdel(i,k)
      faltndc = faloutc(k,m)/pdel(i,k)

      qloc(i,k,wtrc_iatype(m,iwtice)) = qloc(i,k,wtrc_iatype(m,iwtice))-faltndi*dtime/nstep 
      qloc(i,k,wtrc_iatype(m,iwtliq)) = qloc(i,k,wtrc_iatype(m,iwtliq))-faltndc*dtime/nstep 
    end do

!********************************************
!Calculate changes for the rest of the column
!********************************************

    do k=2, pver  

! for cloud liquid and ice, if cloud fraction increases with height
! then add flux from above to both vapor and cloud water of current level
! this means that flux entering clear portion of cell from above evaporates
! instantly

      dum=lcldm(i,k)/lcldm(i,k-1)
      dum=min(dum,1._r8)
      dum1=icldm(i,k)/icldm(i,k-1)
      dum1=min(dum1,1._r8)

      do m=1, wtrc_nwset

        faltndqie=(falouti(k,m)-falouti(k-1,m))/pdel(i,k)
        faltndi=(falouti(k,m)-dum1*falouti(k-1,m))/pdel(i,k)
        faltndqce=(faloutc(k,m)-faloutc(k-1,m))/pdel(i,k)
        faltndc=(faloutc(k,m)-dum*faloutc(k-1,m))/pdel(i,k)

!*********************************
!Add fluxes to water tracer states
!*********************************

        !sedimenting ice sublimation:
        qloc(i,k,wtrc_iatype(m,iwtvap)) = qloc(i,k,wtrc_iatype(m,iwtvap))-&
                                          (faltndqie-faltndi)*dtime/nstep 
        !sedimenting liquid evaporation:
        qloc(i,k,wtrc_iatype(m,iwtvap)) = qloc(i,k,wtrc_iatype(m,iwtvap))-&
                                          (faltndqce-faltndc)*dtime/nstep 

        qloc(i,k,wtrc_iatype(m,iwtice)) = qloc(i,k,wtrc_iatype(m,iwtice))-faltndi*dtime/nstep 
        qloc(i,k,wtrc_iatype(m,iwtliq)) = qloc(i,k,wtrc_iatype(m,iwtliq))-faltndc*dtime/nstep 

!***************************
!Ensure isotopic equilibrium
!***************************

        if(wisotope) then
          if(m .ne. 1) then !not H2O tracer?
            !set fractionation factor
            alpha = wtrc_get_alpha(qloc(i,k,wtrc_iatype(1,iwtvap)), tloc0(i,k), iwspec(wtrc_iatype(m,iwtvap)), &
                                         iwtvap, iwtliq, .false.,1._r8,.false.)

            !NOTE: Equilibrates over CFL sub-step.  This may not be the same as wtrc_niter, and thus it may not converge
            !as cleanly as other cloud physics processes. - JN
            call wtrc_liqvap_equil(alpha, 1._r8, qloc(i,k,wtrc_iatype(1,iwtvap)), qloc(i,k,wtrc_iatype(1,iwtliq)),&
                                  qloc(i,k,wtrc_iatype(m,iwtvap)), qloc(i,k,wtrc_iatype(m,iwtliq)), dliqiso)

          end if
        end if

!****************************
!calculate temperature change
!****************************

        if(m .eq. 1) then !only do for H2O
          tloc(i,k)=tloc(i,k)+(faltndqie-faltndi)*(latvap+latice)*dtime/cpair/nstep
          tloc(i,k)=tloc(i,k)+(faltndqce-faltndc)*latvap*dtime/cpair/nstep
        end if

      end do !water tracers

      !update fractionation temperature:
      tloc0(i,k) = tloc(i,k)       
        
!*********************
!reset fall velocities
!*********************

      fi(k)=max(fi(k)/pdel(i,k),fi(k-1)/pdel(i,k-1))*pdel(i,k)
      fc(k)=max(fc(k)/pdel(i,k),fc(k-1)/pdel(i,k-1))*pdel(i,k)

    end do  !pver (column) loop

!*********************************************
!Add up sedimentation flux out of bottom level
!*********************************************

    do m=1, wtrc_nwset

     !units below are m/s:
      precr(i,m) = precr(i,m)+(faloutc(pver,m)) &
                 /gravit/nstep/1000._r8
      preci(i,m) = preci(i,m)+(falouti(pver,m)) &
                 /gravit/nstep/1000._r8
    end do

!*********
!End loops
!*********

  end do !nstep loop
end do   !ncol loop

!************************************************
!Add summed bottom sediment flux to precipitation
!************************************************

do m=1, wtrc_nwset

  !Point to water tracer rain
  call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstrain,m), srfpcp)
  srfpcp(:) = srfpcp(:) + precr(:,m) !add rain to precipitation variable

  !Point to water tracer snow 
  call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstsnow,m), srfpcp)
  srfpcp(:) = srfpcp(:) + preci(:,m) !add snow to precipitation variable

end do

end subroutine wtrc_sediment

!============================================================================
subroutine stewart_isoevap(ispec, Rr0, rmass0, rmass, irmass0, vmass0, vmass, &
                           ivmass0, airtemp, airpres, pdel, dz, dtime, irmass, ivmass)

!--------------------------------------------------
!Purpose:
!
!To calculate the final isotopic rain ratio
!after the re-evaporation of rain droplets as they 
!fall through sub-cloud air.
!
!Written by:  Jesse Nusbaumer <nusbaume@colorado.edu> - June, 2015
!
!Equations from Stewart, 1975.
!
!------------------------------------------------

!**************
!Use statements
!**************

!use wv_saturation,  only: qsat_water
use wv_saturation,  only: qsat
use physconst,      only: gravit

!********************
!Variable declaration
!********************

!Input variables:
integer, intent(in)   :: ispec   !isotope species
real(r8), intent(in)  :: Rr0     !rain ratio pre-evaporation [unitless]
real(r8), intent(in)  :: rmass0  !bulk water mass pre-evaporation [kg/m/s2]
real(r8), intent(in)  :: rmass   !bulk water mass post-evaporationa [kg/m/s2]
real(r8), intent(in)  :: irmass0 !isotopic water mass pre-evaporation [kg/kg]
real(r8), intent(in)  :: vmass0  !bulk vapor mass pre-evaporation [kg/kg]
real(r8), intent(in)  :: vmass   !bulk vapor mass post-evaporation [kg/kg]
real(r8), intent(in)  :: ivmass0 !isotopic vapor mass pre-evaporation [kg/kig]
real(r8), intent(in)  :: airtemp !air temperature [Kelvin?]
real(r8), intent(in)  :: airpres !air pressure [Pa]
real(r8), intent(in)  :: pdel    !layer thickness [Pa]
real(r8), intent(in)  :: dz      !layer thickness [m]
real(r8), intent(in)  :: dtime   !model time step [s]

!Input/Output variables:
real(r8), intent(inout) :: irmass  !isotopic water mass post-evaporation [kg/kg]
real(r8), intent(inout) :: ivmass  !isotopic vapor mass post-evaporation [kg/kg]

!Local variables:
real(r8) :: fr      !fraction of precipitation mass remaining after evaporation.
real(r8) :: gam     !gamma paramter in Stewart equation
real(r8) :: bet     !beta parameter in stewart equation
real(r8) :: vsat    !saturated vapor mass [kg/kg]
real(r8) :: es      !saturated vapor pressure
real(r8) :: rh      !relative humidity
real(r8) :: ae      !equilibrium fractionation factor
real(r8) :: ak      !kinetic fractionation factor
real(r8) :: Rstw    !Stewart rain ratio
real(r8) :: Rr      !rain ratio post-evaporation    
real(r8) :: Rv0     !vapor ratio pre-evaporation
real(r8) :: tmpv    !temporary bulk vapor variable for equilibration
real(r8) :: itmpv   !temporary vapor variable for equilibration
real(r8) :: itmpr   !temporary rain variable for equilibration
real(r8) :: dliqiso !change in liquid mass due to equilibration [kg/m/s2]
real(r8) :: phi     !local enchancement of relative humidity (tuning parameter)
real(r8) :: heff    !local "effective" relative humidity
real(r8) :: radius  !assumed radius of raindrop [m] (tuning parameter)
real(r8) :: fequil  !fraction of rain the experiences stewart/equilibration

!For kinetic fractionation:
real(r8) :: difrm(pwtspec)
real(r8) :: dkfac    

!*****************
!Tuning parameters
!*****************

!phi value given in Bony et. al., 2008:
phi = 0.9_r8

!Assumed raindrop radius (in meters):
!NOTE:  Now calculated based off rain rate -JN
!radius = 0.00035_r8

!*****************************
!Diffusivity ratios and powers (should be called from water_isotopes.F90 eventually):
!*****************************

difrm = (/ 1._r8, 1._r8, 0.9757_r8, 0.9727_r8 /)         ! Merlivat 1978 (direct from paper)
dkfac = 0.58_r8                                          ! From Stewart, 1975?

!***************************
!Calculate relative humidity
!***************************

!calculate saturation vapor pressure over water:
!call qsat_water(airtemp, airpres, es, vsat)
call qsat(airtemp, airpres, es, vsat)

!Calculate RH:
rh = vmass0/vsat

!calculate effective RH:
heff = phi+((1_r8-phi)*rh) !From Bony et. al., 2008

!*******************************
!Calculate fractionation factors
!*******************************

!equilibrium fractionation factor:
ae = wtrc_get_alpha(vmass0,airtemp,ispec,iwtvap,iwtliq,.false.,1._r8,.false.)

!kinetic fractionation factor:
ak = (1._r8/difrm(ispec))**dkfac

!***********************************************
!Calculate mass-weighted average raindrop radius
!***********************************************

!NOTE: mean diameter in mm = 4/(4.1*R^-0.21), with R in mm/hr, assuming a
!Marshall-Palmer distribution. -JN

radius = 2._r8/(4.1_r8*((rmass0/gravit/dtime)*60._r8*60._r8)**-0.21_r8)

!convert radius to meters:
radius = radius/1000._r8

!***************************************************
!Calculate fraction of precip that will eqiuilibrate
!***************************************************

 call wtrc_equil_time(ispec,airtemp,airpres,radius,dz,ae,difrm(ispec),fequil)

 !Constrain possible values:
 fequil = min(1._r8,max(0._r8,fequil))

!*******************************************************
!Calculate fraction of precipitation lost to evaporation
!*******************************************************

 !Remaining mass fraction:
  fr = rmass/rmass0

 !Constrain possible values:
  fr = min(1._r8,max(0._r8,fr))

!********************************
!Apply Stewart equation if RH < 1
!********************************

if((heff .lt. 1._r8) .and. (fr .lt. 1._r8)) then

!*************************************
!Calculate Stewart equation parameters
!*************************************

  gam = (ae*heff)/(1._r8-ae*ak*(1._r8-heff))
  bet = (1._r8-ae*ak*(1._r8-heff))/(ae*ak*(1._r8-heff))

!*********************************
!Calculate new isotopic rain ratio
!*********************************

  !pre-evaporation vapor ratio:
  Rv0 = wtrc_ratio(ispec,ivmass0,vmass0)

  !Apply stewart equation:
  Rstw = ((Rr0-gam*Rv0)*(fr**bet)) + gam*Rv0

  !Apply equilibration fraction:
  if((fr .eq. 0._r8) .and. (bet .lt. 0._r8)) then !Prevent NaNs (plus there shouldn't be any evaporation anyways)
    Rr = Rr0  !just conserve ratio
  else
    Rr = fequil*Rstw + (1._r8-fequil)*Rr0  !perform fractionation
  end if

!*******************************
!Calculate isotopic mass changes
!*******************************

 dliqiso = irmass0 - Rr*rmass

!*************************************************************
!Calculate new isotopic mass values for current vertical level
!*************************************************************

  !NOTE:  1/pdel needed to convert quantity back
  !       to units of kg/kg.

  !isotopic rain value:
  irmass = irmass0 - dliqiso

  !isotopic vapor value:
  ivmass = ivmass0 + (dliqiso/pdel)

!**********************
else !RH < 1
!**********************
!Equilibrate rain/vapor
!**********************

   !apply evaporation fluxes
    dliqiso = irmass0 - Rr0*rmass
    irmass  = irmass0 - dliqiso
    ivmass  = ivmass0 + dliqiso/pdel

   !save copies of rain and vapor masses
    tmpv = vmass*pdel   !temporary storage variable for bulk vapor
    itmpv = ivmass*pdel !temporary storage variable for isotopic vapor
    itmpr = irmass      !temporary storage variable for isotopic rain 

  if((rmass .ge. 0._r8) .and. (vmass .ge. 0._r8) .and. &
     (itmpv .ge. 0._r8) .and. (itmpr .ge. 0._r8)) then !Ensure quantities are positive  (NOTE: may not be needed - JN).

    call wtrc_liqvap_equil(ae, fequil, tmpv, rmass, itmpv, itmpr, dliqiso)

    irmass = irmass + dliqiso
    ivmass = ivmass - dliqiso/pdel

  end if

!**********************
end if !RH < 1
!**********************
end subroutine stewart_isoevap


!=======================================================================
  subroutine wtrc_clear_precip(pstate, pbuf, itype)
!-----------------------------------------------------------------------
!
! Purpose: Clear the surface precipitation fields.
!
! Method:
!   Clear the accumulated amount of surface precipitation
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use physics_types,  only: physics_state, physics_ptend
  use water_types,    only: pwtype
  use physconst,      only: gravit, rhoh2o
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use constituents,   only: cnst_name
  
  type(physics_state), intent(in)    :: pstate                        ! State of the atmosphere
  type(physics_buffer_desc), pointer :: pbuf(:) !physics buffer

  integer, intent(in)                :: itype                         ! precipitation water type index
  
  integer                            :: iwset       ! water set index           
  integer                            :: lchnk
  integer                            :: ncol
  integer                            :: icol
  integer                            :: srfpcidx    ! Physics Buffer index
  real(r8), pointer, dimension(:)    :: srfpcp      ! Surface precipitation (m/s)

!-----------------------------------------------------------------------
!
  if (trace_water) then
    
    ! Use an internal local ptend, so we can do the mass checking.
    ncol  = pstate%ncol
    lchnk = pstate%lchnk
    
    ! Iterate over the water sets.
    do iwset = 1, wtrc_nwset
    
      ! The surface precipitation is stored in the physics buffer, so
      ! get a pointer to the data.
      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(itype,iwset), srfpcp)

      ! Calculate surface total.
      srfpcp(:ncol) = 0._r8
    end do
  end if

  return
end subroutine wtrc_clear_precip


!=======================================================================
  subroutine wtrc_collect_precip(pstate, pbuf, isrctype, idsttype, sediment_rates, dtime)
!-----------------------------------------------------------------------
!
! Purpose: Calculate the surface precipitation fields from the
!   corresponding 3D precipitation fields.
!
! Method:
!   Collects sedimenting cloud condensate from the bottom model level
!   into the surface precipitation fields. The sediment_rates indicates
!   the sedimentation rate from the bottom model level of the source
!   water type.
!
! Author: Chuck Bardeen
!
!  NOTE:  Subroutine is currently not used. - JN
!
!-----------------------------------------------------------------------
  use physics_types,  only: physics_state, physics_ptend
  use water_types,    only: pwtype
  use physconst,      only: gravit, rhoh2o
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use constituents,   only: cnst_name
  
  type(physics_state), intent(in)    :: pstate                        ! State of the atmosphere
  type(physics_buffer_desc), pointer :: pbuf(:) !physics buffer
  integer, intent(in)                :: isrctype                      ! source water type index
  integer, intent(in)                :: idsttype                      ! destination water type index
  real(r8), intent(in)               :: sediment_rates(pcols,pver,pwtype) ! Sedimentation rates (kg/kg/sec)
  real(r8), intent(in)               :: dtime                         ! length of timestep (s)
  
  integer                            :: iwset       ! water set index           
  integer                            :: lchnk
  integer                            :: ncol
  integer                            :: icol
  integer                            :: srfpcidx    ! Physics Buffer index
  real(r8), pointer, dimension(:)    :: srfpcp      ! Surface precipitation (m/s)

!-----------------------------------------------------------------------
!
  if (trace_water) then
    
    ! Use an internal local ptend, so we can do the mass checking.
    ncol  = pstate%ncol
    lchnk = pstate%lchnk
    
    ! Iterate over the water sets.
    do iwset = 1, wtrc_nwset
    
      ! The surface precipitation is stored in the physics buffer, so
      ! get a pointer to the data.
      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(idsttype,iwset), srfpcp) 
      
      ! The 3D data is in the state structure and/or the tendency. Calculate the mass and
      ! convert from kg/m2 to m/s.
      do icol = 1, ncol
      
        ! Calculate surface total.
        srfpcp(icol) = srfpcp(icol) + sediment_rates(icol, pver, isrctype) * &
                       pstate%pdel(icol,pver) / gravit / rhoh2o / dtime
      end do
    end do
  end if

  return
end subroutine wtrc_collect_precip


!=======================================================================
  subroutine wtrc_diagnose_precip(pstate, pmass, pbuf, iwtype, dtime)
!-----------------------------------------------------------------------
!
! Purpose: Calculate the surface precipitation fields from the
!   corresponding 3D precipitation fields.
!
! Method:
!   Currently, precipitation is a diagnostic field in CAM. While some
!   precipitation may occur from sedimentation of cloud liquid and ice
!   to surface, most of it comes from 3-D rain and snow fields that are
!   the result of microphysical processes (e.g. autoconversion) in the
!   atmosphere. CAM assumes that all snow and rain in the atmosphere
!   reaches the surface at the end of the time step. Thus 3-D fields for
!   water tracer rain and snow need to be collected up and placed into
!   the corresponding surface precipitation fields.
!
!   NOTE: This will just do one water type at a time, so that the timing
!   for the collection of convective precipitation can be different than
!   for stratiform precipitation. However precipitation from all
!   isotopic species will be collected separately.
!
! Author: Chuck Bardeen
!
!  NOTE:  It has been found that the sedimenting ice and cloud liquid that falls
!         through the bottom layer are added to the bulk stratiform precipitation, and thus
!         must be added here.  It does not appear that the sedimentation tendencies 
!         can be used directly in a simple way, as they involve the net flux of condensate
!         into the grid box, while the sediment that counts as precipitation is simply
!          the sediment flux out of the bottom grid box (Flux_out versus Flux_out-Flux_in)
!         There should be no fractionation involved with sedimentation, so only a ratio needs 
!         to be counted.  The routine below uses the pre-sedimentation bottom-level values to 
!         calculate the ratio, but this may not be ideal.  If the stratiform precipitation 
!         appears to have incorrect values, then this ratio may be the culprit. - JN
!
!
!-----------------------------------------------------------------------
  use physics_types,  only: physics_state, physics_ptend
  use water_types,    only: pwtype
  use physconst,      only: gravit, rhoh2o
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use constituents,   only: cnst_name  

   use time_manager,  only: get_nstep

  type(physics_state), intent(in)    :: pstate                        ! State of the atmosphere
  real(r8), intent(inout)            :: pmass(pcols,wtrc_nwset)       ! local constituents
  type(physics_buffer_desc), pointer :: pbuf(:)                       ! physics buffer
  integer,  intent(in)               :: iwtype                        ! water type inde
  real(r8), intent(in)               :: dtime                         ! length of timestep (s)

  integer                            :: iwset
  integer                            :: lchnk
  integer                            :: ncol
  integer                            :: icol
  integer                            :: srfpcidx      ! physics Buffer index 
  real(r8)                           :: R             ! water tracer precipitation ratio (unitless)
  real(r8), pointer, dimension(:)    :: srfpcp        ! surface precipitation (m/s)
  real(r8)                           :: stdpcp(pcols) ! bulk water precipitation <-used for mass fixer - JN

!-----------------------------------------------------------------------
!
  if (trace_water) then
    
    ! Use an internal local ptend, so we can do the mass checking.
    ncol  = pstate%ncol
    lchnk = pstate%lchnk
    
    !initalize variable:
    !stdpcp(:) = 0._r8 !<-used for mass fixer - JN

    ! Iterate over the water sets.
    do iwset = 1, wtrc_nwset
    
      ! The surface precipitation is stored in the physics buffer, so
      ! get a pointer to the data.
      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtype,iwset), srfpcp)
  
      ! The 3D data is in the state structure and/or the tendency. Sum the mass and
      ! convert from kg/m2 to m/s.
      do icol = 1, ncol
      
        ! Calculate surface total.
        srfpcp(icol) = srfpcp(icol) + pmass(icol,iwset)/gravit/rhoh2o/dtime

        !**********************************
        !Stratiform isotopic snow error fix -JN
        !**********************************
        !NOTE:  May not be needed. - JN
        !Save bulk (H2O) water for possible error-correcting.
        if(iwset .eq. 1) then
          stdpcp(icol) = srfpcp(icol)  
        end if

        !Prevent stratiform snow error (the random snow events that preciptate WAY too much isotopic mass).
        !NOTE:  Eventually the root cause of this issue should be found and fixed, but for now this will do. - JN
        if((iwset .ne. 1) .and. (srfpcp(icol) .gt. 2._r8 * stdpcp(icol))) then !isotopic precip too large, possible numerical error.
          !Write out warning to log file:
          if(stdpcp(icol) .gt. wtrc_qmin) write(*,*) 'ERROR: isotopic stratiform precipitation mass error!',srfpcp(icol),stdpcp(icol),iwtype,iwset,icol
          !Adjust water tracers back to standard:
          srfpcp(icol) = stdpcp(icol) !set isotopic mass to be equal to bulk water (aka destroy numerically produced mass).
        end if
        !***********************************

        ! Clear out the field.
        pmass(icol,iwset) = 0._r8     

      end do
    end do
  end if

  return
end subroutine wtrc_diagnose_precip


!=======================================================================
  subroutine wtrc_diagnose_bulk_precip(pstate, ptend, iwtype, dtime)
!-----------------------------------------------------------------------
!
! Purpose: 
!
! Method:
!   Currently, precipitation is a diagnostic field in CAM. If a consituent
!   was added in water tracer to mirror the bulk water, then clear this out
!   now. The actual diagnosis of the bulk precipitation is done in the 
!   CAM microphysics.
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use physics_types,  only: physics_state, physics_ptend
  use water_types,    only: pwtype
  
  type(physics_state), intent(in)    :: pstate                        ! State of the atmosphere
  type(physics_ptend), intent(inout) :: ptend                         ! State tendencies
  integer, intent(in)                :: iwtype                        ! water type inde
  real(r8), intent(in)               :: dtime                         ! length of timestep (s)
  
  integer            :: ncol

!-----------------------------------------------------------------------
!
  if (trace_water) then

    ! Does wtrc own some of the bulk fields?
    if (wtrc_add_stprecip) then
      ncol  = pstate%ncol

      ! Clear out the 3-D bulk precipitation field.
      ptend%lq(wtrc_bulk_indices(iwtype))        = .false.
      ptend%q(:ncol,:,wtrc_bulk_indices(iwtype)) = 0._r8
    end if
  end if
    
  return
end subroutine wtrc_diagnose_bulk_precip


  subroutine wtrc_output_precip(pstate, pbuf)
!-----------------------------------------------------------------------
!
! Purpose: Output the surface precipitation fields.
!
! Method:
!   Clear the accumulated amount of surface precipitation
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use constituents,   only: cnst_name
  use cam_history,    only: outfld
  
  type(physics_state), intent(in)    :: pstate                        ! State of the atmosphere
  type(physics_buffer_desc), pointer :: pbuf(:)     !physics buffer
  
  integer                            :: iwset       ! water set index           
  integer                            :: lchnk
  integer                            :: ncol
  integer                            :: srfpcidx    ! Physics Buffer index
  real(r8), pointer, dimension(:)    :: srfpcp      ! Surface precipitation (m/s)

!-----------------------------------------------------------------------
!
  ! Output fields for the isotopes.
  if (trace_water) then
    ncol  = pstate%ncol
    lchnk = pstate%lchnk
  
    ! Iterate over the water sets.
    do iwset = 1, wtrc_nwset
    
      ! Output the convective rain and snow.
      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvrain,iwset), srfpcp)
      call outfld('PRECRC_' // trim(cnst_name(wtrc_iawset(iwtcvrain, iwset))), srfpcp, pcols, lchnk)      

      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvsnow,iwset), srfpcp)
      call outfld('PRECSC_' // trim(cnst_name(wtrc_iawset(iwtcvsnow, iwset))), srfpcp, pcols, lchnk)      

      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstrain,iwset), srfpcp)   
      call outfld('PRECRL_' // trim(cnst_name(wtrc_iawset(iwtstrain, iwset))), srfpcp, pcols, lchnk)      

      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstsnow,iwset), srfpcp)
      call outfld('PRECSL_' // trim(cnst_name(wtrc_iawset(iwtstsnow, iwset))), srfpcp, pcols, lchnk)      

    end do
  end if

  return
end subroutine wtrc_output_precip


!=======================================================================
subroutine wtrc_mass_fixer(state)
!--------------------------------
!
!Purpose:  To reset the standard (H2O) water tracer back to the bulk water value, and adjust
!all of the proceeding water tracers bay an amount proportional to the H2O reset value.
!
!Written by:  Jesse Nusbaumer <nusbaume@colorado.edu> - April, 2014
!
!--------------------------------

!***************************
!use statements/header files
!***************************

use physics_types, only: physics_state

!*****************
!Delcare variables
!*****************

implicit none

!Input/Output Variables:
type(physics_state), intent(inout) :: state ! Physics state variables

!local variables:
integer  :: i,k,p,m !loop control variables
integer  :: ncol    !number of horizontal variables
real(r8) :: oval    !original value (used to calculate ratio).
real(r8) :: diff    !mass difference between standard tracer and bulk water
real(r8) :: R       !water tracer ratio

!**************
!set ncol value
!**************

ncol = state%ncol

!***************************
!Correct water tracer values
!***************************
do i = 1,ncol
  do k = 1,pver
    do p = 1,pwtype
      if(state%q(i,k,wtrc_iatype(1,p)) .ne. state%q(i,k,wtrc_bulk_indices(p))) then
        !calculate difference:
        diff = state%q(i,k,wtrc_iatype(1,p))-state%q(i,k,wtrc_bulk_indices(p))
        !save original value:
        oval = state%q(i,k,wtrc_iatype(1,p))
        do m=1,wtrc_nwset
          !calculate ratio:
          R = wtrc_ratio(iwspec(wtrc_iatype(m,p)),state%q(i,k,wtrc_iatype(m,p)),&
                         oval)
          !Adjust value:
          state%q(i,k,wtrc_iatype(m,p)) = state%q(i,k,wtrc_iatype(m,p)) - R*diff
        end do
      end if
    end do
  end do
end do

!****************************
!Correct water isotope values
!NOTE:  Some parts of the isotopic physics appear to be causing non-physical numerical errors.
!In order to fix these errors, H216O is corrected to match the bulk water, and all other water
!isotopes are adjusted accordingly. -JN
!****************************
if(wisotope) then
  do i = 1,ncol
    do k = 1,pver
      do p = 1,pwtype
        if(state%q(i,k,wtrc_iatype(2,p)) .ne. state%q(i,k,wtrc_bulk_indices(p))) then
          !calculate difference:
          diff = state%q(i,k,wtrc_iatype(2,p))-state%q(i,k,wtrc_bulk_indices(p))
          !save original value:
          oval = state%q(i,k,wtrc_iatype(2,p))
          ! RPF - restructured here to apply this second mass fixer loop
          ! to *ONLY* the isotopic species, not any H2O species.
          do m=2,wtrc_nwset
            if(iwspec(wtrc_iatype(m,p)) .ne. isph2o) then
              !calculate ratio:
              R = wtrc_ratio(iwspec(wtrc_iatype(m,p)),state%q(i,k,wtrc_iatype(m,p)),&
                           oval)
              !adjust value:
              state%q(i,k,wtrc_iatype(m,p)) = state%q(i,k,wtrc_iatype(m,p)) - R*diff
            end if
            ! end RPF changes
          end do
        end if
      end do
    end do
  end do
end if

end subroutine wtrc_mass_fixer


!=======================================================================
  function wtrc_check_h2o(testname, pstate, qloc, dtime, ptend)
!-----------------------------------------------------------------------
!
! Purpose: Check water tracers for conservation of mass.
!
! Method:
!   The tests are :
!
!     - (wtrc_check_total_h2o) - Checks whether the total mass of the first
!       water set in the model (which should be a total water set) is the
!       same as regular water.
!
!   Errors are reported if a negative total mass is
!   detected or if the relative change in mass is greater than
!   wtrc_qchkmin.

! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use physics_types,  only: physics_state, physics_ptend 
  use water_types,    only: pwtype, wtype_names
  use constituents,   only: cnst_name, pcnst

  
  character(len=*),intent(in)        :: testname                  ! name of the test
  type(physics_state), intent(in)    :: pstate                    ! initial state
  real(r8), intent(in)               :: qloc(pcols,pver,pcnst)    ! current state
  real(r8), intent(in)               :: dtime                     ! length of time step (s)
  type(physics_ptend), intent(inout), optional :: ptend           ! initial tendency
  
  logical                            :: wtrc_check_h2o
  
  integer            :: iwset
  integer            :: icol
  integer            :: ncol
  integer            :: itype
  integer            :: icnst
  integer            :: nwset
  real(r8)           :: bmass(pwtype)
  real(r8)           :: tmass(pwtype)
  real(r8)           :: sbmass
  real(r8)           :: stmass
  real(r8)           :: scale
 
!-----------------------------------------------------------------------
!
  wtrc_check_h2o = .true.
  
  if (trace_water) then

    ncol = pstate%ncol    

    ! Check to see if the total bulk mass is the same as the total wset mass.    
    if (wtrc_check_total_h2o) then
    
      ! When isotopes are used, only the first water set should have a
      ! mass comparable to bulk water.
      if (wisotope) then
        nwset = 1
      else
        nwset = wtrc_nwset
      end if
            
      do icol = 1, ncol
        bmass(:) = 0._r8
        
        do itype = 1, pwtype
          bmass(itype) = bmass(itype) + sum(pstate%q(icol, :, wtrc_bulk_indices(itype)) * pstate%pdel(icol, :))
          if (present(ptend)) then
            bmass(itype) = bmass(itype) + sum(ptend%q(icol,  :, wtrc_bulk_indices(itype)) * pstate%pdel(icol, :)) * dtime
          end if
        end do
        sbmass = sum(bmass)
        
        do iwset = 1, nwset
          tmass(:) = 0._r8
        
          ! Non-wisotope checks, used the fixed rstd to scale the total mass for
          ! the check.
          scale = 1._r8 / wtrc_get_rstd(iwspec(wtrc_iawset(1, iwset)))

          do itype = 1, pwtype
            tmass(itype) = tmass(itype) + sum(qloc(icol, :, wtrc_iawset(itype, iwset)) * pstate%pdel(icol, :))
          end do
          stmass = sum(tmass)
  
          if ((stmass < 0._r8) .or. ((sbmass) < 0._r8)) then
            wtrc_check_h2o = .false.
  
            write(iulog, *) trim(testname) // ": Warning - Negative total mass", iwset, icol, stmass, sbmass
          else if (sbmass > 0._r8) then  
              wtrc_check_h2o = .false.
              
              if (abs(((stmass * scale) - sbmass) / sbmass) > wtrc_qchkmin) then
  
              if (wtrc_warn_only) then
                write(iulog, *) trim(testname) // ": Warning - Total mass conservation limit exceeded ", iwset, icol, scale, stmass*scale, sbmass,(stmass*scale - sbmass) / sbmass
              else
                write(iulog, *) trim(testname) // ": Error - Total mass conservation limit exceeded ", iwset, icol, scale, stmass*scale, sbmass, (stmass*scale - sbmass) / sbmass
              end if
              
              if (wtrc_check_show_types) then
                do itype = 1, pwtype
                  write(iulog, *)  "  " // wtype_names(itype) // " :  ", tmass(itype)*scale, bmass(itype)
                end do
              end if
              
              if (.not. wtrc_warn_only) then
                call endrun()
              end if
            end if
          end if
        end do
      end do
    end if
  end if
  
  return
end function wtrc_check_h2o


!=======================================================================
  subroutine wtrc_adjust_h2o(testname, iwset, ncol, qloc, qmin)
!-----------------------------------------------------------------------
!
! Purpose: Adjust water values to avoid negative and/or small numbers
!
! Method:
!   The numerical techniques used here can create small rounding errors that
!   could generate QNEG3 errors. To compensate for this, mass with a wset
!   is adjust to make sure that all of the tracers are positive definite.
!
!   If there is not enough mass at the grid location for the wset to
!   become positive definite, then the run is stopped.
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use water_types,    only: pwtype, iwtice
  use constituents,   only: cnst_name, pcnst

  
  character(len=*),intent(in)        :: testname                  ! name of adjuster
  integer, intent(in)                :: iwset                     ! wset index
  integer, intent(in)                :: ncol                      ! number of columns
  real(r8), intent(inout)            :: qloc(pcols,pver,pcnst)    ! current state
  real(r8), optional, intent(in)     :: qmin                      ! a minimum value to use instead of wtrc_qmin
  
  integer            :: icol
  integer            :: itype
  integer            :: icnst
  integer            :: i
  integer            :: k
  integer            :: m
  real(r8)           :: qmin_int
 
!-----------------------------------------------------------------------

  if (trace_water) then

    do itype = 1, pwtype
      m = wtrc_iawset(itype, iwset)

      if (present(qmin)) then
        qmin_int = qmin
        
      ! The water types that represent particles have a different threshold
      ! than the gases.
      else if (itype /= iwtvap) then
        qmin_int = wtrc_qmin
      else
        qmin_int = 0._r8
      end if

      do k = 1, pver
        do i = 1, ncol
          if ((qloc(i, k, m) < qmin_int)) then
            qloc(i, k, m) = 0._r8
          end if
        end do
      end do
    end do
  end if
  
  return
end subroutine wtrc_adjust_h2o


!=======================================================================
  subroutine wtrc_adjust_h2o_bulk(testname, ncol, qloc)
!-----------------------------------------------------------------------
!
! Purpose: Adjust water values to avoid negative and/or small numbers
!
! Method:
!   The numerical techniques used here can create small rounding errors that
!   could generate QNEG3 errors. To compensate for this, mass with a wset
!   is adjust to make sure that all of the tracers are positive definite.
!
!   If there is not enough mass at the grid location for the wset to
!   become positive definite, then the run is stopped.
!
! Author: Chuck Bardeen
!
!-----------------------------------------------------------------------
  use water_types,    only: pwtype
  use constituents,   only: cnst_name, pcnst

  
  character(len=*),intent(in)        :: testname                  ! name of adjuster
  integer, intent(in)                :: ncol                      ! number of columns
  real(r8), intent(inout)            :: qloc(pcols,pver,pcnst)    ! current state
  
  integer            :: icol
  integer            :: itype
  integer            :: i
  integer            :: k
  integer            :: m
 
!-----------------------------------------------------------------------

  if (trace_water) then

    do itype = 1, pwtype
      m = wtrc_bulk_indices(itype)

      do k = 1, pver
        do i = 1, ncol

          if ((qloc(i, k, m) < 0._r8)) then
            qloc(i, k, m) = 0._r8
          end if
        end do
      end do
    end do
  end if
  
  return
end subroutine wtrc_adjust_h2o_bulk

!=======================================================================
subroutine wtrc_dicm( niso   , &
                      isp    , feq0   , told   , tnew   , &
                      vapold , liqold , iceold , vapent , &
                      vapnew , liqnew , icenew , vapdet , &
                      rainpr , snowpr , dliqmt , dicefz )
!-----------------------------------------------------------------------
!
! Differential (David's) Isotope Cloud Model (Mirophysics).
!
! This follows the type of model proposed by Merlivat and jouzel, with
! extension of mixed phase by Gendzelman and Arnold, and Ciais and
! Jouzel. Used in a plume model, the results should look very much like
! the Federer et al model. 
!
! Assume changes to total happen linearly.
! This routine is quite generic given some thought about the inputs.
!
!       dvap(m) + dliq(m) + dice(m) + prain(m) + psnow(m) = pvent - pvdet
!
!       dvap(m) = -pliq_vap(m) - pice_vap(m)                             - pvdet(m) + pvent(m)
!       dliq(m) =  pliq_vap(m)               - pice_liq(m) + pliq_ice(m) - prain(m) 
!       dice(m) =                pice_vap(m) + pice_liq(m) - pliq_ice(m) - psnow(m)
!
! For isotopes, there is also and exchange between rain and vapour)
!
! All exchanges are without fractionation EXCEPT:
!  o deposition of vapour onto ice as disillation (with kinetic
!  fractionation)
!  o vapour and liquid are held in equilibrim (equilibrium
!    fractionation,  unless evaporation then kinetic effects)
!  o vapour and rain tend toward equilibrium at some slow rate
!
! Assumes no rain or snow initially. If this is not the case, 
! you might judiciously dump all initial rain into cloud water, 
! and specify an effective precipitation rate (?).
! If this is REALLY a problem, email me, it can be fixed.
! [note to self: snow sublimates no frac, rain as distillation is
! easiest] Similarly, detrainment and entrainment is easily include...
! bot it all costs flops!
!
! In CAM, this module is suitable for shallow and deep convection, 
! and stratiform cloud condensation processes (although the latter
! needs some preprocessing to treat snow melt to rain, and
! evaporation of rain and snow).
!
! After 10 iterations the vapour and liquid delta values are good to the
! 5th significant figure of delta values (8 sig figs for mass). After
! 100 iterations the snow (which is the least accurate) is good to
! almost 3 sig figs delta. I recommend at least 10, and more if you can
! afford it. 20 or 30 seems to be a good choice.
! However, there is a better analytic shortcut that may have better
! convergence properties.

! REMBER: FEQ WILL BE APPLIED DIFFERENTIALLY, so ends up something strange
! (in fact it has an exponential form which can be composed to give the
! correct result...)
!
! Note, magnitude of  accumulated numerical imprecision increases with number of
! iterations. Mass imbalance of 1.e-13 after 1000000 iterations - not
! bad.
!
! David Noone <dcn@colorado.edu> - Wed Jul 21 18:07:18 MDT 2004
!
!-----------------------------------------------------------------------
!  use water_types,    only: wtype_get_alpha 
!-----------------------------------------------------------------------
  implicit none
!------------------------- Input Arguments -----------------------------
  integer , intent(in)    :: niso           ! number of isotope species
  integer , intent(in)    :: isp(niso)      ! species index

  real(r8), intent(in)    :: feq0           ! fraction fractionated
  real(r8), intent(in)    :: told           ! initial temperature
  real(r8), intent(in)    :: tnew           ! final temperature

  real(r8), intent(inout) :: vapold(niso)   ! initial vapour mass
  real(r8), intent(inout) :: liqold(niso)   ! initial cloud liquid mass
  real(r8), intent(inout) :: iceold(niso)   ! initial cloud ice mass
  real(r8), intent(inout) :: vapent(niso)   ! vapour entrainment

  real(r8), intent(inout) :: vapnew(niso)   ! final vapour mass
  real(r8), intent(inout) :: liqnew(niso)   ! final cloud liquid mass
  real(r8), intent(inout) :: icenew(niso)   ! final cloud ice mass
  real(r8), intent(inout) :: vapdet(niso)   ! vapour detrainment
  real(r8), intent(inout) :: rainpr(niso)   ! rain production 
  real(r8), intent(inout) :: snowpr(niso)   ! snow production
!
  real(r8), intent(in)    :: dliqmt         ! change in cloud liquid due to ice melt
  real(r8), intent(in)    :: dicefz         ! change in cloud ice due to liquid freeze

!------------------------- Local Variables -----------------------------

  integer m                         ! isotope number count
  integer itr                       ! iteration count
  integer nitr                      ! iterations for finite integration

  logical ldistdif                  ! differential/integral distillation

  real(r8) pliq_vap(niso)           ! prod. cloud liquid from vapour condensation
  real(r8) pice_vap(niso)           ! prod. cloud ice from vapour condensation
  real(r8) pice_liq(niso)           ! prod. cloud ice from cloud liquid freeze
  real(r8) pliq_ice(niso)           ! prod. cloud liquid from cloud ice melt
  real(r8) prnw_vap(niso)           ! prod. rain water from vapour condensation exchange

  real(r8) pvent(niso)              ! prod. vapour mass by entrainment
  real(r8) pvdet(niso)              ! prod. detrained vapour mass
  real(r8) psnow(niso)              ! prod. snow from conversion of cloud ice
  real(r8) prain(niso)              ! prod. rain from conversion of cloud liq

  real(r8) dvap(niso)               ! change in vapour
  real(r8) dliq(niso)               ! change in cloud liquid
  real(r8) dice(niso)               ! change in cloud ice

  real(r8) qvap(niso)               ! vapour integrand
  real(r8) qliq(niso)               ! cloud liquid integrand
  real(r8) qice(niso)               ! cloud ice integrand
  real(r8) vdet(niso)               ! vapour detrainment accumulator
  real(r8) rain(niso)               ! rain accumulator
  real(r8) snow(niso)               ! snow accumulator

  real(r8) alpice(niso)             ! fractionation factor for ice (kinetic)
  real(r8) alpliq(niso)             ! fractionation factor for liquid

  real(r8) totalold, totalnew       ! totals for conservation checks
  real(r8) vtotnew, visonew         ! updates from distillation to ice
  real(r8) dvapiso                  ! change in vapour from distillation
  real(r8) tk                       ! temperature (kelvin)
  real(r8) fint                     ! fraction through integration

!-----------------------------------------------------------------------
  real(r8)  :: feq_rn = 0.0_r8       ! fraction rain vapour equulibration
!  real(r8)  :: feq_rn = 0.1       ! fraction rain vapour equulibration [VERY SENSITIVE]
!  real(r8)  :: fieq = 0.          ! fraction of initial equilibration (slow convergence)
!  real(r8)  :: fieq = 1.          ! fraction of initial equilibration (better)
  real(r8)  :: fieq = 0.5_r8         ! fraction of initial equilibration (like extra iteration)
!-----------------------------------------------------------------------
  real(r8)  :: qtolerr = 1.e-12   ! small q for error checks
  real(r8)  :: qtiny   = 1.e-22   ! essentially zero q
!-----------------------------------------------------------------------
   alpliq(1) = 1._r8
   alpice(1) = 1._r8
!
! Check input budget
!
   totalold = vapold(1) + liqold(1) + iceold(1) + vapent(1)
   totalnew = vapnew(1) + liqnew(1) + icenew(1) + vapdet(1) + rainpr(1) + snowpr(1)
   if (abs(totalold - totalnew) > qtolerr) then
      write(iulog,*) '(wtrc_dicm) WARNING: total budget does not balance: old /= new'
      write(iulog,2) 'old',vapold(1), liqold(1), iceold(1), vapent(1)
      write(iulog,2) 'new',vapnew(1), liqnew(1), icenew(1), vapdet(1), rainpr(1), snowpr(1)
      write(iulog,2) 'total',totalold, totalnew, totalold - totalnew
!!      call endrun('(wiso_dicm) ABORTED.')
   end if
!
! Check things we assume (all quantities positive definite, and temperature reasonable)
!  
   do m = 1, niso
     if (vapold(m) < 0._r8 .or. liqold(m) < 0._r8 .or. iceold(m) < 0._r8) then
       write(iulog,*) '(wtrc_dicm) old values < 0.: m=',m
       write(iulog,2) 'old:',vapold(m),liqold(m),iceold(m)
       write(iulog,2) 'new:',vapnew(m),liqnew(m),icenew(m)
       vapold(1) = max(vapold(1),0._r8)
       liqold(1) = max(liqold(1),0._r8)
       iceold(1) = max(iceold(1),0._r8)
       call endrun('(wtrc_dicm) Abort')
     end if
     if (vapent(m) < 0._r8) then
        write(iulog,*) '(wtrc_dicm) entrainment < 0. : m=',m
        write(iulog,*) 'vapent:',vapent(m)
        vapent(m) = max(vapent(m),0._r8)
     end if
  end do
  if (vapnew(1) < 0._r8 .or. liqnew(1) < 0. .or. icenew(1) < 0.) then
       write(iulog,*) '(wtrc_dicm) new values < 0.'
       write(iulog,2) 'old:',vapold(1),liqold(1),iceold(1)
       write(iulog,2) 'new:',vapnew(1),liqnew(1),icenew(1)
       vapnew(1) = max(vapnew(1),0._r8)
       liqnew(1) = max(liqnew(1),0._r8)
       icenew(1) = max(icenew(1),0._r8)
     call endrun('(wtrc_dicm) Abort')
  end if
  if (vapdet(1) < 0._r8) then
     write(iulog,*) '(wtrc_dicm) detrainment < 0. : m=',m
     write(iulog,*) 'vapdet:',vapdet(1)
     vapdet(1) = max(vapdet(1),0._r8)
  end if
  if (rainpr(1) < 0._r8 .or. snowpr(1) < 0._r8) then
     write(*,*) '(wtrc_dicm) precipitation production < 0.'
     write(*,*) 'rain, snow:',rainpr(1),snowpr(1)
     call endrun('(wtrc_dicm) Abort...')
     rainpr(1) = max(rainpr(1), 0._r8)
     snowpr(1) = max(snowpr(1), 0._r8)
  end if

  if (told > 330_r8 .or. tnew > 330._r8  .or. &
      told < 130_r8 .or. tnew < 130._r8 ) then
     write(*,*) '(wtrc_dicm) input temperature unreasonable:',told,tnew
!!     call endrun('(wiso_dicm) Abort')
  end if
 2 format(a6,10e16.6)
!
! Given input net changes to total waters, compose small finte differences.
! These are held constant over all iterations.
!  
   dvap(1)  = (vapnew(1) - vapold(1))
   dliq(1)  = (liqnew(1) - liqold(1))
   dice(1)  = (icenew(1) - iceold(1))
!
! Work out what's going to happen, and decide if we need to iterate.
! If iterating use differential form for distillation, as it is slightly faster.
!
  nitr = wtrc_citer             ! default full iterations
  ldistdif = .true.             ! default full differential 
!
  if (abs(dvap(1))   < qtiny .and. &
      abs(dliq(1))   < qtiny .and. &
      abs(dice(1))   < qtiny .and. &
      abs(rainpr(1)) < qtiny .and. &
      abs(snowpr(1)) < qtiny .and. &
      abs(dliqmt)    < qtiny .and. &
      abs(dicefz)    < qtiny) then              ! nothing to do, just one pass
    nitr = 1
    ldistdif = .false.

  else if (iceold(1) < qtiny .and. icenew(1) < qtiny .and.  &
           snowpr(1) < qtiny .and. rainpr(1) < qtiny .and.  &
           abs(dliqmt)<qtiny .and. abs(dicefz)<qtiny ) then     ! equilibrate vap/liq
    nitr = 1

  else if (liqold(1) < qtiny .and. liqnew(1) < qtiny .and.  &
           snowpr(1) < qtiny .and. rainpr(1) < qtiny .and.  &
           abs(dliqmt)<qtiny .and. abs(dicefz)<qtiny ) then     ! distil to ice
    nitr = 1
    ldistdif = .false.
  end if
!
! Convert changes to iterative increments
!
   pvent(:niso) = vapent(:niso) / real(nitr)
!
   dvap(1)  = dvap(1)   / real(nitr)  !<-This variable isn't used? - JN
   dliq(1)  = dliq(1)   / real(nitr)
   dice(1)  = dice(1)   / real(nitr)

   pvdet(1) = vapdet(1) / real(nitr)
   prain(1) = rainpr(1) / real(nitr)
   psnow(1) = snowpr(1) / real(nitr)

   pliq_ice(1) = dliqmt / real(nitr)
   pice_liq(1) = dicefz / real(nitr)
!
! Given the amount of cloud liquid/ice melt freeze, 
! we can back out the other terms
!

pliq_vap(1) = dliq(1) + prain(1) + pice_liq(1) - pliq_ice(1)
pice_vap(1) = dice(1) + psnow(1) - pice_liq(1) + pliq_ice(1)
prnw_vap(1) = 0._r8  

!
! Copy inputs to local variables to iterate on
!  
   qvap(:niso) = vapold(:niso)
   qliq(:niso) = liqold(:niso)
   qice(:niso) = iceold(:niso)
   vdet(:niso) = 0._r8
   rain(:niso) = 0._r8
   snow(:niso) = 0._r8

!
! Integrate over small changes, for all isotopes
! Here order is important... if number of iterations are small.
! The iteration is done in two substeps to improve accuracy 
! (although I'm not convinced it's much better than doubling the number
!  of iterations!)
!
   do itr = 1, nitr
!
! Compute temperature, and assign fractionation factors
!
     fint = (real(itr)-0.5_r8) / real(nitr)
     tk = fint*tnew + (1._r8-fint)*told
!
     !Assume 100% RH:
     do m = 2, niso
         alpliq(m) = wtrc_get_alpha(qvap(m),tk,isp(m),iwtvap,iwtliq,.false.,1._r8,.false.)
         alpice(m) = wtrc_get_alpha(qvap(m),tk,isp(m),iwtvap,iwtice,.false.,1._r8,.true.)
     end do
!
! Start by entraining vapour from the environment
!
     do m = 1, niso
       qvap(m) = qvap(m) + pvent(m)
     end do
!
! Equilibrating liquid and vapour, for a fraction of iteration vapour/liquid production
!

     if (fieq > 0._r8) then
       qliq(1) = qliq(1) + fieq*pliq_vap(1)
       qvap(1) = qvap(1) - fieq*pliq_vap(1)
       do m = 2, niso
         call wtrc_liqvap_equil(alpliq(m),feq0,qvap(1),qliq(1), &
                      qvap(m),qliq(m),pliq_vap(m))
       end do
     end if

!
! Compute the ice melt and liquid freeze, and move it
!
     do m = 2, niso
        pice_liq(m) = pice_liq(1)*wtrc_ratio(isph2o,qliq(m),qliq(1))
        pice_liq(m) = min(pice_liq(m), qliq(m))

        pliq_ice(m) = pliq_ice(1)*wtrc_ratio(isph2o,qice(m),qice(1))
        pliq_ice(m) = min(pliq_ice(m), qice(m))
     end do
!
     do m = 1, niso
       qice(m) = qice(m) + pice_liq(m) - pliq_ice(m)
       qliq(m) = qliq(m) - pice_liq(m) + pliq_ice(m)
     end do
!
! Compute the ice/vapour deposition/sublimation
! Both integral and differential schemes have about the same convergence
! characteristics.
     if (ldistdif) then         !  use differential form for distillation 
       do m = 2, niso
         if (pice_vap(1) > 0._r8) then ! ice deposition by distillation, kinetic
            pice_vap(m) = alpice(m)*pice_vap(1)*wtrc_ratio(isph2o,qvap(m), qvap(1))
            pice_vap(m) = min(pice_vap(m),  qvap(m))
         else                       ! sublimation, no fractionation
            pice_vap(m) = pice_vap(1)*wtrc_ratio(isph2o,qice(m), qice(1))
            pice_vap(m) = max(pice_vap(m), -qice(m))
         end if
       end do
     else                       ! use integral form for distillation
       if (pice_vap(1) > 0._r8) then   ! ice deposition by distillation, kinetic
         vtotnew = qvap(1) - pice_vap(1)
         vtotnew = max(vtotnew, 0._r8)
         do m = 2, niso
           call wtrc_vap_distil(alpice(m),qvap(1),vtotnew,qvap(m),visonew,dvapiso)
           pice_vap(m) = -dvapiso
           pice_vap(m) = min(pice_vap(m),  qvap(m))
         end do
       else                         ! sublimation, no fractionation
         do m = 2, niso
           pice_vap(m) = pice_vap(1)*wtrc_ratio(isph2o,qice(m), qice(1))
           pice_vap(m) = max(pice_vap(m), -qice(m))
         end do
       end if
     end if

     do m = 1, niso
       qice(m) = qice(m) + pice_vap(m)
       qvap(m) = qvap(m) - pice_vap(m)
     end do

!
! equilibrate liquid/vapour equilibrium with second increment
!
     qliq(1) = qliq(1) + (1._r8-fieq)*pliq_vap(1)
     qvap(1) = qvap(1) - (1._r8-fieq)*pliq_vap(1)
     do m = 2, niso
       call wtrc_liqvap_equil(alpliq(m),feq0,qvap(1),qliq(1), &
                    qvap(m),qliq(m),pliq_vap(m))
     end do
!
! Detrain vapour
!
     do m = 2, niso
       pvdet(m) = pvdet(1)*wtrc_ratio(isph2o,qvap(m),qvap(1))
       pvdet(m) = min(pvdet(m), qvap(m))
     end do
!
     do m = 1, niso
       qvap(m) = qvap(m) - pvdet(m)
       vdet(m) = vdet(m) + pvdet(m)
     end do
!
! Complete the iteration by removing rain and snow production
!
     do m = 2, niso             ! must do m=1 first
       prain(m) = prain(1)*wtrc_ratio(isph2o,qliq(m),qliq(1))
       prain(m) = min(prain(m), qliq(m))
       psnow(m) = psnow(1)*wtrc_ratio(isph2o,qice(m),qice(1))
       psnow(m) = min(psnow(m), qice(m))
     end do
!
     do m = 1, niso
       qliq(m) = qliq(m) - prain(m)
       rain(m) = rain(m) + prain(m)
       qice(m) = qice(m) - psnow(m)
       snow(m) = snow(m) + psnow(m)
     end do
!
! Allows some fraction of the rain to equilibrate with cloud vapour
! (isotope only exchange effect)
!
     if (feq_rn > 0._r8) then
       do m = 2, niso
         call wtrc_liqvap_equil(alpliq(m),feq_rn,qvap(1),rain(1), &
                      qvap(m),rain(m),prnw_vap(m))
       end do
     end if

!
  end do                ! itr, iteration loop
!
! For final state, ensure liquid and vapour are in equilibrium
! (might be needed in calling code)
!

  do m = 2, niso
    call wtrc_liqvap_equil(alpliq(m),feq0,qvap(1),qliq(1), &
                 qvap(m),qliq(m),pliq_vap(m))
  end do

!
! Assin final values to output, and check budget was done correctly
!

!write(iulog,*) 'bulk water vapor',vapnew(1),qvap(1)
!write(iulog,*) 'bulk water liquid',liqnew(1),qliq(1)
!write(iulog,*) 'bulk water ice', icenew(1),qice(1)

  do m = 2, niso        ! dont assign to total
!
    vapnew(m) = qvap(m)
    liqnew(m) = qliq(m)
    icenew(m) = qice(m)
    vapdet(m) = vdet(m)
    rainpr(m) = rain(m)
    snowpr(m) = snow(m)
!
! Final check for budgets 
! If not balanced, try more iterations
!  Also, could do an adjustment to enforce mass conservation while
!  preserving ratios...
!
    totalold = vapold(m) + liqold(m) + iceold(m) + vapent(m)
    totalnew = vapnew(m) + liqnew(m) + icenew(m) + vapdet(m) + rainpr(m) + snowpr(m)
    if (abs(totalold - totalnew) > qtolerr) then
       write(*,*) '(wtrc_dicm) WARNING - isotope budget not balanced.'
       write(*,*) m,totalold, totalnew, totalold - totalnew
!!       call endrun('(wiso_dicm) ABORTED.')
    end if
!
  end do

  return
end subroutine wtrc_dicm

!=======================================================================
subroutine wtrc_liqvap_equil(alpha, feq0, vaptot, liqtot, vapiso, liqiso, dliqiso)
!-----------------------------------------------------------------------
! Equilibrate vapour and liquid, assuming the input mass is right
! but distribution is not. liquid vapour partitioning is given by input
! total water. 
!
! Given budgets are unchanged   : v + l = q, and vi + li = qi
! And equilibrium of final state: li/l = alpha vi/v
! Solve for final isotope vapour (thus liquid from budget):
!  vi = efac*qi
!    efac = 1/F; F = alpha*(l/v) + 1
!   
! Nothing fancy here.
!   
! David Noone <dcn@colorado.edu> - Wed Jul 21 19:08:54 MDT 2004
!
!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
  real(r8), intent(in)      :: alpha        ! fractionation fatcor 
  real(r8), intent(in)      :: feq0         ! fraction fractionated
  real(r8), intent(in)      :: vaptot       ! total vapour
  real(r8), intent(in)      :: liqtot       ! total liquid
  real(r8), intent(inout)   :: vapiso       ! isotopic vapour
  real(r8), intent(inout)   :: liqiso       ! isotopic liquid
  real(r8), intent(out)     :: dliqiso      ! change in liquid
!------------------------- Local Variables -----------------------------
  real(r8) dviso                            ! change in vapour
  real(r8) qtot, qiso                       ! total mass of total and isotope
  real(r8) efac                             ! equilibration factor
  real(r8) ratio
!-----------------------------------------------------------------------
!  real(r8) :: qtiny = 1.e-22
  real(r8) :: qtiny = 1.e-36
!-----------------------------------------------------------------------
!
  dliqiso = 0._r8
  qtot = vaptot + liqtot                ! not used
  qiso = vapiso + liqiso

!
! Check for trivial cases - while these fall out in the algebra, 
! the numericas can be a little picky... so treat them explictly.
!
! NOTE:  This code can produce incorrect values if the input is a moisture
!flux with negative values.  Thus it must be guarranteed that every input
!into this subroutine is positive definite - JN.
!
  if (qtot < qtiny) then                ! makes no sense, do nothing
!     write(iulog,*) '(wiso_liqvap_equil) no (trivial) water - doing nothing'
     return
  end if
  if (qiso < qtiny) then
!     write(iulog,*) '(wiso_liqvap_equil) no (trivial) isotope - doing nothing'
     return
  end if
!
  if (liqtot < qtiny) then
!     write(iulog,*) '(wiso_liqvap_equil) no liquid - dump all isotope to vapor'
     dliqiso = -liqiso
     vapiso = vapiso - dliqiso
     liqiso = 0._r8
     return
  end if
  if (vaptot < qtiny) then
!     write(iulog,*) '(wiso_liqvap_equil) no vapour - dump all isotope to liquid'
     dliqiso = vapiso
     vapiso = 0._r8
     liqiso = liqiso + dliqiso
     return
  end if
!
! Compute mass exchange from Rc = a Rv and v + l + q, in the normal way
!
  dviso = wtrc_dqequil(alpha,feq0,vaptot,liqtot,vapiso,liqiso)
  dliqiso = -dviso

!
! Update for outpout
!
  liqiso = liqiso + dliqiso
  vapiso = vapiso - dliqiso
!
  return
end subroutine wtrc_liqvap_equil

!=======================================================================
function wtrc_dqequil(alpha,feq0,vtotnew,ltotnew,visoold,lisoold)
!-----------------------------------------------------------------------
! Computes changes in vapour isotope due to equilibration (optionally
! partial), in a two-phase system. Nothing fancy here:
! Solve the mass balance under constraints dv+dl = 0; Rl = aRv
! Author: David Noone <dcn@colorado.edu> - Wed Aug 11 13:32:00 MDT 2004
!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
  real(r8), intent(in)  :: alpha      ! fractionation factor
  real(r8), intent(in)  :: feq0       ! fraction equilibrated
  real(r8) ,intent(in)  :: vtotnew    ! new vapour
  real(r8) ,intent(in)  :: ltotnew    ! new liquid
  real(r8) ,intent(in)  :: visoold    ! old isotope vapour
  real(r8) ,intent(in)  :: lisoold    ! old isotope liquid
!
  real(r8) wtrc_dqequil               ! return value
!-----------------------------------------------------------------------
  real(r8) qiso                       ! total isotope
  real(r8) vieql                      ! new isotope vapour with equilibration
  real(r8) vinof                      ! new isotope vapour no fractionation
  real(r8) visonew                    ! new isotope vapour
  real(r8) dviso                      ! change in isotope vapour
!-----------------------------------------------------------------------
!
  qiso = visoold + lisoold

! fractionating: Rc = alpha Rv
  vieql = qiso * wtrc_efac(alpha, vtotnew, ltotnew)

! non-fractionating: Rc = Rv
  vinof = qiso * wtrc_efac(1.0_r8  , vtotnew, ltotnew)

! Merge for partial equilibration
  visonew = feq0*vieql + (1._r8-feq0)*vinof

! Compute tendency
  dviso = visonew - visoold

! Check for numerical overflow
  if (dviso < 0._r8) then
    dviso = max(dviso , -visoold )   ! don't move more than available vapour
  else
    dviso = min(dviso ,  lisoold )   ! don't move more than available liquid
  end if
! Assign to output
  wtrc_dqequil = dviso
end function wtrc_dqequil

!=======================================================================
subroutine wtrc_vap_distil(alpha,vtotold,vtotnew,visoold,visonew,dvapiso)
!-----------------------------------------------------------------------
! Performs a rayleigh distillation on some vapour increment assuming
! the parcel is in isolation (this is in integral form). 
!
! Given, dvi/dv = alpha (vi/v), integration from vold to vnew gives:
!     vinew = viold*(vnew/vold)**alpha
!
! Nothing fancy here.
!
! David Noone <dcn@colorado.edu> - Thu Jul 22 11:18:29 MDT 2004
!
!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
  real(r8), intent(in)      :: alpha        ! fractionation factor
  real(r8), intent(in)      :: vtotold      ! initial total vapour
  real(r8), intent(in)      :: vtotnew      ! final total vapour
  real(r8), intent(in)      :: visoold      ! initial isotope vapour
  real(r8), intent(out)     :: visonew      ! initial isotope vapour
  real(r8), intent(out)     :: dvapiso      ! change in isotope vapour (diagnostic)
!-----------------------------------------------------------------------
  real(r8) :: qtiny = 1.e-22
!-----------------------------------------------------------------------
  dvapiso = 0._r8
!
! Check for trivial cases
!
  if (vtotold < vtotnew) then
    write(*,*) '(wtrc_vap_distill) vapour increase, while expecting decrease'
!   call endrun('ABORT')
  endif
  if (vtotold < qtiny) then
    write(*,*) '(wtrc_vap_distil) no (trivial) vapour - doing nothing.'
    visonew = visoold
    return
  end if
!
! Solve the logarithm, and back out change
!
  visonew = visoold*wtrc_ratio(isph2o,vtotnew,vtotold)**alpha
  dvapiso = visonew - visoold
  return
end subroutine wtrc_vap_distil

!=======================================================================
subroutine wtrc_precip_evap(state,rprd,deltat,evpbulk,subbulk,dqdt,prec,snow)
!
!--------------------------------------
!
!Purpose:  To calculate the influence of
!precipitation phase changes (evaporation,
!sublimation, melting, and freezing) on water 
!tracers and water isotopes.
!
!Written by:  Jesse Nusbaumer <nusbaume@colorado.edu> - April 2012
!
!--------------------------------------


use shr_kind_mod,          only: r8 => shr_kind_r8
use physics_types,         only: physics_state
use ppgrid,                only: pcols, pver, pverp
use physconst,             only: gravit, tmelt
use cloud_fraction,        only: cldfrc_fice
use wv_saturation,         only: qsat
use water_types,           only: iwtvap

!*****************
!Declare variables 
!*****************

implicit none

!Input/Output:

type(physics_state), intent(in) :: state          ! Physics state variabiles

real(r8), intent(in)    :: rprd(pcols,pver,pcnst) ! precip production (kg/kg/s)
real(r8), intent(in)    :: deltat                 ! time step

real(r8), intent(in)    :: evpbulk(pcols,pver)    ! bulk evaporation rate (kg/kg/s)
real(r8), intent(in)    :: subbulk(pcols,pver)    ! bulk sublimation rate (kg/kg/s)

real(r8), intent(inout) :: dqdt(pcols,pver,pcnst) ! Change in water vapor (kg/kg/s)
real(r8), intent(inout) :: prec(pcols,pcnst)      ! Conv.-scale precip. rate (m/s)
real(r8), intent(inout) :: snow(pcols,pcnst)      ! Conv.-scale snowfall rate (m/s)

!Local variables:

real(r8) :: rnbulk(pcols,pver)                  !bulk rain evaporation rate (prec-snow) (kg/kg/s)

real(r8) :: evp(pcols,pver,wtrc_ntype(iwtvap))  !Water tracer evaporation (kg/kg/s)
real(r8) :: mlt(pcols,pver,wtrc_ntype(iwtvap))  !Water tracer snow melt   (kg/kg/s)
real(r8) :: frz(pcols,pver,wtrc_ntype(iwtvap))  !Water tracer rain freeze (kg/kg/s)
real(r8) :: sub(pcols,pver,wtrc_ntype(iwtvap))  !Water tracer snow sublimation (kg/kg/s)

real(r8) :: flxpr(pcols,pverp,wtrc_ntype(iwtvap)) !precip flux  (kg/m2/s)
real(r8) :: flxsn(pcols,pverp,wtrc_ntype(iwtvap)) !Snow flux  (kg/m2/s)
real(r8) :: dldt(pcols,pver,wtrc_ntype(iwtvap))   !Rain tendency (kg/kg/s)
real(r8) :: dsdt(pcols,pver,wtrc_ntype(iwtvap))   !Snow tendency

real(r8) :: Rp(pcols,pver)      !Water tracer precip ratio
real(r8) :: Rr(pcols,pver)      !Water tracer rain ratio
real(r8) :: Rv(pcols,pver)      !Water tracer vapor ratio
real(r8) :: Rs(pcols,pver)      !Water tracer snow ratio <-Not used?
real(r8) :: fice(pcols,pver)    !Fraction of ice in layer
real(r8) :: fsnow(pcols,pver)   !Fraction that becomes snow
real(r8) :: est(pcols,pver)     !Saturation vapor pressure
real(r8) :: qst(pcols,pver)     !Saturation specific humidity
real(r8) :: rh(pcols,pver)      !Relative humidity
real(r8) :: heff(pcols,pver)    !Effective humidity

real(r8) :: alpliq              !Equilibrium liquid-vapor fractionation factor
real(r8) :: alpkin              !Kinetic fractionation factor
real(r8) :: phi                 !Tunable parameter to determine heff
real(r8) :: ivtmp               !Temporary water tracer vapor value
real(r8) :: iltmp               !Temporary water tracer rain value
real(r8) :: vtmp                !Temporary standard vapor value
real(r8) :: ltmp                !Temporary standard liquid value
real(r8) :: dliqiso             !Change in isotopic liquid to due equilibration
real(r8) :: Re                  !Rain evaporation ratio (includes fractionation)

!For Stewart rain re-evaporation:
real(r8) :: totrnfx(pcols,pverp) !Total precipitation flux
real(r8) :: gam,bet              !Stewart equation parameters
real(r8) :: fr                   !Fraction of precip remaining after rain evaporation
real(r8) :: Rstw                 !precipitation ratio after application of Stewart equation

!These next two variables are present in water_isotopes.F90, and thus
!ideally should be pulled from there. - JN

real(r8) :: difrm(pwtspec)       !Diffusivity ratio array
real(r8) :: dkfac                !Diffusivity power law factor

!DEBUGGING:
integer     :: ispec                    !water isotope species 
real(r8)    :: dz(pcols,pver)           !layer thickness in height
real(r8)    :: fequil                   !fraction equilibrated (unitless)
!real(r8), parameter :: radius=0.001_r8  !assumed diamaeter of raindrop (m) 
real(r8)    :: radius                   !assumed radius of raindrop (m)

!Index and loop variables:

integer :: ncol               !Number of atmospheric columns

integer :: i,k,m              !loop variables

!For precip mass fixer:
real(r8) :: rdiff
real(r8) :: sdiff
real(r8) :: rmass0
real(r8) :: smass0
real(r8) :: Rd

!********************
!Initialize variables
!********************

!Set variables to zero:
evp(:,:,:)   = 0._r8
mlt(:,:,:)   = 0._r8
frz(:,:,:)   = 0._r8
sub(:,:,:)   = 0._r8
fice(:,:)    = 0._r8
fsnow(:,:)   = 0._r8
rnbulk(:,:)  = 0._r8
flxpr(:,:,:) = 0._r8
flxsn(:,:,:) = 0._r8
totrnfx(:,:) = 0._r8
dz(:,:)      = 0._r8

ncol = state%ncol

!Determine ice fraction in layer
call cldfrc_fice(ncol, state%t, fice, fsnow)

!Determine saturation vapor pressure:
call qsat(state%t(1:ncol, 1:pver), state%pmid(1:ncol, 1:pver), &
           est(1:ncol, 1:pver), qst(1:ncol, 1:pver))

!calculate relative humidity:
rh(:ncol,:) = state%q(:ncol,:,1)/qst(:ncol,:)

!Diffusivity ratios and powers (should be called from water_isotopes.F90 eventually):
difrm = (/ 1._r8, 1._r8, 0.9757_r8, 0.9727_r8 /)         ! Merlivat 1978 (direct from paper)
dkfac = 0.58_r8                                          ! From Stewart, 1975?

!phi value given in Bony et. al., 2008:
phi = 0.9_r8

!calculate layer thickness:
do i = 1,ncol
  do k = 1,pver
    dz(i,k) = state%zi(i,k) - state%zi(i,k+1)
  end do
end do

!**********************************
!Calculate bulk rain re-evaporation
!**********************************

rnbulk(:,:) = evpbulk(:,:)-subbulk(:,:)

!**************************
!loop through water tracers
!**************************

!do m=1,wtrc_ntype(iwtvap)

Rp(:,:) = 0._r8 !initalize ratios
Rv(:,:) = 0._r8

!************************
!Loop through grid points
!************************

do i=1,ncol
  do k=1,pver

   do m=1,wtrc_ntype(iwtvap)

   Rp(:,:) = 0._r8 !initalize ratios
   Rr(:,:) = 0._r8
   Rs(:,:) = 0._r8
   Rv(:,:) = 0._r8

!**********
!Set ratios
!**********

!NOTE:  Evaporation and sublimation are functions of the precip flux, not the precip
!production.  Thus the ratios should be based off of the precip flux. - JN

!Water tracer precipitation ratio:
if(k .eq. 1) then !at top of column?
  Rp(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(1,iwtvap)))
else
  if(flxpr(i,k,1) .ne. 0._r8) then !is there actual precipitation here?
    Rp(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),flxpr(i,k,m),flxpr(i,k,1))
  else !If not, just use local production.
    Rp(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(1,iwtvap)))
  end if
end if

!Water tracer rain ratio:
if(k .eq. 1) then !at top of column?
  Rr(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(1,iwtvap)))
else
  if((flxpr(i,k,1)-flxsn(i,k,1)) .gt. 0._r8) then !is there actual precipitation here (if negative, then assume all snow)?
    Rr(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),(flxpr(i,k,m)-flxsn(i,k,m)),(flxpr(i,k,1)-flxsn(i,k,1)))
  else !If not, just use local production.
    Rr(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(1,iwtvap)))
  end if
end if

!Water tracer snow ratio:
if(k .eq. 1) then !at top of column?
  Rs(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(1,iwtvap)))
else
  if(flxsn(i,k,1) .ne. 0._r8) then !is there actual precipitation here?
    Rs(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),flxsn(i,k,m),flxsn(i,k,1))
  else !If not, just use local production.
    Rs(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(m,iwtvap)),rprd(i,k,wtrc_iatype(1,iwtvap)))
  end if
end if

!Water tracer vapor ratio:
Rv(i,k) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),state%q(i,k,wtrc_iatype(m,iwtvap)),state%q(i,k,wtrc_iatype(1,iwtvap)))

!**************************
!Calculate rain evaporation
!**************************

!Stewart, 1975 isotopic evaporation formulation:

if(wisotope .and. (m .gt. 1) .and. (totrnfx(i,k) .gt. 0._r8)) then !Are water isotopes being used and not H2O or H216O, and is precip. present?
  !calculate equilibrium alpha:
  ispec = iwspec(wtrc_iatype(m,iwtvap)) !water isotope species
  !fractionation factor:
  alpliq = wtrc_get_alpha(state%q(i,k,wtrc_bulk_indices(iwtvap)),state%t(i,k),ispec,&
                          iwtvap,iwtliq,.false.,1._r8,.false.)

  !calculate kinetic alpha:
   alpkin = (1._r8/difrm(ispec))**dkfac
!   alpkin = 1._r8 !turn kinetic fractionation off...

  !calculate effective humidity:
  heff(i,k) = phi+((1_r8-phi)*rh(i,k)) !From Bony et. al., 2008

  !remove evaporated mass in specific grid level:
  fr = (totrnfx(i,k)-(rnbulk(i,k)*state%pdel(i,k)/gravit))/totrnfx(i,k)

  !Constrain ratio fraction:
  fr = min(1._r8,max(0._r8,fr))

  !calculate radius based off rain rate:
  !NOTE: mean diameter in mm = 4/(4.1*R^-0.21), with R in mm/hr, assuming
  !Marshall-Palmer distribution. -JN
  radius = 2._r8/(4.1_r8*(totrnfx(i,k)*60._r8*60._r8)**-0.21_r8)

  !convert radius to meters:
  radius = radius/1000._r8

  !Calculate fraction equilibrated:
  call wtrc_equil_time(ispec,state%t(i,k),state%pmid(i,k),radius,dz(i,k),alpliq,difrm(ispec),fequil) 

 !Constrain equilibration fraction:
  fequil = min(1._r8,max(0._r8,fequil))

  !calculate evaporation rate:
  !if((heff(i,k) < 1.0_r8) .and. (fr .lt. 1._r8)) then  !Is the air un-saturated and rain is evaporating?
  if(.false.) then !Only equilibrate when using wtrc_equil_time. -JN

    !calculate parameters
    gam = (alpliq*heff(i,k))/(1._r8-alpliq*alpkin*(1._r8-heff(i,k)))
    bet = (1._r8-alpliq*alpkin*(1._r8-heff(i,k)))/(alpliq*alpkin*(1._r8-heff(i,k)))

    !Calculate new rain ratio applying Stewart, 1975 equation:
    Rstw = gam*Rv(i,k) + (Rr(i,k) - gam*Rv(i,k)) * fr**bet

    !apply equilibration fraction:
     Re = fequil*Rstw + (1._r8-fequil)*Rr(i,k)

    !Calculate rain re-evaporation:
    if((fr .eq. 0._r8) .and. (bet .lt. 0._r8)) then !Prevent NaNs (plus there shouldn't be any evaporation anyways)
      evp(i,k,m) = rnbulk(i,k)*Rr(i,k)
    else 
     ! evp(i,k,m) = ((flxpr(i,k,m)-flxsn(i,k,m)) - Rstw*((flxpr(i,k,1)-flxsn(i,k,1))-(evp(i,k,1)*state%pdel(i,k)/gravit)))*&
     !              gravit/state%pdel(i,k)
      evp(i,k,m) = ((flxpr(i,k,m)-flxsn(i,k,m)) - Re*((flxpr(i,k,1)-flxsn(i,k,1))-(evp(i,k,1)*state%pdel(i,k)/gravit)))*&
                   gravit/state%pdel(i,k)        
    end if
  else
    !NOTE:  This evp below should always be zero.  If not, then at least it will match the bulk water (even if its physically incorrect) - JN
    evp(i,k,m) = rnbulk(i,k)*Rr(i,k)

    !calculate liquid-vapor equilibrium:
    ivtmp = state%q(i,k,wtrc_iatype(m,iwtvap))+evp(i,k,m)*deltat     !set temporary storage variables
    iltmp = (((flxpr(i,k,m)-flxsn(i,k,m))*gravit/state%pdel(i,k))-evp(i,k,m))*deltat
    vtmp  = state%q(i,k,wtrc_iatype(1,iwtvap))+evp(i,k,1)*deltat
    ltmp  = (((flxpr(i,k,1)-flxsn(i,k,1))*gravit/state%pdel(i,k))-evp(i,k,1))*deltat
    if((iltmp .gt. 0._r8) .and. (ltmp .gt. 0._r8)) then
      call wtrc_liqvap_equil(alpliq, fequil, vtmp, ltmp, ivtmp, iltmp, dliqiso)
    else
      dliqiso = 0._r8
    end if
    evp(i,k,m) = evp(i,k,m) - dliqiso/deltat                    !Add tendency to evaporation

  end if !RH
else
  evp(i,k,m) = rnbulk(i,k)*Rr(i,k)
end if

!*******************
!Calculate snow melt
!*******************

if(state%t(i,k) > tmelt) then ! melt all snow
  mlt(i,k,m) = flxsn(i,k,m) * gravit/state%pdel(i,k) !also convert back to kg/kg/s
else                           ! retain all snow
  mlt(i,k,m) = 0._r8
end if

!***********************
!Calculate rain freezing
!***********************

if(flxpr(i,k,1) > 0._r8) then
  frz(i,k,m) = min(max(0._r8,flxsn(i,k,1)/flxpr(i,k,1)),1._r8) ! guess snow fraction
else
  frz(i,k,m) = 0._r8
end if
frz(i,k,m) = max(fsnow(i,k),frz(i,k,m))  !  guess again 

if(mlt(i,k,1) > 0._r8) frz(i,k,m) = 0._r8

!**************************
!Calculate snow sublimation
!**************************

!No fractionation occurs during sublimation, so water tracers and water isotopes are equivalent. - JN

sub(i,k,m) = subbulk(i,k)*Rs(i,k)

!********************
!Calculate tendencies
!********************

!NOTE:  dqdt seems to work if evpbulk is used - JN

dqdt(i,k,wtrc_iatype(m,iwtvap)) = evp(i,k,m)+sub(i,k,m)       !dqdt must remain pcnst
dldt(i,k,m) = -(evp(i,k,m)+sub(i,k,m))                        !mlt(i,k,m)-frz(i,k,m)
dsdt(i,k,m) = frz(i,k,m)*rprd(i,k,wtrc_iatype(m,iwtvap))-sub(i,k,m)-mlt(i,k,m)

!********************
!Modify precipitation
!********************

!NOTE:  The k+1 allows precipitation to be added to the next level down, which is why the variables
!can be called in the freezing and melting calculations even though the variables are first created
!here. -JN

!pdel/gravit converts kg/kg/s to kg/m2/s 

flxpr(i,k+1,m) = flxpr(i,k,m)+(rprd(i,k,wtrc_iatype(m,iwtvap))+dldt(i,k,m))*state%pdel(i,k)/gravit !modify rain
flxsn(i,k+1,m) = flxsn(i,k,m)+dsdt(i,k,m)*state%pdel(i,k)/gravit                                   !modify snow

!prevent negative numbers:
flxpr(i,k+1,m) = max(flxpr(i,k+1,m), 0._r8)
flxsn(i,k+1,m) = max(flxsn(i,k+1,m), 0._r8)

!Compute total precipitation without any rain re-evaporation
if(wisotope .and. m .eq. 1) then !Are isotope on, and is it H2O?
  totrnfx(i,k+1) = totrnfx(i,k) + ((1._r8-frz(i,k,m))*rprd(i,k,wtrc_iatype(m,iwtvap))+&
                                   mlt(i,k,m)-evp(i,k,m))*state%pdel(i,k)/gravit
  
  !prevent negative numbers:
  totrnfx(i,k+1) = max(totrnfx(i,k+1), 0._r8)
end if

!*********
!End loops
!*********

end do !water tracers

    end do !pver
  end do   !pcols

!Dividing by 1000 converts kg/m2/s to m/s for H2O (aka density of liquid water is 1000 kg/m3) -JN

do m=1,wtrc_nwset

  !Set surface Rain (Precip.) values:
  prec(:ncol,wtrc_iatype(m,iwtvap)) = flxpr(:ncol,pverp,m)/1000_r8

  !Set surface Snow values:
  snow(:ncol,wtrc_iatype(m,iwtvap)) = flxsn(:ncol,pverp,m)/1000_r8

end do     !water tracers

!-----------------------
!correct for mass errors:
!-----------------------
do i=1,ncol
 !Calculate differences:
  rmass0 = (prec(i,wtrc_iatype(2,iwtvap))-snow(i,wtrc_iatype(2,iwtvap)))
  smass0 = snow(i,wtrc_iatype(2,iwtvap))
  rdiff  = rmass0 - (prec(i,wtrc_iatype(1,iwtvap))-snow(i,wtrc_iatype(1,iwtvap)))
  sdiff  = smass0-snow(i,wtrc_iatype(1,iwtvap))
  do m=2,wtrc_nwset
   !Rain errors:
    Rd = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),(prec(i,wtrc_iatype(m,iwtvap))-snow(i,wtrc_iatype(m,iwtvap))),rmass0)
    prec(i,wtrc_iatype(m,iwtvap)) = prec(i,wtrc_iatype(m,iwtvap))-Rd*rdiff
   !Snow errors:
    Rd = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),snow(i,wtrc_iatype(m,iwtvap)),smass0)
    snow(i,wtrc_iatype(m,iwtvap)) = snow(i,wtrc_iatype(m,iwtvap))-Rd*sdiff
   !Giant error check: 
   !NOTE:  !Seems to occur about once every seven years. -JN
   ! if(prec(i,wtrc_iatype(m,iwtvap)) .gt. 10._r8*prec(i,wtrc_iatype(1,iwtvap))) then
    if(prec(i,wtrc_iatype(m,iwtvap)) .gt. 2._r8*prec(i,wtrc_iatype(1,iwtvap))) then
      if(prec(i,wtrc_iatype(1,iwtvap)) .gt. wtrc_qmin) &
      write(*,*) 'ERROR:  Isotopic deep-conv precip error!',prec(i,wtrc_iatype(m,iwtvap)),prec(i,wtrc_iatype(1,iwtvap)),&
                  snow(i,wtrc_iatype(m,iwtvap)),snow(i,wtrc_iatype(1,iwtvap)),i,m
     !Set the tracer precip back to bulk water (it violates mass-conservation,
     !but is probably due to a numerical/non-physical error anyways).
      prec(i,wtrc_iatype(m,iwtvap)) = prec(i,wtrc_iatype(1,iwtvap))
      snow(i,wtrc_iatype(m,iwtvap)) = snow(i,wtrc_iatype(1,iwtvap))
    end if
  end do
end do
!----------------------

!**************
!End subroutine
!**************
end subroutine wtrc_precip_evap


!=======================================================================
subroutine wtrc_equil_time(ispec,temp,pres,rdrop,zdel,alpha,difrm,fequil)
!-----------------------------------------------------------------------
!
! Function to compute adjustment timescale for isotopic equilibration.
! This routine could be bolted into CAM.
!
! Final adjustment fo isotope ratio of rain drops woiuld be:
!
!
!  R_rain = (1.-fequil)*R_original * fequil*R_stewart
!
! This code outputs fequil
!
!
! David Noone <dcn@coas.oregonstate.edu> - Thu Jun  4 13:31:25 MDT 2015
!
! Modified for use in CAM5 by Jesse Nusbaumer <nusbaume@colorado.edu> - June, 2015
!
!-----------------------------------------------------------------------

use physconst,     only: rair, rh2o, gravit, rhoh2o
use wv_saturation, only: qsat

!-----------------------------------------------------------------------
!Declare variables
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
!------------------------ Output Arguments -----------------------------

  real(r8) , intent(out) :: fequil ! adustment fraction

!------------------------- Input Arguments -----------------------------

  integer    , intent(in) :: ispec  ! isotopologue species flag
  real(r8)   , intent(in) :: temp   ! temperature (K)
  real(r8)   , intent(in) :: pres   ! pressure (Pa)
  real(r8)   , intent(in) :: rdrop  ! drop RADIUS (m)
  real(r8)   , intent(in) :: zdel   ! layer thickness (m)
  real(r8)   , intent(in) :: alpha  ! equilibrium fractionation factor
  real(r8)   , intent(in) :: difrm  ! ratio of isotope diffusivities 

!-----------------------------------------------------------------------

  real(r8) :: tadjust                  ! adjustment time (sec)
  real(r8) :: texposure                ! time in layer (sec)

  real(r8) :: rhoa                     ! density of air (kg/m^3)

  real(r8) :: Re                       ! Reynolds nunber
  real(r8) :: Sc                       ! Schemidt number
  real(r8) :: xx                       ! Re**1/2 * Sc** 1/3

  real(r8) :: difa                     ! diffusivity of H2O in air (m^2/s)
  real(r8) :: mu                       ! dynamic viscosity (kg/m/s)
  real(r8) :: vterm                    ! terminal velocity (m/s)
  real(r8) :: cddrop                   ! droplet drag coefficient
  real(r8) :: fvent                    ! ventillation factor (unitless)

 !saturation vapor pressure variables:
  real(r8) :: esat                     ! saturation vapor pressure (Pa)
  real(r8) :: qst                      ! saturation mixing ratio (kg/kg)

!-----------------------------------------------------------------------

!
! Compute density of air on input (equation of state):
!
  rhoa = pres/(Rair*temp)       ! density of air (kg/m3)

!
! Calculate saturation vapor pressure:
!
  call qsat(temp, pres, esat, qst)

!
! Diffusivity of H2O in air (Hall and Pruppacher 1976)
!
  difa = 2.11e-5_r8*difrm*(temp/273.15_r8)**1.94_r8 * (101325._r8/pres)      ! m2/sec

! No fractionation:
!  difa = 2.11e-5_r8*(temp/273.15_r8)**1.94_r8 * (101325._r8/pres) ! m2/sec

!
! Viscosity of air (Sutherland's formula)
!
!NOTE:  use equation from MOZART aerosol subroutines (Rodgers and Yau, pg. 103):
  mu = 1.72e-5_r8 * ((temp/273.0_r8)**1.5_r8) * 393.0_r8 / (temp+120.0_r8)   

!
! Drop terminal fall velocity (from Foote and Du Toit, 1969):
!
  cddrop = 0.6_r8 !From Straka, 2009
  vterm = sqrt((4._r8/3._r8)*gravit*(2._r8*rdrop)*rhoh2o/(cddrop*rhoa) )

!
! Reynolds number (note factor of 2, since length scale should be diamater):
!
  Re = 2._r8*rdrop*rhoa*vterm/mu !equation 10-10 in P&K, 1997

!
! Schmidt number (Pruppacher and Klett):
!
 Sc = mu/(rhoa*difa)

!
! Ventillation factor (Pruppacher and Klett, 1997):
!
 xx = Re**(1._r8/2._r8) * Sc**(1._r8/3._r8)

 if(xx .ge. 1.4_r8) then                      !cut off from Pruppacher and Rasmussen, 1979
   fvent = 0.78_r8 + 0.308_r8*xx              !equation 13-61 in P&K, 1997
 else
   fvent = 1.0_r8 + 0.108_r8*xx*xx            !equation 13-60 in P&K, 1997
 endif

!
! Compute adjustment time (Stewart, Lee and Fung)
! (Note, diffusivity here SHOULD be for isotopes IF you were interested
!  in doing the kinetic component explictly along side the total
!  fractionation)
!
!  Equation from Lee and Fung, 2008:

  tadjust = alpha*(rdrop*rdrop)*rhoh2o*rh2o*temp/(3._r8*fvent*difa*esat)

!No fractionation:
!  tadjust = (rdrop*rdrop)*rhoh2o*rh2o*temp/(3._r8*fvent*difa*esat)

!
! Compute the exposure time as the time falling throught this layer
!
  texposure = zdel/vterm

!
! Finally the adjustment fraction
!
  fequil = 1._r8 - exp(-texposure/tadjust)

!
  return
end subroutine wtrc_equil_time


!=======================================================================
function wtrc_efac(alpha, vapnew, liqnew)
!-----------------------------------------------------------------------
! Computes isotopic equilibration factor - there are different ways to
! compute this which have slightly different numerical properties.
! Author: David Noone <dcn@colorado.edu> - Fri Aug 13 11:26:19 MDT 2004
!-----------------------------------------------------------------------
  real(r8)  , intent(in) :: alpha       ! fractionation factor
  real(r8)  , intent(in) :: vapnew      ! new total vapour
  real(r8)  , intent(in) :: liqnew      ! new total liquid
  real(r8) wtrc_efac                    ! return value
!-----------------------------------------------------------------------
  real(r8) efac                         ! equilibration factor
  real(r8) alov                         ! alpha times l on v
  real(r8) qtot                         ! total water
!-----------------------------------------------------------------------

!#define DIRECTWAY
!#ifdef DIRECTWAY         /* This, most obvious way is less precise */
!  alov = alpha*wtrc_ratio(isph2o, liqnew, vapnew)
!#else                   /* this one can have more precise numerics */
  alov = wtrc_ratio(isph2o, vapnew, vapnew+liqnew)

  alov = alpha*(1.0_r8/alov - 1.0_r8)   ! this is l/v
!#endif 
 
  efac = 1.0_r8/(alov + 1.0_r8)

! 
  efac = max(efac, 0._r8)  !Ensure that equilibration factor is between zero and one
  efac = min(efac, 1._r8)
! 
  wtrc_efac = efac
  return
end function wtrc_efac


!=======================================================================
!Version 3
!=======================================================================
!=======================================================================
subroutine wtrc_q1q2_pjr(dqdt        , ideep, lengath, &
                    q       ,qub     , done   ,lel     ,lcl     , &
                    dupc    ,du      ,dp      ,qstb    ,qsthatb , &
                    hmn     ,hsatb   ,hsthatb ,eu      , &
                    mu      ,md      ,qdb     ,tu      ,td      , &
                    dsubcld ,mx      ,jt      ,il2g    , &
                    msg     ,ed      ,rppe    ,mupc    ,mdpc    , &
                    evpc    ,evp     ,cupc    ,cu      ,pap     , &
                    dz      ,qds     ,rpdpc   ,c0mask  ,eps0    , &
                    jd      ,wtrprd  ,wtdlf   ,dtime)


!----------------------------------------------------------------------- 
! 
! Purpose: This subroutine is designed to 
!calculate the water tracer/isotope tendencies due to
!the ZM deep convective scheme.
! 
! Author: Jesse Nusbaumer <nusbaume@colorado.edu> - August, 2012
!
!
! NOTE: The temperatures used in this code are temperatures that exist
!       after (before?) condensation/evaporation of the model (bulk) water has
!       occured.  Ideally, the temperature before and after the latent
!       heating or cooling would be known, and the phase change tendencies
!       would be solved iterively over a changing temperature, similar
!       to what's done in the cloud physics.  For now,  using the final 
!       temperature should be fine. - JN
! 
!-----------------------------------------------------------------------

!**************                                                                           
!Use statements
!**************

 use physconst,       only: cpair, latvap, rair, epsilo
 use ppgrid,          only: pverp

!********************
!Variable decleration
!********************

 implicit none
!
! Input fields:
!
   logical, intent(in) :: done(pcols,pver) !updraft loop endpoint 
   integer, intent(in) :: lel(pcols)       !Equilibrium level (EL)
   integer, intent(in) :: lcl(pcols)       !lifted condensation level (LCL)    
   integer, intent(in) :: mx(pcols)        !Bottom of updraft
   integer, intent(in) :: jd(pcols)        !Top of downdraft
   integer, intent(in) :: jt(pcols)        !Top of updraft 
   integer, intent(in) :: il2g
   integer, intent(in) :: msg
   integer, intent(in) :: ideep(pcols)         !location of deep convection
   integer, intent(in) :: lengath              !number of grid spaces with deep convection

   real(r8), intent(in) :: q(pcols,pver,pcnst) !Environmental water vapor mixing ratio (state%q) [kg/kg]
   real(r8), intent(in) :: qstb(pcols,pver)    !bulk saturated env. water vapor mixing ratio [kg/kg]
   real(r8), intent(in) :: qsthatb(pcols,pver) !qst at interfaces
   real(r8), intent(in) :: qub(pcols,pver)     !bulk water vapor mixing ratio in updraft [kg/kg]
   real(r8), intent(in) :: hmn(pcols,pver)     !bulk environmental moist static energy
   real(r8), intent(in) :: hsatb(pcols,pver)   !bulk saturated env. moist static energy
   real(r8), intent(in) :: hsthatb(pcols,pver) !hsat at interfaces
   real(r8), intent(in) :: dp(pcols,pver)      !layer thickness [mb]
   real(r8), intent(in) :: mdpc(pcols,pver)    !md before unit-change
   real(r8), intent(in) :: dupc(pcols,pver)    !detrainment rate in updraft pre-unit change
   real(r8), intent(in) :: du(pcols,pver)      !input detrainment rate in updraft 
   real(r8), intent(in) :: mupc(pcols,pver)    !mu before unit-change
   real(r8), intent(in) :: eu(pcols,pver)      !input updraft entraiment rate
   real(r8), intent(in) :: mu(pcols,pver)      !mass flux in updraft 
   real(r8), intent(in) :: md(pcols,pver)      !mass flux in downdraft
   real(r8), intent(in) :: qdb(pcols,pver)     !bulk water vapor mixing ratio in downdraft [kg/kg]
   real(r8), intent(in) :: tu(pcols,pver)      !temperature in updraft [K?]
   real(r8), intent(in) :: td(pcols,pver)      !temperature in downdraft [K?]
   real(r8), intent(in) :: ed(pcols,pver)     !downdraft entrainment rate [?]
   real(r8), intent(in) :: evpc(pcols,pver)    !evaporation rate pre-unit change [kg/m3/s]
   real(r8), intent(in) :: evp(pcols,pver)     !evaporation rate in downdraft
   real(r8), intent(in) :: cupc(pcols,pver)    !condensation rate pre-unit change [kg/m3/s]
   real(r8), intent(in) :: cu(pcols,pver)      !condensation rate in updraft
   real(r8), intent(in) :: rppe(pcols,pver)    !g: rain production pre-evaporation [?]
   real(r8), intent(in) :: dsubcld(pcols)      !sub-cloud layer thickness [mb]
   real(r8), intent(in) :: pap(pcols,pver)     !Pressure at interfaces [Pa]
   real(r8), intent(in) :: dz(pcols,pver)      !Layer thickness in height [m] 
   real(r8), intent(in) :: qds(pcols,pver)     !bulk water saturation value in downdraft [kg/kg]
   real(r8), intent(in) :: rpdpc(pcols,pver)   !precipitation production pre-unit change [?]
   real(r8), intent(in) :: c0mask(pcols)       !Autoconversion rates [1/m]
   real(r8), intent(in) :: eps0(pcols)         !Not quite sure, should look up... 
   real(r8), intent(in) :: dtime               !2xdt (model timestep)
!
! Output fields:
!
   real(r8),intent(out) :: dqdt(pcols,pver,pcnst)     !water tracer tendency [kg/kg/s]
   real(r8),intent(out) :: wtrprd(pcols,pver,pcnst)     !precipitation production rate [kg/kg/s]
   real(r8),intent(out) :: wtdlf(pcols,pver,wtrc_nwset) !tendency due to condensate detrainment [kg/kg/s] 
!
! Work fields:
!
   integer i
   integer k
   integer m
   integer kbm
   integer ktm
   integer kount
   integer blm
   integer tlm

  !NOTE:  These variables below may not need to be 3-D arrays, and could simply be reset
  !after every tracer loop.  However, for now it should be left as is (because it works) - JN

   real(r8) qhat(pcols,pver,wtrc_ntype(iwtvap))  !water vapor mixing ratios at level interfaces [kg/kg]
   real(r8) qu(pcols,pver,wtrc_ntype(iwtvap))    !water vapor mixing ratio in updraft [kg/kg]
   real(r8) qd(pcols,pver,wtrc_ntype(iwtvap))    !water vapor mixing ratio in downdraft [kg/kg]
   real(r8) wtcu(pcols,pver,wtrc_ntype(iwtvap))  !water tracer condensation rate [?]
   real(r8) wtevp(pcols,pver,wtrc_ntype(iwtvap)) !water tracer evaporation rate [?]
   real(r8) hu(pcols,pver,wtrc_ntype(iwtvap))    !tracer moist static energy in updraft
   real(r8) hd(pcols,pver,wtrc_ntype(iwtvap))    !tracer moist static energy in downdraft
   real(r8) wthmn(pcols,pver,wtrc_ntype(iwtvap)) !tracer environmental moist static energy
   real(r8) hsat(pcols,pver,wtrc_ntype(iwtvap))  !tracer saturated env. moist static energy
   real(r8) qst(pcols,pver,wtrc_ntype(iwtvap))   !tracer saturated env. vapor mixing ratio
   real(r8) gamma(pcols,pver,wtrc_ntype(iwtvap)) !No clue, just need it for parameterization
   real(r8) hsthat(pcols,pver,wtrc_ntype(iwtvap))!hsat at interfaces 
   real(r8) qsthat(pcols,pver,wtrc_ntype(iwtvap))!qsat at interfaces 
   real(r8) gamhat(pcols,pver,wtrc_ntype(iwtvap))!Gamma at interfaces 
   real(r8) Ru(pcols,pver,wtrc_ntype(iwtvap))    !water tracer updraft ratio [unitless]
   real(r8) Rd(pcols,pver,wtrc_ntype(iwtvap))    !water tracer downdraft ratio [unitless]
   real(r8) huct(pcols,pver)                     !MSE storage variable for cloud top height

   !For isotopic fractionation:
   real(r8) alpha                                !fractionation factor [unitless]
!   real(r8) alpcn                                !condensation fractionation factor [unitless] 
!   real(r8) alpev                                !evaporation fractionation factor [unitless]
   real(r8) dliqiso                              !needed for isotopic equilbration subroutine [kg/kg]
   real(r8) nu                                   !ice/liquid fraction [unitless]
   real(r8) que                                  !equilibrated updraft vapor [kg/kg] 
   real(r8) qle                                  !equilibrated cloud condensate [kg/kg]   
   real(r8) qur                                  !Rayleigh-distilled updraft vapor [kg/kg]
   real(r8) qlr                                  !Rayleigh-distilled updraft vapor [kg/kg] 
   real(r8) cue                                  !equilibrated updraft condensate [kg/kg]
   real(r8) cur                                  !Rayleigh-distilled updraft condensate [kg/kg]
   real(r8) cuh                                  !H2O updraft condensate [kg/kg]
   real(r8) fr                                   !Fraction of moisture remaining [unitless] 
   real(r8) ltmp                                 !temporary standard liquid amount [kg/kg]
   real(r8) iltmp                                !temporary isotopic liquid amount [kg/kg]

   !For precip production:
   real(r8) totpcp(pcols,wtrc_ntype(iwtvap))     !Total precipitation [?] 
   real(r8) totevp(pcols,wtrc_ntype(iwtvap))     !Total evaporation [?] 
   real(r8) rprd(pcols,pver)                     !Precipitation production pre-unit change [?]
   !NOTE:  Might need to make "pevp" the evaporation term used in the actual
   !tendency calculation! - JN
   real(r8) pevp(pcols,pver)                     !precipitation evaporation 
   real(r8) ql(pcols,pver,wtrc_nwset)            !precipitable liquid [kg/kg]
   real(r8) ql1                                  !temporary precipitable liquid storage [kg/kg]
   real(r8) Rr                                   !precipitation production ratio [unitless] 

   !For mass fixer:
   real(r8) uqdiff                               !correction factor for updraft humidity [kg/kg]
   real(r8) dqdiff                               !correction factor for downdraft humidity [kg/kg]  
   real(r8) oval                                 !original H2O water mass before mass-fixing [kg/kg]
   real(r8) Rfix                                 !Water tracer/isotope ratio used to apply mass fix [unitless]

   !Other variables:
   real(r8) Rc                                   !water tracer condensation ratio [unitless] 
   real(r8) Re                                   !water tracer evaporation ratio [unitless]
   real(r8) qdifr
   real(r8) emc
   real(r8) mdt

!***************************************
!Exit if no deep convection is occurring
!***************************************

if(lengath.eq.0) then

   !Set tendencies to zero:
   do m=1,wtrc_nwset
     dqdt(:,:,wtrc_iatype(m,:)) = 0.0_r8
   end do

   return

else

!***********************
!Loop over water tracers
!***********************

  do m=1,wtrc_nwset 

    !Initalize variables:
    dqdt(:,:,wtrc_iatype(m,:))        = 0.0_r8
    wtrprd(:,:,wtrc_iatype(m,iwtvap)) = 0.0_r8
    wtdlf(:,:,m)         = 0.0_r8 
    qhat(:,:,m)          = q(:,:,wtrc_iatype(m,iwtvap)) !Done to match zm_convr
    qu(:,:,m)            = 0.0_r8                       !<-possibly not needed...
    qd(:,:,m)            = 0.0_r8
    qu(1:lengath,:,m)    = q(ideep(:),:,wtrc_iatype(m,iwtvap))
    qd(1:lengath,:,m)    = q(ideep(:),:,wtrc_iatype(m,iwtvap))
    ql(:,:,m)            = 0.0_r8
    wtcu(:,:,m)          = 0.0_r8
    wtevp(:,:,m)         = 0.0_r8
    wthmn(:,:,m)         = 0.0_r8
    hsat(:,:,m)          = 0.0_r8
    qst(:,:,m)           = 0.0_r8
    gamma(:,:,m)         = 0.0_r8
    hsthat(:,:,m)        = 0.0_r8
    qsthat(:,:,m)        = 0.0_r8
    gamhat(:,:,m)        = 0.0_r8
    hu(:,:,m)            = 0.0_r8
    hd(:,:,m)            = 0.0_r8
    totpcp(:,m)          = 0.0_r8
    totevp(:,m)          = 0.0_r8
    Ru(:,:,m)            = 0.0_r8
    Rd(:,:,m)            = 0.0_r8
    uqdiff               = 0.0_r8
    dqdiff               = 0.0_r8

  end do

!**************************
!Calculate interface values (qhat)
!**************************

!NOTE:  The natural logarithm (log) is NOT a linear operator,
!thus the ratio will not neccesarily be carried through, and
!could result in a loss of tracer/tracer ratio conservation. - JN

!-------------
!water  vapor:
!-------------

  do m=1,wtrc_nwset !water tracers
    do k = msg + 2,pver
      do i = 1,lengath
         qdifr = 0.0_r8
         if (q(ideep(i),k,wtrc_iatype(m,iwtvap)) > 0._r8 .or. q(ideep(i),k-1,wtrc_iatype(m,iwtvap)) > 0._r8) &
           qdifr = abs((q(ideep(i),k,wtrc_iatype(m,iwtvap))- &
                   q(ideep(i),k-1,wtrc_iatype(m,iwtvap)))/max(q(ideep(i),k-1,wtrc_iatype(m,iwtvap)),&
                   q(ideep(i),k,wtrc_iatype(m,iwtvap))))
         if ((qdifr > 1.E-6_r8) .and. (qdifr /= 1._r8)) then !qdifr /= 1 used to prevent divison by zero. - JN
           qhat(i,k,m) = log(q(ideep(i),k-1,wtrc_iatype(m,iwtvap))/q(ideep(i),k,wtrc_iatype(m,iwtvap))) * &
                                                q(ideep(i),k-1,wtrc_iatype(m,iwtvap)) * &
                                                q(ideep(i),k,wtrc_iatype(m,iwtvap))/(q(ideep(i),k-1,wtrc_iatype(m,iwtvap)) - &
                                                q(ideep(i),k,wtrc_iatype(m,iwtvap)))
         else
           qhat(i,k,m) = 0.5_r8* (q(ideep(i),k,wtrc_iatype(m,iwtvap))+q(ideep(i),k-1,wtrc_iatype(m,iwtvap)))
         end if
      end do
    end do
  end do

!***********************************************************
!calculate updraft humidity with cloud condensate production
!***********************************************************

!NOTE:  Although according to the CAM5 documentation, the equation for qu should
!be complete, there still appear to be small errors [O(10e-16)] that occur between the equation
!and the standard bulk updraft water.  Thus a "qdiff" variable is included to compensate for
!these differences.  As this qdiff appears to be directly proportional to the original tracer ratio,
!it should only need to be calculated once. - JN

  do k = pver,msg + 2,-1
    do i = 1,lengath
      if ((k > 1 .and. k < mx(i)) .and. eps0(i) > 0._r8) then !Not sure why this is needed, but it errors without it - JN
        if(mupc(i,k) > 0._r8) then  !Do only if mass flux is positive 

          !++++++++++++++++++++++++++++++++++++++++++
          !determine liquid/ice fraction for isotopes
          !++++++++++++++++++++++++++++++++++++++++++
          !NOTE:  Using same temperature cut-off as Macrophysics for
          !convective detrainment. - JN 
          if(tu(ideep(i),k) > 268.15_r8) then
            nu = 0.0_r8
          elseif(tu(ideep(i),k) < 238.15_r8) then
            nu = 1.0_r8
          else
            nu = (268.15_r8 - tu(ideep(i),k)) / 30._r8
          endif
          !+++++++++++++++++++++++++++++++++++++++++++++

          do m=1,wtrc_nwset      

            !Calculate ratio:
            Ru(i,k,m) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),qu(i,k+1,m),qu(i,k+1,1)) !calculate ratio

           !Calculate condensation:
            wtcu(i,k,m) = Ru(i,k,m)*cupc(i,k)

            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !Calculate isotopic fractionation in the updraft condensation
            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(wisotope .and. (m .gt. 1)) then !only do for isotopes:
!            if(.false.) then
              !Equilibrium fractionation factor:
              alpha = wtrc_get_alpha(qub(i,k+1),tu(ideep(i),k+1),iwspec(wtrc_iatype(m,iwtvap)),&
                                iwtvap,iwtliq,.false.,1._r8,.false.)

             !---------------------------
             !Equilibrate liquid fraction
             !---------------------------
              cue = wtcu(i,k,m)*dz(i,k)/mupc(i,k) !convert to kg/kg
              cuh = wtcu(i,k,1)*dz(i,k)/mupc(i,k) !H2O condensation
              !Equilibrate:
              call wtrc_liqvap_equil(alpha, 1._r8, qu(i,k+1,1), cuh,&
                                     qu(i,k+1,m), cue, dliqiso)

             !determine value of original condensation amount if equilibrated:
              cue = wtcu(i,k,m) + (dliqiso*mupc(i,k)/dz(i,k))
             !multiply by liquid fraction:
              cue = cue * (1._r8 - nu) 
             !---------------------------
             !Distill ice fraction
             !---------------------------
             !calculate new fractionation factor:
              alpha = wtrc_get_alpha(qub(i,k+1),tu(ideep(i),k+1),iwspec(wtrc_iatype(m,iwtvap)),&
                              iwtvap,iwtice,.true.,pap(i,k),.true.)
             !fraction of mass that's condensed:
              fr = (qu(i,k+1,1)-cuh)/qu(i,k+1,1)
             !limit fraction:
              fr = min(1._r8,max(0._r8,fr))
             !distill isotope mass:
              qur = qu(i,k+1,m)*fr**alpha
             !determine amount of mass condensed:
              cur = qu(i,k+1,m)-qur
             !determine value of original condensation amount if distilled:
              cur = cur*mupc(i,k)/dz(i,k)
             !multiply by ice fraction:
              cur = cur * nu
             !--------------------------------------
             !combine into final condensation amount
             !--------------------------------------
              wtcu(i,k,m) = cue + cur
             !--------------------------------------
            end if !wisotope
            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            !Calculate updraft humidity:
            qu(i,k,m) = mupc(i,k+1)/mupc(i,k)*qu(i,k+1,m) + dz(i,k)/mupc(i,k)* (eu(i,k)*q(ideep(i),k,wtrc_iatype(m,iwtvap))- &
                              dupc(i,k)*Ru(i,k,m)*qstb(i,k) - wtcu(i,k,m))

            !Apply mass fix:
            !--------------
            if(m .eq. 1) then !only do for H2O tracer 
              oval   = qu(i,k,m) !save original value
              uqdiff = qub(i,k)-qu(i,k,m)
            end if
            Rfix      = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),qu(i,k,m),oval)
            qu(i,k,m) = qu(i,k,m)+Rfix*uqdiff !add back difference
            !-------------- 

            !Produce cloud condensate:
            if (k >= jt(i)) then
              ql1 = 1._r8/mupc(i,k)* (mupc(i,k+1)*ql(i,k+1,m)- &
              dz(i,k)*dupc(i,k)*ql(i,k+1,m)+dz(i,k)*wtcu(i,k,m))

              ql(i,k,m) = ql1/ (1._r8+dz(i,k)*c0mask(i))
            else
              ql(i,k,m) = 0._r8
            end if

          end do !water tracers
     
        else 

          do m=1,wtrc_nwset
            qu(i,k,m) = q(ideep(i),k,wtrc_iatype(m,iwtvap)) !If not, then assume updraft vapor is env. vapor
          end do

        end if
      end if
      !complain if uqdiff is too large: 
      !----------------
      if((uqdiff/qub(i,k)) .gt. 10e-10_r8) then !relative error greater than 1e-9?
        write(iulog,*) 'wtrc q1q2 uqdiff error here!',uqdiff/qub(i,k),uqdiff,&
                        qub(i,k),i,k
      end if
      !---------------
    end do
  end do

!********************************************
!Calculate downdraft humidity and evaporation
!********************************************

!NOTE:  These two quantities must be computed at the same time. - JN

!NOTE:  Assume initial downdraft vapor is equal to the ratio of the updraft vapor
!at the starting height (assumes that the vapor from the drowndraft comes directly
!from the updraft). - JN

!NOTE:  Also perform unphysical "qdiff" process here, just to make sure
!ratios are conserved. - JN.

!NOTE:  The evaporation term appears to be precipitation re-evaporation.  For
!now, it might be wise to assume that its all rain or all snow, until a suitable
!temperature cutoff can be determined.  NEW NOTE:  Given that rain re-evaporation
!occurs again in the wtrc_precip_evap subroutine, and given that the downdraft should
!be close to saturated (at least if I am remembering the ZM scheme correctly), and
!given that rain is falling in a column of air (the downdraft) that is also
!moving downward, and thus the rain might be in contact with the same air for
!longer than normal, it is assumed that the precipitation is always rain, and
!that it equilibrates isotopically with the downdraft air. - JN 

   do i = 1,lengath
     do m=1,wtrc_nwset
       Rd(i,jd(i),m) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),qu(i,jd(i),m),qu(i,jd(i),1))
       qd(i,jd(i),m) = Rd(i,jd(i),m)*qds(i,jd(i))
     end do
   end do

   do k = msg+2,pver
     do i=1,lengath
       if (k >= jd(i) .and. k < mx(i) .and. eps0(i) > 0._r8) then

         do m=1,wtrc_nwset

           Rd(i,k,m) = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),qd(i,k,m),qd(i,k,1)) !calculate tracer ratio

           wtevp(i,k,m) = Rd(i,k,m)*evpc(i,k)

           qd(i,k+1,m) = (dz(i,k)*(wtevp(i,k,m)+ed(i,k)*q(ideep(i),k,wtrc_iatype(m,iwtvap))) - mdpc(i,k)*qd(i,k,m))/(-mdpc(i,k+1))

           !Apply mass fix:
           !---------------
           if(m .eq. 1) then !only do for H2O tracer 
             oval   = qd(i,k+1,m) !save original value
             dqdiff = qdb(i,k+1)-qd(i,k+1,m)
           end if
           Rfix        = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),qd(i,k+1,m),oval)
           qd(i,k+1,m) = qd(i,k+1,m)+Rfix*dqdiff
           !---------------

           totevp(i,m) = totevp(i,m) - dz(i,k)*ed(i,k)*q(ideep(i),k,wtrc_iatype(m,iwtvap)) 

         end do

       end if
      !complain if dqdiff is too large: 
      !----------------
       if((dqdiff/qdb(i,k)) .gt. 10e-10_r8) then !relative error greater than 1e-9?
         write(iulog,*) 'wtrc q1q2 dqdiff error here!',dqdiff/qdb(i,k),dqdiff,&
                         qdb(i,k),i,k
       end if
      !---------------
     end do
   end do

!***************************
!calculate total evaporation
!***************************

do m=1,wtrc_nwset
 do i= 1,lengath
   totevp(i,m) = totevp(i,m) + mdpc(i,jd(i))*qd(i,jd(i),m) - mdpc(i,mx(i))*qd(i,mx(i),m)
 end do
end do

!*******************************
!Calculate rain production rates (rprd)
!*******************************

 do m=1,wtrc_nwset !wtrc_ntype(iwtvap)

   pevp(:,:) = 0._r8  !inialize precipitable liquid
   rprd(:,:) = 0._r8

!Produce precipitation:
   do k = pver,msg + 2,-1
     do i = 1,lengath
       if (k >= jt(i) .and. k < mx(i) .and. eps0(i) > 0._r8 .and. mupc(i,k) >= 0.0_r8) then
         totpcp(i,m) = totpcp(i,m) + dz(i,k)*(wtcu(i,k,m)-dupc(i,k)*ql(i,k+1,m))
         rprd(i,k) = c0mask(i)*mupc(i,k)*ql(i,k,m) !finalize rain production
       end if
     end do
   end do

!Evaporate precipitation:
    do k = msg + 2,pver
      do i = 1,lengath
        totpcp(i,m) = max(totpcp(i,m),0._r8) !Remove negative values
        totevp(i,m) = max(totevp(i,m),0._r8)
        if (totevp(i,m) > 0._r8 .and. totpcp(i,m) > 0._r8) then
          pevp(i,k) = wtevp(i,k,m)*min(1._r8, totpcp(i,m)/(totevp(i,m)+totpcp(i,m)))
        else
          pevp(i,k) = 0._r8
        end if
        rprd(i,k) = rprd(i,k)-pevp(i,k) !evaporate precipitation
        Rr = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),rprd(i,k),rpdpc(i,k)) !calculate precip ratio
        wtrprd(ideep(i),k,wtrc_iatype(m,iwtvap)) = Rr*wtrprd(ideep(i),k,1)  !calculate precip production (done to avoid unit conversions)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !Isotopically Equilibrate precipitation with downdraft vapor
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !NOTE:  This sections appears to produce regions where the preciptation
        !d-excess is -40 or less, which is way outside the range of observations
        !or theoretical expectations for those locations.  This combined with
        !the fact that the precipitation goes through a second round of possible
        !equilibration in wtrc_precip_evap means that for now it would be better
        !to simply not equilibrate here.  However, one could argue physically
        !that this should still be present, so if later it is found that
        !something generates an abnormally high d-excess over certain land
        !regions, turning this on may be a good way to counterbalance that. - JN
        if(.false.) then 
!        if(wisotope .and. (m .gt. 1)) then !only do for isotopes
          alpha = wtrc_get_alpha(qdb(i,k),td(ideep(i),k),iwspec(wtrc_iatype(m,iwtvap)),&
                                 iwtvap,iwtliq,.false.,1._r8,.false.)

          !create temporary liquid values
          ltmp = wtrprd(ideep(i),k,wtrc_iatype(1,iwtvap))*dtime
          iltmp= wtrprd(ideep(i),k,wtrc_iatype(m,iwtvap))*dtime

          call wtrc_liqvap_equil(alpha, 1._r8, qd(i,k,1), ltmp, qd(i,k,m), iltmp, dliqiso)

          !determine value of original precipitation amount if equilibrated:
          wtrprd(ideep(i),k,wtrc_iatype(m,iwtvap)) = wtrprd(ideep(i),k,wtrc_iatype(m,iwtvap)) + dliqiso/dtime
        end if
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      end do
    end do

 end do !water tracers (m)
!------------

!************************
! find the highest level top and bottom levels of convection
!************************

   ktm = pver
   kbm = pver

   do i=1,lengath
     ktm = min(ktm,jt(i))
     kbm = min(kbm,mx(i))
   end do

!********************************
!Calculate cloud-level tendencies
!********************************

  do m=1,wtrc_nwset
    do k = ktm,pver-1
      do i = 1,lengath
        Rc = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),wtcu(i,k,m),wtcu(i,k,1)) !calculate tracer ratio
        Re = wtrc_ratio(iwspec(wtrc_iatype(m,iwtvap)),wtevp(i,k,m),wtevp(i,k,1))

        !NOTE:  Ratios are used in order to avoid unit conversions. - JN

        emc = -Rc*cu(ideep(i),k)  &    ! condensation in updraft
              +Re*evp(ideep(i),k)      ! evaporating rain in downdraft

         !Calculate tendencies:
         if(dp(i,k) .gt. 0._r8) then !the layer has a finite thickness?
           dqdt(ideep(i),k,wtrc_iatype(m,iwtvap)) = emc + &
                      (+mu(i,k+1)* (qu(i,k+1,m)-qhat(i,k+1,m)) &
                       -mu(i,k)*   (qu(i,k,m)-qhat(i,k,m)) &
                       +md(i,k+1)* (qd(i,k+1,m)-qhat(i,k+1,m)) &
                       -md(i,k)*   (qd(i,k,m)-qhat(i,k,m)) &
                      )/dp(i,k)
         else
           dqdt(ideep(i),k,wtrc_iatype(m,iwtvap)) = emc !no mass fluxes can occur...
         end if

         !calculate detrained condensate:
         wtdlf(ideep(i),k,m) = du(i,k)*ql(i,k+1,m) 

      end do
    end do

!******************************
!Calculate sub-cloud tendencies
!******************************

    do k = kbm,pver
      do i=1,lengath
        if (k == mx(i)) then

           dqdt(ideep(i),k,wtrc_iatype(m,iwtvap)) = (1._r8/dsubcld(i))* &
                                                    (-mu(i,k)*(qu(i,k,m)-qhat(i,k,m)) &
                                                     -md(i,k)*(qd(i,k,m)-qhat(i,k,m)) &
                                                    )
        else if (k > mx(i)) then
           dqdt(ideep(i),k,wtrc_iatype(m,iwtvap)) = dqdt(ideep(i),k-1,wtrc_iatype(m,iwtvap))
        end if
      end do
    end do

  end do  !Water tracers

!******************
!End deep-con check
!******************

end if !Lengath > 0

!*********************

   return
end subroutine wtrc_q1q2_pjr


!=======================================================================
  subroutine  wtrc_setup_diag
!-----------------------------------------------------------------------
! Purpose: Writes configuration of water tracer scheme to standard output.
!-----------------------------------------------------------------------
    use constituents,   only: cnst_name
    use water_isotopes, only: wiso_get_fisub
    use physics_buffer, only: pbuf_get_field_name
    
!------------------------- Local Variables -----------------------------
    integer m, icnst, mp
    real(r8) rstd
!-----------------------------------------------------------------------
1         format(a8,' ', i3,i3,l5,i3,e16.5)
2         format(i4,' ',a8)
3         format(i4,' ',a8,' ',i4,a10)

    write(iulog,*) ' ' 
    write(iulog,*) '---- Water isotopes tracer configuration ----'
    write(iulog,*) 'name      W  S   tag f     Rstd'

    do icnst = 1, wtrc_ncnst
      m    = wtrc_indices(icnst)
      rstd = wtrc_get_rstd(iwspec(m))
        
      write(iulog,1) cnst_name(m),iwater(m), iwspec(m),  &
        iwistag(m), int(wiso_get_fisub(iwspec(m))), rstd
    end do

    write(iulog,*) ' ' 
    write(iulog,*) '---- Water isotopes - surface vapor ----'
    write(iulog,*) '  m  name'
    do icnst = 1, wtrc_nsrfvap
      m = wtrc_indices(wtrc_iasrfvap(icnst))
      write(iulog,2) m, cnst_name(m)
    end do
    
    write(iulog,*) ' ' 
    write(iulog,*) '---- Water isotopes - surface precip ----'
    write(iulog,*) '  m  name       p   name'
    do icnst = 1, wtrc_nsrfpcp
      m  = wtrc_indices(wtrc_iasrfpcp(icnst))
      mp = wtrc_srfpcp_indices(wtrc_types(wtrc_iasrfpcp(icnst)), int((wtrc_iasrfpcp(icnst) - 1) / pwtype) + 1)
      
      if (mp == -1) then
        write(iulog,2) m, cnst_name(m)
      else
        write(iulog,3) m, cnst_name(m), mp, trim(pbuf_get_field_name(mp))
      end if
    end do

    write(iulog,*) '---------------------------------------------'
    write(iulog,*) ' ' 

    write(iulog,*) 'water_tracer_model     :  ', water_tracer_model
    write(iulog,*) 'trace_water            :  ', trace_water
    write(iulog,*) 'wisotope               :  ', wisotope
    write(iulog,*) 'wtrc_detrain_in_macrop :  ', wtrc_detrain_in_macrop
    write(iulog,*) 'wtrc_check_total_h2o   :  ', wtrc_check_total_h2o
    write(iulog,*) 'wtrc_check_show_types  :  ', wtrc_check_show_types
    write(iulog,*) 'wtrc_alpha_kinetic     :  ', wtrc_alpha_kinetic
    write(iulog,*) 'wtrc_add_cvprecip      :  ', wtrc_add_cvprecip
    write(iulog,*) 'wtrc_add_stprecip      :  ', wtrc_add_stprecip
    write(iulog,*) 'wtrc_warn_only         :  ', wtrc_warn_only
    write(iulog,*) 'wtrc_lh2oadj           :  ', wtrc_lh2oadj
    write(iulog,*) 'wtrc_niter             :  ', wtrc_niter
    write(iulog,*) 'wtrc_qchkmin           :  ', wtrc_qchkmin
    write(iulog,*) 'wtrc_qmin              :  ', wtrc_qmin
    write(iulog,*) 'wtrc_fixed_alpha       :  ', wtrc_fixed_alpha
    write(iulog,*) 'wtrc_fixed_rstd        :  ', wtrc_fixed_rstd

    write(iulog,*) ' ' 
    write(iulog,*) '---------------------------------------------'
    write(iulog,*) ' ' 
 
!
    return
  end subroutine wtrc_setup_diag


!=======================================================================
  subroutine wtrc_rescale(q,ncol)
!-----------------------------------------------------------------------
! Purpose: Ensures tracer water mass is exactly the same as prognostic
! Author: David Noone <dcn@caltech.edu> - Mon Mar  8 16:22:30 PST 2004
!-----------------------------------------------------------------------
    use constituents,   only: cnst_get_ind
!---------------------------- Arguments --------------------------------
  integer, intent(in)     :: ncol
  real(r8), intent(inout) :: q(pcols,pver,pcnst)
!-----------------------------------------------------------------------
  real(r8) qerr(pcols,pver)
  real(r8) rat
  integer i,k,m
  integer ixvap                 ! prognostic water species
  integer mbase     ! tracer base for tracer m
  integer mprog     ! prognostic base for tracer m
  integer ivap, itype
!-----------------------------------------------------------------------
!
! Apply appropriate scaling
!
    do itype = 1, pwtype
      qerr(:,:)  = 0._r8

      m          = wtrc_bulk_indices(itype)
    
      do ivap = 1, wtrc_ntype(itype)
        mbase     = wtrc_indices(wtrc_iatype(ivap, itype))
        qerr(:,:) = qerr(:,:) + q(:,:,mbase)
      end do
    
!
! Compute the error
!
      qerr(:,:) = qerr(:,:) - q(:,:,m)
    
!
! Compute tracer ratio, consistent with tracers, then
! apply to total (prognostic) mass
!
      do ivap = 1, wtrc_ntype(itype)
        mbase = wtrc_indices(wtrc_iatype(ivap, itype))

        do k = 1, pver
          do i = 1, ncol
            rat = wtrc_ratio(iwspec(mbase), q(i,k,m), q(i,k,mbase))
!!            q(i,k,m) =  rat*q(i,k,mprog)     ! direct scale
            q(i,k,m) =  q(i,k,m) - rat*qerr(i,k)   ! scale error 
          end do
        end do
      end do
    end do
!
    return
  end subroutine wtrc_rescale

!=======================================================================
  function wtrc_ratio(ispec,qtrc,qtot)
!-----------------------------------------------------------------------
! Purpose: Compute tracer ratio from masses, with numerical checks
! Author David Noone <dcn@colorado.edu> - Sat Jul  3 18:52:40 MDT 2004
!-----------------------------------------------------------------------
    integer,intent(in)   :: ispec       ! isotopic species
    real(r8),intent(in)  :: qtrc        ! tracer water mass
    real(r8),intent(in)  :: qtot        ! "base" water mass
    real(r8)             :: wtrc_ratio  ! return value
!-----------------------------------------------------------------------
!    real(r8) :: qtiny = 1.e-16_r8    ! bigger makes scheme more stable
    real(r8) :: qtiny = 1.e-22_r8   ! smaller makes scheme more accurate
!-----------------------------------------------------------------------
 
!    if (qtot < qtiny) then
!      wtrc_ratio = 1._r8
!    else if (qtot > 0._r8) then
!      wtrc_ratio = qtrc / (qtot+qtiny)
!    else
!      wtrc_ratio = qtrc / (qtot-qtiny)
!    end if

    if (abs(qtot) < wtrc_qmin) then
      wtrc_ratio = wtrc_get_rstd(ispec)
    else
      wtrc_ratio = qtrc / qtot
    end if

    return
  end function wtrc_ratio


!=======================================================================
  subroutine wtrc_ratio_all(ncol,q,rat)
!-----------------------------------------------------------------------
! Computes ratios for all water tracers
!-----------------------------------------------------------------------
    use constituents,   only: cnst_get_ind
!---------------------------- Arguments --------------------------------
   integer , intent(in)   :: ncol
   real(r8), intent(in)   :: q(pcols,pver,pcnst)
   real(r8), intent(out)  :: rat(pcols,pver,pcnst)
!-----------------------------------------------------------------------
   integer                :: i          ! column index
   integer                :: k          ! vertical index
   integer                :: m          ! constituent index for isotopic water
   integer                :: mbase      ! constituent index for "regular" water
   integer                :: iwset      ! water set index
   integer                :: itype      ! water type index
!-----------------------------------------------------------------------
    rat(:,:,:) = 0._r8

    ! Compute ratios based on total water set in the model. By definition,
    ! this should be the first water set.
    do iwset = 1, wtrc_nwset
      do itype = 1, pwtype
        
        m     = wtrc_iawset(itype, iwset)
        mbase = wtrc_iawset(itype, 1)
  
        do k = 1, pver
          do i = 1, ncol
  
            ! TBD: Should the strc routine or the wiso routine be used? Both would seem
            ! to need to be modified to handle smaller mixing ratios.
            rat(i, k, m) = wtrc_ratio(iwspec(m), q(i, k, m), q(i, k, mbase))
          end do
        end do
      end do
    end do

    return
  end subroutine wtrc_ratio_all


!=======================================================================
  function wtrc_get_alpha(q, tk, ispec, isrctype, idsttype, rhclc, porqh, kin)
!-----------------------------------------------------------------------
! Purpose: Retrieve the fractionation for a process that goes from
!          the source water type to the destination water type.
!
! Author: Chuck Bardeen
!
! Modified for RH use by Jesse Nusbaumer
!
!-----------------------------------------------------------------------
    use water_isotopes, only: isph2o
    use water_types,    only: wtype_get_alpha
    use wv_saturation,  only: qsat_water

    real(r8), intent(in)            :: q              ! water vapor mixing ratio (kg/kg)
    real(r8), intent(in)            :: tk             ! temperature (K)
    integer, intent(in)             :: ispec          ! isotope species index
    integer, intent(in)             :: isrctype       ! source water type index
    integer, intent(in)             :: idsttype       ! destination water type index
    logical, intent(in)             :: rhclc          ! Does RH need to be calculated?
    real(r8), intent(in)            :: porqh          ! Either pressure (Pa) or RH, depending on rhclc
    logical, intent(in), optional   :: kin            ! Use kinetic fractionation? (over-rides namelist) 
    real(r8)                        :: wtrc_get_alpha ! return fractionation factor
    
    real(r8)                        :: pmid1(1)       ! pressure (Pa)
    real(r8)                        :: tk1(1)         ! temperature (K)
    real(r8)                        :: es(1)          ! sat. vapor pressure
    real(r8)                        :: qs(1)          ! saturation specific humidity
    real(r8)                        :: gam(1)         ! (l/cp)*(d(qs)/dt) !<-Probably not needed... - JN
    logical                         :: alpkin         ! Kinetic fractionation logical
!-----------------------------------------------------------------------

    if(present(kin)) then         !Is kinetic fractionation being specified?
      alpkin = kin       
    else                          !If not, then let namelist determine
      alpkin = wtrc_alpha_kinetic
    end if

    if(rhclc) then !If relative humidity is being calculated:
      pmid1(1) = porqh
      tk1(1)   = tk

      call qsat_water(tk1, pmid1, es, qs, gam)

      if (wisotope) then
        wtrc_get_alpha = wtype_get_alpha(ispec, isrctype, idsttype, tk, q / qs(1), alpkin)
      else
        wtrc_get_alpha = wtrc_fixed_alpha(ispec)
      end if

    else   !If relative humidity is already known:

      if (wisotope) then
        wtrc_get_alpha = wtype_get_alpha(ispec, isrctype, idsttype, tk, porqh, alpkin)
      else
        wtrc_get_alpha = wtrc_fixed_alpha(ispec)
      end if

    end if
    
    return
  end function wtrc_get_alpha

  
!=======================================================================
  function wtrc_get_rstd(ispec)
!-----------------------------------------------------------------------
! Get the standard isotopic ratio, but adjust the H216O fraction so that
! the result for the model configuration sum to 1.
!
! NOTE: If wisotope is false, then a ratio of 1 is returned for h216o and
! 0 for all other species.
!-----------------------------------------------------------------------
    use water_isotopes, only: wiso_get_rstd, isph2o

    integer, intent(in)     :: ispec
    real(r8)                :: wtrc_get_rstd

!-----------------------------------------------------------------------
    
    ! Either return the true standard isotopic ratio or the default value
    ! if isotopes are not enabled.
    if (wisotope) then
      wtrc_get_rstd = wiso_get_rstd(ispec)
    else
      wtrc_get_rstd = wtrc_fixed_rstd(ispec)
    endif

    return
  end function wtrc_get_rstd






! NOTE: The following routines are to generate test cases and/or perform
! checks. These routines have been brought over from the previous version
! of the CAM isotope code, but are not currently being used. As
! development progresses, they may be removed if they are not found to
! be needed. The are similar to some checks that are incorporated into\
! the new code.
#undef NOCHECK          /* DEBUG: define to bypass conservation checks */
#undef QCHKMESS         /* DEBUG: define to give "sucess" message for checks */
#undef QCHKTERM         /* DEBUG: define to terminate when "qchk" failed */


!=======================================================================
   subroutine wtrc_init_qpert(qpert)
!-----------------------------------------------------------------------
! Initialize constituent perturbation to something (smow?)
!-----------------------------------------------------------------------
    use water_isotopes, only: wiso_get_rstd

    real(r8), intent(inout) :: qpert(pcols,pcnst)

    real(r8) rat
    
    integer  ivap
!-----------------------------------------------------------------------
    qpert(:,:) = 0._r8
    
    if (trace_water) then

      do ivap = 1, wtrc_ntype(iwtvap)

        rat = wiso_get_rstd(wtrc_species(wtrc_iatype(ivap, iwtvap)))
!
        qpert(:,wtrc_indices(wtrc_iatype(ivap,iwtvap))) = rat * qpert(:,wtrc_bulk_indices(iwtvap))
      end do
 
    end if
!   
     return
   end subroutine wtrc_init_qpert


!=======================================================================
  subroutine wtrc_qchk3(subr, vname, ncol, q, qmag0)
!-----------------------------------------------------------------------
! Checks that all tracers are Rstd*prognostic 
! (used for debugging with no fractionation)
!-----------------------------------------------------------------------
    use water_isotopes, only: wiso_get_rstd
    use water_types,    only: wtype_suffix
    use constituents,   only: pcnst, cnst_get_ind

!---------------------------- Arguments --------------------------------
    character(len=*),intent(in) :: subr   ! name of calling subroutine
    character(len=*),intent(in) :: vname  ! name of variable
    integer , intent(in) :: ncol          ! number of columns to scan
    real(r8), intent(in) :: q(pcols,pver,pcnst)  ! tracers
    real(r8), intent(in), optional :: qmag0      ! minimum magnitude of qprg
!------------------------- Local Variables -----------------------------
    real(r8) rstd
    integer m, mbase, ispec, itype, icnst
!-----------------------------------------------------------------------

    do ispec = 1, pwtspec
      rstd = wiso_get_rstd(ispec)
      
      write(*,'(a40,5i6,g16.6)') 'WTRC_QCHK3 ('//trim(subr)//') - tracers:',wtrc_iaspec(:, ispec),rstd
      
      do itype = 1, wtrc_nspec(ispec)
        icnst = wtrc_iaspec(itype, ispec)
        
        m     = wtrc_indices(icnst)
        mbase = wtrc_bulk_indices(icnst)

        call wtrc_qchk2(subr, trim(vname)//trim(wtype_suffix(wtrc_types(icnst))), ncol, q(:,:,m),  rstd*q(:,:,mbase),  qmag0)
      end do
    end do
!
    return
  end subroutine wtrc_qchk3


!=======================================================================
  subroutine wtrc_qchk2(subr,vname,ncol,qtrc,qprg,qmag0)
!-----------------------------------------------------------------------
! Purpose: Check the tracer water mass equal the prognostic
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 19:00:15 MDT 2003
!-----------------------------------------------------------------------

!---------------------------- Arguments --------------------------------
    character(len=*),intent(in) :: subr   ! name of calling subroutine
    character(len=*),intent(in) :: vname  ! name of variable
    integer , intent(in) :: ncol          ! number of columns to scan
    real(r8), intent(in) :: qtrc(pcols,pver)   ! tracer water
    real(r8), intent(in) :: qprg(pcols,pver)   ! prognostic water
    real(r8), intent(in), optional :: qmag0      ! minimum magnitude of qprg
!------------------------- Local Variables -----------------------------
    real(r8) etest                        ! test variable
    real(r8) qmag
    real(r8) qdw, etw                   ! worst values
    integer nbad                        ! number of bad values found
    integer i,k
!-----------------------------------------------------------------------
    nbad = 0
!    qmag = 0._r8
    qmag = wtrc_qchkmin
    qdw = 0._r8
    etw = 0._r8
    if (present(qmag0)) qmag = qmag0
!
    do k = 1, pver
      do i = 1, ncol
       if (wtrc_qchk_one(qtrc(i,k),qprg(i,k),etest,qmag) > 0) then
#ifdef QCHKTERM
          write(iulog,1) 'WTRC_QCHK2: '//'('//trim(subr)//'.'//trim(vname)//')'// &
          ' q(m,1):',i,k,qtrc(i,k),qprg(i,k),etest
1       format(a36,2i4,2e12.3,e12.4)
#endif
         etw = max(etw, abs(etest))
         qdw = max(qdw, abs(qtrc(i,k)-qprg(i,k)))
         nbad = nbad + 1
       end if
      end do
    end do
!
    if (nbad /= 0) then
        write(iulog,*) 'WTRC_QCHK2: '//'('//trim(subr)//'.'//trim(vname)//')',&
              ' *** WARNING - chunk tracers /= Q =',nbad,etw,qdw
#ifdef QCHKTERM 
        call endrun('QCHK2 failed.')
#endif
    else

#ifdef QCHKMESS   /* print a sucess message */
        write(iulog,*) 'WTRC_QCHK2: '//'('//trim(subr)//'.'//trim(vname)//')',&
              ' All OK.'
#endif
    end if
!
    return
  end subroutine wtrc_qchk2

!=======================================================================
  subroutine wtrc_qchk1(subr,vname,ncol,qtrc,qprg,qmag0)
!-----------------------------------------------------------------------
! Purpose: Check the tracer water mass equal the prognostic
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 19:00:15 MDT 2003
!-----------------------------------------------------------------------
!---------------------------- Arguments --------------------------------
    character(len=*),intent(in) :: subr   ! name of calling subroutine
    character(len=*),intent(in) :: vname  ! name of variable
    integer , intent(in) :: ncol    ! number of columns to scan
    real(r8), intent(in) :: qtrc(pcols)   ! tracer water
    real(r8), intent(in) :: qprg(pcols)   ! prognostic water
    real(r8), intent(in),optional :: qmag0 ! minimum q for fail test
!------------------------- Local Variables -----------------------------
    real(r8) etest                        ! test variable
    real(r8) qmag
    real(r8) qdw, etw                   ! worst values
    integer nbad            ! number of bad values found
    integer i
!-----------------------------------------------------------------------
    nbad = 0
!    qmag = 0._r8
    qmag = wtrc_qchkmin
    qdw = 0._r8
    etw = 0._r8
    if (present(qmag0)) qmag = qmag0
!
    do i = 1, ncol
       if (wtrc_qchk_one(qtrc(i),qprg(i),etest,qmag) > 0) then
#ifdef QCHKTERM
          write(iulog,1) 'WTRC_QCHK1: '//'('//trim(subr)//'.'//trim(vname)//')'// &
          ' q(m,1):',i,qtrc(i),qprg(i),etest
1       format(a40,i4,'    ',2e12.4,e12.4)
#endif
         etw = max(etw, abs(etest))
         qdw = max(qdw, abs(qtrc(i)-qprg(i)))
         nbad = nbad + 1
       end if
    end do
!
    if (nbad /= 0) then
        write(iulog,*) 'WTRC_QCHK1: '//'('//trim(subr)//'.'//trim(vname)//')',&
              ' *** WARNING - chunk tracers /= Q =',nbad,etw,qdw
#ifdef QCHKTERM               /* terminate */
        call endrun('QCHK1 failed.')
#endif
    else
#ifdef QCHKMESS   /* print a sucess message */
        write(iulog,*) 'WTRC_QCHK1: '//'('//trim(subr)//'.'//trim(vname)//')',&
              ' All OK.'
#endif
    end if
!
    return
  end subroutine wtrc_qchk1


!=======================================================================
  function wtrc_qchk_one(qtrc,qprg,etest,qmag)
!-----------------------------------------------------------------------
! Purpose: Check the one tracer water mass equal the prognostic
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 19:00:15 MDT 2003
!-----------------------------------------------------------------------
!    real(r8), parameter :: elimit = 1.0e-16 ! precision required
!    real(r8), parameter :: elimit = 1.0e-14 ! precision required (q1q2 fails)
!    real(r8), parameter :: elimit = 1.0e-12 ! precision required
    real(r8), parameter :: elimit = 1.0e-10_r8 ! precision required
!    real(r8), parameter :: qmin = 1.0e-18 ! precision required
    real(r8), parameter :: qmin = 1.e-26_r8 ! precision required
!---------------------------- Arguments --------------------------------
    real(r8), intent(in) :: qtrc          ! tracer water
    real(r8), intent(in) :: qprg          ! prognostic water
    real(r8), intent(in) :: qmag          ! difference limit
    real(r8), intent(out) :: etest        ! test variable
    real(r8)              qdiff,qmabs
    integer wtrc_qchk_one     ! return value
!-----------------------------------------------------------------------
    wtrc_qchk_one = 0

#ifdef NOCHECK
!
! By-pass all checking if running with fractionation
!
#else
    qmabs = max(abs(qprg),abs(qtrc))
    qdiff = abs(qtrc - qprg)
    if (qmabs > qmin .and. qdiff > qmag) then
      etest = qdiff / qmabs
      if (etest > elimit) then
         wtrc_qchk_one = 1
#ifdef QCHKTERM   /* if going to fail, write all diagnostics */
         write(iulog,'(a36,4e16.9)') '(WTRC_QCHK) FAILED TEST:',qtrc,qprg,qdiff,etest
#endif
      end if
    end if
#endif
!
    return
  end function wtrc_qchk_one

!=======================================================================
  subroutine wtrc_chkdelta(subr, ncol, q)
!-----------------------------------------------------------------------
! Checks the delta values of a 2d array (lon,lev)
!-----------------------------------------------------------------------
    use constituents, only: cnst_get_ind
    use water_isotopes, only: wiso_delta
!---------------------------- Arguments --------------------------------
    character(len=*), intent(in) :: subr  ! name of calling routine/message
    integer , intent(in) :: ncol      ! number of columns to scan
    real(r8), intent(in) :: q(pcols,pver,pcnst) ! tracer quantity
    real(r8) del, delbad
    integer i,k,m
    integer mbase     ! prognostic base for tracer m
    integer nbad
    real(r8) qbad
    integer  icnst
!-----------------------------------------------------------------------
!
! Apply appropriate scaling
!
    do icnst = 1, wtrc_ncnst
      m     = wtrc_indices(icnst)
      mbase = wtrc_bulk_indices(wtrc_types(icnst))
!
      delbad = 0._r8
      qbad   = 0._r8
      nbad   = 0._r8
    
      do k = 1, pver
        do i = 1, ncol
          del = wiso_delta(iwspec(m), q(i,k,m), q(i,k,mbase))
        
          if (abs(del) > 1001._r8) then
            nbad = nbad + 1
            if (abs(del) > abs(delbad)) then
               qbad = q(i,k,mbase)
               delbad = del
            endif
!!            call endrun('(wtrc_chkdelta) Stopped.')
          end if
        end do
      end do

!!      if (nbad > 0) then
!!        write(*,*) trim(subr)//' Bad delta values for m=',m
!!        write(*,*) 'nbad = ',nbad, '  worst=',delbad,qbad
!!      end if
    end do
    
    return
  end subroutine wtrc_chkdelta

!=======================================================================
  subroutine wtrc_check(subr, ncol, q)
!-----------------------------------------------------------------------
! Checks H2O tracer (ice, liquid and vapour) is the same as the prognostic
! (optionallly adjust)
!-----------------------------------------------------------------------
    use constituents,   only: cnst_get_ind, qmin
    use water_isotopes, only: wiso_delta
    use water_types,    only: wtype_names
!---------------------------- Arguments --------------------------------
    character(len=*), intent(in) :: subr        ! name of calling routine/message
    integer , intent(in) :: ncol                ! number of columns to scan
    real(r8), intent(inout) :: q(pcols,pver,pcnst) ! tracer quantity (optionally scaled)
    real(r8) dwatr(wtrc_nspec(isph2o),pcols,pver)
    integer i,j,k,m
    integer iw(2)     ! indices of worst values
    integer mbase                       ! prognostic base for tracer m
    logical lerrors
    integer icnst, ixh2oq, itype
!-----------------------------------------------------------------------
    real(r8) :: qtol = 1.e-12_r8   ! qmin(Q) = 1.e-12, qmin(L,I) = 0.)
!    real(r8) :: qtol = 1.e-15    ! qmin(Q) = 1.e-12, qmin(L,I) = 0.)
!    real(r8) :: qtol = 1.e-17    ! Too strict for zm_evap (review?)
!-----------------------------------------------------------------------
!
    lerrors = .false.
!
! Get all differences
!
    do j = 1, wtrc_nspec(isph2o)
      icnst = wtrc_iaspec(j, isph2o)
      
      m     = wtrc_indices(icnst)
      mbase = wtrc_bulk_indices(wtrc_types(icnst))

! There should be one and only one H2O vapor that is not a tagged tracer.      
      if (wtrc_is_vap(m) .and. .not. wtrc_is_tagged(m)) then
        ixh2oq = m
      end if
      
      do k = 1, pver
        do i = 1, ncol
          dwatr(j,i,k) = q(i,k,mbase) - q(i,k,m)
        end do
      end do
    end do
!
! Send reports
!
    if (.not. wtrc_warn_only) then
      do k = 1, pver
        do i = 1, ncol
          if (abs(dwatr(j,i,k)) > max(qtol,qmin(wtrc_bulk_indices(iwtvap)))) then
            write(iulog,1) i,k,q(i,k,wtrc_bulk_indices(iwtvap)), q(i,k,ixh2oq), dwatr(j,i,k) 
1           format (2i3,3e20.10)
          end if
        end do
      end do
    end if
!
    do i = 1, wtrc_nspec(isph2o)
      icnst = wtrc_iaspec(i, isph2o)
      itype = wtrc_types(icnst)
      
      if (count(abs(dwatr(i,:ncol,:))>max(qtol,qmin(wtrc_bulk_indices(itype)))) > 0) then
         write(*,*) '(wtrc_check) ' // wtype_names(itype) // ' differences: '//trim(subr)
         iw = maxloc(abs(dwatr(i,:ncol,:)))
         write(*,2) count(abs(dwatr(i,:ncol,:))>qmin(wtrc_bulk_indices(itype))),iw(1),iw(2),dwatr(i,iw(1),iw(2)), q(iw(1),iw(2),wtrc_indices(icnst))
         if (.not. wtrc_warn_only) &
             call endrun('wtrc_check: ' // wtrc_names(icnst) // ' check failed.')
      endif
    end do
 2  format(i4,' point(s), worst (i,k):',2i5,2e16.6)
!
! Do any adjustments
!
    if (wtrc_lh2oadj) then
!!        write(*,*) 'Applying rescaling to tracers to prohibit drift.'
        call wtrc_rescale(q,ncol)
    end if
!
#ifdef QCHKMESS
   write(iulog,*) '(wtrc_check) all OK: '//trim(subr)
#endif
!
    return
  end subroutine wtrc_check


end module water_tracers
