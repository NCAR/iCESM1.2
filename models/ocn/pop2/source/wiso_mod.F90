!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module wiso_mod

!BOP
! !MODULE: wiso_mod
!
!  Module for Water ISOtopes (WISOs)
!  
!  The delta values of water isotopes are treated as passive tracers.
!     e.g., delta = (R/Rstd - 1) * 1000
!  The units of concentration  for these tracers are per mil
!  The units of surface fluxes for these tracers are per mil * cm/s
!
! !DESCRIPTION:
! ECB adding the case 'coupler' fields
! Add marginal sea balancing for surface tracer fluxes
!     Jiang Zhu (jzhu47@wisc.edu), 20150802
!
! !REVISION HISTORY:
!  SVN:$Id: $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks, only: nx_block, ny_block, block, get_block
   use domain_size, only: max_blocks_clinic, km, nx_global, ny_global
   use domain, only: nblocks_clinic, distrb_clinic, blocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use constants, only: c0, c1, c1000, p5, p001, eps
   use io
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field
   use passive_tracer_tools, only: forcing_monthly_every_ts, &
       ind_name_pair, tracer_read, read_field
   use broadcast
   use netcdf
   use ms_balance_wiso, only: ms_balancing_wiso
   use forcing_fields,  only: PREC_16O_F_nf_ind, PREC_18O_F_nf_ind,     &
                              PREC_HDO_F_nf_ind, EVAP_16O_F_nf_ind,     & 
                              EVAP_18O_F_nf_ind, EVAP_HDO_F_nf_ind,     & 
                              MELT_16O_F_nf_ind, MELT_18O_F_nf_ind,     & 
                              MELT_HDO_F_nf_ind, ROFF_16O_F_nf_ind,     & 
                              ROFF_18O_F_nf_ind, ROFF_HDO_F_nf_ind,     & 
                              IOFF_16O_F_nf_ind, IOFF_18O_F_nf_ind,     & 
                              IOFF_HDO_F_nf_ind, ROCE_16O_nf_ind,       &
                              ROCE_18O_nf_ind,   ROCE_HDO_nf_ind

   implicit none
   save

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: &
       wiso_tracer_cnt,    &
       wiso_init,          &
       wiso_set_sflux,     &
       wiso_tavg_forcing,  &
       wiso_write_restart, &
       wiso_tracer_ref_val

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracers
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
       wiso_tracer_cnt = 2

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      d18o_ind =  1,  & ! Delta18O
      dD_ind   =  2     ! DeltaD

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(wiso_tracer_cnt) :: &
      ind_name_table = (/ &
      ind_name_pair(d18o_ind, 'Delta18O'), &
      ind_name_pair(dD_ind, 'DeltaD') /)

!-----------------------------------------------------------------------
!  mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

   logical (log_kind), dimension(:,:,:), allocatable :: &
      LAND_MASK

!-----------------------------------------------------------------------
!  restoring climatology for surface values
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:), allocatable, target :: &
      SURF_D18O         ! observed annual mean surface D18O

!-----------------------------------------------------------------------
!  forcing related variables
!-----------------------------------------------------------------------

   character(char_len) :: &
      wiso_formulation        ! how to calculate flux (ocmip or coupler)

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK             ! temp array for interpolate_forcing output

   type(forcing_monthly_every_ts) :: &
      prec_file,            & ! wiso flux from precip, if read from file
      dair_file,            & ! delta value of marine air, if read from file
      rh_file,              & ! near-surface relative humidity, if read from file
      d18or_file              ! runoff d18o

!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::   &
      WISO_SFLUX_TAVG

   integer (int_kind) :: &
      tavg_WISO_PREC,         & ! tavg id for wiso flux from PREC
      tavg_WISO_DAIR,         & ! tavg id for delta val of marine air
      tavg_WISO_d18OR,        & ! tavg id for runoff d18o
      tavg_WISO_dHDOR,        & ! tavg id for runoff dHDo
      tavg_WISO_RH,           & ! tavg id for relative humidity
      tavg_PREC_16O_F,        & ! tavg id for 16O PREC_F
      tavg_PREC_18O_F,        & ! tavg id for 18O PREC_F
      tavg_PREC_HDO_F,        & ! tavg id for HDO PREC_F
      tavg_EVAP_16O_F,        & ! tavg id for 16O EVAP_F
      tavg_EVAP_18O_F,        & ! tavg id for 18O EVAP_F
      tavg_EVAP_HDO_F,        & ! tavg id for HDO EVAP_F
      tavg_MELT_16O_F,        & ! tavg id for 16O MELT_F
      tavg_MELT_18O_F,        & ! tavg id for 18O MELT_F
      tavg_MELT_HDO_F,        & ! tavg id for HDO MELT_F
      tavg_ROFF_16O_F,        & ! tavg id for 16O ROFF_F
      tavg_ROFF_18O_F,        & ! tavg id for 18O ROFF_F
      tavg_ROFF_HDO_F,        & ! tavg id for HDO ROFF_F
      tavg_IOFF_16O_F,        & ! tavg id for 16O IOFF_F
      tavg_IOFF_18O_F,        & ! tavg id for 18O IOFF_F
      tavg_IOFF_HDO_F,        & ! tavg id for HDO IOFF_F
      tavg_Roce_16O,          & ! tavg id for Roce/Rstd 16O 
      tavg_Roce_18O,          & ! tavg id for Roce/Rstd 18O 
      tavg_Roce_HDO,          & ! tavg id for Roce/Rstd HDO 
      tavg_STF_PREC_18O,      & ! tavg id for STF flux from PREC
      tavg_STF_PREC_HDO,      & ! tavg id for STF flux from PREC
      tavg_STF_EVAP_18O,      & ! tavg id for STF flux from EVAP
      tavg_STF_EVAP_HDO,      & ! tavg id for STF flux from EVAP
      tavg_STF_MELT_18O,      & ! tavg id for STF flux from MELT
      tavg_STF_MELT_HDO,      & ! tavg id for STF flux from MELT
      tavg_STF_ROFF_18O,      & ! tavg id for STF flux from ROFF
      tavg_STF_ROFF_HDO,      & ! tavg id for STF flux from ROFF
      tavg_STF_IOFF_18O,      & ! tavg id for STF flux from IOFF
      tavg_STF_IOFF_HDO,      & ! tavg id for STF flux from IOFF
      tavg_STF_18O,           & ! tavg id for total STF flux
      tavg_STF_HDO,           & ! tavg id for total STF flux
      tavg_Delta18O_EVAP,     & ! tavg id for d18o flux from EVAP
      tavg_DeltaD_EVAP,       & ! tavg id for dD flux from EVAP
      tavg_Delta18O_ROFF,     & ! tavg id for d18o flux from ROFF
      tavg_DeltaD_ROFF,       & ! tavg id for dD flux from ROFF
      tavg_Delta18O_IOFF,     & ! tavg id for d18o flux from IOFF
      tavg_DeltaD_IOFF          ! tavg id for dD flux from IOFF

!-----------------------------------------------------------------------
!  average surface tracer value related variables
!  used as reference value for virtual flux computations
!-----------------------------------------------------------------------

   logical (log_kind), dimension(wiso_tracer_cnt) :: &
      vflux_flag                ! which tracers get virtual fluxes applied

   integer (int_kind) :: &
      comp_surf_avg_flag        ! time flag id for computing average
                                ! surface tracer values

   real (r8), dimension(wiso_tracer_cnt) :: &
      surf_avg                  ! average surface tracer values

   real (r8) :: &
      wiso_sflx_correction_18o, & ! surface flux correction constant for R18O (kg/m2/s)
      wiso_sflx_correction_hdo    ! surface flux correction constant for RHDO (kg/m2/s)


!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: wiso_sflux_timer

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: wiso_init
! !INTERFACE:

 subroutine wiso_init(init_ts_file_fmt, read_restart_filename, &
                     tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize wiso tracer module. This involves setting metadata, reading
!  the modules namelist and setting initial conditions.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: char_blank, delim_fmt
   use prognostic, only: curtime, oldtime
   use grid, only: KMT, n_topo_smooth, fill_points
   use prognostic, only: tracer_field
   use timers, only: get_timer
   use passive_tracer_tools, only: init_forcing_monthly_every_ts, &
       rest_read_tracer_block, file_read_tracer_block
   use time_management, only: eval_time_flag, check_time_flag, init_time_flag, &
       freq_opt_never, freq_opt_nyear, freq_opt_nmonth

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(wiso_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(nx_block,ny_block,km,wiso_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
character(*), parameter :: sub_name = 'wiso_mod:wiso_init'

   character(char_len) :: &
      init_wiso_option,        & ! option for initialization of wiso 
      init_wiso_init_file,     & ! filename for option 'file'
      init_wiso_init_file_fmt, & ! file format for option 'file'
      wiso_comp_surf_avg_freq_opt     ! choice for freq of comp_surf_avg

   integer (int_kind) :: &
      n,                      & ! index for looping over tracers
      k,                      & ! index for looping over depth levels
      iblock,                 & ! index for looping over blocks
      nml_error                 ! namelist i/o error flag

   type(tracer_read), dimension(wiso_tracer_cnt) :: &
      tracer_init_ext           ! namelist variable for initializing tracers

   type(tracer_read) :: &
      fwiso_flux_prec,        & ! wiso flux from precipitation
      fwiso_flux_dair,        & ! delta value of marine air
      fwiso_flux_d18or,       & ! runoff d18o
      fwiso_flux_rh             ! near-surface relative humidity

   integer (int_kind) :: &
      freq_opt, freq,         & ! args for init_time_flag
      wiso_comp_surf_avg_freq_iopt,& ! choice for freq of comp_surf_avg
      wiso_comp_surf_avg_freq        ! choice for freq of comp_surf_avg

   logical (log_kind) :: &
      wiso_use_nml_surf_vals         ! do namelist surf values override values from restart file

   real (r8) :: &
      surf_avg_d18o_const, surf_avg_dD_const

   namelist /wiso_nml/ &
      init_wiso_option, init_wiso_init_file, init_wiso_init_file_fmt, &
      tracer_init_ext,  &
      wiso_comp_surf_avg_freq_opt, wiso_comp_surf_avg_freq,  &
      wiso_use_nml_surf_vals, surf_avg_d18o_const, surf_avg_dD_const, &
      wiso_formulation, fwiso_flux_prec, fwiso_flux_dair, fwiso_flux_d18or, fwiso_flux_rh, &
      wiso_sflx_correction_18o, wiso_sflx_correction_hdo

   character (char_len) ::  &
      wiso_restart_filename      ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_forcing_monthly_every_ts(prec_file)
   call init_forcing_monthly_every_ts(dair_file)
   call init_forcing_monthly_every_ts(rh_file)
   call init_forcing_monthly_every_ts(d18or_file)

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   do n = 1, wiso_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = 'per mil'
      tracer_d_module(n)%tend_units = 'per mil/s'
      tracer_d_module(n)%flux_units = 'per mil*cm/s'
   end do

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_wiso_option        = 'unknown'
   init_wiso_init_file     = 'unknown'
   init_wiso_init_file_fmt = 'bin'

   do n = 1, wiso_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   wiso_formulation = 'model'

   fwiso_flux_prec%filename     = '/glade/p/work/jiangzhu/data/inputdata/CAM3.PERinput.gx3v7.12mon.nc'
   fwiso_flux_prec%file_varname = 'Pd18O'
   fwiso_flux_prec%scale_factor = c1
   fwiso_flux_prec%default_val  = c0
   fwiso_flux_prec%file_fmt     = 'nc'

   fwiso_flux_dair%filename     = fwiso_flux_prec%filename
   fwiso_flux_dair%file_varname = 'd18OA'
   fwiso_flux_dair%scale_factor = c1
   fwiso_flux_dair%default_val  = -12.
   fwiso_flux_dair%file_fmt     = 'nc'

   fwiso_flux_rh%filename     = fwiso_flux_prec%filename
   fwiso_flux_rh%file_varname = 'RELHUM'
   fwiso_flux_rh%scale_factor = c1
   fwiso_flux_rh%default_val  = 0.75
   fwiso_flux_rh%file_fmt     = 'nc'

   fwiso_flux_d18or%filename     = fwiso_flux_prec%filename
   fwiso_flux_d18or%file_varname = 'd18OR'
   fwiso_flux_d18or%scale_factor = c1
   fwiso_flux_d18or%default_val  = -12.
   fwiso_flux_d18or%file_fmt     = 'nc'

   wiso_comp_surf_avg_freq_opt        = 'never'
   wiso_comp_surf_avg_freq            = 1
   wiso_use_nml_surf_vals             = .false.
   surf_avg_d18o_const         = 0.173845_r8
   surf_avg_dD_const         = 0.173845_r8

   wiso_sflx_correction_18o = c0
   wiso_sflx_correction_hdo = c0

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=wiso_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(sub_name, 'wiso_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ sub_name)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_wiso_option, master_task)
   call broadcast_scalar(init_wiso_init_file, master_task)
   call broadcast_scalar(init_wiso_init_file_fmt, master_task)

   do n = 1, wiso_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(wiso_formulation, master_task)

   call broadcast_scalar(fwiso_flux_prec%filename, master_task)
   call broadcast_scalar(fwiso_flux_prec%file_varname, master_task)
   call broadcast_scalar(fwiso_flux_prec%scale_factor, master_task)
   call broadcast_scalar(fwiso_flux_prec%default_val, master_task)
   call broadcast_scalar(fwiso_flux_prec%file_fmt, master_task)

   prec_file%input = fwiso_flux_prec

   call broadcast_scalar(fwiso_flux_dair%filename, master_task)
   call broadcast_scalar(fwiso_flux_dair%file_varname, master_task)
   call broadcast_scalar(fwiso_flux_dair%scale_factor, master_task)
   call broadcast_scalar(fwiso_flux_dair%default_val, master_task)
   call broadcast_scalar(fwiso_flux_dair%file_fmt, master_task)

   dair_file%input = fwiso_flux_dair

   call broadcast_scalar(fwiso_flux_rh%filename, master_task)
   call broadcast_scalar(fwiso_flux_rh%file_varname, master_task)
   call broadcast_scalar(fwiso_flux_rh%scale_factor, master_task)
   call broadcast_scalar(fwiso_flux_rh%default_val, master_task)
   call broadcast_scalar(fwiso_flux_rh%file_fmt, master_task)

   rh_file%input = fwiso_flux_rh

   call broadcast_scalar(fwiso_flux_d18or%filename, master_task)
   call broadcast_scalar(fwiso_flux_d18or%file_varname, master_task)
   call broadcast_scalar(fwiso_flux_d18or%scale_factor, master_task)
   call broadcast_scalar(fwiso_flux_d18or%default_val, master_task)
   call broadcast_scalar(fwiso_flux_d18or%file_fmt, master_task)

   d18or_file%input = fwiso_flux_d18or

   call broadcast_scalar(wiso_comp_surf_avg_freq_opt, master_task)
   call broadcast_scalar(wiso_comp_surf_avg_freq, master_task)
   call broadcast_scalar(wiso_use_nml_surf_vals, master_task)
   call broadcast_scalar(surf_avg_d18o_const, master_task)
   call broadcast_scalar(surf_avg_dD_const, master_task)

   call broadcast_scalar(wiso_sflx_correction_18o, master_task)
   call broadcast_scalar(wiso_sflx_correction_hdo, master_task)

!-----------------------------------------------------------------------
!  set variables immediately dependent on namelist variables
!-----------------------------------------------------------------------

   select case (wiso_comp_surf_avg_freq_opt)
   case ('never')
      wiso_comp_surf_avg_freq_iopt = freq_opt_never
   case ('nyear')
      wiso_comp_surf_avg_freq_iopt = freq_opt_nyear
   case ('nmonth')
      wiso_comp_surf_avg_freq_iopt = freq_opt_nmonth
   case default
      call document(sub_name, 'wiso_comp_surf_avg_freq_opt', wiso_comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'unknown wiso_comp_surf_avg_freq_opt')
   end select

  call init_time_flag('wiso_comp_surf_avg', comp_surf_avg_flag, &
     default=.false., freq_opt=wiso_comp_surf_avg_freq_iopt,  &
     freq=wiso_comp_surf_avg_freq, owner='wiso_init')

!-----------------------------------------------------------------------
!  namelist consistency checking
!-----------------------------------------------------------------------

   if (wiso_use_nml_surf_vals .and. wiso_comp_surf_avg_freq_iopt /= freq_opt_never) then
      call document(sub_name, 'wiso_use_nml_surf_vals', wiso_use_nml_surf_vals)
      call document(sub_name, 'wiso_comp_surf_avg_freq_opt', wiso_comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'wiso_use_nml_surf_vals can only be .true. if ' /&
                           &/ ' wiso_comp_surf_avg_freq_opt is never')
   endif

!-----------------------------------------------------------------------
!  initialize virtual flux flag array
!-----------------------------------------------------------------------

   vflux_flag = .false.
   vflux_flag(d18o_ind) = .true.
   vflux_flag(dD_ind) = .true.

!-----------------------------------------------------------------------
!   initialize tracers
!-----------------------------------------------------------------------

   select case (init_wiso_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d WISO ratios set to all zeros'
          write(stdout,delim_fmt)
      endif

      if (wiso_use_nml_surf_vals) then
         surf_avg = c0
         surf_avg(d18o_ind) = surf_avg_d18o_const
         surf_avg(dD_ind) = surf_avg_dD_const
      else
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                            TRACER_MODULE(:,:,1,:,curtime,:))
      endif

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      wiso_restart_filename = char_blank

      if (init_wiso_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(sub_name, 'no restart file to read WISOs from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ sub_name)
         endif
         wiso_restart_filename = read_restart_filename
         init_wiso_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         wiso_restart_filename = trim(init_wiso_init_file)

      endif

      call rest_read_tracer_block(init_wiso_init_file_fmt, &
                                  wiso_restart_filename,   &
                                  tracer_d_module,        &
                                  TRACER_MODULE)

      if (wiso_use_nml_surf_vals) then
         surf_avg = c0
         surf_avg(d18o_ind) = surf_avg_d18o_const
         surf_avg(dD_ind) = surf_avg_dD_const
      else
         call extract_surf_avg(init_wiso_init_file_fmt, &
                               wiso_restart_filename)
      endif

      call eval_time_flag(comp_surf_avg_flag) ! evaluates time_flag(comp_surf_avg_flag)%value via time_to_do

      if (check_time_flag(comp_surf_avg_flag)) &
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                            TRACER_MODULE(:,:,1,:,curtime,:))

   case ('file')

      call document(sub_name, 'WISOs being read from separate file')

      call file_read_tracer_block(init_wiso_init_file_fmt, &
                                  init_wiso_init_file,     &
                                  tracer_d_module,        &
                                  ind_name_table,         &
                                  tracer_init_ext,        &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, wiso_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'wiso_init: error in fill_points')
                  return
               endif
            end do
         end do
      endif

      if (wiso_use_nml_surf_vals) then
         surf_avg = c0
         surf_avg(d18o_ind) = surf_avg_d18o_const
         surf_avg(dD_ind) = surf_avg_dD_const
      else
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                            TRACER_MODULE(:,:,1,:,curtime,:))
      endif

   case default
      call document(sub_name, 'init_wiso_option', init_wiso_option)
      call exit_POP(sigAbort, 'unknown init_wiso_option')

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
   do n = 1, wiso_tracer_cnt
      do k = 1, km
         where (k > KMT(:,:,iblock))
            TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
            TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
         end where
      end do
   end do
   end do

!-----------------------------------------------------------------------
!  allocate and initialize LAND_MASK (true for ocean points)
!-----------------------------------------------------------------------

   allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
   LAND_MASK = merge(.true., .false., KMT > 0)

   call get_timer(wiso_sflux_timer, 'WISO_SFLUX', 1, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call wiso_init_tavg
   call wiso_init_sflux

!-----------------------------------------------------------------------
!EOC

 end subroutine wiso_init

!***********************************************************************
!BOP
! !IROUTINE: extract_surf_avg
! !INTERFACE:

 subroutine extract_surf_avg(init_wiso_init_file_fmt, &
                             wiso_restart_filename)

! !DESCRIPTION:
!  Extract average surface values from restart file.
!
! !REVISION HISTORY:
!  by Jiaxu Zhang 1/6/2012

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_wiso_init_file_fmt, & ! file format (bin or nc)
      wiso_restart_filename      ! file name for restart file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   type (datafile) ::&
      restart_file    ! io file descriptor

   integer (int_kind) :: &
      n               ! tracer index

   character (char_len) :: &
      short_name      ! tracer name temporaries

!-----------------------------------------------------------------------

   surf_avg = c0

   restart_file = construct_file(init_wiso_init_file_fmt, &
                                 full_name=trim(wiso_restart_filename), &
                                 record_length=rec_type_dbl, &
                                 recl_words=nx_global*ny_global)

   do n = 1, wiso_tracer_cnt
      if (vflux_flag(n)) then
         short_name = 'surf_avg_' /&
                   &/ ind_name_table(n)%name
         call add_attrib_file(restart_file, trim(short_name), surf_avg(n))
      endif
   end do

   call data_set(restart_file, 'open_read')

   do n = 1, wiso_tracer_cnt
      if (vflux_flag(n)) then
         short_name = 'surf_avg_' /&
                   &/ ind_name_table(n)%name
         call extract_attrib_file(restart_file, trim(short_name), surf_avg(n))
      endif
   end do

   call data_set (restart_file, 'close')

   call destroy_file (restart_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine extract_surf_avg

!***********************************************************************
!BOP
! !IROUTINE: wiso_init_tavg
! !INTERFACE:

 subroutine wiso_init_tavg

! !DESCRIPTION:
!  Define tavg fields not automatically handled by the base model.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      var_cnt             ! how many tavg variables are defined

!-----------------------------------------------------------------------

   var_cnt = 0

   call define_tavg_field(tavg_WISO_PREC,'WISO_PREC',2,           &
                          long_name='WISO fluxes from precip',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_WISO_DAIR,'WISO_DAIR',2,           &
                          long_name='Delta value of marine air',&
                          units='per mil', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_WISO_RH,'WISO_RH',2,           &
                          long_name='Near-surface relative humidity',&
                          units='percent', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_WISO_d18OR,'WISO_d18OR',2,           &
                          long_name='d18O in precip.',&
                          units='per mil', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_Delta18O_EVAP,'Delta18O_EVAP',2,           &
                          long_name='Delta18O fluxes from evapor',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DeltaD_EVAP,'DeltaD_EVAP',2,           &
                          long_name='DeltaD fluxes from evapor',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_Delta18O_ROFF,'Delta18O_ROFF',2, &
                          long_name='Delta18O fluxes from liq runoff',           &
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DeltaD_ROFF,'DeltaD_ROFF',2, &
                          long_name='DeltaD fluxes from liq runoff',           &
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_Delta18O_IOFF,'Delta18O_IOFF',2, &
                          long_name='Delta18O fluxes from ice runoff',           &
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DeltaD_IOFF,'DeltaD_IOFF',2, &
                          long_name='DeltaD fluxes from ice runoff',           &
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

! passed in from cpl in x2o:
   call define_tavg_field(tavg_PREC_16O_F,'PREC_16O_F',2,           &
                          long_name='H216O Precipitation Flux from Cpl (rain+snow)',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_PREC_18O_F,'PREC_18O_F',2,           &
                          long_name='H218O Precipitation Flux from Cpl (rain+snow)',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_PREC_HDO_F,'PREC_HDO_F',2,           &
                          long_name='HDO Precipitation Flux from Cpl (rain+snow)',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_EVAP_16O_F,'EVAP_16O_F',2,           &
                          long_name='H216O Evaporation Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_EVAP_18O_F,'EVAP_18O_F',2,           &
                          long_name='H218O Evaporation Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_EVAP_HDO_F,'EVAP_HDO_F',2,           &
                          long_name='HDO Evaporation Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_MELT_16O_F,'MELT_16O_F',2,           &
                          long_name='H216O Melt Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_MELT_18O_F,'MELT_18O_F',2,           &
                          long_name='H218O Melt Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_MELT_HDO_F,'MELT_HDO_F',2,           &
                          long_name='HDO Melt Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ROFF_16O_F,'ROFF_16O_F',2,           &
                          long_name='H216O LIQ Runoff Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ROFF_18O_F,'ROFF_18O_F',2,           &
                          long_name='H218O LIQ Runoff Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ROFF_HDO_F,'ROFF_HDO_F',2,           &
                          long_name='HDO LIQ Runoff Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_IOFF_16O_F,'IOFF_16O_F',2,           &
                          long_name='H216O ICE Runoff Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_IOFF_18O_F,'IOFF_18O_F',2,           &
                          long_name='H218O ICE Runoff Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_IOFF_HDO_F,'IOFF_HDO_F',2,           &
                          long_name='HDO ICE Runoff Flux from Coupler',&
                          units='kg m-2 s-1', grid_loc='2110')
   var_cnt = var_cnt+1

! passed to cpl in o2x:
   call define_tavg_field(tavg_Roce_16O,'Roce_16O',2,           &
                          long_name='Surface R/Rstd for 16O',&
                          units='mol/mol', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_Roce_18O,'Roce_18O',2,           &
                          long_name='Surface R/Rstd for 18O',&
                          units='mol/mol', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_Roce_HDO,'Roce_HDO',2,           &
                          long_name='Surface R/Rstd for HDO',&
                          units='mol/mol', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_PREC_18O,'STF_PREC_18O',2,           &
                          long_name='STF fluxes from 18O precip from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_EVAP_18O,'STF_EVAP_18O',2,           &
                          long_name='STF fluxes from 18O evap from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_MELT_18O,'STF_MELT_18O',2,           &
                          long_name='STF fluxes from 18O melt from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_ROFF_18O,'STF_ROFF_18O',2,           &
                          long_name='STF fluxes from 18O liq runoff from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_IOFF_18O,'STF_IOFF_18O',2,           &
                          long_name='STF fluxes from 18O ice runoff from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_PREC_HDO,'STF_PREC_HDO',2,           &
                          long_name='STF fluxes from HDO precip from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_EVAP_HDO,'STF_EVAP_HDO',2,           &
                          long_name='STF fluxes from HDO evap from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_MELT_HDO,'STF_MELT_HDO',2,           &
                          long_name='STF fluxes from HDO melt from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_ROFF_HDO,'STF_ROFF_HDO',2,           &
                          long_name='STF fluxes from HDO liq runoff from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_IOFF_HDO,'STF_IOFF_HDO',2,           &
                          long_name='STF fluxes from HDO ice runoff from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_18O,'STF_18O',2,           &
                          long_name='total STF fluxe for 18O from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_STF_HDO,'STF_HDO',2,           &
                          long_name='total STF fluxe for HDO from driver',&
                          units='per mil*cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------

   allocate(WISO_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   WISO_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine wiso_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: wiso_init_sflux
! !INTERFACE:

 subroutine wiso_init_sflux

! !USES:

   use forcing_tools, only: find_forcing_times
   use named_field_mod, only: named_field_get_index, named_field_register, &
       named_field_set

! !DESCRIPTION:
!  Initialize surface flux computations for wiso tracer module.
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: sub_name = 'wiso_mod:wiso_init_sflux'

   integer (int_kind) :: &
      n,                 & ! index for looping over tracers
      iblock               ! index for looping over blocks

   real (r8), dimension (nx_block,ny_block) :: WORK

   real (r8), dimension (nx_block,ny_block,12,max_blocks_clinic), target :: &
      WORK_READ            ! temporary space to read in fields


!-----------------------------------------------------------------------
!  read forcing (if required)
!  otherwise, use values passed in
!-----------------------------------------------------------------------

   select case (wiso_formulation)

   case ('restoring')

!-----------------------------------------------------------------------
!  allocate various wiso allocatable module variables and read in
!-----------------------------------------------------------------------
   
   allocate( SURF_D18O(nx_block,ny_block,max_blocks_clinic) )

      call read_field('nc', &
                      '/glade/u/home/brady/CESM1/iCESM/inputdata/calculated_d18O_gx3v7.nc',   &
                      'd18o', SURF_D18O)
!                     '/glade/home/jizhang/scripts/wiso.test01/ForcingFiles/calculated_d18O_gx3v7.nc',   &

   case ('ocmip')

!-----------------------------------------------------------------------
!  allocate space for interpolate_forcing
!-----------------------------------------------------------------------

      allocate(INTERP_WORK(nx_block,ny_block,max_blocks_clinic,1))

!-----------------------------------------------------------------------
!  first, read file of wiso flux from precipitation
!-----------------------------------------------------------------------

      allocate(prec_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(prec_file%input%file_fmt, &
                      prec_file%input%filename, &
                      prec_file%input%file_varname, &
                      WORK_READ)
      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         prec_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            prec_file%DATA(:,:,iblock,1,n) = c0
         prec_file%DATA(:,:,iblock,1,n) = &
            prec_file%DATA(:,:,iblock,1,n) * prec_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(prec_file%data_time, &
                              prec_file%data_inc, prec_file%interp_type, &
                              prec_file%data_next, prec_file%data_time_min_loc, &
                              prec_file%data_update, prec_file%data_type)

!-----------------------------------------------------------------------
!  next, read file of delta values of marine air
!-----------------------------------------------------------------------

      allocate(dair_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(dair_file%input%file_fmt, &
                      dair_file%input%filename, &
                      dair_file%input%file_varname, &
                      WORK_READ)
      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         dair_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            dair_file%DATA(:,:,iblock,1,n) = c0
         dair_file%DATA(:,:,iblock,1,n) = &
            dair_file%DATA(:,:,iblock,1,n) * dair_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(dair_file%data_time, &
                              dair_file%data_inc, dair_file%interp_type, &
                              dair_file%data_next, dair_file%data_time_min_loc, &
                              dair_file%data_update, dair_file%data_type)

!-----------------------------------------------------------------------
!  third, read relative humidity file
!-----------------------------------------------------------------------

      allocate(rh_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(rh_file%input%file_fmt, &
                      rh_file%input%filename, &
                      rh_file%input%file_varname, &
                      WORK_READ)
      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         rh_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            rh_file%DATA(:,:,iblock,1,n) = c0
         rh_file%DATA(:,:,iblock,1,n) = &
            rh_file%DATA(:,:,iblock,1,n) * rh_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(rh_file%data_time, &
                              rh_file%data_inc, rh_file%interp_type, &
                              rh_file%data_next, rh_file%data_time_min_loc, &
                              rh_file%data_update, rh_file%data_type)

!-----------------------------------------------------------------------
!  finally, read d18OR file
!-----------------------------------------------------------------------

      allocate(d18or_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(d18or_file%input%file_fmt, &
                      d18or_file%input%filename, &
                      d18or_file%input%file_varname, &
                      WORK_READ)
      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         d18or_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            d18or_file%DATA(:,:,iblock,1,n) = c0
         d18or_file%DATA(:,:,iblock,1,n) = &
            d18or_file%DATA(:,:,iblock,1,n) * d18or_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(d18or_file%data_time, &
                              d18or_file%data_inc, d18or_file%interp_type, &
                              d18or_file%data_next, d18or_file%data_time_min_loc, &
                              d18or_file%data_update, d18or_file%data_type)

!-----------------------------------------------------------------------
! register and set surface ratio fields to pass to coupler:  will have to move this...
!-----------------------------------------------------------------------

   call named_field_register('Roce_16O', ROCE_16O_nf_ind)
   call named_field_register('Roce_18O', ROCE_18O_nf_ind)
   call named_field_register('Roce_HDO', ROCE_HDO_nf_ind)
   !$OMP PARALLEL DO PRIVATE(iblock,WORK)
   do iblock=1,nblocks_clinic
      WORK = 1.0_r8 
      call named_field_set(ROCE_16O_nf_ind, iblock, WORK)
      call named_field_set(ROCE_18O_nf_ind, iblock, WORK)
      call named_field_set(ROCE_HDO_nf_ind, iblock, WORK)
   end do
   !$OMP END PARALLEL DO

   case ('model')

      if (my_task == master_task) then
         write(stdout,*)  &
            ' Using fields from coupler forcing for calculating WISO flux'
      endif

!-----------------------------------------------------------------------
!  get indices for surface forcing terms
!-----------------------------------------------------------------------

      call named_field_get_index('PREC_16O_F', PREC_16O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('PREC_18O_F', PREC_18O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('PREC_HDO_F', PREC_HDO_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('EVAP_16O_F', EVAP_16O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('EVAP_18O_F', EVAP_18O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('EVAP_HDO_F', EVAP_HDO_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('MELT_16O_F', MELT_16O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('MELT_18O_F', MELT_18O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('MELT_HDO_F', MELT_HDO_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('ROFF_16O_F', ROFF_16O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('ROFF_18O_F', ROFF_18O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('ROFF_HDO_F', ROFF_HDO_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('IOFF_16O_F', IOFF_16O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('IOFF_18O_F', IOFF_18O_F_nf_ind, &
                               exit_on_err=.true.)
      call named_field_get_index('IOFF_HDO_F', IOFF_HDO_F_nf_ind, &
                               exit_on_err=.true.)

!-----------------------------------------------------------------------
!  register and set surface ratio fields to pass to coupler:  
!-----------------------------------------------------------------------

      call named_field_register('Roce_16O', ROCE_16O_nf_ind)
      call named_field_register('Roce_18O', ROCE_18O_nf_ind)
      call named_field_register('Roce_HDO', ROCE_HDO_nf_ind)
   !$OMP PARALLEL DO PRIVATE(iblock,WORK)
      do iblock=1,nblocks_clinic
         WORK = c1
         call named_field_set(ROCE_16O_nf_ind, iblock, WORK)
         call named_field_set(ROCE_18O_nf_ind, iblock, WORK)
         call named_field_set(ROCE_HDO_nf_ind, iblock, WORK)
      end do
   !$OMP END PARALLEL DO

   case default
      call document(sub_name, 'wiso_formulation', wiso_formulation)

      call exit_POP(sigAbort, &
                    'wiso_init_sflux: Unknown value for wiso_formulation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine wiso_init_sflux

!!***********************************************************************
!!BOP
!! !IROUTINE: read_1dvar_cdf
!! !INTERFACE:
!
! subroutine read_1dvar_cdf(ncid, data_dimid, varname, data, stat)
!
!! !DESCRIPTION:
!!  Subroutine to read in a single 1D variable from a netCDF file
!!  that is supposed to be on a particular dimension
!!
!! !REVISION HISTORY:
!!  same as module
!
!! !USES:
!
!! !INPUT PARAMETERS:
!
!   integer (int_kind), intent(in) :: &
!      ncid              ! netCDF file id
!
!   integer (int_kind), dimension(1), intent(in) :: &
!      data_dimid        ! netCDF dimension id that all data should have
!
!   character (len=*), intent(in) :: &
!      varname           ! name of variable being read
!
!! !OUTPUT PARAMETERS:
!
!   real (r8), dimension(:), intent(out) :: &
!      data              ! where data is going
!
!   integer (int_kind), intent(out) :: &
!      stat              ! status of netCDF call
!!EOP
!!BOC
!!-----------------------------------------------------------------------
!!  local variables
!!-----------------------------------------------------------------------
!
!   integer (int_kind) :: &
!      varid,          & ! netCDF variable id
!      ndims             ! number of dimensions for varid
!
!   integer (int_kind), dimension(1) :: &
!      dimid             ! netCDF dimension id
!
!!-----------------------------------------------------------------------
!
!   stat = nf90_inq_varid(ncid, varname, varid)
!   if (stat /= 0) then
!      write(stdout,*) 'nf_inq_varid for ', trim(varname), ' : ', nf90_strerror(stat)
!      return
!   endif
!
!   stat = nf90_inquire_variable(ncid, varid, ndims=ndims)
!   if (stat /= 0) then
!      write(stdout,*) 'nf_inq_varndims for ', trim(varname), ' : ', nf90_strerror(stat)
!      return
!   endif
!   if (ndims /= 1) then
!      write(stdout,*) 'ndims /= 1 for ', trim(varname)
!      return
!   endif
!
!   stat = nf90_inquire_variable(ncid, varid, dimids=dimid)
!   if (stat /= 0) then
!      write(stdout,*) 'nf_inq_vardimid for ', trim(varname), ' : ', nf90_strerror(stat)
!      return
!   endif
!   if (dimid(1) /= data_dimid(1)) then
!      write(stdout,*) 'dimid mismatch for ', trim(varname)
!      return
!   endif
!
!   stat = nf90_get_var(ncid, varid, data)
!   if (stat /= 0) then
!      write(stdout,*) 'nf_get_var_double for ', trim(varname), ' : ', nf90_strerror(stat)
!      return
!   endif
!
!!-----------------------------------------------------------------------
!!EOC
!
! end subroutine read_1dvar_cdf

!***********************************************************************
!BOP
! !IROUTINE: wiso_set_sflux
! !INTERFACE:

 subroutine wiso_set_sflux(SST,SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Set surface fluxes for wiso tracers
!  Restoring surface value for WISO

! !REVISION HISTORY:
!  by Jiaxu Zhang (jzhang76@wisc.edu), 11/23/2011
! tests for PREC_*_F(ecb, brady@ucar.edu; 4-5-20123)
! adds in case ('model')

! !USES:

   use shr_const_mod
   use constants, only: field_loc_center, field_type_scalar, p5, fwmass_to_fwflux
   use grid,   only: dz
   use forcing_tools, only: update_forcing_data, interpolate_forcing
   use timers, only: timer_start, timer_stop
   use time_management, only: seconds_in_day, check_time_flag, thour00
   use forcing_fields, only: EVAP_F, ROFF_F, PREC_F, IOFF_F, MELT_F
   use named_field_mod, only: named_field_get, named_field_set

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      SST          ! sea surface temperature (C)

   real (r8), dimension(nx_block,ny_block,wiso_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,wiso_tracer_cnt,max_blocks_clinic), &
         intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: sub_name = 'wiso_mod:wiso_set_sflux'

   integer (int_kind) :: &
      iblock             ! block index

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      STF_PREC,        & ! computed wiso flux from prec
      STF_EVAP,        & ! computed wiso flux from evap
      STF_MELT,        & ! computed wiso flux from melt
      STF_ROFF,        & ! computed wiso flux from liq runoff
      STF_IOFF,        & ! computed wiso flux from ice runoff
      alpha_wv,        & ! water to vapor fractionation factor
      ATM_VALS,        & ! marine air tracer values
      RH,              & ! near-surface relative humidity
      d18OP,           & ! d18o in prec flux
      dHDOP,           & ! dhdo in prec flux
      d18OE,           & ! d18o in evap flux
      dHDOE,           & ! dhdo in evap flux
      d18OM,           & ! d18o in melt flux
      dHDOM,           & ! dhdo in melt flux
      d18OR,           & ! d18o in roff flux
      dHDOR,           & ! dhdo in roff flux
      d18OI,           & ! d18o in ioff flux
      dHDOI              ! dhdo in roff flux

   real (r8), dimension(nx_block,ny_block) :: &
      SURF_VALS,       & ! filtered surface tracer values
      dEVAP,           & ! delta values in EVAP
      PREC_16O_F,      & ! Isotopic Precip from coupler 
      PREC_18O_F,      & ! Isotopic Precip 
      PREC_HDO_F,      & ! Isotopic Precip 
      EVAP_16O_F,      & ! Isotopic Evap from coupler 
      EVAP_18O_F,      & ! Isotopic Evap 
      EVAP_HDO_F,      & ! Isotopic Evap 
      MELT_16O_F,      & ! Isotopic melt from coupler 
      MELT_18O_F,      & ! Isotopic melt 
      MELT_HDO_F,      & ! Isotopic melt 
      ROFF_16O_F,      & ! Isotopic ROFF from coupler 
      ROFF_18O_F,      & ! Isotopic ROFF 
      ROFF_HDO_F,      & ! Isotopic ROFF 
      iOFF_16O_F,      & ! Isotopic IOFF from coupler 
      iOFF_18O_F,      & ! Isotopic IOFF 
      iOFF_HDO_F,      & ! Isotopic IOFF
      ROCE_16O,        & ! Isotopic Surface Ratio 
      ROCE_18O,        & ! Isotopic Surface Ratio 
      ROCE_HDO,        & ! Isotopic Surface Ratio 
      STF_PREC_18O,    & ! Delta Prec flux from driver  
      STF_PREC_HDO,    & ! Delta Prec flux from driver 
      STF_EVAP_18O,    & ! Delta Evap flux from driver 
      STF_EVAP_HDO,    & ! Delta Evap flux from driver 
      STF_MELT_18O,    & ! Delta Melt flux from driver 
      STF_MELT_HDO,    & ! Delta Melt flux from driver 
      STF_ROFF_18O,    & ! Delta Runoff flux from driver 
      STF_ROFF_HDO,    & ! Delta Runoff flux from driver 
      STF_IOFF_18O,    & ! Delta Runoff flux from driver 
      STF_IOFF_HDO,    & ! Delta Runoff flux from driver 
      STF_18O,         & ! Total delta flux from driver
      STF_HDO,         & ! Total delta flux from driver
      WORK,            & ! for work computation  
      WISO_surf_sat      ! WISO surface saturation (either d18o or dD) (per mil)

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
      WISO_SFLUX_16O,  &
      WISO_SFLUX_18O,  &
      WISO_SFLUX_HDO

   character (char_len) :: &
      tracer_data_label          ! label for what is being updated

   character (char_len), dimension(1) :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(1) :: &
      tracer_bndy_loc,          &! location and field type for ghost
      tracer_bndy_type           !    cell updates

   real (r8) ::          &
      wiso_restore_tau,         &! restoring time scale
      wiso_restore_rtau          ! reciprocal of restoring time scale

   real (r8) ::        &
      avg_16O_F,       &
      avg_18O_F,       &
      avg_HDO_F

!-----------------------------------------------------------------------
!  local parameters
!-----------------------------------------------------------------------

   real (r8), parameter :: &
!     Rstd_16O = SHR_CONST_VSMOW_16O,      & 
!     Rstd_18O = SHR_CONST_VSMOW_18O,      & 
!     Rstd_HDO = SHR_CONST_VSMOW_D,      & 
!     Rstd_HHO = SHR_CONST_VSMOW_H,      & 
      Rstd_16O = c1,        & 
      Rstd_18O = c1,        &
      Rstd_HDO = c1,        &
      Rstd_HHO = c1,        &
      K  = 6.0e-3,          & ! kinetic fractionation parameter
      a1 = 0.9884,          & ! coefficient for alpha_wv
      a2 = 1.025e-4,        & ! coefficient for alpha_wv
      a3 = -3.57e-7           ! coefficient for alpha_wv

!-----------------------------------------------------------------------

   call timer_start(wiso_sflux_timer)

   if (check_time_flag(comp_surf_avg_flag))  &
      call comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR)

   select case (wiso_formulation)

   case ('restoring')

      wiso_restore_tau  = 30.
      wiso_restore_rtau = c1 / (seconds_in_day * wiso_restore_tau)

      do iblock = 1, nblocks_clinic

         STF_MODULE(:,:,d18o_ind,iblock) = &
                                  (SURF_D18O(:,:,iblock) - &
                                  SURF_VALS_CUR(:,:,d18o_ind,iblock))*        &
                                  wiso_restore_rtau*dz(1)

         STF_MODULE(:,:,dD_ind,iblock) = &
                                  (SURF_D18O(:,:,iblock) - &
                                  SURF_VALS_CUR(:,:,dD_ind,iblock))*        &
                                  wiso_restore_rtau*dz(1)

      end do

   case ('ocmip')

      do iblock = 1, nblocks_clinic
         STF_PREC(:,:,iblock) = c0
         ATM_VALS(:,:,iblock) = c0
         RH(:,:,iblock) = c0
         d18OP(:,:,iblock) = c0 
      end do

!-----------------------------------------------------------------------
!  Interpolate forcing data if necessary
!-----------------------------------------------------------------------

!       if (thour00 >= prec_file%data_update) then
!          tracer_data_names = prec_file%input%file_varname
!          tracer_bndy_loc   = field_loc_center
!          tracer_bndy_type  = field_type_scalar
!          tracer_data_label = 'WISO flux from prec'
!          call update_forcing_data(          prec_file%data_time,   &
!               prec_file%data_time_min_loc,  prec_file%interp_type, &
!               prec_file%data_next,          prec_file%data_update, &
!               prec_file%data_type,          prec_file%data_inc,    &
!               prec_file%DATA(:,:,:,:,1:12), prec_file%data_renorm, &
!               tracer_data_label,            tracer_data_names,     &
!               tracer_bndy_loc,              tracer_bndy_type,      &
!               prec_file%filename,           prec_file%input%file_fmt)
!       endif
!       call interpolate_forcing(INTERP_WORK, &
!            prec_file%DATA(:,:,:,:,1:12), &
!            prec_file%data_time,         prec_file%interp_type, &
!            prec_file%data_time_min_loc, prec_file%interp_freq, &
!            prec_file%interp_inc,        prec_file%interp_next, &
!            prec_file%interp_last,       0)
!       STF_PREC = INTERP_WORK(:,:,:,1)

      if (thour00 >= dair_file%data_update) then
         tracer_data_names = dair_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Delta value of marine air'
         call update_forcing_data(          dair_file%data_time,   &
              dair_file%data_time_min_loc,  dair_file%interp_type, &
              dair_file%data_next,          dair_file%data_update, &
              dair_file%data_type,          dair_file%data_inc,    &
              dair_file%DATA(:,:,:,:,1:12), dair_file%data_renorm, &
              tracer_data_label,            tracer_data_names,     &
              tracer_bndy_loc,              tracer_bndy_type,      &
              dair_file%filename,           dair_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK, &
           dair_file%DATA(:,:,:,:,1:12), &
           dair_file%data_time,         dair_file%interp_type, &
           dair_file%data_time_min_loc, dair_file%interp_freq, &
           dair_file%interp_inc,        dair_file%interp_next, &
           dair_file%interp_last,       0)
      ATM_VALS = INTERP_WORK(:,:,:,1)
 
      if (thour00 >= rh_file%data_update) then
         tracer_data_names = rh_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Near-surface relative humidity'
         call update_forcing_data(        rh_file%data_time,   &
              rh_file%data_time_min_loc,  rh_file%interp_type, &
              rh_file%data_next,          rh_file%data_update, &
              rh_file%data_type,          rh_file%data_inc,    &
              rh_file%DATA(:,:,:,:,1:12), rh_file%data_renorm, &
              tracer_data_label,          tracer_data_names,     &
              tracer_bndy_loc,            tracer_bndy_type,      &
              rh_file%filename,           rh_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK, &
           rh_file%DATA(:,:,:,:,1:12), &
           rh_file%data_time,         rh_file%interp_type, &
           rh_file%data_time_min_loc, rh_file%interp_freq, &
           rh_file%interp_inc,        rh_file%interp_next, &
           rh_file%interp_last,       0)
      RH = INTERP_WORK(:,:,:,1)
      RH = RH * .01                  ! Readin RH is in percent

      if (thour00 >= d18or_file%data_update) then
         tracer_data_names = d18or_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Delta value of marine air'
         call update_forcing_data(          d18or_file%data_time,   &
              d18or_file%data_time_min_loc,  d18or_file%interp_type, &
              d18or_file%data_next,          d18or_file%data_update, &
              d18or_file%data_type,          d18or_file%data_inc,    &
              d18or_file%DATA(:,:,:,:,1:12), d18or_file%data_renorm, &
              tracer_data_label,            tracer_data_names,     &
              tracer_bndy_loc,              tracer_bndy_type,      &
              d18or_file%filename,           d18or_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK, &
           d18or_file%DATA(:,:,:,:,1:12), &
           d18or_file%data_time,         d18or_file%interp_type, &
           d18or_file%data_time_min_loc, d18or_file%interp_freq, &
           d18or_file%interp_inc,        d18or_file%interp_next, &
           d18or_file%interp_last,       0)
      d18OP = INTERP_WORK(:,:,:,1)


 !$OMP PARALLEL DO PRIVATE(iblock,SURF_VALS,dEVAP,                   &
 !$OMP                      ROCE_16O,ROCE_18O,ROCE_HDO,              &
 !$OMP        EVAP_16O_F,EVAP_18O_F,EVAP_HDO_F,STF_EVAP_18O,STF_PREC_18O,  &
 !$OMP                      PREC_16O_F,PREC_18O_F,PREC_HDO_F)
   do iblock = 1, nblocks_clinic
          call named_field_get(PREC_16O_F_nf_ind, iblock, PREC_16O_F)
          call named_field_get(PREC_18O_F_nf_ind, iblock, PREC_18O_F)
          call named_field_get(PREC_HDO_F_nf_ind, iblock, PREC_HDO_F)
          call named_field_get(EVAP_16O_F_nf_ind, iblock, EVAP_16O_F)
          call named_field_get(EVAP_18O_F_nf_ind, iblock, EVAP_18O_F)
          call named_field_get(EVAP_HDO_F_nf_ind, iblock, EVAP_HDO_F)

         where (ATM_VALS(:,:,iblock) < -16.000_r8) &
            ATM_VALS(:,:,iblock) = -16.000_r8
         where (RH(:,:,iblock) > 0.8_r8) &
            RH(:,:,iblock) = 0.8_r8

         where (LAND_MASK(:,:,iblock))
!note this is set up such that dHDO tracer = d18O tracer, fluxes the same.
            STF_PREC(:,:,iblock) = PREC_F(:,:,iblock)*fwmass_to_fwflux*d18OP(:,:,iblock)
            
!let's still compute it, but not use it.. goes to history output now for
!testing and comparison...
            where(PREC_F(:,:,iblock) /= c0 ) 
             STF_PREC_18O(:,:) = 1000.0_r8*(PREC_18O_F(:,:)/(Rstd_16O*Rstd_18O) - PREC_F(:,:,iblock))*fwmass_to_fwflux
             STF_PREC_HDO(:,:) = 1000.0_r8*(PREC_HDO_F(:,:)/(Rstd_HHO*Rstd_HDO) - PREC_F(:,:,iblock))*fwmass_to_fwflux
            elsewhere 
             STF_PREC_18O(:,:) = c0 
             STF_PREC_HDO(:,:) = c0 
            end where

            WISO_SFLUX_TAVG(:,:,1,iblock) = STF_PREC(:,:,iblock)
            WISO_SFLUX_TAVG(:,:,2,iblock) = ATM_VALS(:,:,iblock)
            WISO_SFLUX_TAVG(:,:,3,iblock) = RH(:,:,iblock)
!           where(PREC_F(:,:,iblock) /= c0) 
!             d18OP(:,:,iblock) = PREC_18O_F(:,:)/PREC_F(:,:,iblock) 
!           elsewhere
!             d18OP(:,:,iblock) = c0 
!           end where
   
            WISO_SFLUX_TAVG(:,:,4,iblock) = d18OP(:,:,iblock)

            alpha_wv(:,:,iblock) = a1 + a2*SST(:,:,iblock) + a3*(SST(:,:,iblock)**2)

            SURF_VALS = p5*(SURF_VALS_OLD(:,:,d18o_ind,iblock) + &
                            SURF_VALS_CUR(:,:,d18o_ind,iblock))
            dEVAP = (1-K) / (1 - RH(:,:,iblock)) * & 
                    (alpha_wv(:,:,iblock) * (SURF_VALS + 1.e3) - &
                     RH(:,:,iblock) * (ATM_VALS(:,:,iblock) + 1.e3)) - 1.e3
             STF_EVAP(:,:,iblock) = EVAP_F(:,:,iblock)*fwmass_to_fwflux * dEVAP

!from driver, for comparison...
            where(EVAP_F(:,:,iblock) /= c0 ) 
! note that EVAP_18O_F comes from drver not scaled by Rstd_16O, so not divided by
! here:
             STF_EVAP_18O(:,:) = 1000.0_r8*(EVAP_18O_F(:,:)/Rstd_18O - EVAP_F(:,:,iblock))*fwmass_to_fwflux
             STF_EVAP_HDO(:,:) = 1000.0_r8*(EVAP_HDO_F(:,:)/Rstd_HDO - EVAP_F(:,:,iblock))*fwmass_to_fwflux
            elsewhere 
             STF_EVAP_18O(:,:) = c0 
            end where

            WISO_SFLUX_TAVG(:,:,5,iblock) = STF_EVAP(:,:,iblock)
            STF_ROFF(:,:,iblock) = ROFF_F(:,:,iblock)*fwmass_to_fwflux * d18OP(:,:,iblock)
            WISO_SFLUX_TAVG(:,:,7,iblock) = STF_ROFF(:,:,iblock)
!redundant
            WISO_SFLUX_TAVG(:,:,20,iblock) = STF_ROFF(:,:,iblock)

            STF_MODULE(:,:,d18o_ind,iblock) = STF_EVAP(:,:,iblock) + STF_PREC(:,:,iblock) & 
                                            + STF_ROFF(:,:,iblock)

!-----------------------------------------------
! compute Roce_{16O,18O} from SURF_VALS (delta units to R/Rstd):
!      to be scaled by Rstd in ocn_comp_mod.F90
!-----------------------------------------------
            ROCE_16O = 1.0_r8 
            ROCE_18O = 1._r8 + 1.e-3_r8 * SURF_VALS 

!-----------------------------------------------
! Now do for HDO tracer:
! for testing, because this is equal to the 18O tracer
! now due to the file forcing...
!-----------------------------------------------
           
            SURF_VALS = p5*(SURF_VALS_OLD(:,:,dD_ind,iblock) + &
                            SURF_VALS_CUR(:,:,dD_ind,iblock))
            dEVAP = (1-K) / (1 - RH(:,:,iblock)) * &
                    (alpha_wv(:,:,iblock) * (SURF_VALS + 1.e3) - &
                     RH(:,:,iblock) * (ATM_VALS(:,:,iblock) + 1.e3)) - 1.e3
            STF_EVAP(:,:,iblock) = EVAP_F(:,:,iblock)*fwmass_to_fwflux * dEVAP
            WISO_SFLUX_TAVG(:,:,6,iblock) = STF_EVAP(:,:,iblock)
            STF_ROFF(:,:,iblock) = ROFF_F(:,:,iblock)*fwmass_to_fwflux * d18OP(:,:,iblock)
            WISO_SFLUX_TAVG(:,:,8,iblock) = STF_ROFF(:,:,iblock)
            WISO_SFLUX_TAVG(:,:,24,iblock) = STF_ROFF(:,:,iblock)
            STF_MODULE(:,:,dD_ind,iblock) = STF_EVAP(:,:,iblock) + STF_PREC(:,:,iblock) &
                                          + STF_ROFF(:,:,iblock)

!-----------------------------------------------
! compute Roce_HDO from SURF_VALS (delta units to R/Rstd):
!      to be scaled by Rstd in ocn_comp_mod.F90
!-----------------------------------------------
            ROCE_HDO = 1.0_r8 + 1.0e-3_r8 * SURF_VALS 

!Still have the new iso precip/evap from x2o, though we don't use them here now:
            WISO_SFLUX_TAVG(:,:,9,iblock)  = PREC_16O_F
            WISO_SFLUX_TAVG(:,:,10,iblock) = PREC_18O_F
            WISO_SFLUX_TAVG(:,:,11,iblock) = PREC_HDO_F
            WISO_SFLUX_TAVG(:,:,12,iblock) = EVAP_16O_F
! note that EVAP_*_F comes from drver not scaled by Rstd_16O, so scale here 
            WISO_SFLUX_TAVG(:,:,13,iblock) = Rstd_16O*EVAP_18O_F
            WISO_SFLUX_TAVG(:,:,14,iblock) = Rstd_16O*EVAP_HDO_F
            WISO_SFLUX_TAVG(:,:,15,iblock) = ROCE_16O 
            WISO_SFLUX_TAVG(:,:,16,iblock) = ROCE_18O 
            WISO_SFLUX_TAVG(:,:,17,iblock) = ROCE_HDO 
!but will replace these with the ones computed from read-in files:
            WISO_SFLUX_TAVG(:,:,18,iblock) = STF_PREC_18O(:,:)
            WISO_SFLUX_TAVG(:,:,19,iblock) = STF_EVAP_18O(:,:)
            WISO_SFLUX_TAVG(:,:,21,iblock) = d18OP(:,:,iblock)
            WISO_SFLUX_TAVG(:,:,22,iblock) = STF_PREC_HDO(:,:)
            WISO_SFLUX_TAVG(:,:,23,iblock) = STF_EVAP_HDO(:,:)
         elsewhere
            STF_MODULE(:,:,d18o_ind,iblock) = c0
            STF_MODULE(:,:,dD_ind,iblock) = c0
            ROCE_16O(:,:) = c0
            ROCE_18O(:,:) = c0
            ROCE_HDO(:,:) = c0
         endwhere

! for sending up to coupler...
         call named_field_set(ROCE_16O_nf_ind, iblock, ROCE_16O)
         call named_field_set(ROCE_18O_nf_ind, iblock, ROCE_18O)
         call named_field_set(ROCE_HDO_nf_ind, iblock, ROCE_HDO)

   end do
 !$OMP END PARALLEL DO

   case ('model')

 !$OMP PARALLEL DO PRIVATE(iblock,SURF_VALS,                         &
 !$OMP        ROCE_16O,ROCE_18O,ROCE_HDO,                            &
 !$OMP        PREC_16O_F,PREC_18O_F,PREC_HDO_F,                      &
 !$OMP        EVAP_16O_F,EVAP_18O_F,EVAP_HDO_F,                      &
 !$OMP        MELT_16O_F,MELT_18O_F,MELT_HDO_F,                      &
 !$OMP        ROFF_16O_F,ROFF_18O_F,ROFF_HDO_F,                      &
 !$OMP        IOFF_16O_F,IOFF_18O_F,IOFF_HDO_F,                      &
 !$OMP        STF_PREC_18O,STF_PREC_HDO,STF_EVAP_18O,STF_EVAP_HDO,   &
 !$OMP        STF_MELT_18O,STF_MELT_HDO,STF_ROFF_18O,STF_ROFF_HDO,   &
 !$OMP        STF_IOFF_18O,STF_IOFF_HDO,STF_18O,STF_HDO)

   do iblock = 1, nblocks_clinic

!-----------------------------------------------------------------------
!  get fields from coupler
!-----------------------------------------------------------------------

         call named_field_get(PREC_16O_F_nf_ind, iblock, PREC_16O_F)
         call named_field_get(PREC_18O_F_nf_ind, iblock, PREC_18O_F)
         call named_field_get(PREC_HDO_F_nf_ind, iblock, PREC_HDO_F)
         call named_field_get(EVAP_16O_F_nf_ind, iblock, EVAP_16O_F)
         call named_field_get(EVAP_18O_F_nf_ind, iblock, EVAP_18O_F)
         call named_field_get(EVAP_HDO_F_nf_ind, iblock, EVAP_HDO_F)
         call named_field_get(MELT_16O_F_nf_ind, iblock, MELT_16O_F)
         call named_field_get(MELT_18O_F_nf_ind, iblock, MELT_18O_F)
         call named_field_get(MELT_HDO_F_nf_ind, iblock, MELT_HDO_F)
         call named_field_get(ROFF_16O_F_nf_ind, iblock, ROFF_16O_F)
         call named_field_get(ROFF_18O_F_nf_ind, iblock, ROFF_18O_F)
         call named_field_get(ROFF_HDO_F_nf_ind, iblock, ROFF_HDO_F)
         call named_field_get(IOFF_16O_F_nf_ind, iblock, IOFF_16O_F)
         call named_field_get(IOFF_18O_F_nf_ind, iblock, IOFF_18O_F)
         call named_field_get(IOFF_HDO_F_nf_ind, iblock, IOFF_HDO_F)

         where (LAND_MASK(:,:,iblock))

!-----------------------------------------------------------------------
!  calculate surface fluxes for delta tracers
!-----------------------------------------------------------------------

            STF_PREC_18O = (PREC_18O_F - PREC_16O_F) * fwmass_to_fwflux * c1000
            STF_EVAP_18O = (EVAP_18O_F - EVAP_16O_F) * fwmass_to_fwflux * c1000
            STF_MELT_18O = (MELT_18O_F - MELT_16O_F) * fwmass_to_fwflux * c1000
            STF_ROFF_18O = (ROFF_18O_F - ROFF_16O_F) * fwmass_to_fwflux * c1000
            STF_IOFF_18O = (IOFF_18O_F - IOFF_16O_F) * fwmass_to_fwflux * c1000

            STF_PREC_HDO = (PREC_HDO_F - PREC_16O_F) * fwmass_to_fwflux * c1000
            STF_EVAP_HDO = (EVAP_HDO_F - EVAP_16O_F) * fwmass_to_fwflux * c1000
            STF_MELT_HDO = (MELT_HDO_F - MELT_16O_F) * fwmass_to_fwflux * c1000
            STF_ROFF_HDO = (ROFF_HDO_F - ROFF_16O_F) * fwmass_to_fwflux * c1000
            STF_IOFF_HDO = (IOFF_HDO_F - IOFF_16O_F) * fwmass_to_fwflux * c1000

!-----------------------------------------------------------------------
!  sum up forcing terms
!-----------------------------------------------------------------------

            STF_18O = STF_PREC_18O + STF_EVAP_18O + STF_MELT_18O +  &
                      STF_ROFF_18O + STF_IOFF_18O +                 &
                      wiso_sflx_correction_18o * fwmass_to_fwflux * c1000
            STF_HDO = STF_PREC_HDO + STF_EVAP_HDO + STF_MELT_HDO +  &
                      STF_ROFF_HDO + STF_IOFF_HDO +                 & 
                      wiso_sflx_correction_hdo * fwmass_to_fwflux * c1000

!-----------------------------------------------------------------------
!  compute surface ratios: Roce_{16O,18O}
!-----------------------------------------------------------------------

            ROCE_16O = c1
            ROCE_18O = p001*p5*(SURF_VALS_OLD(:,:,d18o_ind,iblock) + &
                            SURF_VALS_CUR(:,:,d18o_ind,iblock)) + c1
            ROCE_HDO = p001*p5*(SURF_VALS_OLD(:,:,dD_ind,iblock) + &
                            SURF_VALS_CUR(:,:,dD_ind,iblock)) + c1

!-----------------------------------------------------------------------
!  history field updates:
!-----------------------------------------------------------------------

            WISO_SFLUX_TAVG(:,:, 4,iblock) = d18OP(:,:,iblock)

            WISO_SFLUX_TAVG(:,:, 9,iblock) = PREC_16O_F
            WISO_SFLUX_TAVG(:,:,10,iblock) = PREC_18O_F
            WISO_SFLUX_TAVG(:,:,11,iblock) = PREC_HDO_F
            WISO_SFLUX_TAVG(:,:,12,iblock) = EVAP_16O_F
            WISO_SFLUX_TAVG(:,:,13,iblock) = EVAP_18O_F*Rstd_16O
            WISO_SFLUX_TAVG(:,:,14,iblock) = EVAP_HDO_F*Rstd_16O
            WISO_SFLUX_TAVG(:,:,15,iblock) = MELT_16O_F 
            WISO_SFLUX_TAVG(:,:,16,iblock) = MELT_18O_F 
            WISO_SFLUX_TAVG(:,:,17,iblock) = MELT_HDO_F 
            WISO_SFLUX_TAVG(:,:,18,iblock) = ROFF_16O_F
            WISO_SFLUX_TAVG(:,:,19,iblock) = ROFF_18O_F
            WISO_SFLUX_TAVG(:,:,20,iblock) = ROFF_HDO_F
            WISO_SFLUX_TAVG(:,:,21,iblock) = IOFF_16O_F
            WISO_SFLUX_TAVG(:,:,22,iblock) = IOFF_18O_F
            WISO_SFLUX_TAVG(:,:,23,iblock) = IOFF_HDO_F
            WISO_SFLUX_TAVG(:,:,24,iblock) = ROCE_16O
            WISO_SFLUX_TAVG(:,:,25,iblock) = ROCE_18O 
            WISO_SFLUX_TAVG(:,:,26,iblock) = ROCE_HDO
            WISO_SFLUX_TAVG(:,:,27,iblock) = STF_PREC_18O
            WISO_SFLUX_TAVG(:,:,28,iblock) = STF_PREC_HDO
            WISO_SFLUX_TAVG(:,:,29,iblock) = STF_EVAP_18O
            WISO_SFLUX_TAVG(:,:,30,iblock) = STF_EVAP_HDO
            WISO_SFLUX_TAVG(:,:,31,iblock) = STF_MELT_18O
            WISO_SFLUX_TAVG(:,:,32,iblock) = STF_MELT_HDO
            WISO_SFLUX_TAVG(:,:,33,iblock) = STF_ROFF_18O
            WISO_SFLUX_TAVG(:,:,34,iblock) = STF_ROFF_HDO
            WISO_SFLUX_TAVG(:,:,35,iblock) = STF_IOFF_18O
            WISO_SFLUX_TAVG(:,:,36,iblock) = STF_IOFF_HDO
            WISO_SFLUX_TAVG(:,:,37,iblock) = STF_18O
            WISO_SFLUX_TAVG(:,:,38,iblock) = STF_HDO
         elsewhere
            ROCE_16O(:,:) = c0
            ROCE_18O(:,:) = c0
            ROCE_HDO(:,:) = c0
         endwhere

         call named_field_set(ROCE_16O_nf_ind, iblock, ROCE_16O)
         call named_field_set(ROCE_18O_nf_ind, iblock, ROCE_18O)
         call named_field_set(ROCE_HDO_nf_ind, iblock, ROCE_HDO)
   end do
 !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  adjust 18O and HDO surface flux
!-----------------------------------------------------------------------
   WISO_SFLUX_16O(:,:,:) = WISO_SFLUX_TAVG(:,:, 9,:) +  &
                           WISO_SFLUX_TAVG(:,:,12,:) +  &
                           WISO_SFLUX_TAVG(:,:,15,:) +  &
                           WISO_SFLUX_TAVG(:,:,18,:) +  &
                           WISO_SFLUX_TAVG(:,:,21,:)

   WISO_SFLUX_18O(:,:,:) = WISO_SFLUX_TAVG(:,:,10,:) +  &
                           WISO_SFLUX_TAVG(:,:,13,:) +  &
                           WISO_SFLUX_TAVG(:,:,16,:) +  &
                           WISO_SFLUX_TAVG(:,:,19,:) +  &
                           WISO_SFLUX_TAVG(:,:,22,:)

   WISO_SFLUX_HDO(:,:,:) = WISO_SFLUX_TAVG(:,:,11,:) +  &
                           WISO_SFLUX_TAVG(:,:,14,:) +  &
                           WISO_SFLUX_TAVG(:,:,17,:) +  &
                           WISO_SFLUX_TAVG(:,:,20,:) +  &
                           WISO_SFLUX_TAVG(:,:,23,:)

   call comp_tarea_avg(WISO_SFLUX_16O, avg_16O_F)
   call comp_tarea_avg(WISO_SFLUX_18O, avg_18O_F)
   call comp_tarea_avg(WISO_SFLUX_HDO, avg_HDO_F)

   WISO_SFLUX_18O(:,:,:) = WISO_SFLUX_18O(:,:,:) - avg_18O_F + avg_16O_F
   WISO_SFLUX_HDO(:,:,:) = WISO_SFLUX_HDO(:,:,:) - avg_HDO_F + avg_16O_F

   WISO_SFLUX_TAVG(:,:,37,:) = (WISO_SFLUX_18O - WISO_SFLUX_16O) * fwmass_to_fwflux * c1000
   WISO_SFLUX_TAVG(:,:,38,:) = (WISO_SFLUX_HDO - WISO_SFLUX_16O) * fwmass_to_fwflux * c1000

!-----------------------------------------------------------------------
!  do marginal sea balancing before setting the STF array
!-----------------------------------------------------------------------

   call ms_balancing_wiso (WISO_SFLUX_TAVG(:,:,37,:))
   call ms_balancing_wiso (WISO_SFLUX_TAVG(:,:,38,:))

   where (LAND_MASK(:,:,:))
      STF_MODULE(:,:,d18o_ind,:) = WISO_SFLUX_TAVG(:,:,37,:)
      STF_MODULE(:,:,dD_ind,:)   = WISO_SFLUX_TAVG(:,:,38,:)
   elsewhere
      STF_MODULE(:,:,d18o_ind,:) = c0
      STF_MODULE(:,:,dD_ind,:)   = c0
   endwhere

   case default
      call document(sub_name, 'wiso_formulation', wiso_formulation)

      call exit_POP(sigAbort, &
                    'wiso_set_sflux: Unknown value for wiso_formulation')

   end select

   call timer_stop(wiso_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine wiso_set_sflux

!***********************************************************************
!BOP
! !IROUTINE: wiso_tavg_forcing
! !INTERFACE:

 subroutine wiso_tavg_forcing

! !DESCRIPTION:
!  Make accumulation calls for forcing related tavg fields. This is
!  necessary because the forcing routines are called before tavg flags
!  are set.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1, nblocks_clinic
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 1,iblock),tavg_WISO_PREC,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 2,iblock),tavg_WISO_DAIR,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 3,iblock),tavg_WISO_RH,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 4,iblock),tavg_WISO_d18OR,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 5,iblock),tavg_Delta18O_EVAP,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 6,iblock),tavg_DeltaD_EVAP,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 7,iblock),tavg_Delta18O_ROFF,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 8,iblock),tavg_DeltaD_ROFF,iblock,1)

         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:, 9,iblock),tavg_PREC_16O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,10,iblock),tavg_PREC_18O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,11,iblock),tavg_PREC_HDO_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,12,iblock),tavg_EVAP_16O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,13,iblock),tavg_EVAP_18O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,14,iblock),tavg_EVAP_HDO_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,15,iblock),tavg_MELT_16O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,16,iblock),tavg_MELT_18O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,17,iblock),tavg_MELT_HDO_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,18,iblock),tavg_ROFF_16O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,19,iblock),tavg_ROFF_18O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,20,iblock),tavg_ROFF_HDO_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,21,iblock),tavg_IOFF_16O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,22,iblock),tavg_IOFF_18O_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,23,iblock),tavg_IOFF_HDO_F,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,24,iblock),tavg_Roce_16O,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,25,iblock),tavg_Roce_18O,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,26,iblock),tavg_Roce_HDO,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,27,iblock),tavg_STF_PREC_18O,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,28,iblock),tavg_STF_PREC_HDO,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,29,iblock),tavg_STF_EVAP_18O,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,30,iblock),tavg_STF_EVAP_HDO,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,31,iblock),tavg_STF_MELT_18O,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,32,iblock),tavg_STF_MELT_HDO,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,33,iblock),tavg_STF_ROFF_18O,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,34,iblock),tavg_STF_ROFF_HDO,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,35,iblock),tavg_STF_IOFF_18O,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,36,iblock),tavg_STF_IOFF_HDO,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,37,iblock),tavg_STF_18O,iblock,1)
         call accumulate_tavg_field(WISO_SFLUX_TAVG(:,:,38,iblock),tavg_STF_HDO,iblock,1)
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine wiso_tavg_forcing

!***********************************************************************
!BOP
! !IROUTINE: comp_surf_avg
! !INTERFACE:

 subroutine comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR)

! !DESCRIPTION:
!  compute average surface tracer values
!
!  avg = sum(SURF_VAL*TAREA) / sum(TAREA)
!  with the sum taken over ocean points only
!
! !REVISION HISTORY:
!  by Jiaxu Zhang 1/6/2012

   use grid, only: TAREA, RCALCT, area_t
   use global_reductions, only: global_sum

! !INPUT PARAMETERS:

  real (r8), dimension(nx_block,ny_block,wiso_tracer_cnt,max_blocks_clinic), &
      intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,       & ! tracer index
      iblock,  & ! block index
      dcount,  & ! diag counter
      ib,ie,jb,je

   real (r8), dimension(max_blocks_clinic,wiso_tracer_cnt) :: &
      local_sums ! array for holding block sums of each diagnostic

   real (r8) :: &
      sum_tmp ! temp for local sum

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1, &! local work space
      TFACT   ! factor for normalizing sums

   type (block) :: &
      this_block ! block information for current block

!-----------------------------------------------------------------------

   local_sums = c0

!jw   !$OMP PARALLEL DO PRIVATE(iblock,this_block,ib,ie,jb,je,TFACT,n,WORK1)
   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      ib = this_block%ib
      ie = this_block%ie
      jb = this_block%jb
      je = this_block%je
      TFACT = TAREA(:,:,iblock)*RCALCT(:,:,iblock)

      do n = 1, wiso_tracer_cnt
         if (vflux_flag(n)) then
            WORK1 = p5*(SURF_VALS_OLD(:,:,n,iblock) + &
                        SURF_VALS_CUR(:,:,n,iblock))*TFACT
            local_sums(iblock,n) = sum(WORK1(ib:ie,jb:je))
         endif
      end do
   end do
!jw   !$OMP END PARALLEL DO

   do n = 1, wiso_tracer_cnt
      if (vflux_flag(n)) then
         sum_tmp = sum(local_sums(:,n))
         surf_avg(n) = global_sum(sum_tmp,distrb_clinic)/area_t
      endif
   end do

   if(my_task == master_task) then
      write(stdout,*)' Calculating surface tracer averages'
      do n = 1, wiso_tracer_cnt
         if (vflux_flag(n)) then
            write(stdout,*) n, surf_avg(n)
         endif
      end do
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_surf_avg

!***********************************************************************
!BOP
! !IROUTINE: comp_tarea_avg
! !INTERFACE:

 subroutine comp_tarea_avg(XY, XYave)

! !DESCRIPTION:
!  compute average surface everage
!
!  XYave = sum(XY*TAREA) / sum(TAREA)
!  with the sum taken over ocean points only
!
! !REVISION HISTORY:
!  by Jiang Zhu, Sep. 7, 2015

   use grid, only: TAREA, RCALCT, area_t
   use global_reductions, only: global_sum

! !INPUT PARAMETERS:

  real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(in) :: XY
  real (r8), intent(out)  :: &
      XYave ! resulting global ave

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,       & ! tracer index
      iblock,  & ! block index
      dcount,  & ! diag counter
      ib,ie,jb,je

   real (r8), dimension(max_blocks_clinic) :: &
      local_sums ! array for holding block sums of each diagnostic

   real (r8) :: &
      sum_tmp ! temp for local sum

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1, &! local work space
      TFACT   ! factor for normalizing sums

   type (block) :: &
      this_block ! block information for current block

!-----------------------------------------------------------------------

   local_sums = c0

!jw   !$OMP PARALLEL DO PRIVATE(iblock,this_block,ib,ie,jb,je,TFACT,n,WORK1)
   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      ib = this_block%ib
      ie = this_block%ie
      jb = this_block%jb
      je = this_block%je
      TFACT = TAREA(:,:,iblock)*RCALCT(:,:,iblock)

      WORK1 = XY(:,:,iblock) * TFACT
      local_sums(iblock) = sum(WORK1(ib:ie,jb:je))
   end do
!jw   !$OMP END PARALLEL DO

   sum_tmp = sum(local_sums)
   XYave   = global_sum(sum_tmp,distrb_clinic) / area_t

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_tarea_avg

!***********************************************************************
!BOP
! !IROUTINE: wiso_write_restart
! !INTERFACE:

 subroutine wiso_write_restart(restart_file, action)

! !DESCRIPTION:
!  write auxiliary fields & scalars to restart files
!
! !REVISION HISTORY:
!  same as module

   use constants, only: char_blank

! !INPUT PARAMETERS:

   character(*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character (char_len) :: &
      short_name  ! tracer name temporaries

   integer (int_kind) :: n

!-----------------------------------------------------------------------

   if (trim(action) == 'add_attrib_file') then
      short_name = char_blank
      do n=1,wiso_tracer_cnt
         if (vflux_flag(n)) then
            short_name = 'surf_avg_' /&
                      &/ ind_name_table(n)%name
            call add_attrib_file(restart_file,trim(short_name),surf_avg(n))
         endif
      end do
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine wiso_write_restart

!***********************************************************************
!BOP
! !IROUTINE: wiso_tracer_ref_val
! !INTERFACE:

 function wiso_tracer_ref_val(ind)

! !DESCRIPTION:
!  return reference value for tracers using virtual fluxes
!
! !REVISION HISTORY:
!  by Jiaxu Zhang 1/6/2012

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: ind

! !OUTPUT PARAMETERS:

   real (r8) :: wiso_tracer_ref_val

!EOP
!BOC
!-----------------------------------------------------------------------

   if (vflux_flag(ind)) then
      wiso_tracer_ref_val = surf_avg(ind)
   else
      wiso_tracer_ref_val = c0
   endif

!-----------------------------------------------------------------------
!EOC

 end function wiso_tracer_ref_val

!***********************************************************************

end module wiso_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
