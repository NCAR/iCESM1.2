!===============================================================================
! SVN $Id: seq_diag_mct.F90 54832 2013-11-01 20:41:44Z brady@ucar.edu $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/branch_tags/geotrace_tags/geotrace_atm_pop_drvseq4_2_33/driver/seq_diag_mct.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_diag_wiso_mod -- computes spatial \& time averages of fluxed quatities
!
! !DESCRIPTION:
!    The coupler is required to do certain diagnostics, those calculations are
!    located in this module.
!
! !REMARKS:
!    CESM sign convention for fluxes is positive downward with hierarchy being
!       atm/glc/lnd/rof/ice/ocn
!    Sign convention:
!       positive value <=> the model is gaining water, heat, momentum, etc.
!    Unit convention:
!       heat flux     ~ W/m^2
!       momentum flux ~ N/m^2
!       water flux    ~ (kg/s)/m^2
!       salt  flux    ~ (kg/s)/m^2
!
! !REVISION HISTORY:
!    2014-sep-04 - J. Zhu      - modified for water isotopic fluxes
!    2012-aug-20 - T. Craig    - add rof component
!    2008-jul-10 - T. Craig    - updated budget implementation
!    2007-may-07 - B. Kauffman - initial port to cpl7. 
!    2002-nov-21 - R. Jacob    - initial port to cpl6. 
!    199x-mmm-dd - B. Kauffman - original version in cpl4.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_diag_wiso_mct
  
! !USES:

   use shr_kind_mod, only: r8 => shr_kind_r8, in=>shr_kind_in
   use shr_kind_mod, only: i8 => shr_kind_i8,  cl=>shr_kind_cl
   use shr_sys_mod       ! system calls
   use shr_mpi_mod       ! mpi wrappers
   use shr_const_mod     ! shared constants
   use mct_mod           ! mct wrappers
   use esmf

   use seq_comm_mct  ! mpi comm groups & related
   use seq_cdata_mod
   use seq_timemgr_mod

   implicit none
   save
   private

! !PUBLIC TYPES:

   ! none

!PUBLIC MEMBER FUNCTIONS:

   public seq_diag_zero_wiso_mct

   public seq_diag_atm_wiso_mct
   public seq_diag_lnd_wiso_mct
   public seq_diag_rtm_wiso_mct
   public seq_diag_ocn_wiso_mct
   public seq_diag_ice_wiso_mct
   public seq_diag_accum_wiso_mct
   public seq_diag_sum0_wiso_mct
   public seq_diag_print_wiso_mct

!EOP

   !----------------------------------------------------------------------------
   ! Local data
   !----------------------------------------------------------------------------

   !----- local constants -----
   real(r8),parameter :: HFLXtoWFLX = & ! water flux implied by latent heat of fusion
   &  - (shr_const_ocn_ref_sal-shr_const_ice_ref_sal) / &
   &    (shr_const_ocn_ref_sal*shr_const_latice)

   !--- C for component ---
   !--- "r" is recieve in the coupler, "s" is send from the coupler

   integer(in),parameter :: c_size = 20

   integer(in),parameter :: c_atm_as   = 1 ! model index: atm
   integer(in),parameter :: c_atm_ar   = 2 ! model index: atm
   integer(in),parameter :: c_inh_is   = 3 ! model index: ice, northern
   integer(in),parameter :: c_inh_ir   = 4 ! model index: ice, northern
   integer(in),parameter :: c_ish_is   = 5 ! model index: ice, southern
   integer(in),parameter :: c_ish_ir   = 6 ! model index: ice, southern
   integer(in),parameter :: c_lnd_ls   = 7 ! model index: lnd
   integer(in),parameter :: c_lnd_lr   = 8 ! model index: lnd
   integer(in),parameter :: c_ocn_os   = 9 ! model index: ocn
   integer(in),parameter :: c_ocn_or   =10 ! model index: ocn
   integer(in),parameter :: c_rof_rs   =11 ! model index: lnd
   integer(in),parameter :: c_rof_rr   =12 ! model index: lnd
   ! --- on atm grid ---
   integer(in),parameter :: c_inh_as   =13 ! model index: ice, northern
   integer(in),parameter :: c_inh_ar   =14 ! model index: ice, northern
   integer(in),parameter :: c_ish_as   =15 ! model index: ice, southern
   integer(in),parameter :: c_ish_ar   =16 ! model index: ice, southern
   integer(in),parameter :: c_lnd_as   =17 ! model index: lnd
   integer(in),parameter :: c_lnd_ar   =18 ! model index: lnd
   integer(in),parameter :: c_ocn_as   =19 ! model index: ocn
   integer(in),parameter :: c_ocn_ar   =20 ! model index: ocn

   character(len=8),parameter :: cname(c_size) = &
      (/' c2a_atm',' a2c_atm',' c2i_inh',' i2c_inh',' c2i_ish',' i2c_ish', &
        ' c2l_lnd',' l2c_lnd',' c2o_ocn',' o2c_ocn',' c2r_rof',' r2c_rof', &
        ' c2a_inh',' a2c_inh',' c2a_ish',' a2c_ish', &
        ' c2a_lnd',' a2c_lnd',' c2a_ocn',' a2c_ocn' /)

   !--- F for field ---

   integer(in),parameter :: f_size      = 21
   integer(in),parameter :: f_18O       = 8
   integer(in),parameter :: f_HDO       = 15

   integer(in),parameter :: f_wfrz_16O  = 1     ! water: freezing
   integer(in),parameter :: f_wmelt_16O = 2     ! water: melting
   integer(in),parameter :: f_wrain_16O = 3     ! water: precip, liquid
   integer(in),parameter :: f_wsnow_16O = 4     ! water: precip, frozen
   integer(in),parameter :: f_wevap_16O = 5     ! water: evaporation
   integer(in),parameter :: f_wroff_16O = 6     ! water: runoff/flood
   integer(in),parameter :: f_wioff_16O = 7     ! water: frozen runoff
   integer(in),parameter :: f_wfrz_18O  = 8     ! water: freezing
   integer(in),parameter :: f_wmelt_18O = 9     ! water: melting
   integer(in),parameter :: f_wrain_18O = 10    ! water: precip, liquid
   integer(in),parameter :: f_wsnow_18O = 11    ! water: precip, frozen
   integer(in),parameter :: f_wevap_18O = 12    ! water: evaporation
   integer(in),parameter :: f_wroff_18O = 13    ! water: runoff/flood
   integer(in),parameter :: f_wioff_18O = 14    ! water: frozen runoff
   integer(in),parameter :: f_wfrz_HDO  = 15    ! water: freezing
   integer(in),parameter :: f_wmelt_HDO = 16    ! water: melting
   integer(in),parameter :: f_wrain_HDO = 17    ! water: precip, liquid
   integer(in),parameter :: f_wsnow_HDO = 18    ! water: precip, frozen
   integer(in),parameter :: f_wevap_HDO = 19    ! water: evaporation
   integer(in),parameter :: f_wroff_HDO = 20    ! water: runoff/flood
   integer(in),parameter :: f_wioff_HDO = 21    ! water: frozen runoff

   character(len=12),parameter :: fname(f_size) = &
      (/' wfreeze_16O','   wmelt_16O','   wrain_16O','   wsnow_16O', &
        '   wevap_16O',' wrunoff_16O',' wfrzrof_16O', &
        ' wfreeze_18O','   wmelt_18O','   wrain_18O','   wsnow_18O', &
        '   wevap_18O',' wrunoff_18O',' wfrzrof_18O', &
        ' wfreeze_HDO','   wmelt_HDO','   wrain_HDO','   wsnow_HDO', &
        '   wevap_HDO',' wrunoff_HDO',' wfrzrof_HDO'/)

   !--- P for period --- 
   integer(in),parameter :: p_size = 5

   integer(in),parameter :: p_inst = 1
   integer(in),parameter :: p_day  = 2
   integer(in),parameter :: p_mon  = 3
   integer(in),parameter :: p_ann  = 4
   integer(in),parameter :: p_inf  = 5

   character(len=8),parameter :: pname(p_size) = &
      (/'    inst','   daily',' monthly','  annual','all_time' /)

! !PUBLIC DATA MEMBERS

   !--- time-averaged (annual?) global budge diagnostics ---
   !--- note: call sum0 then save budg_dataG and budg_ns on restart from/to root pe ---
   real(r8),public :: budg_dataL_wiso(f_size,c_size,p_size) ! local sum, valid on all pes
   real(r8),public :: budg_dataG_wiso(f_size,c_size,p_size) ! global sum, valid only on root pe
   real(r8),public :: budg_ns_wiso   (f_size,c_size,p_size) ! counter, valid only on root pe

   character(len=*),parameter :: afldname  = 'aream'
   character(len=*),parameter :: latname   = 'lat'
   character(len=*),parameter :: afracname = 'afrac'
   character(len=*),parameter :: lfracname = 'lfrac'
   character(len=*),parameter :: ofracname = 'ofrac'
   character(len=*),parameter :: ifracname = 'ifrac'

   character(*),parameter :: modName = "(seq_diag_wiso_mct) "

   integer(in),parameter :: debug = 0 ! internal debug level

! !PRIVATE DATA MEMBERS

   integer :: index_a2x_Faxa_rainc_16O
   integer :: index_a2x_Faxa_rainc_18O
   integer :: index_a2x_Faxa_rainc_HDO
   integer :: index_a2x_Faxa_rainl_16O
   integer :: index_a2x_Faxa_rainl_18O
   integer :: index_a2x_Faxa_rainl_HDO
   integer :: index_a2x_Faxa_snowc_16O
   integer :: index_a2x_Faxa_snowc_18O
   integer :: index_a2x_Faxa_snowc_HDO
   integer :: index_a2x_Faxa_snowl_16O
   integer :: index_a2x_Faxa_snowl_18O
   integer :: index_a2x_Faxa_snowl_HDO

   integer :: index_x2a_Faxx_evap_16O
   integer :: index_x2a_Faxx_evap_18O
   integer :: index_x2a_Faxx_evap_HDO

   integer :: index_l2x_Fall_evap_16O
   integer :: index_l2x_Fall_evap_18O
   integer :: index_l2x_Fall_evap_HDO

   integer :: index_l2x_Flrl_rofliq_16O
   integer :: index_l2x_Flrl_rofliq_18O
   integer :: index_l2x_Flrl_rofliq_HDO
   integer :: index_l2x_Flrl_rofice_16O
   integer :: index_l2x_Flrl_rofice_18O
   integer :: index_l2x_Flrl_rofice_HDO

   integer :: index_x2l_Faxa_rainc_16O
   integer :: index_x2l_Faxa_rainc_18O
   integer :: index_x2l_Faxa_rainc_HDO
   integer :: index_x2l_Faxa_rainl_16O
   integer :: index_x2l_Faxa_rainl_18O
   integer :: index_x2l_Faxa_rainl_HDO
   integer :: index_x2l_Faxa_snowc_16O
   integer :: index_x2l_Faxa_snowc_18O
   integer :: index_x2l_Faxa_snowc_HDO
   integer :: index_x2l_Faxa_snowl_16O
   integer :: index_x2l_Faxa_snowl_18O
   integer :: index_x2l_Faxa_snowl_HDO
   integer :: index_x2l_Flrr_flood_16O
   integer :: index_x2l_Flrr_flood_18O
   integer :: index_x2l_Flrr_flood_HDO

   integer :: index_r2x_Forr_roff_16O
   integer :: index_r2x_Forr_roff_18O
   integer :: index_r2x_Forr_roff_HDO
   integer :: index_r2x_Forr_ioff_16O
   integer :: index_r2x_Forr_ioff_18O
   integer :: index_r2x_Forr_ioff_HDO
   integer :: index_r2x_Flrr_flood_16O
   integer :: index_r2x_Flrr_flood_18O
   integer :: index_r2x_Flrr_flood_HDO

   integer :: index_x2r_Flrl_rofliq_16O
   integer :: index_x2r_Flrl_rofliq_18O
   integer :: index_x2r_Flrl_rofliq_HDO
   integer :: index_x2r_Flrl_rofice_16O
   integer :: index_x2r_Flrl_rofice_18O
   integer :: index_x2r_Flrl_rofice_HDO

   integer :: index_xao_Faox_evap_16O
   integer :: index_xao_Faox_evap_18O
   integer :: index_xao_Faox_evap_HDO

   integer :: index_x2o_Fioi_meltw_16O
   integer :: index_x2o_Fioi_meltw_18O
   integer :: index_x2o_Fioi_meltw_HDO
   integer :: index_x2o_Faxa_rain_16O
   integer :: index_x2o_Faxa_rain_18O
   integer :: index_x2o_Faxa_rain_HDO
   integer :: index_x2o_Faxa_snow_16O
   integer :: index_x2o_Faxa_snow_18O
   integer :: index_x2o_Faxa_snow_HDO

   integer :: index_i2x_Fioi_meltw_16O
   integer :: index_i2x_Fioi_meltw_18O
   integer :: index_i2x_Fioi_meltw_HDO
   integer :: index_i2x_Faii_evap_16O
   integer :: index_i2x_Faii_evap_18O
   integer :: index_i2x_Faii_evap_HDO

   integer :: index_x2i_Faxa_rain_16O
   integer :: index_x2i_Faxa_rain_18O
   integer :: index_x2i_Faxa_rain_HDO
   integer :: index_x2i_Faxa_snow_16O
   integer :: index_x2i_Faxa_snow_18O
   integer :: index_x2i_Faxa_snow_HDO

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_zero_mct - zero out global budget diagnostic data.
!
! !DESCRIPTION:
!    Zero out global budget diagnostic data.
!
! !REVISION HISTORY:
!    2008-jul-11 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_zero_wiso_mct(EClock,mode)

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock), intent(in),optional :: EClock
   character(len=*), intent(in),optional :: mode

!EOP

   integer(IN) :: ip,yr,mon,day,sec
   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_zero_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. present(EClock) .and. .not. present(mode)) then
      call shr_sys_abort(subName//' ERROR EClock or mode should be present')
   endif

   if (present(EClock)) then
      call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
         curr_mon=mon,curr_day=day,curr_tod=sec)

      do ip = 1,p_size
         if (ip == p_inst) then
            budg_dataL_wiso(:,:,ip) = 0.0_r8
            budg_dataG_wiso(:,:,ip) = 0.0_r8
            budg_ns_wiso(:,:,ip)    = 0.0_r8
         endif
         if (ip==p_day .and. sec==0) then
            budg_dataL_wiso(:,:,ip) = 0.0_r8
            budg_dataG_wiso(:,:,ip) = 0.0_r8
            budg_ns_wiso(:,:,ip)    = 0.0_r8
         endif
         if (ip==p_mon .and. day==1 .and. sec==0) then
            budg_dataL_wiso(:,:,ip) = 0.0_r8
            budg_dataG_wiso(:,:,ip) = 0.0_r8
            budg_ns_wiso(:,:,ip)    = 0.0_r8
         endif
         if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
            budg_dataL_wiso(:,:,ip) = 0.0_r8
            budg_dataG_wiso(:,:,ip) = 0.0_r8
            budg_ns_wiso(:,:,ip)    = 0.0_r8
         endif
      enddo
   endif

   if (present(mode)) then
      if (trim(mode) == 'inst') then
         budg_dataL_wiso(:,:,p_inst)= 0.0_r8
         budg_dataG_wiso(:,:,p_inst)= 0.0_r8
         budg_ns_wiso(:,:,p_inst)   = 0.0_r8
      elseif (trim(mode) == 'day') then
         budg_dataL_wiso(:,:,p_day) = 0.0_r8
         budg_dataG_wiso(:,:,p_day) = 0.0_r8
         budg_ns_wiso(:,:,p_day)    = 0.0_r8
      elseif (trim(mode) == 'mon') then
         budg_dataL_wiso(:,:,p_mon) = 0.0_r8
         budg_dataG_wiso(:,:,p_mon) = 0.0_r8
         budg_ns_wiso(:,:,p_mon)    = 0.0_r8
      elseif (trim(mode) == 'ann') then
         budg_dataL_wiso(:,:,p_ann) = 0.0_r8
         budg_dataG_wiso(:,:,p_ann) = 0.0_r8
         budg_ns_wiso(:,:,p_ann)    = 0.0_r8
      elseif (trim(mode) == 'inf') then
         budg_dataL_wiso(:,:,p_inf) = 0.0_r8
         budg_dataG_wiso(:,:,p_inf) = 0.0_r8
         budg_ns_wiso(:,:,p_inf)    = 0.0_r8
      elseif (trim(mode) == 'all') then
         budg_dataL_wiso(:,:,:)     = 0.0_r8
         budg_dataG_wiso(:,:,:)     = 0.0_r8
         budg_ns_wiso(:,:,:)         = 0.0_r8
      else
         call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
      endif
   endif

end subroutine seq_diag_zero_wiso_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_accum_mct - accum out global budget diagnostic data.
!
! !DESCRIPTION:
!    Accum out global budget diagnostic data.
!
! !REVISION HISTORY:
!    2008-jul-11 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_accum_wiso_mct()

! !INPUT/OUTPUT PARAMETERS:

!EOP

   integer(in) :: ip

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_accum_wiso_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   do ip = p_inst+1,p_size
      budg_dataL_wiso(:,:,ip) = budg_dataL_wiso(:,:,ip) + budg_dataL_wiso(:,:,p_inst)
   enddo
   budg_ns_wiso(:,:,:) = budg_ns_wiso(:,:,:) + 1.0_r8

end subroutine seq_diag_accum_wiso_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_sum0_mct - sum local to global on root
!
! !DESCRIPTION:
!    Sum local values to global on root
!
! !REVISION HISTORY:
!    2008-jul-19 - T. Craig - update
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_sum0_wiso_mct()

! !INPUT/OUTPUT PARAMETERS:

!EOP

   real(r8) :: budg_dataGtmp(f_size,c_size,p_size) ! temporary sum
   integer(in)      :: mpicom      ! mpi comm
   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_sum0_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call seq_comm_setptrs(CPLID,mpicom=mpicom)
   budg_dataGtmp = 0.0_r8
   call shr_mpi_sum(budg_dataL_wiso,budg_dataGtmp,mpicom,subName)
   budg_dataG_wiso = budg_dataG_wiso + budg_dataGtmp
   budg_dataL_wiso = 0.0_r8

end subroutine seq_diag_sum0_wiso_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_atm_wiso_mct - compute global atm input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global atm input/output flux diagnostics for water isotopes
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!    2014-jul-23 - J. Zhu - modified
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_atm_wiso_mct( dom_a, frac_a, a2x_a, x2a_a )

! !INPUT/OUTPUT PARAMETERS:

   type(mct_gGrid),intent(in)          ::  dom_a ! model domain
   type(mct_aVect),intent(in)          :: frac_a ! domain fractions
   type(mct_aVect),intent(in),optional ::  a2x_a ! model to drv bundle
   type(mct_aVect),intent(in),optional ::  x2a_a ! drv to model bundle

!EOP

   !----- local -----
   integer(in)      :: k,n,ic,if,ip      ! generic index
   integer(in)      :: kArea             ! index of area field in aVect
   integer(in)      :: kLat              ! index of lat field in aVect
   integer(in)      :: kl,ka,ko,ki       ! fraction indices
   integer(in)      :: lSize             ! size of aVect
   real(r8)         :: da,di,do,dl       ! area of a grid cell
   logical,save     :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_atm_wiso_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. present(a2x_a) .and. .not. present(x2a_a)) then
      call shr_sys_abort(subName//"ERROR: must input a bundle")
   end if

   kArea = mct_aVect_indexRA(dom_a%data,afldname)
   kLat  = mct_aVect_indexRA(dom_a%data,latname)
   ka    = mct_aVect_indexRA(frac_a,afracname)
   kl    = mct_aVect_indexRA(frac_a,lfracname)
   ko    = mct_aVect_indexRA(frac_a,ofracname)
   ki    = mct_aVect_indexRA(frac_a,ifracname)

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   ip = p_inst

   if (present(a2x_a)) then
      if (first_time) then
         index_a2x_Faxa_rainc_16O   = mct_aVect_indexRA(a2x_a,'Faxa_rainc_16O')
         index_a2x_Faxa_rainc_18O   = mct_aVect_indexRA(a2x_a,'Faxa_rainc_18O')
         index_a2x_Faxa_rainc_HDO   = mct_aVect_indexRA(a2x_a,'Faxa_rainc_HDO')
         index_a2x_Faxa_rainl_16O   = mct_aVect_indexRA(a2x_a,'Faxa_rainl_16O')
         index_a2x_Faxa_rainl_18O   = mct_aVect_indexRA(a2x_a,'Faxa_rainl_18O')
         index_a2x_Faxa_rainl_HDO   = mct_aVect_indexRA(a2x_a,'Faxa_rainl_HDO')
         index_a2x_Faxa_snowc_16O   = mct_aVect_indexRA(a2x_a,'Faxa_snowc_16O')
         index_a2x_Faxa_snowc_18O   = mct_aVect_indexRA(a2x_a,'Faxa_snowc_18O')
         index_a2x_Faxa_snowc_HDO   = mct_aVect_indexRA(a2x_a,'Faxa_snowc_HDO')
         index_a2x_Faxa_snowl_16O   = mct_aVect_indexRA(a2x_a,'Faxa_snowl_16O')
         index_a2x_Faxa_snowl_18O   = mct_aVect_indexRA(a2x_a,'Faxa_snowl_18O')
         index_a2x_Faxa_snowl_HDO   = mct_aVect_indexRA(a2x_a,'Faxa_snowl_HDO')
      end if

      lSize = mct_avect_lSize(a2x_a)
      do n=1,lSize
      do k=1,4

         if (k == 1) then
            ic = c_atm_ar
            da = -dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ka,n)
         elseif (k == 2) then
            ic = c_lnd_ar
            da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(kl,n)
         elseif (k == 3) then
            ic = c_ocn_ar
            da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ko,n)
         elseif (k == 4) then
            if (dom_a%data%rAttr(kLat,n) > 0.0_r8) then
               ic = c_inh_ar
            else
               ic = c_ish_ar
            endif
            da = dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ki,n)
         endif

         if = f_wrain_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_rainc_16O,n) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_rainl_16O,n)
         if = f_wrain_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_rainc_18O,n) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_rainl_18O,n)
         if = f_wrain_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_rainc_HDO,n) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_rainl_HDO,n)
         if = f_wsnow_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_snowc_16O,n) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_snowl_16O,n)
         if = f_wsnow_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_snowc_18O,n) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_snowl_18O,n)
         if = f_wsnow_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_snowc_HDO,n) + &
                                     da*a2x_a%rAttr(index_a2x_Faxa_snowl_HDO,n)
      enddo
      enddo
   end if

   if (present(x2a_a)) then
      if (first_time) then
         index_x2a_Faxx_evap_16O = mct_aVect_indexRA(x2a_a,'Faxx_evap_16O')
         index_x2a_Faxx_evap_18O = mct_aVect_indexRA(x2a_a,'Faxx_evap_18O')
         index_x2a_Faxx_evap_HDO = mct_aVect_indexRA(x2a_a,'Faxx_evap_HDO')
      end if

      lSize = mct_avect_lSize(x2a_a)
      do n=1,lSize
      do k=1,4

         if (k == 1) then
            ic = c_atm_as
            da = -dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ka,n)
         elseif (k == 2) then
            ic = c_lnd_as
            da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(kl,n)
         elseif (k == 3) then
            ic = c_ocn_as
            da =  dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ko,n)
         elseif (k == 4) then
            if (dom_a%data%rAttr(kLat,n) > 0.0_r8) then
               ic = c_inh_as
            else
               ic = c_ish_as
            endif
            da = dom_a%data%rAttr(kArea,n) * frac_a%rAttr(ki,n)
         endif

         if = f_wevap_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*x2a_a%rAttr(index_x2a_Faxx_evap_16O,n)
         if = f_wevap_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*x2a_a%rAttr(index_x2a_Faxx_evap_18O,n)
         if = f_wevap_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     da*x2a_a%rAttr(index_x2a_Faxx_evap_HDO,n)

      enddo
      enddo
   end if

   first_time = .false.

end subroutine seq_diag_atm_wiso_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_lnd_wiso_mct - compute global lnd input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global lnd input/output flux diagnostics for water isotopes
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!    2014-jul-23 - J. Zhu - modified
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_lnd_wiso_mct( dom_l, frac_l, l2x_l, x2l_l)

   type(mct_gGrid),intent(in)          ::  dom_l ! model domain
   type(mct_aVect),intent(in)          :: frac_l ! frac bundle
   type(mct_aVect),intent(in),optional ::  l2x_l ! model to drv bundle
   type(mct_aVect),intent(in),optional ::  x2l_l ! drv to model bundle

!EOP

   !----- local -----
   integer(in)      :: k,n,ic,if,ip      ! generic index
   integer(in)      :: kArea             ! index of area field in aVect
   integer(in)      :: kLat              ! index of lat field in aVect
   integer(in)      :: kl,ka,ko,ki       ! fraction indices
   integer(in)      :: lSize             ! size of aVect
   real(r8)         :: da,di,do,dl       ! area of a grid cell
   logical,save     :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_lnd_wiso_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. present(l2x_l) .and. .not. present(x2l_l)) then
      call shr_sys_abort(subName//"ERROR: must input a bundle")
   end if

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   ip = p_inst

   kArea = mct_aVect_indexRA(dom_l%data,afldname)
   kl    = mct_aVect_indexRA(frac_l,lfracname)

   if (present(l2x_l)) then
      if (first_time) then
         index_l2x_Fall_evap_16O    = mct_aVect_indexRA(l2x_l,'Fall_evap_16O')
         index_l2x_Fall_evap_18O    = mct_aVect_indexRA(l2x_l,'Fall_evap_18O')
         index_l2x_Fall_evap_HDO    = mct_aVect_indexRA(l2x_l,'Fall_evap_HDO')
         index_l2x_Flrl_rofliq_16O  = mct_aVect_indexRA(l2x_l,'Flrl_rofliq_16O')
         index_l2x_Flrl_rofliq_18O  = mct_aVect_indexRA(l2x_l,'Flrl_rofliq_18O')
         index_l2x_Flrl_rofliq_HDO  = mct_aVect_indexRA(l2x_l,'Flrl_rofliq_HDO')
         index_l2x_Flrl_rofice_16O  = mct_aVect_indexRA(l2x_l,'Flrl_rofice_16O')
         index_l2x_Flrl_rofice_18O  = mct_aVect_indexRA(l2x_l,'Flrl_rofice_18O')
         index_l2x_Flrl_rofice_HDO  = mct_aVect_indexRA(l2x_l,'Flrl_rofice_HDO')
      end if

      lSize = mct_avect_lSize(l2x_l)
      ic = c_lnd_lr
      do n=1,lSize
         dl =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
         if = f_wevap_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*l2x_l%rAttr(index_l2x_Fall_evap_16O,n)
         if = f_wevap_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*l2x_l%rAttr(index_l2x_Fall_evap_18O,n)
         if = f_wevap_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*l2x_l%rAttr(index_l2x_Fall_evap_HDO,n)

         if = f_wroff_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*l2x_l%rAttr(index_l2x_Flrl_rofliq_16O,n)
         if = f_wroff_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*l2x_l%rAttr(index_l2x_Flrl_rofliq_18O,n)
         if = f_wroff_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*l2x_l%rAttr(index_l2x_Flrl_rofliq_HDO,n)

         if = f_wioff_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*l2x_l%rAttr(index_l2x_Flrl_rofice_16O,n)
         if = f_wioff_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*l2x_l%rAttr(index_l2x_Flrl_rofice_18O,n)
         if = f_wioff_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*l2x_l%rAttr(index_l2x_Flrl_rofice_HDO,n)
      end do
   end if

   if (present(x2l_l)) then
      if (first_time) then
         index_x2l_Faxa_rainc_16O = mct_aVect_indexRA(x2l_l,'Faxa_rainc_16O')
         index_x2l_Faxa_rainc_18O = mct_aVect_indexRA(x2l_l,'Faxa_rainc_18O')
         index_x2l_Faxa_rainc_HDO = mct_aVect_indexRA(x2l_l,'Faxa_rainc_HDO')
         index_x2l_Faxa_rainl_16O = mct_aVect_indexRA(x2l_l,'Faxa_rainl_16O')
         index_x2l_Faxa_rainl_18O = mct_aVect_indexRA(x2l_l,'Faxa_rainl_18O')
         index_x2l_Faxa_rainl_HDO = mct_aVect_indexRA(x2l_l,'Faxa_rainl_HDO')
         index_x2l_Faxa_snowc_16O = mct_aVect_indexRA(x2l_l,'Faxa_snowc_16O')
         index_x2l_Faxa_snowc_18O = mct_aVect_indexRA(x2l_l,'Faxa_snowc_18O')
         index_x2l_Faxa_snowc_HDO = mct_aVect_indexRA(x2l_l,'Faxa_snowc_HDO')
         index_x2l_Faxa_snowl_16O = mct_aVect_indexRA(x2l_l,'Faxa_snowl_16O')
         index_x2l_Faxa_snowl_18O = mct_aVect_indexRA(x2l_l,'Faxa_snowl_18O')
         index_x2l_Faxa_snowl_HDO = mct_aVect_indexRA(x2l_l,'Faxa_snowl_HDO')
         index_x2l_Flrr_flood_16O = mct_aVect_indexRA(x2l_l,'Flrr_flood_16O')
         index_x2l_Flrr_flood_18O = mct_aVect_indexRA(x2l_l,'Flrr_flood_18O')
         index_x2l_Flrr_flood_HDO = mct_aVect_indexRA(x2l_l,'Flrr_flood_HDO')
      end if

      lSize = mct_avect_lSize(x2l_l)
      ic = c_lnd_ls
      do n=1,lSize
         dl =  dom_l%data%rAttr(kArea,n) * frac_l%rAttr(kl,n)
         if = f_wrain_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_rainc_16O,n) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_rainl_16O,n)
         if = f_wrain_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_rainc_18O,n) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_rainl_18O,n)
         if = f_wrain_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_rainc_HDO,n) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_rainl_HDO,n)

         if = f_wsnow_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_snowc_16O,n) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_snowl_16O,n)
         if = f_wsnow_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_snowc_18O,n) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_snowl_18O,n)
         if = f_wsnow_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_snowc_HDO,n) + &
                                     dl*x2l_l%rAttr(index_x2l_Faxa_snowl_HDO,n)

         if = f_wroff_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*x2l_l%rAttr(index_x2l_Flrr_flood_16O,n)
         if = f_wroff_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*x2l_l%rAttr(index_x2l_Flrr_flood_18O,n)
         if = f_wroff_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     dl*x2l_l%rAttr(index_x2l_Flrr_flood_HDO,n)
      end do
   end if

   first_time = .false.

end subroutine seq_diag_lnd_wiso_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_rtm_wiso_mct - compute global rtm input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global rtm input/output flux diagnostics for water isotopes
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!    2014-jul-23 - J. Zhu - modified
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_rtm_wiso_mct( dom_r, frac_r, r2x_r, x2r_r)

   type(mct_gGrid),intent(in)          ::  dom_r ! model domain
   type(mct_aVect),intent(in)          ::  frac_r ! rtm fractions
   type(mct_aVect),intent(in)          ::  r2x_r ! model to drv bundle
   type(mct_aVect),intent(in)          ::  x2r_r ! drv to model bundle

!EOP

   !----- local -----
   integer(in)      :: k,n,ic,if,ip      ! generic index
   integer(in)      :: kArea             ! index of area field in aVect
   integer(in)      :: kLat              ! index of lat field in aVect
   integer(in)      :: kl,ka,ko,ki,kr    ! fraction indices
   integer(in)      :: lSize             ! size of aVect
   real(r8)         :: da,di,do,dl,dr    ! area of a grid cell
   logical,save     :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_rtm_wiso_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   if (first_time) then
      index_x2r_Flrl_rofliq_16O = mct_aVect_indexRA(x2r_r,'Flrl_rofliq_16O')
      index_x2r_Flrl_rofliq_18O = mct_aVect_indexRA(x2r_r,'Flrl_rofliq_18O')
      index_x2r_Flrl_rofliq_HDO = mct_aVect_indexRA(x2r_r,'Flrl_rofliq_HDO')
      index_x2r_Flrl_rofice_16O = mct_aVect_indexRA(x2r_r,'Flrl_rofice_16O')
      index_x2r_Flrl_rofice_18O = mct_aVect_indexRA(x2r_r,'Flrl_rofice_18O')
      index_x2r_Flrl_rofice_HDO = mct_aVect_indexRA(x2r_r,'Flrl_rofice_HDO')
   end if

   ip = p_inst
   ic = c_rof_rr
   kArea = mct_aVect_indexRA(dom_r%data,afldname)
   lSize = mct_avect_lSize(x2r_r)
   do n=1,lSize
      dr =  dom_r%data%rAttr(kArea,n)
      if = f_wroff_16O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*x2r_r%rAttr(index_x2r_Flrl_rofliq_16O,n)
      if = f_wroff_18O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*x2r_r%rAttr(index_x2r_Flrl_rofliq_18O,n)
      if = f_wroff_HDO;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*x2r_r%rAttr(index_x2r_Flrl_rofliq_HDO,n)

      if = f_wioff_16O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*x2r_r%rAttr(index_x2r_Flrl_rofice_16O,n)
      if = f_wioff_18O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*x2r_r%rAttr(index_x2r_Flrl_rofice_18O,n)
      if = f_wioff_HDO;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*x2r_r%rAttr(index_x2r_Flrl_rofice_HDO,n)
   end do

   if (first_time) then
      index_r2x_Forr_roff_16O   = mct_aVect_indexRA(r2x_r,'Forr_roff_16O')
      index_r2x_Forr_roff_18O   = mct_aVect_indexRA(r2x_r,'Forr_roff_18O')
      index_r2x_Forr_roff_HDO   = mct_aVect_indexRA(r2x_r,'Forr_roff_HDO')
      index_r2x_Forr_ioff_16O   = mct_aVect_indexRA(r2x_r,'Forr_ioff_16O')
      index_r2x_Forr_ioff_18O   = mct_aVect_indexRA(r2x_r,'Forr_ioff_18O')
      index_r2x_Forr_ioff_HDO   = mct_aVect_indexRA(r2x_r,'Forr_ioff_HDO')
      index_r2x_Flrr_flood_16O  = mct_aVect_indexRA(r2x_r,'Flrr_flood_16O')
      index_r2x_Flrr_flood_18O  = mct_aVect_indexRA(r2x_r,'Flrr_flood_18O')
      index_r2x_Flrr_flood_HDO  = mct_aVect_indexRA(r2x_r,'Flrr_flood_HDO')
   end if

   ip = p_inst
   ic = c_rof_rs
   kArea = mct_aVect_indexRA(dom_r%data,afldname)
   lSize = mct_avect_lSize(r2x_r)
   do n=1,lSize
      dr =  dom_r%data%rAttr(kArea,n)
      if = f_wroff_16O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                  dr*r2x_r%rAttr(index_r2x_Forr_roff_16O,n)
      if = f_wroff_18O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                  dr*r2x_r%rAttr(index_r2x_Forr_roff_18O,n)
      if = f_wroff_HDO;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                  dr*r2x_r%rAttr(index_r2x_Forr_roff_HDO,n)

      if = f_wioff_16O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                  dr*r2x_r%rAttr(index_r2x_Forr_ioff_16O,n)
      if = f_wioff_18O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                  dr*r2x_r%rAttr(index_r2x_Forr_ioff_18O,n)
      if = f_wioff_HDO;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                  dr*r2x_r%rAttr(index_r2x_Forr_ioff_HDO,n)

      if = f_wroff_16O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*r2x_r%rAttr(index_r2x_Flrr_flood_16O,n)
      if = f_wroff_18O;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*r2x_r%rAttr(index_r2x_Flrr_flood_18O,n)
      if = f_wroff_HDO;
      budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                  dr*r2x_r%rAttr(index_r2x_Flrr_flood_HDO,n)
   end do

   first_time = .false.

end subroutine seq_diag_rtm_wiso_mct

!BOP ===========================================================================
!
! !IROUTINE: seq_diag_ocn_wiso_mct - compute global ocn input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global ocn input/output flux diagnostics for water isotopes
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!    2014-jul-23 - J. Zhu - modified
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_ocn_wiso_mct( dom_o, frac_o, o2x_o, x2o_o, xao_o, r2x_o)

   type(mct_gGrid),intent(in)          ::  dom_o ! model domain
   type(mct_aVect),intent(in)          :: frac_o ! frac bundle
   type(mct_aVect),intent(in),optional ::  o2x_o ! model to drv bundle
   type(mct_aVect),intent(in),optional ::  x2o_o ! drv to model bundle
   type(mct_aVect),intent(in),optional ::  xao_o ! drv to model bundle
   type(mct_aVect),intent(in),optional ::  r2x_o ! roff to drv bundle

!EOP

   !----- local -----
   integer(in)      :: k,n,if,ic,ip      ! generic index
   integer(in)      :: kArea             ! index of area field in aVect
   integer(in)      :: kLat              ! index of lat field in aVect
   integer(in)      :: kl,ka,ko,ki       ! fraction indices
   integer(in)      :: lSize             ! size of aVect
   real(r8)         :: da,di,do,dl       ! area of a grid cell
   logical,save     :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_ocn_wiso_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. present(o2x_o) .and. .not. present(x2o_o) .and. .not. present(xao_o)) then
      call shr_sys_abort(subName//"ERROR: must input a bundle")
   end if

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   ip = p_inst

   kArea = mct_aVect_indexRA(dom_o%data,afldname)
   ko    = mct_aVect_indexRA(frac_o,ofracname)
   ki    = mct_aVect_indexRA(frac_o,ifracname)

   if (present(xao_o)) then
      if (first_time) then
         index_xao_Faox_evap_16O = mct_aVect_indexRA(xao_o,'Faox_evap_16O')
         index_xao_Faox_evap_18O = mct_aVect_indexRA(xao_o,'Faox_evap_18O')
         index_xao_Faox_evap_HDO = mct_aVect_indexRA(xao_o,'Faox_evap_HDO')
      end if

      lSize = mct_avect_lSize(xao_o)
      ic = c_ocn_or
      do n=1,lSize
         do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
         if = f_wevap_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     do*xao_o%rAttr(index_xao_Faox_evap_16O,n)
         if = f_wevap_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     do*xao_o%rAttr(index_xao_Faox_evap_18O,n)
         if = f_wevap_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     do*xao_o%rAttr(index_xao_Faox_evap_HDO,n)
      end do
   end if

   if (present(x2o_o)) then
      if (first_time) then
         index_x2o_Fioi_meltw_16O = mct_aVect_indexRA(x2o_o,'Fioi_meltw_16O')
         index_x2o_Fioi_meltw_18O = mct_aVect_indexRA(x2o_o,'Fioi_meltw_18O')
         index_x2o_Fioi_meltw_HDO = mct_aVect_indexRA(x2o_o,'Fioi_meltw_HDO')
         index_x2o_Faxa_rain_16O  = mct_aVect_indexRA(x2o_o,'Faxa_rain_16O')
         index_x2o_Faxa_rain_18O  = mct_aVect_indexRA(x2o_o,'Faxa_rain_18O')
         index_x2o_Faxa_rain_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_rain_HDO')
         index_x2o_Faxa_snow_16O  = mct_aVect_indexRA(x2o_o,'Faxa_snow_16O')
         index_x2o_Faxa_snow_18O  = mct_aVect_indexRA(x2o_o,'Faxa_snow_18O')
         index_x2o_Faxa_snow_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_snow_HDO')
      end if

      lSize = mct_avect_lSize(x2o_o)
      ic = c_ocn_os
      do n=1,lSize
         do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
         di =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
         if = f_wmelt_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Fioi_meltw_16O,n)
         if = f_wmelt_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Fioi_meltw_18O,n)
         if = f_wmelt_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Fioi_meltw_HDO,n)

         if = f_wrain_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Faxa_rain_16O,n)
         if = f_wrain_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Faxa_rain_18O,n)
         if = f_wrain_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Faxa_rain_HDO,n)

         if = f_wsnow_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Faxa_snow_16O,n)
         if = f_wsnow_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Faxa_snow_18O,n)
         if = f_wsnow_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*x2o_o%rAttr(index_x2o_Faxa_snow_HDO,n)
      end do
   end if

   if (present(r2x_o)) then
      if (first_time) then
         index_r2x_Forr_roff_16O = mct_aVect_indexRA(r2x_o,'Forr_roff_16O')
         index_r2x_Forr_roff_18O = mct_aVect_indexRA(r2x_o,'Forr_roff_18O')
         index_r2x_Forr_roff_HDO = mct_aVect_indexRA(r2x_o,'Forr_roff_HDO')
         index_r2x_Forr_ioff_16O = mct_aVect_indexRA(r2x_o,'Forr_ioff_16O')
         index_r2x_Forr_ioff_18O = mct_aVect_indexRA(r2x_o,'Forr_ioff_18O')
         index_r2x_Forr_ioff_HDO = mct_aVect_indexRA(r2x_o,'Forr_ioff_HDO')
      end if

      lSize = mct_avect_lSize(r2x_o)
      ic = c_ocn_os
      do n=1,lSize
         do =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ko,n)
         di =  dom_o%data%rAttr(kArea,n) * frac_o%rAttr(ki,n)
         if = f_wroff_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*r2x_o%rAttr(index_r2x_Forr_roff_16O,n)
         if = f_wroff_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*r2x_o%rAttr(index_r2x_Forr_roff_18O,n)
         if = f_wroff_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*r2x_o%rAttr(index_r2x_Forr_roff_HDO,n)

         if = f_wioff_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*r2x_o%rAttr(index_r2x_Forr_ioff_16O,n)
         if = f_wioff_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*r2x_o%rAttr(index_r2x_Forr_ioff_18O,n)
         if = f_wioff_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     (do+di)*r2x_o%rAttr(index_r2x_Forr_ioff_HDO,n)
      end do
   end if

   first_time = .false.

end subroutine seq_diag_ocn_wiso_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_ice_wiso_mct - compute global ice input/output flux diagnostics
!
! !DESCRIPTION:
!     Compute global ice input/output flux diagnostics for water isotopes
!
! !REVISION HISTORY:
!    2008-jul-10 - T. Craig - update
!    2014-jul-23 - J. Zhu - modified
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_diag_ice_wiso_mct( dom_i, frac_i, i2x_i, x2i_i)

   type(mct_gGrid),intent(in)          ::  dom_i ! model domain
   type(mct_aVect),intent(in)          :: frac_i ! frac bundle
   type(mct_aVect),intent(in),optional ::  i2x_i ! model to drv bundle
   type(mct_aVect),intent(in),optional ::  x2i_i ! drv to model bundle

!EOP

   !----- local -----
   integer(in)      :: k,n,ic,if,ip      ! generic index
   integer(in)      :: kArea             ! index of area field in aVect
   integer(in)      :: kLat              ! index of lat field in aVect
   integer(in)      :: kl,ka,ko,ki       ! fraction indices
   integer(in)      :: lSize             ! size of aVect
   real(r8)         :: da,di,do,dl       ! area of a grid cell
   logical,save     :: first_time = .true.

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_ice_wiso_mct) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. present(i2x_i) .and. .not. present(x2i_i)) then
      call shr_sys_abort(subName//"ERROR: must input a bundle")
   end if

   !---------------------------------------------------------------------------
   ! add values found in this bundle to the budget table
   !---------------------------------------------------------------------------

   ip = p_inst

   kArea = mct_aVect_indexRA(dom_i%data,afldname)
   kLat  = mct_aVect_indexRA(dom_i%data,latname)
   ki    = mct_aVect_indexRA(frac_i,ifracname)
   ko    = mct_aVect_indexRA(frac_i,ofracname)

   if (present(i2x_i)) then
         index_i2x_Fioi_meltw_16O   = mct_aVect_indexRA(i2x_i,'Fioi_meltw_16O')
         index_i2x_Fioi_meltw_18O   = mct_aVect_indexRA(i2x_i,'Fioi_meltw_18O')
         index_i2x_Fioi_meltw_HDO   = mct_aVect_indexRA(i2x_i,'Fioi_meltw_HDO')
         index_i2x_Faii_evap_16O    = mct_aVect_indexRA(i2x_i,'Faii_evap_16O')
         index_i2x_Faii_evap_18O    = mct_aVect_indexRA(i2x_i,'Faii_evap_18O')
         index_i2x_Faii_evap_HDO    = mct_aVect_indexRA(i2x_i,'Faii_evap_HDO')

      lSize = mct_avect_lSize(i2x_i)
      do n=1,lSize
         if (dom_i%data%rAttr(kLat,n) > 0.0_r8) then
            ic = c_inh_ir
         else
            ic = c_ish_ir
         endif
         di =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ki,n)
         if = f_wmelt_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     di*i2x_i%rAttr(index_i2x_Fioi_meltw_16O,n)
         if = f_wmelt_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     di*i2x_i%rAttr(index_i2x_Fioi_meltw_18O,n)
         if = f_wmelt_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) - &
                                     di*i2x_i%rAttr(index_i2x_Fioi_meltw_HDO,n)

         if = f_wevap_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*i2x_i%rAttr(index_i2x_Faii_evap_16O,n)
         if = f_wevap_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*i2x_i%rAttr(index_i2x_Faii_evap_18O,n)
         if = f_wevap_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*i2x_i%rAttr(index_i2x_Faii_evap_HDO,n)
      end do
   end if

   if (present(x2i_i)) then
      if (first_time) then
         index_x2i_Faxa_rain_16O   = mct_aVect_indexRA(x2i_i,'Faxa_rain_16O')  
         index_x2i_Faxa_rain_18O   = mct_aVect_indexRA(x2i_i,'Faxa_rain_18O')  
         index_x2i_Faxa_rain_HDO   = mct_aVect_indexRA(x2i_i,'Faxa_rain_HDO')  
         index_x2i_Faxa_snow_16O   = mct_aVect_indexRA(x2i_i,'Faxa_snow_16O')  
         index_x2i_Faxa_snow_18O   = mct_aVect_indexRA(x2i_i,'Faxa_snow_18O')  
         index_x2i_Faxa_snow_HDO   = mct_aVect_indexRA(x2i_i,'Faxa_snow_HDO')  
      end if

      lSize = mct_avect_lSize(x2i_i)
      do n=1,lSize
         if (dom_i%data%rAttr(kLat,n) > 0.0_r8) then
            ic = c_inh_is
         else
            ic = c_ish_is
         endif
         do =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ko,n)
         di =  dom_i%data%rAttr(kArea,n) * frac_i%rAttr(ki,n)
         if  = f_wrain_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*x2i_i%rAttr(index_x2i_Faxa_rain_16O,n)
         if  = f_wrain_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*x2i_i%rAttr(index_x2i_Faxa_rain_18O,n)
         if  = f_wrain_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*x2i_i%rAttr(index_x2i_Faxa_rain_HDO,n)

         if  = f_wsnow_16O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*x2i_i%rAttr(index_x2i_Faxa_snow_16O,n)
         if  = f_wsnow_18O;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*x2i_i%rAttr(index_x2i_Faxa_snow_18O,n)
         if  = f_wsnow_HDO;
         budg_dataL_wiso(if,ic,ip) = budg_dataL_wiso(if,ic,ip) + &
                                     di*x2i_i%rAttr(index_x2i_Faxa_snow_HDO,n)

      end do
   end if

   first_time = .false.

end subroutine seq_diag_ice_wiso_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_diag_print_wiso_mct - print global budget diagnostics
!
! !DESCRIPTION:
!   Print global budget diagnostics for water isotopes
!
! !REVISION HISTORY:
!    2014-jul-23 - J. Zhu - modified
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_diag_print_wiso_mct(EClock,stop_alarm, &
   budg_print_inst, budg_print_daily, budg_print_month, &
   budg_print_ann, budg_print_ltann, budg_print_ltend)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock), intent(in) :: EClock
   logical,intent(in)           :: stop_alarm
   integer,intent(in)           :: budg_print_inst
   integer,intent(in)           :: budg_print_daily
   integer,intent(in)           :: budg_print_month
   integer,intent(in)           :: budg_print_ann
   integer,intent(in)           :: budg_print_ltann
   integer,intent(in)           :: budg_print_ltend

!EOP

   !--- local ---
   integer(in)      :: ic,if,ip    ! data array indicies
   integer(in)      :: ica,icl,icn,ics,ico
   integer(in)      :: icar,icxs,icxr,icas
   integer(in)      :: n           ! loop counter
   integer(in)      :: nday        ! number of days in time avg
   integer(in)      :: cdate,sec   ! coded date, seconds
   integer(in)      :: yr,mon,day  ! date
   integer(in)      :: iam         ! pe number
   integer(in)      :: plev        ! print level
   logical          :: sumdone     ! has a sum been computed yet
   character(len=40):: str         ! string
   real(r8) :: dataGpr (f_size,c_size,p_size) ! values to print, scaled and such

   !----- formats -----
   character(*),parameter :: subName = '(seq_diag_print_wiso_mct) '
   character(*),parameter :: F00   = "('(seq_diag_print_wiso_mct) ',4a)"

   !----- formats -----
   character(*),parameter :: FAH="(4a,i9,i6)"
   character(*),parameter :: FA0="('    ',8x,6(6x,a8,1x))"
   character(*),parameter :: FA1="('    ',a8,6f15.8)"
   character(*),parameter :: FA0r="('    ',8x,7(6x,a8,1x))"
   character(*),parameter :: FA1r="('    ',a8,7f15.8)"

!-------------------------------------------------------------------------------
! print instantaneous budget data
!-------------------------------------------------------------------------------

   sumdone = .false.
   call seq_comm_setptrs(CPLID,iam=iam)
   call seq_timemgr_EClockGetData(EClock,curr_yr=yr, &
      curr_mon=mon,curr_day=day,curr_tod=sec)
   cdate = yr*10000+mon*100+day

   do ip = 1,p_size
      plev = 0
      if (ip == p_inst) then
         plev = max(plev,budg_print_inst)
      endif
      if (ip==p_day .and. sec==0) then
         plev = max(plev,budg_print_daily)
      endif
      if (ip==p_mon .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_month)
      endif
      if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_ann)
      endif
      if (ip==p_inf .and. mon==1 .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_ltann)
      endif
      if (ip==p_inf .and. stop_alarm) then
         plev = max(plev,budg_print_ltend)
      endif

   if (plev > 0) then
! ---- doprint ---- doprint ---- doprint ----

   if (.not.sumdone) then
      call seq_diag_sum0_wiso_mct()
      dataGpr = budg_dataG_wiso
      sumdone = .true.

   !  old budget normalizations (global area and 1e6 for water)
      dataGpr = dataGpr/(4.0_r8*shr_const_pi)
      dataGpr = dataGpr * 1.0e6_r8
      dataGpr = dataGpr/budg_ns_wiso

      if (iam /= 0) return
   endif

   ! ---------------------------------------------------------
   ! ---- detail atm budgets and breakdown into components ---
   ! ---------------------------------------------------------

   if (plev >= 3) then
   do ic = 1,2
      if (ic == 1) then
         ica = c_atm_ar
         icl = c_lnd_ar
         icn = c_inh_ar
         ics = c_ish_ar
         ico = c_ocn_ar
         str = "ATM_to_CPL"
      elseif (ic == 2) then
         ica = c_atm_as
         icl = c_lnd_as
         icn = c_inh_as
         ics = c_ish_as
         ico = c_ocn_as
         str = "CPL_TO_ATM"
      else
         call shr_sys_abort(subname//' ERROR in ic index code 411')
      endif

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' H216O WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
      do if = 1,f_18O-1
         write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
                      dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
                                         dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
                      dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
      enddo
      write(logunit,FA1)    '   *SUM*', sum(dataGpr(1:f_18O-1,ica,ip)),sum(dataGpr(1:f_18O-1,icl,ip)), &
         sum(dataGpr(1:f_18O-1,icn,ip)),sum(dataGpr(1:f_18O-1,ics,ip)),sum(dataGpr(1:f_18O-1,ico,ip)), &
                                        sum(dataGpr(1:f_18O-1,ica,ip))+sum(dataGpr(1:f_18O-1,icl,ip))+ &
         sum(dataGpr(1:f_18O-1,icn,ip))+sum(dataGpr(1:f_18O-1,ics,ip))+sum(dataGpr(1:f_18O-1,ico,ip)) 

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' H218O WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
      do if = f_18O,f_HDO-1
         write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
                      dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
                                         dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
                      dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
      enddo
      write(logunit,FA1)    '   *SUM*',     sum(dataGpr(f_18O:f_HDO-1,ica,ip)),sum(dataGpr(f_18O:f_HDO-1,icl,ip)), &
         sum(dataGpr(F_18O:f_HDO-1,icn,ip)),sum(dataGpr(f_18O:f_HDO-1,ics,ip)),sum(dataGpr(f_18O:f_HDO-1,ico,ip)), &
                                            sum(dataGpr(f_18O:f_HDO-1,ica,ip))+sum(dataGpr(f_18O:f_HDO-1,icl,ip))+ &
         sum(dataGpr(F_18O:f_HDO-2,icn,ip))+sum(dataGpr(f_18O:f_HDO-1,ics,ip))+sum(dataGpr(f_18O:f_HDO-1,ico,ip)) 

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' HDO WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
      do if = f_HDO,f_size
         write(logunit,FA1)    fname(if),dataGpr(if,ica,ip),dataGpr(if,icl,ip), &
                      dataGpr(if,icn,ip),dataGpr(if,ics,ip),dataGpr(if,ico,ip), &
                                         dataGpr(if,ica,ip)+dataGpr(if,icl,ip)+ &
                      dataGpr(if,icn,ip)+dataGpr(if,ics,ip)+dataGpr(if,ico,ip) 
      enddo
      write(logunit,FA1)    '   *SUM*',    sum(dataGpr(f_HDO:f_size,ica,ip)),sum(dataGpr(f_HDO:f_size,icl,ip)), &
         sum(dataGpr(f_HDO:f_size,icn,ip)),sum(dataGpr(f_HDO:f_size,ics,ip)),sum(dataGpr(f_HDO:f_size,ico,ip)), &
                                           sum(dataGpr(f_HDO:f_size,ica,ip))+sum(dataGpr(f_HDO:f_size,icl,ip))+ &
         sum(dataGpr(f_HDO:f_size,icn,ip))+sum(dataGpr(f_HDO:f_size,ics,ip))+sum(dataGpr(f_HDO:f_size,ico,ip)) 

   enddo
   endif   ! plev

   ! ---------------------------------------------------------
   ! ---- detail lnd/ocn/ice component budgets ----
   ! ---------------------------------------------------------

   if (plev >= 2) then
   do ic = 1,4
      if (ic == 1) then
         icar = c_lnd_ar
         icxs = c_lnd_ls
         icxr = c_lnd_lr
         icas = c_lnd_as
         str = "LND"
      elseif (ic == 2) then
         icar = c_ocn_ar
         icxs = c_ocn_os
         icxr = c_ocn_or
         icas = c_ocn_as
         str = "OCN"
      elseif (ic == 3) then
         icar = c_inh_ar
         icxs = c_inh_is
         icxr = c_inh_ir
         icas = c_inh_as
         str = "ICE_NH"
      elseif (ic == 4) then
         icar = c_ish_ar
         icxs = c_ish_is
         icxr = c_ish_ir
         icas = c_ish_as
         str = "ICE_SH"
      else
         call shr_sys_abort(subname//' ERROR in ic index code 412')
      endif

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' H216O WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
      do if = 1, f_18O-1
         write(logunit,FA1)    fname(if),-dataGpr(if,icar,ip),dataGpr(if,icxs,ip), &
                                          dataGpr(if,icxr,ip),-dataGpr(if,icas,ip), &
                                         -dataGpr(if,icar,ip)+dataGpr(if,icxs,ip)+ &
                                          dataGpr(if,icxr,ip)-dataGpr(if,icas,ip)
      enddo
      write(logunit,FA1)    '   *SUM*',-sum(dataGpr(1:f_18O-1,icar,ip)), sum(dataGpr(1:f_18O-1,icxs,ip)), &
                                        sum(dataGpr(1:f_18O-1,icxr,ip)),-sum(dataGpr(1:f_18O-1,icas,ip)), &
                                       -sum(dataGpr(1:f_18O-1,icar,ip))+ sum(dataGpr(1:f_18O-1,icxs,ip))+ &
                                        sum(dataGpr(1:f_18O-1,icxr,ip))- sum(dataGpr(1:f_18O-1,icas,ip))

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' H218O WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
      do if = f_18O, f_HDO-1
         write(logunit,FA1)    fname(if),-dataGpr(if,icar,ip),dataGpr(if,icxs,ip), &
                                          dataGpr(if,icxr,ip),-dataGpr(if,icas,ip), &
                                         -dataGpr(if,icar,ip)+dataGpr(if,icxs,ip)+ &
                                          dataGpr(if,icxr,ip)-dataGpr(if,icas,ip)
      enddo
      write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_18O:f_HDO-1,icar,ip)), sum(dataGpr(f_18O:f_HDO-1,icxs,ip)), &
                                        sum(dataGpr(f_18O:f_HDO-1,icxr,ip)),-sum(dataGpr(f_18O:f_HDO-1,icas,ip)), &
                                       -sum(dataGpr(f_18O:f_HDO-1,icar,ip))+ sum(dataGpr(f_18O:f_HDO-1,icxs,ip))+ &
                                        sum(dataGpr(f_18O:f_HDO-1,icxr,ip))- sum(dataGpr(f_18O:f_HDO-1,icas,ip))

      write(logunit,*) ' '
      write(logunit,FAH) subname,trim(str)//' HDO WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
      do if = f_HDO, f_size
         write(logunit,FA1)    fname(if),-dataGpr(if,icar,ip),dataGpr(if,icxs,ip), &
                                          dataGpr(if,icxr,ip),-dataGpr(if,icas,ip), &
                                         -dataGpr(if,icar,ip)+dataGpr(if,icxs,ip)+ &
                                          dataGpr(if,icxr,ip)-dataGpr(if,icas,ip)
      enddo
      write(logunit,FA1)    '   *SUM*',-sum(dataGpr(f_HDO:f_size,icar,ip)), sum(dataGpr(f_HDO:f_size,icxs,ip)), &
                                        sum(dataGpr(f_HDO:f_size,icxr,ip)),-sum(dataGpr(f_HDO:f_size,icas,ip)), &
                                       -sum(dataGpr(f_HDO:f_size,icar,ip))+ sum(dataGpr(f_HDO:f_size,icxs,ip))+ &
                                        sum(dataGpr(f_HDO:f_size,icxr,ip))- sum(dataGpr(f_HDO:f_size,icas,ip))

   enddo
   endif   ! plev

   ! ---------------------------------------------------------
   ! ---- net summary budgets ----
   ! ---------------------------------------------------------

   if (plev >= 1) then

      write(logunit,*) ' '
      write(logunit,FAH) subname,'NET H216O WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh',' *SUM*  '
      do if = 1, f_18O-1
         write(logunit,FA1r)    fname(if),dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip), &
                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip), &
                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip), &
                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip), &
                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip), &
                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip), &
                                          dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip)+ &
                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip)+ &
                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip)+ &
                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip)+ &
                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip)+ &
                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip)
      enddo
      write(logunit,FA1r) '   *SUM*',sum(dataGpr(1:f_18O-1,c_atm_ar,ip))+sum(dataGpr(1:f_18O-1,c_atm_as,ip)), &
                                     sum(dataGpr(1:f_18O-1,c_lnd_lr,ip))+sum(dataGpr(1:f_18O-1,c_lnd_ls,ip)), &
                                     sum(dataGpr(1:f_18O-1,c_rof_rr,ip))+sum(dataGpr(1:f_18O-1,c_rof_rs,ip)), &
                                     sum(dataGpr(1:f_18O-1,c_ocn_or,ip))+sum(dataGpr(1:f_18O-1,c_ocn_os,ip)), &
                                     sum(dataGpr(1:f_18O-1,c_inh_ir,ip))+sum(dataGpr(1:f_18O-1,c_inh_is,ip)), &
                                     sum(dataGpr(1:f_18O-1,c_ish_ir,ip))+sum(dataGpr(1:f_18O-1,c_ish_is,ip)), &
                                     sum(dataGpr(1:f_18O-1,c_atm_ar,ip))+sum(dataGpr(1:f_18O-1,c_atm_as,ip))+ &
                                     sum(dataGpr(1:f_18O-1,c_lnd_lr,ip))+sum(dataGpr(1:f_18O-1,c_lnd_ls,ip))+ &
                                     sum(dataGpr(1:f_18O-1,c_rof_rr,ip))+sum(dataGpr(1:f_18O-1,c_rof_rs,ip))+ &
                                     sum(dataGpr(1:f_18O-1,c_ocn_or,ip))+sum(dataGpr(1:f_18O-1,c_ocn_os,ip))+ &
                                     sum(dataGpr(1:f_18O-1,c_inh_ir,ip))+sum(dataGpr(1:f_18O-1,c_inh_is,ip))+ &
                                     sum(dataGpr(1:f_18O-1,c_ish_ir,ip))+sum(dataGpr(1:f_18O-1,c_ish_is,ip))

      write(logunit,*) ' '
      write(logunit,FAH) subname,'NET H218O WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh',' *SUM*  '
      do if = f_18O, f_HDO-1
         write(logunit,FA1r)    fname(if),dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip), &
                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip), &
                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip), &
                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip), &
                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip), &
                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip), &
                                          dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip)+ &
                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip)+ &
                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip)+ &
                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip)+ &
                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip)+ &
                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip)
      enddo
      write(logunit,FA1r) '   *SUM*',sum(dataGpr(f_18O:f_HDO-1,c_atm_ar,ip))+sum(dataGpr(f_18O:f_HDO-1,c_atm_as,ip)), &
                                     sum(dataGpr(f_18O:f_HDO-1,c_lnd_lr,ip))+sum(dataGpr(f_18O:f_HDO-1,c_lnd_ls,ip)), &
                                     sum(dataGpr(f_18O:f_HDO-1,c_rof_rr,ip))+sum(dataGpr(f_18O:f_HDO-1,c_rof_rs,ip)), &
                                     sum(dataGpr(f_18O:f_HDO-1,c_ocn_or,ip))+sum(dataGpr(f_18O:f_HDO-1,c_ocn_os,ip)), &
                                     sum(dataGpr(f_18O:f_HDO-1,c_inh_ir,ip))+sum(dataGpr(f_18O:f_HDO-1,c_inh_is,ip)), &
                                     sum(dataGpr(f_18O:f_HDO-1,c_ish_ir,ip))+sum(dataGpr(f_18O:f_HDO-1,c_ish_is,ip)), &
                                     sum(dataGpr(f_18O:f_HDO-1,c_atm_ar,ip))+sum(dataGpr(f_18O:f_HDO-1,c_atm_as,ip))+ &
                                     sum(dataGpr(f_18O:f_HDO-1,c_lnd_lr,ip))+sum(dataGpr(f_18O:f_HDO-1,c_lnd_ls,ip))+ &
                                     sum(dataGpr(f_18O:f_HDO-1,c_rof_rr,ip))+sum(dataGpr(f_18O:f_HDO-1,c_rof_rs,ip))+ &
                                     sum(dataGpr(f_18O:f_HDO-1,c_ocn_or,ip))+sum(dataGpr(f_18O:f_HDO-1,c_ocn_os,ip))+ &
                                     sum(dataGpr(f_18O:f_HDO-1,c_inh_ir,ip))+sum(dataGpr(f_18O:f_HDO-1,c_inh_is,ip))+ &
                                     sum(dataGpr(f_18O:f_HDO-1,c_ish_ir,ip))+sum(dataGpr(f_18O:f_HDO-1,c_ish_is,ip))

      write(logunit,*) ' '
      write(logunit,FAH) subname,'NET HDO WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,sec
      write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh',' *SUM*  '
      do if = f_HDO, f_size
         write(logunit,FA1r)    fname(if),dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip), &
                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip), &
                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip), &
                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip), &
                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip), &
                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip), &
                                          dataGpr(if,c_atm_ar,ip)+dataGpr(if,c_atm_as,ip)+ &
                                          dataGpr(if,c_lnd_lr,ip)+dataGpr(if,c_lnd_ls,ip)+ &
                                          dataGpr(if,c_rof_rr,ip)+dataGpr(if,c_rof_rs,ip)+ &
                                          dataGpr(if,c_ocn_or,ip)+dataGpr(if,c_ocn_os,ip)+ &
                                          dataGpr(if,c_inh_ir,ip)+dataGpr(if,c_inh_is,ip)+ &
                                          dataGpr(if,c_ish_ir,ip)+dataGpr(if,c_ish_is,ip)
      enddo
      write(logunit,FA1r) '   *SUM*',sum(dataGpr(f_HDO:f_size,c_atm_ar,ip))+sum(dataGpr(f_HDO:f_size,c_atm_as,ip)), &
                                     sum(dataGpr(f_HDO:f_size,c_lnd_lr,ip))+sum(dataGpr(f_HDO:f_size,c_lnd_ls,ip)), &
                                     sum(dataGpr(f_HDO:f_size,c_rof_rr,ip))+sum(dataGpr(f_HDO:f_size,c_rof_rs,ip)), &
                                     sum(dataGpr(f_HDO:f_size,c_ocn_or,ip))+sum(dataGpr(f_HDO:f_size,c_ocn_os,ip)), &
                                     sum(dataGpr(f_HDO:f_size,c_inh_ir,ip))+sum(dataGpr(f_HDO:f_size,c_inh_is,ip)), &
                                     sum(dataGpr(f_HDO:f_size,c_ish_ir,ip))+sum(dataGpr(f_HDO:f_size,c_ish_is,ip)), &
                                     sum(dataGpr(f_HDO:f_size,c_atm_ar,ip))+sum(dataGpr(f_HDO:f_size,c_atm_as,ip))+ &
                                     sum(dataGpr(f_HDO:f_size,c_lnd_lr,ip))+sum(dataGpr(f_HDO:f_size,c_lnd_ls,ip))+ &
                                     sum(dataGpr(f_HDO:f_size,c_rof_rr,ip))+sum(dataGpr(f_HDO:f_size,c_rof_rs,ip))+ &
                                     sum(dataGpr(f_HDO:f_size,c_ocn_or,ip))+sum(dataGpr(f_HDO:f_size,c_ocn_os,ip))+ &
                                     sum(dataGpr(f_HDO:f_size,c_inh_ir,ip))+sum(dataGpr(f_HDO:f_size,c_inh_is,ip))+ &
                                     sum(dataGpr(f_HDO:f_size,c_ish_ir,ip))+sum(dataGpr(f_HDO:f_size,c_ish_is,ip))

   endif

   write(logunit,*) ' '
! ---- doprint ---- doprint ---- doprint ----
   endif  ! plev > 0
   enddo  ! ip = 1,p_size

end subroutine seq_diag_print_wiso_mct

!===============================================================================
end module seq_diag_wiso_mct
