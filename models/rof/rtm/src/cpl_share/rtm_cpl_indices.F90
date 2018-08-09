module rtm_cpl_indices
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: rtm_cpl_indices
!
! !DESCRIPTION:
!    Module containing the indices for the fields passed between RTM and
!    the driver. 
!
! !USES:

  use shr_sys_mod,    only : shr_sys_abort
  use RtmVar,         only : wiso_runoff
  implicit none

  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:

  public :: rtm_cpl_indices_set        ! Set the coupler indices

!
! !PUBLIC DATA MEMBERS:
!
  integer, public :: index_x2r_Flrl_rofliq  = 0   ! lnd->rtm liquid runoff forcing from land
  integer, public :: index_x2r_Flrl_rofice  = 0   ! lnd->rtm ice runoff forcing from land

  integer, public :: index_x2r_Flrl_rofliq_16O  = 0     ! lnd->rtm H216O water flux in liquid runoff
  integer, public :: index_x2r_Flrl_rofice_16O  = 0     ! lnd->rtm H216O water flux in ice runoff
  integer, public :: index_x2r_Flrl_rofliq_18O  = 0     ! lnd->rtm H218O water flux in liquid runoff
  integer, public :: index_x2r_Flrl_rofice_18O  = 0     ! lnd->rtm H218O water flux in ice runoff
  integer, public :: index_x2r_Flrl_rofliq_HDO  = 0     ! lnd->rtm HDO water flux in liquid runoff
  integer, public :: index_x2r_Flrl_rofice_HDO  = 0     ! lnd->rtm HDO water flux in ice runoff

  integer, public :: nflds_x2r = 0

  !TODO - nt_rtm and rtm_tracers need to be removed and set by access to the index array
  integer, public :: nt_rtm = 8                         ! number of tracers
  character(len=7), parameter, public :: rtm_tracers(8) = (/character(len=7) :: 'LIQ','ICE', &
                                'LIQ_16O', 'ICE_16O','LIQ_18O','ICE_18O','LIQ_HDO','ICE_HDO'/)

  ! roff to driver (part of land for now) (optional if RTM is off)

  integer, public :: index_r2x_Forr_roff  = 0   ! rtm->ocn liquid runoff to ocean
  integer, public :: index_r2x_Forr_ioff  = 0   ! rtm->ocn ice runoff to ocean
  integer, public :: index_r2x_Flrr_flood = 0   ! rtm->lnd flood runoff (>fthresh) back to land
  integer, public :: index_r2x_Slrr_volr = 0   ! rtm->lnd volr back to land
  integer, public :: nflds_r2x = 0

  integer, public :: index_r2x_Forr_roff_16O    = 0     ! rtm->ocn H216O in liquid runoff to ocean
  integer, public :: index_r2x_Forr_ioff_16O    = 0     ! rtm->ocn H216O in ice runoff to ocean
  integer, public :: index_r2x_Flrr_flood_16O   = 0     ! rtm->lnd H216O in flood runoff back to land
  integer, public :: index_r2x_Slrr_volr_16O    = 0     ! rtm->lnd H216O in volr back to land
  
  integer, public :: index_r2x_Forr_roff_18O    = 0     ! rtm->ocn H218O in liquid runoff to ocean
  integer, public :: index_r2x_Forr_ioff_18O    = 0     ! rtm->ocn H218O in ice runoff to ocean
  integer, public :: index_r2x_Flrr_flood_18O   = 0     ! rtm->lnd H218O in flood runoff back to land
  integer, public :: index_r2x_Slrr_volr_18O    = 0     ! rtm->lnd H218O in volr back to land

  integer, public :: index_r2x_Forr_roff_HDO    = 0     ! rtm->ocn HDO in liquid runoff to ocean
  integer, public :: index_r2x_Forr_ioff_HDO    = 0     ! rtm->ocn HDO in ice runoff to ocean
  integer, public :: index_r2x_Flrr_flood_HDO   = 0     ! rtm->lnd HDO in flood runoff back to land
  integer, public :: index_r2x_Slrr_volr_HDO    = 0     ! rtm->lnd HDO in volr back to land

!=======================================================================
contains

!=======================================================================


  subroutine rtm_cpl_indices_set( )


    !-----------------------------------------------------------------------
    ! !DESCRIPTION: 
    ! Set the coupler indices needed by the rof model coupler interface.
    ! runoff - (rtm -> ocn) and (rtm->lnd)
    !
    ! !USES:
    use seq_flds_mod  , only: seq_flds_r2x_fields, seq_flds_x2r_fields
    use mct_mod       , only: mct_aVect, mct_aVect_init, mct_avect_indexra, &
                              mct_aVect_clean, mct_avect_nRattr
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !REVISION HISTORY:
    ! Author: Mariana Vertenstein
    !
    ! !LOCAL VARIABLES:
    type(mct_aVect)   :: avtmp      ! temporary av
    character(len=32) :: subname = 'rtm_cpl_indices_set'  ! subroutine name
    !-----------------------------------------------------------------------

    ! x2r

    call mct_aVect_init(avtmp, rList=seq_flds_x2r_fields, lsize=1)

    index_x2r_Flrl_rofliq = mct_avect_indexra(avtmp,'Flrl_rofliq')
    index_x2r_Flrl_rofice = mct_avect_indexra(avtmp,'Flrl_rofice')

    ! adding water isotopes
    if ( wiso_runoff ) then
       index_x2r_Flrl_rofliq_16O   = mct_avect_indexra(avtmp,'Flrl_rofliq_16O')
       index_x2r_Flrl_rofice_16O   = mct_avect_indexra(avtmp,'Flrl_rofice_16O')
       index_x2r_Flrl_rofliq_18O   = mct_avect_indexra(avtmp,'Flrl_rofliq_18O')
       index_x2r_Flrl_rofice_18O   = mct_avect_indexra(avtmp,'Flrl_rofice_18O')
       index_x2r_Flrl_rofliq_HDO   = mct_avect_indexra(avtmp,'Flrl_rofliq_HDO')
       index_x2r_Flrl_rofice_HDO   = mct_avect_indexra(avtmp,'Flrl_rofice_HDO')

       ! index_x2r_Flrl_rofliq_16O   = mct_avect_indexra(avtmp,'Flrl_rofliq')
       ! index_x2r_Flrl_rofice_16O   = mct_avect_indexra(avtmp,'Flrl_rofice')
       ! index_x2r_Flrl_rofliq_18O   = mct_avect_indexra(avtmp,'Flrl_rofliq')
       ! index_x2r_Flrl_rofice_18O   = mct_avect_indexra(avtmp,'Flrl_rofice')
       ! index_x2r_Flrl_rofliq_HDO   = mct_avect_indexra(avtmp,'Flrl_rofliq')
       ! index_x2r_Flrl_rofice_HDO   = mct_avect_indexra(avtmp,'Flrl_rofice')
    end if

    nflds_x2r = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

    ! r2x

    call mct_aVect_init(avtmp, rList=seq_flds_r2x_fields, lsize=1)

    index_r2x_Forr_roff  = mct_avect_indexra(avtmp,'Forr_roff')
    index_r2x_Forr_ioff  = mct_avect_indexra(avtmp,'Forr_ioff')
    index_r2x_Flrr_flood = mct_avect_indexra(avtmp,'Flrr_flood')
    index_r2x_Slrr_volr = mct_avect_indexra(avtmp,'Slrr_volr')

    ! adding water isotopes
    if ( wiso_runoff ) then
       index_r2x_Forr_roff_16O     = mct_avect_indexra(avtmp,'Forr_roff_16O')
       index_r2x_Forr_ioff_16O     = mct_avect_indexra(avtmp,'Forr_ioff_16O')
       index_r2x_Forr_roff_18O     = mct_avect_indexra(avtmp,'Forr_roff_18O')
       index_r2x_Forr_ioff_18O     = mct_avect_indexra(avtmp,'Forr_ioff_18O')
       index_r2x_Forr_roff_HDO     = mct_avect_indexra(avtmp,'Forr_roff_HDO')
       index_r2x_Forr_ioff_HDO     = mct_avect_indexra(avtmp,'Forr_ioff_HDO')

       index_r2x_Flrr_flood_16O    = mct_avect_indexra(avtmp,'Flrr_flood_16O')
       index_r2x_Flrr_flood_18O    = mct_avect_indexra(avtmp,'Flrr_flood_18O')
       index_r2x_Flrr_flood_HDO    = mct_avect_indexra(avtmp,'Flrr_flood_HDO')
       index_r2x_Slrr_volr_16O     = mct_avect_indexra(avtmp,'Slrr_volr_16O')
       index_r2x_Slrr_volr_18O     = mct_avect_indexra(avtmp,'Slrr_volr_18O')
       index_r2x_Slrr_volr_HDO     = mct_avect_indexra(avtmp,'Slrr_volr_HDO')
    end if

    nflds_r2x = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

  end subroutine rtm_cpl_indices_set

end module rtm_cpl_indices
