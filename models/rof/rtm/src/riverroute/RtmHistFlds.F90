module RtmHistFlds

!-----------------------------------------------------------------------
! !DESCRIPTION:
! Module containing initialization of RTM history fields and files
! This is the module that the user must modify in order to add new
! history fields or modify defaults associated with existing history
! fields.
!
! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use RunoffMod      , only : runoff
  use RtmHistFile    , only : RtmHistAddfld, RtmHistPrintflds
  use rtm_cpl_indices, only : nt_rtm, rtm_tracers  
  use RtmVar         , only : iulog, wiso_runoff

  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: RtmHistFldsInit 
  public :: RtmHistFldsSet
!
!------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

  subroutine RtmHistFldsInit()

    !-------------------------------------------------------
    ! DESCRIPTION:
    ! Build master field list of all possible fields in a history file.
    ! Each field has associated with it a ``long\_name'' netcdf attribute that
    ! describes what the field is, and a ``units'' attribute. A subroutine is
    ! called to add each field to the masterlist.
    !
    ! !USES:
    ! ARGUMENTS:
    implicit none

    !-------------------------------------------------------

      !! call RtmHistAddfld (fname='QCHANR', units='m3/s',  &
      call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(1)), units='m3/s',  &
           avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(1)), &
           ptr_rof=runoff%runofflnd_nt1)
  
      call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
           avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(2)), &
           ptr_rof=runoff%runofflnd_nt2)
  
      !! call RtmHistAddfld (fname='QCHOCNR', units='m3/s', &
      call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(1)), units='m3/s', &
           avgflag='A', long_name='RTM river discharge into ocean: '//trim(rtm_tracers(1)), &
           ptr_rof=runoff%runoffocn_nt1)
  
      call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(2)), units='m3/s', &
           avgflag='A', long_name='RTM river discharge into ocean: '//trim(rtm_tracers(2)), &
           ptr_rof=runoff%runoffocn_nt2)
  
      call RtmHistAddfld (fname='VOLR'//'_'//trim(rtm_tracers(1)), units='m3',  &
           avgflag='A', long_name='RTM storage: '//trim(rtm_tracers(1)), &
           ptr_rof=runoff%volr_nt1, default='inactive')
  
      call RtmHistAddfld (fname='VOLR'//'_'//trim(rtm_tracers(2)), units='m3',  &
           avgflag='A', long_name='RTM storage: '//trim(rtm_tracers(2)), &
           ptr_rof=runoff%volr_nt2, default='inactive')
  
      call RtmHistAddfld (fname='DVOLRDT_LND', units='mm/s',  &
           avgflag='A', long_name='RTM land change in storage: '//trim(rtm_tracers(1)), &
           ptr_rof=runoff%dvolrdtlnd_nt1, default='inactive')
  
      call RtmHistAddfld (fname='DVOLRDT_LND'//'_'//trim(rtm_tracers(2)), units='mm/s',  &
           avgflag='A', long_name='RTM land change in storage: '//trim(rtm_tracers(2)), &
           ptr_rof=runoff%dvolrdtlnd_nt2, default='inactive')
  
      call RtmHistAddfld (fname='DVOLRDT_OCN', units='mm/s',  &
           avgflag='A', long_name='RTM ocean change of storage: '//trim(rtm_tracers(1)), &
           ptr_rof=runoff%dvolrdtocn_nt1, default='inactive')
  
      call RtmHistAddfld (fname='DVOLRDT_OCN'//'_'//trim(rtm_tracers(2)), units='mm/s',  &
           avgflag='A', long_name='RTM ocean change of storage: '//trim(rtm_tracers(2)), &
           ptr_rof=runoff%dvolrdtocn_nt2, default='inactive')
  
      call RtmHistAddfld (fname='RTMFLOOD', units='m3/s',  &
           avgflag='A', long_name='RTM flooding flux', &
           ptr_rof=runoff%flood, default='inactive')

      ! adding water isotopes
      if ( wiso_runoff ) then
         call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(3)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(3)), &
              ptr_rof=runoff%runofflnd_nt3)
  
         call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(4)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(4)), &
              ptr_rof=runoff%runofflnd_nt4)
  
         call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(5)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(5)), &
              ptr_rof=runoff%runofflnd_nt5)
  
         call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(6)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(6)), &
              ptr_rof=runoff%runofflnd_nt6)
  
         call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(7)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(7)), &
              ptr_rof=runoff%runofflnd_nt7)
  
         call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(8)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(8)), &
              ptr_rof=runoff%runofflnd_nt8)

  
         call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(3)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(3)), &
              ptr_rof=runoff%runoffocn_nt3)
  
         call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(4)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(4)), &
              ptr_rof=runoff%runoffocn_nt4)
  
         call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(5)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(5)), &
              ptr_rof=runoff%runoffocn_nt5)
  
         call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(6)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(6)), &
              ptr_rof=runoff%runoffocn_nt6)
  
         call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(7)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(7)), &
              ptr_rof=runoff%runoffocn_nt7)
  
         call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(8)), units='m3/s',  &
              avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(8)), &
              ptr_rof=runoff%runoffocn_nt8)
      end if

    ! Print masterlist of history fields
    call RtmHistPrintflds()

  end subroutine RtmHistFldsInit

!-----------------------------------------------------------------------

  subroutine RtmHistFldsSet()

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Set rtm history fields as 1d poitner arrays
    !
    implicit none

    !-----------------------------------------------------------------------

    ! Currently only have two tracers

      runoff%runofflnd_nt1(:)  = runoff%runofflnd(:,1)
      runoff%runofflnd_nt2(:)  = runoff%runofflnd(:,2)
  
      runoff%runoffocn_nt1(:)  = runoff%runoffocn(:,1)
      runoff%runoffocn_nt2(:)  = runoff%runoffocn(:,2)
  
      runoff%dvolrdtlnd_nt1(:) = runoff%dvolrdtlnd(:,1)
      runoff%dvolrdtlnd_nt2(:) = runoff%dvolrdtlnd(:,2)
  
      runoff%dvolrdtocn_nt1(:) = runoff%dvolrdtocn(:,1)
      runoff%dvolrdtocn_nt2(:) = runoff%dvolrdtocn(:,2)
  
      runoff%volr_nt1(:)       = runoff%volrlnd(:,1)
      runoff%volr_nt2(:)       = runoff%volrlnd(:,2)

      ! adding water isotopes
      if ( wiso_runoff ) then
         runoff%runofflnd_nt3(:)  = runoff%runofflnd(:,3)
         runoff%runofflnd_nt4(:)  = runoff%runofflnd(:,4)
         runoff%runofflnd_nt5(:)  = runoff%runofflnd(:,5)
         runoff%runofflnd_nt6(:)  = runoff%runofflnd(:,6)
         runoff%runofflnd_nt7(:)  = runoff%runofflnd(:,7)
         runoff%runofflnd_nt8(:)  = runoff%runofflnd(:,8)

         runoff%runoffocn_nt3(:)  = runoff%runoffocn(:,3)
         runoff%runoffocn_nt4(:)  = runoff%runoffocn(:,4)
         runoff%runoffocn_nt5(:)  = runoff%runoffocn(:,5)
         runoff%runoffocn_nt6(:)  = runoff%runoffocn(:,6)
         runoff%runoffocn_nt7(:)  = runoff%runoffocn(:,7)
         runoff%runoffocn_nt8(:)  = runoff%runoffocn(:,8)

         runoff%volr_nt3(:)       = runoff%volrlnd(:,3)
         runoff%volr_nt4(:)       = runoff%volrlnd(:,4)
         runoff%volr_nt5(:)       = runoff%volrlnd(:,5)
         runoff%volr_nt6(:)       = runoff%volrlnd(:,6)
         runoff%volr_nt7(:)       = runoff%volrlnd(:,7)
         runoff%volr_nt8(:)       = runoff%volrlnd(:,8)
      end if

  end subroutine RtmHistFldsSet


end module RtmHistFlds
