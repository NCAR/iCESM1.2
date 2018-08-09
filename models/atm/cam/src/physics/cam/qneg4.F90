
subroutine qneg4 (subnam  ,lchnk   ,ncol    ,ztodt   ,        &
                  qbot    ,srfrpdel,shflx   ,lhflx   ,qflx    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check if moisture flux into the ground is exceeding the total
! moisture content of the lowest model layer (creating negative moisture
! values).  If so, then subtract the excess from the moisture and
! latent heat fluxes and add it to the sensible heat flux.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Olson
! 
! Water isotopes added by J. Nusbaumer - Mar 2011
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,    only: get_lat_p, get_lon_p
   use physconst,    only: gravit, latvap
   use constituents, only: qmin, pcnst
   use cam_logfile,  only: iulog

   !water isotopes:
   use water_types,   only: iwtvap
   use water_tracer_vars, only: trace_water, wtrc_iatype, wtrc_ntype, iwspec
   use water_tracers, only: wtrc_ratio

   implicit none

!
! Input arguments
!
   character*8, intent(in) :: subnam         ! name of calling routine
!
   integer, intent(in) :: lchnk              ! chunk index
   integer, intent(in) :: ncol               ! number of atmospheric columns
!
   real(r8), intent(in) :: ztodt             ! two times model timestep (2 delta-t)
   real(r8), intent(in) :: qbot(pcols,pcnst) ! moisture at lowest model level
   real(r8), intent(in) :: srfrpdel(pcols)   ! 1./(pint(K+1)-pint(K))
!
! Input/Output arguments
!
   real(r8), intent(inout) :: shflx(pcols)   ! Surface sensible heat flux (J/m2/s)
   real(r8), intent(inout) :: lhflx(pcols)   ! Surface latent   heat flux (J/m2/s)
   real(r8), intent(inout) :: qflx (pcols,pcnst)   ! surface water flux (kg/m^2/s)
!
!---------------------------Local workspace-----------------------------
!
   integer :: i,ii              ! longitude indices
   integer :: iw                ! i index of worst violator
   integer :: indxexc(pcols)    ! index array of points with excess flux
   integer :: nptsexc           ! number of points with excess flux
   integer :: m                 ! loop control variable for water isotopes
   integer :: ivap              ! isotope index
!
   real(r8):: worst             ! biggest violator
   real(r8):: excess(pcols)     ! Excess downward sfc latent heat flux

!water isotopes:
   real(r8):: qfxo(pcols,pcnst)   ! initial tracer flux
   real(r8):: rat                 ! tracer ratio

!
!-----------------------------------------------------------------------
!

! Store old value to input for water tracers

   do m = 1, pcnst
     do i = 1, ncol
       qfxo(i,m) = qflx(i,m)
     end do
   end do

! Compute excess downward (negative) q flux compared to a theoretical
! maximum downward q flux.  The theoretical max is based upon the
! given moisture content of lowest level of the model atmosphere.
!
   nptsexc = 0
   do i = 1,ncol
      excess(i) = qflx(i,1) - (qmin(1) - qbot(i,1))/(ztodt*gravit*srfrpdel(i))
!
! If there is an excess downward (negative) q flux, then subtract
! excess from "qflx" and "lhflx" and add to "shflx".
!
      if (excess(i) < 0._r8) then
         nptsexc = nptsexc + 1
         indxexc(nptsexc) = i
         qflx (i,1) = qflx (i,1) - excess(i)
         lhflx(i) = lhflx(i) - excess(i)*latvap
         shflx(i) = shflx(i) + excess(i)*latvap
      end if
   end do
!
! Write out worst value if excess
!
   if (nptsexc.gt.0) then
      worst = 0._r8
      do ii=1,nptsexc
         i = indxexc(ii)
         if (excess(i) < worst) then
            worst = excess(i)
            iw = i
         end if
      end do
      write(iulog,9000) subnam,nptsexc,worst, lchnk, iw, get_lat_p(lchnk,iw),get_lon_p(lchnk,iw)
   end if
!
! Water tracers: where total has change, modify tracers to conserve ratios
!

   if (trace_water) then
     !NOTE:  qfxo may not be needed, as ratio is against H2O tracer, not q. - JN
     do ivap = 1, wtrc_ntype(iwtvap)
       m = wtrc_iatype(ivap, iwtvap)
      
       do ii = 1, nptsexc
         i = indxexc(ii)
!         rat = wtrc_ratio(iwspec(m), qfxo(i,m),qfxo(i,1))
         rat = wtrc_ratio(iwspec(m),qfxo(i,m),qfxo(i,wtrc_iatype(1,iwtvap)))
         qflx(i,m) = qflx(i,m) - rat*excess(i)
       end do
     end do
   end if
!

   return
9000 format(' QNEG4 WARNING from ',a8 &
            ,' Max possible LH flx exceeded at ',i4,' points. ' &
            ,', Worst excess = ',1pe12.4 &
            ,', lchnk = ',i3 &
            ,', i = ',i4 &
            ,', same as indices lat =', i4 &
            ,', lon =', i4 &
           )
end subroutine qneg4
