!=======================================================================
!
!BOP
!
! !MODULE: ice_isotope - Isotope tracer within sea ice
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors: David Bailey, NCAR
!          Jiang Zhu, UW-Madison
!
! 2014: Added i2x evaporation flux
!       Added fractionation options
!       Fixed bugs
!
! !INTERFACE:
!
      module ice_isotope
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_fileunits
      use ice_restart, only: lenstr, restart_dir, restart_file, &
                             pointer_file, runtype
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
!
!EOP
!
      implicit none

      logical (kind=log_kind) :: & 
         restart_iso      ! if .true., read isotope tracer restart file

      character(len=5), parameter ::    &
         frac = 'cfrac'   ! fractionation coefficient calculation method
                          !  cfrac, constant fractionation
                          !  nfrac, nonfractionation
                          !  gfrac, growth-rate dependent for H2_18O

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_isotope
!
! !DESCRIPTION:
!
!  Initialize ice isotope tracer (call prior to reading restart data)
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_isotope
!
! !USES:
!
      use ice_state, only: filename_iso 
!
!EOP
!

      if (trim(filename_iso) /= 'none') restart_iso = .true.

      if (restart_iso) then
         if (trim(runtype) == 'continue') then
            call read_restart_iso
         else
            call read_restart_iso(filename_iso)
         endif
      endif

      end subroutine init_isotope

!=======================================================================

!BOP
!
! !ROUTINE: update_isotope
!
! !DESCRIPTION:
!
!  Increase isotope in ice or snow surface due to deposition
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine update_isotope (nx_block, ny_block,  &
                                dt,       icells,     &
                                indxi,    indxj,      &
                                meltt,    melts,      &
                                meltb,    congel,     &
                                snoice,               &
                                fsnow,                &
                                fiso_Qref,            &
                                trcrn,                &
                                aice_old,             &
                                vice_old, vsno_old,   &
                                vicen, vsnon, aicen,  &
                                fiso_atm, fiso_evapn, &
                                fiso_ocnn, HDO_ocn,   &
                                H2_16O_ocn, H2_18O_ocn)
!
! !USES:
!
      use water_isotopes, only: wiso_alpi

      use ice_domain_size, only: max_ntrcr, nilyr, nslyr, n_iso, n_isomx
      use ice_state, only: nt_iso, nt_Tsfc

!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icells                 ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj           ! compressed indices for cells with ice

      real (kind=dbl_kind), intent(in) :: &
         dt                     ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         meltt,    &
         melts,    &
         meltb,    &
         congel,   &            ! congelation ice growth    (m/step)
         snoice,   &            ! ice thickness increase    (m/step)
         fsnow,    &            ! snowfall       (kg/m^2/s of water)
         vicen,    &            ! volume per unit area of ice    (m)
         vsnon,    &
         aicen,    &
         aice_old, &
         vice_old, &
         vsno_old, &
         HDO_ocn,    &
         H2_16O_ocn, &
         H2_18O_ocn 

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_isomx), &
         intent(in) ::  &
         fiso_atm,      &       ! isotopic snowfall (kg/m^2/s of water)
         fiso_Qref              ! isotope reference humidity

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_isomx), &
         intent(inout) :: &
         fiso_ocnn,     &       ! isotopic freshwater (kg/m^2/s)
         fiso_evapn             ! evaporative water flux (kg/m^2/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr), &
         intent(inout) :: &
         trcrn                  ! tracer concentration (kg/m^3 for isotopes)

!
!  local variables
!
      integer (kind=int_kind) :: i, j, ij, k
      integer (kind=int_kind) :: n  ! print_points
!
      real (kind=dbl_kind), dimension(icells) :: &
         dzssl,     &
         dzint,     &
         dzssli,    &
         dzinti,    &
         evaps,     &           ! evaporation over snow     (m/step)
         evapi,     &           ! evaporation over ice      (m/step)
         dhs_snoice,&           ! snow thickness reduction  (m/step)
         hi,        &
         hs

      real (kind=dbl_kind), dimension(icells,n_isomx) :: &
        isotot, isotot0         ! for diagnostics 

      real (kind=dbl_kind) :: &
         dzssl_new, &
         dzint_new, &
         dzssli_new, &
         dzinti_new, &
         dznew

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_isomx,2) :: &
         isosno, isoice, &      ! mass of isotopes  (kg)
         isosno0, isoice0       ! for diagnostic prints

      real (kind=dbl_kind) :: &
         hs_old, hi_old, hslyr_old, hilyr_old, dhs, dhi,        &
         hslyr, hilyr, sloss1, sloss2,                          &
         TsfK,      &           ! snow/ice temperature (K)
         alphai,    &           ! ice/vapour fractionation coefficient
         ratio,     &           ! isotopic ratio
         work,      &           ! temporary variable
         alpha

! These need to be the same as in the DE code. Put in a common place?
      real (kind=dbl_kind), parameter :: &
        hi_ssl = .050_dbl_kind, &
        hs_ssl = .040_dbl_kind

! initialize
      isosno(:,:,:,:)   = c0
      isoice(:,:,:,:)   = c0
      isosno0(:,:,:,:)  = c0
      isoice0(:,:,:,:)  = c0
      fiso_evapn(:,:,:) = c0
      fiso_ocnn(:,:,:)  = c0

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hs_old=vsno_old(i,j)/aice_old(i,j)
         hi_old=vice_old(i,j)/aice_old(i,j)
         hslyr_old=hs_old/real(nslyr,kind=dbl_kind)
         hilyr_old=hi_old/real(nilyr,kind=dbl_kind)

         dzssl(ij)=min(hslyr_old/c2,hs_ssl)
         dzint(ij)=hs_old-dzssl(ij)
         dzssli(ij)=min(hilyr_old/c2,hi_ssl)
         dzinti(ij)=hi_old-dzssli(ij)

         if (aicen(i,j) > puny) then
            hs(ij) = vsnon(i,j)/aicen(i,j)
            hi(ij) = vicen(i,j)/aicen(i,j)
            dhs_snoice(ij) = snoice(i,j)*rhoi/rhos
         elseif (aice_old(i,j) > puny) then
            hs(ij) = vsnon(i,j)/aice_old(i,j)
            hi(ij) = vicen(i,j)/aice_old(i,j)
            dhs_snoice(ij) = snoice(i,j)*rhoi/rhos
         endif

         evaps(ij) = hs(ij)-(hs_old-melts(i,j)-dhs_snoice(ij)+&
             fsnow(i,j)/rhos*dt)
         evapi(ij) = hi(ij)-(hi_old-meltt(i,j)-meltb(i,j)+ &
             congel(i,j)+snoice(i,j))
      enddo

      do ij = 1, icells
       i = indxi(ij)
       j = indxj(ij)

       do k=1,n_iso
         isosno(i,j,k,:)=&                      ! isotope in snow
          trcrn(i,j,nt_iso+(k-1)*4  :nt_iso+(k-1)*4+1)*vsno_old(i,j)
         isoice(i,j,k,:)=&                      ! isotope in ice
          trcrn(i,j,nt_iso+(k-1)*4+2:nt_iso+(k-1)*4+3)*vice_old(i,j)
         isosno0(i,j,k,:)=isosno(i,j,k,:)
         isoice0(i,j,k,:)=isoice(i,j,k,:)
         isotot0(ij,k)=isosno(i,j,k,2)+isosno(i,j,k,1) &
           +isoice(i,j,k,2)+isoice(i,j,k,1)
       enddo
      enddo

! condensation of vapor onto snow and ice
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         TsfK = trcrn(i,j,nt_Tsfc) + Tffresh

         if (evaps(ij) > c0) then   ! condensation to snow
            do k = 1, n_iso         
               ratio = c1   ! ratio between 18O(HDO) and 16O in humidity
               alphai = c1  ! fractionation coefficient
               if (frac.ne.'nfrac' .and. fiso_Qref(i,j,2)>puny) &
                  ratio = fiso_Qref(i,j,k)/fiso_Qref(i,j,2)
               if (frac.ne.'nfrac' .and. k==1) alphai = wiso_alpi(3,TsfK)
               if (frac.ne.'nfrac' .and. k==2) alphai = wiso_alpi(2,TsfK)
               if (frac.ne.'nfrac' .and. k==3) alphai = wiso_alpi(4,TsfK)
               work = alphai*ratio*rhos*evaps(ij)*aicen(i,j)
               fiso_evapn(i,j,k) = fiso_evapn(i,j,k)+work/dt
               isosno(i,j,k,1) = isosno(i,j,k,1)+work
            enddo
            dzssl(ij) = dzssl(ij)+evaps(ij)
         endif

         if (evapi(ij) > c0) then   ! condensation to ice
            do k = 1, n_iso         
               ratio = c1 ! ratio between 18O(HDO) and 16O in ref humidity
               alphai = c1  ! fractionation coefficient
               if (frac.ne.'nfrac' .and. fiso_Qref(i,j,2)>puny) &
                  ratio = fiso_Qref(i,j,k)/fiso_Qref(i,j,2)
               if (frac.ne.'nfrac' .and. k==1) alphai = wiso_alpi(3,TsfK)
               if (frac.ne.'nfrac' .and. k==2) alphai = wiso_alpi(2,TsfK)
               if (frac.ne.'nfrac' .and. k==3) alphai = wiso_alpi(4,TsfK)
               work = alphai*ratio*rhoi*evapi(ij)*aicen(i,j)
               fiso_evapn(i,j,k) = fiso_evapn(i,j,k)+work/dt
               isoice(i,j,k,1) = isoice(i,j,k,1)+work
            enddo
            dzssli(ij) = dzssli(ij)+evapi(ij)
         endif
      enddo                 ! icells

!     basal ice growth and isotope uptake
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (congel(i,j) > c0) then
            alpha = isoice_alpha(congel(i,j)/dt,'HDO',frac)
            work = alpha*HDO_ocn(i,j)*rhoi*congel(i,j)*aicen(i,j)
            isoice(i,j,1,2) = isoice(i,j,1,2)+work
            fiso_ocnn(i,j,1) = fiso_ocnn(i,j,1)-work/dt

            alpha = isoice_alpha(congel(i,j)/dt,'H2_16O',frac)
            work = alpha*H2_16O_ocn(i,j)*rhoi*congel(i,j)*aicen(i,j)
            isoice(i,j,2,2) = isoice(i,j,2,2)+work
            fiso_ocnn(i,j,2) = fiso_ocnn(i,j,2)-work/dt

            alpha = isoice_alpha(congel(i,j)/dt,'H2_18O',frac)
            work = alpha*H2_18O_ocn(i,j)*rhoi*congel(i,j)*aicen(i,j)
            isoice(i,j,3,2) = isoice(i,j,3,2)+work
            fiso_ocnn(i,j,3) = fiso_ocnn(i,j,3)-work/dt

            dzinti(ij) = dzinti(ij)+congel(i,j)
         endif
     enddo

! sublimation of snow and ice
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (evaps(ij) < c0) then   ! snow sublimation (no fractionation)
            do k = 1, n_iso         
               !ratio = c1 ! ratio between 18O(HDO) and 16O in snow ssl
               !if (isosno(i,j,2,1) > puny)      &
               !   ratio = isosno(i,j,k,1)/isosno(i,j,2,1)
               !if (ratio > c5) ratio = c1   !! remove latter?
               !work = ratio*rhos*evaps(ij)*aicen(i,j)
               !fiso_evapn(i,j,k) = fiso_evapn(i,j,k)+work/dt
               
               sloss1 = c0
               sloss2 = c0
               if (dzssl(ij) > puny)                &
                  sloss1 = isosno(i,j,k,1)*         &
                     min(-evaps(ij),dzssl(ij))/dzssl(ij)
               if (isosno(i,j,k,1) >= sloss1) then
                  isosno(i,j,k,1) = isosno(i,j,k,1)-sloss1
               else
                  sloss1 = isosno(i,j,k,1)
                  isosno(i,j,k,1) = c0
               endif
               if (dzint(ij) > puny)                &
                  sloss2 = isosno(i,j,k,2)*         &
                     max(-evaps(ij)-dzssl(ij),c0)/dzint(ij)
               if (isosno(i,j,k,2) >= sloss2) then
                  isosno(i,j,k,2) = isosno(i,j,k,2)-sloss2
               else
                  sloss2 = isosno(i,j,k,2)
                  isosno(i,j,k,2) = c0
               endif
               fiso_evapn(i,j,k) = fiso_evapn(i,j,k)-(sloss1+sloss2)/dt
            enddo

            dzint(ij) = dzint(ij)+min(dzssl(ij)+evaps(ij),c0)
            dzssl(ij) = max(dzssl(ij)+evaps(ij),c0)
            if ( dzssl(ij) <= c0) then  ! ssl goes away
               fiso_evapn(i,j,:) = fiso_evapn(i,j,:)-isosno(i,j,:,1)/dt !CRT
               isosno(i,j,:,1) = c0
               dzssl(ij) = c0
            endif
            if (dzint(ij) <= c0) then   ! int goes away
               fiso_evapn(i,j,:) = fiso_evapn(i,j,:)-isosno(i,j,:,2)/dt !CRT
               isosno(i,j,:,2) = c0
               dzint(ij) = c0
            endif
         endif

         if (evapi(ij) < c0) then   ! ice sublimation (no fractionation)
            do k = 1, n_iso         
               !!ratio = c1 ! ratio between 18O(HDO) and 16O in ice ssl
               !!if (isoice(i,j,2,1) > puny)      &
               !!   ratio = isoice(i,j,k,1)/isoice(i,j,2,1)
               !!if (ratio > c5) ratio = c1   ! remove latter?
               !!work = ratio*rhoi*evapi(ij)*aicen(i,j)
               !!fiso_evapn(i,j,k) = fiso_evapn(i,j,k)+work/dt

               sloss1 = c0
               sloss2 = c0
               if (dzssli(ij) > puny)               &
                  sloss1 = isoice(i,j,k,1)*         &
                     min(-evapi(ij),dzssli(ij))/dzssli(ij)
               if (isoice(i,j,k,1) >= sloss1) then
                  isoice(i,j,k,1) = isoice(i,j,k,1)-sloss1
               else
                  sloss1 = isoice(i,j,k,1)
                  isoice(i,j,k,1) = c0
               endif
               if (dzinti(ij) > puny)               &
                  sloss2 = isoice(i,j,k,2)*         &
                     max(-evapi(ij)-dzssli(ij),c0)/dzinti(ij)
               if (isoice(i,j,k,2) >= sloss2) then
                  isoice(i,j,k,2) = isoice(i,j,k,2)-sloss2
               else
                  sloss2 = isoice(i,j,k,2)
                  isoice(i,j,k,2) = c0
               endif
               fiso_evapn(i,j,k) = fiso_evapn(i,j,k)-(sloss1+sloss2)/dt
            enddo

            dzinti(ij) = dzinti(ij)+min(dzssli(ij)+evapi(ij),c0)
            dzssli(ij) = max(dzssli(ij)+evapi(ij),c0)
            if ( dzssli(ij) <= c0) then ! ssl goes away
               fiso_evapn(i,j,:) = fiso_evapn(i,j,:)-isoice(i,j,:,1)/dt !CRT
               isoice(i,j,:,1) = c0
               dzssli(ij) = c0
            endif
            if (dzinti(ij) <= c0) then  ! int goes away
               fiso_evapn(i,j,:) = fiso_evapn(i,j,:)-isoice(i,j,:,2)/dt !CRT
               isoice(i,j,:,2) = c0
               dzinti(ij) = c0
            endif
         endif
      enddo             ! icells

!     surface snow melt
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)

        if (melts(i,j) > c0) then
         do k=1,n_iso
          sloss1=c0
          sloss2=c0
          if (dzssl(ij) > puny)         &
           sloss1 = isosno(i,j,k,1)     &
                 *min(melts(i,j),dzssl(ij))/dzssl(ij)
          if (isosno(i,j,k,1) >= sloss1) then
             isosno(i,j,k,1) = isosno(i,j,k,1)-sloss1
          else
             sloss1 = isosno(i,j,k,1)
             isosno(i,j,k,1) = c0
          endif
          if (dzint(ij) > puny)             &
             sloss2=isosno(i,j,k,2)         &
                   *max(melts(i,j)-dzssl(ij),c0)/dzint(ij)
          if (isosno(i,j,k,2) >= sloss2) then
             isosno(i,j,k,2) = isosno(i,j,k,2)-sloss2
          else
             sloss2 = isosno(i,j,k,2)
             isosno(i,j,k,2) = c0
          endif
          fiso_ocnn(i,j,k)=fiso_ocnn(i,j,k)+(sloss1+sloss2)/dt
         enddo  ! n_iso

         dzint(ij)=dzint(ij)+min(dzssl(ij)-melts(i,j),c0)
         dzssl(ij)=max(dzssl(ij)-melts(i,j),c0)
         if ( dzssl(ij) <= c0) then ! ssl melts away
          fiso_ocnn(i,j,:)= fiso_ocnn(i,j,:)+isosno(i,j,:,1)/dt
          isosno(i,j,:,1) = c0
          dzssl(ij) = c0
         endif
         if (dzint(ij) <= c0) then  ! int melts away
          fiso_ocnn(i,j,:) = fiso_ocnn(i,j,:)+isosno(i,j,:,2)/dt
          isosno(i,j,:,2) = c0
          dzint(ij) = c0
         endif
        endif
      enddo         ! icells

!     surface ice melt
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)

        if (meltt(i,j) > c0) then
         do k=1,n_iso
          sloss1=c0
          sloss2=c0
          if (dzssli(ij) > puny)    &
           sloss1=isoice(i,j,k,1)   &
                 *min(meltt(i,j),dzssli(ij))/dzssli(ij)
          if (isoice(i,j,k,1) >= sloss1) then
             isoice(i,j,k,1) = isoice(i,j,k,1)-sloss1
          else
             sloss1 = isoice(i,j,k,1)
             isoice(i,j,k,1) = c0
          endif
          if (dzinti(ij) > puny)    &
           sloss2=isoice(i,j,k,2)   &
                 *max(meltt(i,j)-dzssli(ij),c0)/dzinti(ij)
          if (isoice(i,j,k,2) >= sloss2) then
             isoice(i,j,k,2) = isoice(i,j,k,2)-sloss2
          else
             sloss2 = isoice(i,j,k,2)
             isoice(i,j,k,2) = c0
          endif
          fiso_ocnn(i,j,k)=fiso_ocnn(i,j,k)+(sloss1+sloss2)/dt
         enddo

         dzinti(ij)=dzinti(ij)+min(dzssli(ij)-meltt(i,j),c0)
         dzssli(ij)=max(dzssli(ij)-meltt(i,j),c0)
         if (dzssli(ij) <= c0) then   ! ssl ice melts away
          fiso_ocnn(i,j,:) = fiso_ocnn(i,j,:)+isoice(i,j,:,1)/dt
          isoice(i,j,:,1) = c0
          dzssli(ij) = c0
         endif
         if (dzinti(ij) <= c0) then   ! int ice melts away
          fiso_ocnn(i,j,:) = fiso_ocnn(i,j,:)+isoice(i,j,:,2)/dt
          isoice(i,j,:,2) = c0
          dzinti(ij) = c0
         endif
        endif
       enddo        ! icells

!      basal ice melt.  Assume all isotopes lost in basal melt
       do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)

        if (meltb(i,j) > c0) then
         do k=1,n_iso
          sloss1=c0
          sloss2=c0
          if (dzssli(ij) > puny)  &
           sloss1=max(meltb(i,j)-dzinti(ij),c0)  &
                 *isoice(i,j,k,1)/dzssli(ij)
          if (isoice(i,j,k,1) >= sloss1) then
             isoice(i,j,k,1) = isoice(i,j,k,1)-sloss1
          else
             sloss1 = isoice(i,j,k,1)
             isoice(i,j,k,1) = c0
          endif
          if (dzinti(ij) > puny)  &
           sloss2=min(meltb(i,j),dzinti(ij))  &
                 *isoice(i,j,k,2)/dzinti(ij)
          if (isoice(i,j,k,2) >= sloss2) then
             isoice(i,j,k,2) = isoice(i,j,k,2)-sloss2
          else
             sloss2 = isoice(i,j,k,2)
             isoice(i,j,k,2) = c0
          endif
          fiso_ocnn(i,j,k)=fiso_ocnn(i,j,k)+(sloss1+sloss2)/dt
         enddo
 
         dzssli(ij) = dzssli(ij)+min(dzinti(ij)-meltb(i,j), c0)
         dzinti(ij) = max(dzinti(ij)-meltb(i,j), c0)           
         if (dzssli(ij) <= c0) then   ! ssl ice melts away
          fiso_ocnn(i,j,:) = fiso_ocnn(i,j,:)+isoice(i,j,:,1)/dt
          isoice(i,j,:,1) = c0
          dzssli(ij) = c0
         endif
         if (dzinti(ij) <= c0) then   ! int ice melts away
          fiso_ocnn(i,j,:) = fiso_ocnn(i,j,:)+isoice(i,j,:,2)/dt
          isoice(i,j,:,2) = c0
          dzinti(ij) = c0
         endif
        endif
       enddo        ! icells

!     snowfall and isotope deposition
       do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (fsnow(i,j) > c0) then
           isosno(i,j,:,1) = isosno(i,j,:,1)+               &
                             fiso_atm(i,j,:)*aicen(i,j)*dt
           dzssl(ij) = dzssl(ij)+fsnow(i,j)/rhos*dt
         endif
       enddo

!     snoice formation
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)

        if (dhs_snoice(ij) > c0) then
         do k=1,n_iso
          sloss1=c0
          sloss2=c0
          if (dzint(ij) > puny)                         &
           sloss2 = min(dhs_snoice(ij),dzint(ij))      &
                  *isosno(i,j,k,2)/dzint(ij)
          if (isosno(i,j,k,2) >= sloss2) then
             isosno(i,j,k,2) = isosno(i,j,k,2)-sloss2
          else
             sloss2 = isosno(i,j,k,2)
             isosno(i,j,k,2) = c0
          endif
          if (dzssl(ij) > puny)                         &
           sloss1 = max(dhs_snoice(ij)-dzint(ij),c0)   &
                  *isosno(i,j,k,1)/dzssl(ij)
          if (isosno(i,j,k,1) >= sloss1) then
             isosno(i,j,k,1) = isosno(i,j,k,1)-sloss1
          else
             sloss1 = isosno(i,j,k,1)
             isosno(i,j,k,1) = c0
          endif
          isoice(i,j,k,1) = isoice(i,j,k,1)+sloss1+sloss2
         enddo

         dzssl(ij)=dzssl(ij)-max(dhs_snoice(ij)-dzint(ij),c0)
         dzint(ij)=max(dzint(ij)-dhs_snoice(ij),c0)
         dzssli(ij)=dzssli(ij)+snoice(i,j)
         if ( dzssl(ij) <= c0) then ! ssl goes away
          fiso_ocnn(i,j,:)= fiso_ocnn(i,j,:)+isosno(i,j,:,1)/dt
          isosno(i,j,:,1) = c0
          dzssl(ij) = c0
         endif
         if (dzint(ij) <= c0) then  ! int goes away
          fiso_ocnn(i,j,:) = fiso_ocnn(i,j,:)+isosno(i,j,:,2)/dt
          isosno(i,j,:,2) = c0
          dzint(ij) = c0
         endif
        endif
      enddo         ! icells

!     redistribute isotope within vertical layers
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hslyr = hs(ij)/real(nslyr,kind=dbl_kind)
         hilyr = hi(ij)/real(nilyr,kind=dbl_kind)
         dzssl_new = min(hslyr/c2,hs_ssl)       ! new ssl for snow
         dzint_new = hs(ij)-dzssl_new
         dzssli_new = min(hilyr/c2,hi_ssl)      ! new ssl for ice
         dzinti_new = hi(ij)-dzssli_new

         do k=1,n_iso

          dznew=min(dzssl_new-dzssl(ij),c0)
          sloss1=c0
          if (dzssl(ij) > puny) &
           sloss1=dznew*isosno(i,j,k,1)/dzssl(ij) ! not neccesarily a loss term
          dznew=max(dzssl_new-dzssl(ij),c0)
          if (dzint(ij) > puny) &
           sloss1=sloss1+isosno(i,j,k,2)*dznew/dzint(ij) ! not really a loss term
          isosno(i,j,k,1) =isosno(i,j,k,1)+sloss1 
          isosno(i,j,k,2) =isosno(i,j,k,2)-sloss1

          sloss2=c0
          dznew=min(dzssli_new-dzssli(ij),c0)
          if (dzssli(ij) > puny) & 
           sloss2=dznew*isoice(i,j,k,1)/dzssli(ij)
          dznew=max(dzssli_new-dzssli(ij),c0)
          if (dzinti(ij) > puny) & 
           sloss2=sloss2+isoice(i,j,k,2)*dznew/dzinti(ij)
          isoice(i,j,k,1) = isoice(i,j,k,1)+sloss2
          isoice(i,j,k,2) = isoice(i,j,k,2)-sloss2

          isotot(ij,k)=isosno(i,j,k,2)+isosno(i,j,k,1)      &
           +isoice(i,j,k,2)+isoice(i,j,k,1)
          if ( (isotot(ij,k)-isotot0(ij,k))                 &
             - fiso_atm(i,j,k)*aicen(i,j)*dt                &
             - fiso_evapn(i,j,k)*dt                         &
             + fiso_ocnn(i,j,k)*dt > 1e-3) then
            write(nu_diag,*) 'isotope tracer:      ',k
            write(nu_diag,*) 'isotot-isotot0     ',isotot(ij,k)-isotot0(ij,k) 
            write(nu_diag,*) 'fiso_atm-fiso_ocnn          ',fiso_atm(i,j,k)*aicen(i,j)*dt &
                                                         + fiso_evapn(i,j,k)*dt &
                                                         - fiso_ocnn(i,j,k)*dt
          endif

         enddo          ! n_iso
      enddo             ! icells

!     reload tracers
      do ij = 1, icells
       i = indxi(ij)
       j = indxj(ij)

       do k=1,n_iso
          ! Update tracers only when vsnon/vicen is large enough.
          ! Otherwise, they are unchanged from last timestep. 
          if (vsnon(i,j) > puny) then
             trcrn(i,j,nt_iso+(k-1)*4  ) = isosno(i,j,k,1)/vsnon(i,j)
             trcrn(i,j,nt_iso+(k-1)*4+1) = isosno(i,j,k,2)/vsnon(i,j)
          endif
          if (vicen(i,j) > puny) then
             trcrn(i,j,nt_iso+(k-1)*4+2) = isoice(i,j,k,1)/vicen(i,j)
             trcrn(i,j,nt_iso+(k-1)*4+3) = isoice(i,j,k,2)/vicen(i,j)
          endif

         !do n = 1,2
         ! limit the trcrn to be positive
         !  if (trcrn(i,j,nt_iso+(k-1)*4+n-1) < puny) then
         !     trcrn(i,j,nt_iso+(k-1)*4+n-1) = c0
         !     fiso_ocnn(i,j,k) = fiso_ocnn(i,j,k) + &
         !       trcrn(i,j,nt_iso+(k-1)*4+n-1)*vsnon(i,j)/dt
         !  endif
         !  if (trcrn(i,j,nt_iso+(k-1)*4+n+1) < puny) then
         !     trcrn(i,j,nt_iso+(k-1)*4+n+1) = c0
         !     fiso_ocnn(i,j,k) = fiso_ocnn(i,j,k) + &
         !       trcrn(i,j,nt_iso+(k-1)*4+n+1)*vicen(i,j)/dt
         !  endif
         !enddo

       enddo        ! n_iso

!     scale fiso_ocnn. It will be re-scaled by aicen latter in merge_fluxes
      if (aicen(i,j) > puny) then
         fiso_ocnn(i,j,:) = fiso_ocnn(i,j,:)/aicen(i,j)
         fiso_evapn(i,j,:) = fiso_evapn(i,j,:)/aicen(i,j)
      endif

      enddo             ! icells

      do ij = 1, icells
       i = indxi(ij)
       j = indxj(ij)
       if (trcrn(i,j,nt_iso) < -puny .or. trcrn(i,j,nt_iso+1) < -puny    &
       .or. trcrn(i,j,nt_iso+2) < -puny .or. trcrn(i,j,nt_iso+3) < -puny) then
           write(nu_diag,*) 'isotope negative in isotope code'
           write(nu_diag,*) 'INT neg in isotope my_task = ',&
                               my_task &
                               ,' printing point = ',n &
                               ,' i and j = ',i,j
           write(nu_diag,*) 'Int Neg iso snowssl    = ',isosno0(i,j,1,1)
           write(nu_diag,*) 'Int Neg iso new snowssl= ',isosno(i,j,1,1)
           write(nu_diag,*) 'Int Neg iso snowint    = ',isosno0(i,j,1,2)
           write(nu_diag,*) 'Int Neg iso new snowint= ',isosno(i,j,1,2)
           write(nu_diag,*) 'Int Neg iso icessl     = ',isoice0(i,j,1,1)
           write(nu_diag,*) 'Int Neg iso new ice ssl= ',isoice(i,j,1,1)
           write(nu_diag,*) 'Int Neg iso ice int    = ',isoice0(i,j,1,2)
           write(nu_diag,*) 'Int Neg iso new ice int= ',isoice(i,j,1,2)
           write(nu_diag,*) 'Int Neg iso vicen      = ',vice_old(i,j)
           write(nu_diag,*) 'Int Neg iso vsnon      = ',vsno_old(i,j)
           write(nu_diag,*) 'Int Neg iso aicen      = ',aicen(i,j)
           write(nu_diag,*) 'Int Neg iso new vicen  = ',vicen(i,j)
           write(nu_diag,*) 'Int Neg iso new vsnon  = ',vsnon(i,j)
           write(nu_diag,*) 'Int Neg iso melts      = ',melts(i,j)
           write(nu_diag,*) 'Int Neg iso meltt      = ',meltt(i,j)
           write(nu_diag,*) 'Int Neg iso meltb      = ',meltb(i,j)
           write(nu_diag,*) 'Int Neg iso congel     = ',congel(i,j)
           write(nu_diag,*) 'Int Neg iso snoice     = ',snoice(i,j)
           write(nu_diag,*) 'Int Neg iso evap sno   = ',evaps(ij)
           write(nu_diag,*) 'Int Neg iso evap ice   = ',evapi(ij)
           write(nu_diag,*) 'Int Neg iso fsnow      = ',fsnow(i,j)
           write(nu_diag,*) 'Int Neg iso fiso_atm   = ',fiso_atm(i,j,1)
           write(nu_diag,*) 'Int Neg iso fiso_ocnn  = ',fiso_ocnn(i,j,1)
          endif
      enddo

      end subroutine update_isotope



!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_iso - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_iso(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for restarting
!
! !REVISION HISTORY:
!
! authors Elizabeth Hunke, LANL (original version)
!         David Bailey, NCAR
!         Marika Holland, NCAR
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
         
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.iso.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_iso,filename,0)

      if (my_task == master_task) then
        write(nu_dump_iso) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do k = 1, n_iso
      do n = 1, ncat
       call ice_write(nu_dump_iso,0,trcrn(:,:,nt_iso  +(k-1)*4,n,:),'ruf8',diag)
       call ice_write(nu_dump_iso,0,trcrn(:,:,nt_iso+1+(k-1)*4,n,:),'ruf8',diag)
       call ice_write(nu_dump_iso,0,trcrn(:,:,nt_iso+2+(k-1)*4,n,:),'ruf8',diag)
       call ice_write(nu_dump_iso,0,trcrn(:,:,nt_iso+3+(k-1)*4,n,:),'ruf8',diag)
      enddo
      enddo

      if (my_task == master_task) close(nu_dump_iso)

      end subroutine write_restart_iso

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_iso - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_iso(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for an ice isotope restart
!
! !REVISION HISTORY:
!
! authors Elizabeth Hunke, LANL (original version)
!         David Bailey, NCAR
!         Marika Holland, NCAR
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag

      if (my_task == master_task) then
         ! reconstruct path/file
         if (present(filename_spec)) then
            filename = filename_spec
         else
            open(nu_rst_pointer,file=pointer_file)
            read(nu_rst_pointer,'(a)') filename0
            filename = trim(filename0)
            close(nu_rst_pointer)

            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice('isotope restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.iso', &
               string2(1:lenstr(string2))
         endif
      endif ! master_task

      call ice_open(nu_restart_iso,filename,0)

      if (my_task == master_task) then
        read(nu_restart_iso) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do k = 1, n_iso
      do n = 1, ncat
       call ice_read(nu_restart_iso,0,trcrn(:,:,nt_iso  +(k-1)*4,n,:),'ruf8',&
            diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call ice_read(nu_restart_iso,0,trcrn(:,:,nt_iso+1+(k-1)*4,n,:),'ruf8',&
            diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call ice_read(nu_restart_iso,0,trcrn(:,:,nt_iso+2+(k-1)*4,n,:),'ruf8',&
            diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call ice_read(nu_restart_iso,0,trcrn(:,:,nt_iso+3+(k-1)*4,n,:),'ruf8',&
            diag,field_type=field_type_scalar,field_loc=field_loc_center)
      enddo
      enddo

      if (my_task == master_task) close(nu_restart_iso)

      end subroutine read_restart_iso

!=======================================================================
      function isoice_alpha(growth_rate, sp, frac)

! !DESCRIPTION:
!
! calculate the fractionation coefficient for sea-ice formation
!
! !REVISION HISTORY:
!
! authors: Jiang Zhu, UW-Madison 
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::   &
         growth_rate                     ! sea-ice formation rate (m/s)
      character(*), intent(in) ::   &
         sp,frac                         ! species: H2_16O, H2_18O, HDO
                                         ! calculation methods:
                                         !  cfrac, constant fractionation
                                         !  nfrac, nonfractionation
                                         !  gfrac, growth-rate dependent
      real (kind=dbl_kind) ::   &
         isoice_alpha                    ! return fractionation

      if (frac == 'nfrac') isoice_alpha = c1
      if (sp == 'H2_16O')  isoice_alpha = c1

      ! Lehmann and Siegenthaler, 1991
      !--------------------------------------------------
      if (frac == 'cfrac' .and. sp == 'HDO')            &
         isoice_alpha = 1.02120_dbl_kind
      if (frac == 'cfrac' .and. sp == 'H2_18O')         &
         isoice_alpha = 1.00291_dbl_kind
         
      ! Eq.9, Toyota et al., 2013
      ! For HDO, 7.2852 = 0.2120/0.00291
      !--------------------------------------------------
      if (frac == 'gfrac' .and. sp == 'HDO')                        &
         isoice_alpha = c1+7.2852_dbl_kind*1.2280E-3_dbl_kind+      &
            0.7311E-3_dbl_kind*exp(-growth_rate/8.0100E8_dbl_kind)+ &
            0.8441E-3_dbl_kind*exp(-growth_rate/0.7800E6_dbl_kind)
      if (frac == 'gfrac' .and. sp == 'H2_18O')                     &
         isoice_alpha = c1+1.2280E-3_dbl_kind+                      &
            0.7311E-3_dbl_kind*exp(-growth_rate/8.0100E8_dbl_kind)+ &
            0.8441E-3_dbl_kind*exp(-growth_rate/0.7800E6_dbl_kind)
      return

      end function isoice_alpha

!=======================================================================

      end module ice_isotope

!=======================================================================
