!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module nd_mod

!BOP
! !MODULE: nd_mod
!
! !DESCRIPTION:Module for Nd 
!  The concentration of Nd143 and Nd144 are treated as tracers.
!  Unit for Nd143 and Nd144 is pmol/m^3
!
! !REVISION HISTORY:
!  SVN:$Id: nd_mod.F90 26603 2011-01-28 23:09:02Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use blocks, only: nx_block, ny_block, block
   use domain_size, only: max_blocks_clinic, km, nx_global, ny_global
   use domain, only: nblocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use prognostic, only: tracer_field
   use kinds_mod
   use constants, only: c0, c1, p5, char_blank, delim_fmt, field_type_scalar
   use io, only: data_set
   use io_types, only: stdout, nml_in, nml_filename, datafile, io_dim,  &
       io_field_desc, rec_type_dbl, construct_file, construct_io_dim,   &
       construct_io_field, destroy_file, destroy_io_field
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
   use passive_tracer_tools, only: ind_name_pair, tracer_read, &
       rest_read_tracer_block, file_read_tracer_block, read_field
   
   
   
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: nd_tracer_cnt,        &
             nd_init,              &
             nd_set_interior,      &
             nd_tavg_forcing

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracer
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
       nd_tracer_cnt = 2

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
       nd143_ind = 1,      &     ! nd143 index
       nd144_ind = 2          ! nd144 index

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(nd_tracer_cnt) :: &
       ind_name_table = (/ ind_name_pair(nd143_ind, 'ND143'),ind_name_pair(nd144_ind, 'ND144') /)

!-----------------------------------------------------------------------
!  tavg ids for non-standard tavg variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      
       tavg_ND143_D,                &
       tavg_ND143_P,                &
       tavg_ND144_D,                &
       tavg_ND144_P,                &
       tavg_ND143_SOURCE,           &
       tavg_ND143_SINK,             &
       tavg_ND144_SOURCE,           &
       tavg_ND144_SINK,             &
       tavg_ND143_RESET_TEND,        &
       tavg_ND144_RESET_TEND
     
      

!-----------------------------------------------------------------------
!  mask the continental margin area
!-----------------------------------------------------------------------
   logical(log_kind),dimension(:,:,:,:),allocatable :: &
       MARGIN_MASK

!-----------------------------------------------------------------------
!  Surface dust deposition surf_dust: read from dust.nc unit: gm^-2s^-1
!-----------------------------------------------------------------------
   real(r8),dimension(:,:,:),allocatable :: &
       surf_dust,           & !surface dust deposition 
       ir_dust,             &
       river_con,           &
       ir_river,            &
       river_source,        &
       river_source143,     &
       river_source144,     &
       ir_margin,           &
       dust_source,        &
       dust_source144,  &
       dust_source143
 
    
   real(r8),dimension(:,:,:,:),allocatable :: &
       calcite_r,                             & ! calcite mass concentration ratio
       poc_r,                                 &  ! poc mass concentration ratio	
       opal_r,                                &  ! opal mass concentration ratio
       dust_r,                                & ! dust mass concentration ratio
       margin_source,                         &
       margin_source144,                      &
       margin_source143,                      &
       prodk_cndk,                            &  ! denominator in equation 14(rempfer2011)
       nd143_d,                               &
       nd143_p,                               &
       nd144_d,                               &
       nd144_p,                               &
       nd143_source,                          &
       nd143_sink,                            &
       nd144_source,                          &
       nd144_sink

  
   real(r8)          :: &
       nd_w 


!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: nd_init
! !INTERFACE:
subroutine nd_init(init_ts_file_fmt, read_restart_filename, &
                      tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize nd tracer module. This involves setting metadata, reading
!  the module namelist and setting initial conditions.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use broadcast, only: broadcast_scalar
   use prognostic, only: curtime, oldtime
   use grid, only: KMT, n_topo_smooth, fill_points, dz
   use time_management, only: seconds_in_year

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
       init_ts_file_fmt,       &   ! format (bin or nc) for input file
       read_restart_filename       ! file name for restart file

    
! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(nd_tracer_cnt), intent(inout) :: &
       tracer_d_module   ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,nd_tracer_cnt,3,max_blocks_clinic), &
       intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
       errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'nd_mod:Nd_init'

   character(char_len) :: &
       init_nd_option,           & ! option for initialization of nd
       init_nd_init_file,        & ! filename for option 'file'
       init_nd_init_file_fmt,    & ! file format for option 'file'
  
       prescribed_filename,       & ! fielname for dust deposition
       surf_dust_varname,        & ! variale name for dust deposition
       dust_ir_varname,          &
       river_con_varname,        & ! variable name for Nd concentration in river runoff
       river_ir_varname,         &
       margin_ir_varname,        & ! variable name for river ir
       calcite_r_varname,        & ! variable name for calcite R
       poc_r_varname,            & ! variable name for poc R
       opal_r_varname,           & ! variable name for opal R
       dust_r_varname              ! variable name for dust R

   real(r8)          :: &
       cnd_dust,                 & ! global concentration of Nd in dust
       beta_dust,                & ! fraction of Nd released
       nd_mass,                  & ! Nd mass g/pmol
       a_boundary_tot,           & ! total area for boundary
       f_boundary,               &  ! total source from boundary	
       k_poc,                    &
       k_calcite,                &
       k_opal,                   &
       k_dust
  
   logical(log_kind) :: &
       lnml_found             ! Was nd_nml found ?

   integer(int_kind) :: &
       n,                   & ! index for looping over trlsacers
       k,                   & ! index for looping over depth levels
       nx,                  & ! index for looping over x direction
       ny,                  & ! index for looping over y direction
       iblock,              & ! index for looping over blocks
       nml_error              ! namelist i/o error flag

!     l,                   & ! index for looping over time levels
   type(tracer_read), dimension(nd_tracer_cnt) :: &
       nd_init_ext        ! namelist variable for initializing tracers

   namelist /nd_nml/ &
       init_nd_option, init_nd_init_file, nd_init_ext, &
       init_nd_init_file_fmt, prescribed_filename, surf_dust_varname, &
        dust_ir_varname, river_con_varname,&
       river_ir_varname, margin_ir_varname,&
       calcite_r_varname,poc_r_varname,       &
       opal_r_varname, dust_r_varname,        &
       cnd_dust, beta_dust, nd_mass,a_boundary_tot, f_boundary, k_poc,         &
       k_calcite,k_opal,k_dust,nd_w


   character (char_len) ::  &
       nd_restart_filename  ! modified file name for restart file
      
  
!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success
  
   tracer_d_module(nd143_ind)%short_name = 'ND143'
   tracer_d_module(nd143_ind)%long_name  = 'ND143'
   tracer_d_module(nd143_ind)%units      = 'pmol/m^3'
   tracer_d_module(nd143_ind)%tend_units = 'pmol/(m^3s)'
!   tracer_d_module(nd143_ind)%flux_units = 'cm years/s'


   tracer_d_module(nd144_ind)%short_name = 'ND144'
   tracer_d_module(nd144_ind)%long_name  = 'ND144'
   tracer_d_module(nd144_ind)%units      = 'pmol/m^3'
   tracer_d_module(nd144_ind)%tend_units = 'pmol/(m^3s)'
!  tracer_d_module(nd144_ind)%flux_units = 'cm years/s'
!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_nd_option = 'unknown'
   init_nd_init_file = 'unknown'
   init_nd_init_file_fmt = 'bin'

   do n = 1,nd_tracer_cnt
       nd_init_ext(n)%mod_varname  = 'unknown'
       nd_init_ext(n)%filename     = 'unknown'
       nd_init_ext(n)%file_varname = 'unknown'
       nd_init_ext(n)%scale_factor = c1
       nd_init_ext(n)%default_val  = c0
       nd_init_ext(n)%file_fmt     = 'nc'
   end do
   
   prescribed_filename = 'unknown'
   surf_dust_varname = 'unknown'
   dust_ir_varname = 'unknown'
   river_con_varname = 'unknown'
   river_ir_varname = 'unknown'
   margin_ir_varname = 'unknown'
   calcite_r_varname = 'unknown'
   poc_r_varname =     'unknown'
   opal_r_varname  =   'unknown'
   dust_r_varname  =   'unknown'
  
   cnd_dust = c0
   beta_dust = c0
   nd_mass = c0
   a_boundary_tot = c0
   f_boundary = c0
   k_poc = c0
   k_calcite = c0
   k_opal = c0
   k_dust = c0
   nd_w = c0
  
   if (my_task == master_task) then
       open (nml_in, file=nml_filename, status='old',iostat=nml_error)
       if (nml_error /= 0) then  
           nml_error = -1
       else
           nml_error =  1      
       endif
       do while (nml_error > 0)
           read(nml_in, nml=nd_nml,iostat=nml_error)
       end do
       if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
       call document(subname, 'nd_nml not found')
       call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_nd_option , master_task)
   call broadcast_scalar(init_nd_init_file, master_task)
   call broadcast_scalar(init_nd_init_file_fmt, master_task)

   do n = 1,nd_tracer_cnt
      call broadcast_scalar(nd_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(nd_init_ext(n)%filename, master_task)
      call broadcast_scalar(nd_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(nd_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(nd_init_ext(n)%default_val, master_task)
      call broadcast_scalar(nd_init_ext(n)%file_fmt, master_task)
   end do


   call broadcast_scalar(prescribed_filename, master_task)
   call broadcast_scalar(surf_dust_varname, master_task)
   call broadcast_scalar(dust_ir_varname, master_task)
   call broadcast_scalar(river_con_varname, master_task)
   call broadcast_scalar(river_ir_varname, master_task)
   call broadcast_scalar(margin_ir_varname, master_task)
   call broadcast_scalar(calcite_r_varname, master_task)
   call broadcast_scalar(poc_r_varname, master_task)
   call broadcast_scalar(opal_r_varname, master_task)
   call broadcast_scalar(dust_r_varname, master_task)
   call broadcast_scalar(cnd_dust, master_task)
   call broadcast_scalar(beta_dust, master_task)
   call broadcast_scalar(nd_mass, master_task)
   call broadcast_scalar(a_boundary_tot, master_task)
   call broadcast_scalar(f_boundary, master_task)
   call broadcast_scalar(k_poc, master_task)
   call broadcast_scalar(k_calcite, master_task)
   call broadcast_scalar(k_opal, master_task)
   call broadcast_scalar(k_dust, master_task)
   call broadcast_scalar(nd_w, master_task)
   
!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   select case (init_nd_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d  isotope ratio set to all zeros' 
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
       
   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )
      nd_restart_filename = char_blank
      if (init_nd_init_file == 'same_as_TS') then
        if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read Nd from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
        endif
        nd_restart_filename = read_restart_filename
        init_nd_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file
        nd_restart_filename = trim(init_nd_init_file)

      endif
      call rest_read_tracer_block(init_nd_init_file_fmt, &
                                  nd_restart_filename,   &
                                  tracer_d_module,         &
                                  TRACER_MODULE)

   case ('file')
      call document(subname, 'Nd being read from separate file')

      call file_read_tracer_block(init_nd_init_file_fmt, &
                                  init_nd_init_file,     &
                                  tracer_d_module,         &
                                  ind_name_table,          &
                                  nd_init_ext,         &
                                  TRACER_MODULE)
 
      if (n_topo_smooth > 0) then
        do n = 1,nd_tracer_cnt
         do k=1,km
            call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'nd_init: error in fill_points for Nd')
               return
            endif
          end do
         end do
         
      endif
    case default
      call document(subname, 'init_nd_option = ', init_nd_option)
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
      do n = 1,nd_tracer_cnt
         do k = 1,km
            where (k > KMT(:,:,iblock))
                 TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
                 TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do
   enddo

!-----------------------------------------------------------------------

   call define_tavg_field(tavg_ND143_RESET_TEND, 'ND143_RESET_TEND',2,  &
                          long_name='surface reset tendency of ND', &
                          units='years/s', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   call define_tavg_field(tavg_ND144_RESET_TEND, 'ND144_RESET_TEND',2,  &
                          long_name='surface reset tendency of ND', &
                          units='years/s', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   
   
   
  
    
   call define_tavg_field(tavg_ND143_D, 'ND143_D',3,&
                          long_name='ND143 disolved',&
                          units='pmol/m^3', grid_loc='3111')
     
     
   call define_tavg_field(tavg_ND143_P, 'ND143_P',3,&
                          long_name='ND143 particle related',&
                          units='pmol/m^3', grid_loc='3111')
     
     
   call define_tavg_field(tavg_ND144_D, 'ND144_D',3,&
                          long_name='ND144 disolved',&
                          units='pmol/m^3', grid_loc='3111')
     
   call define_tavg_field(tavg_ND144_P, 'ND144_P',3,&
                          long_name='ND144 particle related',&
                           units='pmol/m^3', grid_loc='3111')
     
  
     
   call define_tavg_field(tavg_ND143_SOURCE, 'ND143_SOURCE',3,&
                          long_name='source term: dust and boundary',&
                          grid_loc='3111',coordinates='TLONG TLAT z_t time')
  
   call define_tavg_field(tavg_ND143_SINK, 'ND143_SINK',3,&
                          long_name='sink term',&
                          grid_loc='3111',coordinates='TLONG TLAT z_t time')
     
   call define_tavg_field(tavg_ND144_SOURCE, 'ND144_SOURCE',3,&
                          long_name='source term: dust and boundary',&
                          grid_loc='3111',coordinates='TLONG TLAT z_t time')
  
   call define_tavg_field(tavg_ND144_SINK, 'ND144_SINK',3,&
                          long_name='sink term',&
                          grid_loc='3111',coordinates='TLONG TLAT z_t time')
    
   write(stdout,delim_fmt)
   write(stdout,*) ' successful define tavg_SURFACE_DUST' 
   write(stdout,delim_fmt)

!-----------------------------------------------------------------------
!  allocate and initialize MARGIN_MASK (true for margin points)
!-----------------------------------------------------------------------
   allocate( MARGIN_MASK(nx_block,ny_block,km,max_blocks_clinic))

   MARGIN_MASK(:,:,:,:) = .false.

   do k = 1,51
       where (KMT .eq. k)
            MARGIN_MASK(:,:,k,:) = .true.
       end where
   end do


!-----------------------------------------------------------------------
!  read in surface dust deposition
!-----------------------------------------------------------------------

   allocate(surf_dust(nx_block,ny_block,max_blocks_clinic))
   surf_dust = c0
   call read_field('nc',prescribed_filename,surf_dust_varname,surf_dust)

   allocate(ir_dust(nx_block,ny_block,max_blocks_clinic))
   ir_dust = c0
   call read_field('nc',prescribed_filename,dust_ir_varname,ir_dust)

   allocate(river_con(nx_block,ny_block,max_blocks_clinic))
   river_con = c0
   call read_field('nc',prescribed_filename,river_con_varname,river_con)

   allocate(ir_river(nx_block,ny_block,max_blocks_clinic))
   ir_river = c0
   call read_field('nc',prescribed_filename,river_ir_varname,ir_river)

   allocate(river_source(nx_block,ny_block,max_blocks_clinic))
   river_source = c0

   allocate(river_source143(nx_block,ny_block,max_blocks_clinic))
   river_source143 = c0
   allocate(river_source144(nx_block,ny_block,max_blocks_clinic))
   river_source144 = c0

   allocate(ir_margin(nx_block,ny_block,max_blocks_clinic))
   ir_margin = c0
   call read_field('nc',prescribed_filename,margin_ir_varname,ir_margin)




!-----------------------------------------------------------------------
!  read in particle concentration ratio 
!-----------------------------------------------------------------------

   allocate(calcite_r(nx_block,ny_block,km,max_blocks_clinic))
   calcite_r = c0
   call read_field_3D('nc',prescribed_filename,calcite_r_varname,calcite_r)

   allocate(poc_r(nx_block,ny_block,km,max_blocks_clinic))
   poc_r = c0
   call read_field_3D('nc',prescribed_filename,poc_r_varname,poc_r)

   allocate(opal_r(nx_block,ny_block,km,max_blocks_clinic))
   opal_r = c0
   call read_field_3D('nc',prescribed_filename,opal_r_varname,opal_r)

   allocate(dust_r(nx_block,ny_block,km,max_blocks_clinic))
   dust_r = c0
   call read_field_3D('nc',prescribed_filename,dust_r_varname,dust_r)

!-----------------------------------------------------------------------
!  calculate prodk_cndk
!-----------------------------------------------------------------------
   allocate(prodk_cndk(nx_block,ny_block,km,max_blocks_clinic))



!-----------------------------------------------------------------------
!  calculate dust source term for total Nd, use ir_dust to calculate 
!  dust source for Nd143 and Nd144
!-----------------------------------------------------------------------

   allocate(dust_source(nx_block,ny_block,max_blocks_clinic))
   allocate(dust_source144(nx_block,ny_block,max_blocks_clinic))
   allocate(dust_source143(nx_block,ny_block,max_blocks_clinic))

   dust_source = surf_dust*cnd_dust*beta_dust/(dz(1)*0.01_r8*nd_mass)


!-----------------------------------------------------------------------
!  calculate margin source term for total Nd, use ir_dust to calculate 
!  dust source for Nd143 and Nd144
!-----------------------------------------------------------------------
   allocate(margin_source(nx_block,ny_block,km,max_blocks_clinic))
   allocate(margin_source144(nx_block,ny_block,km,max_blocks_clinic))
   allocate(margin_source143(nx_block,ny_block,km,max_blocks_clinic))


   do iblock = 1,nblocks_clinic
     do nx = 1,nx_block
       do ny = 1,ny_block 
            dust_source144(nx,ny,iblock) = 0.36_r8*dust_source(nx,ny,iblock)/(ir_dust(nx,ny,iblock)+c1)
            dust_source143(nx,ny,iblock) = 0.36_r8*dust_source(nx,ny,iblock)-dust_source144(nx,ny,iblock)
            do k = 1,km
                margin_source(nx,ny,k,iblock) = f_boundary/(a_boundary_tot*dz(k)*0.01_r8*seconds_in_year*nd_mass)
                margin_source144(nx,ny,k,iblock) = 0.36_r8*margin_source(nx,ny,k,iblock)/(ir_margin(nx,ny,iblock)+c1)
                margin_source143(nx,ny,k,iblock) = 0.36_r8*margin_source(nx,ny,k,iblock)-margin_source144(nx,ny,k,iblock)
                prodk_cndk(nx,ny,k,iblock) = c1 + k_poc*poc_r(nx,ny,k,iblock) +k_calcite*calcite_r(nx,ny,k,iblock) + &
                                             k_opal*opal_r(nx,ny,k,iblock)+ k_dust*dust_r(nx,ny,k,iblock)
            end do
      end do
     end do
   end do

   where(MARGIN_MASK)
   elsewhere
      margin_source = c0
      margin_source144 = c0
      margin_source143 = c0
   end where


!-----------------------------------------------------------------------
!  allocate space for nd143_d,nd143_p,nd144_d,nd144_p
!-----------------------------------------------------------------------
allocate(nd143_d(nx_block,ny_block,km,max_blocks_clinic))
allocate(nd143_p(nx_block,ny_block,km,max_blocks_clinic))
allocate(nd144_d(nx_block,ny_block,km,max_blocks_clinic))
allocate(nd144_p(nx_block,ny_block,km,max_blocks_clinic))

allocate(nd143_source(nx_block,ny_block,km,max_blocks_clinic))
allocate(nd143_sink(nx_block,ny_block,km,max_blocks_clinic))
allocate(nd144_source(nx_block,ny_block,km,max_blocks_clinic))
allocate(nd144_sink(nx_block,ny_block,km,max_blocks_clinic))





!EOC

 end subroutine nd_init

!***********************************************************************
!BOP
! !IROUTINE: nd_set_interior
! !INTERFACE:

 
subroutine nd_set_interior(k, TRACER_MODULE_OLD, TRACER_MODULE_CUR, &
                           DTRACER_MODULE,this_block)

! !DESCRIPTION:
!  set interior source/sink term for nd isotope ratio tracer
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use time_management, only: seconds_in_year
   use grid, only : dz,KMT,REGION_MASK
   
   use forcing_fields, only : ROFF_F
! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: &
      k                   ! vertical level index
      
   real (r8), dimension(:,:,:,:), intent(in) :: &
      TRACER_MODULE_OLD, &! old tracer values
      TRACER_MODULE_CUR   ! current tracer values
      
   type (block), intent(in) :: &
      this_block          ! block info for the current block

! !OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,nd_tracer_cnt), intent(out) :: &
      DTRACER_MODULE      ! computed source/sink term
    
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   
      
   integer (int_kind)       :: &
      bid  ,        &               ! local_block id
      nx,           &
      ny,           &
      ntracer
   
   logical(log_kind),dimension(nx_block,ny_block) :: &
      mask
  
   real(r8), dimension(nx_block,ny_block) :: &
      p143_k_upper,        &
      p143_k,              &
      p143_k_lower,        &
      p144_k_upper,        &
      p144_k,              &
      p144_k_lower,        &
      p143_top,            &
      p143_bot,            &
      p144_top,            &
      p144_bot,            &
      nd143_cur,           &
      nd144_cur,           &
      nd143_upper,         &
      nd144_upper,         &
      nd143_lower,         &
      nd144_lower
   
!EOP
!BOC

    bid = this_block%local_id
    
    river_source(:,:,bid) = ROFF_F(:,:,bid)* river_con(:,:,bid)/(dz(1)*0.01_r8)
    
!-----------------------------------------------------------------------
    DTRACER_MODULE = c0
!-----------------------------------------------------------------------
!  apply dust source and river source
!-----------------------------------------------------------------------
    if (k==1) then
        river_source(:,:,bid) = 0.3_r8*ROFF_F(:,:,bid)* river_con(:,:,bid)/(dz(1)*0.01_r8)
        river_source144(:,:,bid) = 0.36_r8*river_source(:,:,bid)/(ir_river(:,:,bid)+c1)
        river_source143(:,:,bid) = 0.36_r8*river_source(:,:,bid)-river_source144(:,:,bid)
        DTRACER_MODULE(:,:,nd143_ind) = DTRACER_MODULE(:,:,nd143_ind) + dust_source143(:,:,bid) +river_source143(:,:,bid)
        DTRACER_MODULE(:,:,nd144_ind) = DTRACER_MODULE(:,:,nd144_ind) + dust_source144(:,:,bid) +river_source144(:,:,bid)
    endif
!-----------------------------------------------------------------------
!  apply boundary source
!-----------------------------------------------------------------------   
    mask= MARGIN_MASK(:,:,k,bid)
  
    where(mask)
        DTRACER_MODULE(:,:,nd143_ind) = DTRACER_MODULE(:,:,nd143_ind) + margin_source143(:,:,k,bid)
        DTRACER_MODULE(:,:,nd144_ind) = DTRACER_MODULE(:,:,nd144_ind) + margin_source144(:,:,k,bid) 
    end where
    
    
!-----------------------------------------------------------------------
!  sink term
!-----------------------------------------------------------------------  

!  calculate nd143_d,nd143_p,nd144_d,nd144p and p143_top, p143_bot,p144_top,p144_bot
    p143_top = c0
    p143_bot = c0
    p144_top = c0
    p144_bot = c0
    nd143_cur = p5*(TRACER_MODULE_OLD(:,:,k,nd143_ind)+TRACER_MODULE_CUR(:,:,k,nd143_ind))
    nd144_cur = p5*(TRACER_MODULE_OLD(:,:,k,nd144_ind)+TRACER_MODULE_CUR(:,:,k,nd144_ind))

    where(KMT(:,:,bid)>=k)
        nd143_d(:,:,k,bid) = nd143_cur/prodk_cndk(:,:,k,bid)
        nd143_p(:,:,k,bid) = nd143_cur - nd143_d(:,:,k,bid)
    
        nd144_d(:,:,k,bid) = nd144_cur/prodk_cndk(:,:,k,bid)
        nd144_p(:,:,k,bid) = nd144_cur - nd144_d(:,:,k,bid)
    end where  
    
    
    
    where(nd143_d(:,:,k,bid).lt.c0)
        nd143_d(:,:,k,bid) = c0
    end where
    
    where(nd143_p(:,:,k,bid).lt.c0)
        nd143_p(:,:,k,bid) = c0
    end where
    
    where(nd144_d(:,:,k,bid).lt.c0)
        nd144_d(:,:,k,bid) = c0
    end where
    
    where(nd144_p(:,:,k,bid).lt.c0)
        nd144_p(:,:,k,bid) = c0
    end where
                  
    if(k==1) then
        p143_k = nd143_p(:,:,k,bid)
        nd143_lower =  p5*(TRACER_MODULE_OLD(:,:,k+1,nd143_ind)+ TRACER_MODULE_CUR(:,:,k+1,nd143_ind))
        nd144_lower =  p5*(TRACER_MODULE_OLD(:,:,k+1,nd144_ind)+ TRACER_MODULE_CUR(:,:,k+1,nd144_ind))              
        p143_k_lower = nd143_lower*(1-c1/prodk_cndk(:,:,k+1,bid))
    
        p143_bot     = p5*(p143_k + p143_k_lower)                      
        
        p144_k       = nd144_p(:,:,k,bid)        
        p144_k_lower = nd144_lower*(1-c1/prodk_cndk(:,:,k+1,bid))
    
        p144_bot     = p5*(p144_k + p144_k_lower) 
                        
    
    
    else if(k >1) then
       where((KMT(:,:,bid))>k)
    
            nd143_upper  = p5*(TRACER_MODULE_OLD(:,:,k-1,nd143_ind)+ TRACER_MODULE_CUR(:,:,k-1,nd143_ind))
            nd143_lower =  p5*(TRACER_MODULE_OLD(:,:,k+1,nd143_ind)+ TRACER_MODULE_CUR(:,:,k+1,nd143_ind))
     
            p143_k_upper = nd143_upper*(1-c1/prodk_cndk(:,:,k-1,bid))                         
   
            p143_k       = nd143_p(:,:,k,bid)
    
            p143_k_lower = nd143_lower*(1-c1/prodk_cndk(:,:,k+1,bid))
    
            p143_top     = p5*(p143_k + p143_k_upper) 
            p143_bot     = p5*(p143_k + p143_k_lower) 
    
            nd144_upper  = p5*(TRACER_MODULE_OLD(:,:,k-1,nd144_ind)+ TRACER_MODULE_CUR(:,:,k-1,nd144_ind))
            nd144_lower =  p5*(TRACER_MODULE_OLD(:,:,k+1,nd144_ind)+ TRACER_MODULE_CUR(:,:,k+1,nd144_ind))
    
            p144_k_upper = nd144_upper*(1-c1/prodk_cndk(:,:,k-1,bid))                       
        
            p144_k       = nd144_p(:,:,k,bid)
    
            p144_k_lower = nd144_lower*(1-c1/prodk_cndk(:,:,k+1,bid))
    
            p144_top     = p5*(p144_k + p144_k_upper) 
            p144_bot     = p5*(p144_k + p144_k_lower)
    
    
       end where
       where((KMT(:,:,bid)).eq.k)
            nd143_upper  = p5*(TRACER_MODULE_OLD(:,:,k-1,nd143_ind)+ TRACER_MODULE_CUR(:,:,k-1,nd143_ind))
            p143_k_upper = nd143_upper*(1-c1/prodk_cndk(:,:,k-1,bid))
            p143_k       = nd143_p(:,:,k,bid)
            p143_k_lower = c0
            p143_top     = p5*p143_k_upper
            p143_bot     = p5*p143_k
    
            nd144_upper  = p5*(TRACER_MODULE_OLD(:,:,k-1,nd144_ind)+ TRACER_MODULE_CUR(:,:,k-1,nd144_ind))
            p144_k_upper = nd144_upper*(1-c1/prodk_cndk(:,:,k-1,bid)) 
            p144_k       = nd144_p(:,:,k,bid)
            p144_k_lower = c0
    
            p144_top     = p5*p144_k_upper
            p144_bot     = p5*p144_k
    
    
       end where
    

    end if
    
    nd143_source(:,:,k,bid) =  DTRACER_MODULE(:,:,nd143_ind)
    nd144_source(:,:,k,bid) =  DTRACER_MODULE(:,:,nd143_ind)
    
    nd143_sink(:,:,k,bid) = nd_w*(p143_bot - p143_top)/dz(k)
    nd144_sink(:,:,k,bid) = nd_w*(p144_bot - p144_top)/dz(k)
    
    DTRACER_MODULE(:,:,nd143_ind) = DTRACER_MODULE(:,:,nd143_ind) - nd143_sink(:,:,k,bid)
    DTRACER_MODULE(:,:,nd144_ind) = DTRACER_MODULE(:,:,nd144_ind) - nd144_sink(:,:,k,bid)
 
    where(REGION_MASK(:,:,bid).eq.-14 .or. REGION_MASK(:,:,bid).eq.-13 .or. REGION_MASK(:,:,bid).eq.-12 .or.    &
           REGION_MASK(:,:,bid).eq.-5 .or. REGION_MASK(:,:,bid).eq.4 .or. REGION_MASK(:,:,bid).eq.11) 
              DTRACER_MODULE(:,:,nd143_ind) = c0
              DTRACER_MODULE(:,:,nd144_ind) = c0
    end where
      
    
!-----------------------------------------------------------------------
!EOC

 end subroutine nd_set_interior

!***********************************************************************

!***********************************************************************
!BOP
! !IROUTINE: nd_tavg_forcing
! !INTERFACE:

 subroutine nd_tavg_forcing
 


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
      iblock,   &      ! block loop index
      k

!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1, nblocks_clinic
       do k = 1,km
          call accumulate_tavg_field(nd143_d(:,:,k,iblock),tavg_ND143_D,iblock,k)
          call accumulate_tavg_field(nd143_p(:,:,k,iblock),tavg_ND143_P,iblock,k)
          call accumulate_tavg_field(nd144_d(:,:,k,iblock),tavg_ND144_D,iblock,k)
          call accumulate_tavg_field(nd144_p(:,:,k,iblock),tavg_ND144_P,iblock,k)
          call accumulate_tavg_field(nd143_source(:,:,k,iblock),tavg_ND143_SOURCE,iblock,k)
          call accumulate_tavg_field(nd143_sink(:,:,k,iblock),tavg_ND143_SINK,iblock,k)
          call accumulate_tavg_field(nd144_source(:,:,k,iblock),tavg_ND144_SOURCE,iblock,k)
          call accumulate_tavg_field(nd144_sink(:,:,k,iblock),tavg_ND144_SINK,iblock,k)
   
      end do  
   
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine nd_tavg_forcing

!***********************************************************************

!***********************************************************************
!BOP
! !IROUTINE: read_field_3D
! !INTERFACE:

 subroutine read_field_3D(fmt, filename, fieldname, FIELD, record_length)

! !DESCRIPTION:
!  read 3D field from a file
!  Assumes the field is (nx_global,ny_global), cell centered, and scalar.
!  The length of the 3rd dimension is determined by the dimension of FIELD.
!  For binary files, the default external precision is double precision.
!  This can be overridden by passing the desired precision into record_length.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      fmt,                 & ! format (bin or nc)
      filename,            & ! file to read from
      fieldname              ! field to be read

   integer(int_kind), intent(in), optional :: &
      record_length          ! record length type for binary files

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(:,:,:,:), intent(inout), target :: &
      FIELD                  ! field to be read in

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'passive_tracer_tools:read_field_3D'

   integer(int_kind) :: &
      record_length_loc    ! record length type for binary files

   type (io_field_desc) :: &
      FIELD_DESC           ! IO field descriptors for FIELD

   type (datafile) :: &
      restart_file         ! io file descriptor

   type (io_dim) :: &
      i_dim, j_dim, k_dim  ! dimension descriptors

!-----------------------------------------------------------------------

   call document(subname, 'reading ' /&
                       &/ trim(fieldname) /&
                       &/ ' from ' /&
                       &/ trim(filename))

   if (present(record_length)) then
      record_length_loc = record_length
   else
      record_length_loc = rec_type_dbl
   endif

   restart_file =                                     &
      construct_file(fmt,                             &
                     full_name=trim(filename),        &
                     record_length=record_length_loc, &
                     recl_words=nx_global*ny_global)

   call data_set(restart_file, 'open_read')

   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)
   k_dim = construct_io_dim('k', size(FIELD,3))

   FIELD_DESC =                                       &
      construct_io_field(trim(fieldname),             &
                         dim1=i_dim,                  &
                         dim2=j_dim,                  &
                         dim3=k_dim,                  &
                         grid_loc ='3111',            &
                         d3d_array = FIELD)

   call data_set (restart_file, 'define', FIELD_DESC)

   call data_set (restart_file, 'read', FIELD_DESC)

   call destroy_io_field (FIELD_DESC)

   call data_set (restart_file, 'close')

   call destroy_file (restart_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine read_field_3D


end module nd_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
