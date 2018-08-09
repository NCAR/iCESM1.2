module POP_CplIndices
  
  use seq_flds_mod
  use mct_mod

  implicit none

  SAVE
  public                               ! By default make data private

  ! ocn -> drv

  integer :: index_o2x_So_t      
  integer :: index_o2x_So_u
  integer :: index_o2x_So_v
  integer :: index_o2x_So_s
  integer :: index_o2x_So_dhdx
  integer :: index_o2x_So_dhdy
  integer :: index_o2x_Fioo_q
  integer :: index_o2x_Faoo_fco2_ocn
  integer :: index_o2x_Faoo_fdms_ocn
  integer :: index_o2x_So_roce_16O
  integer :: index_o2x_So_roce_18O
  integer :: index_o2x_So_roce_HDO

  ! drv -> ocn

  integer :: index_x2o_Si_ifrac        ! fractional ice wrt ocean
  integer :: index_x2o_So_duu10n       ! 10m wind speed squared           (m^2/s^2)
  integer :: index_x2o_Sa_pslv         ! sea-level pressure               (Pa)
  integer :: index_x2o_Sa_co2prog      ! bottom atm level prognostic CO2
  integer :: index_x2o_Sa_co2diag      ! bottom atm level diagnostic CO2
  integer :: index_x2o_Foxx_taux       ! zonal wind stress (taux)         (W/m2   )
  integer :: index_x2o_Foxx_tauy       ! meridonal wind stress (tauy)     (W/m2   )
  integer :: index_x2o_Foxx_swnet      ! net short-wave heat flux         (W/m2   )
  integer :: index_x2o_Foxx_sen        ! sensible heat flux               (W/m2   )
  integer :: index_x2o_Foxx_lat        
  integer :: index_x2o_Foxx_lwup       ! longwave radiation (up)          (W/m2   )
  integer :: index_x2o_Faxa_lwdn       ! longwave radiation (down)        (W/m2   )
  integer :: index_x2o_Fioi_melth      ! heat flux from snow & ice melt   (W/m2   )
  integer :: index_x2o_Fioi_meltw      ! snow melt flux                   (kg/m2/s)
  integer :: index_x2o_Fioi_salt       ! salt                             (kg(salt)/m2/s)
  integer :: index_x2o_Foxx_evap       ! evaporation flux                 (kg/m2/s)
  integer :: index_x2o_Faxa_prec         
  integer :: index_x2o_Faxa_snow       ! water flux due to snow           (kg/m2/s)
  integer :: index_x2o_Faxa_rain       ! water flux due to rain           (kg/m2/s)
  integer :: index_x2o_Faxa_bcphidry   ! flux: Black   Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_bcphodry   ! flux: Black   Carbon hydrophobic dry deposition
  integer :: index_x2o_Faxa_bcphiwet   ! flux: Black   Carbon hydrophilic wet deposition
  integer :: index_x2o_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2o_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_x2o_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition
  integer :: index_x2o_Forr_roff       ! river runoff flux                (kg/m2/s)
  integer :: index_x2o_Forr_ioff       ! ice runoff flux                  (kg/m2/s)

  integer :: index_x2o_Faxa_prec_16O         
  integer :: index_x2o_Faxa_snow_16O       ! H216O water flux due to snow           (kg/m2/s)
  integer :: index_x2o_Faxa_rain_16O       !  "    water flux due to rain           (kg/m2/s)
  integer :: index_x2o_Faxa_prec_18O         
  integer :: index_x2o_Faxa_snow_18O       ! H218O water flux due to snow           (kg/m2/s)
  integer :: index_x2o_Faxa_rain_18O       !   "   water flux due to rain           (kg/m2/s)
  integer :: index_x2o_Faxa_prec_HDO         
  integer :: index_x2o_Faxa_snow_HDO       ! HDO   water flux due to snow           (kg/m2/s)
  integer :: index_x2o_Faxa_rain_HDO       !  "    water flux due to rain           (kg/m2/s)
  integer :: index_x2o_Foxx_evap_16O       ! H216O evaporation flux                 (kg/m2/s)
  integer :: index_x2o_Foxx_evap_18O       ! H218O evaporation flux                 (kg/m2/s)
  integer :: index_x2o_Foxx_evap_HDO       ! H218O evaporation flux                 (kg/m2/s)
  integer :: index_x2o_Forr_roff_16O       ! H216O liq river runoff flux            (kg/m2/s)
  integer :: index_x2o_Forr_ioff_16O       ! H216O ice river runoff flux            (kg/m2/s)
  integer :: index_x2o_Forr_roff_18O       ! H218O liq river runoff flux            (kg/m2/s)
  integer :: index_x2o_Forr_ioff_18O       ! H218O ice river runoff flux            (kg/m2/s)
  integer :: index_x2o_Forr_roff_HDO       ! HDO   liq river runoff flux            (kg/m2/s)
  integer :: index_x2o_Forr_ioff_HDO       ! HDO   ice river runoff flux            (kg/m2/s)
  integer :: index_x2o_Fioi_meltw_16O      ! snow melt flux                   (kg/m2/s)
  integer :: index_x2o_Fioi_meltw_18O      ! snow melt flux                   (kg/m2/s)
  integer :: index_x2o_Fioi_meltw_HDO      ! snow melt flux                   (kg/m2/s)

contains

  subroutine POP_CplIndicesSet( )

    type(mct_aVect) :: o2x      ! temporary
    type(mct_aVect) :: x2o      ! temporary

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
    call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

    index_o2x_So_t          = mct_avect_indexra(o2x,'So_t')
    index_o2x_So_u          = mct_avect_indexra(o2x,'So_u')
    index_o2x_So_v          = mct_avect_indexra(o2x,'So_v')
    index_o2x_So_s          = mct_avect_indexra(o2x,'So_s')
    index_o2x_So_dhdx       = mct_avect_indexra(o2x,'So_dhdx')
    index_o2x_So_dhdy       = mct_avect_indexra(o2x,'So_dhdy')
    index_o2x_Fioo_q        = mct_avect_indexra(o2x,'Fioo_q')
    index_o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2_ocn',perrWith='quiet')
    index_o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')

    index_o2x_So_roce_16O = mct_avect_indexra(o2x,'So_roce_16O',perrWith='quiet')
    index_o2x_So_roce_18O = mct_avect_indexra(o2x,'So_roce_18O',perrWith='quiet')
    index_o2x_So_roce_HDO = mct_avect_indexra(o2x,'So_roce_HDO',perrWith='quiet')

    index_x2o_Si_ifrac      = mct_avect_indexra(x2o,'Si_ifrac')
    index_x2o_Sa_pslv       = mct_avect_indexra(x2o,'Sa_pslv')
    index_x2o_So_duu10n     = mct_avect_indexra(x2o,'So_duu10n')

    index_x2o_Foxx_tauy     = mct_avect_indexra(x2o,'Foxx_tauy')
    index_x2o_Foxx_taux     = mct_avect_indexra(x2o,'Foxx_taux')
    index_x2o_Foxx_swnet    = mct_avect_indexra(x2o,'Foxx_swnet')
    index_x2o_Foxx_lat      = mct_avect_indexra(x2o,'Foxx_lat')
    index_x2o_Foxx_sen      = mct_avect_indexra(x2o,'Foxx_sen')
    index_x2o_Foxx_lwup     = mct_avect_indexra(x2o,'Foxx_lwup')
    index_x2o_Faxa_lwdn     = mct_avect_indexra(x2o,'Faxa_lwdn')
    index_x2o_Fioi_melth    = mct_avect_indexra(x2o,'Fioi_melth')   
    index_x2o_Fioi_meltw    = mct_avect_indexra(x2o,'Fioi_meltw')
    index_x2o_Fioi_salt     = mct_avect_indexra(x2o,'Fioi_salt')   
    index_x2o_Faxa_prec     = mct_avect_indexra(x2o,'Faxa_prec')   
    index_x2o_Faxa_snow     = mct_avect_indexra(x2o,'Faxa_snow')   
    index_x2o_Faxa_rain     = mct_avect_indexra(x2o,'Faxa_rain')   
    index_x2o_Foxx_evap     = mct_avect_indexra(x2o,'Foxx_evap')
    index_x2o_Forr_roff     = mct_avect_indexra(x2o,'Forr_roff')
    index_x2o_Forr_ioff     = mct_avect_indexra(x2o,'Forr_ioff')
    index_x2o_Faxa_bcphidry = mct_avect_indexra(x2o,'Faxa_bcphidry')
    index_x2o_Faxa_bcphodry = mct_avect_indexra(x2o,'Faxa_bcphodry')
    index_x2o_Faxa_bcphiwet = mct_avect_indexra(x2o,'Faxa_bcphiwet')
    index_x2o_Faxa_ocphidry = mct_avect_indexra(x2o,'Faxa_ocphidry')
    index_x2o_Faxa_ocphodry = mct_avect_indexra(x2o,'Faxa_ocphodry')
    index_x2o_Faxa_ocphiwet = mct_avect_indexra(x2o,'Faxa_ocphiwet')
    index_x2o_Faxa_dstdry1  = mct_avect_indexra(x2o,'Faxa_dstdry1')
    index_x2o_Faxa_dstdry2  = mct_avect_indexra(x2o,'Faxa_dstdry2')
    index_x2o_Faxa_dstdry3  = mct_avect_indexra(x2o,'Faxa_dstdry3')
    index_x2o_Faxa_dstdry4  = mct_avect_indexra(x2o,'Faxa_dstdry4')
    index_x2o_Faxa_dstwet1  = mct_avect_indexra(x2o,'Faxa_dstwet1')
    index_x2o_Faxa_dstwet2  = mct_avect_indexra(x2o,'Faxa_dstwet2')
    index_x2o_Faxa_dstwet3  = mct_avect_indexra(x2o,'Faxa_dstwet3')
    index_x2o_Faxa_dstwet4  = mct_avect_indexra(x2o,'Faxa_dstwet4')
    index_x2o_Sa_co2prog    = mct_avect_indexra(x2o,'Sa_co2prog',perrWith='quiet')
    index_x2o_Sa_co2diag    = mct_avect_indexra(x2o,'Sa_co2diag',perrWith='quiet')
    index_x2o_Faxa_prec_16O = mct_avect_indexra(x2o,'Faxa_prec_16O',perrWith='quiet')   
    index_x2o_Faxa_snow_16O = mct_avect_indexra(x2o,'Faxa_snow_16O',perrWith='quiet')   
    index_x2o_Faxa_rain_16O = mct_avect_indexra(x2o,'Faxa_rain_16O',perrWith='quiet')   
    index_x2o_Faxa_prec_18O = mct_avect_indexra(x2o,'Faxa_prec_18O',perrWith='quiet')   
    index_x2o_Faxa_snow_18O = mct_avect_indexra(x2o,'Faxa_snow_18O',perrWith='quiet')   
    index_x2o_Faxa_rain_18O = mct_avect_indexra(x2o,'Faxa_rain_18O',perrWith='quiet')   
    index_x2o_Faxa_prec_HDO = mct_avect_indexra(x2o,'Faxa_prec_HDO',perrWith='quiet')   
    index_x2o_Faxa_snow_HDO = mct_avect_indexra(x2o,'Faxa_snow_HDO',perrWith='quiet')   
    index_x2o_Faxa_rain_HDO = mct_avect_indexra(x2o,'Faxa_rain_HDO',perrWith='quiet')   
    index_x2o_Foxx_evap_16O = mct_avect_indexra(x2o,'Foxx_evap_16O',perrWith='quiet')
    index_x2o_Foxx_evap_18O = mct_avect_indexra(x2o,'Foxx_evap_18O',perrWith='quiet')
    index_x2o_Foxx_evap_HDO = mct_avect_indexra(x2o,'Foxx_evap_HDO',perrWith='quiet')
    index_x2o_Forr_roff_16O = mct_avect_indexra(x2o,'Forr_roff_16O',perrWith='quiet')
    index_x2o_Forr_ioff_16O = mct_avect_indexra(x2o,'Forr_ioff_16O',perrWith='quiet')
    index_x2o_Forr_roff_18O = mct_avect_indexra(x2o,'Forr_roff_18O',perrWith='quiet')
    index_x2o_Forr_ioff_18O = mct_avect_indexra(x2o,'Forr_ioff_18O',perrWith='quiet')
    index_x2o_Forr_roff_HDO = mct_avect_indexra(x2o,'Forr_roff_HDO',perrWith='quiet')
    index_x2o_Forr_ioff_HDO = mct_avect_indexra(x2o,'Forr_ioff_HDO',perrWith='quiet')
    index_x2o_Fioi_meltw_16O= mct_avect_indexra(x2o,'Fioi_meltw_16O')
    index_x2o_Fioi_meltw_18O= mct_avect_indexra(x2o,'Fioi_meltw_18O')
    index_x2o_Fioi_meltw_HDO= mct_avect_indexra(x2o,'Fioi_meltw_HDO')

    call mct_aVect_clean(x2o)
    call mct_aVect_clean(o2x)

  end subroutine POP_CplIndicesSet

end module POP_CplIndices
