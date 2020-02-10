;iclm4_IC_create.ncl
;
;Goal:  To create an isotopic Initial Conditions (IC) file for iCLM4.  The file is needed not only to allow the model to run
;properly, but it can also grealy help the spin-up of deep snow pack and soil water.
;
;Written by:  Jesse Nusbaumer <nusbaume@colorado.edu> - June, 2014.
;
;**************
;Load Libraries
;**************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;*************
;Begin Program
;*************
begin
;***********************
;Constants/Control panel
;***********************

Rstd = 1.0  ;Standard ratio used to convert from delta values to water masses.

;isotope delta values:
ifilo = "/glade/p/work/nusbaume/isotope_obs/GNIP/d18o.2x2.ltm.nc"
ifild = "/glade/p/work/nusbaume/isotope_obs/GNIP/hdo.2x2.ltm.nc"

;CLM4 IC file:

;data directory:
cdir = "/glade/p/work/nusbaume/inputdata/clm4/"

;0.5x.0.5:
cfil = "clmi.BCN.2000-01-01_0.47x0.63_gx1v6_simyr2000_c121030.nc"

;1x1:
;cfil = "clmi.BCN.2000-01-01_0.9x1.25_gx1v6_simyr2000_c100303.nc"

;2x2:
;cfil = "clmi.BCN.2000-01-01_1.9x2.5_gx1v6_simyr2000_c100309.nc"

;1x1 - 1850:
;cfil = "clmi.Bi1850CN.0863-01-01_0.9x1.25_gx1v6_simyr1850_c150724.nc"

;2x2 - 1850:
;cfil = "clmi.BCN.1949-01-01_1.9x2.5_gx1v6_simyr1850_c121113.nc"

;4x5 - 1850:
;cfil = "clmi.BCN.1949-01-01_4x5_gx3v7_simyr1850_c100322.nc"

;2x2 - 850:
;cfil = "b.e11.B1850C5CN.f19_g16.0850cntl.001.clm2.r.0850-01-01-00000.nc"

;********************
;Read in CLM4 IC file
;********************

cin = addfile(cdir+cfil,"r")   ;open netCDF file

h2osno = cin->H2OSNO      ;read in variables
h2oliq = cin->H2OSOI_LIQ
h2oice = cin->H2OSOI_ICE
h2ocan = cin->H2OCAN
wa     = cin->WA
wt     = cin->WT

;lat/lon data needed to apply delta values:
latc = cin->cols1d_lat ;latitudes per column
lonc = cin->cols1d_lon ;longitudes per column
latp = cin->pfts1d_lat ;latitudes per pft
lonp = cin->pfts1d_lon ;longitudes per pft

;**************************************
;Read in delta-values from obs or model
;**************************************

oin = addfile(ifilo,"r") ;open netCDF files
din = addfile(ifild,"r") 

d18o = oin->d18o         ;read in variables
dDo  = din->hdo

;lat/long data needed to apply delta values:
lati = oin->lat
loni = oin->lon

;For Nik's GNIP data set, one must calculate lat/lons,
;the actual values aren't present in the netCDF file.
;----------------------------------
latn = fspan(-90,90,dimsizes(lati))
lonn = fspan(0,360,dimsizes(loni))

delete(lati)
delete(loni)

lati = latn
loni = lonn
;--------------------------------

;Average over time (if needed):
d18 = dim_avg(d18o(lat|:,lon|:,time|:))
dD  = dim_avg(dDo(lat|:,lon|:,time|:))

;******************************
;Convert delta values to ratios
;******************************

R18 = ((d18/1000.)+1.0)*Rstd
RD  = ((dD/1000.)+1.0)*Rstd

;print checks:
print(max(RD))
print(min(RD))

print(max(R18))
print(min(R18))

;*********************************
;Find closest column/pft locations
;*********************************

;NOTE: Using the nearest coordinate method may result 
;in the loss of data, given that an individual point
;is used to represent the grid cell, instead of an
;interpolated average.  However, for now this is probably
;an ok alternative. -JN

;columns:
ilat_col = ind_nearest_coord(latc,lati,0)
ilon_col = ind_nearest_coord(lonc,loni,0)

;pfts:
ilat_pft = ind_nearest_coord(latp,lati,0)
ilon_pft = ind_nearest_coord(lonp,loni,0)

;*************************************
;Create new water mass/ratio variables
;*************************************

;1-D variables:

;snow water:
h2osno_h2o = h2osno ;copy dimension sizes
h2osno_o18 = h2osno 
h2osno_hdo = h2osno
do c = 0,dimsizes(h2osno)-1 ;loop over columns
  h2osno_o18(c) = h2osno(c)*R18(ilat_col(c),ilon_col(c))
  h2osno_hdo(c) = h2osno(c)*RD(ilat_col(c),ilon_col(c))
end do

;unconfined aquifer:
wa_h2o = wa
wa_o18 = wa
wa_hdo = wa

;total water storage:
wt_h2o = wt
wt_o18 = wt
wt_hdo = wt

;GNIP isotopic ratio:
Rgnip_h2o = wt
Rgnip_o18 = wt
Rgnip_hdo = wt

Rgnip_h2o(:)= 1.0 ;set all values to one

;Canopy water:
h2ocan_h2o = h2ocan
h2ocan_o18 = h2ocan
h2ocan_hdo = h2ocan

;Canopy isotopic vapor ratio:
Rcan_h2o   = h2ocan
Rcan_o18   = h2ocan
Rcan_hdo   = h2ocan

Rcan_h2o(:)= 1.0 ;set all values to one.

;Sun isotopic leaf water ratio:
Rsun_h2o   = h2ocan
Rsun_o18   = h2ocan
Rsun_hdo   = h2ocan

Rsun_h2o(:)= 1.0 ;set all values to one.

;Shade isotopic leaf water ratio:
Rsha_h2o   = h2ocan
Rsha_o18   = h2ocan
Rsha_hdo   = h2ocan

Rsha_h2o(:)= 1.0 ;set all values to one

;Isotopic ratio of sun leaf produced cellulose:
Rsncl_h2o   = h2ocan
Rsncl_o18   = h2ocan
Rsncl_hdo   = h2ocan

Rsncl_h2o(:) = 1.0 ;set all values to one

;Isotopic ratio of shade leaf produced cellulose:
Rshcl_h2o   = h2ocan
Rshcl_o18   = h2ocan
Rshcl_hdo   = h2ocan

Rshcl_h2o(:) = 1.0 ;set all values to one

;Isotopic ratio of sun cellulose * photosynthesis:
Rsnclp_h2o   = h2ocan
Rsnclp_o18   = h2ocan
Rsnclp_hdo   = h2ocan

Rsnclp_h2o(:) = 1.0 ;set all values to one

;Isotopic ratio of shade cellulose * photosynthesis:
Rshclp_h2o   = h2ocan
Rshclp_o18   = h2ocan
Rshclp_hdo   = h2ocan

Rshclp_h2o(:) = 1.0 ;set all values to one

;Peclet effect - included sunlit leaf ratios (diagnostic only):
Rsndg_h2o   = h2ocan
Rsndg_o18   = h2ocan
Rsndg_hdo   = h2ocan

Rsndg_h2o(:) = 1.0 ;set all values to one

;Peclet effect - included shaded leaf ratios (diagnostic only):
Rshdg_h2o   = h2ocan
Rshdg_o18   = h2ocan
Rshdg_hdo   = h2ocan

Rshdg_h2o(:) = 1.0 ;set all values to one

;steady-state sunlit leaf ratios:
Rsnss_h2o   = h2ocan
Rsnss_o18   = h2ocan
Rsnss_hdo   = h2ocan

Rsnss_h2o(:) = 1.0 ;set all values to one

;steady-state shaded leaf ratios:
Rshss_h2o   = h2ocan
Rshss_o18   = h2ocan
Rshss_hdo   = h2ocan

Rshss_h2o(:) = 1.0 ;set all values to one

;steady-state sunlit leaf ratios * photosynthesis:
Rsnssp_h2o   = h2ocan
Rsnssp_o18   = h2ocan
Rsnssp_hdo   = h2ocan

Rsnssp_h2o(:) = 1.0 ;set all values to one

;steady-state shaded leaf ratios * photosynthesis:
Rshssp_h2o   = h2ocan
Rshssp_o18   = h2ocan
Rshssp_hdo   = h2ocan

Rshssp_h2o(:) = 1.0 ;set all values to one

;Isotopic ratio of Xylem water:
Rxym_h2o    = h2ocan
Rxym_o18    = h2ocan
Rxym_hdo    = h2ocan

Rxym_h2o(:) = 1.0

;Isotopic canopy liquid water ratio:
Rcnwt_h2o   = h2ocan
Rcnwt_o18   = h2ocan
Rcnwt_hdo   = h2ocan

Rcnwt_h2o   = 1.0

;2-D variables:

;soil water:
h2oliq_h2o = h2oliq ;copy dimensions
h2oliq_o18 = h2oliq
h2oliq_hdo = h2oliq

;soil ice:
h2oice_h2o = h2oice
h2oice_o18 = h2oice
h2oice_hdo = h2oice

;*******************************************
;Calculate isotopic values for new variables
;*******************************************

do c = 0,dimsizes(h2osno)-1 ;loop over columns

  ;1-D variables:

  ;Snow water:
  h2osno_o18(c) = h2osno(c)*R18(ilat_col(c),ilon_col(c))
  h2osno_hdo(c) = h2osno(c)*RD(ilat_col(c),ilon_col(c))

  ;unconfined aquifer:
  ;Gvien that wa appears to be constanly reset in CLM4,
  ;it is easier just to set it equal to bulk water. - JN

;  wa_o18(c) = wa(c)*R18(ilat_col(c),ilon_col(c))
;  wa_hdo(c) = wa(c)*RD(ilat_col(c),ilon_col(c))

  ;total water storage:
  wt_o18(c) = wt(c)*R18(ilat_col(c),ilon_col(c))
  wt_hdo(c) = wt(c)*RD(ilat_col(c),ilon_col(c))

  ;GNIP ratios:
  Rgnip_o18(c) = R18(ilat_col(c),ilon_col(c))
  Rgnip_hdo(c) = RD(ilat_col(c),ilon_col(c)) 

  ;2-D variables:
  
  do k = 0,dimsizes(h2oliq(0,:))-1 ;loop over vertical levels:
    ;soil water:
    h2oliq_o18(c,k) = h2oliq(c,k)*R18(ilat_col(c),ilon_col(c))
    h2oliq_hdo(c,k) = h2oliq(c,k)*RD(ilat_col(c),ilon_col(c))

    ;soil ice:
    h2oice_o18(c,k) = h2oice(c,k)*R18(ilat_col(c),ilon_col(c))
    h2oice_hdo(c,k) = h2oice(c,k)*RD(ilat_col(c),ilon_col(c))
  end do

end do

do p = 0,dimsizes(h2ocan)-1 ;loop over PFTs

  ;Canopy water:
  h2ocan_o18(p) = h2ocan(p)*R18(ilat_pft(p),ilon_pft(p))
  h2ocan_hdo(p) = h2ocan(p)*RD(ilat_pft(p),ilon_pft(p))

  ;Canopy isotopic vapor ratio:
  Rcan_o18(p)   = R18(ilat_pft(p),ilon_pft(p))
  Rcan_hdo(p)   = RD(ilat_pft(p),ilon_pft(p))

  ;Sun isotopic leaf water ratio:
  Rsun_o18(p)   = R18(ilat_pft(p),ilon_pft(p))
  Rsun_hdo(p)   = RD(ilat_pft(p),ilon_pft(p))

  ;Shade isotopic leaf water ratio:
  Rsha_o18(p)   = R18(ilat_pft(p),ilon_pft(p))
  Rsha_hdo(p)   = RD(ilat_pft(p),ilon_pft(p))

  ;Isotopic ratio of sun leaf produced cellulose:
  Rsncl_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rsncl_hdo(p)  = RD(ilat_pft(p),ilon_pft(p))

  ;Isotopic ratio of shade leaf produced cellulose:
  Rshcl_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rshcl_hdo(p)  = RD(ilat_pft(p),ilon_pft(p)) 

  ;Isotopic ratio of sun cellulose * photosynthesis:
  Rsnclp_o18(p) = R18(ilat_pft(p),ilon_pft(p))
  Rsnclp_hdo(p) = RD(ilat_pft(p),ilon_pft(p))

  ;Isotopic ratio of shade cellulose * photosynthesis:
  Rshclp_o18(p) = R18(ilat_pft(p),ilon_pft(p))
  Rshclp_hdo(p) = RD(ilat_pft(p),ilon_pft(p))

  ;Peclet effect - included sunlit leaf ratios (diagnostic only):
  Rsndg_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rsndg_hdo(p)  = RD(ilat_pft(p),ilon_pft(p))

  ;Peclet effect - included shaded leaf ratios (diagnostic only):
  Rshdg_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rshdg_hdo(p)  = RD(ilat_pft(p),ilon_pft(p))   

  ;steady-state sunlit leaf ratios:
  Rsnss_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rsnss_hdo(p)  = RD(ilat_pft(p),ilon_pft(p))

  ;steady-state shaded leaf ratios:
  Rshss_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rshss_hdo(p)  = RD(ilat_pft(p),ilon_pft(p))

  ;steady-state sunlit leaf ratios * photosynthesis:
  Rsnssp_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rsnssp_hdo(p)  = RD(ilat_pft(p),ilon_pft(p))

  ;steady-state shaded leaf ratios * photosynthesis:
  Rshssp_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rshssp_hdo(p)  = RD(ilat_pft(p),ilon_pft(p))

  ;Isotopic ratio of Xylem water:
  Rxym_o18(p)   = R18(ilat_pft(p),ilon_pft(p))
  Rxym_hdo(p)   = RD(ilat_pft(p),ilon_pft(p))

  ;Isotopic canopy liquid water ratio:
  Rcnwt_o18(p)  = R18(ilat_pft(p),ilon_pft(p))
  Rcnwt_hdo(p)  = RD(ilat_pft(p),ilon_pft(p))

end do

;******************
;generate meta-data
;******************

;----------
;dimensions:
;----------

;1-D variables:

h2osno_h2o!0 = "column"
h2osno_o18!0 = "column"
h2osno_hdo!0 = "column"

wa_h2o!0     = "column"
wa_o18!0     = "column"
wa_hdo!0     = "column"

wt_h2o!0     = "column"
wt_o18!0     = "column"
wt_hdo!0     = "column"

Rgnip_h2o!0  = "column"
Rgnip_o18!0  = "column"
Rgnip_hdo!0  = "column"

h2ocan_h2o!0 = "pft"
h2ocan_o18!0 = "pft"
h2ocan_hdo!0 = "pft"

Rcan_h2o!0   = "pft"
Rcan_o18!0   = "pft"
Rcan_hdo!0   = "pft"

Rsun_h2o!0   = "pft"
Rsun_o18!0   = "pft"
Rsun_hdo!0   = "pft"

Rsha_h2o!0   = "pft"
Rsha_o18!0   = "pft"
Rsha_hdo!0   = "pft"

Rsncl_h2o!0   = "pft"
Rsncl_o18!0   = "pft"
Rsncl_hdo!0   = "pft"

Rshcl_h2o!0   = "pft"
Rshcl_o18!0   = "pft"
Rshcl_hdo!0   = "pft"

Rsnclp_h2o!0   = "pft"
Rsnclp_o18!0   = "pft"
Rsnclp_hdo!0   = "pft"

Rshclp_h2o!0   = "pft"
Rshclp_o18!0   = "pft"
Rshclp_hdo!0   = "pft"

Rsndg_h2o!0   = "pft"
Rsndg_o18!0   = "pft"
Rsndg_hdo!0   = "pft"

Rshdg_h2o!0   = "pft"
Rshdg_o18!0   = "pft"
Rshdg_hdo!0   = "pft"

Rsnss_h2o!0   = "pft"
Rsnss_o18!0   = "pft"
Rsnss_hdo!0   = "pft"

Rshss_h2o!0   = "pft"
Rshss_o18!0   = "pft"
Rshss_hdo!0   = "pft"

Rsnssp_h2o!0  = "pft"
Rsnssp_o18!0  = "pft"
Rsnssp_hdo!0  = "pft"

Rshssp_h2o!0  = "pft"
Rshssp_o18!0  = "pft"
Rshssp_hdo!0  = "pft"

Rxym_h2o!0    = "pft"
Rxym_o18!0    = "pft"
Rxym_hdo!0    = "pft"

Rcnwt_h2o!0   = "pft"
Rcnwt_o18!0   = "pft"
Rcnwt_hdo!0   = "pft"

;2-D variables:
h2oliq_h2o!0 = "column"
h2oliq_o18!0 = "column"
h2oliq_hdo!0 = "column"

h2oliq_h2o!1 = "levtot"
h2oliq_o18!1 = "levtot"
h2oliq_hdo!1 = "levtot"

h2oice_h2o!0 = "column"
h2oice_o18!0 = "column"
h2oice_hdo!0 = "column"

h2oice_h2o!1 = "levtot"
h2oice_o18!1 = "levtot"
h2oice_hdo!1 = "levtot"
;-----------

;----------
;Attributes:
;----------

;1-D variables:

h2osno_h2o@long_name = "H2OTR snow water"
h2osno_h2o@units     = "kg/m2"
h2osno_o18@long_name = "H218O snow water"
h2osno_o18@units     = "kg/m2"
h2osno_hdo@long_name = "HDO snow water"
h2osno_hdo@units     = "kg/m2"

wa_h2o@long_name     = "H2OTR water in the unconfined aquifer"
wa_h2o@units         = "mm"
wa_o18@long_name     = "H218O water in the unconfined aquifer"
wa_o18@units         = "mm"
wa_hdo@long_name     = "HDO water in the unconfined aquifer"
wa_hdo@units         = "mm"

wt_h2o@long_name     = "H2OTR total water storage"
wt_h2o@units         = "mm"
wt_o18@long_name     = "H218O total water storage"
wt_o18@units         = "mm"
wt_hdo@long_name     = "HDO total water storage"
wt_hdo@units         = "mm"

Rgnip_h2o@long_name  = "H2OTR GNIP isotope ratio"
Rgnip_h2o@units      = "unitless"
Rgnip_o18@long_name  = "H218O GNIP isotope ratio"
Rgnip_o18@units      = "unitless"
Rgnip_hdo@long_name  = "HDO GNIP isotope ratio"
Rgnip_hdo@units      = "unitless"  

h2ocan_h2o@long_name = "H2OTR canopy water"
h2ocan_h2o@units     = "kg/m2"
h2ocan_o18@long_name = "H218O canopy water"
h2ocan_o18@units     = "kg/m2"
h2ocan_hdo@long_name = "HDO canopy water"
h2ocan_hdo@units     = "kg/m2"

Rcan_h2o@long_name   = "H2OTR canopy vapor ratio"
Rcan_h2o@units       = "unitless"
Rcan_o18@long_name   = "H218O canopy vapor ratio"
Rcan_o18@units       = "unitless"
Rcan_hdo@long_name   = "HDO canopy vapor ratio"
Rcan_hdo@units       = "unitless"

Rsun_h2o@long_name   = "H2OTR sun leaf water ratio"
Rsun_h2o@units       = "unitless"
Rsun_o18@long_name   = "H218O sun leaf water ratio"
Rsun_o18@units       = "unitless"
Rsun_hdo@long_name   = "HDO sun leaf water ratio"
Rsun_hdo@units       = "unitless"

Rsha_h2o@long_name   = "H2OTR shade leaf water ratio"
Rsha_h2o@units       = "unitless"
Rsha_o18@long_name   = "H218O shade leaf water ratio"
Rsha_o18@units       = "unitless"
Rsha_hdo@long_name   = "HDO shade leaf water ratio"
Rsha_hdo@units       = "unitless"

Rsncl_h2o@long_name   = "H2OTR to H2O ratio of sunlit cellulose (incl Peclet)"
Rsncl_h2o@units       = "unitless"
Rsncl_o18@long_name   = "H218O to H2O ratio of sunlit cellulose (incl Peclet)"
Rsncl_o18@units       = "unitless"
Rsncl_hdo@long_name   = "HDO to H2O ratio of sunlit cellulose (incl Peclet)"
Rsncl_hdo@units       = "unitless"

Rshcl_h2o@long_name   = "H2OTR to H2O ratio of shaded cellulose (incl Peclet)"
Rshcl_h2o@units       = "unitless"
Rshcl_o18@long_name   = "H218O to H2O ratio of shaded cellulose (incl Peclet)"
Rshcl_o18@units       = "unitless"
Rshcl_hdo@long_name   = "HDO to H2O ratio of shaded cellulose (incl Peclet)"
Rshcl_hdo@units       = "unitless"

Rsnclp_h2o@long_name   = "H2OTR to H2O ratio of sunlit cellulose (incl Peclet), times sunlit photosynthesis"
Rsnclp_h2o@units       = "unitless"
Rsnclp_o18@long_name   = "H218O to H2O ratio of sunlit cellulose (incl Peclet), times sunlit photosynthesis"
Rsnclp_o18@units       = "unitless"
Rsnclp_hdo@long_name   = "HDO to H2O ratio of sunlit cellulose (incl Peclet), times sunlit photosynthesis"
Rsnclp_hdo@units       = "unitless"

Rshclp_h2o@long_name   = "H2OTR to H2O ratio of shaded cellulose (incl Peclet), times sunlit photosynthesis"
Rshclp_h2o@units       = "unitless"
Rshclp_o18@long_name   = "H218O to H2O ratio of shaded cellulose (incl Peclet), times sunlit photosynthesis"
Rshclp_o18@units       = "unitless"
Rshclp_hdo@long_name   = "HDO to H2O ratio of shaded cellulose (incl Peclet), times sunlit photosynthesis"
Rshclp_hdo@units       = "unitless"

Rsndg_h2o@long_name    = "H2OTR to H2O ratio of sunlit leaf water (incl Peclet)"
Rsndg_h2o@units        = "unitless"
Rsndg_o18@long_name    = "H218O to H2O ratio of sunlit leaf water (incl Peclet)"
Rsndg_o18@units        = "unitless"
Rsndg_hdo@long_name    = "HDO to H2O ratio of sunlit leaf water (incl Peclet)"
Rsndg_hdo@units        = "unitless"

Rshdg_h2o@long_name    = "H2OTR to H2O ratio of shaded leaf water (incl Peclet)"
Rshdg_h2o@units        = "unitless"
Rshdg_o18@long_name    = "H218O to H2O ratio of shaded leaf water (incl Peclet)"
Rshdg_o18@units        = "unitless"
Rshdg_hdo@long_name    = "HDO to H2O ratio of shaded leaf water (incl Peclet)"
Rshdg_hdo@units        = "unitless"

Rsnss_h2o@long_name    = "H2OTR to H2O ratio of sunlit leaf water (steady state)"
Rsnss_h2o@units        = "unitless"
Rsnss_o18@long_name    = "H218O to H2O ratio of sunlit leaf water (steady state)"
Rsnss_o18@units        = "unitless"
Rsnss_hdo@long_name    = "HDO to H2O ratio of sunlit leaf water (steady state)"
Rsnss_hdo@units        = "unitless"

Rshss_h2o@long_name    = "H2OTR to H2O ratio of shaded leaf water (steady state)"
Rshss_h2o@units        = "unitless"
Rshss_o18@long_name    = "H218O to H2O ratio of shaded leaf water (steady state)"
Rshss_o18@units        = "unitless"
Rshss_hdo@long_name    = "HDO to H2O ratio of shaded leaf water (steady state)"
Rshss_hdo@units        = "unitless"

Rsnssp_h2o@long_name    = "H2OTR to H2O ratio of sunlit leaf water (incl Peclet), times sunlit photosynthesis"
Rsnssp_h2o@units        = "unitless"
Rsnssp_o18@long_name    = "H218O to H2O ratio of sunlit leaf water (incl Peclet), times sunlit photosynthesis"
Rsnssp_o18@units        = "unitless"
Rsnssp_hdo@long_name    = "HDO to H2O ratio of sunlit leaf water (incl Peclet), times sunlit photosynthesis"
Rsnssp_hdo@units        = "unitless"

Rshssp_h2o@long_name    = "H2OTR to H2O ratio of shaded leaf water (incl Peclet), times shaded photosynthesis"
Rshssp_h2o@units        = "unitless"
Rshssp_o18@long_name    = "H218O to H2O ratio of shaded leaf water (incl Peclet), times shaded photosynthesis"
Rshssp_o18@units        = "unitless"
Rshssp_hdo@long_name    = "HDO to H2O ratio of shaded leaf water (incl Peclet), times shaded photosynthesis"
Rshssp_hdo@units        = "unitless"

Rxym_h2o@long_name    = "H2OTR to H2O ratio of xylem water"
Rxym_h2o@units        = "unitless"
Rxym_o18@long_name    = "H218O to H2O ratio of xylem water"
Rxym_o18@units        = "unitless"
Rxym_hdo@long_name    = "HDO to H2O ratio of xylem water"
Rxym_hdo@units        = "unitless"

Rcnwt_h2o@long_name    = "H2OTR to H2O ratio of Canopy liquid water"
Rcnwt_h2o@units        = "unitless"
Rcnwt_o18@long_name    = "H218O to H2O ratio of Canopy liquid water"
Rcnwt_o18@units        = "unitless"
Rcnwt_hdo@long_name    = "HDO to H2O ratio of Canopy liquid water"
Rcnwt_hdo@units        = "unitless"

;2-D variables:

h2oliq_h2o@long_name = "H2OTR liquid water"
h2oliq_h2o@units     = "kg/m2"
h2oliq_o18@long_name = "H218O liquid water"
h2oliq_o18@units     = "kg/m2"
h2oliq_hdo@long_name = "HDO liquid water"
h2oliq_hdo@units     = "kg/m2"

h2oice_h2o@long_name = "H2OTR ice lens"
h2oice_h2o@units     = "kg/m2"
h2oice_o18@long_name = "H218O ice lens"
h2oice_o18@units     = "kg/m2"
h2oice_hdo@long_name = "HDO ice lens"
h2oice_hdo@units     = "kg/m2"

;********************
;Generate netCDF file
;********************

;temporary file:
tfil = "/glade/p/work/nusbaume/iso_clm4_IC.nc"

system("rm -f "+tfil)      ;remove old file

out = addfile("/glade/p/work/nusbaume/iso_clm4_IC.nc","c") ;create new netCDF file

;Snow water:
out->H2OSNO_H2OTR = h2osno_h2o
out->H2OSNO_HDO   = h2osno_hdo
out->H2OSNO_H218O = h2osno_o18

;Unconfined aquifer:
out->WA_H2OTR = wa_h2o
out->WA_HDO   = wa_hdo
out->WA_H218O = wa_o18

;Total water storage:
out->WT_H2OTR = wt_h2o
out->WT_HDO   = wt_hdo
out->WT_H218O = wt_o18

;GNIP isotope ratios:
out->RGNIP_H2OTR = Rgnip_h2o
out->RGNIP_HDO   = Rgnip_hdo
out->RGNIP_H218O = Rgnip_o18

;Canopy water:
out->H2OCAN_H2OTR = h2ocan_h2o
out->H2OCAN_HDO   = h2ocan_hdo
out->H2OCAN_H218O = h2ocan_o18

;Isotopic ratio of Canopy vapor:
out->Rcanopy_H2OTR = Rcan_h2o
out->Rcanopy_HDO   = Rcan_hdo
out->Rcanopy_H218O = Rcan_o18

;Isotopic ratio of sun leaf water:
out->Rleafsun_H2OTR = Rsun_h2o
out->Rleafsun_HDO   = Rsun_hdo
out->Rleafsun_H218O = Rsun_o18

;Isotopic ratio of shade leaf water:
out->Rleafsha_H2OTR = Rsha_h2o
out->Rleafsha_HDO   = Rsha_hdo
out->Rleafsha_H218O = Rsha_o18

;Soil water:
out->H2OSOI_LIQ_H2OTR = h2oliq_h2o
out->H2OSOI_LIQ_HDO   = h2oliq_hdo
out->H2OSOI_LIQ_H218O = h2oliq_o18

;Soil ice:
out->H2OSOI_ICE_H2OTR = h2oice_h2o
out->H2OSOI_ICE_HDO   = h2oice_hdo
out->H2OSOI_ICE_H218O = h2oice_o18

;Isotopic ratio of sun leaf produced cellulose:
out->Rcellsun_H2OTR = Rsncl_h2o
out->Rcellsun_HDO   = Rsncl_hdo
out->Rcellsun_H218O = Rsncl_o18

;Isotopic ratio of shade leaf produced cellulose:
out->Rcellsha_H2OTR = Rshcl_h2o
out->Rcellsha_HDO   = Rshcl_hdo
out->Rcellsha_H218O = Rshcl_o18

;Isotopic ratio of sun cellulose * photosynthesis:
out->Rcellsunpsn_H2OTR = Rsnclp_h2o
out->Rcellsunpsn_HDO   = Rsnclp_hdo
out->Rcellsunpsn_H218O = Rsnclp_o18

;Isotopic ratio of shade cellulose * photosynthesis:
out->Rcellshapsn_H2OTR = Rshclp_h2o
out->Rcellshapsn_HDO   = Rshclp_hdo
out->Rcellshapsn_H218O = Rshclp_o18

;Peclet effect - included sunlit leaf ratios (diagnostic only):
out->Rleafsundiag_H2OTR = Rsndg_h2o
out->Rleafsundiag_HDO   = Rsndg_hdo
out->Rleafsundiag_H218O = Rsndg_o18

;Peclet effect - included shaded leaf ratios (diagnostic only):
out->Rleafshadiag_H2OTR = Rshdg_h2o
out->Rleafshadiag_HDO   = Rshdg_hdo
out->Rleafshadiag_H218O = Rshdg_o18

;steady-state sunlit leaf ratios:
out->Rleafsunstdy_H2OTR = Rsnss_h2o
out->Rleafsunstdy_HDO   = Rsnss_hdo
out->Rleafsunstdy_H218O = Rsnss_o18

;steady-state shaded leaf ratios:
out->Rleafshastdy_H2OTR = Rshss_h2o
out->Rleafshastdy_HDO   = Rshss_hdo
out->Rleafshastdy_H218O = Rshss_o18

;steady-state sunlit leaf ratios * photosynthesis:
out->Rleafsunpsn_H2OTR = Rsnssp_h2o
out->Rleafsunpsn_HDO   = Rsnssp_hdo
out->Rleafsunpsn_H218O = Rsnssp_o18

;steady-state sunlit leaf ratios * photosynthesis:
out->Rleafshapsn_H2OTR = Rshssp_h2o
out->Rleafshapsn_HDO   = Rshssp_hdo
out->Rleafshapsn_H218O = Rshssp_o18

;Isotopic ratio of Xylem water:
out->Rxylem_H2OTR = Rxym_h2o
out->Rxylem_HDO   = Rxym_hdo
out->Rxylem_H218O = Rxym_o18

;Isotopic ratio of Canopy liquid water:
out->Rcanopywater_H2OTR = Rcnwt_h2o
out->Rcanopywater_HDO   = Rcnwt_hdo
out->Rcanopywater_H218O = Rcnwt_o18

;*****************************
;Combine netCDF files together
;*****************************

;file output directory:
odir = "/glade/p/work/nusbaume/inputdata/isotopic_inputdata/clm4/GNIP_based/"

;Copy over base CLM4 file:
system("cp "+cdir+cfil+" "+odir+cfil)

;Ensure file can be written:
system("chmod u+w "+odir+cfil)

;Add isotopic variables to file:
system("ncks -A "+tfil+" "+odir+cfil)

;***********
;End program
;***********
end
