#!/bin/csh -f

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------
@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number  ($my_stream)
   exit 5
endif

@ s1 = 1   # use base-model stream 1
   
cat >! $CASEROOT/Buildconf/pop2conf/abio_dic_dic14_tavg_contents << EOF
$s1  ABIO_DIC
$s1  ABIO_DIC14 
$s1  ABIO_D14Cocn
$s1  FvICE_ABIO_DIC
$s1  FvPER_ABIO_DIC
$s1  FvICE_ABIO_DIC14
$s1  FvPER_ABIO_DIC14
$s1  FG_ABIO_DIC14
$s1  FG_ABIO_DIC
$s1  ABIO_CO2_IFRAC
$s1  ABIO_CO2_XKW
$s1  ABIO_CO2_ATM_PRESS
$s1  ABIO_CO2_SCHMIDT
$s1  ABIO_CO2_PV
$s1  ABIO_pCO2SURF
$s1  ABIO_DCO2STAR
$s1  ABIO_CO2STAR
$s1  ABIO_DpCO2
$s1  ABIO_pCO2
$s1  ABIO_D14Catm
$s1  ABIO_ALK_SURF
$s1  ABIO_PH_SURF
$s1  ABIO_DIC_SURF
$s1  ABIO_DIC14_SURF
$s1  STF_ABIO_DIC
$s1  STF_ABIO_DIC14
$s1  J_ABIO_DIC
$s1  J_ABIO_DIC14
$s1  Jint_ABIO_DIC
$s1  Jint_ABIO_DIC14
$s1  RESID_ABIO_DIC 
$s1  RESID_ABIO_DIC14 
EOF

if ($OCN_TAVG_TRACER_BUDGET == TRUE) then
cat >> $CASEROOT/Buildconf/pop2conf/abio_dic_dic14_tavg_contents << EOF
$s1  KPP_SRC_ABIO_DIC14
$s1  KPP_SRC_ABIO_DIC
$s1  DIA_IMPVF_ABIO_DIC14
$s1  DIA_IMPVF_ABIO_DIC
$s1  HDIFE_ABIO_DIC14
$s1  HDIFE_ABIO_DIC
$s1  HDIFN_ABIO_DIC14
$s1  HDIFN_ABIO_DIC
$s1  HDIFB_ABIO_DIC14
$s1  HDIFB_ABIO_DIC
$s1  UE_ABIO_DIC14
$s1  UE_ABIO_DIC
$s1  VN_ABIO_DIC14
$s1  VN_ABIO_DIC
$s1  WT_ABIO_DIC14
$s1  WT_ABIO_DIC
EOF
endif

#===============================================================================
# The following are fields computed by the CFC modules that are not placed in
# the tavg file by default.
# $s1 $s1  Jint_100m_DIC
# $s1 $s1  Jint_100m_DIC14
#===============================================================================
