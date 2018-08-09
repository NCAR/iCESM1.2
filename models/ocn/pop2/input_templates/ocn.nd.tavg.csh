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
   echo invalid my_stream number $my_stream
   exit 5
endif

@ s1 = 1   # use base-model stream 1

cat >! $CASEROOT/Buildconf/pop2conf/nd_tavg_contents << EOF
$s1  ND143
$s1  ND144

EOF

#$s1  ND143_P
#$s1  ND144_P
#$s1  ND143_D
#$s1  ND144_D
#$s1  ND143_SOURCE
#$s1  ND143_SINK
#$s1  ND144_SOURCE
#$s1  ND144_SINK


#-------------------------------------------------------------------------------------
# Add optional tracer budget terms
#-------------------------------------------------------------------------------------
if ($OCN_TAVG_TRACER_BUDGET == TRUE) then
cat >> $CASEROOT/Buildconf/pop2conf/nd_tavg_contents << EOF
$s1  ND143_RESET_TEND
$s1  ND144_RESET_TEND
$s1  DIA_IMPVF_ND
$s1  HDIFE_ND
$s1  HDIFN_ND
$s1  HDIFB_ND
$s1  UE_ND
$s1  VN_ND
$s1  WT_ND
EOF
endif

#  disable the following until they are computed correctly
#  IAGE_SQR 
#  UE_IAGE
#  VN_IAGE
#  WT_IAGE
#  ADV_IAGE
#  J_IAGE
#  Jint_IAGE
#  STF_IAGE
#  RESID_IAGE
#  FvPER_IAGE
#  FvICE_IAGE

