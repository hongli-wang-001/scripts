#!/bin/ksh 


set -euax

pmfile=/scratch2/BMC/wrfruc/rli/WF1/PM/Data
fname=pm25.airnow

yy=2019
mm=08

cd ${pmfile}/anowpm_hourly

for dd in 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
do

  cyc_day=${yy}${mm}${dd}
  cd ${cyc_day}

  for hh in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23
  do
     mv  pm25.bufr.${hh} $fname.${cyc_day}${hh}.bufr 
  done 
  
  dd_next=$((dd + 1))
  if [ ${dd_next} -lt 10  ]; then
    dd_next=0${dd_next}
  fi

  mv pm25.bufr.24 ../${yy}${mm}${dd_next}/${fname}.${yy}${mm}${dd_next}00.bufr
  
  cd ..

done
