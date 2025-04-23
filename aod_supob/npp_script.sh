#!/bin/sh -l

#set -ux

source ./app_modules.sh

# This script is to convert viirs obs to ioda format using gdas_obsprovider2ioda.x.

year=2020
months="09"
days="01" #"02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30"
hours="00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23"

sat=npp  #npp j01
inpath=/scratch1/BMC/wrfruc/hwang/wf1/data/case2020sep/aod_viis/${sat}/001/
outpath=./Data/${sat}/stride16_minob36_qc0  

mkdir -p $outpath
mkdir -p YAMLs

#thinning="0.999"
#thinning="0.99"
#thinning="0.9"
thinning="0.0"   # superobbing, no thinning

start_time=$(date +%s.%N)

for month in $months; do
   for day in $days; do
      if [[ -d $inpath ]];then
        for hh in $hours; do
           echo  ${year} ${month} ${day}  ${hh}

           # N20
           # run run_viirs2ioda_n20.py and retrun all the obs files within the time window
           inputfiles=$(python run_viirs2ioda_npp.py ${year} ${month} ${day} ${hh} ${inpath} ${outpath})
           #echo $inputfiles

           # the processed output filename
           #outputfile=$outpath/${year}${month}${day}"/gdas.t"$hh"z.viirs_n20."$year$month$day$hh".nc4"
           outputfile=$outpath/${year}${month}${day}/"ioda.viirs_npp."$year$month$day$hh".nc4"
           python Edit_yaml_weight.py $year$month$day$hh $outputfile $thinning "4" "npp" $inputfiles


           # NPP
           # run run_viirs2ioda_npp.py and retrun all the obs files within the time window
           #inputfiles=$(python run_viirs2ioda_npp.py ${year} ${month} ${day} ${hh} ${inpath} ${outpath})
           #echo $inputfiles

           # the processed output filename
           #outputfile=$outpath/${year}${month}${day}"/gdas.t"$hh"z.viirs_npp."$year$month$day$hh".nc4"
           #python Edit_yaml_weight.py $year$month$day$hh $outputfile $thinning "4" "npp" $inputfiles

        done
      fi
   done
done


end_time=$(date +%s.%N)
elapsed_time=$(echo "$end_time - $start_time" | bc)

echo $elapsed_time
