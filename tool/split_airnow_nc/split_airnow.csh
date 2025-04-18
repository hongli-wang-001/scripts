#!/bin/ksh -f

set -x

export std=2020090100
export etd=2020093100

export sth=0
while [ $std -le $etd ]
do
   export year=`echo $std | cut -c1-4`
   export month=`echo $std | cut -c5-6`
   export dd=`echo $std | cut -c7-8`
   export hh=`echo $std | cut -c9-10`
   export cyc=`echo $std | cut -c9-10`

   ncks -d time,${sth},${sth} AIRNOW_20200901_20200930.nc AIRNOW_${std}.nc

   let "sth=sth+1"

   export std=`/home/Hongli.Wang/bin_theia/da_advance_time.exe $std 1h`
done

