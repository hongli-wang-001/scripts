#!/bin/ksh -f

# module purge
# module load python/3.6.3
# ON Hera, 20220218
#	module purge
module use -a /contrib/anaconda/modulefiles
module load intel
module load anaconda/latest

set -x

#ln -s /scratch2/BMC/wrfruc/hwang/wf1/data_monet/AIRNOW_20190801_20190831_b.nc . 
# check ncl file
export std=2020092606
export etd=2020093118

#exp="dapm25_w2h_conus13_wsep2020"
#exp="daaod_pm25_w2h_conus13_wsep2020"
exp="dapapm25_w2h_conus13_wsep2020"
#exp="daaod2h_conus13_wsep2020"
export out=/scratch2/BMC/wrfruc/hwang/wf1/post/out_case2020sep_v2/${exp}
export obsdir=/scratch2/BMC/wrfruc/hwang/wf1/data/pm2p5_airnow_case2020sep/decode_bufr/hour 
#/scratch2/BMC/wrfruc/hwang/wf1/data_monet/hour
# ./vrf_rrfcmaq_mypm2p5_v4.ncl 
while [ $std -le $etd ]
do
   
   export outdir=${out}/$std/pm25
   mkdir -p $outdir 
   export sth=01 
   export eth=24
   export year=`echo $std | cut -c1-4`
   export month=`echo $std | cut -c5-6`
   export dd=`echo $std | cut -c7-8`
   export cyc=`echo $std | cut -c9-10`

   mkdir -p vrf_${exp}/$std

   while [ $sth -le $eth ] ; do
   hr=$sth
   typeset -Z3 sth
#   infile=${indir}/$std/dynf${sth}.nc
   outfile1=${outdir}/rcmaq.t${cyc}z.pm25tot_f${sth}.nc

   export vrf_date=`/home/Hongli.Wang/bin_theia/da_advance_time.exe $std ${sth}h` 
   echo "vrf_date= $vrf_date"  
   export vrf_date_pre=`/home/Hongli.Wang/bin_theia/da_advance_time.exe $vrf_date -1h` 
   ln -sf $obsdir/pm2p5.${vrf_date}.txt pm2p5.txt
   export vrf_day=`echo $vrf_date | cut -c1-8`
#   ln -sf $obsdir/pm2p5.${vrf_date_pre}.txt pm2p5_pre.txt
#20190806/pm25.airnow.2019080618
   ln -sf $obsdir/$vrf_day/pm25.airnow.$vrf_date pm2p5.txt
   ln -sf $outfile1 rcmaq.pm25tot_fcst.nc 
   ncl vrf_rrfcmaq_mypm2p5_v4.ncl > vrf_${exp}/$std/log.vrfpm2p5.f${sth}.txt 
   mv vobs.txt  vrf_${exp}/$std/vobs.f${sth}.txt

   let "sth=sth+1" 
    
   done

#   export std=`$NDATE +24 $std` 
#   export std=`/gpfs/dell1/nco/ops/nwprod/prod_util.v1.1.0/exec/ips/ndate +24 $std`
#   export std=$(/home/Ruifang.Li/bin/bumpidx ${std} 6)
   export std=`/home/Hongli.Wang/bin_theia/da_advance_time.exe $std 6h` 
done
  

