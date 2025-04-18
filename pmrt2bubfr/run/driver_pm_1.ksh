#!/bin/ksh -l

#SBATCH --partition=bigmem
#SBATCH --mem=250g
#SBATCH --ntasks=1
#SBATCH --qos=batch  # partition bigmem: batch debug urgent 
#SBATCH --time=0:30:00
#SBATCH --account=zrtrr
#SBATCH --job-name=pm_prepbufr
##SBATCH --output log.pm

set -euax

cwf pmbufr

file_dir=/mnt/lfs1/BMC/wrfruc/hwang/fire2/pmrt2bubfr/Data
#/mnt/lfs4/BMC/wrfruc/Ruifang.Li/WF2/PM/Data
here=$(pwd)
pm_exec=${here}/pmbufr.x

mkdir -p ${file_dir}/BUFR
cd ${file_dir}/BUFR

yy=2022
mm=10

for dd in 18; do 

  cyc_day=${yy}${mm}${dd}
  mkdir -p ${cyc_day}
  cd ${cyc_day}
  rm -f *

  ln -sf ${here}/pm.table . 
  ln -sf ${file_dir}/${cyc_day}/Monitoring_Site_Locations_V2.dat .
  site=Monitoring_Site_Locations_V2.dat

  # remove | in site file and get only PM2.5 site
  awk '-F|' '{print $1,$4,$8,$12,$13}' ${site} | grep PM2.5 > site_loc
  cnt_site=`wc  -l < site_loc`
  
  for k in $(seq -w  00 23); do 

    echo 'cycle='${cyc_day}${k}
    ln -sf ${file_dir}/${cyc_day}/HourlyData_${cyc_day}${k}.dat .
    fname=HourlyData_${cyc_day}${k}.dat

    # remove | and / from hourly data
    awk '-F|' '{print $1,$2,$3,$6,$8}' ${fname} | grep PM2.5 > fvar_$k
    awk '-F/' '{print $1,$2,$3,$4,$5,$6}' fvar_$k > file_var_$k
    cnt_rec=`wc  -l < file_var_$k`
    cp file_var_$k file_var
    ${pm_exec} "site_loc" ${cnt_site} "file_var" ${cnt_rec} ${cyc_day}${k}
  done


  rm fvar_*
  cd ..

done


exit


#sed 's/|/  /g' $site_loc > site_loc1                                            # replace | with empty space
#awk '/PM2.5/{++cnt} END {print "Count = ", cnt}' HourlyData_${cyc_day}${k}.dat  # count number which match 
