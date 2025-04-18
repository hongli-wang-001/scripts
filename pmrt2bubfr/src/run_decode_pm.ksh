#!/bin/ksh
set -euax

cwf pmbufr_decode

cyc_day=20221018
here=$(pwd)
pm_exec=${here}/pmbufr_decode.x

cd Data/BUFR/$cyc_day

rm -f log.*

for i in  $(seq -w  00 23);
do

#cat pm.bufr.01 pm.bufr.02 pm.bufr.03 pm.bufr.04 pm.bufr.05 pm.bufr.06 \
#    pm.bufr.07 pm.bufr.08 pm.bufr.09 pm.bufr.10 pm.bufr.11 pm.bufr.12 \
#    pm.bufr.13 pm.bufr.14 pm.bufr.15 pm.bufr.16 pm.bufr.17 pm.bufr.18 \
#    pm.bufr.19 pm.bufr.20 pm.bufr.21 pm.bufr.22 pm.bufr.23 pm.bufr.24 > pm.bufr.decode

ln -fs HourlyData_${cyc_day}$i.bufr pm.bufr
${pm_exec} > log.${cyc_day}$i 2>&1
rm pm.bufr

#grep subset= log.pm.$i | awk '{print $8}' > log.tm00.$i
#diff ../../log.tm00.$i .   # sequence number should equal to that from one day data

done
