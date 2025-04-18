#!/bin/csh

foreach dir (vrf_dapapm25_w2h_conus13_wsep2020)
foreach hz ( 00 12 )

if ( -d ${dir} ) then
 mkdir -p ${hz}z
 grep avg $dir/202009*${hz}/log* | awk '{print $3 }' > ${hz}z/log.obs.ts.${dir}
 grep avg $dir/202009*${hz}/log* | awk '{print $4 }' > ${hz}z/log.fcst.ts.${dir}
 grep rms $dir/202009*${hz}/log* | awk '{print $3 }' > ${hz}z/log.rmsfmo.ts.${dir}
 grep avefmo $dir/202009*${hz}/log* | awk '{print $3 }' > ${hz}z/log.avefmo.ts.${dir}
 grep stdfmo $dir/202009*${hz}/log* | awk '{print $3 }' > ${hz}z/log.stdfmo.ts.${dir}
 grep "Num of valid obs" $dir/202009*${hz}/log* | awk '{print $9 }' > ${hz}z/log.numobs.ts.${dir}
endif
end
end

