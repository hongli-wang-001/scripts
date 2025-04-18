#https://stackoverflow.com/questions/70622457/chanaging-standard-name-and-long-name-with-nco
#ncatted -a standard_name,pt,o,c,air_potential_temperature -a long_name,pt,o,c,Potential_temperature pt_19891020-19891022.nc
rm -fr numatkn.nc pm25at_out.nc 
ncks -v numatkn $1 numatkn.nc
ncap2 -F -s "pm25co=(0.0*numatkn+0.05)" numatkn.nc pm25at_out.nc
#ncrename -v flash_extent_density,fed_old $1
ncks -h -A -v pm25co pm25at_out.nc $1
ncatted -h -a long_name,pm25co,o,c,pm25co $1
ncdump -hs $1
##ncatted -a history_of_appended_files,global,d,, ori.20200903.180000.fv_tracer.res.tile1.nc
#ncdump -hs ori.20200903.180000.fv_tracer.res.tile1.nc
#ncatted -a history,global,d,, ori.20200903.180000.fv_tracer.res.tile1.nc

