#https://stackoverflow.com/questions/70622457/chanaging-standard-name-and-long-name-with-nco
#ncatted -a standard_name,pt,o,c,air_potential_temperature -a long_name,pt,o,c,Potential_temperature pt_19891020-19891022.nc
rm -fr numatkn.nc pm25at_out.nc 
ncks -v numatkn $1 numatkn.nc
ncap2 -F -s "pm25ac=(0.0*numatkn+1.0)" numatkn.nc pm25at_out.nc
#ncrename -v flash_extent_density,fed_old $1
ncks -h -A -v pm25ac pm25at_out.nc $1
ncatted -h -a long_name,pm25ac,o,c,pm25ac $1
ncdump -hs $1
