#rm -fr pm25at.nc 
ncks -v numatkn $1 numatkn.nc
ncap2 -F -s "pm25at=(0.0*numatkn+1.0)" numatkn.nc pm25at_out.nc
#ncrename -v flash_extent_density,fed_old $1
ncks -h -A -v pm25at pm25at_out.nc $1
ncatted -h -a long_name,pm25at,o,c,pm25at $1
ncdump -hs $1
