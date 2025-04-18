#https://stackoverflow.com/questions/70622457/chanaging-standard-name-and-long-name-with-nco
#ncatted -a standard_name,pt,o,c,air_potential_temperature -a long_name,pt,o,c,Potential_temperature pt_19891020-19891022.nc
#https://sourceforge.net/p/nco/discussion/9829/thread/13174fed94/
#ncap2 -s 'defdim("lengthd",10);date[lengthd,time]=AN_ARRAY_OF_THIS_SIZE' input.nc output.nc

rm -fr temp.nc
ncks -v specific_humidity $1 temp.nc
ncap2 -F -s "mixing_ratio_of_smoke_wrt_dry_air=(0*specific_humidity+0.2)" temp.nc temp_smoke.nc
ncap2 -F -s "mixing_ratio_of_dust_wrt_dry_air=(0*specific_humidity+0.1)" temp.nc temp_dust.nc
ncap2 -F -s "mass_density_of_particulate_matter_2p5_in_airr=(0*specific_humidity+0.3)" temp.nc temp_pm2p5.nc
#ncap2 -F -s "dry_air_density_levels_minus_one=(specific_humidity/0.018)" temp.nc temp_air_dens.nc
#ncrename -v flash_extent_density,fed_old $1
ncks -A -v mixing_ratio_of_smoke_wrt_dry_air temp_smoke.nc $1
ncks -A -v mixing_ratio_of_dust_wrt_dry_air temp_dust.nc $1
ncks -A -v mass_density_of_particulate_matter_2p5_in_airr temp_pm2p5.nc $1
ncap2 -s 'defdim("nz_mixing_ratio_of_smoke_wrt_dry_air",70)' $1
ncap2 -s 'defdim("nz_mixing_ratio_of_dust_wrt_dry_air",70)' $1 
ncap2 -s 'defdim("nz_mass_density_of_particulate_matter_2p5_in_airr",70)' $1 
#ncatted -a long_name,coarsepm,o,c,coarsepm $1
ncdump -hs $1
