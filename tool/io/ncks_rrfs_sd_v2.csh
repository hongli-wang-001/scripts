#https://stackoverflow.com/questions/70622457/chanaging-standard-name-and-long-name-with-nco
#ncatted -a standard_name,pt,o,c,air_potential_temperature -a long_name,pt,o,c,Potential_temperature pt_19891020-19891022.nc
#https://sourceforge.net/p/nco/discussion/9829/thread/13174fed94/
#ncap2 -s 'defdim("lengthd",10);date[lengthd,time]=AN_ARRAY_OF_THIS_SIZE' input.nc output.nc

ncap2 -s 'defdim("nz_mixing_ratio_of_smoke_wrt_dry_air",70)' aero_gauss_state_F12.nc
ncap2 -s 'defdim("nz_mixing_ratio_of_dust_wrt_dry_air",70)' aero_gauss_state_F12.nc
ncap2 -s 'defdim("nz_mass_density_of_particulate_matter_2p5_in_air",70)' aero_gauss_state_F12.nc 
ncap2 -s "mixing_ratio_of_smoke_wrt_dry_air[nz_mixing_ratio_of_smoke_wrt_dry_air,ny,nx]=0.1" aero_gauss_state_F12.nc
ncap2 -s "mixing_ratio_of_dust_wrt_dry_air[nz_mixing_ratio_of_dust_wrt_dry_air,ny,nx]=0.2" aero_gauss_state_F12.nc
ncap2 -s "mass_density_of_particulate_matter_2p5_in_air[nz_mass_density_of_particulate_matter_2p5_in_air,ny,nx]=0.3" aero_gauss_state_F12.nc
ncrename -v dry_air_density_levels_minus_one,ori_dry_air_density_levels_minus_one aero_gauss_state_F12.nc
ncap2 -s "dry_air_density_levels_minus_one[nz_dry_air_density_levels_minus_one,ny,nx]=1.0" aero_gauss_state_F12.nc
#ncap2 -s 'mixing_ratio_of_smoke_wrt_dry_air[nz_mixing_ratio_of_smoke_wrt_dry_air,ny,nx]=0.1'
#ncap2 -s 'mixing_ratio_of_dust_wrt_dry_air[nz_mixing_ratio_of_dust_wrt_dry_air,ny,nx]=0.2'
dry_air_density_levels_minus_one

Exit
rm -fr temp.nc
ncks -v specific_humidity $1 temp.nc
ncap2 -F -s "mixing_ratio_of_smoke_wrt_dry_air=(0*specific_humidity+0.2)" temp.nc temp_smoke.nc
ncap2 -F -s "mixing_ratio_of_dust_wrt_dry_air=(0*specific_humidity+0.1)" temp.nc temp_dust.nc
ncap2 -F -s "mass_density_of_particulate_matter_2p5_in_air=(0*specific_humidity+0.3)" temp.nc temp_pm2p5.nc
#ncap2 -F -s "dry_air_density_levels_minus_one=(specific_humidity/0.018)" temp.nc temp_air_dens.nc
#ncrename -v flash_extent_density,fed_old $1
ncks -A -v mixing_ratio_of_smoke_wrt_dry_air temp_smoke.nc $1
ncks -A -v mixing_ratio_of_dust_wrt_dry_air temp_dust.nc $1
ncks -A -v mass_density_of_particulate_matter_2p5_in_air temp_pm2p5.nc $1

#ncatted -a long_name,coarsepm,o,c,coarsepm $1
ncdump -hs $1
