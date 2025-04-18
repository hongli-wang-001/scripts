ncap2 -h -s 'defdim("nz_mixing_ratio_of_smoke_wrt_dry_air",70)' aero_gauss_state_F12.nc
ncap2 -h -s 'defdim("nz_mixing_ratio_of_dust_wrt_dry_air",70)' aero_gauss_state_F12.nc
ncap2 -h -s 'defdim("nz_mass_density_of_particulate_matter_2p5_in_air",70)' aero_gauss_state_F12.nc 
ncap2 -h -s "mixing_ratio_of_smoke_wrt_dry_air[nz_mixing_ratio_of_smoke_wrt_dry_air,ny,nx]=0.1" aero_gauss_state_F12.nc
ncap2 -h -s "mixing_ratio_of_dust_wrt_dry_air[nz_mixing_ratio_of_dust_wrt_dry_air,ny,nx]=0.2" aero_gauss_state_F12.nc
ncap2 -h -s "mass_density_of_particulate_matter_2p5_in_air[nz_mass_density_of_particulate_matter_2p5_in_air,ny,nx]=0.3" aero_gauss_state_F12.nc
ncrename -h -v dry_air_density_levels_minus_one,ori_dry_air_density_levels_minus_one aero_gauss_state_F12.nc
ncap2 -h -s "dry_air_density_levels_minus_one[nz_dry_air_density_levels_minus_one,ny,nx]=1.0" aero_gauss_state_F12.nc
#ncdump -hs $1
