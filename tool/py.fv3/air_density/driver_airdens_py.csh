rm airdens.nc
ncks -v pm25_tot gfsic.nc airdens.nc
ncrename -v pm25_tot,dry_air_density airdens.nc
ncks -v pm25_tot gfsic.nc airdens_res.nc
python3 compute_airdens_v5.py

