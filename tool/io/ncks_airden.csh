	ncks -v dry_air_density /scratch2/BMC/wrfruc/hwang/fire2_jedi/test/JEDI_BEC_TEST/BUMP_02/Data/inputs/lam_cmaq/bkg/20200912.120000.fv_core.res.tile1_with_airden.nc airden.nc
  1847	0:13	history | grep REST
  1848	0:14	ncks -h -A -v dry_air_density airden.nc /scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm/expt_dirs/../nco_dirs/com/aqm/v7.0/aqm.20200901/12/RESTART/20200902.120000.fv_core.res.tile1.nc
  1849	0:15	ncdump -hs /scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm/expt_dirs/../nco_dirs/com/aqm/v7.0/aqm.20200901/12/RESTART/20200902.120000.fv_core.res.tile1.nc

