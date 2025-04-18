ln -sf /scratch2/BMC/wrfruc/hwang/wf1/ncl_script/*v3*ncl .
ln -sf  /scratch2/BMC/wrfruc/hwang/wf1/ncl_script/inc-group_aero_pm2p5_13km.ncl . 

#ncl inc-group_aero_pm2p5_13km.ncl  > log.group_aero_fv3lam_ugm2
ncl obs_v3_aod.ncl > log.obs.aod
ncl cmaq_v3_aod.ncl > log.bkg.aod
ncl anl_v3_aod.ncl > log.anl.aod
ncl omb_v3_aod.ncl > log.omb.aod
ncl oma_v3_aod.ncl > log.oma.aod

ncl obs_v3_pm25.ncl > log.obs.pm25 
ncl cmaq_v3_pm25.ncl > log.bkg.pm25 
ncl anl_v3_pm25.ncl > log.anl.pm25
ncl omb_v3_pm25.ncl > log.omb.pm25
ncl oma_v3_pm25.ncl > log.oma.pm25
