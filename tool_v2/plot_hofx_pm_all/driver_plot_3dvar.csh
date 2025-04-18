python3 obserr-pm25.py $1 OBS_ERROR OBSERR-PM2p5
python3 simulated-amb-pm25.py $1 A-B AMB-PM2p5
python3 obs-pm25_v2.py     $1 AirNOW_PM2.5 OBS_PM2p5
python3 obs_pm25_outlier_v3.py $1 AirNOW_PM2.5 OBS_PM2p5_outlier
python3 simulated-bkg-pm25.py  $1 bkg_PM2.5 BKG_PM2p5
python3 simulated-ana-pm25.py  $1 ana_PM2.5 ANA_PM2p5
python3 omb_3dvar_pm25_v2.py $1 O-B OMB_PM2p5_VADER >! log.omb.$1
python3 oma_3dvar_pm25_v2.py $1 O-A OMA_PM2p5_VADER >! log.oma.$1
mkdir -p figure_hofx_$1
mv *.png figure_hofx_$1
