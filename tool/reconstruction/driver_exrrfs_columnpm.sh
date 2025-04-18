#!/bin/bash

 export nwges_dir=/scratch2/BMC/wrfruc/hwang/fire2_jedi/tool/reconstruction
 export restart_prefix=20200912.180000
 ncap2 -O -v -s 'column_pm2p5_total=pm25_tot' ${nwges_dir}/RESTART/${restart_prefix}.fv_tracer.res.tile1.nc ${nwges_dir}/RESTART/${restart_prefix}.tmp.nc
 ncap2 -O -v -s 'pm2p5_total=pm25_tot' ${nwges_dir}/RESTART/${restart_prefix}.fv_tracer.res.tile1.nc ${nwges_dir}/RESTART/${restart_prefix}.tmp.nc
 python -u exrrfs_process_columnpm.py
