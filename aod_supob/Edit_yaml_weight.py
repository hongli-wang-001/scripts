import sys
import yaml
from datetime import datetime,timedelta
import os
import subprocess as sp

def process_arguments(date, outputfile, thinning, channel,sensor, inputfiles):
    # Your code to process the arguments
    print("date: ", date, " outputfile: ", outputfile, " thinning thresould: ",thinning, " channle no.: ",channel, " sensor: ",sensor)
    print("List of inputfiles:")
    for string in inputfiles:
        print(string)


date = sys.argv[1]
outputfile = sys.argv[2]
thinning = float(sys.argv[3])
channel = int(sys.argv[4])
sensor = sys.argv[5]
inputfiles = sys.argv[6:]

process_arguments(date,outputfile, thinning, channel, sensor, inputfiles)

date_anl = datetime.strptime(str(date),'%Y%m%d%H')
window_begin = date_anl - timedelta(hours=1)
window_end = date_anl + timedelta(hours=0)


with open('./gdas_viirsaod2ioda_weightfunc_org_C96.yaml') as fp:
#with open('./gdas_viirsaod2ioda_org.yaml') as fp:
#with open('./gdas_viirsaod2ioda_weightfunc_org_C192.yaml') as fp:
#with open('./gdas_viirsaod2ioda_weightfunc_org_C96_thin.yaml') as fp:  # thin the obs to the superobbing grids without superobbing
#with open('./gdas_viirsaod2ioda_weightfunc_org_C192_thin.yaml') as fp:  # thin the obs to the superobbing grids without superobbing
    data = yaml.safe_load(fp)

#print(data['input files'] )
data['window begin'] = window_begin.strftime('%Y-%m-%dT%H:%M:%S')+'Z'
data['window end'] = window_end.strftime('%Y-%m-%dT%H:%M:%S')+'Z'
data['output file'] = outputfile
data['input files'] = inputfiles 

#print(data['input files'] )
#for i, file in enumerate(inputfiles):
#  data['input files'][i] = file

data['thinning']['threshold'] = thinning
data['channel'] = channel
#print(data['input files'] )

#gdas_yaml_temp = '/scratch1/NCEPDEV/da/Yaping.Wang/JEDI_work/python_wrappers/Viirs_x/YAMLs/gdas_viirsaod2ioda_'+sensor+'_'+date+'.yaml'
gdas_yaml_temp = './YAMLs/gdas_viirsaod2ioda_'+sensor+'_'+date+'.yaml'
with open(gdas_yaml_temp,'w') as f:
     yaml.dump(data, f)


#cmd='/scratch2/BMC/wrfruc/rli/WF1/GDASApp_pm_20240726/build/bin/gdas_obsprovider2ioda.x '+gdas_yaml_temp
#cmd='/scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm/GDASApp_pm_20240726/build/bin/gdas_obsprovider2ioda.x '+gdas_yaml_temp
cmd='/scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm/GDASapp_update/GDASApp/build/bin/gdas_obsprovider2ioda.x '+gdas_yaml_temp
my_env = os.environ.copy()
my_env['OMP_NUM_THREADS'] = '4'

proc = sp.Popen(cmd,env=my_env,shell=True)
proc.communicate()  # wait for this job to complete

