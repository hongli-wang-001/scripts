#!/usr/bin/env python3
# read airnow data and convert to netcdf
import netCDF4 as nc
import numpy as np
import inspect, os, sys, argparse
import pandas as pd
from datetime import datetime
from pathlib import Path

IODA_CONV_PATH = Path(__file__).parent/"../lib/pyiodaconv"
if not IODA_CONV_PATH.is_dir():
    IODA_CONV_PATH = Path(__file__).parent/'..'/'lib-python'
sys.path.append(str(IODA_CONV_PATH.resolve()))
import meteo_utils
import ioda_conv_ncio as iconv
from collections import defaultdict, OrderedDict
from orddicts import DefaultOrderedDict
import pyresample
from math import radians, degrees, sin, cos, asin, acos, sqrt
def great_circle(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return 6371 * (
        acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
    )

def assign_grid(gridlon,gridlat,lon2,lat2,reso):
    ngrid=np.nan
    for i in range(len(gridlon)):
        if great_circle(gridlon[i],gridlat[i],lon2,lat2) <= reso:
            ngrid=i
            return ngrid
    return ngrid
      
def read_monitor_file(sitefile=None):

    colsinuse = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    airnow = pd.read_csv(sitefile, delimiter='|', header=None,
                         usecols=colsinuse, dtype={0: str}, encoding="ISO-8859-1")
    airnow.columns = [
        'siteid', 'Site_Code', 'Site_Name', 'Status', 'Agency',
        'Agency_Name', 'EPA_region', 'latitude', 'longitude', 'Elevation',
        'GMT_Offset', 'Country_Code', 'CMSA_Code', 'CMSA_Name', 'MSA_Code',
        'MSA_Name', 'state_Code', 'state_Name', 'County_Code',
        'County_Name', 'City_Code']
    airnow['airnow_flag'] = 'AIRNOW'
    airnow.columns = [i.lower() for i in airnow.columns]
    return airnow


def filter_bad_values(df):
    """Short summary.

    Returns
    -------
    type
        Description of returned object.

    """

    df.loc[(df.obs > 3000) | (df.obs < 0), 'obs'] = np.NaN
    return df


def long_to_wide(df):
    from pandas import Series, merge
    w = df.pivot_table(
        values='obs', index=['time', 'siteid'],
        columns='variable').reset_index()
    cols = Series(df.columns)
    g = df.groupby('variable')
    for name, group in g:
        w[name + '_unit'] = group.units.unique()[0]

    return merge(w, df, on=['siteid', 'time'])


def add_data(infile, sitefile):
    df = pd.read_csv(infile, delimiter='|',
                     header=None,
                     error_bad_lines=False,
                     encoding='ISO-8859-1')
    cols = ['date', 'time', 'siteid', 'site', 'utcoffset', 'variable', 'units',
            'obs', 'source']
    df.columns = cols
    df['obs'] = df.obs.astype(float)
    df['siteid'] = df.siteid.str.zfill(9)
    df['utcoffset'] = df.utcoffset.astype(int)
    df['time'] = pd.to_datetime(df.date + ' ' + df.time,
                                format='%m/%d/%y %H:%M',
                                exact=True)
    df.drop(['date'], axis=1, inplace=True)
    df['time_local'] = df.time + pd.to_timedelta(df.utcoffset, unit='H')

    monitor_df = read_monitor_file(sitefile)
    df = pd.merge(df, monitor_df, on='siteid')
    df.drop_duplicates(inplace=True)
    df = filter_bad_values(df)
    return long_to_wide(df)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=(
            'Reads single AIRNow text file '
            ' and converts into IODA formatted output files.')
    )

    required = parser.add_argument_group(title='arguments')
    required.add_argument(
        '-i', '--input',
        help="path of AIRNow text input file",
        type=str, required=True)
    required.add_argument(
        '-s', '--sitefile',
        help="path of AIRNow site list file",
        type=str, required=True)
    required.add_argument(
        '-grid', '--gridfile',
        help="geolat/geolon file of the grid",
        type=str, required=True)
    required.add_argument(
        '-o', '--output',
        help="path of IODA output file",
        type=str, required=True)

    args = parser.parse_args()
    print('infile=', args.input, args.sitefile)
    f = add_data(args.input, args.sitefile).drop_duplicates(subset=['PM2.5', 'OZONE', 'siteid', 'latitude', 'longitude'])
    
    f1= f.dropna(subset=['latitude','longitude'], how='any')
    f2 = f1.dropna(subset=['PM2.5','OZONE'], how='all').reset_index()
    
    aswath = pyresample.geometry.SwathDefinition(lons=f2.longitude,lats=f2.latitude)
    
    fgrid=nc.Dataset(args.gridfile,'r')
    geolon=fgrid['geolon'][:]
    geolon[geolon>180]=geolon[geolon>180]-360
    geolat=fgrid['geolat'][:]
    grid_radius=great_circle(geolon[0,0],geolat[0,0],geolon[0,1],geolat[0,1])*1000/2
    gshape=len(geolat.flatten())-1 
    grid1=pyresample.geometry.GridDefinition(lats=geolat, lons=geolon)
    
    mask_grid1, mask_aswath, index_a, distance_a = pyresample.kd_tree.get_neighbour_info(
      source_geo_def=grid1, target_geo_def=aswath, reduce_data=False, radius_of_influence=grid_radius,
      neighbours=1)
     
    grid_ind1=index_a.astype(float) 
    grid_ind1[grid_ind1>gshape]=np.nan
#    print('len=',len(f2.latitude),len(grid_ind1))
    f2['grid_indx']=grid_ind1

#    f2['grid_indx']=np.array([ assign_grid(geolon,geolat,f2['longitude'][i],f2['latitude'][i],grid_radius) 
#      for i in range(len(f2['latitude'])) ])
      
    f2a = f2.dropna(subset=['grid_indx'], how='any').reset_index()
    
    f3 = f2a.groupby('grid_indx').agg([np.size,np.mean,np.std])
    f3site = f2a.groupby('grid_indx')['siteid'].transform(lambda x: '-'.join(x)).drop_duplicates()
    
    nlocs, columns = f3.shape

    obsvars = {'pm25': 'pm25', 'o3': 'o3', }
    AttrData = {'converter': os.path.basename(__file__), }

    locationKeyList = [("latitude", "float"), ("longitude", "float"),
                       ("station_elevation", "float"), ("height", "float"), ("station_id", "string"),
                       ("datetime", "string"), ("nobs_o3","int"), ("nobs_pm25", "int")]

    writer = iconv.NcWriter(args.output, locationKeyList)

    varDict = defaultdict(lambda: defaultdict(dict))
    outdata = defaultdict(lambda: DefaultOrderedDict(OrderedDict))
    loc_mdata = defaultdict(lambda: DefaultOrderedDict(OrderedDict))
    var_mdata = defaultdict(lambda: DefaultOrderedDict(OrderedDict))
    units = {}
    units['pm25'] = 'microgram/m3'
    units['o3'] = 'ppmV'

    for i in ['pm25', 'o3']:
        varDict[i]['valKey'] = i, writer.OvalName()
        varDict[i]['errKey'] = i, writer.OerrName()
        varDict[i]['qcKey'] = i, writer.OqcName()

    d = np.empty([nlocs], 'S20')
    print(f2a.time)
    d[:] = np.array(f2a.time[1].strftime('%Y-%m-%dT%H:%M:%SZ'))
    loc_mdata['datetime'] = writer.FillNcVector(d, 'datetime')
    loc_mdata['latitude'] = f3['latitude']['mean'].values
    loc_mdata['longitude'] = f3['longitude']['mean'].values
    loc_mdata['height'] = np.full((nlocs), 10.)
    loc_mdata['station_elevation'] = f3['elevation']['mean'].values
    loc_mdata['nobs_o3'] = f3['OZONE']['size'].values.astype(int)
    loc_mdata['nobs_pm25'] = f3['PM2.5']['size'].values.astype(int)
     
    c = np.empty([nlocs], dtype='S20')
    c[:] = f3site.values
    loc_mdata['station_id'] = writer.FillNcVector(c, 'string')

    outdata[varDict['pm25']['valKey']] = np.array(f3['PM2.5']['mean'].fillna(nc.default_fillvals['f4']))
    outdata[varDict['o3']['valKey']] = np.array((f3['OZONE']['mean']/1000).fillna(nc.default_fillvals['f4']))
    
    o3err = np.fmax(f3['OZONE']['std'].values,f3['OZONE']['mean'].values/10)/1000
    o3err[np.isnan(o3err)] = 0.001
    outdata[varDict['o3']['errKey']] = np.clip(o3err, 0.001, 0.1)
    pmerr = np.fmax(f3['PM2.5']['std'].values,f3['PM2.5']['mean'].values/10)
    pmerr[np.isnan(pmerr)] = 0.1
    outdata[varDict['pm25']['errKey']] = np.clip(pmerr, 0.1, 50)
    
    for i in ['pm25', 'o3']:      
        outdata[varDict[i]['qcKey']] = np.full((nlocs), 0)

    writer._nvars = 2
    writer._nlocs = nlocs
    writer.BuildNetcdf(outdata, loc_mdata, var_mdata, AttrData, units)
