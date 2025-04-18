#!/usr/bin/env python3

# Read airnow text data file and convert to IODA netcdf
import os
from datetime import datetime
import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import pyiodaconv.ioda_conv_engines as iconv
import pyiodaconv.ioda_conv_ncio as iconio
from collections import defaultdict
from pyiodaconv.orddicts import DefaultOrderedDict

os.environ["TZ"] = "UTC"

# Dictionary of output variables (ObsVal, ObsError, and PreQC).
varDict = {
    'PM2.5': ['particulatematter2p5Surface', 'ug m-3'],
    'OZONE': ['ozoneSurface', 'ppmV']
}

locationKeyList = [
    ("latitude", "float", "degrees_north"),
    ("longitude", "float", "degrees_east"),
    ("dateTime", "long", "seconds since 1970-01-01T00:00:00Z"),
    ("stationElevation", "float", "m"),
    ("height", "float", "m"),
    ("pressure", "float", "hpa"),
    ("stationIdentification", "string", "")
]

meta_keys = [m_item[0] for m_item in locationKeyList]

GlobalAttrs = {
    'converter': os.path.basename(__file__),
    'ioda_version': 2,
    'description': 'AIRNow data (converted from text/csv to IODA',
    'source': 'Unknown (ftp)'
}

iso8601_string = locationKeyList[meta_keys.index('dateTime')][2]
epoch = datetime.fromisoformat(iso8601_string[14:-1])

metaDataName = iconv.MetaDataName()
obsValName = iconv.OvalName()
obsErrName = iconv.OerrName()
qcName = iconv.OqcName()

float_missing_value = nc.default_fillvals['f4']
int_missing_value = nc.default_fillvals['i4']
double_missing_value = nc.default_fillvals['f8']
long_missing_value = nc.default_fillvals['i8']
string_missing_value = '_'

missing_vals = {
    'string': string_missing_value,
    'integer': int_missing_value,
    'long': long_missing_value,
    'float': float_missing_value,
    'double': double_missing_value
}

dtypes = {
    'string': object,
    'integer': np.int32,
    'long': np.int64,
    'float': np.float32,
    'double': np.float64
}

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
    if sitefile == 'True':
        df.loc[(df.obs > 3000) | (df.obs < 0), 'obs'] = np.NaN
    else:
        df.loc[(df.pm25 > 3000) | (df.pm25 < 0), 'pm25'] = np.NaN
        df.loc[(df.ozone > 3000) | (df.ozone < 0), 'ozone'] = np.NaN
    return df

def long_to_wide(df):
    w = df.pivot_table(values='obs', index=['time', 'siteid'],
                       columns='variable').reset_index()
    g = df.groupby('variable')
    for name, group in g:
        w[name + '_unit'] = group.units.unique()[0]
    return w

def add_data(infile, sitefile):
    df = pd.read_csv(infile, delimiter='|',
                     header=None,
                     on_bad_lines='warn',
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

def average_data(df):
    # Bin latitude and longitude
    df['lat_bin'] = (df['latitude'] // 0.1) * 0.1
    df['lon_bin'] = (df['longitude'] // 0.1) * 0.1
    return df.groupby(['lat_bin', 'lon_bin']).agg(
        PM2_5=('PM2.5', 'mean'),
        OZONE=('OZONE', 'mean'),
        latitude=('latitude', 'mean'),
        longitude=('longitude', 'mean')
    ).reset_index()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Reads single AIRNow text file and converts into IODA formatted output files.')
    
    sitefile = os.environ['sitefile']
    required = parser.add_argument_group(title='required arguments')
    
    if sitefile == 'True':
        required.add_argument(
            '-s', '--sitefile',
            help="path of AIRNow site list file",
            type=str, required=True)

    required.add_argument(
        '-i', '--input',
        help="path of AIRNow text input file",
        type=str, required=True)
    
    required.add_argument(
        '-o', '--output',
        help="path of IODA output file",
        type=str, required=True)

    args = parser.parse_args()
    
    if sitefile == 'True':
        print('infile=', args.input, args.sitefile)
        f = add_data(args.input, args.sitefile).drop_duplicates(subset=['PM2.5', 'OZONE', 'siteid', 'latitude', 'longitude'])
    else:
        print('infile=', args.input)
        f = add_data(args.input)

    f3 = f.dropna(subset=['PM2.5'], how='any').reset_index()
    f3 = average_data(f3)  # Apply averaging function
    nlocs, columns = f3.shape

    dt = f3.time[1].to_pydatetime()
    time_offset = round((dt - epoch).total_seconds())

    ioda_data = {}
    data = {}
    for key in varDict.keys():
        data[key] = []
    for key in meta_keys:
        data[key] = []

    # Load CMAQ lat lon 
    ll_file = 'lonlat.nc'
    with nc.Dataset(ll_file, 'r') as ll_nc:
        lon = ll_nc.variables['grid_lont'][:]
        lat = ll_nc.variables['grid_latt'][:]

    hgt_file = 'phis.nc'
    with nc.Dataset(hgt_file, 'r') as hgt_nc:
        hgt = hgt_nc.variables['phis'][0, :, :]
        hgt = hgt / 9.8

    latitude = np.array(f3['latitude'])
    longitude = np.array(f3['longitude'])
    longitude = np.where(longitude < 0, longitude + 360, longitude)
    
    # Interpolation to observation locations
    print('longitude:', longitude.size)
    interpolated_hgt = np.full_like(longitude, np.nan)
    if longitude.size > 0:
        interpolated_hgt = griddata(
            (lon.flatten(), lat.flatten()),
            hgt.flatten(),
            (longitude, latitude),
            method='linear',
            fill_value=np.nan
        )

    # Fill the temporary data arrays from input file column data
    data['stationIdentification'] = np.full(nlocs, f3.siteid, dtype='S20')
    data['dateTime'] = np.full(nlocs, np.int64(time_offset))
    data['latitude'] = np.array(f3['latitude'])
    data['longitude'] = np.array(f3['longitude'])
    data['stationElevation'] = interpolated_hgt
    data['height'] = np.full(nlocs, float_missing_value, dtype=dtypes['float'])

    for key in varDict.keys():
        data[varDict[key][0]] = np.array(f3[key])

    # Finalize IODA data for writing
    ioda_data = {}
    for key in varDict.keys():
        ioda_data[varDict[key][0]] = data[varDict[key][0]]

    # Setup the IODA writer and write everything out.
    writer = iconv.IodaWriter(args.output, locationKeyList)
    writer.BuildIoda(ioda_data, None, None, GlobalAttrs)


