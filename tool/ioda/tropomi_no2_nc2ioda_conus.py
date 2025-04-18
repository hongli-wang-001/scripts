#!/usr/bin/env python3

#
# (C) Copyright 2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

import argparse
import netCDF4 as nc
import numpy as np
import os

import pyiodaconv.ioda_conv_engines as iconv
from collections import defaultdict, OrderedDict
from pyiodaconv.orddicts import DefaultOrderedDict

#[HChoi, 6/23/2023]
import os, sys
import math as m

locationKeyList = [
    ("latitude", "float"),
    ("longitude", "float"),
    ("dateTime", "string"),
]

AttrData = {
    'converter': os.path.basename(__file__),
    'nvars': np.int32(1),
}

DimDict = {
}


class tropomi(object):
    def __init__(self, filenames, varname, columnType, qa_flg, thin, obsVar):
        self.filenames = filenames
        self.varname = varname
        self.columnType = columnType
        self.qa_flg = qa_flg
        self.thin = thin
        self.obsVar = obsVar
        self.varDict = defaultdict(lambda: defaultdict(dict))
        self.outdata = defaultdict(lambda: DefaultOrderedDict(OrderedDict))
        self.varAttrs = DefaultOrderedDict(lambda: DefaultOrderedDict(dict))
        self._read()

    # Open input file and read relevant info
    def _read(self):
        # set up variable names for IODA
        varname_str = list(self.obsVar.keys())[0]
        print('Processing variable: %s' % (varname_str), flush=1)
        iodavar = self.obsVar[varname_str]
        self.varDict[iodavar]['valKey'] = iodavar, iconv.OvalName()
        self.varDict[iodavar]['errKey'] = iodavar, iconv.OerrName()
        self.varDict[iodavar]['qcKey'] = iodavar, iconv.OqcName()
        self.varAttrs[iodavar, iconv.OvalName()]['coordinates'] = 'longitude latitude'
        self.varAttrs[iodavar, iconv.OerrName()]['coordinates'] = 'longitude latitude'
        self.varAttrs[iodavar, iconv.OqcName()]['coordinates'] = 'longitude latitude'
        self.varAttrs[iodavar, iconv.OvalName()]['units'] = 'mol m-2'
        self.varAttrs[iodavar, iconv.OerrName()]['units'] = 'mol m-2'
        # loop through input filenames
        first = True
        for f in self.filenames:
            ncd = nc.Dataset(f, 'r')

            # get global attributes
            AttrData['date_time_string'] = ncd.getncattr('time_reference')[0:19]+'Z'
            AttrData['sensor'] = ncd.getncattr('sensor')
            AttrData['platform'] = ncd.getncattr('platform')

            # many variables are time, scanline, ground_pixel
            # but others are just time, scanline
            lats = ncd.groups['PRODUCT'].variables['latitude'][:].ravel()
            nlocs = len(lats)
            lons = ncd.groups['PRODUCT'].variables['longitude'][:].ravel()
            qa_value = ncd.groups['PRODUCT'].variables['qa_value'][:]  # 2D
            times = np.empty_like(qa_value, dtype=object)
            qa_value = qa_value.ravel()
            nlevs = ncd.groups['PRODUCT'].dimensions['layer'].size

            # adding ability to pre filter the data using the qa value
            # and also perform thinning using random uniform draw
            qaf = qa_value > self.qa_flg
#[HChoi 6/23/2023]
            print ('Checking if Obs falls inside the domain... it takes a few minutes.')
            chk = np.zeros(len(qaf), dtype=bool)
            ict = 0
            for idx in range(len(lats)):
                if (lats[idx] > 20.0 and lats[idx] < 55.0 and lons[idx] > -135.0 and lons[idx] < -60.0):
                   TorF = rrfs_domain_check(lats[idx], lons[idx])
                   chk[idx] = not TorF     #reversed so "True => inside"  "False => outside"

            thi = np.random.uniform(size=len(lons)) > self.thin
#[HChoi 6/23/2023]
           #flg = np.logical_and(qaf, thi)
            flg = np.logical_and(np.logical_and(qaf, thi), chk)
            qc_flag = ncd.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS']\
                .variables['processing_quality_flags'][:]
            qc_flag = qc_flag.ravel().astype('int32')
            time1 = ncd.groups['PRODUCT'].variables['time_utc'][:]
            for t in range(len(time1[0])):
                times[0, t, :] = time1[0, t][0:19]+'Z'
            times = times.ravel()

            if self.varname == 'no2':
                # grab the averaging kernel and reshape it
                avg_kernel = ncd.groups['PRODUCT'].variables['averaging_kernel'][:]
                avg_kernel = np.flip(np.reshape(avg_kernel, (nlocs, nlevs)), axis=1)

                if self.columnType == 'tropo':
                    trop_layer = ncd.groups['PRODUCT'].variables['tm5_tropopause_layer_index'][:].ravel()
                    total_airmass = ncd.groups['PRODUCT'].variables['air_mass_factor_total'][:].ravel()
                    trop_airmass = ncd.groups['PRODUCT'].variables['air_mass_factor_troposphere'][:].ravel()
                    # do not loop over nlocs here this makes the execution very slow
                    for k in range(nlevs):
                        avg_kernel[..., k][np.full((nlocs), k, dtype=int) <= trop_layer] = 0
                        avg_kernel[..., k] *= total_airmass / trop_airmass

                # construct the pressure vertices array
                ps = ncd.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['INPUT_DATA'].\
                    variables['surface_pressure'][:]
                # bottom of layer is vertice 0, very top layer is TOA (0hPa)
                ak = ncd.groups['PRODUCT'].variables['tm5_constant_a'][:, :]
                bk = ncd.groups['PRODUCT'].variables['tm5_constant_b'][:, :]
                preslv = np.flip(np.transpose(ak[..., 0][:, np.newaxis] + np.outer(bk[..., 0],
                                 ps[...].ravel())), axis=1)
                top = ak[nlevs-1, 1] + bk[nlevs-1, 1]*ps[...].ravel()

            elif self.varname == 'co':
                # grab the averaging kernel and reshape it
                avg_kernel = ncd.groups['PRODUCT'].groups['SUPPORT_DATA'].\
                    groups['DETAILED_RESULTS'].variables['column_averaging_kernel'][:]
                avg_kernel = np.reshape(avg_kernel, (nlocs, nlevs))

                # construct the pressure vertices array
                preslv = ncd.groups['PRODUCT'].groups['SUPPORT_DATA'].\
                    groups['DETAILED_RESULTS'].variables['pressure_levels'][:]
                preslv = np.reshape(preslv, (nlocs, nlevs))
                top = np.zeros(nlocs)

            # assemble presvertices with top vertice
            preslv = np.append(top[:, np.newaxis], preslv, axis=1)

            # scale the avk using AMF ratio and tropopause level for tropo column
            nlocf = len(lats[flg])
            scaleAK = np.ones((nlocf, nlevs), dtype=np.float32)

            if first:
                # add metadata variables
                self.outdata[('dateTime', 'MetaData')] = times[flg]
                self.outdata[('latitude', 'MetaData')] = lats[flg]
                self.outdata[('longitude', 'MetaData')] = lons[flg]
                self.outdata[('quality_assurance_value', 'MetaData')] = qa_value[flg]

                self.outdata[('averagingKernel', 'RetrievalAncillaryData')] = avg_kernel[flg]
                self.outdata[('pressureVertice', 'RetrievalAncillaryData')] = preslv[flg]

            else:
                self.outdata[('dateTime', 'MetaData')] = np.concatenate((
                    self.outdata[('dateTime', 'MetaData')], times[flg]))
                self.outdata[('latitude', 'MetaData')] = np.concatenate((
                    self.outdata[('latitude', 'MetaData')], lats[flg]))
                self.outdata[('longitude', 'MetaData')] = np.concatenate((
                    self.outdata[('longitude', 'MetaData')], lons[flg]))
                self.outdata[('quality_assurance_value', 'MetaData')] = np.concatenate((
                    self.outdata[('quality_assurance_value', 'MetaData')], qa_value[flg]))

                self.outdata[('averagingKernel', 'RetrievalAncillaryData')] = np.concatenate((
                    self.outdata[('averagingKernel', 'RetrievalAncillaryData')], avg_kernel[flg]))
                self.outdata[('pressureVertice', 'RetrievalAncillaryData')] = np.concatenate((
                    self.outdata[('pressureVertice', 'RetrievalAncillaryData')], preslv[flg]))

            for ncvar, iodavar in self.obsVar.items():

                if ncvar in ['nitrogendioxide_tropospheric_column',
                             'carbonmonoxide_total_column']:
                    data = ncd.groups['PRODUCT'].variables[ncvar][:].ravel()[flg]
                    err = ncd.groups['PRODUCT'].variables[ncvar+'_precision'][:].ravel()[flg]
                else:
                    data = ncd.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables[ncvar][:].ravel()[flg]
                    err = ncd.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables[ncvar+'_precision'][:].ravel()[flg]
                if first:
                    self.outdata[self.varDict[iodavar]['valKey']] = data
                    self.outdata[self.varDict[iodavar]['errKey']] = err
                    self.outdata[self.varDict[iodavar]['qcKey']] = qc_flag[flg]
                else:
                    self.outdata[self.varDict[iodavar]['valKey']] = np.concatenate(
                        (self.outdata[self.varDict[iodavar]['valKey']], data))
                    self.outdata[self.varDict[iodavar]['errKey']] = np.concatenate(
                        (self.outdata[self.varDict[iodavar]['errKey']], err))
                    self.outdata[self.varDict[iodavar]['qcKey']] = np.concatenate(
                        (self.outdata[self.varDict[iodavar]['qcKey']], qc_flag[flg]))

            first = False

        DimDict['Location'] = len(self.outdata[('dateTime', 'MetaData')])
        AttrData['Location'] = np.int32(DimDict['Location'])
        DimDict['Layer'] = nlevs
        AttrData['Layer'] = np.int32(DimDict['Layer'])
        DimDict['Vertice'] = nlevs + 1
        AttrData['Vertice'] = np.int32(DimDict['Vertice'])

        varname = 'pressureVertice'
        vkey = (varname, 'RetrievalAncillaryData')
        self.varAttrs[vkey]['coordinates'] = 'longitude latitude'
        self.varAttrs[vkey]['units'] = 'Pa'

        varname = 'averagingKernel'
        vkey = (varname, 'RetrievalAncillaryData')
        self.varAttrs[vkey]['coordinates'] = 'longitude latitude'
        self.varAttrs[vkey]['units'] = ''

# Check to determine if a point falls inside CONUS-13km NA [HChoi 6/23/2023]
# Code from lam_domaincheck_esg_c subroutine in LAMDomainCheck.interface.F90
def rrfs_domain_check(lat,lon):
  alpha = 0.112806427227431     # alpha parameter for ESG grid definition
  kappa = -0.350523595226368    # kappa paramater for ESG grid definition
  plat = 38.5                   # center point latitude of ESG grid (degrees)
  plon = -97.5                  # center point longitude of ESG grid (degrees)
  pazi = 0.0                    # azimuth angle for ESG grid definition
  npx = 396                     # number of grid points in x direction
  npy = 232                     # number of grid points in y direction
#  dx = 0.01349                 # grid spacing in degrees (3 km grid spacing)
#  dy = 0.01349                 # grid spacing in degrees (3 km grid spacing)
  delx = 0.0010127              # grid spacing of supergrid in map units
  dely = 0.0010176              # grid spacing of supergrid in map units
  failure = False               # failure code

  plat = plat * (m.pi / 180.0)  # convert to radians
  plon = plon * (m.pi / 180.0)  # convert to radians
  delx = delx * 2               # multiply by 2 to get actual grid resolution
  dely = dely * 2               # multiply by 2 to get actual grid resolution
  lat = lat * (m.pi / 180.0)    # convert to radians
  lon = lon * (m.pi / 180.0)    # convert to radians
#  print (delx)
#  sys.exit('STOP!!')
# gtoxm_ak_rr subroutine from Jim Purser's pesg.f90 code

  clat = m.cos(plat)
  slat = m.sin(plat)
  clon = m.cos(plon)
  slon = m.sin(plon)
  cazi = m.cos(pazi)
  sazi = m.sin(pazi)

  azirot = [[cazi,  sazi, 0],
            [-sazi, cazi, 0],
            [0,     0,    1]]
  prot = [[-slon,      clon,       0],
          [-slat*clon, -slat*slon, clat],
          [clat*clon,  clat*slon,  slat]]

  azirot = np.array(azirot)
  prot = np.array(prot)
  prot = np.matmul(prot,azirot)

  sla = m.sin(lat)
  cla = m.cos(lat)
  slo = m.sin(lon)
  clo = m.cos(lon)

  xe = [[cla*clo],
        [cla*slo],
        [sla]]

  xe = np.array(xe)
  xc = np.matmul(prot,xe)       # Do NOT use the transpose of prot here

  xc = np.array(xc)
  zp = float(xc[2]) + 1.0
  xs = xc[0:2]/zp

  xs = np.array(xs)

  s = kappa * (float(xs[0])*float(xs[0]) + float(xs[1])*float(xs[1]))
  sc = 1.0 - s
  if (abs(s) >= 1.0):
    failure = True
  xt = (2.0 * xs) / sc

  xm = xt       # Define xm, will set correct values below

  if (alpha > 0):
    ra = m.sqrt(alpha)
    razt = ra * xt[0]
    xm[0] = m.atan(razt) / ra
  elif (alpha < 0):
    ra = m.sqrt(-alpha)
    razt = ra * xt[0]
    if (abs(razt) >= 1.0):
      failure = True
    xm[0] = m.atanh(razt) / ra
  else:
    xm[0] = xt[0]

  if (alpha > 0):
    ra = m.sqrt(alpha)
    razt = ra * xt[1]
    xm[1] = m.atan(razt) / ra
  elif (alpha < 0):
    ra = m.sqrt(-alpha)
    razt = ra * xt[1]
    if (abs(razt) >= 1.0):
      failure = True
    xm[1] = m.atanh(razt) / ra
  else:
    xm[1] = xt[1]

  xm[0] = xm[0]/delx
  xm[1] = xm[1]/dely
 #print((xm))

# use xm to determine if point is good or bad
  if ((abs(xm[0])) < (npx/2)) and ((abs(xm[1])) < (npy/2)) and (failure == False):
  # Point is inside the ESG grid
    check = False
  else:
  # Point is outside the ESG grid
    check = True

  return(check)

def main():

    # get command line arguments
    parser = argparse.ArgumentParser(
        description=(
            'Reads TROPOMI NO2/CO netCDF files: official Copernicus product'
            'and converts into IODA formatted output files. Multiple'
            'files are able to be concatenated.')
    )

    required = parser.add_argument_group(title='required arguments')
    required.add_argument(
        '-i', '--input',
        help="path of TROPOMI L2 NO2/CO observation netCDF input file(s)",
        type=str, nargs='+', required=True)
    required.add_argument(
        '-o', '--output',
        help="path of IODA output file",
        type=str, required=True)
    required.add_argument(
        '-v', '--variable',
        help="name of varibale, available list: [no2, co]",
        type=str, required=True)
    required.add_argument(
        '-c', '--column',
        help="type of column: total or tropo",
        type=str, required=True)
    optional = parser.add_argument_group(title='optional arguments')
    optional.add_argument(
        '-q', '--qa_value',
        help="qa value used to preflag data that goes into file before QC"
        "default at 0.75 (no2) as suggested in the documentation. See:"
        "https://sentinel.esa.int/documents/247904/2474726/"
        "Sentinel-5P-Level-2-Product-User-Manual-Nitrogen-Dioxide.pdf section 8.6"
        "0.5 is suggested for co. See:"
        "https://sentinel.esa.int/documents/247904/2474726/"
        "Sentinel-5P-Level-2-Product-User-Manual-Carbon-Monoxide.pdf section 8.3",
        type=float, default=0.75)
    optional.add_argument(
        '-n', '--thin',
        help="percentage of random thinning from 0.0 to 1.0. Zero indicates"
        " no thinning is performed. (default: %(default)s)",
        type=float, default=0.0)

    args = parser.parse_args()

    if args.variable == "co":
        var_name = 'carbonmonoxide'
        if args.column == "tropo":
            print('CO is only available for total column, reset column to total', flush=1)
            args.column = 'total'
    elif args.variable == "no2":
        var_name = 'nitrogendioxide'

    if args.column == "tropo":

        obsVar = {
            var_name+'_tropospheric_column': var_name+'Column'
        }

        varDims = {
            var_name+'Column': ['Location']
        }

    elif args.column == "total":

        obsVar = {
            var_name+'_total_column': var_name+'Total'
        }

        varDims = {
            var_name+'Total': ['Location']
        }

    varDims['averagingKernel'] = ['Location', 'Layer']
    varDims['pressureVertice'] = ['Location', 'Vertice']

    # Read in the NO2 data
    var = tropomi(args.input, args.variable, args.column, args.qa_value, args.thin, obsVar)

    # setup the IODA writer
    writer = iconv.IodaWriter(args.output, locationKeyList, DimDict)

    # write everything out
    writer.BuildIoda(var.outdata, varDims, var.varAttrs, AttrData)


if __name__ == '__main__':
    main()
