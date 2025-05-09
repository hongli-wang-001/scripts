#!/usr/bin/env python3
#
# (C) Copyright 2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

import sys
import argparse
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
from pathlib import Path

IODA_CONV_PATH = Path(__file__).parent/"../lib/pyiodaconv"
if not IODA_CONV_PATH.is_dir():
    IODA_CONV_PATH = Path(__file__).parent/'..'/'lib-python'
sys.path.append(str(IODA_CONV_PATH.resolve()))

import ioda_conv_ncio as iconv
from orddicts import DefaultOrderedDict

vName = {
    'A': "aerosol_optical_depth_4",
}

locationKeyList = [
    ("latitude", "float"),
    ("longitude", "float"),
    ("datetime", "string")
]

AttrData = {}


class AOD(object):

    def __init__(self, filename, method, mask, thin, writer):
        self.filename = filename
        self.method = method
        self.mask = mask
        self.thin = thin
        self.data = DefaultOrderedDict(lambda: DefaultOrderedDict(dict))
        self.writer = writer
        self._read()

    def _read(self):
        ncd = nc.Dataset(self.filename)
        lons = ncd.variables['Longitude'][:].ravel()
        lats = ncd.variables['Latitude'][:].ravel()
        vals = ncd.variables['AOD550'][:].ravel()
        errs = ncd.variables['Residual'][:].ravel()
        qcpath = ncd.variables['QCPath'][:].ravel()
        qcall = ncd.variables['QCAll'][:].ravel()
        if self.mask == "maskout":
            mask = np.logical_not(vals.mask)
            vals = vals[mask]
            lons = lons[mask]
            lats = lats[mask]
            errs = errs[mask]
            qcpath = qcpath[mask]
            qcall = qcall[mask]
        # get global attributes
        gatts = {attr: getattr(ncd, attr) for attr in ncd.ncattrs()}
        base_datetime = gatts["time_coverage_end"]
        self.satellite = gatts["satellite_name"]
        self.sensor = gatts["instrument_name"]
        ncd.close()

        valKey = vName['A'], self.writer.OvalName()
        errKey = vName['A'], self.writer.OerrName()
        qcKey = vName['A'], self.writer.OqcName()

        # apply thinning mask
        if self.thin > 0.0:
            mask_thin = np.random.uniform(size=len(lons)) > self.thin
            lons = lons[mask_thin]
            lats = lats[mask_thin]
            vals = vals[mask_thin]
            errs = errs[mask_thin]
            qcpath = qcpath[mask_thin]
            qcall = qcall[mask_thin]

        # set variable keys
        sfcKey = "surface_type", "MetaData"
        szaKey = "sol_zenith_angle", "MetaData"
        saaKey = "sol_azimuth_angle", "MetaData"
        mduKey = "modis_deep_blue_flag", "MetaData"
        biaKey = "aerosol_optical_depth_4", "KnownObsBias"

        # defined surface type, uncertainty, and bias array
        sfctyp = 0*qcall
        uncertainty = 0.0*errs
        uncertainty1 = 0.0*errs
        uncertainty2 = 0.0*errs
        bias = 0.0*errs
        bias1 = 0.0*errs
        bias2 = 0.0*errs

        if self.method == "nesdis":
            # Case of water high quality
            uncertainty = 0.00784394 + 0.219923*vals
            bias = 0.0151799 + 0.0767385*vals
            # case of bright land high quality
            uncertainty1 = 0.0550472 + 0.299558*vals
            bias1 = -0.0107621 + 0.150480*vals
            # case of dark land high quality
            uncertainty2 = 0.111431 + 0.128699*vals
            bias2 = -0.0138969 + 0.157877*vals

        for i in range(len(lons)):

            # convert byte to integer
            sfctyp[i] = int.from_bytes(qcpath[i], byteorder='big')
            if self.method == "nesdis":
                if sfctyp[i] == 1:   # case of bright land high quality
                    bias[i] = bias1[i]
                    uncertainty[i] = uncertainty1[i]
                else:   # case of dark land high quality
                    bias[i] = bias2[i]
                    uncertainty[i] = uncertainty2[i]
                errs[i] = uncertainty[i]

            locKey = lats[i], lons[i], base_datetime
            self.data[0][locKey][valKey] = vals[i]
            self.data[0][locKey][errKey] = errs[i]
            self.data[0][locKey][qcKey] = qcall[i]
            self.data[0][locKey][biaKey] = bias[i]

            # solar zenith angle (sza) is set all 0 for test
            # solar azimuth angle (saa) is set  all 0 for test
            # modis_deep_blue_flag (mdu)is set all 0 for test
            self.data[0][locKey][szaKey] = 0.0
            self.data[0][locKey][saaKey] = 0.0
            self.data[0][locKey][sfcKey] = sfctyp[i]
            self.data[0][locKey][mduKey] = 0

            # write global attributes out
            if self.satellite == 'NPP':
                self.satellite = "suomi_npp"
            if self.sensor == 'VIIRS':
                self.sensor = "v.viirs-m_npp"

            AttrData["observation_type"] = "AOD"
            AttrData["satellite"] = self.satellite
            AttrData["sensor"] = self.sensor
            AttrData['date_time_string'] = base_datetime


def main():

    parser = argparse.ArgumentParser(
        description=('Read VIIRS aerosol optical depth file(s) and Converter'
                     ' of native NetCDF format for observations of optical'
                     ' depth from VIIRS AOD550 to IODA netCDF format.')
    )
    parser.add_argument('-i', '--input',
                        help="name of viirs aod input file(s)",
                        type=str, required=True)
    parser.add_argument('-o', '--output',
                        help="name of ioda output file",
                        type=str, required=True)
    optional = parser.add_argument_group(title='optional arguments')
    optional.add_argument(
        '-m', '--method',
        help="calculation bias/error method: nesdis/default, default=none",
        type=str, required=True)
    optional.add_argument(
        '-k', '--mask',
        help="maskout missing values: maskout/default, default=none",
        type=str, required=True)
    optional.add_argument(
        '-t', '--thin',
        help="percentage of random thinning, from 0.0 to 1.0. Zero indicates"
             " no thinning is performed. (default: %(default)s)",
        type=float, default=0.0)

    args = parser.parse_args()

    writer = iconv.NcWriter(args.output, locationKeyList)

    # Read in the profiles
    aod = AOD(args.input, args.method, args.mask, args.thin, writer)

    (ObsVars, LocMdata, VarMdata) = writer.ExtractObsData(aod.data)

    # set constants for the four variables
    VarMdata['frequency'] = np.full(1, 5.401666e+14, dtype='f4')
    VarMdata['polarization'] = np.full(1, 1, dtype='i4')
    VarMdata['wavenumber'] = np.full(1, 18161.61, dtype='f4')
    VarMdata['sensor_channel'] = np.full(1, 4, dtype='i4')

    writer.BuildNetcdf(ObsVars, LocMdata, VarMdata, AttrData)


if __name__ == '__main__':
    main()
