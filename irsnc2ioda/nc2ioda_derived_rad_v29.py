# nc2ioda_derived_rad_v29.py

import sys
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
import netCDF4 as nc

from pyiodaconv.ioda_conv_engines import MetaDataName, OvalName, OerrName
from pyiodaconv.def_jedi_utils import set_metadata_attributes, set_obspace_attributes
import pyiodaconv.ioda_conv_engines as iconv
from pyiodaconv.orddicts import DefaultOrderedDict


# Constants used globally
from netCDF4 import default_fillvals

missing_vals_float = default_fillvals['f4']
missing_vals_int   = default_fillvals['i4']

WMO_SAT_ID = 72
WMO_SEN_ID = 212

def read_out_channel_info(nc_path, spectral_group_name='lwir'):
    """
    Read n_channels and out_channels from the netCDF file at nc_path,
    using the derived group under data → spectral_group_name → derived.

    Args:
        nc_path (str): path to the netCDF file.
        spectral_group_name (str): name of the spectral group, e.g. 'lwir' or 'mwir'.

    Returns:
        n_channels (int): size of the n_channels_out dimension.
        out_channels (np.ndarray of int32): array of required_channels (1‑based).

    Raises:
        KeyError: if expected groups/variables/dimensions are not found.
    """
    with Dataset(nc_path, mode='r') as ncfile:
        # Navigate to the derived group
        try:
            derived = ncfile.groups['data'].groups[spectral_group_name].groups['derived']
        except KeyError as e:
            raise KeyError(f"Could not find group data/{spectral_group_name}/derived in file {nc_path}: {e}")

        # Read the dimension
        if 'n_channels_out' not in derived.dimensions:
            raise KeyError(f"Dimension 'n_channels_out' not found in derived group for spectral group '{spectral_group_name}'")
        n_channels = derived.dimensions['n_channels_out'].size

        # Read the required_channels variable
        if 'required_channels' not in derived.variables:
            raise KeyError(f"Variable 'required_channels' not found in derived group for spectral group '{spectral_group_name}'")
        required_channels_var = derived.variables['required_channels']
        out_channels = required_channels_var[:].astype(np.int32)

        # Optional consistency check
        if out_channels.shape[0] != n_channels:
            print(f"[WARN] out_channels length {out_channels.shape[0]} != n_channels {n_channels}")

        return n_channels, out_channels

def compute_sensor_angles_general(
    sat_lat_deg, sat_lon_deg,
    site_lat_deg, site_lon_deg,
    Re=6371.0, H=35786.0
):
    """
    Compute sensor zenith and azimuth angles from satellite to ground points.

    Parameters:
    - sat_lat_deg: scalar, satellite latitude in degrees
    - sat_lon_deg: scalar, satellite longitude in degrees
    - site_lat_deg: array-like, ground point latitude(s) in degrees
    - site_lon_deg: array-like, ground point longitude(s) in degrees
    - Re: Earth radius in km (default 6371 km)
    - H: Satellite height above Earth's surface in km (default GEO 35786 km)

    Returns:
    - zenith_angle_deg: array-like, sensor zenith angle(s) in degrees
    - azimuth_deg: array-like, sensor azimuth angle(s) in degrees (clockwise from North)
    """

    # Convert all angles to radians
    sat_lat = np.radians(sat_lat_deg)
    sat_lon = np.radians(sat_lon_deg)
    site_lat = np.radians(np.asarray(site_lat_deg))
    site_lon = np.radians(np.asarray(site_lon_deg))

    # Satellite position in ECEF
    r_sat = Re + H
    x_sat = r_sat * np.cos(sat_lat) * np.cos(sat_lon)
    y_sat = r_sat * np.cos(sat_lat) * np.sin(sat_lon)
    z_sat = r_sat * np.sin(sat_lat)

    # Ground point position in ECEF
    x_site = Re * np.cos(site_lat) * np.cos(site_lon)
    y_site = Re * np.cos(site_lat) * np.sin(site_lon)
    z_site = Re * np.sin(site_lat)

    # Vector from site to satellite
    dx = x_sat - x_site
    dy = y_sat - y_site
    dz = z_sat - z_site

    # Distance and norms
    r_norm = np.sqrt(dx**2 + dy**2 + dz**2)
    site_norm = np.sqrt(x_site**2 + y_site**2 + z_site**2)

    # Zenith angle: angle between vector to satellite and local zenith (site vector)
    dot = dx * x_site + dy * y_site + dz * z_site
    cos_zenith = dot / (r_norm * site_norm)
    cos_zenith = np.clip(cos_zenith, -1, 1)  # numerical stability
    zenith_angle_rad = np.arccos(cos_zenith)
    zenith_angle_deg = np.degrees(zenith_angle_rad)

    # Local ENU vectors at ground site
    sin_lat = np.sin(site_lat)
    cos_lat = np.cos(site_lat)
    sin_lon = np.sin(site_lon)
    cos_lon = np.cos(site_lon)

    east = np.array([-sin_lon, cos_lon, np.zeros_like(site_lat)])
    north = np.array([
        -sin_lat * cos_lon,
        -sin_lat * sin_lon,
        cos_lat
    ])

    # Reshape for broadcasting if needed
    # east and north have shape (3, N) where N is number of points
    # dx, dy, dz have shape (N, )
    # So dot products computed element-wise over N

    # Project vector to satellite onto east and north
    e = dx * east[0] + dy * east[1] + dz * east[2]
    n = dx * north[0] + dy * north[1] + dz * north[2]

    azimuth_rad = np.arctan2(e, n)
    azimuth_deg = (np.degrees(azimuth_rad) + 360) % 360

    return zenith_angle_deg, azimuth_deg

def compute_sensor_angles_array(sat_lon_deg, site_lat_deg, site_lon_deg, H=35786.0, Re=6371.0):
    # H=35786.0
    # Convert inputs to numpy arrays
    site_lat = np.radians(np.asarray(site_lat_deg))
    site_lon = np.radians(np.asarray(site_lon_deg))
    sat_lon = np.radians(sat_lon_deg)

    # Satellite position (fixed, scalar)
    x_sat = (Re + H) * np.cos(sat_lon)
    y_sat = (Re + H) * np.sin(sat_lon)
    z_sat = 0.0

    # Ground point positions in ECEF
    x_site = Re * np.cos(site_lat) * np.cos(site_lon)
    y_site = Re * np.cos(site_lat) * np.sin(site_lon)
    z_site = Re * np.sin(site_lat)

    # Vector from site to satellite
    dx = x_sat - x_site
    dy = y_sat - y_site
    dz = z_sat - z_site

    # Range vector magnitude
    r_norm = np.sqrt(dx**2 + dy**2 + dz**2)
    zenith_norm = np.sqrt(x_site**2 + y_site**2 + z_site**2)

    # Zenith angle
    dot_product = dx * x_site + dy * y_site + dz * z_site
    cos_zenith_angle = dot_product / (r_norm * zenith_norm)
    zenith_angle_rad = np.arccos(np.clip(cos_zenith_angle, -1.0, 1.0))
    zenith_angle_deg = np.degrees(zenith_angle_rad)

    # Compute ENU (east, north) vectors for azimuth
    sin_lat = np.sin(site_lat)
    cos_lat = np.cos(site_lat)
    sin_lon = np.sin(site_lon)
    cos_lon = np.cos(site_lon)

    # ENU basis
    east_x = -sin_lon
    east_y = cos_lon
    east_z = 0

    north_x = -sin_lat * cos_lon
    north_y = -sin_lat * sin_lon
    north_z = cos_lat

    # Project satellite vector into ENU
    e = dx * east_x + dy * east_y + dz * east_z
    n = dx * north_x + dy * north_y + dz * north_z

    azimuth_rad = np.arctan2(e, n)
    azimuth_deg = (np.degrees(azimuth_rad) + 360) % 360

    return zenith_angle_deg, azimuth_deg

def compute_sensor_angles(sat_lon_deg, site_lat_deg, site_lon_deg, Re=6371.0, H=35786.0):
    # Convert degrees to radians
    sat_lon = np.radians(sat_lon_deg)
    site_lat = np.radians(site_lat_deg)
    site_lon = np.radians(site_lon_deg)

    # Satellite position in ECEF (Earth-centered Earth-fixed), assuming it's above equator
    x_sat = (Re + H) * np.cos(sat_lon)
    y_sat = (Re + H) * np.sin(sat_lon)
    z_sat = 0  # Equator

    # Ground point position in ECEF
    x_site = Re * np.cos(site_lat) * np.cos(site_lon)
    y_site = Re * np.cos(site_lat) * np.sin(site_lon)
    z_site = Re * np.sin(site_lat)

    # Vector from site to satellite
    dx = x_sat - x_site
    dy = y_sat - y_site
    dz = z_sat - z_site

    # Range vector (satellite - ground)
    r = np.array([dx, dy, dz])
    r_norm = np.linalg.norm(r)

    # Local zenith vector (same as position vector for spherical Earth)
    zenith = np.array([x_site, y_site, z_site])
    zenith_norm = np.linalg.norm(zenith)

    # Zenith angle = angle between r and zenith
    cos_zenith_angle = np.dot(r, zenith) / (r_norm * zenith_norm)
    zenith_angle_rad = np.arccos(cos_zenith_angle)
    zenith_angle_deg = np.degrees(zenith_angle_rad)

    # Azimuth angle: angle between projection of r onto local horizontal plane and north
    # Compute ENU (east-north-up) basis at the site
    east = np.array([-np.sin(site_lon), np.cos(site_lon), 0])
    north = np.array([-np.sin(site_lat) * np.cos(site_lon),
                      -np.sin(site_lat) * np.sin(site_lon),
                       np.cos(site_lat)])
    up = np.array([np.cos(site_lat) * np.cos(site_lon),
                   np.cos(site_lat) * np.sin(site_lon),
                   np.sin(site_lat)])

    # Project r into ENU frame
    e = np.dot(r, east)
    n = np.dot(r, north)

    azimuth_rad = np.arctan2(e, n)
    azimuth_deg = (np.degrees(azimuth_rad) + 360) % 360

    return zenith_angle_deg, azimuth_deg

def apply_scaling(var):
    """Apply scale_factor and add_offset if they exist; return float32 array."""
    data = var[:].astype(np.float32)
    scale = getattr(var, 'scale_factor', 1.0)
    offset = getattr(var, 'add_offset', 0.0)
    arr = data.reshape(-1)
    print(f"[INFO] scale_factor '{scale}'")
    print(f"[INFO] add_offset '{offset}'")
    # Every 25th values
    print(f"  Sample every 25th of var: {arr[::25]}")
    return data * scale + offset


def read_derived_rad(file_path):
    """Read data & metadata from derived_rad netCDF, compute datetime array, select channels."""
    ncfile = Dataset(file_path, mode='r')
    data_group = ncfile.groups['data']

    dwell_row = data_group.dimensions['dwell_row'].size
    dwell_column = data_group.dimensions['dwell_column'].size
    nloc = dwell_row * dwell_column
    print(f"[INFO] dwell_row={dwell_row}, dwell_column={dwell_column}, nloc={nloc}")

    # Parse global attribute time_coverage_start into datetime
    tcs = getattr(ncfile, 'time_coverage_start', None)
    if tcs is None:
        raise ValueError("Global attribute 'time_coverage_start' not found")
    if len(tcs) == 14:
        fmt = "%Y%m%d%H%M%S"
    elif len(tcs) == 12:
        fmt = "%Y%m%d%H%M"
    else:
        raise ValueError(f"Unrecognized format for time_coverage_start: {tcs}")
    dt0 = datetime.strptime(tcs, fmt)
    print(f"[INFO] time_coverage_start as datetime: {dt0}")

    # Seconds since UNIX epoch (1970‑01‑01)
    epoch = datetime(1970, 1, 1)
    seconds_since_epoch = int((dt0 - epoch).total_seconds())
    datetime_array = np.full((nloc,), seconds_since_epoch, dtype=np.int64)

    # Load spatial vars
    spatial_vars = [
        'latitude', 'longitude',
        'satellite_azimuth_angle', 'satellite_zenith_angle',
        'solar_azimuth_angle', 'solar_zenith_angle',
        'cloud_signal', 'cloud_fraction'
    ]
    data = {'data': {}}
    for var_name in spatial_vars:
        if var_name in data_group.variables:
            var = data_group.variables[var_name]
            ##var_raw = var[:].reshape(-1)
            ##arr = apply_scaling(var).reshape(-1)
            #arr = var[:].reshape(-1)
            # Thinning data
            thinned = var[0::1, 0::1]  # Start at 3rd row/col, step by 5
            arr = thinned.reshape(-1)
            data['data'][var_name] = arr
            print(f"[INFO] Scaled + reshaped '{var_name}'")
            # Print a small sample of values
            # First 5 values
            print(f"  Sample first 5 of {var_name}: {arr[:5]}")
            # Every 25th values
            #print(f"  Sample every 25th of raw {var_name}: {var_raw[::25]}")
            #print(f"  Sample every 25th of scaled {var_name}: {arr[::25]}")
            # Or random sample
            # import numpy as np
            # idx = np.random.choice(arr.size, size=5, replace=False)
            # print(f"  Random sample of {var_name}: {arr[idx]}")
            # Store shape once to calculate nloc after loop
            if var_name == 'latitude':
                thinned_shape = thinned.shape  # (rows, cols)
            print(f"[INFO] Thinned + reshaped '{var_name}' to shape {thinned.shape}")

        else:
            print(f"[WARN] Spatial var '{var_name}' not found in file")

    nloc = thinned_shape[0] * thinned_shape[1]
    datetime_array = np.full((nloc,), seconds_since_epoch, dtype=np.int64)
    print(f"[INFO] Updated nloc after thinning: {nloc}")

    # Load metadata vars
    meta_vars = ['dwell_type', 'dwell_number', 'stroke_direction', 'dust_warning']
    for var_name in meta_vars:
        if var_name in data_group.variables:
            arr = data_group.variables[var_name][:].reshape(-1)
            data['data'][var_name] = arr
            print(f"[INFO] Loaded metadata '{var_name}'")
        else:
            print(f"[WARN] Metadata var '{var_name}' not found")

    # LWIR wavenumber and related info
    lwir_group = data_group.groups['lwir']
    wavenumber = lwir_group.variables['wavenumber'][:] * 100.0  # conversion to m⁻¹
    start_wavenumber = lwir_group.variables['start_wavenumber'][:]
    end_wavenumber   = lwir_group.variables['end_wavenumber'][:]
    wavenumber_step  = lwir_group.variables['wavenumber_step'][:]
    spatial_sampling_distance = lwir_group.variables['spatial_sampling_distance'][:]
    print(f"[INFO] LWIR wavenumber shape: {wavenumber.shape}")
    print(f"[INFO] LWIR start_wavenumber: {start_wavenumber}")
    print(f"[INFO] LWIR end_wavenumber: {end_wavenumber}")
    print(f"[INFO] LWIR wavenumber_step: {wavenumber_step}")
    print(f"[INFO] LWIR :spatial_sampling_distance {spatial_sampling_distance}")

    # Determine valid channel indices (0‑based)
    req_zero = out_channels_lw - 1
    valid_indices = req_zero[(req_zero >= 0) & (req_zero < wavenumber.size)]
    selected_wavenumbers = wavenumber[valid_indices]
    print(f"[INFO] Selected {len(selected_wavenumbers)} LWIR wavenumbers")

    data['lwir'] = {
        'selected_wavenumbers': selected_wavenumbers,
        'selected_channel_numbers': out_channels_lw
    }

    # Read reconstructed radiance and mask/fix invalid values
    derived = lwir_group.groups['derived']
    #rec_rad = derived.variables['reconstructed_radiance'][:]
    #Thinning data
    rec_rad = derived.variables['reconstructed_radiance'][0::1, 0::1, :]
    print(f"[INFO] rec_rad thinning raw min value: {np.min(rec_rad)}")
    print(f"[INFO] rec_rad thinning raw max value: {np.max(rec_rad)}")
    rec_rad[rec_rad < 0.0] = missing_vals_float
    print(f"[INFO] rec_rad filter min value: {np.min(rec_rad)}")
    print(f"[INFO] rec_rad filter max value: {np.max(rec_rad)}")
    # Mask invalid values
    #valid_mask = (rec_rad >= 0.0) & (rec_rad != missing_vals_float)
    valid_mask = (rec_rad >= 0.0) & (rec_rad != missing_vals_float) & ~np.isnan(rec_rad)

    # Use masked array so that invalid values are ignored
    masked = np.ma.array(rec_rad, mask=~valid_mask)
    # Compute stats along the third dimension (axis = 2)
    mean_per_channel = masked.mean(axis=(0,1))       # shape → (dim3,)
    median_per_channel = np.ma.median(masked, axis=(0,1))
    std_per_channel = masked.std(axis=(0,1))         # default ddof=0; use ddof=1 for sample std
    min_per_channel    = masked.min(axis=(0, 1))
    max_per_channel    = masked.max(axis=(0, 1))

    # Print / inspect
    for i in range(rec_rad.shape[2]):
        # First 5 values
        valid_obs_num = masked[:, :, i].count()
        print(f"  Sample first 5 of reconstructed rad: {masked[::1,::1,i]}")
        print(
            f"Channel {i} {out_channels_lw[i]}: "
            f"valid_obs_num={valid_obs_num},"
            f"mean={mean_per_channel[i]:.6f}, "
            f"median={median_per_channel[i]:.6f}, "
            f"std={std_per_channel[i]:.6f}, "
            f"min={min_per_channel[i]:.6f}, "
            f"max={max_per_channel[i]:.6f}"
    )

    #for i in range(rec_rad.shape[2]):
    #    print(f"Channel {i} {out_channels_lw[i]}: mean={mean_per_channel[i]:.4f}, "
    #          f"median={median_per_channel[i]:.4f}, std={std_per_channel[i]:.4f}")

    valid_values = rec_rad[valid_mask]

    num_valid = np.count_nonzero(valid_mask)
    total_values = rec_rad.size
    print(f"[INFO] Total radiance values: {total_values}, valid: {num_valid}, fraction: {num_valid/total_values:.4f}")

    if num_valid > 0:
        mean_val = np.mean(valid_values)
        median_val = np.median(valid_values)
        std_val = np.std(valid_values)   # default population std; use ddof=1 for sample

        print(f"[INFO] Mean of valid radiance: {mean_val:.4f}")
        print(f"[INFO] Median of valid radiance: {median_val:.4f}")
        print(f"[INFO] Std dev of valid radiance: {std_val:.4f}")
    else:
        print("[WARN] No valid radiance values to compute statistics")


    # Apply unit conversion and mask again
    # rec_rad = rec_rad * wavenumber_step
    #rec_rad[rec_rad < 0.0] = missing_vals_float
    #rec_rad[rec_rad < 0.0] = 0.5600E-03
    data['lwir_derived'] = {
        'reconstructed_radiance': rec_rad
    }
    print(f"[INFO] Loaded reconstructed_radiance with shape {rec_rad.shape}")

    ncfile.close()
    return data, nloc, datetime_array


def write_ioda(data, nloc, datetime_array, output_file):
    """Build and write IODA file using pyiodaconv, using your original structure."""
    channel_numbers = np.arange(1, n_channels + 1, dtype=np.int32)

    # Flatten radiance
    rec_rad = data['lwir_derived']['reconstructed_radiance']
    radiance_flat = rec_rad.reshape(nloc, n_channels).astype(np.float32)
    print(f"[INFO] radiance_flat shape: {radiance_flat.shape}")

    md = MetaDataName()
    oval = OvalName()
    oerr = OerrName()

    # Prepare IODA data dict
    ioda_data = {}

    # Metadata variables
    ioda_data[('latitude', md)] = data['data']['latitude'].astype(np.float32)
    ioda_data[('longitude', md)] = data['data']['longitude'].astype(np.float32)
    ioda_data[('dateTime', md)] = datetime_array

    # Channel metadata
    ioda_data[('sensorChannelNumber', md)] = out_channels_lw
    ioda_data[('sensorCentralWavenumber', md)] = data['lwir']['selected_wavenumbers'].astype(np.float32)

    # Compute angles for MTG-S1
    zenith, azimuth = compute_sensor_angles_array(
        sat_lon_deg=-3.4,  
        site_lat_deg=data['data']['latitude'].astype(np.float32),
        site_lon_deg=data['data']['longitude'].astype(np.float32)
    )
    # Angular & other metadata per location
    #ioda_data[('sensorAzimuthAngle', md)] = data['data']['satellite_azimuth_angle'].astype(np.float32)
    #ioda_data[('sensorZenithAngle', md)] = data['data']['satellite_zenith_angle'].astype(np.float32)
    ioda_data[('sensorAzimuthAngle', md)] = azimuth
    ioda_data[('sensorZenithAngle', md)] = zenith

    ioda_data[('solarAzimuthAngle', md)] = data['data']['solar_azimuth_angle'].astype(np.float32)
    ioda_data[('solarZenithAngle', md)] = data['data']['solar_zenith_angle'].astype(np.float32)
    ioda_data[('cloudAmount', md)] = data['data']['cloud_fraction'].astype(np.float32)

    # Identifiers
    ioda_data[('satelliteIdentifier', md)] = np.full((nloc,), WMO_SAT_ID, dtype=np.int32)
    ioda_data[('instrumentIdentifier', md)] = np.full((nloc,), WMO_SEN_ID, dtype=np.int32)

    # Observations and error
    ioda_data[('radiance', oval)] = radiance_flat
    ioda_data[('radiance', oerr)] = np.full((nloc, n_channels), 1.0, dtype=np.float32)

    # Dimension dictionary
    DimDict = {
        'Location': np.int32(nloc),
        'Channel': channel_numbers
    }

    # Variable → dimensions mapping
    varDims = {
        'radiance': ['Location', 'Channel'],
        'sensorChannelNumber': ['Channel'],
        'sensorCentralWavenumber': ['Channel']
    }

    # Prepare variable attributes, location key list, etc.
    varAttrs = DefaultOrderedDict(lambda: DefaultOrderedDict(dict))

    locationKeyList = [
        ("latitude", "float", "degrees_north"),
        ("longitude", "float", "degrees_east"),
        ("dateTime", "long", "seconds since 1970-01-01T00:00:00Z"),
        ("sensorCentralWavenumber", "float", "m-1"),
        ("sensorChannelNumber", "integer", ""),
        ("sensorAzimuthAngle", "float", "degree"),
        ("sensorZenithAngle", "float", "degree"),
        ("solarAzimuthAngle", "float", "degree"),
        ("solarZenithAngle", "float", "degree")
    ]
    #missing_vals_float = np.finfo(np.float32).max
    #missing_vals_float = 9.96921e+36
    #missing_vals_int = np.iinfo(np.int32).max

    # Populate varAttrs for metadata keys
    for (key, dtype, units) in locationKeyList:
        attr = varAttrs[(key, md)]
        attr['units'] = units
        if key == 'sensorChannelNumber':
            attr['_FillValue'] = missing_vals_int
            attr['long_name'] = 'Sensor Channel Number'
        elif key == 'sensorCentralWavenumber':
            attr['_FillValue'] = missing_vals_float
            attr['long_name'] = 'Sensor Central Wavenumber'
        else:
            attr['_FillValue'] = missing_vals_float

    # Channel variable attributes
    varAttrs['Channel'] = {
        'long_name': 'Channel index',
        'units': '',
        '_FillValue': missing_vals_int
    }

    # Radiance variable attributes
    varAttrs[('radiance', iconv.OvalName())] = {
        '_FillValue': missing_vals_float,
        'units': 'W m-2 sr-1 m',
        'long_name': 'LWIR Reconstructed Radiance'
    }

    # Global attributes
    GlobalAttrs = {
        'sensor': WMO_SEN_ID,
        'platform': WMO_SAT_ID,
        'sensorCommonName': 'IRS',
        'platformCommonName': 'MTG-S1',
        'processingLevel': 'Level-1',
        'datetimeReference': datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    }

    # Write IODA
    writer = iconv.IodaWriter(output_file, locationKeyList, DimDict)
    set_metadata_attributes(varAttrs)
    set_obspace_attributes(varAttrs)
    writer.BuildIoda(ioda_data, varDims, varAttrs, GlobalAttrs)

    print(f"[INFO] Wrote IODA file: {output_file}")

def write_ioda_2(data, nloc, datetime_array, output_file):
    """Build and write IODA file using pyiodaconv, filtering out invalid lat/lon."""
    channel_numbers = np.arange(1, n_channels + 1, dtype=np.int32)

    # Extract latitude and longitude
    lat = data['data']['latitude'].astype(np.float32)
    lon = data['data']['longitude'].astype(np.float32)

    # Define valid mask: exclude NaNs and out-of-bounds values
    spatial_valid_mask = np.isfinite(lat) & np.isfinite(lon) & \
                 (lat >= -90.0) & (lat <= 90.0) & \
                 (lon >= -180.0) & (lon <= 180.0)

    # Define valid mask: exclude NaNs and out-of-bounds values
    #spatial_valid_mask = np.isfinite(lat) & np.isfinite(lon) & \
    #             (lat >= 20.0) & (lat <= 75.0) & \
    #             (lon >= -45.0) & (lon <= 50.0)

    # Flatten radiance
    rec_rad = data['lwir_derived']['reconstructed_radiance']
    radiance_flat = rec_rad.reshape(nloc, n_channels).astype(np.float32)

    # Define radiance-valid mask:
    # True if at least one channel in the observation is not missing
    missing_val = missing_vals_float
    radiance_valid_mask = np.any(
        (radiance_flat != missing_val) & np.isfinite(radiance_flat),
        axis=1
    )

    # Combine both masks
    valid_mask = spatial_valid_mask & radiance_valid_mask

    # Filter data using valid mask
    lat = lat[valid_mask]
    lon = lon[valid_mask]
    datetime_array = datetime_array[valid_mask]
    nloc_valid = len(lat)

    if nloc_valid == 0:
        print("[WARNING] No valid observations found — exiting.")
        sys.exit(0)

    # Flatten radiance and filter
    #rec_rad = data['lwir_derived']['reconstructed_radiance']
    #radiance_flat = rec_rad.reshape(nloc, n_channels).astype(np.float32)
    radiance_flat = radiance_flat[valid_mask]

    # Filter other location-dependent variables
    solar_azimuth = data['data']['solar_azimuth_angle'].astype(np.float32)[valid_mask]
    solar_zenith = data['data']['solar_zenith_angle'].astype(np.float32)[valid_mask]
    cloud_fraction = data['data']['cloud_fraction'].astype(np.float32)[valid_mask]

    # Recompute sensor angles with filtered lat/lon
    zenith, azimuth = compute_sensor_angles_array(
        sat_lon_deg=-3.4,
        site_lat_deg=lat,
        site_lon_deg=lon
    )

    md = MetaDataName()
    oval = OvalName()
    oerr = OerrName()

    # Prepare IODA data dict
    ioda_data = {}

    # Metadata variables (filtered)
    ioda_data[('latitude', md)] = lat
    ioda_data[('longitude', md)] = lon
    ioda_data[('dateTime', md)] = datetime_array

    # Channel metadata (unchanged)
    ioda_data[('sensorChannelNumber', md)] = out_channels_lw
    ioda_data[('sensorCentralWavenumber', md)] = data['lwir']['selected_wavenumbers'].astype(np.float32)

    # Angular & other metadata per location (filtered)
    ioda_data[('sensorAzimuthAngle', md)] = azimuth
    ioda_data[('sensorZenithAngle', md)] = zenith
    ioda_data[('solarAzimuthAngle', md)] = solar_azimuth
    ioda_data[('solarZenithAngle', md)] = solar_zenith
    ioda_data[('cloudAmount', md)] = cloud_fraction
    # Add sensorScanPosition
    # First, pre-allocate the array with a fill value (optional if assigning right after)
    sensor_scan_position = np.full((nloc_valid,), 0, dtype=np.int32)  # or use a meaningful default if needed

    # Then assign actual values by converting zenith to int (rounding down)
    sensor_scan_position = zenith.astype(np.int32)

    # Add to IODA data
    ioda_data[('sensorScanPosition', md)] = sensor_scan_position

    # Add sensorViewAngle
    ioda_data[('sensorViewAngle', md)] = ioda_data[('sensorZenithAngle', md)]

    # Identifiers (filtered size)
    ioda_data[('satelliteIdentifier', md)] = np.full((nloc_valid,), WMO_SAT_ID, dtype=np.int32)
    ioda_data[('instrumentIdentifier', md)] = np.full((nloc_valid,), WMO_SEN_ID, dtype=np.int32)

    # Observations and error (filtered)
    ioda_data[('radiance', oval)] = radiance_flat
    ioda_data[('radiance', oerr)] = np.full((nloc_valid, n_channels), 1.0, dtype=np.float32)

    # Dimension dictionary (updated to filtered nloc)
    DimDict = {
        'Location': np.int32(nloc_valid),
        'Channel': channel_numbers
    }

    # Variable → dimensions mapping
    varDims = {
        'radiance': ['Location', 'Channel'],
        'sensorChannelNumber': ['Channel'],
        'sensorCentralWavenumber': ['Channel']
    }

    # Prepare variable attributes, location key list, etc.
    varAttrs = DefaultOrderedDict(lambda: DefaultOrderedDict(dict))

    locationKeyList = [
        ("latitude", "float", "degrees_north"),
        ("longitude", "float", "degrees_east"),
        ("dateTime", "long", "seconds since 1970-01-01T00:00:00Z"),
        ("sensorCentralWavenumber", "float", "m-1"),
        ("sensorChannelNumber", "integer", ""),
        ("sensorAzimuthAngle", "float", "degree"),
        ("sensorZenithAngle", "float", "degree"),
        ("sensorScanPosition", "integer", ""),
        ("sensorViewAngle", "float", "degree"),
        ("solarAzimuthAngle", "float", "degree"),
        ("solarZenithAngle", "float", "degree")
    ]
    #missing_vals_float = np.finfo(np.float32).max
    #missing_vals_float = 9.96921e+36
    #missing_vals_int = np.iinfo(np.int32).max

    # Populate varAttrs for metadata keys
    for (key, dtype, units) in locationKeyList:
        attr = varAttrs[(key, md)]
        attr['units'] = units
        if key == 'sensorChannelNumber':
            attr['_FillValue'] = missing_vals_int
            attr['long_name'] = 'Sensor Channel Number'
        elif key == 'sensorCentralWavenumber':
            attr['_FillValue'] = missing_vals_float
            attr['long_name'] = 'Sensor Central Wavenumber'
        elif key == 'sensorScanPosition':
            attr['_FillValue'] = missing_vals_int
            attr['long_name'] = 'Sensor Scan Position'
        elif key == 'sensorViewAngle':
            attr['_FillValue'] = missing_vals_float
            attr['long_name'] = 'Sensor View Angle'      
        else:
            attr['_FillValue'] = missing_vals_float

    # Channel variable attributes
    varAttrs['Channel'] = {
        'long_name': 'Channel index',
        'units': '',
        '_FillValue': missing_vals_int
    }

    # Radiance variable attributes
    varAttrs[('radiance', iconv.OvalName())] = {
        '_FillValue': missing_vals_float,
        'units': 'W m-2 sr-1 m',
        'long_name': 'LWIR Reconstructed Radiance'
    }

    # Global attributes
    GlobalAttrs = {
        'sensor': WMO_SEN_ID,
        'platform': WMO_SAT_ID,
        'sensorCommonName': 'IRS',
        'platformCommonName': 'MTG-S1',
        'processingLevel': 'Level-1',
        'datetimeReference': datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    }

    # Write IODA
    writer = iconv.IodaWriter(output_file, locationKeyList, DimDict)
    set_metadata_attributes(varAttrs)
    set_obspace_attributes(varAttrs)
    writer.BuildIoda(ioda_data, varDims, varAttrs, GlobalAttrs)

    print(f"[INFO] Wrote IODA file: {output_file} with {nloc_valid} valid locations")


if __name__ == "__main__":
    input_nc = "derived_rad.nc"
    output_ioda = "derived_rad_v29.ioda.nc"

    n_channel_lw, out_channels_lw = read_out_channel_info(input_nc, spectral_group_name='lwir')
    print(f"n_channel_lw = {n_channel_lw}")
    print(f"out_channels_lw = {out_channels_lw}")
    n_channels = n_channel_lw

    data, nloc, datetime_array = read_derived_rad(input_nc)
    write_ioda_2(data, nloc, datetime_array, output_ioda)

