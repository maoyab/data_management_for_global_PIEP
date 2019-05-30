import numpy as np
from pandas import date_range, read_csv, DataFrame
import h5py
import os
import xarray as xr
from datetime import datetime
import calendar

from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib.pyplot as plt


drive = '/Volumes/PIEP_data/'


def get_raster_data(data_filename, transform_filename):
    #from osgeo import gdal
    #dataset = gdal.Open(filename)
    #transform = dataset.GetGeoTransform()
    #x_origin = transform[0]
    #y_origin = transform[3]
    #pixel_width = transform[1]
    #pixel_height = -transform[5]
    #cols = dataset.RasterXSize
    #rows = dataset.RasterYSize
    #band = dataset.GetRasterBand(1)
    #data = band.ReadAsArray(0, 0, cols, rows)
    data = np.load(data_filename)
    transform = np.load(transform_filename)
    x_origin = transform[0]
    y_origin = transform[3]
    pixel_width = transform[1]
    pixel_height = -transform[5]
    return data, x_origin, y_origin, pixel_width, pixel_height


def get_georaster_loc_data(lat, lon, data, x_origin, y_origin, pixel_width, pixel_height):
    col = int((lon - x_origin) / pixel_width)
    row = int((y_origin - lat ) / pixel_height)
    try:
        x = data[row][col]
        if x >= 0:
            return x
        else:
            return np.nan
    except:
        print row, col


def cal_MvG_s(ds, h):
    # h in cm water
    n = ds['n_fit_5cm']
    alpha = ds['alpha_fit_5cm']
    theta_r = ds['s_h'] * ds['n']
    theta_s = ds['n']
    m = 1 - 1 / n
    se = (1 + ( alpha * h) ** n) ** (- m)
    s = (se * (theta_s - theta_r) + theta_r) / theta_s
    return s


def cal_MvG_h(ds, s):
    # h in cm water
    n = ds['n_fit_5cm']
    alpha = ds['alpha_fit_5cm']
    theta_r = ds['s_h'] * ds['n']
    theta_s = ds['n']
    m = 1 - 1 / n
    se = (s * theta_s - theta_r) / (theta_s - theta_r)
    h = ((se ** (-1 / m) - 1) ** (1 / n)) / alpha
    return h


def get_clsm_data(f, name):
    data = f[name][:]
    data = data.astype(float)
    _FillValue = f[name].attrs['_FillValue']
    valid_max = f[name].attrs['valid_max']
    valid_min = f[name].attrs['valid_min']
    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)

    data[invalid] = np.nan
    return np.ma.masked_where(np.isnan(data), data)


def lat_lon_to_36kmEASEGRID2_row_cols(lats, lons):
    # from Brodzik et al., 2012
    # WGS 84 ellipsoid
    a = 6378137.  # equatorial radius [m]
    e = 0.0818191908426  # eccentricity
    c = 36000 # size of grid cell [m]
    phi_1 = 30. * np.pi / 180.  # latitude of true scale
    lamb0 = 0.  # map reference longitude

    x0 = (964 - 1) / 2.
    y0 = (406 - 1) / 2.

    lats[lats < -83.67156] = -83.67156
    lats[lats > 83.67156] = 83.67156
    lons[lons < -179.8133] = -179.8133
    lons[lons > 179.8133] = 179.8133

    phi = lats * np.pi / 180.
    lamb = lons * np.pi / 180.

    k0 = (np.cos(phi_1))\
         / ((1 - (e ** 2) * (np.sin(phi_1) ** 2)) ** 0.5)

    q_lamb = (1 - e ** 2) \
             * (np.sin(phi) / (1 - (e ** 2) * (np.sin(phi) ** 2))\
                - 1 / (2 * e)\
                * np.log((1 - e * np.sin(phi)) / (1 + e * np.sin(phi)))
               )

    x = a * k0 * (lamb - lamb0)
    y = (a * q_lamb) / (2 * k0)

    cols = np.round(x0 + x / c).astype(int)
    rows = np.round(y0 - y / c).astype(int)

    return cols, rows


def merge_att_ds(ds, ds_attributes):
    ds = ds.groupby(['y', 'x']).mean()
    ds = ds.to_xarray()
    ds = ds.assign_coords(col=ds.x)
    ds = ds.assign_coords(row=ds.y)
    return xr.merge([ds_attributes, ds])


def stochastic_rf_char(ds):
    mu = ds.mean(dim='time')['rainfall'].values
    var = ds.var(dim='time')['rainfall'].values
    l = 2 * mu ** 2 / var
    a = mu / l
    return a, l


def make_season_list(n):
    months = range(1, 13)
    seasons = []
    for ii, month in enumerate(months):
        g = [month]
        for ni in range(1, n):
            g.append(months[ii - ni])
        seasons.append(g)
    return seasons


def get_dry_season_rf_attributes(ds_attributes, ds_rf, nyrs, dry_rf_min=0.05):
    ds_monthly_rf = ds_rf.resample('MS', dim='time', how='sum')
    df_annual_rf = ds_rf.sum(dim='time') / nyrs
    ds_monthly_rf['rainfall'] = (['time', 'y', 'x'], ds_monthly_rf['rainfall'].values / df_annual_rf['rainfall'].values)

    (nrows, ncols) = np.shape(ds_attributes['n'].values)
    dry_rf = np.zeros((nrows, ncols))
    dry_rf[:] = np.nan
    dry_season_start = np.zeros((nrows, ncols))
    dry_season_start[:] = np.nan
    n_dry_months = np.zeros((nrows, ncols))
    td = np.zeros((nrows, ncols))
    alpha_rf, lambda_rf = stochastic_rf_char(ds_rf)

    for n in reversed(range(2, 12)):
        for ii, dry_months in enumerate(make_season_list(n)):
            td_i = 0
            dry_rf_i = np.zeros((nrows, ncols))
            for mi in dry_months:
                td_i = td_i + calendar.monthrange(2017, mi)[1]
                m_t = [t for t in ds_monthly_rf['time'].values \
                    if (t.astype('datetime64[M]').astype(int) % 12 + 1) == mi]
                ds_monthly_i = ds_monthly_rf.sel(time=m_t)
                dry_rf_i = dry_rf_i + ds_monthly_i.mean(dim='time')['rainfall'].values

            mask = dry_rf_i < dry_rf_min
            mask[np.isnan(dry_rf) == 0] = False

            wet_t = [t for t in ds_rf['time'].values \
                    if (t.astype('datetime64[M]').astype(int) % 12 + 1) not in dry_months]
            ds_rf_s = ds_rf.sel(time=wet_t)
            alpha_rf_i, lambda_rf_i = stochastic_rf_char(ds_rf_s)

            if ii == 0:
                dry_rf_0 = dry_rf_i
            else:
                mask[dry_rf_i > dry_rf_0] = False

            dry_rf[mask] = dry_rf_i[mask]
            td[mask] = td_i
            dry_season_start[mask] = dry_months[0]
            n_dry_months[mask] = n
            alpha_rf[mask] = alpha_rf_i[mask]
            lambda_rf[mask] = lambda_rf_i[mask]
    mask = np.isnan(ds_attributes['n'])
    dry_rf[mask] = np.nan
    alpha_rf[mask] = np.nan
    lambda_rf[mask] = np.nan
    n_dry_months[mask] = np.nan
    dry_season_start[mask] = np.nan
    td[mask] = np.nan

    return alpha_rf, lambda_rf, n_dry_months, dry_season_start, td, dry_rf


def get_season_et0(ds_attributes, ds_et0):

    (nrows, ncols) = np.shape(ds_attributes['n'].values)
    et0_wet = np.zeros((nrows, ncols))
    et0_wet[:] = np.nan
    et0_dry = np.zeros((nrows, ncols))
    et0_dry[:] = np.nan

    for n in reversed(range(2, 12)):
        for ii, dry_months in enumerate(make_season_list(n)):
            mask = np.logical_and(ds_attributes['n_dry_months'] == n,
                                  ds_attributes['dry_season_start'] == dry_months[0])

            wet_t = [t for t in ds_et0['time'].values \
                    if (t.astype('datetime64[M]').astype(int) % 12 + 1) not in dry_months]
            ds_et0_wet = ds_et0.sel(time=wet_t)
            ds_et0_wet = ds_et0_wet.mean(dim='time')
            et0_wet[mask] = ds_et0_wet['et0'].values[mask]

            dry_t = [t for t in ds_et0['time'].values \
                    if (t.astype('datetime64[M]').astype(int) % 12 + 1) in dry_months]
            ds_et0_dry = ds_et0.sel(time=dry_t)
            ds_et0_dry = ds_et0_dry.mean(dim='time')
            et0_dry[mask] = ds_et0_dry['et0'].values[mask]

    return et0_wet, et0_dry


if __name__ == "__main__":

# resample and save clsm landcover class & soil params nc file.......................................................

    f_clsm = os.path.join(drive, 'Original_Data/n5eil01u.ecs.nsidc.org/SMAP', 'SPL4SMLM.004', 'SMAP_L4_SM_lmc_00000000T000000_Vv4030_001.h5')
    with h5py.File(f_clsm, mode='r') as f:
        clsm_poros = get_clsm_data(f, '/Land-Model-Constants_Data/clsm_poros')
        clsm_wp = get_clsm_data(f, '/Land-Model-Constants_Data/clsm_wp')
        mwrtm_vegcls = get_clsm_data(f, '/Land-Model-Constants_Data/mwrtm_vegcls')

        latitude = f['cell_lat'][:]
        longitude = f['cell_lon'][:]
        n_rows_resampled = np.shape(clsm_poros)[0] // 4
        n_cols_resampled = np.shape(clsm_poros)[1] // 4
        rows_resampled = range(n_rows_resampled)
        cols_resampled = range(n_cols_resampled)

        resampled_latitude = latitude.reshape(n_rows_resampled, 4, n_cols_resampled, 4).mean(-1).mean(1)
        resampled_longitude = longitude.reshape(n_rows_resampled, 4, n_cols_resampled, 4).mean(-1).mean(1)
        resampled_clsm_poros = clsm_poros.reshape(n_rows_resampled, 4, n_cols_resampled, 4).mean(-1).mean(1)
        resampled_clsm_wp = clsm_wp.reshape(n_rows_resampled, 4, n_cols_resampled, 4).mean(-1).mean(1)
        resampled_mwrtm_vegcls = mwrtm_vegcls.reshape(n_rows_resampled, 4, n_cols_resampled, 4).mean(-1).mean(1)

    cols, rows = np.meshgrid(cols_resampled, rows_resampled )
    ds_attributes = xr.Dataset({'clsm_poros': (['y', 'x'], resampled_clsm_poros)})
    ds_attributes['clsm_wp'] = (['y', 'x'], resampled_clsm_wp)
    ds_attributes['vegcls'] = (['y', 'x'], resampled_mwrtm_vegcls)
    ds_attributes['lat'] = (['y', 'x'], resampled_latitude)
    ds_attributes['lon'] = (['y', 'x'], resampled_longitude)
    ds_attributes['col_i'] = (['y', 'x'], cols)
    ds_attributes['row_i'] = (['y', 'x'], rows)

    ds_attributes = ds_attributes.assign_coords(col=ds_attributes.x)
    ds_attributes = ds_attributes.assign_coords(row=ds_attributes.y)
    ds_attributes = ds_attributes.set_index(x='col', inplace=True)
    ds_attributes = ds_attributes.set_index(y='row', inplace=True)

    nc_filename_out = os.path.join(drive, 'Processed_Data/SMAP/attributes', 'clsm_params_EASEGRID_36km.nc')
    ds_attributes.to_netcdf(nc_filename_out)


# add ET0, rainfall AI .......................................................................................................

    print 'adding ET0, rainfall AI'
    nc_attributes = os.path.join(drive, 'Processed_Data/SMAP/attributes', 'clsm_params_EASEGRID_36km.nc')
    ds_attributes = xr.open_dataset(nc_attributes)

    f = os.path.join(drive, 'Processed_Data/SMAP/daily_L4', 'smap_daily_L4_all.nc')
    ds_l4 = xr.open_dataset(f)
    ds_l4 = ds_l4.sel(time=date_range('2015-04-01', '2018-03-31', freq='1D'))
    ds_l4 = ds_l4.mean(dim='time')
    ai = ds_l4['et0'] / ds_l4['rainfall']
    ds_attributes['ai'] = (('y', 'x'), ai)
    ds_attributes['et0_a'] = (('y', 'x'), ds_l4['et0'].values)
    ds_attributes['rainfall'] = (('y', 'x'), ds_l4['rainfall'].values)


# resample Monzka soil water retention parameters and merge ...................................................

    print 'adding soil water retention parameters'
    mvg_file = os.path.join(drive, 'Original_Data/Hydraul_Param_SoilGrids_Schaap_0', 'Hydraul_Param_SoilGrids_Schaap_sl2.nc')
    ds_mvg = xr.open_dataset(mvg_file)

    lats = ds_mvg['latitude'].values
    lons = ds_mvg['longitude'].values
    cols0, rows0 = lat_lon_to_36kmEASEGRID2_row_cols(lats, lons)
    cols, rows = np.meshgrid(cols0, rows0)
    ds_mvg['y'] = (['latitude', 'longitude'], rows)
    ds_mvg['x'] = (['latitude', 'longitude'], cols)
    ds_mvg = ds_mvg.to_dataframe()
    ds_attributes = merge_att_ds(ds_mvg, ds_attributes)


# correct theta_r & theta_s with record min / max calculate s_0.033MPa and s_1.5MPa .......................................................

    print 'adjusting soil parameters'
    f = os.path.join(drive, 'Processed_Data/SMAP/daily_sm_L3', 'smapL3_daily_sm_all.nc')
    ds_sm = xr.open_dataset(f)
    sel_period = date_range('2015-04-01', '2018-03-31', freq='1D')
    ds_sm = ds_sm.sel(time=sel_period)

    ds_attributes['len_s_obs'] = (['y', 'x'], ds_sm.count(dim='time')['soil_moisture'].values)

    c = np.nanmax([ds_attributes['mean_theta_s_5cm'].values, ds_sm.max(dim='time')['soil_moisture'].values], axis=0)
    ds_attributes['mean_theta_s_5cm_corrected'] = (['y', 'x'], c)

    c = np.nanmax([ds_attributes['clsm_poros'].values, ds_sm.max(dim='time')['soil_moisture'].values], axis=0)
    ds_attributes['clsm_poros_corrected'] = (['y', 'x'], c)

    ds_attributes['n'] = (['y', 'x'], ds_attributes['mean_theta_s_5cm_corrected'].values)

    s_h = np.nanmin([ds_attributes['mean_theta_r_5cm'].values / ds_attributes['mean_theta_s_5cm'].values,
                   ds_sm.min(dim='time')['soil_moisture'].values / ds_attributes['n'].values - 0.01], axis=0)

    ds_attributes['s_h'] = (['y', 'x'], s_h)

    ds_attributes['s_1.5MPa'] = (['y', 'x'], cal_MvG_s(ds_attributes, 15295.7))
    ds_attributes['s_0.033MPa'] = (['y', 'x'], cal_MvG_s(ds_attributes, 336.5))
    s_fc = ds_attributes['s_0.033MPa']
    ds_attributes['s_fc'] = (['y', 'x'], s_fc)

    ds_attributes['max_s'] = (['y', 'x'], ds_sm.max(dim='time')['soil_moisture'].values / ds_attributes['n'].values)
    ds_attributes['min_s'] = (['y', 'x'], ds_sm.min(dim='time')['soil_moisture'].values / ds_attributes['n'].values)
    ds_attributes['mean_s'] = (['y', 'x'], ds_sm.mean(dim='time')['soil_moisture'].values / ds_attributes['n'].values)
    ds_attributes['std_s'] = (['y', 'x'], ds_sm.std(dim='time')['soil_moisture'].values / ds_attributes['n'].values)


# resample PFT and merge ...................................................

    pft_file = os.path.join(drive, 'Original_Data/global_maps_of_plant_traits', 'filtered_preds.csv')
    pft_data = read_csv(pft_file)

    lons = pft_data['lon'].values
    lats = pft_data['lat'].values
    resample_locs = []
    steps = [-0.125, 0.125]
    resample_locs = [[lat + step_lat, lon + step_lon]\
                     for step_lon in steps\
                     for step_lat in steps\
                     for lat, lon in zip(lats, lons)]
    latsi, lonsi = zip(*resample_locs)
    pft_df = DataFrame(data={'latitude': latsi, 'longitude': lonsi})
    for i in range(15):
        resample_data = []
        steps = [-0.125, 0.125]
        resample_data = [pfti / 100.\
                         for step_lon in steps\
                         for step_lat in steps\
                         for lat, lon, pfti\
                         in zip(lats, lons, pft_data['PFT%s' % i])]
        pft_df['PFT%s' % i] = resample_data

    lats = pft_df['latitude'].values
    lons = pft_df['longitude'].values
    cols, rows = lat_lon_to_36kmEASEGRID2_row_cols(lats, lons)
    pft_df['y'] = rows
    pft_df['x'] = cols

    ds_attributes = merge_att_ds(pft_df, ds_attributes)
    ds_attributes = ds_attributes.drop(['latitude', 'longitude'])
    for i in range(1, 15):
        mask = np.logical_and(np.isnan(ds_attributes['PFT%s' % i]) == 1, 
                              np.isnan(ds_attributes['clsm_poros']) == 0)
        ds_attributes['PFT%s' % i] = xr.where(mask, 0, ds_attributes['PFT%s' % i])

    mask = np.logical_and(np.isnan(ds_attributes['PFT0']) == 1, 
                              np.isnan(ds_attributes['clsm_poros']) == 0)
    ds_attributes['PFT0'] = xr.where(mask, 1, ds_attributes['PFT0'])


# add RF characteristics annual .......................................

    print 'adding RF characteristics annual'
    f = os.path.join(drive, 'Processed_Data/SMAP/daily_L4/', 'smap_daily_L4_all.nc')
    ds_climate = xr.open_dataset(f)
    ds_rf = ds_climate.drop(['et0', 'specific_humidity'])

    sel_period = date_range('2015-04-01', '2018-03-31', freq='1D')
    ds_rf = ds_rf.sel(time=sel_period)
    a, l = stochastic_rf_char(ds_rf)
    ds_attributes['rf_alpha_a'] = (['y', 'x'], a)
    ds_attributes['rf_lambda_a'] = (['y', 'x'], l)


# add seasonal rainfall characteristics .......................................

    print 'adding RF characteristics seasonal'
    nyrs = 3
    alpha_rf, lambda_rf,\
    n_dry_months, dry_season_start,\
    td, dry_rf = get_dry_season_rf_attributes(ds_attributes, ds_rf, nyrs, dry_rf_min=0.03)

    ds_attributes['rf_alpha_ws'] = (['y', 'x'], alpha_rf)
    ds_attributes['rf_lambda_ws'] = (['y', 'x'], lambda_rf)
    ds_attributes['n_dry_months'] = (['y', 'x'], n_dry_months)
    ds_attributes['dry_season_start'] = (['y', 'x'], dry_season_start)
    ds_attributes['td'] = (['y', 'x'], td)
    ds_attributes['dry_season_rf'] = (['y', 'x'], dry_rf)


# add seasonal et0 characteristics ..............................................

    print 'adding et0 characteristics seasonal'
    ds_et0 = ds_climate.drop(['rainfall', 'specific_humidity'])
    sel_period = date_range('2015-04-01', '2018-03-31', freq='1D')
    ds_et0 = ds_et0.sel(time=sel_period)

    et0_wet, et0_dry = get_season_et0(ds_attributes, ds_et0)

    ds_attributes['et0_wet'] = (['y', 'x'], et0_wet)
    ds_attributes['et0_dry'] = (['y', 'x'], et0_dry)


# save .........................................................................

    print 'saving'
    nc_filename_out = os.path.join(drive, 'Processed_Data/SMAP/attributes', 'sw_model_parameters_for_it.nc')
    ds_attributes.to_netcdf(nc_filename_out)

