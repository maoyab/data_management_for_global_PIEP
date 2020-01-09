import numpy as np
from pandas import date_range
import h5py
import os
import xarray as xr
from datetime import datetime
import calendar


drive = '/Volumes/PIEP_data/'


def get_filename_from_date_L4(date):
    hdf_folder = os.path.join(drive, 'Original_Data/n5eil01u.ecs.nsidc.org/SMAP', 'SPL4SMGP.004')
    f = os.path.join(hdf_folder, '%d.%02d.%02d' % (date.year, date.month, date.day),
        'SMAP_L4_SM_gph_%d%02d%02dT%02d3000_Vv4030_001.h5' \
        % (date.year, date.month, date.day, date.hour))
    return f


def resample_L4_day(dates, day):
    names = ['/Geophysical_Data/precipitation_total_surface_flux',
             '/Geophysical_Data/temp_lowatmmodlay',
             '/Geophysical_Data/net_downward_shortwave_flux',
             '/Geophysical_Data/net_downward_longwave_flux',
             '/Geophysical_Data/heat_flux_ground',
             '/Geophysical_Data/specific_humidity_lowatmmodlay']
    all_data = []
    for name in names:
        resampled_data = []
        for date in dates:
            f_i = get_filename_from_date_L4(date)
            with h5py.File(f_i, mode='r') as f:
                data = f[name][:]
                _FillValue = f[name].attrs['_FillValue']
                valid_max = f[name].attrs['valid_max']
                valid_min = f[name].attrs['valid_min']
                invalid = np.logical_or(data > valid_max,
                                        data < valid_min)
                invalid = np.logical_or(invalid, data == _FillValue)
                data[invalid] = np.nan
                n_rows_rs = np.shape(data)[0] // 4
                n_cols_rs = np.shape(data)[1] // 4
                resampled_data.append(data.reshape(n_rows_rs, 4, n_cols_rs, 4).mean(-1).mean(1))
        all_data.append(resampled_data)
    rainfall, ta, rs, rl, g, sh = all_data
    pet = [et0(tai, rsi, rli, gi) for tai, rsi, rli, gi in zip(ta, rs, rl, g)]

    time = [day for i in dates]
    rows_rs = range(n_rows_rs)
    cols_rs = range(n_cols_rs)
    ds = xr.Dataset({'rainfall': (['time', 'y', 'x'], rainfall)},
                      coords={'row': rows_rs, 'col': cols_rs, 'time': time})
    ds['et0'] = (['time', 'y', 'x'], pet)
    ds['specific_humidity'] = (['time', 'y', 'x'], sh)
    ds = ds.resample(time='D').mean()

    ds['rainfall'] = ds['rainfall'] * 3600 * 24 / 997 * 1000.  # convert kg m-2 s-1 to mm/day
    ds['et0'] = ds['et0'] / 28.94  # convert W/m2 to mm /day
    return ds


def et0(Tair, Rs, Rl, G):
    # Priestley and Taylor, 1972
    gamma = 81.8 * 0.665 * 10 ** (-3)  # kPa / kelvin
    delta = vaporPressureSlope(Tair)
    a_pt = 1.26
    return a_pt * delta / (delta + gamma) * (Rs + Rl - G)


def airVaporPressureSat(Ta):
    # eq 2.14 p.28 Brutsaert, 2005
    a0 = 6984.505294
    a1 = -188.903931
    a2 = 2.133357675
    a3 = -1.288580973 * 10.0 ** (-2)
    a4 = 4.393587233 * 10.0 ** (-5)
    a5 = -8.023923082 * 10.0 ** (-8)
    a6 = 6.136820929 * 10.0 ** (-11)
    eStar = (a0 + Ta * (a1 + Ta * (a2 + Ta * (a3 + Ta * (a4 + Ta * (a5 + Ta * a6)))))) * 0.1  # kPa
    return eStar


def vaporPressureSlope(Ta):
    # eq. 2.2 p. 42 Brutsaert, 1982
    T = 273.16
    eStar = airVaporPressureSat(Ta)
    tra = 1.0 - T / Ta
    delta = (T * eStar / Ta ** 2) * (13.3185 - 3.952 * tra - 1.9335 * tra ** 2 - 0.5196 * tra ** 3)  # kPa/C
    return delta


if __name__ == "__main__":


# resample L4 9km:3h to 36km: 1day ............................................................

    days = date_range('2015-03-31', '2019-03-31', freq='1D')
    for day in days:
        print day
        datetime_day = date_range('%d-%02d-%02d 01:30:00' % (day.year, day.month, day.day), \
                                  '%d-%02d-%02d 22:30:00' % (day.year, day.month, day.day), freq='3H')
        ds = resample_L4_day(datetime_day, day)
        nc_filename_out = os.path.join(drive, 'Processed_Data/SMAP/daily_L4', 'smap_daily_L4_%d%02d%02d.nc' % (day.year, day.month, day.day))
        ds.to_netcdf(nc_filename_out)

# combine L4 daily ts by month nc files ............................................................
    for year in [2015, 2016, 2017, 2018, 2019]:
        for month in range(1, 13):
            print year, month
            n_days = calendar.monthrange(year, month)[1]
            days = date_range('%d-%02d-01' % (year, month), '%d-%02d-%02d' % (year, month, n_days), freq='1D')
            ds = []
            for day in days:
                f = os.path.join(drive, 'Processed_Data/SMAP/daily_L4', 'smap_daily_L4_%d%02d%02d.nc' \
                    % (day.year, day.month, day.day))
                try:
                    dsi = xr.open_dataset(f)
                    ds.append(dsi)
                    dsi.close()
                except:
                    print f
            ds = xr.merge(ds)
            f = os.path.join(drive, 'Processed_Data/SMAP/daily_L4', 'smap_daily_L4_%d%02d_m.nc' % (year, month))
            ds.to_netcdf(f)

# combine month nc files to year nc files............................................................
    for year in [2015, 2016, 2017, 2018, 2019]:
        ds = []
        for month in range(1, 13):
            print year, month
            f = os.path.join(drive, 'Processed_Data/SMAP/daily_L4', 'smap_daily_L4_%d%02d_m.nc' % (year, month))
            dsi = xr.open_dataset(f)
            ds.append(dsi)
            dsi.close()
        ds = xr.merge(ds)
        f = os.path.join(drive, 'Processed_Data/SMAP/daily_L4', 'smap_daily_L4_%d_y.nc' % (year))
        ds.to_netcdf(f)


# combine yearly nc files to all nc file ............................................................
    ds = []
    for year in [2015, 2016, 2017, 2018, 2019]:
        print year
        f = os.path.join(drive, 'Processed_Data/SMAP/daily_L4', 'smap_daily_L4_%d_y.nc' % (year))
        dsi = xr.open_dataset(f)
        ds.append(dsi)
        dsi.close()

    ds = xr.merge(ds)
    f = os.path.join(drive, 'Processed_Data/SMAP/daily_L4', 'smap_daily_L4_all.nc')
    ds.to_netcdf(f)
