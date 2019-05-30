import numpy as np
from pandas import date_range
import os
import xarray as xr
from cPickle import dump


def save_iterable(ds, ds_p, cols_iii, rows_iii, ix, it_name, it_path, seasonal=False):
    iterables = []
    for nit, (row, col) in enumerate(zip(rows_iii, cols_iii)):
        ds_pi = ds_p.isel(y=row, x=col)
        ds_i = ds.isel(y=row, x=col)
        s_obs = ds_i['soil_saturation'].values.flatten()
        s_obs = [np.float(s) for s in s_obs if np.isnan(s) == 0]

        iti = {'lat': np.float(ds_pi['lat'].values),
                 'lon': np.float(ds_pi['lon'].values),
                 'loc_i': cols_iii.index(col) + ix,
                 'row': row,
                 'col': col,
                 'date_0': ds_i['time'].values[0],
                 'date_f': ds_i['time'].values[1],
                 's_obs': s_obs,
                 'len_s_obs': len(s_obs),
                 'min_s': np.min(s_obs),
                 'mean_s': np.mean(s_obs),
                 'max_s': np.max(s_obs),
                 'std_s': np.std(s_obs),
                 'delta': 0,
                 'Zr': 50.,
                 'n': np.float(ds_pi['n'].values),
                 's_fc': np.float(ds_pi['s_fc'].values),
                 's_star': np.float(ds_pi['s_0.033MPa'].values),  # placeholder for init
                 's_wilt': np.float(ds_pi['s_1.5MPa'].values),  # placeholder for init
                 's_h': np.float(ds_pi['s_h'].values),
                 'b': 1 / (np.float(ds_pi['n_fit_5cm'].values) - 1),
                 'Ks': np.float(ds_pi['mean_Ks_5cm'].values) * 10.,
                 'f_w': 0.05,  # placeholder until random init
                 'f_max': 1,  # placeholder until random init
                 'f_max_dry': 1,  # placeholder until random init
                 }
        if seasonal is True:
            iti['t_d'] = np.float(ds_pi['td'].values)
            s = ds_pi['dry_season_start'].values
            n = np.int(ds_pi['n_dry_months'].values)
            iti['dry months'] = [s + mi for mi in range(n)]
            iti['rf_lambda'] = np.float(ds_pi['rf_lambda_ws'].values)
            iti['rf_alpha'] = np.float(ds_pi['rf_alpha_ws'].values)
            iti['rf_lambda'] = np.float(ds_pi['rf_lambda_ws'].values)
            iti['et0'] = np.float(ds_pi['et0_wet'].values)
            iti['et0_dry'] = np.float(ds_pi['et0_dry'].values)

        else:
            iti['t_d'] = 0
            iti['dry months'] = []
            iti['rf_lambda'] = np.float(ds_pi['rf_lambda_a'].values)
            iti['rf_alpha'] = np.float(ds_pi['rf_alpha_a'].values)
            iti['et0'] = np.float(ds_pi['et0_a'].values)
            iti['et0_dry'] = np.float(ds_pi['et0_a'].values)

        iterables.append([nit, iti])

    picklename = '%s_%s.pickle' % (it_name, iii)
    picklename = os.path.join(it_path, picklename)
    with open(picklename, 'wb') as f:
        dump(iterables, f, -1)


if __name__ == "__main__":


# make iterables...........................................................
    it_path = '../../data/iterables/'
    sel_period = date_range('2015-04-01', '2018-03-31', freq='1D')

    f = os.path.join('../../data/Processed_Data', 'sw_model_parameters_for_it.nc')
    ds_p = xr.open_dataset(f)

    f = os.path.join('../../data/Processed_Data', 'smapL3_daily_sm_all.nc')
    ds = xr.open_dataset(f)
    ds = ds.sel(time=sel_period)

    ds['soil_saturation'] = ds['soil_moisture'] / ds_p['n']
    ds_p = ds_p.where(ds_p['len_s_obs'] >= 365)
    ds_p = ds_p.where(np.isnan(ds_p['s_fc']) == 0)
    ### for seasonal select loc where dry season length >0 only
    #ds_p = ds_p.where(ds_p['td'] > 0)
    rows = ds_p['row_i'].values.flatten()
    rows = rows[np.isnan(rows) == 0]
    rows = [np.int(r) for r in rows]

    cols = ds_p['col_i'].values.flatten()
    cols = cols[np.isnan(cols) == 0]
    cols = [np.int(c) for c in cols]

    n_locs = len(rows)
    #52379
    print n_locs

    n_per_it = 4800
    # can make iterables with less locatiosn if seasonal
    n_iterables = np.int(np.floor(n_locs / np.float(n_per_it)))

    for iii in range(n_iterables):
        cols_iii = cols[iii * n_per_it: (1 + iii) * n_per_it]
        rows_iii = rows[iii * n_per_it: (1 + iii) * n_per_it]
        print iii * n_per_it, (1 + iii) * n_per_it
        save_iterable(ds, ds_p, cols_iii, rows_iii, iii * n_per_it,
                      'globe_smap_annual', it_path, seasonal=False)

    iii = n_iterables 
    cols_iii = cols[iii * n_per_it:]
    rows_iii = rows[iii * n_per_it:]
    print iii * n_per_it, n_locs
    save_iterable(ds, ds_p, cols_iii, rows_iii, iii * n_per_it,
                  'globe_smap_annual', it_path, seasonal=False)

