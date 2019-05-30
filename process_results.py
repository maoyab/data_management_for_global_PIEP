import os
import sys
from pandas import HDFStore, DataFrame, merge, read_pickle, date_range
import cPickle as pickle
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime
import statsmodels.api as sm
sys.path.insert(0, '../PIEP')
from sswm import SM_A




def cal_MvG_s(h, alpha, n, theta_s, theta_r):
    m = 1 - 1/ n
    se = (1 + (alpha * h) ** n) ** (-m)
    s = (se * (theta_s - theta_r) + theta_r) / theta_s
    return s



def cal_MvG_h(s, alpha, n, theta_s, theta_r):
    m = 1 - 1/ n
    se = (s * theta_s - theta_r) / (theta_s - theta_r)
    h = ((se ** (-1 / m) - 1) ** (1 / n)) / alpha
    return h


def df_to_ds(df0):
    df0['loc_i'] = df0.index
    print df0['loc_i']

    v = 'loc_i'
    df = df0.set_index(['row', 'col', 'time'])
    df = df[v].unstack()
    df = df.to_panel()

    var_list = df0.keys()
    var_list = [vi for vi in var_list if vi not in ['row', 'col', 'time', 'sobs', 's_obs']]

    ds = xr.Dataset({v: (['time', 'y', 'x'], df.values)},
                      coords={'col': ('x', df.minor_axis), 'row': ('y', df.major_axis), 'time': df.items})
    for vi in var_list:
        print vi
        df = df0.set_index(['row', 'col', 'time'])
        df = df[vi].unstack()
        df = df.to_panel()
        ds[vi] = (['time', 'y', 'x'], df.values)
    ds = ds.mean(dim='time')
    return ds


def apply_theta(theta, s_obs):
    def __nse(obs, mod):
        mo = np.mean(obs)
        a = np.sum([(mi - oi) ** 2 for mi, oi in zip(mod, obs)])
        b = np.sum([(oi - mo) ** 2 for oi in obs])
        return 1 - a / b

    def __evapotranspiration(s):
        #laio et al., 2001
        if s > smpdf.s_star:
            return smpdf.e_max
        elif s > smpdf.s_wilt:
            return (smpdf.e_max - smpdf.e_w) * (s - smpdf.s_wilt) / (smpdf.s_star - smpdf.s_wilt) + smpdf.e_w
        elif s > smpdf.s_h:
            return smpdf.e_w * (s - smpdf.s_h) / (smpdf.s_wilt - smpdf.s_h) + smpdf.e_w
        else:
            return 0

    def __plant_stress(s, q = 2):
        # Porporato et al., 2001
        if s <= smpdf.s_wilt:
            return 1
        elif (s < smpdf.s_wilt and s <= smpdf.s_star):
            return ((smpdf.s_star - s) / (smpdf.s_star - smpdf.s_wilt)) ** q
        else:
            return 0

    smpdf = SM_A(theta)
    p_fitted_norm = smpdf.get_p0()
    cdf = np.cumsum(p_fitted_norm)
    cdf_m_n = cdf / np.max(cdf)

    s_mod = [(np.abs(cdf_m_n - qi / 365.)).argmin() / np.float(len(p_fitted_norm) - 1) for qi in range(1, 365)]
    obs_l = [np.percentile(s_obs, qi / 365. * 100) for qi in range(1, 365)]
    NSE = __nse(obs_l, s_mod)

    mean_plant_stress = np.nanmean([__plant_stress(s, q=2) for s in s_obs])
    et_mod = np.nanmean([__evapotranspiration(s, ) for s in s_obs])
    norm_wu = et_mod / (smpdf.rf_alpha * smpdf.rf_lambda)
    swnwu = norm_wu * (1 - mean_plant_stress)

    return NSE, mean_plant_stress, norm_wu, swnwu


def test_theta(df_ri):
    if np.isnan(df_ri['s_wilt']) == 0:
        theta_ref = [df_ri['s_h'], df_ri['s_1.5MPa'], df_ri['s_0.033MPa'], df_ri['s_fc'],
                 1, 0.0001, df_ri['et0'],
                 df_ri['rf_alpha'], df_ri['rf_lambda'], df_ri['delta'],
                 df_ri['b'], df_ri['Ks'], df_ri['n'], df_ri['Zr'],
                 df_ri['et0_dry'], df_ri['t_d']]

        theta_est = [df_ri['s_h'], df_ri['s_wilt'], df_ri['s_star'], df_ri['s_fc'],
                 df_ri['f_max'], df_ri['f_w'], df_ri['et0'],
                 df_ri['rf_alpha'], df_ri['rf_lambda'], df_ri['delta'],
                 df_ri['b'], df_ri['Ks'], df_ri['n'], df_ri['Zr'],
                 df_ri['et0_dry'], df_ri['t_d']]
        s_obs = ds_sm.isel(y=df_ri['row_i'], x=df_ri['col_i'])['soil_saturation'].values
        s_obs = [s for s in s_obs if np.isnan(s) == 0]

        nse_rc, stress_rc, norm_wu_rc, swnwu_rc = apply_theta(theta_ref, s_obs)
        nse, stress, norm_wu, swnwu = apply_theta(theta_est, s_obs)
        return [nse_rc, stress_rc, norm_wu_rc, swnwu_rc, nse, stress, norm_wu, swnwu]
    else:
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

if __name__ == '__main__':


    data_path = '../../data/Processed_Data'
    results_path = '../../results'

# combine results by iterable....................................................................................................

    li = 0
    picklename = os.path.join(results_path, 'results_by_iterable', 'combined_results_globe_smap_annual_f_%s.pickle' % li)
    with open(picklename, 'rb') as f:
        df0 = pickle.load(f)
        print picklename
        print len(df0.index)
    for li in range(1, 11):
        picklename = os.path.join(results_path, 'results_by_iterable', 'combined_results_globe_smap_annual_f_%s.pickle' % li)
        with open(picklename, 'rb') as f:
            df = pickle.load(f)
            print picklename
            print len(df.index)
        df0 = df0.append(df)
    print len(df0.index)
    picklename = os.path.join(results_path, 'intermediary_files', 'combined_results_smap_global_annual_f.pickle')
    with open(picklename, 'wb') as f:
        pickle.dump(df0, f)
    


# merge attributes, calculate MvG thresholds make .nc...........................................................
    
    picklename = os.path.join(results_path, 'intermediary_files', 'combined_results_smap_global_annual_f.pickle')
    print picklename
    with open(picklename, 'rb') as f:
        df_r = pickle.load(f)

    f = os.path.join(data_path, 'sw_model_parameters_for_it.nc')
    ds_p = xr.open_dataset(f)

    f = os.path.join(data_path, 'smapL3_daily_sm_all.nc')
    ds_sm = xr.open_dataset(f)
    sel_period = date_range('2015-04-01', '2018-03-31', freq='1D')
    ds_sm = ds_sm.sel(time=sel_period)
    ds_sm['soil_saturation'] = ds_sm['soil_moisture'] / ds_p['n']

    df_p = ds_p.to_dataframe()
    sel_vars = [kk for kk in df_p.keys() if kk not in df_r.keys()]
    df_p = df_p[sel_vars]

    df_r['x'] = df_r['col']
    df_r['y'] = df_r['row']
    df_r = merge(df_r, df_p, on=['y', 'x'], how='outer')

    df_r['psi_0'] = cal_MvG_h(df_r['s_wilt'], df_r['alpha_fit_5cm'], df_r['n_fit_5cm'],
                                         df_r['n'], df_r['s_h'] * df_r['n']) * 9.8067 * 10 ** (-5)

    df_r['psi_1'] = cal_MvG_h(df_r['s_star'], df_r['alpha_fit_5cm'], df_r['n_fit_5cm'],
                                         df_r['n'], df_r['s_h'] * df_r['n']) * 9.8067 * 10 ** (-5)

    df_r['time'] = '2015-2018'

    test_0 = [test_theta(df_ri) for index, df_ri in df_r.iterrows()]
    nse_rc, stress_rc, norm_wu_rc, swnwu_rc, nse, stress, norm_wu, swnwu = zip(*test_0)

    df_r['NSE_pdf_rc'] = nse_rc
    df_r['stress_index_rc'] = stress_rc
    df_r['swnwu_rc'] = swnwu_rc
    df_r['norm_wu_rc'] = norm_wu_rc

    df_r['NSE_pdf'] = nse
    df_r['stress_index'] = stress
    df_r['swnwu'] = swnwu
    df_r['norm_wu'] = norm_wu

    df_r = df_r.groupby(['y', 'x']).mean()
    ds = df_r.to_xarray()

    mask = ds['num_pl'] < 3
    piep_var = ['s_star', 's_wilt', 'f_max', 'f_w', 'psi_0', 'psi_1',
                's_wilt_grd', 's_star_grd', 'f_max_grd', 'f_w_grd', 'efficiency',
                'NSE_pdf', 'stress_index', 'swnwu', 'norm_wu']
    for v in piep_var:
        ds[v] = xr.where(mask, np.nan, ds[v])

    nc_filename = os.path.join(results_path, 'intermediary_files', 'sswbinv_results_smap_global_annual_f.pickle')
    print ds
    ds.to_netcdf(nc_filename)
    

# save selected attributes to nc data for publication...............................................................................

    f = os.path.join(results_path, 'intermediary_files', 'sswbinv_results_smap_global_annual_f.pickle')
    ds = xr.open_dataset(f)

    ds['pct Bare'] = ds['PFT0']
    ds['pct C4'] = ds['PFT14']
    ds['pct C3'] = ds['PFT12'] + ds['PFT13']
    ds['pct Tree'] = ds['PFT1'] + ds['PFT2'] + ds['PFT3'] + ds['PFT4'] + ds['PFT5'] + ds['PFT6'] + ds['PFT7'] + ds['PFT8']
    ds['pct Shrub'] = ds['PFT9'] + ds['PFT10'] + ds['PFT11']
    ds['pct Crop'] = 1 - ds['pct C4'] - ds['pct C3'] - ds['pct Shrub'] - ds['pct Tree'] - ds['pct Bare']
    ds['alpha_MvG'] = ds['alpha_fit_5cm']
    ds['n_MvG'] = ds['n_fit_5cm']
    ds['Z'] = ds['Zr']
    ds['E_p'] = ds['et0']
    ds['aridity_index'] = ds['ai']
    ds['latitude'] = ds['lat']
    ds['longitude'] = ds['lon']
    selected_vars = [['latitude', 'degrees', 'latitude of grid centroid'],
                    ['longitude', 'degrees', 'longitude of grid centroid'],
                    ['len_s_obs', 'unitless', 'number of soil moisture observations in L3 SMAP timeseries (April 2015 to March 2018) used in analysis'],
                    ['aridity_index', 'unitless', 'aridity index, ratio of average potential evapotranspiration to rainfall'],
                    ['max_s', 'unitless', 'maximum soil saturation value in L3 SMAP timeseries (April 2015 to March 2018)'],
                    ['min_s', 'unitless', 'minimum soil saturation value in L3 SMAP timeseries (April 2015 to March 2018)'],
                    ['mean_s', 'unitless', 'mean soil saturation value in L3 SMAP timeseries (April 2015 to March 2018)'],
                    ['std_s', 'unitless', 'standard deviation of observed soil saturation values in L3 SMAP timeseries (April 2015 to March 2018)'],
                    ['alpha_MvG', 'unitless', 'empirical parameter used in the Mulalem-van Genuchten equation'],
                    ['n_MvG', 'unitless', 'empirica; parameter used in the Mulalem-van Genuchten equation'],
                    ['Z', 'mm', 'soil depth'],
                    ['n', 'unitless', 'soil porosity'],
                    ['b', 'unitless', 'empirical parameter used in the Clapp and Hornberger soil water retention curve equation'],
                    ['Ks', 'mm/day', 'saturated soil hydraulic conductivity'],
                    ['s_fc', 'unitless', 'soil saturation at field capacity'],
                    ['s_h', 'unitless', 'soil saturation at the hygroscopic point'],
                    ['s_1.5MPa', 'unitless', 'soil saturation at 1.5 MPa soil water potential'],
                    ['s_0.033MPa', 'unitless', 'soil saturation at 0.033 MPa soil water potential'],
                    ['rf_alpha', 'mm/day', 'average daily rainfall depth'],
                    ['rf_lambda', 'unitless', 'average daily rainfall frequency'],
                    ['E_p', 'mm/day', 'average daily potential evapotranspiration'],
                    ['s_star', 'unitless', 'soil saturation at the point of incipient stomatal closure'],
                    ['s_wilt', 'unitless', 'soil saturation at the wilting point'],
                    ['f_max', 'unitless', 'ratio of maximum soil water uptake to potential evapotranspiration'],
                    ['f_w', 'unitless', 'ratio of soil water uptake at the wilting point to potential evapotranspiration'],
                    ['psi_0', 'MPa', 'soil water potential at the point of no soil water uptake'],
                    ['psi_1', 'MPa', 'soil water potential at the point of downregulation of soil water uptake'],
                    ['s_wilt_grd', 'unitless', 'Gelman-Rubin diagnostic for s_wilt'],
                    ['s_star_grd', 'unitless', 'Gelman-Rubin diagnostic for s_star'],
                    ['f_max_grd', 'unitless', 'Gelman-Rubin diagnostic for f_max'],
                    ['f_w_grd', 'unitless', 'Gelman-Rubin diagnostic for f_w'],
                    ['efficiency', 'unitless', 'efficiency of Metropolis-Hastings Markov chain Monte Carlo algorithm'],
                    ['NSE_pdf', 'unitless', 'quantile level Nash Sutcliffe efficiency between theoretical and empirical soil saturation probability distribution using inferred critical ecohydrological thresholds'],
                    ['NSE_pdf_rc', 'unitless', 'quantile level Nash Sutcliffe efficiency between theoretical and empirical soil saturation probability distribution using constant reference critical ecohydrological thresholds'],
                    ['stress_index', 'unitless', 'soil moisture stress index estimated using inferred critical ecohydrological thresholds'],
                    ['stress_index_rc', 'unitless', 'soil moisture stress index estimated using constant reference critical ecohydrological thresholds'],
                    ['norm_wu', 'unitless', 'soil water use normalized by rainfall estimated using inferred critical ecohydrological thresholds'],
                    ['norm_wu_rc', 'unitless', 'soil water use normalized by rainfall estimated using constant reference critical ecohydrological thresholds'],
                    ['swnwu', 'unitless', 'stress weighted normalized water use estimated using inferred critical ecohydrological thresholds'],
                    ['swnwu_rc', 'unitless', 'stress weighted normalized water use estimated using constant reference critical ecohydrological thresholds'], 
                    ['pct Tree', 'unitless', 'abundance of tree plant functional type'],
                    ['pct Shrub', 'unitless', 'abundance of shrub plant functional type'],
                    ['pct C3', 'unitless', 'abundance of C3 grass plant functional type'],
                    ['pct C4', 'unitless', 'abundance of C4 grass plant functional type'],
                    ['pct Bare', 'unitless', 'percentage  bare ground'],
                    ['pct Crop', 'unitless', 'abundance of crop plant functional type'],
                    ['vegcls', 'unitless', 'IGBP land cover class'],
                    ]

    drop_vars = [ii for ii in ds.data_vars.keys() if ii not in zip(*selected_vars)[0]]
    ds = ds.drop(drop_vars)

    ds.attrs = {'title':'Global Dataset of Ecohydrological Parameters',
                'author': 'Maoya Bassiouni', 'contact': 'bassioum@oregonstate.edu',
                'date': '2019',
                'grid': '36 km Equal-Area Scalable Earth Grid, Version 2.0 (EASE-Grid 2.0) in a global cylindrical projection',
                'description': 'input parameters and results associated with analysis in "Global Variation in Critical Soil Water Potentials" (Bassiouni et al., in preparation)',
                'methods': 'for methods and references to original datasets used see "Global Variation in Critical Soil Water Potentials" (Bassiouni et al., in preparation)'}
    for v, unit, descr in selected_vars:
        ds[v].attrs['unit'] = unit
        ds[v].attrs['description'] = descr
    f_out = os.path.join(results_path, 'piep_smap_results.nc')
    ds.to_netcdf(f_out)
    ds_re = xr.open_dataset(f_out)
    print ds_re
