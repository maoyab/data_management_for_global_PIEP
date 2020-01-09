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


def test_theta(df_ri, ds_sm):
    if np.isnan(df_ri['s_wilt']) == 0:
        s_obs = ds_sm.isel(y=df_ri['row_i'], x=df_ri['col_i'])['soil_saturation'].values
        s_obs = [s for s in s_obs if np.isnan(s) == 0]
        if len(s_obs) > 100:
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

            nse_rc, stress_rc, norm_wu_rc, swnwu_rc = apply_theta(theta_ref, s_obs)
            nse, stress, norm_wu, swnwu = apply_theta(theta_est, s_obs)
            return [nse_rc, stress_rc, norm_wu_rc, swnwu_rc, nse, stress, norm_wu, swnwu]
        else:
            return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    else:
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]


if __name__ == '__main__':


    data_path = '../../data/Processed_Data'
    results_path = '../../results'

# combine results by iterable...................................................................................
    '''
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
    '''

# merge attributes, calculate MvG thresholds make .nc...........................................................
    '''
    picklename = os.path.join(results_path, 'intermediary_files', 'combined_results_smap_global_annual_f.pickle')
    print picklename
    with open(picklename, 'rb') as f:
        df_r = pickle.load(f)

    f = os.path.join(data_path, 'sw_model_parameters_for_it_revision.nc')
    ds_p = xr.open_dataset(f)

    f = os.path.join(data_path, 'smapL3_daily_sm_all.nc')
    ds_sm = xr.open_dataset(f)

    sel_period = date_range('2018-04-01', '2019-03-31', freq='1D')
    ds_sm_v = ds_sm.sel(time=sel_period)
    ds_sm_v['soil_saturation'] = ds_sm_v['soil_moisture'] / ds_p['n']

    sel_period = date_range('2015-04-01', '2018-03-31', freq='1D')
    ds_sm_e = ds_sm.sel(time=sel_period)
    ds_sm_e['soil_saturation'] = ds_sm_e['soil_moisture'] / ds_p['n']

    sel_period = date_range('2015-04-01', '2019-03-31', freq='1D')
    ds_sm_o = ds_sm.sel(time=sel_period)
    ds_sm_o['soil_saturation'] = ds_sm_o['soil_moisture'] / ds_p['n']

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

    df_r['time'] = '2015-2019'

    test_0 = [test_theta(df_ri, ds_sm_e) for index, df_ri in df_r.iterrows()]
    nse_rc_e, stress_rc_e, norm_wu_rc_e, swnwu_rc_e, nse_e, stress_e, norm_wu_e, swnwu_e = zip(*test_0)

    df_r['NSE_pdf_rc_e'] = nse_rc_e
    df_r['stress_index_rc_e'] = stress_rc_e
    df_r['norm_wu_rc_e'] = norm_wu_rc_e
    df_r['swnwu_rc_e'] = swnwu_rc_e

    df_r['NSE_pdf_e'] = nse_e
    df_r['stress_index_e'] = stress_e
    df_r['norm_wu_e'] = norm_wu_e
    df_r['swnwu_e'] = swnwu_e

    df_r['rf_alpha'] = df_r['rf_alpha_a_o']
    df_r['rf_lambda'] = df_r['rf_lambda_a_o']
    df_r['et0'] = df_r['et0_a_o']
    test_0 = [test_theta(df_ri, ds_sm_o) for index, df_ri in df_r.iterrows()]
    nse_rc_o, stress_rc_o, norm_wu_rc_o, swnwu_rc_o, nse_o, stress_o, norm_wu_o, swnwu_o = zip(*test_0)

    df_r['NSE_pdf_rc_o'] = nse_rc_o
    df_r['stress_index_rc_o'] = stress_rc_o
    df_r['norm_wu_rc_o'] = norm_wu_rc_o
    df_r['swnwu_rc_o'] = swnwu_rc_o

    df_r['NSE_pdf_o'] = nse_o
    df_r['stress_index_o'] = stress_o
    df_r['norm_wu_o'] = norm_wu_o
    df_r['swnwu_o'] = swnwu_o

    df_r['rf_alpha'] = df_r['rf_alpha_a_v']
    df_r['rf_lambda'] = df_r['rf_lambda_a_v']
    df_r['et0'] = df_r['et0_a_v']
    test_0 = [test_theta(df_ri, ds_sm_v) for index, df_ri in df_r.iterrows()]
    nse_rc_v, stress_rc_v, norm_wu_rc_v, swnwu_rc_v, nse_v, stress_v, norm_wu_v, swnwu_v = zip(*test_0)

    df_r['NSE_pdf_rc_v'] = nse_rc_v
    df_r['stress_index_rc_v'] = stress_rc_v
    df_r['norm_wu_rc_v'] = norm_wu_rc_v
    df_r['swnwu_rc_v'] = swnwu_rc_v

    df_r['NSE_pdf_v'] = nse_v
    df_r['stress_index_v'] = stress_v
    df_r['norm_wu_v'] = norm_wu_v
    df_r['swnwu_v'] = swnwu_v

    df_r = df_r.groupby(['y', 'x']).mean()
    ds = df_r.to_xarray()

    nc_filename = os.path.join(results_path, 'intermediary_files', 'sswbinv_results_smap_global_annual_f_revision.pickle')
    print ds
    ds.to_netcdf(nc_filename)
    '''

# save selected attributes to nc data for publication...........................................................

    f = os.path.join(results_path, 'intermediary_files', 'sswbinv_results_smap_global_annual_f_revision.pickle')
    ds = xr.open_dataset(f)

    mask = ds['num_pl'] < 3
    piep_var = ['s_star', 's_wilt', 'f_max', 'f_w', 'psi_0', 'psi_1',
                's_wilt_grd', 's_star_grd', 'f_max_grd', 'f_w_grd',
                's_wilt_std', 's_star_std', 'f_max_std', 'f_w_std', 'efficiency',
                'NSE_pdf_e', 'stress_index_e', 'swnwu_e', 'norm_wu_e',
                'NSE_pdf_v', 'stress_index_v', 'swnwu_v', 'norm_wu_v',
                'NSE_pdf_o', 'stress_index_o', 'swnwu_o', 'norm_wu_o']
    for v in piep_var:
        ds[v] = xr.where(mask, np.nan, ds[v])

    ds['alpha_MvG'] = ds['alpha_fit_5cm']
    ds['n_MvG'] = ds['n_fit_5cm']
    ds['Z'] = ds['Zr']

    ds['E_p_v'] = ds['et0_a_v']
    ds['E_p_e'] = ds['et0_a']

    ds['rf_alpha_v'] = ds['rf_alpha_a_v']
    ds['rf_lambda_v'] = ds['rf_lambda_a_v']

    ds['rf_alpha_e'] = ds['rf_alpha_a']
    ds['rf_lambda_e'] = ds['rf_lambda_a']

    ds['aridity_index'] = ds['ai_o']
    ds['latitude'] = ds['lat']
    ds['longitude'] = ds['lon']
    ds['stress_index'] =  ds['stress_index_o']
    ds['stress_index_rc'] =  ds['stress_index_rc_o']
    ds['swnwu'] = ds['swnwu_o']
    ds['swnwu_rc'] = ds['swnwu_rc_o']
    ds['norm_wu'] = ds['norm_wu_o']
    ds['norm_wu_rc'] = ds['norm_wu_rc_o']

    ds['max_s'] = ds['max_s_o']
    ds['min_s'] = ds['min_s_o']
    ds['mean_s'] = ds['mean_s_o']
    ds['std_s'] = ds['std_s_o']

    ds['max_s_v'] = ds['max_s']
    ds['min_s_v'] = ds['min_s']
    ds['mean_s_v'] = ds['mean_s']
    ds['std_s_v'] = ds['std_s']

    selected_vars = [['latitude', 'degrees', 'latitude of grid centroid'],
                    ['longitude', 'degrees', 'longitude of grid centroid'],
                    ['len_s_obs', 'unitless', 'number of soil moisture observations in used in estimation (L3 SMAP 04/2015 to 03/2018)'],
                    
                    ['aridity_index', 'unitless', 'ratio of average potential evapotranspiration to rainfall (L4 SMAP 04/2015 to 03/2019)'],
                    ['vegcls', 'unitless', 'IGBP land cover class'],
                    
                    ['alpha_MvG', 'unitless', 'empirical parameter used in the Mulalem-van Genuchten equation'],
                    ['n_MvG', 'unitless', 'empirical parameter used in the Mulalem-van Genuchten equation'],
                    ['Z', 'mm', 'soil depth'],
                    ['n', 'unitless', 'soil porosity'],
                    ['b', 'unitless', 'empirical parameter used in the Clapp and Hornberger soil water retention curve equation'],
                    ['Ks', 'mm/day', 'saturated soil hydraulic conductivity'],
                    ['s_fc', 'unitless', 'soil saturation at field capacity'],
                    ['s_h', 'unitless', 'soil saturation at the hygroscopic point'],
                    ['s_1.5MPa', 'unitless', 'soil saturation at -1.5 MPa soil water potential'],
                    ['s_0.033MPa', 'unitless', 'soil saturation at -0.033 MPa soil water potential'],
                    
                    ['rf_alpha_e', 'mm/day', 'average daily rainfall depth of estimation period (L4 SMAP 04/2015 to 03/2018)'],
                    ['rf_lambda_e', 'unitless', 'average daily rainfall frequency of estimation period (L4 SMAP 04/2015 to 03/2018)'],
                    ['E_p_e', 'mm/day', 'average daily potential evapotranspiration of estimation period (L4 SMAP 04/2015 to 03/2018)'],
                    
                    ['rf_alpha_v', 'mm/day', 'average daily rainfall depth of validation period (L4 SMAP 04/2018 to 03/2019)'],
                    ['rf_lambda_v', 'unitless', 'average daily rainfall frequency of validation period (L4 SMAP 04/2018 to 03/2019)'],
                    ['E_p_v', 'mm/day', 'average daily potential evapotranspiration of validation period (L4 SMAP 04/2018 to 03/2019)'],

                    ['s_star', 'unitless', 'soil saturation at the point of incipient stomatal closure, mean of posterior estimates'],
                    ['s_wilt', 'unitless', 'soil saturation at the wilting point, mean of posterior estimates'],
                    ['f_max', 'unitless', 'ratio of maximum rate of surface soil water losses from ET to E_p, mean of posterior estimates'],
                    ['f_w', 'unitless', 'ratio of surface soil water losses at the wilting point to E_p, mean of posterior estimates'],
                    ['s_wilt_std', 'unitless', 'standard deviation of posterior estimates of s_wilt'],
                    ['s_star_std', 'unitless', 'standard deviation of posterior estimates of s_star'],
                    ['f_max_std', 'unitless', 'standard deviation of posterior estimates of f_max'],
                    ['f_w_std', 'unitless', 'standard deviation of posterior estimates of f_w'],
                    ['s_wilt_grd', 'unitless', 'Gelman-Rubin diagnostic for s_wilt'],
                    ['s_star_grd', 'unitless', 'Gelman-Rubin diagnostic for s_star'],
                    ['f_max_grd', 'unitless', 'Gelman-Rubin diagnostic for f_max'],
                    ['f_w_grd', 'unitless', 'Gelman-Rubin diagnostic for f_w'],
                    ['f_w_grd', 'unitless', 'Gelman-Rubin diagnostic for f_w'],
                    ['efficiency', 'unitless', 'efficiency of Metropolis-Hastings Markov chain Monte Carlo algorithm'],
                    
                    ['psi_0', 'MPa', 'soil water potential when plant water uptake is zero'],
                    ['psi_1', 'MPa', 'soil water potential when plant water uptake is downregulated from its maximum rate'],
                    
                    ['NSE_pdf_e', 'unitless', 'estimation period quantile level Nash Sutcliffe efficiency between best-fit theoretical and empirical soil saturation probability distribution using inferred ecohydrological parameters'],
                    ['NSE_pdf_rc_e', 'unitless', 'estimation period quantile level Nash Sutcliffe efficiency between referent theoretical and empirical soil saturation probability distribution using constant reference ecohydrological parameters'],
                     
                    ['NSE_pdf_v', 'unitless', 'validation period quantile level Nash Sutcliffe efficiency between best-fit theoretical and empirical soil saturation probability distribution using inferred ecohydrological parameters'],
                    ['NSE_pdf_rc_v', 'unitless', 'validation period quantile level Nash Sutcliffe efficiency between reference theoretical and empirical soil saturation probability distribution using constant reference ecohydrological parameters'],
                   
                    ['stress_index', 'unitless', 'average soil moisture stress index estimated using inferred ecohydrological parameters (04/2015 to 03/2019)'],
                    ['stress_index_rc', 'unitless', 'average soil moisture stress index estimated using reference constants (04/2015 to 03/2019)'],
                    ['norm_wu', 'unitless', 'average water uptake normalized by rainfall estimated using inferred ecohydrological parameters (04/2015 to 03/2019)'],
                    ['norm_wu_rc', 'unitless', 'average water uptake normalized by rainfall estimated using reference constants (04/2015 to 03/2019)'],
                    ['swnwu', 'unitless', 'average stress weighted normalized water uptake estimated using inferred ecohydrological parameters (04/2015 to 03/2019)'],
                    ['swnwu_rc', 'unitless', 'average stress weighted normalized water uptake estimated using reference constants (04/2015 to 03/2019)'], 
                    ['max_s', 'unitless', 'maximum soil saturation of estimation data (L3 SMAP 04/2015 to 03/2018)'],
                    ['min_s', 'unitless', 'minimum soil saturation of estimation data(L3 SMAP 04/2015 to 03/2018)'],
                    ['mean_s', 'unitless', 'mean soil saturation of estimation data(L3 SMAP 04/2015 to 03/2018)'],
                    ['std_s', 'unitless', 'standard deviation of estimation data(L3 SMAP 04/2015 to 03/2018)'],
                    
                    ['max_s_v', 'unitless', 'maximum soil saturation value in L3 SMAP timeseries (April 2018 to March 2019)'],
                    ['min_s_v', 'unitless', 'minimum soil saturation value in L3 SMAP timeseries (April 2018 to March 2019)'],
                    ['mean_s_v', 'unitless', 'mean soil saturation value in L3 SMAP timeseries (April 2018 to March 2019)'],
                    ['std_s_v', 'unitless', 'standard deviation of observed soil saturation values in L3 SMAP timeseries (April 2018 to March 2019)'],
                    ]

    drop_vars = [ii for ii in ds.data_vars.keys() if ii not in zip(*selected_vars)[0]]
    ds = ds.drop(drop_vars)

    ds.attrs = {'title':'Global Dataset of Ecohydrological Parameters',
                'author': 'Maoya Bassiouni', 'contact': 'maoya.bassiouni@slu.se',
                'date': '2020',
                'grid': '36 km Equal-Area Scalable Earth Grid, Version 2.0 (EASE-Grid 2.0) in a global cylindrical projection',
                'description': 'input parameters and results associated with analysis in "Plant Water Uptake Thresholds Inferred from Satellite Soil Moisture" (Bassiouni et al., in preparation)',
                'methods': 'for methods and references to original datasets used see "Plant Water Uptake Thresholds Inferred from Satellite Soil Moisture" (Bassiouni et al., in preparation)'}
    for v, unit, descr in selected_vars:
        ds[v].attrs['unit'] = unit
        ds[v].attrs['description'] = descr
    f_out = os.path.join(results_path, 'piep_smap_results_revision.nc')
    ds.to_netcdf(f_out)
    ds_re = xr.open_dataset(f_out)
    print ds_re
    
