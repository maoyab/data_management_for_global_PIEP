import os
import sys
from pandas import DataFrame
import cPickle as pickle
import numpy as np

def combine_by_resultpath(resultpath, li, fi_name='results_globe_smap', model_params_estimate=['f_max', 'f_w', 's_wilt', 's_star']):
    r = []
    resultpath_i = os.path.join(resultpath, '%s_%s' % (fi_name, li))
    result_files = [os.path.join(resultpath_i, f) for f in os.listdir(resultpath_i) 
                    if f.endswith('pickle') 
                        and (f.startswith('x') == 0) 
                        and (f.startswith('combined')==0)]
    x_result_files = [os.path.join(resultpath_i, f) for f in os.listdir(resultpath_i) 
                    if f.endswith('pickle') 
                        and (f.startswith('x') == 1)]
    for f in result_files:
        with open(f, 'rb') as fp:
            loc_res = pickle.load(fp)
        vv = []
        for v in loc_res.keys():
            vv.append(loc_res[v])
        vv.append(0)
        r.append(vv)
    v_keys = loc_res.keys()
    for f in x_result_files:
        with open(f, 'rb') as fp:
            loc_res = pickle.load(fp)
        vv = []
        for v in v_keys:
            if v in loc_res.keys() and v not in model_params_estimate:
                vv.append(loc_res[v])
            else:
                vv.append(np.nan)
        vv.append(1)
        r.append(vv)
    v_keys.append('x_file')
    if len(r) > 0:
        print len(r)
        r = zip(*r)
        df = {}
        for i, v in enumerate(v_keys):
            df[v] = r[i]
        df = DataFrame.from_dict(df)
        picklename = os.path.join(resultpath, 'combined_%s_%s.pickle' % (fi_name, li))
        with open(picklename, 'wb') as f:
            pickle.dump(df, f)


if __name__ == '__main__':

# save df multiple runs................................................................
    fv = 'results'
    resultpath = '../results'
    li = sys.argv[1]
    combine_by_resultpath(resultpath, li, fi_name='results_globe_smap', model_params_estimate=['f_max', 'f_w', 's_wilt', 's_star'])
