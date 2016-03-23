import numpy
from numpy import array

import dadi

import demographic_models

data = dadi.Spectrum.from_file('YRI_CEU.fs')
ns = data.sample_sizes

pts_l = [40,50,60]

func = demographic_models.prior_onegrow_mig
params = array([1.881, 0.0710, 1.845, 0.911, 0.355, 0.111])

func_ex = dadi.Numerics.make_extrap_log_func(func)
model = func_ex(params, ns, pts_l)
ll_model = dadi.Inference.ll_multinom(model, data)
print 'Model log-likelihood:', ll_model
theta = dadi.Inference.optimal_sfs_scaling(model, data)

grid = dadi.Inference.index_exp[1:100:2j,0:3:3j]
popt, fopt, grid, fout, thetas = dadi.Inference.optimize_grid(data, func_ex, pts_l, gridin, full_output=True, fixed_params = (None,params[1],params[2],params[3],None,params[5]), verbose=1)

print 'Optimized parameters', repr(popt)
model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, data)
print 'Optimized log-likelihood:', ll_opt

dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                    pop_ids =('YRI','CEU'))
