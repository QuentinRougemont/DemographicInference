#!/usr/bin/env python3                       
# -*- coding: utf-8 -*-                      

import sys
import dadi
import modeledemo_new_models_folded
#import modeledemo_new_models_folded
import os
import numpy
from numpy import array
#import Godambe                               

#faire lire input                            
try:
    sfs = sys.argv[1]     #sfs file
    model = sys.argv[2]   #name of the model              
    grids = sys.argv[3]   #grid size
except:
    print(__doc__)


assert model in ["SC2N2mG", "SC2N2m", "SC2N","SC2NG", "SC2m", "SC2mG", "SC", "SCG", 
"IM2N2mG", "IM2N2m", "IM2N","IM2NG", "IM2m", "IM2mG", "IM", "IMG",
"AM2N2mG", "AM2N2m", "AM2N","AM2NG", "AM2m", "AM2mG", "AM", "AMG",
"SI2N","SI2NG", "SI", "SIG"]  , "Model must be 'f1' or 'f2'"

if model == "SC2N2mG":
    func = modeledemo_new_models_folded.SC2N2mG
elif model == "IM2N2mG":
    func = modeledemo_new_models_folded.IM2N2mG
elif model == "AM2N2mG":
    func = modeledemo_new_models_folded.AM2N2mG
elif model == "SC2N2m":
    func = modeledemo_new_models_folded.SC2N2m
elif model == "IM2N2m":
    func = modeledemo_new_models_folded.IM2N2m
elif model == "AM2N2m":
    func = modeledemo_new_models_folded.AM2N2m
elif model == "IM2N":
    func = modeledemo_new_models_folded.IM2N
elif model == "SC2N":
    func = modeledemo_new_models_folded.SC2N
elif model == "AM2N":
    func = modeledemo_new_models_folded.AM2N
elif model == "IM2NG":
    func = modeledemo_new_models_folded.IM2NG
elif model == "SC2NG":
    func = modeledemo_new_models_folded.SC2NG
elif model == "AM2NG":
    func = modeledemo_new_models_folded.AM2NG
elif model == "IM2m":
    func = modeledemo_new_models_folded.IM2m
elif model == "SC2m":
    func = modeledemo_new_models_folded.SC2m
elif model == "AM2m":
    func = modeledemo_new_models_folded.AM2m
elif model == "IM2mG":
    func = modeledemo_new_models_folded.IM2mG
elif model == "SC2mG":
    func = modeledemo_new_models_folded.SC2mG
elif model == "AM2mG":
    func = modeledemo_new_models_folded.AM2mG
elif model == "IM":
    func = modeledemo_new_models_folded.IM
elif model == "SC":
    func = modeledemo_new_models_folded.SC
elif model == "AM":
    func = modeledemo_new_models_folded.AM
elif model == "IMG":
    func = modeledemo_new_models_folded.IMG
elif model == "SCG":
    func = modeledemo_new_models_folded.SCG
elif model == "AMG":
    func = modeledemo_new_models_folded.AM2mG
    
else:
    print("Model must be 'SC2N2mG', 'SC2N2m', 'SC2N','SC2NG', 'SC2m', 'SC2mG', 'SC', 'SCG', 'IM2N2mG', 'IM2N2m', 'IM2N','IM2NG', 'IM2m', 'IM2mG', 'IM', 'IMG','AM2N2mG', 'AM2N2m', 'AM2N','AM2NG', 'AM2m', 'AM2mG', 'AM', 'AMG','SI2N','SI2NG', 'SI', 'SIG' ")
    sys.exit()



func_ex = dadi.Numerics.make_extrap_log_func(func)

fs = dadi.Spectrum.from_file(sfs)          
p0 = numpy.fromfile('param.txt',sep=",")     
print(p0)

grid_pts = numpy.fromfile('grids.txt', sep = ",")
#print(grid_pts)
#grid_pts = numpy.array(grid_pts)
print(fs.S())
print(grid_pts)

if not os.path.exists('sample_data/'):
      os.makedirs('sample_data/')

nboots = 100
sfs = []
for i in range(nboots):
    s=fs.sample()
    sfs.append(s)
    s.to_file("sample_data/boot{0}.fs".format(str(i)))


fi = open(model + 'demographic_confidence_intervals.txt','w')
fi.write('Optimized parameters: {0}\n\n'.format(p0))

pts_l = [140,160,180]
print("------------ running uncertainty-----------")
# we want to try a few different step sizes (eps) to see if
# uncertainties very wildly with changes to step size.
for eps in [0.01, 0.001, 0.0001]:
     uncert = dadi.Godambe.GIM_uncert(
     func_ex=func_ex,
     grid_pts=pts_l,
     all_boot=sfs,
     p0=p0,
     data=fs,
     log=False, multinom=True, eps=eps, return_GIM=False)
     fi.write('Estimated 95% uncerts (with step size '+str(eps)+'): {0}\n'.format(1.96*uncert[:-1]))
     fi.write('Lower bounds of 95% confidence interval : {0}\n'.format(p0-1.96*uncert[:-1]))
     fi.write('Upper bounds of 95% confidence interval : {0}\n\n'.format(p0+1.96*uncert[:-1]))
fi.close()

