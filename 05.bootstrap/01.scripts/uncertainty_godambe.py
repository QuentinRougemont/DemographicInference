#! ~/py33/bin                                
# -*- coding: utf-8 -*-                      

import sys                                   
import dadi                                  
import modeledemo_new_models_folded          
import numpy                                 
from numpy import array                      
import Godambe                               

#faire lire input                            
try:                                         
    sfs = sys.argv[1]     #sfs file
    model = sys.argv[2]   # #arg2                 
except:   
    print(__doc__)                           


assert model in ["SC2N2mG", "SC2N2m", "SC2N","SC2NG", "SC2m", "SC2mG", "SC", "SCG", 
"IM2N2mG", "IM2N2m", "IM2N","IM2NG", "IM2m", "IM2mG", "IM", "IMG",
"AM2N2mG", "AM2N2m", "AM2N","AM2NG", "AM2m", "AM2mG", "AM", "AMG",
"SI2N","SI2NG", "SI", "SIG"]  , "Model must be 'f1' or 'f2'"

if model == "SC2N2mG":
    funcex = modeledemo_new_models_folded.SC2N2mG
elif model == "IM2N2mG":
    funcex = modeledemo_new_models_folded.IM2N2mG
elif model == "AM2N2mG":
    funcex = modeledemo_new_models_folded.AM2N2mG
elif model == "SC2N2m":
    funcex = modeledemo_new_models_folded.SC2N2m
elif model == "IM2N2m":
    funcex = modeledemo_new_models_folded.IM2N2m
elif model == "AM2N2m":
    funcex = modeledemo_new_models_folded.AM2N2m
elif model == "IM2N":
    funcex = modeledemo_new_models_folded.IM2N
elif model == "SC2N":
    funcex = modeledemo_new_models_folded.SC2N
elif model == "AM2N":
    funcex = modeledemo_new_models_folded.AM2N
elif model == "IM2NG":
    funcex = modeledemo_new_models_folded.IM2NG
elif model == "SC2NG":
    funcex = modeledemo_new_models_folded.SC2NG
elif model == "AM2NG":
    funcex = modeledemo_new_models_folded.AM2NG
elif model == "IM2m":
    funcex = modeledemo_new_models_folded.IM2m
elif model == "SC2m":
    funcex = modeledemo_new_models_folded.SC2m
elif model == "AM2m":
    funcex = modeledemo_new_models_folded.AM2m
elif model == "IM2mG":
    funcex = modeledemo_new_models_folded.IM2mG
elif model == "SC2mG":
    funcex = modeledemo_new_models_folded.SC2mG
elif model == "AM2mG":
    funcex = modeledemo_new_models_folded.AM2mG
elif model == "IM":
    funcex = modeledemo_new_models_folded.IM
elif model == "SC":
    funcex = modeledemo_new_models_folded.SC
elif model == "AM":
    funcex = modeledemo_new_models_folded.AM
elif model == "IMG":
    funcex = modeledemo_new_models_folded.IMG
elif model == "SCG":
    funcex = modeledemo_new_models_folded.SCG
elif model == "AMG":
    funcex = modeledemo_new_models_folded.AM2mG
    
else:
    print("Model must be 'SC2N2mG', 'SC2N2m', 'SC2N','SC2NG', 'SC2m', 'SC2mG', 'SC', 'SCG', 'IM2N2mG', 'IM2N2m', 'IM2N','IM2NG', 'IM2m', 'IM2mG', 'IM', 'IMG','AM2N2mG', 'AM2N2m', 'AM2N','AM2NG', 'AM2m', 'AM2mG', 'AM', 'AMG','SI2N','SI2NG', 'SI', 'SIG' ")
    sys.exit()

fs = dadi.Spectrum.from_file(sfs)          
p0 = numpy.fromfile('param.txt',sep=",")     

sfs = []             
for i in range(100): 
    s=fs.sample()    
    sfs.append(s)    
                    
uncert0 = Godambe.GIM_uncert(    
     func_ex=funcex,   
     grid_pts=numpy.array([80]), 
     all_boot=sfs,               
     p0=p0,                      
     data=fs,                    
     log=False, multinom=True, eps=0.01, return_GIM=False)
                                                                          
with open('model_uncertainty'+str(model), 'w') as f:   
    for p in uncert0:                                                     
         f.write(str(p) +'\n')                                            

