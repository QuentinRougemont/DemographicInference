#! ~/py33/bin                                                                                                                                                                                                      
# -*- coding: utf-8 -*-                                                                                                                                                                                            
                                                                                                                                                                                                                   
                                                                                                                                                                                                                   
import sys                                                                                                                                                                                                         
import dadi                                                                                                                                                                                                        
import modeledemo_new_models_folded                                                                                                                                                                                
import numpy                                                                                                                                                                                                       
from numpy import array                                                                                                                                                                                            
import Godambe                                                                                                                                                                                                     
                                                                                                                                                                                                                   
# IM2N2mG                                                                                                                                                                                                          
#read external argument                                                                                                                                                                                            
#faire lire input                                                                                                                                                                                                  
try:                                                                                                                                                                                                               
    input = sys.argv[1]    #basename of the ms output file                                                                                                                                                         
except:                                                                                                                                                                                                            
    print(__doc__)                                                                                                                                                                                                 
                                                                                                                                                                                                                   
fs = dadi.Spectrum.from_file(input)                                                                                                                                                                                
p0 = numpy.fromfile('param.txt',sep=",")                                                                                                                                                                           

#param.txt = file containing the parameter values of the best model obtained from dadi

print(model)                                                                                                                                                                                                       
                                                                                                                                                                                                                   
sfs = []                                                                                                                                                                                                           
for i in range(100):                                                                                                                                                                                               
        s=fs.sample()                                                                                                                                                                                              
        sfs.append(s)                                                                                                                                                                                              
                                                                                                                                                                                                                   
uncert0 = Godambe.GIM_uncert(                                                                                                                                                                                      
     func_ex=modeledemo_new_models_folded.IM2N2mG,                                                                                                                                                                 
     grid_pts=numpy.array([80]),                                                                                                                                                                                   
     all_boot=sfs,                                                                                                                                                                                                 
     p0=p0,                                                                                                                                                                                                        
     data=fs,                                                                                                                                                                                                      
     log=False, multinom=True, eps=0.01, return_GIM=False)                                                                                                                                                         
                                                                                                                                                                                                                   
                                                                                                                                                                                                                   
with open('calif1_thompson.SC2N2m_3_masked_2020_2_7_632148', 'w') as f:                                                                                                                                            
    for p in uncert0:                                                                                                                                                                                              
         f.write(str(p) +'\n')     
