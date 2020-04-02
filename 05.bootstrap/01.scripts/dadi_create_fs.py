#!/usr/bin/env python2


import dadi
import os
import pylab
import sys
from os import path

#faire lire input
try: 
    input = sys.argv[1]    #basename of the ms output file
    P = float(sys.argv[2]) #P parameter value 
    Q = float(sys.argv[3]) #Q parameter value
    id = int(sys.argv[4])  #id of the seq
    model = sys.argv[5] 
except:
    print(__doc__)

input1=input +"neutral"
input2=input +"islands"
input3=input +"normal_recomb"
input4=input +"linkedsel"

#Normal region
fsN = dadi.Spectrum.from_ms_file(input1)                                               
#Island region
fsI = dadi.Spectrum.from_ms_file(input2)                                               
#Normal region (no hfr)
fsNr= dadi.Spectrum.from_ms_file(input3)                                               
#region under Hrf
fsLr= dadi.Spectrum.from_ms_file(input4)                                               

#combine the spectrum:
fs = Q*fsLr + (1-Q)*fsNr + P*fsN + (1-P)*fsI
fs=fs.fold()


fs.to_file("02.boot_fs_"+ str(model) +"/data_boot." + str(id) + ".fs")
