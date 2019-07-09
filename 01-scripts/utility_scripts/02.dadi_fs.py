#!/usr/bin/python2.7

import sys
import dadi
import pylab

try:
        proj = int(sys.argv[1])  # size of fs to project must be equal or smaller than the number of individuals
except:
        print(__doc__)

dd = dadi.Misc.make_data_dict ('spectrum_NP.data2')
fs = dadi.Spectrum.from_data_dict ( dd , pop_ids =[ 'pop1', 'pop2'] , projections =[proj,proj] , polarized = False) 
s = fs.S()
fst = fs.Fst()

#python3
#print(f'number of s = {s}')
#print(f'differentiation - Fst = {fst}') 

dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.1)

fs.to_file("pop1_pop2_NP_" + str(proj) + ".fs")

fs . mask [0 ,1] = True
fs . mask [1 ,0] = True
s1 = fs.S()
fst1 = fs.Fst()
fs.to_file("pop1_pop2_NP_masked" + str(proj) + ".fs")

print 'number of site = {0:.5f}.'.format(s)
print 'differentiation (Fst) = {0:.5f}.'.format(fst)

print 'number of site masking singleton = {0:.5f}.'.format(s1)
print 'differentiation (Fst) without singleton = {0:.5f}.'.format(fst1)
