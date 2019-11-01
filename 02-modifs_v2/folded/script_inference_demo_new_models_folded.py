#! ~/py33/bin
# -*- coding: utf-8 -*-

# Library importation
import os
import sys
import getopt
import time
import matplotlib
matplotlib.use('Agg')

from scipy import stats
import pylab
import numpy
from numpy import array
import dadi
import modeledemo_new_models_folded
#from matplotlib import pyplot as plt

#Help function
def usage():
    """ Function for help """
    print("# This script allow you to test different demographic models on your genomic data\n"+
        "# and will give you which one is the best fitted.\n\n"+
        "# This is an exemple of the most complete command line :\n"+
        "# -o pathoutput -y population1 -x population2 -p 10,20,30 -f pathfsfile -m SI,SI2N,SIG,SI2NG,EM,IM,IMG,IM2N,IM2N2m,IM2NG,IM2mG,IM2N2mG, AM,AMG,AM2N,AM2N2m,AM2NG,AM2mG,AM2N2mG,PAM,SC,SC2N,SCG,SC2N2m,SC2NG,SC2mG,SC2N2mG,PSC,EM2M,IM2m,AM2m,SC2m,EM2M2P,IM2M2P,AM2M2P,PAM2M2P,SC2M2P,PSC2M2P -l -a -h -v\n\n"+
        "# This is an exemple of the shortest command line:\n"+
        "# -f pathfsfile\n\n"+
        "# -h --help : Display the help you are looking at.\n"+
        "# -v --verbose : Print steps while the code is running\n"+
        "# -y --population1 : Take the name of the first population in the sfs (y-axis)\n"+
        "# -x --population2 : Take the name of the second population in the sfs (x-axis)\n"+
        "# -o --outputname : Take the path of output file.\n"+
        "# -f --fs_file_name : Take the path of the fs file from thr parent directory.\n"+
        "# -p --grid_points : Take 3 numbers separated by a coma, for the size of grids for extrapolation.\n"+
        "# -m --model_list : Take until 18 name of model (SI,SI2N,SIG,SI2NG,EM,IM,IMG,IM2N,IM2N2m, IM2NG,IM2mG,IM2N2mG, AM,AMG,AM2N,AM2N2m,AM2NG,AM2mG,AM2N2mG,PAM,SC,SC2N,SCG,SC2N2m,SC2NG,SC2mG,SC2N2mG,PSC,EM2M,IM2m,AM2m,SC2m,EM2M2P,IM2M2P,AM2M2P,PAM2M2P,SC2M2P,PSC2M2P) separated by a coma.\n"+
        "# For more information on models see docstrings in the module modeledemo.\n"+
        "# -z : mask the singletons.\n"+
        "# -l : record the final parameters in the output file.\n\n\n"
        "########################## Enjoy ###########################")
    return()

#Argument function
def takearg(argv):
    """ Function which record arguments from the command line."""
    # default values
    masked = False # freq 0,1 and 1,0 masked if masked = 1
    pts_l = None  # Grids sizes for extrapolation
    outputname = "mis_fs_2d_optlog"
    model_list = ["SI", "SI2N", "SIG", "SI2NG", "IM", "IMG", "IM2N", "IM2N2m","IM2NG", "IM2mG", "IM2N2mG", "AM", "AMG", "AM2N","AM2N2m", "AM2NG", "AM2mG", "AM2N2mG", "SC", "SC2N", "SCG" , "SC2N2m", "SC2NG", "SC2mG", "SC2N2mG", "IM2m", "AM2m", "SC2m", "PIM2m", "PIM2N","PSC2m","PSC2N"]
    verbose = False
    logparam = False
    nompop1 = "Pop1"
    nompop2 = "Pop2"

    checkfile = False #initilization. if True fs file needed exists, if False it doesn't

    if len(argv) < 2:
        print("You should give, at least, the name of the fs file !")
        sys.exit(1)
    try:
        opts, args = getopt.getopt(argv[1:], "hvo:y:x:azf:p:m:l", ["help", "verbose", "outputname=", "population1=", "population2=", "masked", "fs_file_name=", "grid_points=", "model_list=", "log"])
    except getopt.GetoptError as err:
        # Affiche l'aide et quitte le programme
        print(err) # Va afficher l'erreur en anglais
        usage() # Fonction à écrire rappelant la syntaxe de la commande
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()                     
            sys.exit()
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-o", "--outputname"):
            outputname = arg
        elif opt in ("-y", "--population1"):
            nompop1 = arg
        elif opt in ("-x", "--population2"):
            nompop2 = arg
        elif opt in ("-z", "--masked"):
            masked = True
        elif opt in ("-f", "--fs_file_name"):
            fs_file_name = arg
            checkfile = True
        elif opt in ("-p", "--grid_points"):
            pts_l = arg.split(",")
        elif opt in ("-m", "--model_list"):
            model_list = arg.split(",")
        elif opt in ("-l", "--log"):
            logparam = True
        else:
            print("Option {} inconnue".format(opt))
            sys.exit(2)
    if not checkfile:
        print("You should give, at least, the name of the fs file !")
        sys.exit(1)
    return(masked, pts_l, outputname, nompop1, nompop2, fs_file_name, model_list, verbose, logparam)

#Inference function
def callmodel(func, data, output_file, modeldemo, ll_opt_dic, nbparam_dic,
          nompop1="Pop1", nompop2="Pop2", params=None, fixed_params=None, lower_bound=None, upper_bound=None,
          pts_l=None, ns=None,outputname=None, verbose=False, maxiter=20, 
          Tini=50, Tfin=0, learn_rate=0.005, schedule= "cauchy"):

    # Make the extrapolating version of our demographic model function.
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    # Calculate the model AFS.
    model = func_ex(params, ns, pts_l)
    # Likelihood of the data given the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, data)
    print 'Model log-likelihood:', ll_model
    # The optimal value of theta (4*No*u) given the model.
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    print 'theta:', theta
    
    # Do the optimization. By default we assume that theta is a free parameter,
    # since it's trivial to find given the other parameters. If you want to fix
    # theta, add a multinom=False to the call.
    # (This is commented out by default, since it takes several minutes.)
    # The maxiter argument restricts how long the optimizer will run. For production
    # runs, you may want to set this value higher, to encourage better convergence.
    # Tini = initial temperature of the chain.
        # Learn rate = decreasing rate in the probability of accepting worse solutions as it explores the solution space. 
    if optimizationstate == "anneal_hot" :
        # Perturb our parameter array before optimization. This does so by taking each
        # parameter a up to a factor of two up or down.
        p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=lower_bound, upper_bound=upper_bound)

        popt = dadi.Inference.optimize_anneal(p0, data, func_ex, pts_l, 
                              lower_bound=lower_bound,
                              upper_bound=upper_bound,
                              verbose=verbose,
                              maxiter=maxiter, Tini=Tini, Tfin=Tfin, 
                              learn_rate=learn_rate, schedule=schedule)
    elif optimizationstate == "anneal_cold" :
        popt = dadi.Inference.optimize_anneal(params, data, func_ex, pts_l, 
                              lower_bound=lower_bound,
                              upper_bound=upper_bound,
                              verbose=verbose,
                              maxiter=maxiter/2, Tini=Tini/2, Tfin=Tfin, 
                              learn_rate=learn_rate*2, schedule=schedule)
 
    else :
        popt = dadi.Inference.optimize_log(params, data, func_ex, pts_l, 
                           lower_bound=lower_bound,
                           upper_bound=upper_bound,
                           verbose=verbose,
                           maxiter=maxiter/2)
    
    # Computation of statistics
    model = func_ex(popt, ns, pts_l)
    ll_opt = dadi.Inference.ll_multinom(model, data)
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    AIC = 2*len(params)-2*ll_opt
    ll_opt_dic[modeldemo] = ll_opt
    nbparam_dic[modeldemo] = len(params)

    # Print results
    print 'Optimized parameters', repr(popt)
    print 'Optimized log-likelihood:', ll_opt
    print 'theta:', theta
    
    # Write results
    line = ("\n" + str(modeldemo) + "\n" + "Model log-likelihood: " + repr(ll_model) + "\n" "Optimization : " + repr(optimizationstate) + "\n"  "Optimized parameters: " + repr(popt) + "\n" + "Optimized log-likelihood: " + repr(ll_opt) + "\n" + "theta: " + repr(theta) + "\n" + "AIC: " + repr(AIC) + "\n")
    output_file.write(line)

    # Plot a comparison of the resulting fs with the data.
    if optimizationstate == "BFGS" :
        import pylab
        pylab.figure()
        dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=0.1, resid_range=3,
                                pop_ids =(nompop1,nompop2),
                                saveplot=True, nomplot=(outputname + "_" + modeldemo), showplot=False)
    done=True
    return(done, ll_opt_dic, nbparam_dic, popt)

##############################
##############################

# Load parameters
masked, pts_l, outputname, nompop1, nompop2, fs_file_name, model_list, verbose, logparam = takearg(sys.argv)

if pts_l != None:
    for i in range(len(pts_l)):
        pts_l[i] = int(pts_l[i])

# Load the data
data = dadi.Spectrum.from_file(fs_file_name)
ns = data.sample_sizes

# Creation of outputname and setting default params if they are not in the args
datastate = "not_masked"
opt_list = ["anneal_hot", "anneal_cold", "BFGS"]
if pts_l == None:
    pts_l = [ns[0]+10,ns[0]+20,ns[0]+30]
if masked:
    data.mask[1,0] = True
    data.mask[0,1] = True
    outputname = outputname + "_masked"
    datastate = "masked"

outputname = outputname + "_" + repr(time.localtime()[0]) + "_" + repr(time.localtime()[1]) + "_" + repr(time.localtime()[2]) + "_" + repr(time.localtime()[3]) + repr(time.localtime()[4]) + repr(time.localtime()[5])

# Create output dir and file
os.mkdir(outputname)
output_file = open((outputname + "/" + outputname + ".txt"), "w")

# Save the parameters
if logparam :
    line = ("Model(s) : " + repr(model_list) + "\n" + "Data state : " + repr(datastate) + "\n" + "Grid points : " + repr(pts_l) + "\n\n\n")
    output_file.write(line)
    
# Create dic for ll to make lrt
ll_opt_dic = {}
nbparam_dic = {}

# ML inference for each model
for namemodel in model_list:
    print namemodel
    time.sleep(1.0)

    if namemodel == "SI":

        # Custom Simple Isolation model: nu1, nu2, Ts, O
        func = modeledemo_new_models_folded.SI

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2])
            else:
                params = (popt[0], popt[1], popt[2])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 15,]
            lower_bound = [0.01, 0.01, 0,]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "SI2N":

        # Custom Simple Isolation model with different recombination along the genome: nu1, nu2, Ts, nr, bf
        func = modeledemo_new_models_folded.SI2N

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 0.5, 1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4])
            else:
                params = (popt[0], popt[1], popt[2], popt[3], popt[4])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 15, 1, 0.999 ]
            lower_bound = [0.01, 0.01, 0, 0, 0.001]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "SIG": 

        # Custom Strict Isolation model: nu1, nu2, b1, b2, Ts, O
        func = modeledemo_new_models_folded.SIG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4])
            else:
                params = (popt[0], popt[1], popt[2], popt[3], popt[4])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [80, 80, 100, 100, 15]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "SI2NG":

        # Custom Simple Isolation model with different recombination along the genome: nu1, nu2, b1, b2, hrf, Ts, Q
        func = modeledemo_new_models_folded.SI2NG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 1, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])
            else:
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 100, 100, 1, 15, 0.5]
            lower_bound = [0.01, 0.01, 0.01, 00.1, 0.01, 0.01, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "IM":

        # Custom Isolation with Migration model: nu1, nu2, m12, m21, Ts 
        func = modeledemo_new_models_folded.IM

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4])
        
            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 60, 60, 15]
            lower_bound = [0.01, 0.01, 0, 0, 0]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "IMG":

        # Custom Isolation with Migration model: nu1, nu2, b1, b2, m12, m21, Ts
        func = modeledemo_new_models_folded.IMG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 1, 1, 1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])
            else:
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 100, 100, 60, 60, 15]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "IM2N":

        # Custom Ancient Migration with 2 Migration rate model: nu1, nu2, hrf, m12, m21, Ts, Q 
        func = modeledemo_new_models_folded.IM2N

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 0.8, 5, 5, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 1, 60, 60, 15, 0.5]
            lower_bound = [0.01, 0.01, 0.1, 0, 0, 0, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "IM2N2m":

        # Custom Ancient Migration with 2 Migration rate model and 2 population size: nu1, nu2, hrf, m12, m21, me12, me21, Ts, P, Q 
        func = modeledemo_new_models_folded.IM2N2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 0.8, 5, 5, 0.5, 0.5, 1, 0.5, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 1, 60, 60, 20, 20, 15, 0.95, 0.5]
            lower_bound = [0.01, 0.01, 0.1, 0, 0, 0, 0, 0, 0.01, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "IM2N2mG":

        # Custom Ancient Migration with 2 Migration rate model and 2 population size: nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, P, Q 
        func = modeledemo_new_models_folded.IM2N2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 1, 1, 0.8, 5, 5, 0.5, 0.5, 1, 0.5, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 150,150,1, 60, 60, 20, 20, 15, 0.95, 0.5]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0.1, 0, 0, 0, 0, 0, 0.01, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "IM2NG":

        # Custom Ancient Migration with 2 Migration rate model: nu1, nu2, b1, b2, hrf, m12, m21, Ts, Q 
        func = modeledemo_new_models_folded.IM2NG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 1, 1, 0.8, 5, 5, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100,150,150, 1, 60, 60, 15, 0.5]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0.1, 0, 0, 0, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "IM2mG":

        # Custom Isolation with 2 Migration rate model: nu1, nu2, b1, b2, m12, m21, me12, me21, Ts, P, O 
        func = modeledemo_new_models_folded.IM2mG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 1, 1, 5, 5, 0.5, 0.5, 1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])
        
            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 200, 200, 40, 40, 30, 30, 10, 0.95]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if (namemodel == "AM") or (namemodel == "PAM"):

        # Custom Ancient Migration Model: nu1, nu2, m12, m21, Ts, Tam
        if namemodel == "AM":
            func = modeledemo_new_models_folded.AM
        else :
            func = modeledemo_new_models_folded.PAM

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
        
            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 60, 60, 15, 2]
            lower_bound = [0.01, 0.01, 0, 0, 0, 0]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))
        
    if namemodel == "AM2N":

        # Custom Ancient Migration with 2 Migration rate model: nu1, nu2, hrf, m12, m21, Tam, Ts, Q
        func = modeledemo_new_models_folded.AM2N

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 0.8, 5, 5, 0.1, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100,100, 1, 60, 60, 12, 15, 0.5]
            lower_bound = [0.01, 0.01, 0.1, 0, 0, 0, 0, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))
        
        
    if (namemodel == "AMG"):

        # Custom Ancient Migration Model: nu1, nu2, b1, b2, m12, m21, Tam, Ts
        func = modeledemo_new_models_folded.AMG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 1, 1, 0.1, 1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])
            else:
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 100, 100, 30, 30, 10, 2]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "AM2N2m":

        # Custom Ancient Migration with 2 Migration rate model: nu1, nu2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q, O
        func = modeledemo_new_models_folded.AM2N2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 0.8, 1, 1, 5, 5, 0.1, 1, 0.5, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 1, 60, 60, 20, 20, 10, 15, 0.95, 0.5]
            lower_bound = [0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "AM2NG":

        # Custom Ancient Migration with 2 Migration rate model: nu1, nu2, b1, b2, hrf, m12, m21, Tam, Ts, Q, O
        func = modeledemo_new_models_folded.AM2NG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 1, 1, 0.8, 5, 5, 0.1, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [90, 90, 60, 60, 1, 50, 50, 10, 10, 0.5]
            lower_bound = [0.01, 0.01, 0.1, 0, 0, 0, 0, 0, 0, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "AM2mG":

        # Custom Ancient Migration with 2 Migration rate model: nu1, nu2, b1, b2, m12, m21, me12, me21, Tam, Ts, P, O
        func = modeledemo_new_models_folded.AM2mG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 1, 1, 5, 5, 0.5, 0.5, 0.1, 1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [90, 90, 40, 40, 30, 30, 20, 20, 2, 10, 0.95]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "AM2N2mG":

        # Custom Ancient Migration with 2 Migration rate model:  nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q
        func = modeledemo_new_models_folded.AM2N2mG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 1, 1, 0.8, 1, 1, 5, 5, 0.1, 1, 0.5, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 50, 50, 1, 60, 60, 20, 20, 7, 15, 0.95, 0.5]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))



    if (namemodel == "SC") or (namemodel == "PSC"):

        # Custom Simple Secondary Contact Model: nu1, nu2, m12, m21, Ts, Tsc
        if namemodel == "SC":
            func = modeledemo_new_models_folded.SC
        else :
            func = modeledemo_new_models_folded.PSC
        
        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])

        
            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [200, 200, 60, 60, 15, 5]
            lower_bound = [0.01, 0.01, 0, 0, 0, 0]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))
        
        
    if namemodel == "SCG":

        # Custom Secondary contact with exponential growth model: nu1, nu2, b1, b2, m12, m21, Ts, Tsc
        func = modeledemo_new_models_folded.SCG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 1, 1, 1, 1, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 150, 150, 60, 60, 15, 5]
            lower_bound = [0.01, 0.01, 0, 0, 0, 0, 0, 0]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "SC2N":

        # Custom Simple Secondary Contact Model: nu1, nu2, hrf, m12, m21, Ts, Tsc, Q
        func = modeledemo_new_models_folded.SC2N

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 0.8, 1, 1, 1, 0.1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 1, 60, 60, 15, 5, 0.5]
            lower_bound = [0.01, 0.01, 0.1, 0, 0, 0, 0, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "SC2N2m":

        # Custom Simple Secondary Contact Model: nu1, nu2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q
        func = modeledemo_new_models_folded.SC2N2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 0.8, 1, 1, 1, 1, 10, 3, 0.5, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])

        
            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 1, 60, 60, 30, 30, 15, 5, 0.95, 0.5]
            lower_bound = [0.01, 0.01, 0.1, 0, 0, 0, 0, 0, 0, 0.5, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "SC2NG":

        # Custom Simple Secondary Contact Model: nu1, nu2, b1, b2, hrf, m12, m21, Ts, Tsc, Q
        func = modeledemo_new_models_folded.SC2NG
        
        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 0.8, 1, 1, 1, 0.1, 0.95)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])

        
            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 150, 150, 1, 60, 60, 15, 8, 0.5]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0.1, 0, 0, 0, 0, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "SC2mG":

        # Custom Secondary contact with 2 Migration rate model: nu1, nu2, b1, b2, m12, m21, me12, me21, Ts, Tsc, P
        func = modeledemo_new_models_folded.SC2mG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 1, 1, 5, 5, 0.5, 0.5, 1, 0.1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 150, 150, 60, 60, 10, 10, 15, 5, 0.95]
            lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "SC2N2mG":

        # Custom Simple Secondary Contact Model: nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q
        func = modeledemo_new_models_folded.SC2N2mG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 0.8, 1, 1, 5, 5, 1, 0.1, 0.5, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12])

        
            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 200, 200, 1, 60, 60, 10, 10, 15, 7, 0.95, 0.5]
            lower_bound = [0.01, 0.01, 0.1, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.5, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "PIM2m":

        # Custom Periodic Isol. w. Migration with 2 Migration rate model: nu1, nu2, mA12, mA21, mAe12, mAe21, m12, m21, me12, me21, Ts, Tam, Tsc, P
        func = modeledemo_mis_new_models.PIM2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
            #nu1, nu2, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts, Tam, Tsc, P
                params = (1, 1, 5, 5, 0.5, 0.5, 5, 5, 0.5, 0.5, 1, 0.3, 0.1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 50, 50, 10, 10, 50, 50, 10, 10, 15, 10, 5, 0.95]
            lower_bound = [0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "PIM2N":

        # Custom Periodic Isol. w. Migration with 2 Migration rate and linked sel.: nu1, nu2, hrf, mA12, mA21, m12, m21, Ts, Tam, Tsc, Q
        func = modeledemo_mis_new_models.PIM2N

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
            #nu1, nu2, hrf, mA12, mA21, m12, m21, Ts, Tam, Tsc, P
                params = (1, 1, 0.8, 5, 5, 5, 5,  1, 0.3, 0.1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 1, 50, 50, 50, 50, 15, 10, 5, 0.95]
            lower_bound = [0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "IM2m":

        # Custom Isolation with 2 Migration rate model: nu1, nu2, m12, m21, me12, me21, Ts, P
        func = modeledemo_new_models_folded.IM2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 5, 5, 0.5, 0.5, 1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])
        
            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [250, 250, 150, 150, 30, 30, 15, 0.95]
            lower_bound = [0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "AM2m":

        # Custom Ancient Migration with 2 Migration rate model: nu1, nu2, m12, m21, me12, me21, Ts, Tam, P
        func = modeledemo_new_models_folded.AM2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 5, 5, 0.5, 0.5, 1, 0.1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [250, 250, 150, 150, 30, 30, 15, 5, 0.95]
            lower_bound = [0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "SC2m":

        # Custom Secondary contact with 2 Migration rate model: nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P
        func = modeledemo_new_models_folded.SC2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":       
                params = (1, 1, 5, 5, 0.5, 0.5, 1, 0.1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [250, 250, 150, 150, 30, 30, 15, 7, 0.95]
            lower_bound = [0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "PSC2m":

        # Custom Periodic SC with 2 Migration rate model: nu1, nu2, mA12, mA21, mAe12, mAe21, m12, m21, me12, me21, Ts, Tsc1, Tam, Tsc, P
        func = modeledemo_mis_new_models.PSC2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
            #nu1, nu2, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts,Tsc1, Tam, Tsc, P
                params = (1, 1, 5, 5, 0.5, 0.5, 5, 5, 0.5, 0.5, 1, 0.6, 0.3, 0.1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13], popt[14])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13], popt[14])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            #nu1, nu2, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts,Tsc1, Tam, Tsc, P
            upper_bound = [100, 100, 50, 50, 10, 10, 50, 50, 10, 10, 20, 15, 10, 5, 0.95]
            lower_bound = [0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "PSC2N":

        # Custom Periodic SC with 2 Migration rate model: nu1, nu2, hrf, mA12, mA21, m12, m21, Ts, Tsc1, Tam, Tsc, P
        func = modeledemo_mis_new_models.PSC2N

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
            #nu1, nu2,hrf,  mA12, mA21, m12, m21, Ts,Tsc1, Tam, Tsc, P
                params = (1, 1, 0.8, 5, 5, 5, 5, 1, 0.6, 0.3, 0.1, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11]
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 1, 50, 50, 50, 50, 20, 15, 10, 5, 0.95]
            lower_bound = [0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))
        
    if namemodel == "PSC2N2m":

        # Custom Periodic SC with 2 Migration rate model: nu1, nu2, hrf, mA12, mA21, mAe12, mAe21, m12, m21, me12, me21, Ts, Tsc1, Tam, Tsc, P, Q
        func = modeledemo_mis_new_models.PSC2N2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
            #nu1, nu2, hrf, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts,Tsc1, Tam, Tsc, P, Q
                params = (1, 1, 0.8, 5, 5, 0.5, 0.5, 5, 5, 0.5, 0.5, 1, 0.6, 0.3, 0.1, 0.5, 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13], popt[14], popt[15], popt[16])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13], popt[14], popt[15], popt[16])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            #nu1, nu2, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts,Tsc1, Tam, Tsc, P
            upper_bound = [100, 100, 1, 50, 50, 10, 10, 50, 50, 10, 10, 20, 15, 10, 5, 0.95, 0.95]
            lower_bound = [0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0.05]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))


    if namemodel == "PIM2N2m":

        # Custom Periodic Isol. w. Migration with 2 Migration rate model: nu1, nu2, mA12, mA21, mAe12, mAe21, m12, m21, me12, me21, Ts, Tam, Tsc, P
        func = modeledemo_mis_new_models.PIM2N2m

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
            #nu1, nu2, hrf, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts, Tam, Tsc, P, Q
                params = (1, 1, 0.8, 5, 5, 0.5, 0.5, 5, 5, 0.5, 0.5, 1, 0.3, 0.1, 0.5 0.5)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13], popt[14], popt[15])
            else :
                params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11], popt[12], popt[13]), popt[14], popt[15]

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [100, 100, 1, 50, 50, 10, 10, 50, 50, 10, 10, 15, 10, 5, 0.95, 0.95]
            lower_bound = [0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0.95]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic, 
                                  nompop1=nompop1, nompop2=nompop2, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))



output_file.close()
