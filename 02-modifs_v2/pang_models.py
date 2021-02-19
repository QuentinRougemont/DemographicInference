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
        "# -o pathoutput -y population1 -p 10,20,30 -f pathfsfile -m PAN, PANG, PANGb, PANG2N, PANGb2N -l -a -h -v\n\n"+
        "# This is an exemple of the shortest command line:\n"+
        "# -f pathfsfile\n\n"+
        "# -h --help : Display the help you are looking at.\n"+
        "# -v --verbose : Print steps while the code is running\n"+
        "# -o --outputname : Take the path of output file.\n"+
        "# -f --fs_file_name : Take the path of the fs file from the parent directory.\n"+
        "# -p --grid_points : Take 3 numbers separated by a coma, for the size of grids for extrapolation.\n"+
        "# -m --model_list : Take until 18 name of model (PAN, PANG, PANGb, PANG2N, PANGb2N)separated by a coma.\n"+
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
    model_list = ["PAN","PANG","PANGb","PANG2N","PANGb2N"]
    verbose = False
    logparam = False
    nompop1 = "Pop1"

    checkfile = False #initilization. if True fs file needed exists, if False it doesn't

    if len(argv) < 2:
        print("You should give, at least, the name of the fs file !")
        sys.exit(1)
    try:
        opts, args = getopt.getopt(argv[1:], "hvo:y:azf:p:m:l", ["help", "verbose", "outputname=", "population1=", "masked", "fs_file_name=", "grid_points=", "model_list=", "log"])
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
    return(masked, pts_l, outputname, nompop1, fs_file_name, model_list, verbose, logparam)

#Inference function
def callmodel(func, data, output_file, modeldemo, ll_opt_dic, nbparam_dic,
          nompop1="Pop1", params=None, fixed_params=None, lower_bound=None, upper_bound=None,
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
        dadi.Plotting.plot_1d_comp_multinom(model, data,fig_num=None, residual="Anscombe" ,
                                  plot_masked=False, saveplot=True, nomplot=(outputname + "_" + modeldemo)) #, showplot=False)
    done=True
    return(done, ll_opt_dic, nbparam_dic, popt)

#plot_1d_comp_multinom(model, data, fig_num=None, residual='Anscombe',
#                          plot_masked=False):

##############################
##############################
# Load parameters
masked, pts_l, outputname, nompop1, fs_file_name, model_list, verbose, logparam = takearg(sys.argv)

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
    data.mask[1] = True
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
    if namemodel == "PAN": 

        # two epoch model: nu1,T
	#warning: this corresponds to the two epoch model of dadi
        func = modeledemo_new_models_folded.PAN

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1,1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1])
            else:
                params = (popt[0], popt[1])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [200, 20]
            lower_bound = [0.01, 0.001]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "PANG": 

        # panmictic model: nu1, T
        func = modeledemo_new_models_folded.PANG

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1])
            else:
                params = (popt[0], popt[1])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [200, 40 ]
            lower_bound = [0.01,0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "PANGb": 

        # panmictic model: nu1, nuF, T
        func = modeledemo_new_models_folded.PANGb

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
            upper_bound = [200, 200, 20 ]
            lower_bound = [0.001, 0.001, 0.0001]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "PANG2N": 

        # panmictic model: nu1, hrf, T, Q
        func = modeledemo_new_models_folded.PANG2N

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1,  1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3])
            else:
                params = (popt[0], popt[1], popt[2], popt[3])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [200,  1, 20, 0.5 ]
            lower_bound = [0.01, 0.01, 0.0001, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

    if namemodel == "PANGb2N": 

        # panmictic model: nuF, nuF,hrf, T, Q
        func = modeledemo_new_models_folded.PANGb2N

        for optimizationstate in opt_list:
            print optimizationstate

            if optimizationstate == "anneal_hot":
                params = (1, 1, 1, 1, 0.1)
            elif optimizationstate == "anneal_cold":
                params = (popt[0], popt[1], popt[2], popt[3], popt[4])
            else:
                params = (popt[0], popt[1], popt[2], popt[3], popt[4])

            # The upper_bound array is for use in optimization. Occasionally the optimizer
            # will try wacky parameter values. We in particular want to exclude values with
            # very long times, as they will take a long time to evaluate.
            upper_bound = [200, 200, 1, 20, 0.5]
            lower_bound = [0.01, 0.01, 0.01, 0.0001, 0.01]

            done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, namemodel, ll_opt_dic, nbparam_dic,
                                  nompop1=nompop1, params=params, fixed_params=None, lower_bound=lower_bound, 
                                  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
                                  outputname=outputname + "/" + outputname, 
                                  verbose=verbose, maxiter=20, Tini=50, Tfin=0, learn_rate=0.005, 
                                  schedule= "cauchy")
        if done: print(("\n" + namemodel + " : done\n"))

