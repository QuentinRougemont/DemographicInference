#! ~/py33/bin
# -*- coding: utf-8 -*- 

"""
module modeledemo with demographic models of divergence used in Rougeux et al.2017
"""

import numpy
import dadi


def SI(params, (n1,n2), pts):
    nu1, nu2, Ts, O = params
    """
    Model with split and complete isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    O: The proportion of accurate orientation
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Finally, calculate the spectrum.
    # oriented
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # mis-oriented
    fsM = dadi.Numerics.reverse_array(fsO)

    ### Sum the two spectra in proportion O
    fs = O*fsO+(1-O)*fsM
    return fs


def SI2N(params, (n1,n2), pts):
    nu1, nu2, Ts, nr, bf, O = params
    """
    Model with split and complete isolation, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: The scaled time of the split
    n1,n2: Size of fs to generate.
    nr: Proportion of non/low-recombining regions
    bf : Background factor (to which extent the effective population size is reduced in the non-recombining regions)
    O: The proportion of accurate orientation
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # Spectrum of non-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2	
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    # Finally, calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    
    
	# Spectrum of recombining regions
    # phi for the equilibrium ancestral population
    phir = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    # We set the population sizes after the split to nu1 and nu2
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Finally, calculate the spectrum.
    # oriented
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    # mis-oriented
    fsrM = dadi.Numerics.reverse_array(fsrO)

    ### Sum the two spectra in proportion O
    fs= O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)

    return fs


def SIG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, Ts, O = params
    
    """ 
    Model with split and complete isolation.

    nu1: Size of population 1 at split.
    nu2: Size of population 2 at split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We start the population size change after the split independantly in each population
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # mis-oriented
    fsM = dadi.Numerics.reverse_array(fsO)

    ### Sum the two spectra in proportion O
    fs = O*fsO+(1-O)*fsM
    return fs


def SI2NG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, Ts, Q, O = params
    """
    Model with split and complete isolation, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    Ts: The scaled time of the split
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We start the population size change after the split independantly in each population
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) and integrate the hrf for low-recombining regions
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)
    
    ### Sum the two spectra in proportion O (and Q)
    fs= O*((1-Q)*fsnrO + Q*fslrO) + (1-O)*((1-Q)*fsnrM + Q*fslrM)
    return fs


def IM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, O = params
    """
    Model with migration during the divergence.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    O: The proportion of accurate orientation
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # Finally, calculate the spectrum.
    # oriented
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # mis-oriented
    fsM = dadi.Numerics.reverse_array(fsO)

    ### Sum the two spectra in proportion O
    fs = O*fsO+(1-O)*fsM
    return fs


def IMG(params, (n1,n2), pts): 
    nu1, nu2, b1, b2, m12, m21, Ts, O = params
    """
    Model with migration during the divergence.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split (in units of 2*Na generations).
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We start the population size change after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # mis-oriented
    fsM = dadi.Numerics.reverse_array(fsO)
    
    ### Sum the two spectra in proportion O
    fs = O*fsO+(1-O)*fsM
    return fs

def IM2N(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, Ts, Q, O = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    #### Sum the two spectra in proportion O (and Q)
    fs= O*((1-Q)*fsnrO + Q*fslrO) + (1-O)*((1-Q)*fsnrM + Q*fslrM)
    return fs


def IM2NG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, Ts, Q, O = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We start the population size change after the split independantly in each population and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) & integrate the hrf for low-recombining regions and set the migration rates to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b1**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    #### Sum the two spectra in proportion O (and Q)
    fs= O*((1-Q)*fsnrO + Q*fslrO) + (1-O)*((1-Q)*fsnrM + Q*fslrM)
    return fs


def IM2mG(params, (n1,n2), pts): 
    nu1, nu2, b1, b2, m12, m21, me12, me21, Ts, P, O = params
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    ### Sum the two spectra in proportion O (and P)
    fs = O*(P*fsNO + (1-P)*fsIO) + (1-O)*(P*fsNM + (1-P)*fsIM)
    return fs


def AM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tam, O = params

    """
    Model with split, ancient migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    O: The proportion of accurate orientation
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21 
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=0, m21=0)

    # Finally, calculate the spectrum.
    # oriented
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # mis-oriented
    fsM = dadi.Numerics.reverse_array(fsO)

    ### Sum the two spectra in proportion O
    fs = O*fsO+(1-O)*fsM
    return fs


def AMG(params, (n1,n2), pts): 
    nu1, nu2, b1, b2, m12, m21, Tam, Ts, O = params
    """
    Model with split, ancient migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
  	b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
    Ts: The scaled time between the end of ancient migration and present.
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# We start the population size change after the split independantly in each population and set the migration rate to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tam)
    bnu2_func = lambda t: nu2 * b2**(t/Tam)
    phi = dadi.Integration.two_pops(phi, xx, Tam, bnu1_func, bnu2_func, m12=m12, m21=m21)
    # We continue the population size change after ancient migration (until present) independantly in each population and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    # oriented
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # mis-oriented
    fsM = dadi.Numerics.reverse_array(fsO)

    ### Sum the two spectra in proportion O
    fs = O*fsO+(1-O)*fsM
    return fs

def AM2N(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, Tam, Ts, Q, O = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration.
    Ts: The scaled time between the end of ancient migration and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
	
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We keep the population sizes after the split to hrf*nu1 and hrf*nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    #### Sum the two spectra in proportion O and 1-O
    fs= O*((1-Q)*fsnrO + Q*fslrO) + (1-O)*((1-Q)*fsnrM + Q*fslrM)
    return fs
    
    
def AM2N2m(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q, O = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rate to zero 
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We keep the population sizes after the split to hrf*nu1 and hrf*nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    #### Sum the spectra
    fs = O*(P*fsNO+(1-P)*fsIO+(1-Q)*fsnrO+Q*fslrO) + (1-O)*(P*fsNM+(1-P)*fsIM+(1-Q)*fsnrM+Q*fslrM)
    return fs


def AM2NG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, Tam, Ts, Q, O = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration.
    Ts: The scaled time between the end of ancient migration and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We start the population size change after the split independantly in each population and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We start the population size change after the split independantly (bnu{1,2}_func) & integrate the hrf for low-recombining regions in each population and set the migration rates to zero
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b1**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    #### Sum the two spectra in proportion O and 1-O
    fs= O*((1-Q)*fsnrO + Q*fslrO) + (1-O)*((1-Q)*fsnrM + Q*fslrM)
    return fs
    

def AM2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, m12, m21, me12, me21, Tam, Ts, P, O = params
    """
    Model of semi permeability with split, ancient migration with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
    Ts: The scaled time between the end of ancient migration and present.
    P: The porportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes nu1 & nu2 after the split and set the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We start the population reduction after the split and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes nu1 & nu2 after the split and set the migration rate to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    # We start the population reduction after the split and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    ### Sum the two spectra in proportion O (and P)
    fs = O*(P*fsNO + (1-P)*fsIO) + (1-O)*(P*fsNM + (1-P)*fsIM)
    return fs
    

def AM2N2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q, O = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We SET the population size nu1 & nu2 after the split and set the migration rates to m12 & m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We set the population sizes change independently in each population after ancient migration to bnu1_func and bnu2_func and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population size after the split in each population and set the migration rates to me12 & me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes change after ancient migration to bnu1_func and bnu2_func and set the migration rates to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
	
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to me12 and me21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We set the population sizes after the split to hrf*bnu1_func and hrf*bnu2_func and set the migration rates to zero
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)
 
 	#### Sum the spectra
    fs = O*(P*fsNO + (1-P)*fsIO + (1-Q)*fsnrO + Q*fslrO) + (1-O)*(P*fsNM + (1-P)*fsIM + (1-Q)*fsnrM + Q*fslrM)
    return fs


def SC(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, O = params
    """
    Model with split, complete isolation, followed by secondary contact

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    n1,n2: Size of fs to generate.
    O: The proportion of accurate orientation
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # Finally, calculate the spectrum.
    # oriented
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # mis-oriented
    fsM = dadi.Numerics.reverse_array(fsO)

    ### Sum the two spectra in proportion O
    fs = O*fsO+(1-O)*fsM
    return fs


def SC2N(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, Ts, Tsc, Q, O = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    ### Sum the spectra 
    fs = O*((1-Q)*fsnrO + Q*fslrO) + (1-O)*((1-Q)*fsnrM + Q*fslrM)
    return fs


def SCG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, m12, m21, Ts, Tsc, O = params
   
    """
    Model with split, complete isolation, followed by secondary contact with exponential growth

    nu1: Size of population 1 at split.
    nu2: Size of population 2 at split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population reduction after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phi = dadi.Integration.two_pops(phi, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Calculate the spectrum
    # Oriented
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # Mis-oriented
    fsM = dadi.Numerics.reverse_array(fsO)
    
    # calculate the spectrum.
    fs = O*fsO+(1-O)*fsM
    return fs
    
    
def SC2N2m(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q, O = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    ### Sum the spectra 
    fs = O*(Q*fslrO+(1-Q)*fsnrO+P*fsNO+(1-P)*fsIO) + (1-O)*(Q*fslrM+(1-Q)*fsnrM+P*fsNM+(1-P)*fsIM)
    return fs


def SC2NG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, Ts, Tsc, Q, O = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We set the population sizes changes independently between pops after the split and isolation to nu1 and nu2 and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)   
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We set the population sizes changes independently between pops after the split and isolation to bnu1 and bnu2 (bnu{1,2}_func) & integrate the hrf and set the migration rates to m12 & m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tsc)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tsc)) * hrf   
    philr = dadi.Integration.two_pops(philr, xx, Tsc, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    ### Sum the spectra 
    fs = O*((1-Q)*fsnrO+Q*fslrO) + (1-O)*((1-Q)*fsnrM+Q*fslrM)
    return fs


def SC2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, m12, m21, me12, me21, Ts, Tsc, P, O = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates and exponential growth

    nu1: Size of pop 1 after split.
    nu2: Size of pop 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P: The porportion of the genome evolving neutrally.
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population size change after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population size change after the split and set the migration rates to me12 and me21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    ### Sum the two spectra in proportion P (and O)
    fs = O*(P*fsNO+(1-P)*fsIO) + (1-O)*(P*fsNM+(1-P)*fsIM)
    return fs


def SC2N2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q, O = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population sizes change independently in each populations after the split to bnu1 and bnu2 and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes change after the split to bnu1 and bnu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes change independently in each populations after the isolation to bnu1 and bnu2 and set the migration rates to m12 & m21
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # mis-oriented
    fsnrM = dadi.Numerics.reverse_array(fsnrO)

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tsc)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tsc)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Tsc, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    # oriented
    fslrO = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    # mis-oriented
    fslrM = dadi.Numerics.reverse_array(fslrO)

    ### Sum the spectra 
    fs = O*(Q*fslrO+(1-Q)*fsnrO+P*fsNO+(1-P)*fsIO) + (1-O)*(Q*fslrM+(1-Q)*fsnrM+P*fsNM+(1-P)*fsIM)
    return fs


def IM2m(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, P, O = params
    
    """
    nu1 = 1.6563
    nu2 = 0.7078
    m12 = 9.9748
    m21 = 0.2344
    me12 = 0.3558
    me21 = 1.3781
    Ts = 8.2590
    P = 0.8261
    O = 0.9312
    n1 = 26
    n2 = 26
    pts = 5
    """
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    ### Sum the two spectra in proportion P (and O)
    fs = O*(P*fsNO+(1-P)*fsIO) + (1-O)*(P*fsNM+(1-P)*fsIM)
    return fs

def AM2m(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tam, P, O = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)

  
    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    ### Sum the two spectra in proportion P (and O)
    fs = O*(P*fsNO+(1-P)*fsIO) + (1-O)*(P*fsNM+(1-P)*fsIM)
    return fs

def SC2m(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P, O = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P: The proportion of the genome evolving neutrally
    O: The proportion of accurate orientation
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    # oriented
    fsNO = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    # mis-oriented
    fsNM = dadi.Numerics.reverse_array(fsNO)


    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsIO = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    # mis-oriented
    fsIM = dadi.Numerics.reverse_array(fsIO)

    ### Sum the two spectra in proportion P (and O)
    fs = O*(P*fsNO+(1-P)*fsIO) + (1-O)*(P*fsNM+(1-P)*fsIM)
    return fs
