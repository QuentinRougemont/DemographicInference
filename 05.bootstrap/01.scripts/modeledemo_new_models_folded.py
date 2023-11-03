#!/usr/bin/env python
# -*- coding: utf-8 -*- 

"""
module modeledemo with demographic models of divergence
"""

import numpy
import dadi

def SC2N2mG(params, ns, pts):
    nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q = params
    n1, n2 = ns

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
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

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
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

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
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

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
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI 
    return fs



def IM2N2mG(params, ns, pts):
    nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, P, Q = params
    n1, n2 = ns

    """
    Model of semi permeability with split, ongoing migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

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
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
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
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    #### Spectrum of low-recombining regions

    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the spectra
    fs = (P*fsN+(1-P)*fsI+(1-Q)*fsnr+Q*fslr) 
    return fs


def AM2mG(params, ns, pts):
    nu1, nu2, b1, b2, m12, m21, me12, me21, Tam, Ts, P = params
    n1, n2 = ns
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
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

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
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion O (and P)
    fs =(P*fsN + (1-P)*fsI) 
    return fs
    

def AM2N2mG(params, ns, pts):
    nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q = params
    n1, n2 = ns

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
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

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
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    
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
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

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
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
 
    #### Sum the spectra
    fs = (P*fsN + (1-P)*fsI + (1-Q)*fsnr + Q*fslr) 
    return fs
