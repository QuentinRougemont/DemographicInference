#! ~/py33/bin
# -*- coding: utf-8 -*-

"""module modeledemo contenant les différents modèles démographiques de divergence"""

import numpy
import dadi


def SI(params, (n1,n2), pts):
    nu1, nu2, Ts= params
    """
    Model with split and complete isolation.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
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
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
 
    return fs

def IM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts = params
    """
    Model with migration during the divergence.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
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
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


def AM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tam = params
    """
    Model with split, ancient migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
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
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def PAM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tam = params
    """
    Model with split, followed by two periods of ancient migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
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
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21 
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=0, m21=0)

    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def SC(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc = params
    """
    Model with split, complete isolation, followed by secondary contact

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
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
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # Finally, calculate the spectrum.
    # oriented
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def PSC(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc = params
    """
    Model with split, complete isolation, followed by two periods of secondary contact

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
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
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # Finally, calculate the spectrum.
    # oriented
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def IM2M(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, P = params
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
  
    ### Sum the two spectra in proportion P (and O)
    fs = P*fsNO+(1-P)*fsIO
    return fs

def AM2M(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tam, P = params
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
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
  
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
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = P*fsN+(1-P)*fsI
    return fs

def SC2M(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P = params
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
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))



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
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = P*fsN+(1-P)*fsI
    return fs

def SC3M(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, mi12, mi21, Ts, Tsc, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 arround genomic islands.
    me21: Effective migration from pop 1 to pop 2 arround genomic islands.
    mi12: Effective migration from pop 2 to pop 1 within genomic islands.
    mi21: Effective migration from pop 1 to pop 2 within genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome influenced by genomic islands
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
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
    
    ### Calculate the spectrum arround genomic islands
    # phi for the equilibrium ancestral population
    phiE = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiE = dadi.PhiManip.phi_1D_to_2D(xx, phiE)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiE = dadi.Integration.two_pops(phiE, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiE = dadi.Integration.two_pops(phiE, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsE = dadi.Spectrum.from_phi(phiE, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to mi12 and mi21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=mi12, m21=mi21)
    # calculate the spectrum.
    # oriented
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))


    ### Sum the three spectra in proportion P Q and (1-Q) (and O)
    fs = P*fsN+Q*fsE+(1-(P+Q))*fsI
    return fs

def IM2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, P1, P2 = params
    """
    Model with migration during the divergence with two type of migration and two proportions

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    P1: The proportion of the genome evolving neutrally in population 1
    P2: The proportion of the genome evolving neutrally in population 2
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    # oriented
    fsN1N2 = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsI1I2 = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))

    ### Calculate the neutral spectrum in population 1 and the genomic island spectrum in population 2
    # phi for the equilibrium ancestral population
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and me21
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=m12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsN1I2 = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum in population 1 and the neutral spectrum in population 2
    # phi for the equilibrium ancestral population
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and m21
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=me12, m21=m21)
    # calculate the spectrum.
    # oriented
    fsI1N2 = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))


    ### Sum the four spectra
    fs = P1*P2*fsN1N2 + (1-P1)*(1-P2)*fsI1I2 + P1*(1-P2)*fsN1I2 + (1-P1)*P2*fsI1N2
    return fs


def AM2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tam, P1, P2 = params
    """
    Model of semi permeability with split, complete isolation, followed by ancient migration with 2 migration rates and two proportions

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P1: The proportion of the genome evolving neutrally in population 1
    P2: The proportion of the genome evolving neutrally in population 2
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsN1N2 = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to me12 and me21
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsI1I2 = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
 
    ### Calculate the neutral spectrum in population 1 and the genomic island spectrum in population 2
    # phi for the equilibrium ancestral population
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and me21
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=m12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsN1I2 = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
 

    ### Calculate the genomic island spectrum in population 1 and the neutral spectrum in population 2
    # phi for the equilibrium ancestral population
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to me12 and m21
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=me12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsI1N2 = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))


    ### Sum the four spectra
    fs = P1*P2*fsN1N2 + (1-P1)*(1-P2)*fsI1I2 + P1*(1-P2)*fsN1I2 + (1-P1)*P2*fsI1N2
    return fs

def PAM2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tam, P1, P2 = params
    """
    Model of semi permeability with split, complete isolation, followed by two periods of ancient migration with 2 migration rates and two proportions

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P1: The proportion of the genome evolving neutrally in population 1
    P2: The proportion of the genome evolving neutrally in population 2
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to me12 and me21
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to me12 and me21
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))


    ### Calculate the neutral spectrum in population 1 and the genomic island spectrum in population 2
    # phi for the equilibrium ancestral population
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and me21
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=m12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and me21
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=m12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    
    ### Calculate the genomic island spectrum in population 1 and the neutral spectrum in population 2
    # phi for the equilibrium ancestral population
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to me12 and m21
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=me12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to me12 and m21
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=me12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    # oriented
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
 
    ### Sum the four spectra
    fs = P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O
    return fs

def SC2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P1, P2 = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates and two proportions

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P1: The proportion of the genome evolving neutrally in population 1
    P2: The proportion of the genome evolving neutrally in population 2
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    # oriented
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))

    ### Calculate the neutral spectrum in population 1 and the genomic island spectrum in population 2
    # phi for the equilibrium ancestral population
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and me21
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum in population 1 and the neutral spectrum in population 2
    # phi for the equilibrium ancestral population
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and m21
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1, nu2, m12=me12, m21=m21)
    # calculate the spectrum.
    # oriented
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))

    ### Sum the four spectra
    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) 
    return fs

def PSC2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P1, P2 = params
    """
    Model of semi permeability with split, complete isolation, followed by two periods of secondary contact with 2 migration rates and two proportions

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P1: The proportion of the genome evolving neutrally in population 1
    P2: The proportion of the genome evolving neutrally in population 2
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    # oriented
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum in population 1 and 2
    # phi for the equilibrium ancestral population
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))

    ### Calculate the neutral spectrum in population 1 and the genomic island spectrum in population 2
    # phi for the equilibrium ancestral population
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and me21
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=me21)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and me21
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum in population 1 and the neutral spectrum in population 2
    # phi for the equilibrium ancestral population
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and m21
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1, nu2, m12=me12, m21=m21)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and m21
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1, nu2, m12=me12, m21=m21) 
    # calculate the spectrum.
    # oriented
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))

    ### Sum the four spectra
    fs = P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O
    return fs

