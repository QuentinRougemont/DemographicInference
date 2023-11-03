# DemographicInference
This repository contains files and scripts from dadi v1.6.3 (Gutenkunst et al. 2009) that were modified to improve optimization, as well as custom models that account for semi-permeability (i.e by including variable effective migration rates across the genome), as introduced by Tine et al 2014. The repository also contains data from Rougemont et al. 2016 who applied these models to European lamprey ecotypes using both Folded and Unfolded JSFS. The approach implemented in dadi relies on the comparison of alternative demographic models of divergence to genome-wide polymorphism data.



[![DOI](https://zenodo.org/badge/54107776.svg)](https://zenodo.org/badge/latestdoi/54107776)


## TO DO: 

make everything compatible with python3 ....


## WARNING:

This pipeline requires `scipy` v0.13.3 or older versions. 

The simulated annealing implemented in  `model_demo.py` is deprecated in newer version of scipy.


### Dependencies:

`python2.7`

`scipy v0.13` or older

`numpy`

`matplotlib`

`pylab`

To extract population SFS using `utility_scripts`:
`vcftools`
`bgzip`
`tabix`

### REFERENCES:
 
Gutenkunst RN, Hernandez RD, Williams SH, Bustamante CD (2009) Inferring the joint demographic history of multiple populations from multidimensional SNP frequency data. PLoS Genet 5: e1000695.

Tine M, Kuhl H, Gagnaire P-A et al. (2014) European sea bass genome and its variation provide insights into adaptation to euryhalinity and speciation. Nature Communications, 5:5770,doi:10.1038/ncomms6770

Rougemont, Q., Gagnaire, P.-A., Perrier, C., Genthon, C., Besnard, A.-L., Launey, S. and Evanno, G. (2016), Inferring the demographic history underlying parallel genomic divergence among pairs of parasitic and non-parasitic lamprey ecotypes. Mol Ecol. Accepted Author Manuscript. https://doi:10.1111/mec.13664

Hugo Cayuela\*, Quentin Rougemont\*, Martin Laporte, Claire Mérot, Eric Normandeau, Yann Dorant, Ole K. Tørresen, Siv Nam Khang Hoff, Sissel Jentoft, Pascal Sirois, Martin Castonguay, Teunis Jansen, Kim Praebel, Marie Clément, Louis Bernatchez
Shared ancestral polymorphism and chromosomal rearrangements as potential drivers of local adaptation in a marine fish
https://doi.org/10.1111/mec.15499
bioRxiv 782201; doi: https://doi.org/10.1101/782201  

Rougemont Q, Moore J-S, Leroy T, Normandeau E, Rondeau EB, Withler RE, et al. (2020) 
Demographic history shaped geographical patterns of deleterious mutation load in a broadly distributed Pacific Salmon. PLoS Genet 16(8):e1008348. 
https://doi.org/10.1371/journal.pgen.1008348

Population genomics and history of speciation reveal fishery management gaps in two related redfish species (Sebastes mentella and Sebastes fasciatus)” (2020)
Laura Benestan, Quentin Rougemont, Caroline Senay, Eric Normandeau, Eric Parent, Rick Rideout, Louis Bernatchez, Yvan Lambert, Céline Audet, Geneviève J. Parent 
Evoluationay application https://doi.org/10.1111/eva.13143


See also: https://popgensealab.wordpress.com/dadi-inference/

### Dadi and Its modification:

The modified version of dadi contain several improvements

1. The possibility to explictly test alternative models of speciation such as models of divergence with continuous gene flow (IM models), models of secondary contact (SC models), models of strict alloptatric divergence (SI) or models of divergence with initial gene flow followed by isolation without gene flow (AM model). The models are:  
  * IM 
  * AM 
  * SI 
  * SC   

   and are implemented in the model_demo.py see [here](https://github.com/QuentinRougemont/DemographicInference/tree/master/02-modifs_v2)

the models are displayed below:
![example_graph](https://github.com/QuentinRougemont/DemographicInference/blob/master/pictures/Fig1.png)


2. The possibility to integrate semipermeability in the demographic models of divergence with gene flow (IM2m, AM2m) and secondary contacts (SC2m).

It is also possible to include various effect of linked selection (i.e. background selection, selective sweeps and hitchhiking effects) locally reducing Ne (suffix 2N).
These models can be combined to integrate both the effect of linked selection and of barriers to gene flow:
  * therefore the following model can be tested:
  * IM2N, AM2N, SC2N
  * IM2m, AM2m, SC2m
  * IM2N2m, AM2N2mn SC2N2m  

Modelling variable migration rates across the genome  is a simple way to capture the effect of genetic barriers to gene flow (due to genetic incompatibilities or local adaptation loci), that locally reduce gene flow and
generate heterogeneous patterns of genomic divergence. 

Here, two categories of loci are considered to take into account these variable rates of gene flow across the genome, those that are "neutrally" exchanged (occurring in proportion P) 
and those with reduced effective migration rates (occurring in proportion 1-P). 

The same idea is used to define loci undergoing higher rates of drift.

Details of the models are available here: `02-modifs_v2/unfolded/modeledemo_new_models.py`

It is also possible to model population expansion and contraction but I ussually do not use these models to avoid overfitting.


3. An improvement of the method that search for the global optimum of the likelihood function. 

Optimization is performed using a Simulated Annealing (SA) optimization procedure (implemented in the inference.py module) before BFGS optimization. 

In practice, one hot and one cold SA optimization are performed before  one additional round of BFGS optimization. 
See the script script_inference_anneal2_newton.py for further details. 

Of course, this does not preclude to perform multiple runs to check for convergence in likelihood estimations among independent runs.

4. The error rate for ancestral allelic state identification in unfolded JSFS is taken into account by a misorientation parameter (O).

These modifications were originally developed by Tine et al. (2014). Yoann Anciaux, Christelle Fraïsse, Pierre-Alexandre Gagnaire, and Nicolas Bierne (see Tine et al. 2014).

5. additional modification:

model with/without growth for all population:
 * IM/IMG, IM2N/IM2NG, IM2m/IM2mG, IM2N2m/IM2N2mG
 * AM/AMG, AM2N/AM2NG, AM2m/AM2mG, AM2N2m/AM2N2mG
 * SC/SCG, SC2N/SC2NG, SC2m/SC2mG, SC2N2m/SC2N2mG
 * SI/SIG, SI2N/SI2NG
 * PANG

![example_graph](https://github.com/QuentinRougemont/DemographicInference/blob/master/pictures/Fig2.png)



model with/without growth IN ANCESTRAL population:
 * IMA/IMAG, IMA2N/IMA2NG, IMA2m/IMA2mG, IMA2N2m/IMA2N2mG
 * AMA/AMAG, AMA2N/AMA2NG, AMA2m/AMA2mG, AMA2N2m/AMA2N2mG
 * SCA/SCAG, SCA2N/SCA2NG, SCA2m/SCA2mG, SCA2N2m/SCA2N2mG
 * SIA/SIAG, SIA2N/SIA2NG

![example_graph](https://github.com/QuentinRougemont/DemographicInference/blob/master/pictures/Fig3.png)


models with multiple phase of isolation contacts are also implemented and can be modified easily

models for one pop:
 * with instantaneous pop change (PAN)
 * with instantaneous pop change followed by exp. growth (PANG)
 * same with bottleneck (PANGb)
 * same with linked selection (suffix 2N)

 These are defined [here](https://github.com/QuentinRougemont/DemographicInference/tree/master/02-modifs_v2/folded) for the folded version  
 & [here](https://github.com/QuentinRougemont/DemographicInference/tree/master/02-modifs_v2/unfolded) for the unfolded version

 The script to define the prior and launch the one pop model is [here](https://github.com/QuentinRougemont/DemographicInference/blob/master/02-modifs_v2/pang_models.py)


## Unfolding and Folding:

### Unfolded JSFS

In dadi, when an outgroup is available, it is possible to determine the ancestral state of the mutations and then to count the number of derived mutations segregating in each population (See Gutenkunst et al. 2009 for more information). In the implementation by Tine et al., it is possible to infer the proportion of misoriented mutations, which is included as a model parameter (O). In Rougemont et al. (2016) the Sea lamprey outgroup (Petromyzon marinus) was a too distant to determine the ancestral state of all mutations, so the majority of mutations remained unoriented motivating the use of the folded JSFS.

### Folded JSFS

The folded (J)SFS includes all SNPs without distinguishing ancestral and derived alleles. The folded SFS is therefore less informative with regards to the directionality of migration and drift effects. However, it is still informative enough to make inferences using simple demographic divergence models as the ones presented before. In particular, the folded JSFS can retain the signal of a secondary contact if the duration of divergence was long enough compared to the duration of secondary gene flow. Moreover, heterogeneous gene flow can be detected from an unfolded JSFS, using models incorporating variation in migration rate across the genome.

To circumvent the loss of information and validate the accuracy of parameters estimates from the folded JSFS, our strategy was to re-integrate fixed values for parameters estimates such as Ne1, Ne2, T, Tsc in the unfolded version of the JSFS and then to re-estimate migration rate parameters (M12,M21,Me12,Me21) as well as the proportion of loci displaying a reduced migration rate (1-P). 

## Details:

Models are defined here:

`02-modifs_v2/unfolded/modeledemo_new_models.py`

and here for the folded version:

`02-modifs_v2/folded/modeledemo_new_models_folded.py`

The unfolded prior are in : 

please check the prior before running the models !!!  

`02-modifs_v2/unfolded/script_inference_demo_new_models.py`

and  the folded priors are in:

`02-modifs_v2/folded/script_inference_demo_new_models_folded.py`

Then the scripts:

`01-scripts/00.run_dadi_parallel_v2.sh`

and `01-run_model_iteration_v2.sh`

provides example for running all models in parallel. These should run easily everywhere if you have parallel module installed.

I used them on the ComputeCanada super computer(s) with examples of launcher provided in `graham_cedar/` folder and `colosse` folder

## To DO:  
add PSC, PIM,  PSC2N, PIM2N, PSC2N2m, PIM2N2m 
