# DemographicInference
This repository contains files and scripts from dadi v1.6.3 (Gutenkunst et al. 2009) that were modified to improve optimization, as well as custom models that account for semi-permeability (i.e by including variable effective migration rates across the genome), as introduced by Tine et al 2014. The repository also contains data from Rougemont et al. 2016 who applied these models to European lamprey ecotypes using both Folded and Unfolded JSFS. The approach implemented in dadi relies on the comparison of alternative demographic models of divergence to genome-wide polymorphism data.

## WARNING:

This pipeline requires `scipy` v0.13.3 or older versions. 

The simulated annealing implemented in  `model_demo.py` is deprecated in newer version of scipy.


Dependencies:
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

Rougemont, Q., Gagnaire, P.-A., Perrier, C., Genthon, C., Besnard, A.-L., Launey, S. and Evanno, G. (2016), Inferring the demographic history underlying parallel genomic divergence among pairs of parasitic and non-parasitic lamprey ecotypes. Mol Ecol. Accepted Author Manuscript. doi:10.1111/mec.13664

See also: https://popgensealab.wordpress.com/dadi-inference/

### Dadi and Its modification:

The modified version of dadi contain 4 improvements:

1. The possibility to explictly test alternative models of speciation such as models of divergence with continuous gene flow (IM models), models of secondary contact (SC models), models of strict alloptatric divergence (SI) or models of divergence with initial gene flow followed by isolation without gene flow (AM model). The models are implemented in the model_demo.py

2. The possibility to integrate semipermeability in the demographic models of divergence with gene flow (IM2m, AM2m) and secondary contacts (SC2m). Modelling variable migration rates across the genome  is a simple way to capture the effect of genetic barriers to gene flow (due to genetic incompatibilities or local adaptation loci), that locally reduce gene flow and generate heterogeneous patterns of genomic divergence. Here, two categories of loci are considered to take into account these variable rates of gene flow across the genome, those that are "neutrally" exchanged (occurring in proportion P) and those with reduced effective migration rates (occurring in proportion 1-P). The models are also implemented in the script model_demo.py

3. An improvement of the method that search for the global optimum of the likelihood function. Optimization is performed using a Simulated Annealing (SA) optimization procedure (implemented in the inference.py module) before BFGS optimization. In practice, one hot and one cold SA optimization are performed before  one additional round of BFGS optimization. See the script script_inference_anneal2_newton.py for further details. Of course, this does not preclude to perform multiple runs to check for convergence in likelihood estimations among independent runs.

4. The error rate for ancestral allelic state identification in unfolded JSFS is taken into account by a misorientation parameter (O).

These modifications were originally developed by Tine et al. (2014). Yoann Anciaux, Christelle Fraïsse, Pierre-Alexandre Gagnaire, and Nicolas Bierne (see Tine et al. 2014).

## Unfolding and Folding:

### Unfolded JSFS

In dadi, when an outgroup is available, it is possible to determine the ancestral state of the mutations and then to count the number of derived mutations segregating in each population (See Gutenkunst et al. 2009 for more information). In the implementation by Tine et al., it is possible to infer the proportion of misoriented mutations, which is included as a model parameter (O). In Rougemont et al. (2016) the Sea lamprey outgroup (Petromyzon marinus) was a too distant to determine the ancestral state of all mutations, so the majority of mutations remained unoriented motivating the use of the folded JSFS.

### Folded JSFS

The folded (J)SFS includes all SNPs without distinguishing ancestral and derived alleles. The folded SFS is therefore less informative with regards to the directionality of migration and drift effects. However, it is still informative enough to make inferences using simple demographic divergence models as the ones presented before. In particular, the folded JSFS can retain the signal of a secondary contact if the duration of divergence was long enough compared to the duration of secondary gene flow. Moreover, heterogeneous gene flow can be detected from an unfolded JSFS, using models incorporating variation in migration rate across the genome.

To circumvent the loss of information and validate the accuracy of parameters estimates from the folded JSFS, our strategy was to re-integrate fixed values for parameters estimates such as Ne1, Ne2, T, Tsc in the unfolded version of the JSFS and then to re-estimate migration rate parameters (M12,M21,Me12,Me21) as well as the proportion of loci displaying a reduced migration rate (1-P). 

