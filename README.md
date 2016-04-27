# DemographicInference
This repository contains files and script from dadi v1.6.3 (Gutenkunst et al. 2009) as well as a modified version of the program introduced by Tine et al 2014 to account for  semipermeability (i.e to include variable levels of gene flow accorss the genome).
The repository also contains data from Rougemont et al. 2016 who applied the model to European lamprey ecotypes using both Folded and Unfolded JSFS. 
The idea behind dadi is to compare demographic models of genomic divergence.

REFERENCE:
Gutenkunst RN, Hernandez RD, Williams SH, Bustamante CD (2009) Inferring the joint demographic history of multiple populations from multidimensional SNP frequency data. PLoS Genet 5: e1000695.

Tine M, Kuhl H, Gagnaire P-A et al. (2014) European sea bass genome and its variation provide insights into adaptation to euryhalinity and speciation. Nature Communications, 5:5770,doi:10.1038/ncomms6770

Rougemont, Q., Gagnaire, P.-A., Perrier, C., Genthon, C., Besnard, A.-L., Launey, S. and Evanno, G. (2016), 
Inferring the demographic history underlying parallel genomic divergence among pairs of parasitic and non-parasitic lamprey ecotypes. Mol Ecol. Accepted Author Manuscript. doi:10.1111/mec.13664

See also:https://popgensealab.wordpress.com/dadi-inference/

#################### Dadi and Its modification ###############################################
The modified version of dadi contain 3 improvements:

1. The possibility to explictly test alternative models of speciation such as models of divergence with continuous gene flow (IM models), models of secondary contact (SC models), models of strict alloptatric divergence (SI) or models of divergence with initial gene flow followed by allopatric divergence or stop of gene flow (AM model). The models are implemented in the model_demo.py

2. The possibility to integrate semipermeability in the demographic model of divergence with gene flow (IM2m, AM2m) and secondary contacts (SC2m).
Such idea was first implemented by Roux et al (2013) in an ABC context. Also it is very accurate, it is more time consuming.
Incorporating varying migration rate is howver important as species boundary are long known to be semipermeable. 
Accordingly, modelling variable rates of gene flow accross the genome can help taking into account the effect of barrier to gene flow (due to selection or other processes), that locally reduced gene flow, 
and generate heterogeneous patterns of genomic divergence. Here, two categories of loci, those "neutrally" exchanged (P) and those with a non neutral rates of exchange (1-P) are considered to take into account
 these variable rates of gene flow accross the genome. The models are also implemented in the model_demo.py

3. An optimisation of the search for the global optimum leading to improved inferences of the parameters and likelihood. 
Optimisation is performed using a Simulated Annealing (SA) optimisation procedure (implemented in the inference.py module).
In the inference, a hot and a cold SA followed by one round of BFGS optimisation. See the script script_inference_anneal2_newton.py for further details. 
Of course, this does not preclude to perform multiple runs to check for convergence in likelihood estimations.

Note importantly that all these modifications were originally developped by Tine et al. (2014).

############################Unfolding and Folding ########################################
#Unfolded JSFS
In dadi when an outgroup is available, it is possible to determine the ancestral state of the mutations and then to count the number of derived mutations segregating in one, two or three derived populations. (See Gutenkunst et al. 2009 for more information).
In the implementation by Tine et al. it is possible to determine the proportion of misoriented mutation, which result in a distortion of the SFS, this is included as a model parameters and 
In Rougemont et al. (2016) the outgroup Sea lamprey (Petromyzon marinus) was a too distant genus so that it was not possible to determine the ancestral state, and the majority of mutations remained unoriented in this case.

#Folded JSFS
As a result, it is possible to use the folded (J)SFS that includes all SNPs. 
The folded SFS was less informative with regards to estimations of the directionnality of parameters such as migration rates. However, it was still informative to largely improves inferences. In particular we showed how incorporating variation in migration rate accross the genome improved inferences.

To circumvent the loss of information and validate the accuracy of parameters estimates from the folded JSFS, our strategy was to re-integrate fixed values for parameters estimates such as Ne1, Ne2, T,Tsc in the unfolded version of the JSFS and then to re-estimate muigration rate (M12,M21) effective migration rate (Me12,Me21) as well as the proportion of loci displaying a reduced migration rate (1-P). 
