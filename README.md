# GWAS-GS
Working code for GWAS assisted GS 
Demo code for GWAS assisted GS in a spring wheat population for five traits. This code was originally written to be ran within an HPC. Replicates can be modified to decrease running times.

GWAS assisted GS is accomplished with GBLUP,CBLUP, and RKHS methods. I have tested inclusion of the top one, three, or five more significant by P value GWAS results to be incorproated as fixed effects in prediction models. GWAS is only done with materials in the model training pool, not model testing pool, to avoid contaminating accuracy comparison.
