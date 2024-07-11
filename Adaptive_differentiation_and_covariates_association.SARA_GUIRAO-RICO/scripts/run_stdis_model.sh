#!/bin/bash                                                                                                                         

# run BayPass (STDis Model)
mkdir results_standard
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./input/hgdp.geno -efile ./input/covariates -scalecov -outprefix ./results_standard/hgdp_stdis
