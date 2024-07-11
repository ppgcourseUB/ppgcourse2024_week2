#!/bin/bash                                                                                                             

# run BayPass (CORE Model) with the 10000 PODs as input
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./results_core/G.hgdp_pods_10000 -efile ./input/covariates -scalecov  -outprefix ./results_standard/hgdp_stdis_10000_pods 
