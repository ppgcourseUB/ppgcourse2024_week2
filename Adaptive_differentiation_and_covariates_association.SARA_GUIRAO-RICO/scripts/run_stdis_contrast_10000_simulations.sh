#!/bin/bash                                                                                                             

# run BayPass (CORE Model) with the 10000 c2 PODs as input
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./results_stdis_contr/G.hgdp_C2_10000_pods -contrastfile ./input/covariates_eu -efile ./input/covariates_eu -outprefix ./results_stdis_contr/hgdp_contrast_10000_pods 
