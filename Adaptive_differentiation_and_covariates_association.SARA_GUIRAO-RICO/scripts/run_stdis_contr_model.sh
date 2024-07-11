#!/bin/bash     

# run BayPass (STDis and Contrast Model)
mkdir results_stdis_contr
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./input/hgdp.geno -contrastfile ./input/covariates_eu -efile ./input/covariates_eu -outprefix ./results_stdis_contr/hgdp_contrast
