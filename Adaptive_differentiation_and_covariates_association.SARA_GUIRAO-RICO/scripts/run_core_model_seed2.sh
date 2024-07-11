#!/bin/bash

# run BayPass (CORE Model) with seed2
mkdir results_core
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./input/hgdp.geno -seed 26847 -outprefix ./results_core/hgdp_core_s2

