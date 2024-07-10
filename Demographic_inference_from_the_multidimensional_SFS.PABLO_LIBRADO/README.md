# Demographic inference from the multidimensional SFS: Momi2

This session will showcase how to run momi2 based on a publicly-available dataset (Allentoft et al. 2024; https://erda.ku.dk/archives/917f1ac64148c3800ab7baa29402d088/published-archive.html), which includes thousands of ancient human genomes. For computational reasons, the dataset has been filtered to include SNPs from chromosome 1 only, and 5 ancient samples of **decent** quality:
 
| Sample | Date (ya) | Depth (X) | Continent |
|--|--|--|--|
| ela01 | 493 | 13.47 | Southern Africa |
| UstIshim | 45,020 |35.19 | Asia (Russia) |
| sf12 | 8,895 | 58.87| Europe (Sweeden)
| Funadomari_23 | 3,755 | 39.44 | Asia (Japan)
| AltaiNeandertal | ~120,000  | 50 | Asia

Download the folder for this practice in your running container.

```bash
ghget https://github.com/ppgcourseUB/ppgcourse2024_week2/tree/main/Demographic_inference_from_the_multidimensional_SFS.PABLO_LIBRADO
```

The input VCF file should contain 644,366 SNPs:

`gunzip -c 1.neo_impute_ph.filteredAA.bi.vcf.gz | grep -cv '#'` 


## Transform this VCF into a mSFS
Momi2 includes two scripts to calculate the multidimensional Site Frequency Spectrum (mSFS) from the VCF, resulting in the required input format:  

**Step1**:

`python -m momi.read_vcf --bed 20141020.strict_mask.whole_genome.chr1.bed 1.neo_impute_ph.filteredAA.bi.vcf.gz ind2pop.txt AlleleCounts.gz` 

Note that the file *ind2pop.txt* associates each sample to its corresponding population group (HG = Hunter Gatherer). The script polarizes the variants as ancestral or derived according to the AA value in the INFO field. The BED file contains the genomic coordinates that are trustworthy (eg. not affected by CNVs, showing good mappability to the hg19 reference assembly), and thus to be included in the subsequent analyses.

Open the AlleleCounts.gz file to explore its format:

`gunzip -c AlleleCounts.gz | less -S` 

**QUESTIONS**: Why there are 24,779 SNPs excluded? What do you think this line represents? 

> [[2, 0], [2, 0], [2, 0], [2, 0], [0, 0]],

**Solutions: reflects an allele configuration, characterised by 2 ancestral alleles in the African population, 2 ancestral alleles in HGAncientEurasia, and so on.* 

If you press *G*, you will be redirected to the end of the file, where each SNP has an associated allele configuration.


**Step2**:
To estimate the standard error of the observed mSFS, one can: 

 - Split chromosome 1 into *10* different blocks (equally-sized)
  - Calculate an independent SFS for each block
  - Bootstrap or jackknife across blocks:

`python -m momi.extract_sfs sfs.gz 10 AlleleCounts.gz` 

See the bottom of the resulting file *sfs.gz*. 

## Specifying our first momi2 model
As momi2 is a python module, we can open a python console, and run one command at a time, to understand what it does. 

First, we import the momi2 module within the python console, but also another module for logging out the messages to a file. We then define some very general and crucial evolutionary parameters (basically, ancient effective population size, generation time, and mutation rate).  

    import momi
    import logging
    logging.basicConfig(level=logging.INFO,
                    filename="tutorial.log")
    model = momi.DemographicModel(N_e=1.2e4, gen_time=29, muts_per_gen=1.25e-8)

**QUESTION**: Why do we need to define a mutation rate? and a Generation time? 

Then, we read the mSFS data, and link the data to a population model (yet to be defined, below):

    sfs = momi.Sfs.load("sfs.gz")
    model.set_data(sfs) 

Now we must add the samples as population representatives, and their associated population parameters. The parameter *N* represents their  effective population size, while *t* represents their respective time of sampling (in years ago), as stipulated in the first table of this tutorial.  

    model.add_leaf("Africa", N="N_Africa", t=493)
	model.add_leaf("HGAncientEurasia", N="N_HGAncientEurasia", t=45000)
	model.add_leaf("HGJapan", N="N_Japan", t=3755)
	model.add_leaf("HGSweeden", N="N_HGSweeden", t=8895)
	model.add_leaf("Neandertal", N="N_Neandertal", t=120000)

While *t* are fixed thanks to external information (eg. radiocarbon dating), we provided a name (string) to *Ns*, because these effective population sizes will need to be estimated.  Therefore, we need to create these parameters.

    model.add_size_param("N_Africa")
    model.add_size_param("N_HGAncientEurasia")
    model.add_size_param("N_Japan")
    model.add_size_param("N_HGSweeden")
    model.add_size_param("N_Neandertal")

After reaching this stage, it is time to revisit our knowledge about  human population history. We know that Neandertals are sister to all Anatomically Modern Humans (AMHs). Within AMHs, we know that Europeans and Asians are closer to each other than to Africans. Therefore, I propose the following initial topology (in Newick format):

    ((((HGJapan, HGSweeden), HGAncientEurasia), Africa), Neandertal);

in momi2, every internal node of this population tree needs to be specified manually. First node joins *HGSweeden* to *HGJapan* at time *t_HGJapan*:

    model.move_lineages("HGSweeden", "HGJapan", t="t_HGJapan", N="Na_HGJapan") 

Their ancestral population will have an effective population size of *Na_HGJapan* (parameter to be estimated as well). As *HGSweeden* joined *HGJapan*, and not viceversa, their ancestral population will be named *HGJapan*.  We proceed with the following split nodes:

    model.move_lineages("HGJapan", "HGAncientEurasia", t="t_HGAncientEurasia", N="Na_HGAncientEurasia") 
    model.move_lineages("HGAncientEurasia", "Africa", t="t_Africa", N="Na_Africa") 
    model.move_lineages("Africa", "Neandertal", t="t_root", N="Na_root")

We next need to create the corresponding parameters, to be estimated:

    model.add_size_param("Na_HGJapan")
    model.add_size_param("Na_HGAncientEurasia")
    model.add_size_param("Na_Africa")
    model.add_size_param("Na_root")
    
    model.add_time_param("t_HGJapan",lower=8895)
    model.add_time_param("t_HGAncientEurasia",lower=45000, lower_constraints=["t_HGJapan"])
    model.add_time_param("t_Africa",lower=450,lower_constraints=["t_HGAncientEurasia"])
    model.add_time_param("t_root", lower_constraints=["t_Africa"])

Note how we set a lowest boundary for *t_HGJapan*, because the split time between *HGSweeden* and *t_HGJapan* cannot be posterior the age of the oldest sample (*HGSweeden*). Likewise, we can also constrain the values for the other parameters.

**QUESTION**: Draw the population tree, and add the names of the parameters to the internal and external branches of the population tree.

Finally, fit the model (population tree) to the data (mSFS), using:

    model.optimize(method="TNC")
    
        message: Max. number of function evaluations reached
        success: False
         status: 3
            fun: 0.020349770964971054
              x: [ 1.154e+01  2.842e+00 ...  1.302e+05  5.680e+04]
            nit: 21
            jac: [ 9.806e-06 -1.831e-04 ...  1.043e-08 -7.782e-09]
           nfev: 131
     parameters:            N_Africa: 103038.4024389235
                  N_HGAncientEurasia: 17.143677483861627
                             N_Japan: 1103.1961644752325
                         N_HGSweeden: 2018.7414059526031
                        N_Neandertal: 695.3271571444569
                          Na_HGJapan: 6288.321217226941
                 Na_HGAncientEurasia: 8725.431559584407
                           Na_Africa: 3760.725889439739
                             Na_root: 18793.645754404668
                           t_HGJapan: 31154.933720726753
                  t_HGAncientEurasia: 45039.10885877785
                            t_Africa: 175224.0809136692
                              t_root: 232019.1424775751
                              
	 kl_divergence: 0.020349770964971054
	 log_likelihood: -1560946.8262688008


To evaluate fit, we can (i) check the likelihood value, and (ii) calculate diagnostic statistics, such as the observed (from the mSFS) vs expected (as predicted from the model fit) genetic distance between each pair of populations:

    model_fit_stats = momi.SfsModelFitStats(model)
    model_fit_stats.all_pairs_ibs()

                    Pop1              Pop2  Expected  Observed         Z
    0            HGJapan        Neandertal  0.552346  0.561391  1.773136
    1   HGAncientEurasia        Neandertal  0.559671  0.568641  1.515447
    2             Africa  HGAncientEurasia  0.629624  0.625349 -1.286339
    3          HGSweeden        Neandertal  0.553259  0.557463  1.022765
    4         Neandertal        Neandertal  0.962486  0.965463  0.731525
    5   HGAncientEurasia           HGJapan  0.707396  0.704530 -0.701925
    6            HGJapan         HGSweeden  0.711731  0.708436 -0.650945
    7             Africa            Africa  0.631729  0.630060 -0.507527
    8             Africa           HGJapan  0.622298  0.623359  0.492626
    9   HGAncientEurasia         HGSweeden  0.708309  0.705801 -0.440942
    10  HGAncientEurasia  HGAncientEurasia  0.725724  0.720317 -0.415958
    11           HGJapan           HGJapan  0.809975  0.808871 -0.121746
    12         HGSweeden         HGSweeden  0.761733  0.760462 -0.111567
    13            Africa         HGSweeden  0.623211  0.623209 -0.000584
    14            Africa        Neandertal  0.551766  0.551764 -0.000372

|Z-scores|>2 or 3 are often considered indicative of poor fit. Fit seems good, except that *HGJapan*, *HGAncientEurasia* and *HGSweeden* want to be closer to Neandertal. Why?

## Including Neandertal introgression into Eurasians
We know that present-day Eurasians carry 2-4%  Neandertal ancestry more than Africans. However, our Neandertal sample from Altai is 120,000 years-old, and predates our out-of-Africa. The Neandertal population that contributed genetic ancestry to AMHs in Eurasia had to be posterior to 120,000 ya. We will create a new ghost population (*GhostNeandertal*). How?

    import momi
    import logging
    logging.basicConfig(level=logging.INFO, filename="tutorial.log")
    model = momi.DemographicModel(N_e=1.2e4, gen_time=29, muts_per_gen=1.25e-8)
    sfs = momi.Sfs.load("sfs.gz")
    
    model.set_data(sfs)
    model.add_leaf("Africa", N="N_Africa", t=450)
    model.add_leaf("HGAncientEurasia", N="N_HGAncientEurasia", t=45000)
    model.add_leaf("HGJapan", N="N_Japan", t=3755)
    model.add_leaf("HGSweeden", N="N_HGSweeden", t=8895)
    model.add_leaf("Neandertal", N="N_Neandertal", t=120000)

    model.add_size_param("N_Africa")
    model.add_size_param("N_HGAncientEurasia")
    model.add_size_param("N_Japan")
    model.add_size_param("N_HGSweeden")
    model.add_size_param("N_Neandertal")
    
    model.move_lineages("HGSweeden", "HGJapan", t="t_HGJapan", N="Na_HGJapan") 
    model.move_lineages("HGJapan", "HGAncientEurasia", t="t_HGAncientEurasia", N="Na_HGAncientEurasia") 
    
    #THIS LINE INCLUDES GENE FLOW
    model.move_lineages("HGAncientEurasia", "GhostNeandertal", t="t_pulse", p="p_pulse") 
    
    model.move_lineages("HGAncientEurasia", "Africa", t="t_Africa", N="Na_Africa") 
    model.move_lineages("GhostNeandertal", "Neandertal", t="t_Neandertal", N="Na_Neandertal") 
    model.move_lineages("Africa", "Neandertal", t="t_root", N="Na_root")

    model.add_size_param("Na_HGJapan")
    model.add_size_param("Na_HGAncientEurasia")
    model.add_size_param("Na_Africa")
    model.add_size_param("Na_Neandertal")
    model.add_size_param("Na_root")
    
    model.add_pulse_param("p_pulse")
    model.add_time_param("t_HGJapan",lower=8895)
    model.add_time_param("t_HGAncientEurasia",lower=45000, lower_constraints=["t_HGJapan"])
    model.add_time_param("t_Africa",lower=2000,lower_constraints=["t_HGAncientEurasia"])
    model.add_time_param("t_pulse", upper_constraints=["t_Africa"],lower_constraints=["t_HGAncientEurasia"])
    model.add_time_param("t_Neandertal",lower=120000,lower_constraints=["t_pulse"])
    model.add_time_param("t_root", lower_constraints=["t_Neandertal","t_Africa"])
    model.optimize(method="TNC")
which returns:

       message: Max. number of function evaluations reached
            success: False
             status: 3
                fun: 0.017635994325904895
                  x: [ 1.763e+01  1.853e+01 ...  6.570e+03  2.761e+05]
                nit: 23
                jac: [-5.014e-07  0.000e+00 ... -4.060e-06 -1.851e-08]
               nfev: 171
         parameters:            N_Africa: 45405854.62299604
                      N_HGAncientEurasia: 111338819.40762144
                                 N_Japan: 931.4199230462771
                             N_HGSweeden: 3252.7269767444545
                            N_Neandertal: 65.84877251696435
                              Na_HGJapan: 4991.976326537487
                     Na_HGAncientEurasia: 2235.2965086541317
                               Na_Africa: 14454.597029029488
                           Na_Neandertal: 9346.947451485023
                                 Na_root: 16293.727657481453
                                 p_pulse: 0.1032362983236373
                               t_HGJapan: 25816.318894580996
                      t_HGAncientEurasia: 45000.0
                                t_Africa: 87251.15114901628
                                 t_pulse: 63620.467322138706
                            t_Neandertal: 126570.42641445273
                                  t_root: 402642.1607573377
      kl_divergence: 0.017635994325904895
     log_likelihood: -1559905.295077387

Once optimised, we can also plot this model:

    yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5]
    fig = momi.DemographyPlot(model, ["HGSweeden","HGJapan","HGAncientEurasia","Africa", "Neandertal","GhostNeandertal"], linthreshy=1e5, figsize=(10,12),major_yticks=yticks)
    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages("momi2.case1.pdf")
    pdf.savefig(fig.draw())
    pdf.close()


Is it the likelihood better? and significantly better? Why the introgression fraction into Eurasia is greater than 10%? Why the population size is so large in N_Africa? 

 
    zless -S 1.neo_impute_ph.filteredAA.bi.vcf.gz | grep -v "##" | awk '{ if ($10 ~ /(0\|1)|(1\|0)/){sum++} } END{print sum}' #Africa
    208981
    zless -S 1.neo_impute_ph.filteredAA.bi.vcf.gz | grep -v "##" | awk '{ if ($12 ~ /(0\|1)|(1\|0)/){sum++} } END{print sum} ' #Ancient Eurasia
    155461
    zless -S 1.neo_impute_ph.filteredAA.bi.vcf.gz | grep -v "##" | awk '{ if ($13 ~ /(0\|1)|(1\|0)/){sum++} } END{print sum} ' #HGSweeden
    135933
    zless -S 1.neo_impute_ph.filteredAA.bi.vcf.gz | grep -v "##" | awk '{ if ($14 ~ /(0\|1)|(1\|0)/){sum++} } END{print sum} ' # HGJapan
    112945

The Africa sample really is more heterozygous.

## Removing nucleotide transitions (TS)

The original VCF file contains nucleotide transitions (TS) and nucleotide transversions (TV) . This may be problematic as postmortem DNA damage, typical of ancient DNA molecules, is often manifested as an excess of nucleotide transitions (spontaneous Cytosine deamination leads to C>T errors). It is unlikely but still possible that African genomes retain some excess of DNA damage, as imputation works worse for ancestries under-represented in the reference panel. 

Try to remove nucleotide transitions (remember to index the resulting file), to see if it fixes the potential problem. For example:
    
    gunzip -c 1.neo_impute_ph.filteredAA.bi.vcf.gz | awk -F '\t' '( !(($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C")) || $0 ~ /^#/ )' | bcftools/bcftools view -m 2 -o 1.neo_impute_ph.filteredAA.bi.TV.vcf.gz -O z -

**TASK**: Format *1.neo_impute_ph.filteredAA.bi.TV.vcf.gz* into the momi2 input file (*sfs.TV.gz*) using the commands above. 

**QUESTION**: Why the mutation rate needs to be scaled after removing nucleotide transitions (see line 4: momi.DemographicModel)?

	import momi
	import logging
	logging.basicConfig(level=logging.INFO, filename="tutorial.log")
	model = momi.DemographicModel(N_e=1.2e4, gen_time=29, muts_per_gen=0.00000000390625)
    sfs = momi.Sfs.load("sfs.TV.gz")
    model.set_data(sfs)
    model.add_leaf("Africa", N="N_Africa", t=450)
    model.add_leaf("HGAncientEurasia", N="N_HGAncientEurasia", t=45000)
    model.add_leaf("HGJapan", N="N_Japan", t=3755)
    model.add_leaf("HGSweeden", N="N_HGSweeden", t=8895)
    model.add_leaf("Neandertal", N="N_Neandertal", t=120000)

    model.add_size_param("N_Africa")
    model.add_size_param("N_HGAncientEurasia")
    model.add_size_param("N_Japan")
    model.add_size_param("N_HGSweeden")
    model.add_size_param("N_Neandertal")
    
    model.move_lineages("HGSweeden", "HGJapan", t="t_HGJapan", N="Na_HGJapan") 
    model.move_lineages("HGJapan", "HGAncientEurasia", t="t_HGAncientEurasia", N="Na_HGAncientEurasia") 
    
    #THIS LINE INCLUDES GENE FLOW
    model.move_lineages("HGAncientEurasia", "GhostNeandertal", t="t_pulse", p="p_pulse") 
    
    model.move_lineages("HGAncientEurasia", "Africa", t="t_Africa", N="Na_Africa") 
    model.move_lineages("GhostNeandertal", "Neandertal", t="t_Neandertal", N="Na_Neandertal") 
    model.move_lineages("Africa", "Neandertal", t="t_root", N="Na_root")

    model.add_size_param("Na_HGJapan")
    model.add_size_param("Na_HGAncientEurasia")
    model.add_size_param("Na_Africa")
    model.add_size_param("Na_Neandertal")
    model.add_size_param("Na_root")
    
    model.add_pulse_param("p_pulse")
    model.add_time_param("t_HGJapan",lower=8895)
    model.add_time_param("t_HGAncientEurasia",lower=45000, lower_constraints=["t_HGJapan"])
    model.add_time_param("t_Africa",lower=2000,lower_constraints=["t_HGAncientEurasia"])
    model.add_time_param("t_pulse", upper_constraints=["t_Africa"],lower_constraints=["t_HGAncientEurasia"])
    model.add_time_param("t_Neandertal",lower=120000,lower_constraints=["t_pulse"])
    model.add_time_param("t_root", lower_constraints=["t_Neandertal","t_Africa"])
    model.optimize(method="TNC")

    yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5]
    fig = momi.DemographyPlot(model, ["HGSweeden","HGJapan","HGAncientEurasia","Africa", "Neandertal","GhostNeandertal"], linthreshy=1e5, figsize=(10,12),major_yticks=yticks)
    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages("momi2.case1.pdf")
    pdf.savefig(fig.draw())
    pdf.close()

Does it solve the issue if large *N_Africa*? One useful QC is to calculate the TS:TV ratio, which in humans is known to be about 2.1-2.2. This could have been done using the following command:

`bcftools stats 1.neo_impute_ph.filteredAA.bi.vcf.gz | grep TSTV` 

which returns a value of 2.21. We assume that the original VCF file was correct.

The fit of this last model is again reasonable, but the estimates make little sense. Be critical. 

	fit = momi.SfsModelFitStats(model)
	
	fit.f4("HGSweeden", "Africa", "Neandertal", None)
	JackknifeGoodnessFitStat(expected=0.007575222600482934, observed=0.00455537508882347, bias=-3.581825600430799e-07, sd=0.0014467483471938651, z_score=-2.0873343436104874)
	
	fit.f4("HGAncientEurasia", "Africa", "Neandertal", None)
	JackknifeGoodnessFitStat(expected=0.007575222600482913, observed=0.0074002639325956715, bias=1.6046377471518047e-06, sd=0.002886309956413501, z_score=-0.06061672880921051)

**QUESTION**: What do these f4 values really tell us?

## Archaic introgression into Iron Age Africans?
Another unexplored possibility is that Iron Age Africans carried ancestry from a yet-uncharacterized archaic hominin, which could explain also its increased heterozygosity, and unrealistic estimate of effective population sizes. 

    import momi
    import logging
    logging.basicConfig(level=logging.INFO, filename="tutorial.log")
    model = momi.DemographicModel(N_e=1.2e4, gen_time=29, muts_per_gen=1.25e-8)
    sfs = momi.Sfs.load("sfs.gz")
    
    model.set_data(sfs)
    model.add_leaf("Africa", N="N_Africa", t=450)
    model.add_leaf("HGAncientEurasia", N="N_HGAncientEurasia", t=45000)
    model.add_leaf("HGJapan", N="N_Japan", t=3755)
    model.add_leaf("HGSweeden", N="N_HGSweeden", t=8895)
    model.add_leaf("Neandertal", N="N_Neandertal", t=120000)

    model.add_size_param("N_Africa")
    model.add_size_param("N_HGAncientEurasia")
    model.add_size_param("N_Japan")
    model.add_size_param("N_HGSweeden")
    model.add_size_param("N_Neandertal")
    
    #THIS LINE INCLUDES GENE FLOW #2
    model.move_lineages("Africa", "GhostNeandertal2", t="t_pulse2", p="p_pulse2") 
        
    model.move_lineages("HGSweeden", "HGJapan", t="t_HGJapan", N="Na_HGJapan") 
    model.move_lineages("HGJapan", "HGAncientEurasia", t="t_HGAncientEurasia", N="Na_HGAncientEurasia") 
    
    #THIS LINE INCLUDES GENE FLOW #1
    model.move_lineages("HGAncientEurasia", "GhostNeandertal", t="t_pulse", p="p_pulse") 

    model.move_lineages("HGAncientEurasia", "Africa", t="t_Africa", N="Na_Africa") 
    model.move_lineages("GhostNeandertal", "Neandertal", t="t_Neandertal", N="Na_Neandertal") 
    model.move_lineages("Africa", "GhostNeandertal2", t="t_Neandertal2", N="Na_Neandertal2")
	model.move_lineages("GhostNeandertal2","Neandertal", t="t_root", N="Na_root")
	 
    model.add_size_param("Na_HGJapan")
    model.add_size_param("Na_HGAncientEurasia")
    model.add_size_param("Na_Africa")
    model.add_size_param("Na_Neandertal")
    model.add_size_param("Na_Neandertal2")
    model.add_size_param("Na_root")
    
    model.add_pulse_param("p_pulse")
    model.add_pulse_param("p_pulse2")
    
    model.add_time_param("t_HGJapan",lower=8895)
    model.add_time_param("t_HGAncientEurasia",lower=45000, lower_constraints=["t_HGJapan"])
    model.add_time_param("t_Africa",lower=2000,lower_constraints=["t_HGAncientEurasia"])
    model.add_time_param("t_pulse", upper_constraints=["t_Africa"],lower_constraints=["t_HGAncientEurasia"])
    model.add_time_param("t_pulse2", upper_constraints=["t_Africa"])
    model.add_time_param("t_Neandertal",lower=120000,lower_constraints=["t_pulse"])
    model.add_time_param("t_Neandertal2",lower_constraints=["t_pulse2","t_Africa"])
    model.add_time_param("t_root", lower_constraints=["t_Neandertal2","t_Neandertal"])
    
    model.optimize(method="TNC")
	    message: Max. number of function evaluations reached
        success: False
         status: 3
            fun: 0.011410006705917162
              x: [ 1.027e+01  1.944e+00 ...  3.186e+05  3.648e+04]
            nit: 24
            jac: [ 1.048e-05 -5.340e-06 ...  4.489e-10  4.233e-10]
           nfev: 210
     parameters:            N_Africa: 28933.737559174813
                  N_HGAncientEurasia: 6.983262096944417
                             N_Japan: 251.9379661276343
                         N_HGSweeden: 95.71592985046232
                        N_Neandertal: 20.257692376459225
                          Na_HGJapan: 18561.752946439115
                 Na_HGAncientEurasia: 894.4945157210198
                           Na_Africa: 16024.473877937186
                       Na_Neandertal: 7311.428244731894
                      Na_Neandertal2: 8217.035798157078
                             Na_root: 16810.864728741944
                             p_pulse: 0.057220616673059686
                            p_pulse2: 0.16006057168517643
                           t_HGJapan: 9933.288008539948
                  t_HGAncientEurasia: 45012.11109635486
                            t_Africa: 58302.11899863937
                             t_pulse: 45012.1371474509
                            t_pulse2: 41225.23299108599
                        t_Neandertal: 122269.03267414535
                       t_Neandertal2: 376869.05149820144
                              t_root: 413353.6422748137
	 kl_divergence: 0.011410006705917162
	 log_likelihood: -1557515.7983847614
 
    yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5]
    fig = momi.DemographyPlot(model, ["HGSweeden","HGJapan","HGAncientEurasia","Africa", "GhostNeandertal2","Neandertal","GhostNeandertal"], linthreshy=1e5, figsize=(10,12),major_yticks=yticks)
    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages("momi2.case2.pdf")
    pdf.savefig(fig.draw())
    pdf.close()

**IMPORTANT**: If you run this population model several times, you may find different results with almost the same likelihood value. Depending on how the two ghost populations are modelled, the parameters change for the rest of the tree.
