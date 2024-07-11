  seed =  -1

       seqfile = bears_bpp.txt
      Imapfile = bears_map.txt
       outfile = bears_out.txt
      mcmcfile = bears_mcmc.txt
      *Threads = 8 1 1

 *speciesdelimitation = 0        # fixed species tree
*speciesdelimitation = 1 0 2    # species delimitation rjMCMC algorithm0 and finetunee
speciesdelimitation = 1 1 2 1  # species delimitation rjMCMC algorithm1 finetune a m
        speciestree = 0        # species-tree by NNI pSlider ExpandRatio ShrinkRatio
*  speciesmodelprior = 1        # 0: uniform LH 1:uniform rooted trees 2: uniformSLH 3: uniformSRooted

speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 10  thibetanus americanus deningeri spelaeusA ingressus spelaeusB arctosA arctosB arctosC maritimus
                    6 4 1 8 1 1 37 7 5 19
                 ((thibetanus, americanus), ((deningeri, (spelaeusA, (ingressus, spelaeusB))), (arctosA, (arctosB, (arctosC, maritimus)))));
        phase =   0  0  0  0  0  0  0  0  0  0

       usedata = 1            # 0: no data prior 1:seq like
         nloci = 12          # number of data sets in seqfile

     cleandata = 0            # remove sites with ambiguity data 1:yes 0:no?

    *thetaprior = 3 0.002    # invgammaa b for theta integrated out by default add E or e to also sample theta
      *tauprior = 3 0.004      # invgammaa b for root tau & Dirichleta for other tau's
      thetaprior = 2 0.002   # invgamma(a, b) for theta
      tauprior = 2 0.01    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

*     phiprior = 0          # Betaa b for root tau & Dirichleta for other tau's

*     heredity = 0
*    locusrate = 0

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 8000
      sampfreq = 2
       nsample = 100000

