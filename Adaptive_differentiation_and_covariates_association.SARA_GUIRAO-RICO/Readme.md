# Adaptive differentiation and covariates association

Instructor: Sara Guirao Rico

In this practical class we will try to identify human **SNPs with evidences of local adaptive diversification** and/or **significantly associated to geographical variables** while taking into account the null correlation of allele frequencies across populations.

* We will use a subset of the data used in Coop et al. (2010). It consists of **genotypes** corresponding to **2335 SNPs** (2333 + 2) for **927 individuals** from **52 human populations** of the Human Genome Diversity Project (HGDP) panel (Conrad et al. 2006). The last two SNPs (rs12913832 and rs1042602) were added manually. 
These **two SNPs** have been previously reported to be **under positive selection in European populations** and are located in genes (ERC2 and TYR) involved in **light skin and eye color** (Wilde et al.2014).
* The covariates (geographical variables) explored are **latitude**, **longitude** and a categorical variable with value = 1 if the population is **European** and -1 if is not.

> [!Warning]
> Note that even we are performing the analysis with a very reduce dataset, the time it takes to run some of the models is above the available time in this practice and therefore, in some cases, we will use some of the results previously obtained (**precomputed results**) to calculate some of the statistics or to plot them using R.

## Working space
We are going to work in a container of the image *ppgcourseub/ppg24:latest*

You can find the executable binary of `BayPass` (g_baypass) in the folder `./software/baypass_public/source. 
 

## Data

### 1. Retrieve the **input data** and the **scripts** for this session:  
 
```bash
ghget https://github.com/ppgcourseUB/ppgcourse2024_week2/tree/main/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO
```
> [!Tip]
> Remember that Mac users may need first to download the folder in the "/tmp" directory and then move it to the ppguser 

* The folder has three main subfolders:

	*  input: genotype and covariate input data needed for the analysis. 
	*  scripts: scripts to execute some of the BayPass models and to perform simulations (PODs).
 	*  precomputed_results: folder with three different subfolders with precomputed results.

* The files in the precomputed_results folder are classified in subfolders according to the model/process (e.g., CORE, STDis,...) and for each, in turn, in simulations or plot subfolders.

### 2. The "baypass_utils.R" file 

This file contains six R functions (R Core Team, 2015) (simulate.baypass(); plot.omega(); fmd.dist(); geno2YN() and simulate.PC.correlated.covariate() that may helpful to interpret some of the results obtained with BAYPASS.  

Now you are ready to run BayPass.

## The CORE Model
The core model allows to perform **genome scans for differentiation** (covariate free) using the **XtX statistics** (\~Fst).
The main advantage of this approach is that explicitly **accounts for the covariance structure in population allele frequencies** (via estimating) resulting from the demographic history of the populations.

To run this model (using allele count data) you will need:
* The **number of populations** in the analysis (```-npop```)
* The **genotype file** (hgdp.geno): the genotypes (allele counts) for each SNP and population. In rows, the SNPs ordered according to their physical position on the chromosomes (if possible), except for the last two SNPs that where "artificially" introduced. In columns: populations. Each population has two columns: one for the number of the reference and the other for the alternative allele counts (```-gfile```). 
* A random number for the **seed** (in case of needed; ```-seed```)
* A **prefix to name the output** (```-outprefix```)

> *  For more specifications see [BayPass manual](https://forgemia.inra.fr/mathieu.gautier/baypass_public/-/blob/master/manual/BayPass_manual.pdf?ref_type=heads) 

### 1. Run BayPass under the **CORE Model**
 
> [!Warning]
> The **CORE model** should be run with **three different seeds** by executing the script "run_core_model_seed1.sh", "run_core_model_seed2.sh" and "run_core_model_seed3.sh" but **here we are going to run it only with one** for a matter of time:

```bash
bash ./scripts/run_core_model_seed1.sh
#bash ./scripts/run_core_model_seed2.sh
#bash  ./scripts/run_core_model_seed3.sh
```
> * This is the code in "run_core_model_seed1.sh" script:

```
#!/bin/bash

# run BayPass (CORE Model) with seed1
mkdir results_core
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./input/hgdp.geno -seed 15263 -outprefix ./results_core/hgdp_core_s1
```

> * It takes ~ 9-15 mins each one
> * This will generate 7 files for each seed.

### 2. **Sanity Check**. 

2.1. **In an R session**, we can compare the omega matrices obtained under the CORE Model when using different seeds to check consistency in the estimation of parameters of the model. 

```R
setwd("./results_core")
source("../baypass_utils.R")
require(corrplot); require(ape); require(geigen); require(mvtnorm)

#Read omegas obtained from running the CORE Model with three different seeds
omega_s1=as.matrix(read.table(file="hgdp_core_s1_mat_omega.out", header=F))
omega_s2=as.matrix(read.table(file="/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_core/hgdp_core_s2_mat_omega.out", header=F))
omega_s3=as.matrix(read.table(file="/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_core/hgdp_core_s3_mat_omega.out", header=F))

#Plot the comparison between omega-seed1 and omega-seed2:
pdf(file="omega_s1_s2_comparison.pdf")
	plot(omega_s1, omega_s2) ; abline(a=0,b=1)
dev.off()
```

```diff
- QUESTION: Are they similar?
```

2.2. Compute the **distances between pairs of omegas** to **check consistency*** in the estimation of parameters of the model.

```R
dist.12=fmd.dist(omega_s1, omega_s2)
dist.13=fmd.dist(omega_s1, omega_s3)
dist.23=fmd.dist(omega_s2, omega_s3)
dist.12
dist.13
dist.23
```

```diff
- QUESTION: Are they similar?
```

> * If the omegas are not significantly different we can assume that there is consistency in the parameters estimation and hence, you should choose one of the omegas to perform the subsequent analyses (e.g., omega 1).

### 3. **Visualization** of the shared history of populations (**after creating the pdf file, move it to your shared folder to be able to visualize the result in your laptop**).

3.1. Explore the **shared history of populations** by transforming the omega covariance matrix into a **correlation matrix** using the R function cov2cor().

```R
#Setup the population names for the omega matrix:
pop.names=c("Papuan","Melanesian","Surui","Pima","Maya","Karitiana","Columbian","Yi","Yaku","Xibo","Uygur","Tujia","Tu","She","Orogen","Naxi","Mongolia","Miao","Lahu","Japanese","Hezhen","Han","Daur",
"Dal","Cambodian","Sindhi","Pathan","Makrani","Kalash","Hazara","Burusho","Brahui","Balochi","Palestinian","Mozabite","Druze","Bedouin","Tuscan","Sardinian","Russian","Orcadian","French","Italian","Basque",
"Adygei","Yoruba","San","MbutiPygmy","Mandenka","BiakaPygmy","BantuSouthAfrica","BantuKenya")
dimnames(omega_s1)=list(pop.names,pop.names)

#Transform the omega covariance matrix into a correlation matrix
cor.mat_s1=cov2cor(omega_s1)

#Plot the omega correlation matrix:
pdf(file="correlation_matrix_core_model.pdf")
corrplot(cor.mat_s1, method="color",mar=c(2,1,2,2)+0.1, 
		main=expression("Correlation map based on"~hat(Omega)), 
		tl.cex=0.5)
dev.off()
```

3.2. Explore the **shared history of populations** by transforming the correlation matrix into a **hierarchical clustering tree** using the R function hclust().

```R
#Transform the correlation matrix into a hierarchical clustering tree
hgdp.tree_s1=as.phylo(hclust(as.dist(1-cor.mat_s1**2)))

#Plot the hierarchical clustering tree:
pdf(file="correlation_matrix_core_model_tree.pdf")
plot(hgdp.tree_s1, type="p", 
	main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"), 
	cex=0.5)
dev.off()
```

3.3. Explore the **shared history of populations** by performing a **heatmap and hierarchical clustering tree** (using the average agglomeration method). 

```R
pdf(file="omega_heatmap.pdf")
hclust.ave <- function(x) hclust(x, method="average")
heatmap(1-cor.mat_s1,hclustfun = hclust.ave,
	main=expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()
```

3.4. Explore the **shared history of populations** by performing an eigen-decomposition of the scaled covariance matrix of the population allele frequencies to allow representation in a two-dimension plot. This actually corresponds to a (between population) **PCA–like analysis** (**R in your laptop**).

```R
pdf(file="omega_PCA-like.pdf")
plot.omega(omega_s1,PC=c(1,2),pop.names=pop.names,
	main=expression("SVD of "*Omega),
	col=rainbow(nrow(omega_s1)), pch=20,pos=3)
dev.off()
```
> * If you want to modify the size of the text you should modify the baypass_utils.R script

### 4. Explore the values of the **XtXst statistic** (~Fst) obtained under the CORE Model (**after creating the pdf file, move it to your shared folder to be able to visualize the result in your laptop**)

```R
#Read the XtX file:
hgdp_s1.snp.res=read.table("hgdp_core_s1_summary_pi_xtx.out",h=T)
head(hgdp_s1.snp.res)

#Plot the XtXst values:
pdf(file="XtXst_core_model.pdf")
plot(hgdp_s1.snp.res$XtXst, xlab="SNP", ylab="XtXst", main="XtXst Seed 1")
	points(x= 2334, y=hgdp_s1.snp.res[hgdp_s1.snp.res[,1] == 2334, ]$XtXst, col = "red", pch=20)
	points(x= 2335, y=hgdp_s1.snp.res[hgdp_s1.snp.res[,1] == 2335, ]$XtXst, col = "blue", pch=20)
dev.off()
```

```diff
- QUESTION: Which are the XtXst outliers? How many there are? Do we need to perform a test to know how many of them are significant?
```

### 5. Check behavior of the *P*-values** associated to the XtXst estimator (**after creating the pdf file, move it to your shared folder to be able to visualize the result in your laptop**).

```R
pdf("omega_XtXst_pvals_hist.pdf")
hist(10**(-1*hgdp_s1.snp.res$log10.1.pval.),freq=F,breaks=50, 
    main=expression('XtXst '*italic(P)*'-value distribution'), 
    xlab=expression(italic(P)*'-value'))
	abline(h=1)
dev.off()

#Plot the XtXst and their respective P-values
pdf("omega_XTXst_and_Pvals.pdf")
layout(matrix(1:2,2,1))
plot(hgdp_s1.snp.res$XtXst,xlab="SNP", ylab="XtXst", main="XtXst" )
	points(x=2334, y=hgdp_s1.snp.res[hgdp_s1.snp.res[,1]==2334, ]$XtXst, col="red", pch=20)
	points(x=2335, y=hgdp_s1.snp.res[hgdp_s1.snp.res[,1]==2335, ]$XtXst, col="blue", pch=20)
plot(hgdp_s1.snp.res$log10.1.pval.,ylab=expression('XtXst '*italic(P)*'-value (-log10 scale)'), 
    	xlab="SNP", main=expression('XtXst '*italic(P)*'-value'))
	points(x=2334, y=hgdp_s1.snp.res[hgdp_s1.snp.res[,1]==2334, ]$log10.1.pval., col="red", pch=20)
	points(x=2335, y=hgdp_s1.snp.res[hgdp_s1.snp.res[,1]==2335, ]$log10.1.pval., col="blue", pch=20)
	abline(h=3,lty=2) #correspods to a P-value theshold = 0.001
dev.off()
```

```diff
- QUESTION: Where are the two putative outliers (red and blue points)?
```

5.1. Inspect the *P*-values associated two the two putative outliers.

```R
#Inspect the P-values associated two the two putative outliers
hgdp_s1.snp.res[hgdp_s1.snp.res[,1] == 2334, ]$log10.1.pval.
hgdp_s1.snp.res[hgdp_s1.snp.res[,1] == 2335, ]$log10.1.pval.
```

```diff
- QUESTION: Which XtXst values are significant?
```

### 6. Correct by False Discovery Rate (FDR)** by transforming the *P*-values into *q*-values (**after creating the pdf file, move it to your shared folder to be able to visualize the result in your laptop**).

```R
library("qvalue")

#Get the XtXst P-values
pvals <- 10^(-hgdp_s1.snp.res$log10.1.pval)

#Create a qobj
qobj <- qvalue(p = pvals)

#Get the qvalues from P-values
qvals <- qobj$qvalues

#Count and compare how many SNPs are significant without and with FDR correction with alpha = 0.001
sum(pvals < 0.001)
sum(qvals < 0.001)
```
> * In case the *P*-values of the XtXst do not behave well, you will need to perform "neutral" simulations (Pseudo Observed Data –PODs–).

### 7. Pseudo Observed Data (PODs) 

Here, we are going to **simulate data (PODs)** using the R function simulate.baypass() in the baypass_utils.R script (provided in the BayPass package). PODs are simulated under the inference model (e.g., using posterior estimates of the covariance matrix and the a and b parameters of the beta prior distribution for the overall (across population) SNP allele frequencies).
Once these PODS are simulated, we need to **run again the CORE Model with the PODs as input** to built the \"expected\" distribution of the XtX values under the inference model in order to find which of the observed XtX values are significantly different from the expected (**calibration process**)

> We want to perform **two different sets of simulations** to inspect how many simulations are needed to retrieve the estimated demographic history: i) **simulating 1,000 PODs**; ii) simulating **100,000 PODs**.

> [!Warning]
> * However, here we are going to run only the first (simu.hgdp_1000) of the two simulation experiments for a matter of time. Instead, we will use the precomputed files with the 100,000 simulations.

7.1. **Simulate 1,000** Pseudo Observed Data (PODs):

In **R**

```R
#Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution obtained when running the CORE Model
pi.beta.coef=read.table("hgdp_core_s1_summary_beta_params.out",h=T)$Mean

#Upload the original data to obtain total allele count (sample size for each population). Do this by using the geno2YN() function in baypass_utils.R script
hgdp.data<-geno2YN("../input/hgdp.geno")

#Read the omega matrix from seed1 obtained from running the CORE Model:
omega_s1=as.matrix(read.table(file="hgdp_core_s1_mat_omega.out", header=F))

#Simulated 1000 PODs
simu.hgdp_1000 <- simulate.baypass(omega.mat=omega_s1, nsnp=1000, 
    sample.size=hgdp.data$NN, beta.pi=pi.beta.coef, pi.maf=0, suffix="hgdp_pods_1000")
q()
```

> * The G.hgdp_pods_1000 file is now the new genotype input file resulting from the simulation process.  

7.2. **Run again, the CORE Model** only with the first set of simulations (G.hgdp\_pods\_1000) using the script "run_core_1000_simulations.sh":

```bash
bash ./scripts/run_core_1000_simulations.sh 
```
> * This is the code to run the "run_core_1000_simulations.sh" script

```bash
#!/bin/bash
  
# run BayPass (CORE Model) with the 1,000 PODs as input
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./results_core/G.hgdp_pods_1000 -outprefix ./results_core/hgdp_pod_1000
```

> [!Warning]
> For the second set of simulations, we are going to **use the precomputed files** resulting from running the CORE Model with the 100,000 simulations as input.
  
7.3. **Sanity Check**.

Here, we are **comparing the simulated data (PODS)** under the inference model **to the observed data** to assess if the inference model (posterior distributions for the covariance matrix and the other hyperparameters) is giving us \"valid\" predictions about the \"reality\". In other words, if the model we have inferred is able to generate data similar results to those observed and in case of yes, how many simulations are needed.

**In R**

```R
setwd("./results_core")
source("../baypass_utils.R")

#Read the omega matrix from seed1 obtained from running the CORE Model: 
omega_s1=as.matrix(read.table(file="hgdp_core_s1_mat_omega.out", header=F))

#Get the omega estimated from the PODs:
pod.omega_1000=as.matrix(read.table("hgdp_pod_1000_mat_omega.out"))

#Plot the observed versus the simulated omegas:
pdf(file="comparison_omega_obs_sim_1000_pods.pdf")
plot(pod.omega_1000, omega_s1, xlab="Omega from PODs", 
    ylab="Omega from real data", main="1000 PODs") 
	abline(a=0,b=1)
dev.off()

#Compute the distance between the observed and simulated omegas
fmd.dist(pod.omega_1000, omega_s1)

#Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution obtained when running the CORE Model
pi.beta.coef=read.table("hgdp_core_s1_summary_beta_params.out",h=T)$Mean

#Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the PODs:
pod.pi.beta.coef_1000=read.table("hgdp_pod_1000_summary_beta_params.out",h=T)$Mean

#Plot the observed versus the simulated pi.beta:
pdf(file="comparison_pi.beta_obs_sim_1000_pods.pdf")
plot(pod.pi.beta.coef_1000, pi.beta.coef, 
	xlab="pi.beta from PODs", ylab="pi.beta from real data", main="1000 PODs") 
	abline(a=0,b=1)
dev.off()

### For 100000 PODS
#Get the omega estimated from the PODs:
pod.omega_100000=as.matrix(read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_core/hgdp_pod_100000_mat_omega.out"))

#Plot the observed versus the simulated omegas:
pdf(file="comparison_omega_obs_sim_100000_pods.pdf")
plot(pod.omega_100000, omega_s1, xlab="Omega from PODs", 
    ylab="Omega from real data", main="100000 PODs") 
	abline(a=0,b=1)
dev.off()

#Compute the distance between the observed and simulated omegas:
fmd.dist(pod.omega_100000, omega_s1)

#Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the PODs:
pod.pi.beta.coef_100000=read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_core/hgdp_pod_100000_summary_beta_params.out",h=T)$Mean

#Plot the observed versus the simulated pi.beta:
pdf(file="comparison_pi.beta_obs_sim_100000_pods.pdf")
plot(pod.pi.beta.coef_100000, pi.beta.coef, xlab="pi.beta from PODs", 
	ylab="pi.beta from real data", main="100000 PODs")
	abline(a=0,b=1)
dev.off()
```

```diff
- QUESTION: Are the omega matrix distances between the observed data and those from the 1,000 and 100,000 PODs similar?
```

```diff
- QUESTION: Where are the dots are falling in both sets of similations? What is the main difference when comparing the two simulation experiments (1,000 and 100,000 PODs) to the observed data?
```
    
## The STANDARD Model (STDis): importance sampling 
This model allows to evaluate to which extent the **population covariables are (linearly) associated to each marker/SNP**.
The estimation of the beta regression coefficients for each SNP and covariable is performed using the importance sampling approach.

* Remember that this model is **recommended when the number of populations is small (e.g.,  8**) and/or **when populations are highly differentiated**.
* The importance sampling algorithm relies on approximations to estimate the regression coefficients (do not sample them from the posterior distribution and hence, **the model should be run 3-5 times with different seeds** to check consistency across runs). 
* The **median of the parameter of insterest** (e.g., beta) among runs can be taken as the estimate (beta) of the parameter.
* Bayes factor importance sampling (BFis) and the approximated Bayesian *P*-value (eBPis) are used to evaluate if a particular SNP is associated with a particular covariable and, as the XtX statistic in the CORE Model, they should be calibrated. 
* Note that unlike the other models (AUX, AUX-LD, STDmcmc), this model calculates the covariance matrix (omega) and the correlation parameter (beta) at the same time.

To run this model (using allele count data) you will need:
* The **number of populations** in the analysis (```-npop```)
* The **genotype file** (hgdp.geno): the genotypes (allele counts) for each SNP and population. In rows, the SNPs ordered according to their physical position on the chromosomes (if possible), except for the last two SNPs that where "artificially" introduced. In columns: populations. Each population has two columns: one for the reference and the other for the alternative allele counts (```-gfile```). 
* The **covariates file** (covariates in the input folder): In rows, the covariates. In columns, populations (one column per population). The order of the populations should be the same as in the genotype file (```-efile```).
* To specify if you want to **scale** covariables (```-scalecov```)
* **A prefix to name the output** (```-outprefix```)

### 1. Run BayPass with the STANDARD model importance sampling** by submit the job script "run_stdis_model.sh" using **the sbatch command:

```bash
bash ./scripts/run_stdis_model.sh
```
> * This is the code to run the "run_stdis_model.sh" script

```bash
#!/bin/bash                                                                                                                         

# run BayPass (STDis Model)
mkdir results_standard
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./input/hgdp.geno -efile ./input/covariates -scalecov -outprefix ./results_standard/hgdp_stdis
```

> * This generates 9 output files.
> * It takes ~ 9-15 mins.

### 2. Inspect the obtained **results** (**after creating the pdf file, move it to your shared folder to be able to visualize the result in your laptop**).

```R
setwd("./results_standard")
source("../baypass_utils.R")

#Read the file with the BF, the eBPis and the correlation coefficients parameters
hgdp_stdis.snp.res=read.table("hgdp_stdis_summary_betai_reg.out",h=T)

#Plot BF, the eBPis and the correlation coefficients for the Latitude covariable
pdf("hgdp_stdis_model_lat.pdf")
layout(matrix(1:3,3,1))
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1,]$BF.dB.,
	xlab="SNP",ylab="BFis (in dB)")
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2334, ]$BF.dB., col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2335, ]$BF.dB., col="blue", pch=20)
	abline(h=10,lty=2, col="green") #strong evidence
	abline(h=15,lty=2, col="blue")  #very strong evidence
	abline(h=20,lty=2, col="red")   #decesive evidence
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1,]$eBPis,
	xlab="SNP",ylab="eBPis")
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2334, ]$eBPis, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2335, ]$eBPis, col="blue", pch=20)
	abline(h=3,lty=2, col="red")   #Pval = 0.001
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1] == 1,]$Beta_is,
	xlab="SNP",ylab=expression(beta~"coefficient"))
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2334, ]$Beta_is, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2335, ]$Beta_is, col="blue", pch=20)
	mtext("STDis MODEL: Latitude",side=3,line=- 2,outer=TRUE)
dev.off()

#Plot BF, the eBPis and the correlation coefficients for the Longitude covariable
pdf("hgdp_stdis_model_lon.pdf")
layout(matrix(1:3,3,1))
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2,]$BF.dB.,
	xlab="SNP",ylab="BFis (in dB)")
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2334, ]$BF.dB., col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2335, ]$BF.dB., col="blue", pch=20)
	abline(h=10,lty=2, col="green") #strong evidence
	abline(h=15,lty=2, col="blue")  #very strong evidence
	abline(h=20,lty=2, col="red")   #decesive evidence
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1] == 2,]$eBPis,
	xlab="SNP",ylab="eBPis")
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2334, ]$eBPis, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2335, ]$eBPis, col="blue", pch=20)
	abline(h=3,lty=2, col="red")   #Pval=0.001
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2,]$Beta_is,
	xlab="SNP",ylab=expression(beta~"coefficient"))
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2334, ]$Beta_is, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2335, ]$Beta_is, col="blue", pch=20)
	mtext("STDis MODEL: Longitude",side=3,line=- 2,outer=TRUE)
dev.off()

#Plot BF, the eBPis and the correlation coefficients for the European origin covariable
pdf("hgdp_stdis_model_eu.pdf")
layout(matrix(1:3,3,1))
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3,]$BF.dB.,
	xlab="SNP",ylab="BFis (in dB)")
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2334, ]$BF.dB., col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2335, ]$BF.dB., col="blue", pch=20)
	abline(h=10,lty=2, col="green") #strong evidence
	abline(h=15,lty=2, col="blue")  #very strong evidence
	abline(h=20,lty=2, col="red")   #decesive evidence
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3,]$eBPis,
	xlab="SNP",ylab="eBPis")
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2334, ]$eBPis, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2335, ]$eBPis, col="blue", pch=20)
	abline(h=3,lty=2, col="red")   #Pval=0.001
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3,]$Beta_is,
xlab="SNP",ylab=expression(beta~"coefficient"))
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2334, ]$Beta_is, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2335, ]$Beta_is, col="blue", pch=20)
	mtext("STDis MODEL: European origin",side=3,line=- 2,outer=TRUE)
dev.off()
```

```diff
- QUESTION: How many significant SNPs are correlating with any of the covariates? based on what criteria, BF or eBPis? Are all of them correlating in the same way?
```

### 3. Calibrate the STDis Parameters (BF, the eBPis and the correlation coefficients) using PODs.

3.1. **Simulate 10,000 PODs** using script "run_10000_simulations.sh".

**In R**

```R
setwd("./results_core")
source("../baypass_utils.R")

#Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution obtained when running the CORE Model
pi.beta.coef=read.table("hgdp_core_s1_summary_beta_params.out",h=T)$Mean

#Upload the original data to obtain total allele count (sample size for each population). Do this by using the geno2YN() function in baypass_utils.R script
hgdp.data<-geno2YN("../input/hgdp.geno")

#Read the omega matrix from seed1 obtained when running the CORE Model:
omega_s1=as.matrix(read.table(file="hgdp_core_s1_mat_omega.out", header=F))

#Generate 10,000 PODs
simu.hgdp_10000 <- simulate.baypass(omega.mat=omega_s1, nsnp=10000,
	sample.size= hgdp.data$NN, beta.pi=pi.beta.coef, pi.maf=0, suffix="hgdp_pods_10000")
```

3.2. Run STDis model with 10,000 PODS as input with the script "run_stdis_10000_simulations.sh" 

```bash
bash ./scripts/run_stdis_10000_simulations.sh
```
> * This is the code to run the "run_stdis_10000_simulations.sh" script:

```bash
#!/bin/bash                                                                                                             

# run BayPass (CORE Model) with the 10000 PODs as input
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./results_core/G.hgdp_pods_10000 -efile ./input/covariates -scalecov  -outprefix ./results_standard/hgdp_stdis_10000_pods
```

3.3. **Sanity check**.

```R
setwd("./results_standard")
source("../baypass_utils.R")

#Read the files resulting from running the CORE Model and from simulating PODs
omega_s1=as.matrix(read.table(file="/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/results_core/hgdp_core_s1_mat_omega.out", header=F))
pi.beta.coef=read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/results_core/hgdp_core_s1_summary_beta_params.out",h=T)$Mean

pod.omega_10000=as.matrix(read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_standard/hgdp_stdis_10000_pods_mat_omega.out"))
pod.pi.beta.coef_10000=read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_standard/hgdp_stdis_10000_pods_summary_beta_params.out",h=T)$Mean

#Plot the observed versus the simulated omegas:
pdf(file="comparison_omega_obs_sim_10000_pods.pdf")
plot(pod.omega_10000, omega_s1, xlab="Omega from PODs", 
    ylab="Omega from real data", main="10000 PODs") 
    abline(a=0,b=1)
dev.off()

#Compute the distance between the observed and simulated omegas
fmd.dist(pod.omega_10000, omega_s1)

#Plot the observed versus the simulated pi.beta:
pdf(file="comparison_pi.beta_obs_sim_10000pods.pdf")
plot(pod.pi.beta.coef_10000, pi.beta.coef, xlab="pi.beta from PODs", 
    ylab="pi.beta from Data", main="10000 PODs") 
    abline(a=0,b=1)
dev.off()
```

3.4. **Calibrate** the BF, the eBPis and the correlation coefficients parameters.

```R
#Read the file with the BF, the eBPis and the correlation coefficients parameters
hgdp_stdis.snp.res=read.table("hgdp_stdis_summary_betai_reg.out",h=T)

#Read the output file with the BF, eBPis and Beta correlation coefficients calculated from pods
hgdp_stdis_10000_pods_param=read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_standard/hgdp_stdis_10000_pods_summary_betai_reg.out",h=T)
 
### Calibrate the STDis parameters for Latitude
#Compute the 1% threshold of BF
hgdp_stdis_10000_pods_BF_thresh_lat=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==1,]$BF.dB.,probs=0.99)
#Compute the 1% threshold of eBPis
hgdp_stdis_10000_pods_eBPis_thresh_lat=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==1,]$eBPis,probs=0.99)
#Compute the 1% right threshold of Beta_is
hgdp_stdis_10000_pods_Beta_is_thresh1_lat=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==1,]$Beta_is,probs=0.99)
#Compute the 1% left threshold of Beta_is
hgdp_stdis_10000_pods_Beta_is_thresh2_lat=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==1,]$Beta_is,probs=0.01)

#Plot the observed parameters with the new empiral thresholds obtained from the PODs
pdf("stdis_model_calibrated_lat.pdf")
layout(matrix(1:3,3,1))
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1,]$BF.dB.,
    ylab="Calibrated BF", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_BF_thresh_lat,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2334, ]$BF.dB., col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2335, ]$BF.dB., col="blue", pch=20)

plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1,]$eBPis,
    ylab="Calibrated eBPis", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_eBPis_thresh_lat,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2334, ]$eBPis, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2335, ]$eBPis, col="blue", pch=20)

plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1,]$Beta_is,
    ylab="Calibrated Beta", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_Beta_is_thresh1_lat,lty=2)
	abline(h=hgdp_stdis_10000_pods_Beta_is_thresh2_lat,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2334, ]$Beta_is, col = "red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==1 & hgdp_stdis.snp.res[,2]==2335, ]$Beta_is, col="blue", pch=20)
	mtext("STDis Model: Latitude",side=3,line=-1.5,outer=TRUE)
dev.off()

### Calibrate the STDis parameters for Longitude
#Compute the 1% threshold of BF
hgdp_stdis_10000_pods_BF_thresh_lon=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==2,]$BF.dB.,probs=0.99)
#Compute the 1% threshold of eBPis
hgdp_stdis_10000_pods_eBPis_thresh_lon=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==2,]$eBPis,probs=0.99)
#Compute the 1% right threshold of Beta_is
hgdp_stdis_10000_pods_Beta_is_thresh1_lon=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==2,]$Beta_is,probs=0.99)
#Compute the 1% left threshold of Beta_is
hgdp_stdis_10000_pods_Beta_is_thresh2_lon=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==2,]$Beta_is,probs=0.01)

#Plot the observed parameters with the new empiral thresholds obtained from the PODs
pdf("stdis_model_calibrated_lon.pdf")
layout(matrix(1:3,3,1))
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2,]$BF.dB.,
    ylab="Calibrated BF", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_BF_thresh_lon,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2334, ]$BF.dB., col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2335, ]$BF.dB., col="blue", pch=20)

plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2,]$eBPis,
    ylab="Calibrated eBPis", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_eBPis_thresh_lon,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2334, ]$eBPis, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2335, ]$eBPis, col="blue", pch=20)

plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2,]$Beta_is,
    ylab="Calibrated Beta", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_Beta_is_thresh1_lon,lty=2)
	abline(h=hgdp_stdis_10000_pods_Beta_is_thresh2_lon,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2334, ]$Beta_is, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==2 & hgdp_stdis.snp.res[,2]==2335, ]$Beta_is, col="blue", pch=20)
	mtext("STDis Model: Longitude",side=3,line=- 1.5,outer=TRUE)
dev.off()

### Calibrate the STDis parameters for European origin
#Compute the 1% threshold of BF
hgdp_stdis_10000_pods_BF_thresh_eu=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==3,]$BF.dB.,probs=0.99)
#Compute the 1% threshold of eBPis
hgdp_stdis_10000_pods_eBPis_thresh_eu=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==3,]$eBPis,probs=0.99)
#Compute the 1% right threshold of Beta_is
hgdp_stdis_10000_pods_Beta_is_thresh1_eu=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==3,]$Beta_is,probs=0.99)
#Compute the 1% left threshold of Beta_is
hgdp_stdis_10000_pods_Beta_is_thresh2_eu=quantile(hgdp_stdis_10000_pods_param[hgdp_stdis_10000_pods_param[, 1]==3,]$Beta_is,probs=0.01)

#Plot the observed parameters with the new empiral thresholds obtained from the PODs
pdf("stdis_model_calibrated_eu.pdf")
layout(matrix(1:3,3,1))
plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3,]$BF.dB.,
    ylab="Calibrated BF", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_BF_thresh_eu,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2334, ]$BF.dB., col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2335, ]$BF.dB., col="blue", pch=20)

plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3,]$eBPis,
    ylab="Calibrated eBPis", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_eBPis_thresh_eu,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2334, ]$eBPis, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2335, ]$eBPis, col="blue", pch=20)

plot(hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3,]$Beta_is,
    ylab="Calibrated Beta", xlab="SNP")
	abline(h=hgdp_stdis_10000_pods_Beta_is_thresh1_eu,lty=2)
	abline(h=hgdp_stdis_10000_pods_Beta_is_thresh2_eu,lty=2)
	points(x=2334, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2334, ]$Beta_is, col="red", pch=20)
	points(x=2335, y=hgdp_stdis.snp.res[hgdp_stdis.snp.res[,1]==3 & hgdp_stdis.snp.res[,2]==2335, ]$Beta_is, col="blue", pch=20)
	mtext("STDis Model: European Origin",side=3,line=- 1.5,outer=TRUE)
dev.off()
```

```diff
- QUESTION: How many significant SNPs are correlating with any of the covariates? based on what criteria, BF or eBPis? Are all of them correlating in the same way?
```

# EXTRA EXERCISES

## The STDis Model and CONTRAST Analysis
This **"combined" analysis** allows to evaluate to which extent the **population binary covariables are associated to specific markers/SNPs** using two different approaches: i) running the **STDis Model** and hence, testing if covariates are (linearly) associated to each marker/SNP or ii) performing a **contrast analysis** to compute contrast of standardized allele frequencies between two groups of populations.

* The contrast anaysis can be interesting if the assumed linear relationship between the covariates and allele frequencies is not entirely satisfactory. 

* The association analyses (STDis, STDmcmc and AUX Models) carried out with categorical or binary covariables can be problematic when dealing with small data sets or if one wishes to disregard some populations. 

* Remember all the recomendations described above for running the STDis Model

To run this analysis (using allele count data) you will need:
* The **number of populations** in the analysis (```-npop```)
* The **genotype file** (hgdp.geno): the genotypes (allele counts) for each SNP and population. In rows, the SNPs ordered according to their physical position on the chromosomes (if possible), except for the last two SNPs that where "artificially" introduced. In columns: populations. Each population has two columns: one for the reference and the other for the alternative allele counts (```-gfile```). 
* The **binary covariates file** (covariates_eu in the input folder): In rows, the covariates. In columns, populations (one column per population). The order of the populations should be the same as in the genotype file (```-efile```).
* The **binary covariates file** (covariates_eu in the input folder): In rows, the binary covariates. In columns, the group membership of each population, 1 for first group, -1 for the alternative group, and possibly 0 if excluded from the contrast computation (```-contrastfile```).
* To specify if you want to **scale covariables** (```-scalecov```)
* A **prefix to name the output** (```-outprefix```)
 
1. Run the **STDis and the Contrast Model** using the script "run_stdis_contr_model.sh". 

```bash
bash  ./scripts/run_stdis_contr_model.sh
```

> * This is the code to run the "run_stdis_contr_model.sh" script

```bash
#!/bin/bash                                                                                                              
# run BayPass (STDis and Contrast Model)
mkdir results_stdis_contr
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./input/hgdp.geno -contrastfile ./input/covariates_eu -efile ./input/covariates_eu -outprefix ./results_stdis_contr/hgdp_contrast
````
> * It generates 9 output files.
> * It takes about 6-15 mins.

2. **Inspect** the obtained **results**.

```R
setwd("./results_stdis_contr")

#Read the files with the BF and the C2
covariates_eu.bf=read.table("hgdp_contrast_summary_betai_reg.out",h=T)$BF.dB.
covariates_eu.C2=read.table("hgdp_contrast_summary_contrast.out",h=T)

#Check the behavior of the P-values associated to the C2
pdf("C2_pvals_hist.pdf")
hist(10**(-1*covariates_eu.C2$log10.1.pval.),freq=F,breaks=50,
    main=expression('C2 '*italic(P)*'-value distribution'), 
    xlab=expression(italic(P)*'-value'))
	abline(h=1)
dev.off()
```

```diff
- QUESTION: Are the C2 *P*-values behaving well?
```

3. **Calibrate** C2 and BF statistics.

3.1. **Simulate 10,000 neutral PODs** using the script "run_10000_c2_simulations.sh". 

**In R**

```R
setwd("./results_stdis_contr")
source("../baypass_utils.R")

#Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
c2.pi.beta.coef=read.table("hgdp_contrast_summary_beta_params.out",h=T)$Mean

#Upload the original data to obtain total allele count
hgdp.data<-geno2YN("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/input/hgdp.geno")

#Read the omega matrix:
omega_contrast=as.matrix(read.table(file="hgdp_contrast_mat_omega.out", header=F))

#Generate 10000 PODs
simu.C2.10000 <- simulate.baypass(omega.mat=omega_contrast,nsnp=10000,
	sample.size=hgdp.data$NN, beta.pi=c2.pi.beta.coef, pi.maf=0, 
	suffix="hgdp_C2_10000_pods")
```

3.2. Run the **STDis and contrast Models with the 10,000 PODs as input** using the script "run_stdis_contrast_10000_simulations.sh". 

```bash
bash ./scripts/run_stdis_contrast_10000_simulations.sh
```

> * This is the code to run the "run_stdis_contrast_10000_simulations.sh" script:

```bash
#!/bin/bash                                                                                                             

# run BayPass (CORE Model) with the 10000 c2 PODs as input
../software/baypass_public/sources/g_baypass -npop 52 -gfile ./results_stdis_contr/G.hgdp_C2_10000_pods -contrastfile ./input/covariates_eu -efile ./input/covariates_eu -outprefix ./results_stdis_contr/hgdp_contrast_10000_pods 
```

> * It generates 9 output files
> * It takes about 20-40 mins

3.3. **Sanity Check**.

```R
setwd("./results_stdis_contr")
source("../baypass_utils.R")

#Read the omega matrix from observed data:
omega_contrast=as.matrix(read.table(file="hgdp_contrast_mat_omega.out", header=F))

#Get estimate of omega from the PODs
pod.c2.omega=as.matrix(read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_stdis_contr/hgdp_contrast_10000_pods_mat_omega.out"))

#Plot the observed versus the simulated omegas:
pdf(file="comparison_omega_obs_sim_10000_pods.pdf")
plot(pod.c2.omega,omega_contrast) 
    abline(a=0,b=1)
dev.off()
    
#Get the distance between the simulated and the real omega    
fmd.dist(pod.c2.omega,omega_contrast)

#Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the observed data
c2.pi.beta.coef=read.table("hgdp_contrast_summary_beta_params.out",h=T)$Mean

#Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the POD analysis
pod.c2.pi.beta.coef=read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_stdis_contr/hgdp_contrast_10000_pods_summary_beta_params.out",h=T)$Mean

#Plot the observed versus the simulated pi.beta:
pdf(file="comparison_pi.beta_obs_sim_10000_pods.pdf")
plot(pod.c2.pi.beta.coef,c2.pi.beta.coef) 
    abline(a=0,b=1)
dev.off()
```

3.4. **Calibrate** the C2 and BF.

```R
setwd("./results_stdis_contr")
source("../baypass_utils.R")

#Read the obtained results from real data
covariates_eu.bf.all <-
	read.table(file="hgdp_contrast_summary_betai_reg.out", h=T)

#Read the files with the BF and the C2
covariates_eu.bf=read.table("hgdp_contrast_summary_betai_reg.out",h=T)$BF.dB.
covariates_eu.C2=read.table("hgdp_contrast_summary_contrast.out",h=T)

#Read the files with the simulated C2 and BF
pod.c2_10000=read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_stdis_contr/hgdp_contrast_10000_pods_summary_contrast.out",h=T)
pod.BF.10000=read.table("/home/ppguser/Adaptive_differentiation_and_covariates_association.SARA_GUIRAO-RICO/precomputed_results/results_stdis_contr/hgdp_contrast_10000_pods_summary_betai_reg.out",h=T)

#Compute the 1% threshold of BF
pod.c2_10000_thresh=quantile(pod.c2_10000$M_C2,probs=0.99)
pod.BF_10000_thresh=quantile(pod.BF.10000$BF.dB.,probs=0.99)

#Plot the observed C2 and BF with the new thresholds obtained from PODs
pdf("Calibrated_C2_and_BF.pdf")
plot(covariates_eu.bf,covariates_eu.C2$M_C2,
	xlab="calibrated BF",ylab="calibrated C2")
	abline(h=pod.c2_10000_thresh,lty=2) #
	abline(v=pod.BF_10000_thresh,lty=2) #
	points(x=covariates_eu.bf.all[covariates_eu.bf.all[,2]==2334, ]$BF.dB., y=covariates_eu.C2[covariates_eu.C2[ ,2 ]==2334, ]$M_C2, col="red", pch=20)
	points(x=covariates_eu.bf.all[covariates_eu.bf.all[,2]==2335, ]$BF.dB., y=covariates_eu.C2[covariates_eu.C2[ ,2 ]==2335, ]$M_C2, col="blue", pch=20)
dev.off()
```

3.5. **Plot** the observed C2 and BF for a matter of comparison (**R in your laptop**).

```R
setwd("./results_stdis_contr")
source("../baypass_utils.R")

#Read the obtained results for the real data
covariates_eu.bf.all <-
	read.table(file="hgdp_contrast_summary_betai_reg.out", h=T)

#Plot
pdf("C2_and_BF_and_Pvals.pdf")
plot(covariates_eu.bf,covariates_eu.C2$log10.1.pval.,
	xlab="BF",ylab="C2 p-value (-log10 scale)")
	abline(h=3,lty=2) #0.001 p--value theshold
	abline(v=10,lty=2) #BF threshold for strong evidence (according to Jeffreys’ rule)
	points(x=covariates_eu.bf.all[covariates_eu.bf.all[,2]==2334, ]$BF.dB., 
		y=covariates_eu.C2[covariates_eu.C2[ ,2 ]==2334, ]$log10.1.pval., col="red", pch=20)
	points(x=covariates_eu.bf.all[covariates_eu.bf.all[,2]==2335, ]$BF.dB., 
		y=covariates_eu.C2[covariates_eu.C2[ ,2 ]==2335, ]$log10.1.pval., col="blue", pch=20)
dev.off()
```

```diff
- QUESTION: How many SNPs are significant before and after calibrating the C2 statistic using PODs? and how many after calibrating the BF?
```

### BIBLIOGRAPHY

* Gautier M. 2015. Genome-Wide Scan for Adaptive Divergence and Association with Population-Specific Covariates. Genetics 201(4): 1555–1579. 
* Günther T and Coop G. 2013. Robust Identification of Local Adaptation from Allele Frequencies. Genetics 195(1): 205–220.
* Coop G. Witonsky D, Di Rienzo A, Pritchard JK. 2010. Using environmental correlations to identify loci underlying local adaptation. Genetics 185: 1411–1423.
* Hoban S, Kelley JL, Lotterhos KE, Antolin MF, Bradburd G, Lowry DB, Poss ML, Reed LK, Storfer A, Whitlock MC. 2016. Finding the Genomic Basis of Local Adaptation: Pitfalls, Practical Solutions, and Future Directions. Am Nat. 188(4): 379–397


