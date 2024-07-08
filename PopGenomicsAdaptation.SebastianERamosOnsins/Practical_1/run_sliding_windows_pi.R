setwd("./")
#setwd("/Users/sramos/Documents/CV/Classes/PopGenomics_Adaptation_Course-2024/2024/Practical_1/")

#install vcfR library
if (!require('vcfR')) install.packages('vcfR'); library('vcfR')

#read files from slim output
slim_files <- system("ls *.practical1.output.VCF",intern=T)

pdf("SelectiveSweep_SlidingWindows_Pi.pdf")
for(model in slim_files) {
  #read vcf file
  vcf <- read.vcfR(file=model)
  vcfg <- vcfR2genind(vcf)
  vcfgs <- summary(vcfg)
  
  #take the positions of SNPs and their heterozygosity
  pos <- as.numeric(vcf@fix[,2])
  #pi.o <- vcfgs$Hobs
  #pi.e <- vcfgs$Hexp
  pi.o <- NULL
  #pi.e <- NULL
  for(l1 in 1:length(vcf@gt[,1])) {
    r1 <- sapply(vcf@gt[l1,-1],strsplit,split="|",fixed=T)
    Ho <- 0; #He <- 0
    for(r in 1:length(r1)) {
      Ho <- Ho + as.numeric(r1[[r]][1]==r1[[r]][2])
      #He <- He + sum(r1[[r]][1]=="1",r1[[r]][2]=="1")
    }
    Ho <- Ho/length(r1); Ho <- min(Ho,1-Ho)
    #He <- He/(2*length(r1)); He <- 2*He*(1-He)
    pi.o <- c(pi.o,Ho)
    #pi.e <- c(pi.e,He)
  }
  pos.pi <- cbind(pos,pi.o)
  
  #define windows
  L <- 1e5
  len.win <- 1e3
  w <- seq(1,L,len.win)
  winpi <- NULL; 
  #calculate sum(Ho) for each window
  pm<-1; 
  for(end in 2:length(w)) {
    sum <- 0; 
    while(pos.pi[pm,1] < w[end]){ 
      sum <- sum + pos.pi[pm,2];
      pm  <- pm + 1
    }; 
    winpi <- c(winpi,sum)
  }
  #plot pi/window
  winpi <- winpi/len.win
  plot(winpi,type="b",pch=16,main=sprintf("Obs Heterozygosity per window (%.0fbp)\n%s",len.win,model),xlab="window",ylab="Pi")
  mu<-2e-7; Ne<-4000; theta <- 4*Ne*mu
  abline(h=theta,col="blue")
}
dev.off()
