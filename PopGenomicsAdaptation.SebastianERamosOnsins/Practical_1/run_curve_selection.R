#read files from slim output
curve_files <- system("ls *.practical1.output_curve.txt",intern=T)

pdf("SelectiveSweep_frequency.pdf")
for(model in curve_files) {
  data <- read.table(file=model)
  freq <- data[-1,2]
  if(length(freq))
    plot(x=c(1:length(freq)),y=freq,xlab="Generations",ylab="Frequency",ylim=c(0,1),main=sprintf("Frequency of Selected Position\n%s",model),type="l")
}
dev.off()
