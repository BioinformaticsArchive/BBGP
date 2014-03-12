args<-commandArgs(TRUE)
alen=length(args)
description=args[alen]

data<-read.table(args[1])
outfile<-args[2]

names(data)<-c('CHR','BP','P')

source("../manhattan_Dmel.R")

quartz(height=6,width=12)
par(bg="white")
par(mar=c(5,5,4,2))

manhattan(data,yLabel="-log10(P-value)", colors=c("black","slategrey"),limitchromosomes=c("X","2L","2R","3L","3R","4"),main="",nplog=TRUE, cex.lab=2, cex.axis=2)

dev.print(png, file=paste(outfile,".png",sep=""),width=1024)

