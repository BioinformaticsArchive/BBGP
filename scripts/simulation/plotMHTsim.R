args<-commandArgs(TRUE)
alen=length(args)
description=args[alen]

data<-read.table(args[1])
sel<-read.table(args[2])
outfile<-args[3]

names(data)<-c('CHR','BP','P')
sel$V3=NULL
sel$V4=NULL
sel$V5=NULL
sel$V6=NULL
names(sel)<-c('CHR','BP')
true<-merge(data,sel,c("CHR","BP"))

source("../manhattan_Dmel.R")

quartz(height=6,width=12)
par(bg="white")
par(mar=c(5,5,4,2))

manhattan(data,yLabel="-log10(P-value)", colors=c("black","slategrey"),limitchromosomes=c("X","2L","2R","3L","3R","4"),main="",nplog=FALSE, cex.lab=2, cex.axis=2)

# Plot truly selected SNPs:

points(true[true$CHR=="2L",]$BP, true[true$CHR=="2L",]$P, pch=20, cex=1, col="red")

points(true[true$CHR=="2R",]$BP+max(data[data$CHR=="2L",]$BP), true[true$CHR=="2R",]$P, pch=20, cex=1, col="red")

points(true[true$CHR=="3L",]$BP+max(data[data$CHR=="2L",]$BP)+ max(data[data$CHR=="2R",]$BP), true[true$CHR=="3L",]$P, pch=20, cex=1, col="red")

points(true[true$CHR=="3R",]$BP+max(data[data$CHR=="2L",]$BP)+ max(data[data$CHR=="2R",]$BP)+ max(data[data$CHR=="3L",]$BP), true[true$CHR=="3R",]$P, pch=20, cex=1, col="red")

dev.print(png, file=paste(outfile,".png",sep=""),width=1024)

