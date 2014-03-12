##############################################################################
###Manhattan plot function 
###from http://gettinggeneticsdone.blogspot.co.at/2011/04/annotated-manhattan-plots-and-qq-plots.html
###Adapted for Dmel chromosomes
##############################################################################

manhattan <- function(dataframe, yLabel=expression(-log[10](italic(P))),dCEX=1, colors=c("gray60","black"), ymax="max", limitchromosomes=c("X","2L","2R","3L","3R","4","2LHet","2RHet","3LHet","3RHet","XHet","YHet"), suggestiveline=NULL, genomewideline=NULL, annotate=NULL, cand_col="black", ymin=0, pointscale="linear", ps_max = 1, ps_min = 0.1, p_ch=20, nplog=TRUE, ...) {
	# to only get a region give limitchromosomes=c("X:10000-200000")
	# to change the scaling of points give pointscale either as "linear", "none" to turn it off 
	# or to a threshold value underneath which cex=ps_min and above cex=ps_max
	# to change the printchar set p_ch to something else
	# to annotate a region or a point give in the following format
	# annotate=c("2L:10000-15000,2R:190000-20000","X:110000"), cand_col=c("red","green")
	# the cand_color vector must be of the same size as the annotate vector
	# if you do not want to have a logarithmic scale give nplog=FALSE and set y min to something low
	d=dataframe
	if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    chromregion = NULL 
    if (length(limitchromosomes)>0) {
    	# check if region for chromosome given
    	if (length(limitchromosomes)==1 && grepl('\\w+\\:\\d+\\-\\d+',limitchromosomes,perl=TRUE)){
    		limitchromosomes = unlist(strsplit(limitchromosomes,"[:-]+"))
    		chromregion = as.numeric(limitchromosomes[2:3])
    		limitchromosomes <- limitchromosomes[1]
    	}
    	d=d[d$CHR %in% limitchromosomes, ]
    	   	
    	}
    ymin<- ifelse(nplog,10^-(ymin),ymin)
	# order by chromosome vector
	d$CHR <- factor(d$CHR,levels<-limitchromosomes)
	if (nplog) { d = subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=ymin))}
	else {d = subset(na.omit(d[order(d$CHR, d$BP), ]), (P > ymin)) }
    # remove nas, sort, and keep only 0<P<=10^ymin for nplog
    # if positions given extract only that region
    if ( length(chromregion) > 0) d = subset(d, (d$BP >= chromregion[1]  & d$BP <= chromregion[2]))
	if(nplog) {d$logp = -log10(d$P)}
	else {d$logp = d$P}
	d$pos=NA
	ticks=NULL
	lastbase=0
	dmax = ceiling(max(d$logp))
	if (pointscale == "linear") { d$CEX= (ps_min + (ps_max-ps_min) * d$logp/dmax) }
	else if (pointscale == "none") {d$CEX = dCEX}
	else {
		threshold <- as.numeric(pointscale)
		d$CEX= ifelse(d$logp >= threshold, ps_max, ps_min)
	}
	numchroms=length(unique(d$CHR))
	colors <- rep(colors,numchroms)[1:numchroms]
	# setting the pointsizes in dependence of p Values
	if (ymax=="max") ymax<-ceiling(max(d$logp))
	#if (ymax<8) ymax<-8
	print("Setting up plot... ")
	if (numchroms==1) {
		d$pos=d$BP
		ticks=floor(length(d$pos))/2+1
		} 
	else {
		cc=0
		for (i in unique(d$CHR)) {
			cc=cc+1
			if (cc==1) {
				d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
			} else {
				lastbase=lastbase+tail(subset(d,CHR==unique(d$CHR)[cc-1])$BP, 1)
				d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
			}
			ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
		}
	}
	print("Plotting points... ")
	if (numchroms==1) {
		with(d, plot(pos, logp, ylim=c(0,ymax), ylab=yLabel, xlab=unique(d$CHR), col=colors, cex.lab=1.1, pch=p_ch, cex=d$CEX, ...))
		}	else {
		with(d, plot(pos, logp, ylim=c(0,ymax), ylab=yLabel, xlab="Chromosome", xaxt="n", type="n", pch=p_ch, cex=d$CEX, ...))
		axis(1, at=ticks, lab=unique(d$CHR), hadj=0.0005, padj=0.0005, ...)
		icol=1
		for (i in unique(d$CHR)) {
			with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], pch=p_ch,cex=d[d$CHR==i, ]$CEX, ...))
			icol=icol+1
			}
		}

    if (!is.null(annotate)) {
    	nn=0
    	for(set in annotate){
    		nn=nn+1
    		print("Annotating candidates... ")
    		cand_regions = unlist(strsplit(set,"[,]+"))
    		for (i in cand_regions) {
    			region = unlist(strsplit(i,"[:-]+"))
    			if (length(region) == 2 ) {region[3] = region[2] }
    			chromosome = region[1]
    			chrom_region = as.numeric(region[2:3])
				d.annotate=d[which(d$CHR == chromosome & chrom_region[1] <= d$BP & d$BP <= chrom_region[2]),]
    			with(d.annotate, points(pos, logp, col=cand_col[nn], cex=CEX, pch=p_ch, ...)) 
    		

    	 	}
    	}
    }	
    if(!is.null(suggestiveline)) abline(h=suggestiveline, col="blue")
    if(!is.null(genomewideline)) abline(h=genomewideline, col="black")
}

