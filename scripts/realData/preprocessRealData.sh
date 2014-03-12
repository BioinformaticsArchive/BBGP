# Here we describe how the real data was prepared for CMH results and BBGP comparison
# Created by Agnes Jonas, 06/03/2014

#----------------------------------------------------------------------------
# Get real data from: http://datadryad.org/resource/doi:10.5061/dryad.60k68/6
#----------------------------------------------------------------------------

curl -L 'http://datadryad.org/bitstream/handle/10255/dryad.39719/BF37.zip' | unzip

# If the link is not working the data can be downloaded from the Dryad data base (http://datadryad.org/) under doi:10.5061/dryad.60k68/6

# unzip BF37.zip

#------------------------------------------
# Visualize CMH results with Manhattan plot
#------------------------------------------

# The file already contains the p-values from the CMH test in the last column.

mkdir mhtCmh

awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$NF}' BF37 > ./mhtCmh/BF37.pvalOnly

R --vanilla --args ./mhtCmh/BF37.pvalOnly ./mhtCmh/mhtCMH_realData < ./plotMHtrealCMH.R

#------------------
# Create BBGP input
#------------------

# Cut last column:

awk 'BEGIN{FS="\t";OFS="\t"}{$NF = ""; print}' BF37 > BF37_counts

# Filter out rows which contain zero counts

python ./filtZerosSync.py --input BF37_counts --output BF37_filtZeros
rm BF37_counts

# Filter out tri- and quad-allelic cases:

python ./filtTriAllele.py --input BF37_filtZeros --output BF37_filt
rm BF37_filtZeros

# The real data contains the following generation-replicate pairs: BR1,BR3,BR2,F15R4,F15R5,F23R1,F27R5,F37R4,F37R5,F37R1
# Append header:

awk 'BEGIN{FS="\t";OFS="\t"; print "chro","pos","ref",0,0,0,15,15,23,27,37,37,37}{print $0}' BF37_filt > BF37_bbgpInput
rm BF37_filt






