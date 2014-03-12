# Here we describe how we use MimicrEE to simulate a data set to test 
# the performance of CMH and BBGP.
# Created by Agnes Jonas, 6/03/2014 


# cd pathTo/simulation

# The scripts we used are parts of the softwares listed below.
# The full softwares can be downloaded form the following websites:
# MimicrEE: https://code.google.com/p/mimicree/
# PoPoolation2: https://code.google.com/p/popoolation2/

mkdir ./data/selectedSNPs
mkdir ./data/simData
mkdir ./data/simData_Poiss
mkdir ./data/cmhInput
mkdir ./data/bbgpInput
mkdir ./data/cmhOutput
mkdir ./data/cmhMHT

# Base haplotypes and the recombination rates can be downloaded from: http://i122server.vu-wien.ac.at/pop/Kosiol_website/GP/

curl 'http://i122server.vu-wien.ac.at/pop/Kosiol_website/GP/N1000-hzc200-lrf.mimhap.gz' | gunzip > ./data/N1000-hzc200-lrf.mimhap
curl 'http://i122server.vu-wien.ac.at/pop/Kosiol_website/GP/dmel.rr.txt' > ./data/dmel.rr.txt

#--------------------------------
# Pick 100 selected SNPs 3 times:
#--------------------------------

# pick-freq-dependent-random-addtive-snps.py choose random SNPs and assign selective advantage to them
# We want to find these truly selected SNPs with the proposed methods.

# Note the amount of produced data can be very large. For shorter the running times and less data, change the following options:
# use the provided sample haplotype file: N1000-hzc200-lrf.mimhap_sample 
# --loci-count 2000000 ==> 100
# number of selected SNPs: -n 100 ==> 10


for i in {1..3}; do python ./pick-freq-dependent-random-addtive-snps.py -n 100 -s 0.1 -e 0.5 --min-frequency 0.12 -m 0.8 --loci-count 2000000 --input ./data/N1000-hzc200-lrf.mimhap > ./data/selectedSNPs/s01-n$i.txt; done

#-------------------------------
# Start simulation: 5 replicates
#-------------------------------

# Input haplotypes: ./data/N1000-hzc200-lrf.mimhap
# Recombination rate: ./data/dmel.rr.txt
# For more informations on the parameters please visit: https://code.google.com/p/mimicree/wiki/ManualMimicrEESummary

# The set of selected SNPs used for the simulations can be downloaded from: http://i122server.vu-wien.ac.at/pop/Kosiol_website/GP/

curl 'http://i122server.vu-wien.ac.at/pop/Kosiol_website/GP/selectedSNPs_original.zip' > data/selectedSNPs_original.zip 
cd ./data
tar -xzf selectedSNPs_original.zip
cd ..

# If you want to use the same set of selected SNPs please modify the following flag:
# --additive ./data/selectedSNPs/s01-n$i.txt ==> ./data/selectedSNPs_original/s01-n$i.txt

for i in {1..3}; do nohup java -Xmx8g -jar  ./MimicrEESummary.jar --haplotypes-g0 ./data/N1000-hzc200-lrf.mimhap --recombination-rate ./data/dmel.rr.txt --output-mode 2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60 --replicate-runs 5 --output-format sync --threads 8 --additive ./data/selectedSNPs/s01-n$i.txt --output-file ./data/simData/allRep_s01-n$i.sync; done

#----------------------------------------
# Add 3-fold Poisson coverage to the data
#----------------------------------------

# poisson-3fold-sample.py Poisson sample the simulated counts to provide coverage information

for i in {1..3}; do python ./poisson-3fold-sample.py --input ./data/simData/allRep_s01-n$i.sync --coverage 45 >> ./data/simData_Poiss/allRep_poissCoverage_s01-n$i.sync; done

#-----------------------------------
# Create BBGP input by adding header
#-----------------------------------

for i in {1..3}; do awk 'BEGIN{FS="\t";OFS="\t"; print "chro","pos","ref",0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60}{print $0}' ./data/simData_Poiss/allRep_poissCoverage_s01-n$i.sync > ./data/bbgpInput/allRep_poissCoverage_withHeader_s01-n$i.sync; done

#------------------------------
# Perform CMH with 5 replicates
#------------------------------

# CMH test can be performed with various number of replicates. For more information on the parameters please visit https://code.google.com/p/popoolation2/wiki/Manual#cmh-test.pl

for i in {1..3}; do awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$34,$35,$65,$66,$96,$97,$127,$128,$158}' ./data/simData_Poiss/allRep_poissCoverage_s01-n$i.sync > ./data/cmhInput/BE_allRep_poissCoverage_s01-n$i.sync; done

for i in {1..3}; do perl ./cmh-test.pl --input ./data/cmhInput/BE_allRep_poissCoverage_s01-n$i.sync --output ./data/cmhOutput/BE_allRep_poissCoverage_s01-n$i.pval --min-count 2 --min-coverage 10 --max-coverage 2000 --population 1-2,3-4,5-6,7-8,9-10; done

for i in {1..3}; do awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$NF}' ./data/cmhOutput/BE_allRep_poissCoverage_s01-n$i.pval > ./data/cmhOutput/BE_allRep_poissCoverage_s01-n$i.pvalOnly; done

#-------------------------------------------
# Visualize CMH results with Manhattan plots
#-------------------------------------------

for i in {1..3}; do R --vanilla --args ./data/cmhOutput/BE_allRep_poissCoverage_s01-n$i.pvalOnly ./data/selectedSNPs/s01-n$i.txt ./data/cmhMHT/mhtCMH_n$i < ./plotMHTsim.R; done


