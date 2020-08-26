#!/bin/bash -l

nInd=$1 # number of individuals in a sample
nLoci=$2 # number of loci analysed
ll=$3 # locus length
cov=$4 # mean coverage per locus analyzed
distrib=$5 # one of "beta", "gamma", "norm", "unif"
reps=$6 #number of simulations to run (had set to 1000)

# so call for 5 indiv, 100 loci, 10000 loci length, 100x coverage, distrib, simulations
# sh ms_simulate_haplotypes.sh 5 100 10000 100 norm 5

# calc ploidy and gene copies for ms
ploidy=2 # ploidy of individuals
let gc=$nInd*$ploidy # gene copies
echo "Gene copies=$gc"

if [ ! -d "results" ]; then
        echo "create results"
        mkdir results  ### All output goes here ###
fi
echo "saving to dir: results"

echo "Starting Job: "
date

echo "Generating distributions..."

# now call R script to generate distributions for sampling
Rscript code/generate_distributions_txt.R $nInd

# these get generated into the "results" folder

# calculate theta
mu=1e-8 # per bp mutation rate
Ne=25000 # coalescent effect size of population
theta="4*$Ne*$ll*$mu"
theta=`perl -E "say ${theta}"`

# sequencing uncertainty
error=0.1 # error rate, 0.1 equals 10%
thresh=2 # number of times haplotype must be seen to count


# loop through and run
x=1
while [ $x -le $reps ] 
do
  # local
	~/Downloads/msdir/ms $gc $nLoci -t $theta > results/ms_${nInd}_${nLoci}_${theta}_${cov}_${x}.ms
	
	# cluster
	#ms $gc $nLoci -t $theta > results/ms_${nInd}_${nLoci}_${theta}_${cov}_${x}.out
	
	perl code/ms_sample_2.pl results/ms_${nInd}_${nLoci}_${theta}_${cov}_${x}.ms $cov $error $thresh results/dist_${distrib}_ind${nInd}.txt > results/sample_${nInd}_${nLoci}_${theta}_${cov}_${error}_${thresh}_${distrib}_${x}.out
	
	x=$(( $x + 1 ))

done

paste results/sample_${nInd}_${nLoci}_${theta}_${cov}_${error}_${thresh}_${distrib}*.out > results/run_${nInd}_${nLoci}_t${theta}_${cov}_${error}_${thresh}_${distrib}_reps${reps}.out

# rm temp files
echo "Removing sample and ms out temp files"
rm seedms
# rm results/sample_*.out
# rm results/ms_*.ms

echo "Done!!"

# to run: sh simulate_haplotypes.sh

