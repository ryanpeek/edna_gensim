#!/bin/bash -l

nInd=$1 # number of individuals in a sample
nLoci=$2 # number of loci analysed
ll=$3 # locus length
cov=$4 # mean coverage per locus analyzed
reps=1000 #number of simulations to run

# calc ploidy and gene copies for ms
ploidy=2 # ploidy of individuals
gc=$nInd*$ploidy # gene copies

if [ ! -d "results" ]; then
        echo "create results"
        mkdir results  ### All output goes here ###
fi
echo "saving to dir: results"

echo "Starting Job: "
date

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

	#LOCAL CALL ~/Downloads/msdir/ms $gc $nLoci -t $theta > results/ms_${nInd}_${nLoci}_${theta}_${x}.out
	ms $gc $nLoci -t $theta > results/ms_${nInd}_${nLoci}_${theta}_${cov}_${x}.out
	perl ms_sample.pl results/ms_${nInd}_${nLoci}_${theta}_${cov}_${x}.out $cov $error $thresh > results/sample_${nInd}_${nLoci}_${theta}_${cov}_${error}_${thresh}_${x}.out
	x=$(( $x + 1 ))

done

paste results/sample_${nInd}_${nLoci}_${theta}_${cov}_${error}_${thresh}*.out > results/run_${nInd}_${nLoci}_t${theta}_${cov}_${error}_${thresh}_reps${reps}.out

# rm temp files
echo "Removing sample and ms out temp files"
rm seedms
rm results/sample_*.out
rm results/ms_*.out

echo "Done!!"

# to run: sh simulate_haplotypes.sh

