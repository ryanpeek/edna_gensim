#!/bin/bash -l

nInd=$1
nLoci=$2
ll=$3
cov=$4


if [ ! -d "results" ]; then
        echo "create results"
        mkdir results  ### All output goes here ###
fi
echo "directory exists"
echo "Starting Job: "
date

mu=1e-8 #per bp mutation rate
#ll=1000 #locus length
Ne=25000 #coalescent effect size of population

theta="4*$Ne*$ll*$mu"
theta=`perl -E "say ${theta}"`

#nInd=10 #number of individuals in water sample
ploidy=2 #ploidy of individuals
gc=$nInd*$ploidy

#cov=10 #mean coverage per locus analyzed
#nLoci=100 #number of loci analysed
reps=1000 #number of simulations to run


x=1
while [ $x -le $reps ] 
do

	~/Downloads/msdir/ms $gc $nLoci -t $theta > results/ms_${nInd}_${nLoci}_${theta}_${x}.out
	perl ms_sample.pl results/ms_${nInd}_${nLoci}_${theta}_${x}.out $cov > results/sample_${nInd}_${nLoci}_${theta}_${x}.out
	x=$(( $x + 1 ))

done

paste results/sample_${nInd}_${nLoci}_${theta}_*.out > results/run_${nInd}_${nLoci}_t${theta}_reps${reps}.out

# rm temp files
echo "Removing sample and ms out temp files"
rm seedms
rm results/sample_*.out
rm results/ms_*.out

echo "Done!!"

# to run: sh simulate_haplotypes.sh

