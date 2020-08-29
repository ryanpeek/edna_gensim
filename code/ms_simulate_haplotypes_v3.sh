#!/bin/bash -l

nInd=$1 # number of individuals in a sample
nLoci=100 # $2 number of loci analysed
ll=10000 # $3 # locus length
cov=100 #$4 # mean coverage per locus analyzed (100x)
distrib=$2 # one of "beta", "gamma", "norm", "unif", "equal", "betaU", "gammaX"
reps=100 #number of simulations to run (had set to 1000)
resDir=$3 # results dir to create

# simplified call
# sh 5 norm sim003
# sh code/ms_simulate_haplotypes_v3.sh 10 norm sim004

# calc ploidy and gene copies for ms
ploidy=2 # ploidy of individuals
gc=$(($nInd * $ploidy)) # gene copies


if [ ! -d "results" ]; then
        echo "create results"
        mkdir results  ### All output goes here ###
fi

if [ ! -d "results/${resDir}" ]; then
        mkdir results/${resDir}  ### All output goes here ###
fi

echo "Starting Job: "
date

echo "Saving to output to: results/${resDir} \n"

echo "Working on $distrib distribution & using gene copies=$gc"

nInd=$(printf "%03g" $nInd)

echo "\nGenerating distributions for $nInd individuals \n"

echo "Calculating theta...\n"

# calculate theta
mu=1e-8 # per bp mutation rate
Ne=25000 # coalescent effect size of population
theta="4*$Ne*$ll*$mu"
theta=`perl -E "say ${theta}"`

# sequencing uncertainty
error=0.1 # error rate, 0.1 equals 10%
thresh=2 # number of times haplotype must be seen to count

echo "Looping through $reps replicates and simulating...\n"
# loop through and run
x=1
while [ $x -le $reps ]	
do
	#printf -v repx "%03g" $x
  
  # CLUSTER: 
  # ms simulations
  #ms $gc $nLoci -t $theta > results/ms_${nInd}_${nLoci}_${theta}_${cov}_$(printf "%03g" $x).ms

  # generate distributions for sampling
  #Rscript --vanilla generate_distributions.R $nInd $distrib $resDir $(printf "%03g" $x)
  
  # sample for haplotypes using distributions above
  # perl ms_sample_v3.pl results/ms_${nInd}_${nLoci}_${theta}_${cov}_$(printf "%03g" $x).ms $cov $error $thresh results/${resDir}/dist_${distrib}_ind${nInd}_$(printf "%03g" $x).txt > results/${resDir}/sample_${nInd}_${nLoci}_${theta}_${cov}_${error}_${thresh}_${distrib}_$(printf "%03g" $x).out
  
  # LOCAL: 
  # ms simulations
  ~/Downloads/msdir/ms $gc $nLoci -t $theta > results/ms_${nInd}_${nLoci}_${theta}_${cov}_$(printf "%03g" $x).ms
	
  # generate distributions for sampling
  Rscript --vanilla code/generate_distributions.R $nInd $distrib $resDir $(printf "%03g" $x)
  
  # sample for haplotypes using distributions above
  perl code/ms_sample_v3.pl results/ms_${nInd}_${nLoci}_${theta}_${cov}_$(printf "%03g" $x).ms $cov $error $thresh results/${resDir}/dist_${distrib}_ind${nInd}_$(printf "%03g" $x).txt > results/${resDir}/sample_${nInd}_${nLoci}_${theta}_${cov}_${error}_${thresh}_${distrib}_$(printf "%03g" $x).out
	
	x=$(( $x + 1 ))

done

echo "Simulations Done! Cleaning up..."

paste results/${resDir}/sample_${nInd}_${nLoci}_${theta}_${cov}_${error}_${thresh}_${distrib}*.out > results/${resDir}/run_${nInd}_${nLoci}_t${theta}_${cov}_${error}_${thresh}_${distrib}_reps${reps}.out


if [ ! -d "results/${resDir}/distribs" ]; then
        mkdir results/${resDir}/distribs  ### All output goes here ###
fi

# move distributions files
mv results/${resDir}/dist_*.txt results/${resDir}/distribs/

# rm temp files
echo "\tRemoving sample and ms out temp files \n"
rm seedms
rm results/ms_*.ms
#rm results/${resDir}/dist*.txt
rm results/${resDir}/sample_*.out

echo "Done!!"
date


