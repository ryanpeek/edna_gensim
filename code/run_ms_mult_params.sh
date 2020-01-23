#!/bin/bash -l
# 2019 Ryan Peek

if [ "$#" -lt 3 ]; then
        echo -e "Runs ms multiple times with mulitple parameter sets.\n\nRequired Arguments:"
        echo -e "\tParameter File"
        echo -e "\tsbatch time argument"
        echo -e "\tOutput Filename"
        echo -e "\nParameter File contains the following information, delimited by a single space. Multiple lines can be given for multiple parameter sets to run."
        echo -e "\tnumber of individuals in sample"
        echo -e "\tnumber of loci analyzed"
        echo -e "\tlocus length"
        echo -e "\tmean coverage per locus analyzed"
        exit
fi

# input parameters
file=$1
stime=$2
output=$3

echo -e "\n\nInput file: ${file}.\nRuntime: ${stime}.\nOutput: ${output}.\n\n\n"

# get end from number of lines in file
end=`wc -l ${file}`
end=${end% *}

echo -e "Submitting:"

x=1
while [ $x -le $end ]
do
  # get variables from current line of input file
  line=`sed "${x}q;d" ${file}`
	#read -a ARRAY <<< ${line} # a here string to make array
	ARRAY="echo ${line}";
	# could also use arr=(`echo ${line}`);
	nInd="${ARRAY[0]}"
  nLoci="${ARRAY[1]}"
  ll="${ARRAY[2]}"
  cov="${ARRAY[3]}"
  mu=1e-8 # per bp mutation rate
  Ne=25000 # coalescent effect size of population
	theta="4*$Ne*$ll*$mu" # calculate theta based on ll, Ne, mu
	theta=`perl -E "say ${theta}"`
	ploidy=2 # ploidy of individuals
	gc=$nInd*$ploidy # calculate gene copies
	
	# CAN CHANGE THIS TO CHANGE # of sims
	reps=1000 #number of simulations to run

  echo -e "Run: ${x}. Parameters:\n\tnIndividuals: ${nInd}.\n\tnLoci: ${nLoci}.\n\tLociLength: ${ll}.\n\tCoverage: ${cov}.\n\tTheta: ${theta}.\n\tIterations: ${reps}.\n" 

        #run sims, assumes ms in bin
	echo "#!/bin/bash -l

	x=1
	while [ $x -le $reps ]
	do

        	ms $gc $nLoci -t $theta > results/ms_${nInd}_${nLoci}_${theta}_${x}.out
        	perl ms_sample.pl results/ms_${nInd}_${nLoci}_${theta}_${x}.out $cov > results/sample_${nInd}_${nLoci}_${theta}_${x}.out
	
        	x=$(( $x + 1 ))

	done

	paste results/sample_${nInd}_${nLoci}_${theta}_*.out > results/run_${nInd}_${nLoci}_t${theta}_reps${reps}.out
	
	# rm temp files
	echo "Removing sample and ms out temp files"
	rm seedms
	rm results/sample_*.out
	rm results/ms_*.out

	echo 'Done!!'"  > ${nInd}_${nLoci}_${theta}_${ll}_${cov}/ms_run.sh
        #sbatch -t $stime -p high --mem=16G ${nInd}_${nLoci}_${theta}_${ll}_${cov}/ms_run.sh

        echo -e "\n"

        x=$(( $x + 1 ))
done

