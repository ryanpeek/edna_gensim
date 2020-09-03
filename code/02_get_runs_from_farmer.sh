#!/bin/bash

# for testing
#ls run*.out
#get run*.out  

sftp -P 2022 rapeek@farm.cse.ucdavis.edu <<EOF
cd projects/edna_pop/results/sim005
lcd results/sim005
get run*.out 

exit
EOF
