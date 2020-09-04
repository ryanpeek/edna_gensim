#!/bin/bash

# move to the correct local folder
cd results/sim005

# this had to be run locally on farm

# 01. FIND THE LAST FILE LOCALLY (BY TIME)
lastfile=$(ls -ltr run* | tail -1 | awk '{print $9}')
echo "last file: $lastfile"


## TO GET ALL FILES
sftp -P 2022 rapeek@farm.cse.ucdavis.edu <<EOF
cd projects/edna_pop/results/sim005
get run_*
exit
EOF

