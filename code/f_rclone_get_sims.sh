#!/bin/bash

# move to the correct local folder
#cd results/sim005

## TO GET ALL FILES ONTO BOX
ssh -p 2022 rapeek@agri.cse.ucdavis.edu <<EOF
module load rclone
cd projects/edna_pop/results/sim005
rclone copy -P . --include "run*" box:PROJECTS/SERDP/sim005
exit
EOF

cd ~/Box\ Sync/PROJECTS/SERDP/sim005

ls -1tr run* | wc -l
