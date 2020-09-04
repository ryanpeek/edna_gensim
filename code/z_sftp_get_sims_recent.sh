#!/bin/bash

# move to the correct local folder
cd results/sim005

# this had to be run locally on farm

# 01. FIND THE LAST FILE LOCALLY (BY TIME)
lastfile=$(ls -ltr run* | tail -1 | awk '{print $9}')
echo "last file: $lastfile"

lasttime=$(date -r $lastfile +"%m/%d/%Y %H:%M:%S")
##lasttime="09/03/2020 09:19:10"
echo "last modified file: $lasttime"

## 02. ON FARM: find files that are newer than a date/time and redirect to list
ssh -p 2022 rapeek@farm.cse.ucdavis.edu /bin/bash <<EOF
cd /home/rapeek/projects/edna_pop/results/sim005
find run_* -type f -newermt "$lasttime" > newest
exit
EOF

# 03. GET THE NEWEST LIST locally
sftp -P 2022 rapeek@farm.cse.ucdavis.edu <<EOF
cd projects/edna_pop/results/sim005
get newest
exit
EOF

# 04. RUN THIS TO THEN PULL THE FILES
file='newest'
totnum=$(wc -l newest)
echo $totnum
#totnum=20
fileNames=$(echo "ls -1rt run*" | sftp -P 2022 rapeek@farm.cse.ucdavis.edu:/home/rapeek/projects/edna_pop/results/sim005 | tail -n $totnum) 

for f in $fileNames
do
  #echo "get $f"
  echo "get $f" | sftp -P 2022 rapeek@farm.cse.ucdavis.edu:/home/rapeek/projects/edna_pop/results/sim005
done
  
## OLD: TO GET ALL FILES
# sftp -P 2022 rapeek@farm.cse.ucdavis.edu <<EOF
# cd projects/edna_pop/results/sim005
# get run_*
# exit
# EOF

