#!/bin/bash
#$ -pe smp-verbose 8 -o ~/Meshfree/preform.out_$JOB_ID
echo Running preform program    
mkdir ~/Meshfree/$JOB_ID
rm -r ~/Meshfree/$JOB_ID/Displacement
rm -r ~/Meshfree/$JOB_ID/History
rm -r ~/Meshfree/$JOB_ID/srRod
mkdir ~/Meshfree/$JOB_ID/History
cp ~/MFREEX-master/build/bin/preform/preform.nodes ~/Meshfree/$JOB_ID/preform.nodes
cp ~/MFREEX-master/build/bin/preform/preform.boundary ~/Meshfree/$JOB_ID/preform.boundary
cd ~/Meshfree/$JOB_ID; ~/MFREEX-master/build/bin/preform/preform
echo Finished running preform program 