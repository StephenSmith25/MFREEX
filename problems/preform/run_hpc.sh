#!/bin/sh
#$ -pe smp-verbose 8 -o ~/Meshfree/meshfree.out.$JOB_ID
#$ -M ssmith54@qub.ac.uk –m bea
START=$(date +%s)
echo Running preform program	
rm -r Displacement
rm -r History
rm -r srRod
mkdir History
./preform
echo Finished running preform program 
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
