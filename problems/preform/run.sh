#!/bin/sh
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
