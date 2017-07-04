#!/bin/bash

# rmd last one

CD=`pwd`
ARRAY=(adp atp)

for (( l=0; l<${#ARRAY[@]}; l++ )); do
rm -f DSSP_averages_${ARRAY[l]}.txt
grep -v \# lava_${ARRAY[l]}.dssp_half1_extend.xvg | grep -v \@ | tail -n 1000 > ${ARRAY[l]}.dssp.half1_extend.xvg #last 100ns only
ave1=`awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' ${ARRAY[l]}.dssp.half1_extend.xvg`
grep -v \# lava_${ARRAY[l]}.dssp_half1.xvg | grep -v \@ | tail -n  1000 > ${ARRAY[l]}.dssp.half1.xvg
ave2=`awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' ${ARRAY[l]}.dssp.half1.xvg`
grep -v \# lava_${ARRAY[l]}.dssp_half2.xvg | grep -v \@ | tail -n 1000 > ${ARRAY[l]}.dssp.half2.xvg
ave3=`awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' ${ARRAY[l]}.dssp.half2.xvg`
echo ${ARRAY[l]} $ave1 $ave2 $ave3 >> DSSP_averages_${ARRAY[l]}.txt	
done
