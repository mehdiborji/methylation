#!/bin/bash

read1=$1
read2=$2
out_dir=$3
ref_dir=$4

echo 'read1 =' $1
echo 'read2 =' $2
echo 'out_dir =' $3
echo 'ref_dir =' $4
    
/home/meb521/Bismark-0.24.2/bismark \
    --pbat \
    --score_min L,0,-0.4 \
    $ref_dir \
    -1 $read1 \
    -2 $read2 \
    -o $out_dir