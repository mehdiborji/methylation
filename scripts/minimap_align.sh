#!/bin/bash

echo 'mode =' $1
echo 'ref =' $2
echo 'dir=' $3
echo 'read1=' $4
echo 'read2=' $5
echo 'outfile=' $6

# make index fist minimap2 -x sr -t 8 -d GRCh38_v44_chrs2.mmi GRCh38_v44_chrs.fasta

minimap2 -aY --eqx -x $1 -F 1000 -t 8 --secondary=no --sam-hit-only $2 $3/$4 $3/$5 > $3/$6.sam
samtools view -@8 -o $3/$6.bam $3/$6.sam
rm $3/$6.sam
#samtools sort -@20 -o $3/$6.bam $3/$6.sam
#samtools index -@20 $3/$6.bam
#samtools view -@20 -c $3/$6.bam
#samtools view -h -b -F 2308 -@16 $3/$5.bam > $3/$5_pri.bam
#samtools view -@16 -c $3/$5_pri.bam
#samtools index -@16 $3/$5_pri.bam
