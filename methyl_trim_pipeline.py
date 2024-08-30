import argparse
import methyl_utils
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
args = parser.parse_args()

indir = args.indir
sample = args.sample

######################################################

parts = methyl_utils.find_sub_fastq_parts(indir, sample)

print(f"found {len(parts)} parts")

parts = methyl_utils.find_sub_fastq_parts(indir, sample)

for p in parts:
    
    fq1 = f"{indir}/{sample}/split/{sample}_R1.part_{p}_clean.fastq"
    fq2 = f"{indir}/{sample}/split/{sample}_R3.part_{p}_clean.fastq"
    
    
    out_fq1 = f"{indir}/{sample}/split/{sample}_R1.part_{p}_clean_trim.fastq"
    out_fq2 = f"{indir}/{sample}/split/{sample}_R3.part_{p}_clean_trim.fastq"
    
    out_json = f"{indir}/{sample}/split/{sample}_part_{p}.json"
    out_html = f"{indir}/{sample}/split/{sample}_part_{p}.html"
    
    options = "--adapter_sequence=CTATCTCTTATACACATCTCCAAACC --adapter_sequence_r2=CTGTCTCTTATACACATCTGACGCTG"
    
    input_output = f"-i {fq1} -I {fq2} -o {out_fq1} -O {out_fq2} -j {out_json} -h {out_html}"
    
    command = f"/home/meb521/fastp -w 10 {options} {input_output}"
    
    print(command)
    
    try:
        subprocess.run(command, shell=True)
    except Exception as e:
        print(f"A {e} occurred")