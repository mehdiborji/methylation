import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-p", "--parts_batch", type=int)
args = parser.parse_args()

indir = args.indir
sample = args.sample
parts_batch = args.parts_batch
######################################################

multiplier_per_task = 15

start = (parts_batch - 1) * multiplier_per_task
end = parts_batch * multiplier_per_task

for b in range(start, end):
    
    p = str(b).zfill(3)
    
    fq1 = f"{indir}/{sample}/split/{sample}_R1.part_{p}_clean.fastq"
    fq2 = f"{indir}/{sample}/split/{sample}_R3.part_{p}_clean.fastq"
    
    
    out_fq1 = f"{indir}/{sample}/split/{sample}_R1.part_{p}_clean_trim.fastq"
    out_fq2 = f"{indir}/{sample}/split/{sample}_R3.part_{p}_clean_trim.fastq"
    
    out_json = f"{indir}/{sample}/split/{sample}_part_{p}.json"
    out_html = f"{indir}/{sample}/split/{sample}_part_{p}.html"
    
    #options = "--adapter_sequence=CTATCTCTTATACACATCTCCAAACC --adapter_sequence_r2=CTGTCTCTTATACACATCTGACGCTG"
    
    options = "--adapter_sequence=CTGTCTCTTATACACATCT --adapter_sequence_r2=CTGTCTCTTATACACATCT"
    
    input_output = f"-i {fq1} -I {fq2} -o {out_fq1} -O {out_fq2} -j {out_json} -h {out_html}"
    
    command = f"/home/meb521/fastp -w 4 {options} {input_output}"
    
    print(command)
    
    try:
        subprocess.run(command, shell=True)
    except Exception as e:
        print(f"A {e} occurred")