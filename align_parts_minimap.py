import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-r", "--ref_index", type=str)
parser.add_argument("-p", "--parts_batch", type=int)

args = parser.parse_args()

indir = args.indir
sample = args.sample
ref_index = args.ref_index
parts_batch = args.parts_batch

######################################################

print("current dir:", os.getcwd())
os.chdir(indir)
print("current dir changed to:", os.getcwd())

command = "/home/meb521/methylation/scripts/minimap_align.sh"

# this many of the fastq parts will be aligned in a `loop` in each task of the array job
# total fastq parts must be <= multiplier_per_task * array_job total tasks
# x153 DNA library with 1.43b total and 10m per part reads: 143 < 9 * 16

multiplier_per_task = 6

start = (parts_batch - 1) * multiplier_per_task
end = parts_batch * multiplier_per_task

for part in range(start, end):
    
    p = str(part).zfill(3)

    outdir = f"{sample}/split"
    fq1 = f"{sample}_R1.part_{p}_clean_trim.fastq"
    fq2 = f"{sample}_R3.part_{p}_clean_trim.fastq"
    file_name = f"{sample}_part_{p}"
    
    bam_index = f"{sample}/split/{file_name}.bam.bai"
    bam_index_exists = os.path.isfile(bam_index)

    if not bam_index_exists:
        submit = f"{command} sr {ref_index} {outdir} {fq1} {fq2} {file_name}"
        print(submit)
        subprocess.run(submit, shell=True)
    else:
        print(bam_index,'exists, alignment seems finished')
