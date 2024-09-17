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

print("current dir:", os.getcwd())
os.chdir(indir)
print("current dir changed to:", os.getcwd())
# ref_dir = "/n/scratch/users/m/meb521/GRCh38_v44/"

p = str(parts_batch - 1).zfill(3)


bam = f"{sample}/split/output_{p}/{sample}_R1.part_{p}_clean_trim_bismark_bt2_pe.bam"
command = "/home/meb521/Bismark-0.24.2/filter_non_conversion"
options = "--paired --threshold 5 --consecutive"

submit = f"{command} {options} {bam}"
print(submit)
subprocess.run(submit, shell=True)
