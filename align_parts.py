import argparse
import methyl_utils
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-r", "--ref_dir", type=str)
parser.add_argument("-p", "--parts_batch", type=int)
args = parser.parse_args()

indir = args.indir
sample = args.sample
ref_dir = args.ref_dir
parts_batch = args.parts_batch

######################################################

print("current dir:", os.getcwd())
os.chdir(indir)
print("current dir changed to:", os.getcwd())
# ref_dir = "/n/scratch/users/m/meb521/GRCh38_v44/"

command = "/home/meb521/methylation/scripts/bismark_align.sh"

p = str(parts_batch - 1).zfill(3)

fq1 = f"{sample}/split/{sample}_R1.part_{p}_clean_trim.fastq"
fq2 = f"{sample}/split/{sample}_R3.part_{p}_clean_trim.fastq"
outdir = f"{sample}/split/output_{p}"

bam = f"{sample}/split/output_{p}/{sample}_R1.part_{p}_clean_trim_bismark_bt2_pe.bam"

mbias = f"{sample}/split/output_{p}/{sample}_R1.part_{p}_clean_trim_bismark_bt2_pe.M-bias.txt"

mbias_exists = os.path.isfile(mbias)
if not mbias_exists:
    submit = f"{command} {fq1} {fq2} {outdir} {ref_dir} {bam}"
    print(submit)
    subprocess.run(submit, shell=True)
