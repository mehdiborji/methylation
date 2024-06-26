import argparse
import methyl_utils
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-r", "--ref_dir", type=str)
parser.add_argument("-j", "--job_submit", default=False, action="store_true")
parser.add_argument("-b", "--sbatch", default=False, action="store_true")
args = parser.parse_args()

indir = args.indir
sample = args.sample
ref_dir = args.ref_dir
job_submit = args.job_submit
sbatch = args.sbatch

######################################################

parts = methyl_utils.find_sub_fastq_parts(indir, sample)

print(f"found {len(parts)} parts")

print("current dir:", os.getcwd())
os.chdir(indir)
print("current dir changed to:", os.getcwd())
# ref_dir = "/n/scratch/users/m/meb521/GRCh38_v44/"

command = "/home/meb521/methylation/scripts/bismark_align_methex.sh"

for p in parts:
    fq1 = f"{sample}/split/{sample}_R1.part_{p}_clean.fastq"
    fq2 = f"{sample}/split/{sample}_R3.part_{p}_clean.fastq"
    outdir = f"{sample}/split/output_{p}"

    bam = f"{sample}/split/output_{p}/{sample}_R1.part_{p}_clean_bismark_bt2_pe.bam"

    mbias = f"{sample}/split/output_{p}/{sample}_R1.part_{p}_clean_bismark_bt2_pe.M-bias.txt"

    mbias_exists = os.path.isfile(mbias)
    if not mbias_exists:
        if sbatch:
            submit = f"sbatch {command} {fq1} {fq2} {outdir} {ref_dir} {bam}"
        else:
            submit = f"{command} {fq1} {fq2} {outdir} {ref_dir} {bam}"
        print(submit)
        if job_submit:
            subprocess.call(submit, shell=True)
