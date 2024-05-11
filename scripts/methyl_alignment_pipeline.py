import argparse
import methyl_utils
import os

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cores', type=str)
parser.add_argument('-i', '--indir', type=str)
parser.add_argument('-s', '--sample', type=str)
parser.add_argument('-r', '--ref_dir', type=str)
parser.add_argument('-j', '--job_submit', default=False, action='store_true')
args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
ref_dir = args.ref_dir
job_submit = args.job_submit

######################################################

parts = methyl_utils.find_sub_fastq_parts(indir,sample)
args = [(indir, sample, part, limit) for part in parts]

os.chdir(indir)

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
        submit = f"sbatch {command} {fq1} {fq2} {outdir} {ref_dir} {bam}"
        print(submit)
        if job_submit:
            subprocess.call(submit, shell=True)
