import argparse
import methyl_utils
import os
import subprocess
import concurrent.futures

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-c", "--cores", type=int)
args = parser.parse_args()

indir = args.indir
sample = args.sample
cores = args.cores

######################################################

parts = methyl_utils.find_sub_fastq_parts(indir, sample)

print(f"found {len(parts)} parts")

print("current dir:", os.getcwd())
os.chdir(f"{indir}/{sample}/split")
print("current dir changed to:", os.getcwd())
# ref_dir = "/n/scratch/users/m/meb521/GRCh38_v44/"

def run_sed_fastq(fastq_file):
    
    backup_file = fastq_file + '.bak'
    
    backup_exists = os.path.isfile(backup_file)
    
    if not backup_exists:
        if 'R1' in fastq_file:
            command  = f"sed -i.bak '/^@/s/_1_/_/g' {fastq_file}"
        else:
            command  = f"sed -i.bak '/^@/s/_2_/_/g' {fastq_file}"

        print(command)
        subprocess.run(command, shell=True)
    else:
        print(backup_file,'exists')

R1_files = [f"{sample}_R1.part_{p}_clean_trim.fastq" for p in parts]
R3_files = [f"{sample}_R3.part_{p}_clean_trim.fastq" for p in parts]

with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:

    executor.map(run_sed_fastq, R1_files)
    
    executor.map(run_sed_fastq, R3_files)
