import argparse
import pandas as pd
from multiprocessing import Pool
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=str)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample

def allc_process(indir, sample, bc):
    
    allcools_dir = f"{indir}/{sample}/allcools"
    
    cg_file = f'{allcools_dir}/{bc}_CpG_context.tsv'
    ch_file = f'{allcools_dir}/{bc}_Non_CpG_context.tsv'
    allc_file = f'{allcools_dir}/{bc}_allc.tsv'
    
    submit = f"cat {cg_file} {ch_file} > {allc_file}"
    subprocess.run(submit, shell=True)
    
    allc_sorted_file  = f'{allcools_dir}/{bc}_allc_sorted.tsv'
    
    submit = f"sort -k1,1 -k2,2n {ch_file} > {allc_sorted_file}"
    subprocess.run(submit, shell=True)

    submit = f"bgzip -f {allc_sorted_file}"
    subprocess.run(submit, shell=True)

    submit = f"tabix -s1 -b2 -e2 -f {allc_sorted_file}.gz"
    subprocess.run(submit, shell=True)
    
def fix_allc(indir, sample, bc):
    
    print(bc)
    allcools_dir = f"{indir}/{sample}/allcools"
    
    allc_sorted_file  = f'{allcools_dir}/{bc}_allc_sorted.tsv.gz'
    
    allc_sorted_mod_file  = f'{allcools_dir}/{bc}_allc_sorted_mod.tsv.gz'
  
    submit = f"zcat {allc_sorted_file} | sed 's/CH/CA/g' | bgzip > {allc_sorted_mod_file}"
    
    subprocess.run(submit, shell=True)

    submit = f"tabix -s1 -b2 -e2 -f {allc_sorted_mod_file}"
    
    subprocess.run(submit, shell=True)
    
whitelist = pd.read_csv(f"{indir}/{sample}/{sample}_whitelist.csv")
all_bcs = sorted(whitelist.tenx_whitelist)
args = [(indir, sample, bc) for bc in all_bcs]

pool = Pool(int(cores))
results = pool.starmap(allc_process, args)
pool.close()
pool.join()