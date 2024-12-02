import argparse
from multiprocessing import Pool
import subprocess
import pandas as pd
import json

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=str)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-p", "--parts_batch", type=int)

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
parts_batch = args.parts_batch

def allc_process(indir, sample, bc):
    allcools_dir = f"{indir}/{sample}/allcools"

    cg_file = f"{allcools_dir}/{bc}_CpG_context.tsv"
    ch_file = f"{allcools_dir}/{bc}_Non_CpG_context.tsv"
    allc_file = f"{allcools_dir}/{bc}_allc.tsv"
    allc_sorted_file = f"{allcools_dir}/{bc}_allc_sorted.tsv"
    
    submit = f"cat {cg_file} {ch_file} > {allc_file}"
    subprocess.run(submit, shell=True)

    submit = f"sort -k1,1 -k2,2n {allc_file} > {allc_sorted_file}"
    subprocess.run(submit, shell=True)

    submit = f"bgzip -f {allc_sorted_file}"
    subprocess.run(submit, shell=True)

    submit = f"tabix -s1 -b2 -e2 -f {allc_sorted_file}.gz"
    subprocess.run(submit, shell=True)
    
def allc_stats_process(indir, sample, bc):
    
    allcools_dir = f"{indir}/{sample}/allcools"
    allc_sorted_file = f"{allcools_dir}/{bc}_allc_sorted.tsv.gz"
    cc = pd.read_table(allc_sorted_file,header=None)
    
    all_cnt = pd.crosstab(cc[3],cc[4])
    all_cnt['total_c'] = all_cnt.sum(axis=1)
    all_cnt['ratio'] = all_cnt[1]/all_cnt['total_c']*100
    c_calls = [bc]
    c_calls.extend(all_cnt.melt().value.to_list())
    
    return c_calls
    

multiplier_per_task = 4
    
start = (parts_batch - 1) * multiplier_per_task
end = parts_batch * multiplier_per_task

with open(f"{indir}/{sample}/{sample}_whitelist_batches.json", "r") as file:
    bc_splits = json.load(file)
    
for b in range(start, end):
    
    print(b)
    args = [(indir, sample, bc) for bc in bc_splits[b]]
    
    print(args)
    
    pool = Pool(int(cores))
    results = pool.starmap(allc_process, args)
    pool.close()
    pool.join()

    
    
    

