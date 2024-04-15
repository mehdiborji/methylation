import argparse
from multiprocessing import Pool
import methyl_utils
import os
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cores', type=str)
parser.add_argument('-i', '--indir', type=str)
parser.add_argument('-s', '--sample', type=str)
parser.add_argument('-w', '--window_size', type=int)
parser.add_argument('-l', '--limit', default=False, action='store_true')

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
window_size = args.window_size
limit = args.limit

######################################################


######################################################

parts = methyl_utils.find_sub_fastq_parts(indir,sample)
args = [(indir, sample, part, limit) for part in parts]

"""
pool = Pool(int(cores))
results = pool.starmap(methyl_utils.extract_bc_from_bam, args)
pool.close()
pool.join()

methyl_utils.aggregate_bc_dicts(indir,sample)


pool = Pool(int(cores))
results = pool.starmap(methyl_utils.save_quad_batch_json, args)
pool.close()
pool.join()
"""

bcs = pd.read_csv(f'{indir}/{sample}/{sample}_whitelist.csv')
sub_batch_N = int(np.sqrt(len(bcs)))+1
args = [(indir, sample, str(j+1).zfill(3)) for j in range(sub_batch_N)]

"""
pool = Pool(10)
results = pool.starmap(methyl_utils.aggregate_quad_parts, args)
pool.close()
pool.join()
"""

chr_idx_dict = methyl_utils.chrom_to_windows(methyl_utils.hg38_chroms, window_size)

print('chr_idx_dict lenght = ',len(chr_idx_dict))

bcs = pd.read_csv(f'{indir}/{sample}/{sample}_whitelist.csv')
sub_batch_N = int(np.sqrt(len(bcs)))+1
args = [(indir, sample, str(j+1).zfill(3), window_size, chr_idx_dict) for j in range(sub_batch_N)]

pool = Pool(int(cores))
results = pool.starmap(methyl_utils.make_count_sparse_mtx_batch_windows, args)
pool.close()
pool.join()
