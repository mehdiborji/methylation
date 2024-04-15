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
parser.add_argument('-w', '--window_size', type=str)
parser.add_argument('-l', '--limit', default=False, action='store_true')

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
window_size = args.window_size
limit = args.limit

hg38_chroms = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}

chr_idx_dict = methyl_utils.chrom_to_windows(hg38_chroms, window_size)

print('chr_idx_dict lenght = ',len(chr_idx_dict))

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

pool = Pool(int(cores))
#results = pool.starmap(methyl_utils.make_count_sparse_mtx_batch, args)
results = pool.starmap(methyl_utils.make_count_sparse_mtx_batch_windows, args[:1])
pool.close()
pool.join()