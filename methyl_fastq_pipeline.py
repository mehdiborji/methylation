import argparse
from multiprocessing import Pool
import methyl_utils

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=str)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-l", "--limit", default=False, action="store_true")

args = parser.parse_args()
cores = args.cores
indir = args.indir
sample = args.sample
limit = args.limit

######################################################

methyl_utils.split_fastq_by_lines(indir, sample, 20e6)

######################################################

parts = methyl_utils.find_sub_fastq_parts(indir, sample)
args = [(indir, sample, part, limit) for part in parts]

pool = Pool(int(cores))
results = pool.starmap(methyl_utils.extract_clean_fastq, args)
pool.close()
pool.join()
"""
methyl_utils.aggregate_bc_dicts(indir,sample)

pool = Pool(int(cores))
results = pool.starmap(methyl_utils.extract_bc_from_bam, args)
pool.close()
pool.join()

"""
