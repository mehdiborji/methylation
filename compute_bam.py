import argparse
from multiprocessing import Pool
import methyl_utils
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=str)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-r", "--reference_genome_index", type=str)
parser.add_argument("-l", "--limit", default=False, action="store_true")
args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
reference_genome_index = args.reference_genome_index
limit = args.limit

chrs = pd.read_table(reference_genome_index, header=None)

chrs = chrs[chrs[0].str.startswith("chr")][0].tolist()

args = [(indir, sample, c, limit) for c in chrs]
pool = Pool(int(cores))
results = pool.starmap(methyl_utils.compute_dup_rate, args)
pool.close()
pool.join()
