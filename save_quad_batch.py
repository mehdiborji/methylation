import argparse
from multiprocessing import Pool
import methyl_utils

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=int)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-l", "--limit", default=False, action="store_true")
parser.add_argument("-p", "--parts_batch", type=int)
args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
limit = args.limit
parts_batch = args.parts_batch

parts = methyl_utils.find_sub_fastq_parts(indir, sample)

start = (parts_batch - 1) * cores
end = parts_batch * cores

args = [(indir, sample, part, limit) for part in parts[start:end]]
pool = Pool(cores)
results = pool.starmap(methyl_utils.save_quad_batch_from_bam, args)
pool.close()
pool.join()
