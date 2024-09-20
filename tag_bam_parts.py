import argparse
from multiprocessing import Pool
import methyl_utils

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=str)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-g", "--genomic", default=False, action="store_true")
parser.add_argument("-l", "--limit", default=False, action="store_true")

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
genomic = args.genomic
limit = args.limit


parts = methyl_utils.find_sub_fastq_parts(indir, sample)
args = [(indir, sample, part, limit) for part in parts]
pool = Pool(int(cores))

if genomic:
    results = pool.starmap(methyl_utils.tag_bam_with_all_barcodes, args)
else:
    results = pool.starmap(methyl_utils.tag_bam_with_barcodes, args)

pool.close()
pool.join()
