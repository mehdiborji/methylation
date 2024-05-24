import argparse
from multiprocessing import Pool
import methyl_utils

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=str)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-w", "--window_size", type=int)
parser.add_argument("-m", "--methylation_context", type=str)
parser.add_argument("-r", "--reference_genome_index", type=str)
parser.add_argument("-l", "--limit", default=False, action="store_true")
args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
window_size = args.window_size
methylation_context = args.methylation_context
reference_genome_index = args.reference_genome_index
limit = args.limit

sub_batch_N = 20
parts = methyl_utils.find_sub_fastq_parts(indir, sample)

######################################################

args = [(indir, sample, part, methylation_context, limit) for part in parts]
pool = Pool(int(cores))
results = pool.starmap(methyl_utils.save_quad_batch_json, args)
pool.close()
pool.join()

######################################################

args = [(indir, sample, str(j+1).zfill(3), methylation_context) for j in range(sub_batch_N)]

pool = Pool(int(cores))
results = pool.starmap(methyl_utils.aggregate_quad_parts, args)
pool.close()
pool.join()

######################################################

chr_idx_dict = methyl_utils.fasta_index_to_windows(reference_genome_index, window_size)

print("chr_idx_dict lenght = ", len(chr_idx_dict))

args = [
    (indir, sample, str(j + 1).zfill(3), window_size, chr_idx_dict, methylation_context)
    for j in range(sub_batch_N)
]

pool = Pool(int(cores))
results = pool.starmap(methyl_utils.make_count_sparse_mtx_batch_windows, args)
pool.close()
pool.join()
######################################################

parts = methyl_utils.find_sub_fastq_parts(indir, sample)
args = [(indir, sample, part, limit) for part in parts]
pool = Pool(int(cores))
results = pool.starmap(methyl_utils.tag_bam_with_barcodes, args)
pool.close()
pool.join()

######################################################
