import argparse
from multiprocessing import Pool
import methyl_utils

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=int)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-w", "--window_size", type=int)
parser.add_argument("-m", "--methylation_context", type=str)
parser.add_argument("-r", "--reference_genome_index", type=str)
args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
window_size = args.window_size
context = args.methylation_context
reference_genome_index = args.reference_genome_index

chr_idx_dict = methyl_utils.fasta_index_to_windows(reference_genome_index, window_size)
print("chr_idx_dict lenght = ", len(chr_idx_dict))

methyl_utils.stack_mtx_windows(indir, sample, window_size, chr_idx_dict, context, cores)