import argparse
import methyl_utils

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-w", "--window_size", type=int)
parser.add_argument("-m", "--methylation_context", type=str)
parser.add_argument("-r", "--reference_genome_index", type=str)
parser.add_argument("-p", "--parts_batch", type=int)
args = parser.parse_args()

indir = args.indir
sample = args.sample
window_size = args.window_size
methylation_context = args.methylation_context
reference_genome_index = args.reference_genome_index
parts_batch = args.parts_batch

chr_idx_dict = methyl_utils.fasta_index_to_windows(reference_genome_index, window_size)

print("chr_idx_dict lenght = ", len(chr_idx_dict))

if methylation_context == "Non_CpG_context":
    multiplier_per_task = 1  # will run this many of batches in each task, could be one by one or many ar once using pool
else:
    multiplier_per_task = 10

start = (parts_batch - 1) * multiplier_per_task
end = parts_batch * multiplier_per_task

for b in range(start, end):
    batch = str(b + 1).zfill(3)
    print(batch)
    methyl_utils.make_count_sparse_mtx_batch_windows(
        indir, sample, batch, window_size, chr_idx_dict, methylation_context
    )
