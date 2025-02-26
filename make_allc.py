import argparse
import methyl_utils

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-m", "--methylation_context", type=str)
parser.add_argument("-p", "--parts_batch", type=int)
args = parser.parse_args()

indir = args.indir
sample = args.sample
methylation_context = args.methylation_context
parts_batch = args.parts_batch

# will run this many of batches in each task
# could be one by one or many ar once using pool

if methylation_context == "Non_CpG_context":
    multiplier_per_task = 2
else:
    multiplier_per_task = 20

start = (parts_batch - 1) * multiplier_per_task
end = parts_batch * multiplier_per_task

for b in range(start, end):
    batch = str(b + 1).zfill(3)
    print(batch)
    methyl_utils.make_allc_tsv(indir, sample, batch, methylation_context)
