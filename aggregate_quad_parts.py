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

multiplier_per_task = 5

start = (parts_batch - 1) * multiplier_per_task
end = parts_batch * multiplier_per_task

for b in range(start, end):
    batch = str(b + 1).zfill(3)
    print(batch)
    methyl_utils.aggregate_quad_parts(indir, sample, batch, methylation_context)
