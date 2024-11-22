import argparse
import methyl_utils

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-l", "--limit", default=False, action="store_true")
parser.add_argument("-p", "--parts_batch", type=int)
args = parser.parse_args()

indir = args.indir
sample = args.sample
limit = args.limit
parts_batch = args.parts_batch

multiplier_per_task = 2

start = (parts_batch - 1) * multiplier_per_task
end = parts_batch * multiplier_per_task


for part in range(start, end):
    print(part)
    methyl_utils.save_quad_batch_from_bam(indir, sample, str(part).zfill(3), limit)