import argparse
from multiprocessing import Pool
import methyl_utils
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=str)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-b", "--barcodes", type=str)
parser.add_argument("-l", "--limit", default=False, action="store_true")

args = parser.parse_args()
cores = args.cores
indir = args.indir
sample = args.sample
barcodes = args.barcodes
limit = args.limit


py_dir = os.path.dirname(os.path.abspath(__file__))

######################################################

methyl_utils.split_fastq_by_lines(indir, sample, 20e6)

######################################################

parts = methyl_utils.find_sub_fastq_parts(indir, sample)
args = [(indir, sample, part, limit) for part in parts]

pool = Pool(int(cores))
results = pool.starmap(methyl_utils.extract_clean_fastq, args)
pool.close()
pool.join()
######################################################

methyl_utils.aggregate_bc_dicts(indir,sample)

methyl_utils.write_bc_raw_reads(indir, sample, 20000)

######################################################
if barcodes is None:
    whitelist = f'{py_dir}/data/737K-arc-v1.txt.gz'

methyl_utils.write_bc_whitelist(indir, sample, whitelist)

subprocess.call([ f'{py_dir}/scripts/barcode_ref.sh', f'{indir}/{sample}/{sample}_bc_whitelist.fasta', f'{indir}/{sample}/{sample}_ref/'])

subprocess.call([ f'{py_dir}/scripts/barcode_align.sh', f'{indir}/{sample}/{sample}_bc_raw_reads.fasta', f'{indir}/{sample}/{sample}_ref/', f'{indir}/{sample}/{sample}_matching', cores, '-1'])
######################################################
    
    
"""
pool = Pool(int(cores))
results = pool.starmap(methyl_utils.extract_bc_from_bam, args)
pool.close()
pool.join()

"""
