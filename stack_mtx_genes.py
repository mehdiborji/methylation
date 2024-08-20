import argparse
import methyl_utils
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=int)
parser.add_argument("-i", "--indir", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-m", "--methylation_context", type=str)
parser.add_argument("-r", "--reference_gencode", type=str)

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
context = args.methylation_context
reference_gencode = args.reference_gencode

df_gtf_genes = pd.read_csv(reference_gencode)
gene_idx_list = df_gtf_genes["gene_id"].tolist()

methyl_utils.stack_mtx_genes(indir, sample, gene_idx_list, context, cores)
