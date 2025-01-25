import pysam
import csv
import os
import numpy as np
import pandas as pd
import edlib
import json
import subprocess
import gzip
import re
import mappy
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io
from scipy.sparse import csr_matrix
from scipy.sparse import vstack
import concurrent.futures
from anndata import AnnData
import scanpy as sc
import time
import pybedtools

N_read_extract = 3e5  # maximum reads for limited moded in testing
N_interval_log = 1e6  # interval for logging

print(N_read_extract)

def seq_counter(seq_dict, seq_instance):
    if seq_instance in seq_dict:
        seq_dict[seq_instance] += 1
    else:
        seq_dict[seq_instance] = 1


def quad_dict_store(quad_dict, quad_key, quad_items):
    if quad_key in quad_dict:
        quad_dict[quad_key].extend([quad_items])
    else:
        quad_dict[quad_key] = [quad_items]


def quality_calc(seq, quals, bases_dict, quals_dict):
    for i in range(len(seq)):
        if i in bases_dict:
            seq_counter(bases_dict[i], seq[i])
        else:
            bases_dict[i] = {}
            seq_counter(bases_dict[i], seq[i])

        if i in quals_dict:
            seq_counter(quals_dict[i], quals[i])
        else:
            quals_dict[i] = {}
            seq_counter(quals_dict[i], quals[i])


def string_position_count(seq, seqs_dict):
    for i in range(len(seq)):
        if i in seqs_dict:
            seq_counter(seqs_dict[i], seq[i])
        else:
            seqs_dict[i] = {}
            seq_counter(seqs_dict[i], seq[i])


def quality_df(quals_dict):
    quals_df = pd.DataFrame(quals_dict)
    quals_df = quals_df.T
    quals_df = quals_df.fillna(0)
    quals_df = quals_df.stack()
    quals_df = quals_df.reset_index()
    # quals_df.columns = ['base', 'quality', 'tot_count']
    # quals_df['mult'] = quals_df.quality * quals_df.tot_count
    # quals_df_grouped = quals_df.groupby('base').sum()
    quals_df.columns = ["position", "quantity", "total_cnt"]
    quals_df.position = quals_df.position.astype("int")
    # quals_df[quals_df.position.isin(np.arange(10))]
    counts_df = quals_df.groupby("position").sum()
    quals_df["position_cnt"] = quals_df.position.apply(
        lambda x: counts_df.loc[x].total_cnt
    )
    quals_df["frequency"] = quals_df.total_cnt / quals_df.position_cnt * 100
    return quals_df


def edit_match(input_seq, target_seq, max_dist):
    if input_seq == target_seq:
        dist = 0
        match = True
    else:
        edit = edlib.align(input_seq, target_seq, "NW", "path", max_dist)
        dist = edit["editDistance"]
        if dist >= 0 and dist <= max_dist:
            cigar = edit["cigar"]
            if "D" in cigar or "I" in cigar:
                match = False
                dist = "indel"
            else:
                match = True
        else:
            match = False

    return (match, dist)


def find_sub_fastq_parts(indir, sample):
    pattern = re.compile(r"_R1.part_(.*?)\.fastq")
    all_files = os.listdir(f"{indir}/{sample}/split/")
    parts = sorted(
        [
            f.split(".part_")[1].split(".fastq")[0]
            for f in all_files
            if pattern.search(f)
        ]
    )
    parts = sorted(
        np.unique([f.split(".part_")[1][:3] for f in all_files if pattern.search(f)])
    )
    # part + 3 digits because we did split suffix with 3 digits

    return parts

def run_command(command):
    subprocess.run(command, shell=True)

def split_fastq_by_lines(indir, sample, lines=4e6):
    splitted_file = f"{indir}/{sample}/split/{sample}_R1.part_000.fastq"

    if os.path.isfile(splitted_file):
        print(splitted_file, " splitted fastq exists, skip splitting")
    else:
        print(splitted_file, " splitted fastq does not exist")

        R1 = f"{indir}/{sample}_R1_001.fastq.gz"

        split_dir = f"{indir}/{sample}/split"
        if not os.path.exists(split_dir):
            os.makedirs(split_dir)
            print(f"{split_dir} created")
        else:
            print(f"{split_dir} already exists")

        split_R1_name = f"{split_dir}/{sample}_R1.part_"
        
        # pigz -p 8 -d --stdout {R1} #potential minor speed gain instead of zcat:
        command_R1 = f"zcat {R1} | split -a 3 -l {int(lines)} -d --additional-suffix=.fastq - {split_R1_name}"
        command_R2 = command_R1.replace("_R1", "_R2")
        command_R3 = command_R1.replace("_R1", "_R3")
        
        from concurrent.futures import ThreadPoolExecutor
        
        commands = [command_R1, command_R2, command_R3]
        with ThreadPoolExecutor() as executor:
            executor.map(run_command, commands)

def extract_clean_fastq(indir, sample, part, limit):
    
    total_reads = 0
    
    start_time = time.time()
    
    R1_fastq = f"{indir}/{sample}/split/{sample}_R1.part_{part}.fastq"
    R2_fastq = R1_fastq.replace("_R1.", "_R2.")
    R3_fastq = R1_fastq.replace("_R1.", "_R3.")

    R1_fastq_clean = f"{indir}/{sample}/split/{sample}_R1.part_{part}_clean.fastq"
    R3_fastq_clean = R1_fastq_clean.replace("_R1.", "_R3.")

    bcs_json = f"{indir}/{sample}/split/{sample}.part_{part}_bcs.json"

    if os.path.isfile(bcs_json):
        print(bcs_json, " exists, skip")
        return

    bcs_dict = {}

    R1_clean = open(R1_fastq_clean, "w")
    R3_clean = open(R3_fastq_clean, "w")

    r1_qual_dict = {}
    r2_qual_dict = {}
    r3_qual_dict = {}
    r1_base_dict = {}
    r2_base_dict = {}
    r3_base_dict = {}

    # store_bases = False
    # store_quals = True
    do_qc = True

    with pysam.FastxFile(R1_fastq) as R1, pysam.FastxFile(
        R2_fastq
    ) as R2, pysam.FastxFile(R3_fastq) as R3:
        for r1, r2, r3 in zip(R1, R2, R3):
            
            total_reads += 1

            seq1 = r1.sequence
            seq2 = r2.sequence
            seq3 = r3.sequence

            len1 = len(seq1)
            len3 = len(seq3)

            # bc = seq2 for sciMET data, matched already
            # seq_counter(bcs_dict,bc)

            if do_qc and total_reads % 500 == 0:
                quals1 = r1.get_quality_array()
                quality_calc(seq1, quals1, r1_base_dict, r1_qual_dict)

                quals2 = r2.get_quality_array()
                quality_calc(seq2, quals2, r2_base_dict, r2_qual_dict)

                quals3 = r3.get_quality_array()
                quality_calc(seq3, quals3, r3_base_dict, r3_qual_dict)

            if len1 >= 40 and len3 >= 40:
                
                bc = seq2[8:24]

                match, dist = edit_match(seq2[:8], "CAGACGCG", 2)

                if match:
                    
                    seq_counter(bcs_dict, bc)

                    R1_clean.write(f"@{r1.name}_{bc} r1\n")
                    R1_clean.write(f"{r1.sequence}\n")
                    R1_clean.write("+\n")
                    R1_clean.write(f"{r1.quality}\n")

                    R3_clean.write(f"@{r3.name}_{bc} r2\n")
                    R3_clean.write(f"{r3.sequence}\n")
                    R3_clean.write("+\n")
                    R3_clean.write(f"{r3.quality}\n")
            
            if total_reads % N_interval_log == 0:
                elapsed_time = time.time() - start_time
                message = f"Processed {total_reads} reads in {elapsed_time:.2f}s"
                print(f"{message} in part {part}", flush=True)
            
            if total_reads > N_read_extract and limit:
                break

    r1_qual_df = quality_df(r1_qual_dict)
    r2_qual_df = quality_df(r2_qual_dict)
    r3_qual_df = quality_df(r3_qual_dict)

    r1_base_df = quality_df(r1_base_dict)
    r2_base_df = quality_df(r2_base_dict)
    r3_base_df = quality_df(r3_base_dict)

    r1_qual_df.to_csv(R1_fastq.replace(".fastq", "_quals.csv"))
    r2_qual_df.to_csv(R2_fastq.replace(".fastq", "_quals.csv"))
    r3_qual_df.to_csv(R3_fastq.replace(".fastq", "_quals.csv"))

    r1_base_df.to_csv(R1_fastq.replace(".fastq", "_bases.csv"))
    r2_base_df.to_csv(R2_fastq.replace(".fastq", "_bases.csv"))
    r3_base_df.to_csv(R3_fastq.replace(".fastq", "_bases.csv"))

    R1_clean.close()
    R3_clean.close()

    with open(bcs_json, "w") as json_file:
        json.dump(bcs_dict, json_file)


def extract_bc_from_bam(indir, sample, part, limit):
    i = 0

    bam = f"{indir}/{sample}/split/output_{part}/{sample}_R1.part_{part}_clean_bismark_bt2_pe.bam"

    bcs_json = f"{indir}/{sample}/split/{sample}.part_{part}_bcs.json"

    if os.path.isfile(bcs_json):
        print(bcs_json, " exists, skip")
        # return

    bcs_dict = {}

    samfile = pysam.AlignmentFile(bam, "r")

    for read in samfile.fetch(until_eof=True):
        i += 1

        bc = read.qname.split("_")[-2]

        seq_counter(bcs_dict, bc)

        if i > N_read_extract and limit:
            break

    with open(bcs_json, "w") as json_file:
        json.dump(bcs_dict, json_file)


def tag_bismark_bam_with_whitelist_barcodes(indir, sample, part, limit):
    
    matching_csv = pd.read_csv(f"{indir}/{sample}/{sample}_raw_to_tenx_whitelist.csv")
    raw_to_tenx = dict(zip(matching_csv.bc, matching_csv.tenx_whitelist))

    sam_tag = f"{indir}/{sample}/split/{sample}.part_{part}_tagged.bam"
    sam = f"{indir}/{sample}/split/output_{part}/{sample}_R1.part_{part}_clean_trim_bismark_bt2_pe.bam"

    #if os.path.isfile(sam_tag):
    #    print(sam_tag, " exists, skip")
    #    return

    samfile = pysam.AlignmentFile(sam, "rb")
    tagged_bam = pysam.AlignmentFile(sam_tag, "wb", template=samfile)

    total_reads = 0
    
    start_time = time.time()
    
    for read in samfile.fetch(until_eof=True):
        
        total_reads += 1
        raw_barcode = read.qname.split("_")[-2]
        if raw_barcode in raw_to_tenx:
            matched_barcode = raw_to_tenx[raw_barcode]

            read.set_tag("CB", matched_barcode)
            tagged_bam.write(read)
            
        if total_reads % N_interval_log == 0:
                elapsed_time = time.time() - start_time
                print(f"Processed {total_reads} reads of part {part} in {elapsed_time:.2f}s", flush=True)
            
        if total_reads > N_read_extract and limit:
            break

    tagged_bam.close()
    samfile.close()
    
def tag_minimap_bam_with_all_barcodes(indir, sample, part, limit):
    
    matching_csv = pd.read_csv(f"{indir}/{sample}/{sample}_matching_raw_to_tenx_passing.csv")
    raw_to_tenx = dict(zip(matching_csv.bc, matching_csv.tenx_whitelist))

    sam_tag = f"{indir}/{sample}/split/{sample}.part_{part}_tagged.bam"
    sam = f"{indir}/{sample}/split/{sample}_part_{part}.bam"

    if os.path.isfile(sam_tag):
        print(sam_tag, " exists, skip")
        return

    samfile = pysam.AlignmentFile(sam, "rb")
    tagged_bam = pysam.AlignmentFile(sam_tag, "wb", template=samfile)

    total_reads = 0
    
    start_time = time.time()
    
    for read in samfile.fetch(until_eof=True):
        total_reads += 1
        
        if read.is_proper_pair and read.mapq > 5:
            
            raw_barcode = read.qname.split("_")[-1]

            if raw_barcode in raw_to_tenx:
                matched_barcode = raw_to_tenx[raw_barcode]

                read.set_tag("CB", matched_barcode)
                tagged_bam.write(read)

        if total_reads % N_interval_log == 0:
                elapsed_time = time.time() - start_time
                print(f"Processed {total_reads} reads of part {part} in {elapsed_time:.2f}s", flush=True)
            
        if total_reads > N_read_extract and limit:
            break
    tagged_bam.close()
    samfile.close()


def compute_dup_rate(indir, sample, chrom, limit):
    N_interval_log = 1e6
    
    sam = f"{indir}/{sample}/{sample}_markdup.bam"
    
    if not os.path.isfile(sam):
        sam = f"{indir}/{sample}/{sample}_markdup_piped.bam"
        
    samfile = pysam.AlignmentFile(sam, "rb")

    BC_dup_count = {}
    BC_unique_count = {}

    start_time = time.time()

    total_reads = 0

    for read in samfile.fetch(chrom):
        total_reads += 1
        bc = read.get_tag("CB")
        if read.is_duplicate:
            seq_counter(BC_dup_count, bc)
        else:
            seq_counter(BC_unique_count, bc)

        if total_reads % N_interval_log == 0:
            elapsed_time = time.time() - start_time
            print(
                f"Processed {total_reads} reads of {chrom} in {elapsed_time:.2f}s",
                flush=True,
            )
        if total_reads > N_read_extract and limit:
            break
    samfile.close()

    BC_dup_df = pd.Series(BC_dup_count)
    BC_uni_df = pd.Series(BC_unique_count)
    BC_dup_df.name = "dup_cnt"
    BC_uni_df.name = "uniq_cnt"
    mrg = pd.merge(BC_uni_df, BC_dup_df, how="outer", left_index=True, right_index=True)
    mrg["total_cnt"] = mrg.sum(axis=1)
    mrg["dup_rate"] = mrg.total_cnt / mrg.uniq_cnt
    mrg["log10cnt"] = np.log10(mrg.total_cnt)
    mrg["log10cnt_uniq"] = np.log10(mrg.uniq_cnt)
    mrg.to_csv(f"{indir}/{sample}/{sample}_dup_rate_{chrom}.csv")


def aggregate_bc_dicts(indir, sample):
    
    dir_split = f"{indir}/{sample}/split"
    files = os.listdir(dir_split)
    jsons = sorted([f for f in files if "_bcs.json" in f])

    agg_read_csv = f"{indir}/{sample}/{sample}_agg_cnt_raw_bcs.csv"

    if os.path.isfile(agg_read_csv):
        print(agg_read_csv, " exists, skip")
        return

    data_agg = {}

    for jsn in jsons:
        
        jsn_to_load = f"{dir_split}/{jsn}"
        
        with open(jsn_to_load, "r") as json_file:
            data_sub = json.load(json_file)
            for k in data_sub:
                if data_agg.get(k) is not None:
                    data_agg[k] += data_sub[k]
                else:
                    data_agg[k] = data_sub[k]

    pd.Series(data_agg).to_csv(agg_read_csv)


def write_bc_whitelist(indir, sample, bc_file):
    
    fasta_file = f"{indir}/{sample}/{sample}_bc_whitelist.fasta"

    if os.path.isfile(fasta_file):
        print(fasta_file, " exists, skip")
        return

    bcs = pd.read_table(bc_file, names=["bc"])
    bcs = pd.DataFrame(bcs.bc.apply(lambda x: x.split("-")[0]))

    bcs = bcs.sort_values(by="bc")
    bcs["rev_bc"] = bcs.bc.apply(lambda x: mappy.revcomp(x))

    with open(fasta_file, "w") as f:
        for bc in bcs.rev_bc:
            f.write(f">{bc}\n")
            f.write(f"{bc}\n")

            
def find_local_min(df):
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks

    bins = 100

    cnt_sorted = df.sort_values(ascending=False)
    sub = cnt_sorted.iloc[100:].copy()

    x = np.histogram(sub, bins)  

    smooth = gaussian_filter1d(x[0], 3)
    # find the local minimum
    peak_idx, _ = find_peaks(-smooth)
    # take the mid point of point before and after
    print(peak_idx, x[1][:-1][peak_idx])
    mean_hist = (x[1][1:][peak_idx] + x[1][:-1][peak_idx]) / 2  

    # take the last value in list of local minima (could be more than one)
    mean_hist = mean_hist[-1]
    
    return mean_hist


def DNA_RNA_exact_shared(indir, sample, dna_10x_wl_file, rna_10x_wl_file, max_expected_barcodes):
    
    atac = pd.read_table(dna_10x_wl_file, names=["atac"])
    rna = pd.read_table(rna_10x_wl_file, names=["rna"])
    all_bcs = atac.join(rna)
    all_bcs['rev_comp_atac'] = all_bcs.atac.apply(lambda x: mappy.revcomp(x))

    
    N_interval_log = 2e6
    
    raw_bcs_DNA_csv = f"{indir}/{sample}/{sample}_agg_cnt_raw_bcs.csv"
    rev_comp_atac_dict = dict.fromkeys(all_bcs['rev_comp_atac'], None)

    with open(raw_bcs_DNA_csv, "r") as f:
        reader = csv.reader(f, delimiter=",")
        for i, line in enumerate(reader):
            
            if line[0] in rev_comp_atac_dict:
                rev_comp_atac_dict[line[0]] = int(line[1])
            if i % N_interval_log == 0:
                print(i, line)
            #if i>10000000:
            #    break

    raw_bcs_RNA_csv = f"{indir}/{sample}/{sample}_agg_cnt_raw_bcs_RNA.csv" 
    rna_dict = dict.fromkeys(all_bcs['rna'], None)

    with open(raw_bcs_RNA_csv, "r") as f:
        reader = csv.reader(f, delimiter=",")
        for i, line in enumerate(reader):
            
            if line[0] in rna_dict:
                rna_dict[line[0]] = int(line[1])
            if i % N_interval_log == 0:
                print(i, line)
            #if i>10000000:
            #    break
        
    all_bcs['RNA_cnt'] = all_bcs.rna.apply(lambda x: rna_dict[x])
    all_bcs['DNA_cnt'] = all_bcs.rev_comp_atac.apply(lambda x: rev_comp_atac_dict[x])

    filt_bcs = all_bcs[ (all_bcs.RNA_cnt>10) & (all_bcs.DNA_cnt>10) ].copy()

    top_RNA = filt_bcs['RNA_cnt'].nlargest(max_expected_barcodes)
    top_DNA = filt_bcs['DNA_cnt'].nlargest(max_expected_barcodes)

    #filt_bcs = filt_bcs[filt_bcs['RNA_cnt'].isin(top_RNA) & filt_bcs['DNA_cnt'].isin(top_DNA)].copy()
    
    filt_bcs = filt_bcs[ (filt_bcs.RNA_cnt>=top_RNA.min()) & (filt_bcs.DNA_cnt>=top_DNA.min()) ].copy()

    filt_bcs['log10_RNA_cnt'] = np.log10(filt_bcs['RNA_cnt'])
    filt_bcs['log10_DNA_cnt'] = np.log10(filt_bcs['DNA_cnt'])

    filt_bcs.to_csv(f"{indir}/{sample}/{sample}_RNA_DNA_share_bcs.csv")
    
    DNA_threshold = find_local_min(filt_bcs.log10_DNA_cnt)
    RNA_threshold = find_local_min(filt_bcs.log10_RNA_cnt)

    wl_RNA_DNA = filt_bcs[ (filt_bcs.log10_DNA_cnt>DNA_threshold) & (filt_bcs.log10_RNA_cnt>RNA_threshold) ].copy()

    wl_RNA_DNA.to_csv(f"{indir}/{sample}/{sample}_RNA_DNA_whitelist.csv")
    
    
def generate_c_to_t_variants(seq):
    variants = [[None, 0, seq, seq]]

    for i, nucleotide in enumerate(seq):
        if nucleotide == 'C':
            variant = seq[:i] + 'T' + seq[i+1:]
            variants.append([i, 1, seq, variant])
            
    for i, nucleotide in enumerate(seq):
        if nucleotide == 'C':
            for j in range(i + 1, len(seq)):
                if seq[j] == 'C':
                    variant = seq[:i] + 'T' + seq[i+1:j] + 'T' + seq[j+1:]
                    variants.append([i, j, 2, seq, variant])
                    
    for i, nucleotide in enumerate(seq):
        if nucleotide == 'C':
            for j in range(i + 1, len(seq)):
                if seq[j] == 'C':
                    for k in range(j + 1, len(seq)):
                        if seq[k] == 'C':
                            variant = seq[:i] + 'T' + seq[i+1:j] + 'T' + seq[j+1:k] + 'T' + seq[k+1:]
                            variants.append([i, j, k, 3, seq, variant])

    return variants


def c_to_t_correction(indir, sample):
    
    white = pd.read_csv(f"{indir}/{sample}/{sample}_RNA_DNA_whitelist.csv", index_col=0)

    all_hammings = []
    for kk, seq in enumerate(white.rev_comp_atac):
        variants = generate_c_to_t_variants(seq)
        all_hammings.extend([c[-3:] for c in variants])
    
    all_df = pd.DataFrame(all_hammings)
    N_interval_log = 2e6

    raw_bcs_DNA_csv = f"{indir}/{sample}/{sample}_agg_cnt_raw_bcs.csv"
    all_ham_dict = dict.fromkeys(all_df[2], None)

    with open(raw_bcs_DNA_csv, "r") as f:
        reader = csv.reader(f, delimiter=",")

        non_matched = {}
        total_reads = 0
        for i, line in enumerate(reader):

            raw_bc = line[0]
            read_cnt = int(line[1])
            total_reads += read_cnt
            if read_cnt>=10:
                if raw_bc in all_ham_dict:
                    all_ham_dict[raw_bc] = read_cnt
                elif read_cnt>100:
                    non_matched[raw_bc] = read_cnt
            if i % N_interval_log == 0:
                print(i, line)
            if >1000000:
                break
    
    non_matched_csv = f"{indir}/{sample}/{sample}_raw_bcs_non_matched.csv"
    pd.Series(non_matched).sort_values(ascending=False).to_csv(non_matched_csv)
    
    all_df[3] = all_df[2].apply(lambda x: all_ham_dict[x])
    all_df_sub = all_df[all_df[3]>0]
    ham_reps = all_df_sub[2].value_counts()
    ham_reps = ham_reps[ham_reps>1].index

    all_df_reps = all_df_sub[all_df_sub[2].isin(ham_reps)]
    
def write_bc_raw_reads(indir, sample, threshold):
    fasta_file = f"{indir}/{sample}/{sample}_bc_raw_reads.fasta"

    if os.path.isfile(fasta_file):
        print(fasta_file, " exists, skip")
        return

    bc_file = f"{indir}/{sample}/{sample}_agg_cnt_raw_bcs.csv"
    bcs = pd.read_csv(bc_file)
    bcs.columns = ["bc", "read_cnt"]
    
    bcs = bcs[bcs.read_cnt > threshold].copy()
    bcs = bcs.sort_values(by="bc", ascending=False)
    with open(fasta_file, "w") as f:
        for bc in bcs.bc:
            f.write(f">{bc}\n")
            f.write(f"{bc}\n")


def processing_matching(indir, sample, AS_min=12):
    all_AS = []
    all_pairs = []

    samfile = pysam.AlignmentFile(f"{indir}/{sample}/{sample}_matching.sam", "rb")
    
    matched_csv = f"{indir}/{sample}/{sample}_matching_raw_to_tenx_passing.csv" 
    if os.path.isfile(matched_csv):
        print(matched_csv, " exists, skip")
        return
    
    for read in samfile.fetch():
        AS = read.get_tag("AS")
        all_AS.append([AS, read.flag])
        if read.flag == 0:
            bc = read.reference_name
            seq = read.query
            all_pairs.append([AS, bc, seq])

    print("making all_pairs DF")
    all_pairs = pd.DataFrame(all_pairs)
    bc = pd.read_csv(f"{indir}/{sample}/{sample}_agg_cnt_raw_bcs.csv")
    bc.columns = ["bc", "read_cnt"]
    all_pairs.columns = ["AS", "tenx_whitelist", "bc"]
    matched = bc.merge(all_pairs, how="left", on="bc")
    print("saving matching_raw_to_tenx DF")
    matched.to_csv(f"{indir}/{sample}/{sample}_matching_raw_to_tenx.csv", index=None)

    matched = matched[matched.AS >= AS_min].copy()

    print("saving tenx_passin")
    matched.to_csv(matched_csv, index=None
    )

def filered_barcodes(indir, sample, max_expected_barcodes=50000, bins=100):
    
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
    from matplotlib.backends.backend_pdf import PdfPages

    print("opening raw_to_tenx_passing_AS DF")
    matched = pd.read_csv(f"{indir}/{sample}/{sample}_matching_raw_to_tenx_passing.csv")
    matched_sub = matched[["read_cnt", "tenx_whitelist"]]
    matched_sub_grouped = matched_sub.groupby("tenx_whitelist").sum()

    qc_pdf_file = f"{indir}/{sample}/{sample}_QC.pdf"

    #if os.path.isfile(qc_pdf_file):
    #    print(qc_pdf_file, " exists, skip")
    #else:
    qc_pdfs = PdfPages(qc_pdf_file)

    agg_bcs = matched_sub_grouped.sort_values(by="read_cnt", ascending=False)
    agg_bcs["log10_read_cnt"] = np.log10(agg_bcs.read_cnt)

    # select top max_bc except first 100
    sub = agg_bcs.iloc[20:max_expected_barcodes].copy()
    # fit a histogram
    x = np.histogram(sub.log10_read_cnt, bins)  
    # smooth histogram
    smooth = gaussian_filter1d(x[0], 3)
    # find the local minimum
    peak_idx, _ = find_peaks(-smooth)
    # take the mid point of point before and after
    print(peak_idx, x[1][:-1][peak_idx])
    mean_hist = (x[1][1:][peak_idx] + x[1][:-1][peak_idx]) / 2  

    # take the last value in list of local minima (could be more than one)
    mean_hist = mean_hist[-1]  

    print(f"filtering for final whitelist with min log10={round(mean_hist,4)} counts")
    wl_df = agg_bcs[agg_bcs.log10_read_cnt >= mean_hist].copy()
    wl_df.to_csv(f"{indir}/{sample}/{sample}_whitelist.csv")
    # wl_reads=wl_df.read_cnt.sum()

    matched_filtered = matched[matched.tenx_whitelist.isin(wl_df.index)].copy()
    print("saving raw_to_tenx_whitelist")
    matched_filtered.to_csv(
        f"{indir}/{sample}/{sample}_raw_to_tenx_whitelist.csv", index=None
    )

    white_list_size = wl_df.shape[0]

    print("plotting rankplot and histogram of 10x matched read counts")
    plt.figure(figsize=(4, 3))
    log10_ranks = np.log10(np.arange(1, len(agg_bcs) + 1))
    log10_cnts = agg_bcs.log10_read_cnt
    plt.plot(log10_ranks, log10_cnts)  # ,label='Rank Plot of Reads')
    plt.xlabel("Log10 Ranks")
    plt.ylabel("Log10 Read Counts")

    reads_in_wh = round(wl_df.read_cnt.sum()/1e6,1)
    matche_reads = round(agg_bcs.read_cnt.sum()/1e6,1)

    print(reads_in_wh,matche_reads,reads_in_wh/matche_reads*100)
    plt.title(f"{sample}\n {white_list_size} white listed")

    plt.plot(
        [0, log10_ranks[-1]],
        [mean_hist, mean_hist],
        linewidth=1,
        label="log10 threshold",
        c="tab:green",
    )
    log10_wl = np.log10(white_list_size)
    plt.plot(
        [log10_wl, log10_wl],
        [log10_cnts.min(), log10_cnts.max()],
        linewidth=1,
        label="log10 size",
        c="tab:orange",
    )
    plt.legend(loc="best")

    qc_pdfs.savefig(bbox_inches="tight")

    plt.figure(figsize=(4, 3))
    plt.plot(x[1][:-1], x[0], label="Raw Histogram")
    plt.plot(x[1][:-1], smooth, label="Gaussian Smoothed")
    plt.xlabel("Log10 Read Counts")
    plt.ylabel("Bin Height")
    plt.title(f"{sample}\n in_whitelist={reads_in_wh}m \n bc_matched={matche_reads}m")
    plt.plot(
        [mean_hist, mean_hist],
        [0, np.max(x[0])],
        linewidth=2,
        label="Whitelist Threshold",
    )
    plt.legend(loc="best")
    qc_pdfs.savefig(bbox_inches="tight")
    
    
    plt.figure(figsize=(4, 3))
    #plt.hist(sub.log10_read_cnt, 100);
    sns.histplot(sub.log10_read_cnt,bins=bins);
    plt.title(f"{sample}\n")
    qc_pdfs.savefig(bbox_inches="tight")

    qc_pdfs.close()


def split_bcs_to_batches(indir, sample):
    bcs = pd.read_csv(f"{indir}/{sample}/{sample}_whitelist.csv", index_col=0)
    
    sub_batch_N = int(bcs.read_cnt.sum()/1e6/15)
    
    print(f'total batches of 15m reads {sub_batch_N}')
    # Initialize N batches and their sums
    batches = [[] for _ in range(sub_batch_N)]
    batch_sums = [0] * sub_batch_N

    # Distribute numbers to batches based on value, but store original indices
    for index, number in enumerate(bcs["read_cnt"]):
        # Find the batch with the minimum sum
        min_batch_index = batch_sums.index(min(batch_sums))
        # Add the index to this batch
        batches[min_batch_index].append(index)
        # Update the batch sum with the number's value
        batch_sums[min_batch_index] += number

    bc_batches = []
    for b in batches:
        print(len(b), bcs.iloc[b]["read_cnt"].sum())
        bc_batches.append(bcs.iloc[b].index.tolist())

    with open(f"{indir}/{sample}/{sample}_whitelist_batches.json", "w") as file:
        json.dump(bc_batches, file)


def save_quad_batch_json(indir, sample, part, context, limit):
    dir_split = f"{indir}/{sample}/split"
    meth_file = f"{dir_split}/output_{part}/{context}_{sample}_R1.part_{part}_clean_bismark_bt2_pe.txt.gz"

    batch = str(1).zfill(3)
    batch_json = f"{dir_split}/quads_part_{part}_batch_{batch}_{context}.json"

    if os.path.isfile(batch_json):
        print(batch_json, " exists, skip")
        # return

    matching_csv = pd.read_csv(f"{indir}/{sample}/{sample}_raw_to_tenx_whitelist.csv")
    raw_to_tenx = dict(zip(matching_csv.bc, matching_csv.tenx_whitelist))

    with open(f"{indir}/{sample}/{sample}_whitelist_batches.json", "r") as file:
        bc_splits = json.load(file)
    sub_batch_N = len(bc_splits)

    quad_dict = {}

    conversion = False
    if context == "Non_CpG_context":
        conversion = True

    Non_CpG_to_CpG_dict = {"x": "z", "h": "z", "X": "Z", "H": "Z"}

    all_failed_bc = []

    i = 0

    with gzip.open(meth_file, "rt") as f:
        for line in f:
            i += 1

            split_line = line.strip().split("\t")
            # print(split_line)

            if conversion:
                split_line[-1] = Non_CpG_to_CpG_dict[split_line[-1]]

            # print('converted', split_line)
            # raw_bc = split_line[0].split(':')[0] # sciMET fastqs
            raw_bc = split_line[0].split("_")[-1]  # raw bc added to name fastqs

            if raw_bc in raw_to_tenx:
                matched_bc = raw_to_tenx[raw_bc]

                quad_dict_store(quad_dict, matched_bc, "_".join(split_line[-3:]))
            else:
                all_failed_bc.append(raw_bc)

            if i > N_read_extract and limit:
                break

    print(len(all_failed_bc))

    for j in range(sub_batch_N):
        batch = str(j + 1).zfill(3)
        batch_json = f"{dir_split}/quads_part_{part}_batch_{batch}_{context}.json"
        sub_agg = {}
        for a in bc_splits[j]:
            if a in quad_dict:
                sub_agg[a] = quad_dict[a]

        with open(batch_json, "w") as json_file:
            json.dump(sub_agg, json_file)


def save_quad_batch_from_bam(indir, sample, part, limit):
    dir_split = f"{indir}/{sample}/split"
    
    suffix = "clean_trim_bismark_bt2_pe.nonCG_filtered.bam"
    #suffix = "clean_trim_bismark_bt2_pe.bam"
    
    bam_file = f"{dir_split}/output_{part}/{sample}_R1.part_{part}_{suffix}"

    bias_meth_file = (
        f"{dir_split}/output_{part}/{sample}_part_{part}_bias_methylation.csv"
    )

    if os.path.isfile(bias_meth_file):
        print(bias_meth_file, " exists, skip")
        # return

    matching_csv = pd.read_csv(f"{indir}/{sample}/{sample}_raw_to_tenx_whitelist.csv")
    raw_to_tenx = dict(zip(matching_csv.bc, matching_csv.tenx_whitelist))

    with open(f"{indir}/{sample}/{sample}_whitelist_batches.json", "r") as file:
        bc_splits = json.load(file)
    sub_batch_N = len(bc_splits)

    R1_meth_dict = {}
    R2_meth_dict = {}

    context_conversion = {"x": "z", "h": "z", "X": "Z", "H": "Z"}

    quad_dict_CpG = {}
    quad_dict_Non_CpG = {}

    total_failed_reads = 0
    diffs = []
    # frags = {}

    bam = pysam.AlignmentFile(bam_file, "r")

    start_time = time.time()

    total_reads = 0

    for read in bam.fetch(until_eof=True):
        total_reads += 1

        # frag_size = read.template_length
        chrom = read.reference_name
        if read.is_reverse:
            meth = read.get_tag("XM")[::-1]
            chrom_pos = read.get_reference_positions()[::-1]
        else:
            meth = read.get_tag("XM")
            chrom_pos = read.get_reference_positions()

        r1_left_clip = 15
        r1_right_clip = 2
        r2_left_clip = 2
        r2_right_clip = 15
        r1_trim_after = 150
        r2_trim_after = 1000
        
        if read.is_read1:  # mean it's R2
            meth = meth[r2_left_clip:-r2_right_clip][:r2_trim_after]
            chrom_pos = chrom_pos[r2_left_clip:-r2_right_clip][:r2_trim_after]
            string_position_count(meth, R2_meth_dict)
        else:
            meth = meth[r1_left_clip:-r1_right_clip][:r1_trim_after]
            chrom_pos = chrom_pos[r1_left_clip:-r1_right_clip][:r1_trim_after]

            string_position_count(meth, R1_meth_dict)

        if len(chrom_pos) != len(meth):
            diff = len(chrom_pos) - len(meth)
            diffs.append(diff)
            # print()
            # break
            # if diff<-10:
            #    print(read.cigar,diff,len(read.get_tag('XM')),len(read.get_reference_positions()))
            #    break
        else:
            raw_bc = read.qname.split("_")[-2]

            if raw_bc in raw_to_tenx:
                matched_bc = raw_to_tenx[raw_bc]

                for i, char in enumerate(meth):
                    if char in ["H", "X", "x", "h"]:
                        quad_dict_store(
                            quad_dict_Non_CpG,
                            matched_bc,
                            f"{chrom}_{chrom_pos[i]+1}_{context_conversion[char]}",
                        )
                    elif char in ["z", "Z"]:
                        quad_dict_store(
                            quad_dict_CpG,
                            matched_bc,
                            f"{chrom}_{chrom_pos[i]+1}_{char}",
                        )
            else:
                total_failed_reads += 1

        if total_reads % N_interval_log == 0:
            elapsed_time = time.time() - start_time
            print(f"Processed {total_reads} reads in {elapsed_time:.2f} seconds.", flush=True)

        if total_reads > N_read_extract and limit:
            break

    print("total bc fail reads = ", total_failed_reads, "total reads = ", total_reads)

    context = "Non_CpG_context"
    for j in range(sub_batch_N):
        batch = str(j + 1).zfill(3)
        batch_json = f"{dir_split}/quads_part_{part}_batch_{batch}_{context}.json"
        sub_agg = {}
        for a in bc_splits[j]:
            if a in quad_dict_Non_CpG:
                sub_agg[a] = quad_dict_Non_CpG[a]

        with open(batch_json, "w") as json_file:
            json.dump(sub_agg, json_file)

    context = "CpG_context"
    for j in range(sub_batch_N):
        batch = str(j + 1).zfill(3)
        batch_json = f"{dir_split}/quads_part_{part}_batch_{batch}_{context}.json"
        sub_agg = {}
        for a in bc_splits[j]:
            if a in quad_dict_CpG:
                sub_agg[a] = quad_dict_CpG[a]

        with open(batch_json, "w") as json_file:
            json.dump(sub_agg, json_file)

    R1_meth_df = quality_df(R1_meth_dict)
    R2_meth_df = quality_df(R2_meth_dict)

    R1_meth_df.to_csv(bias_meth_file.replace("_bias_", "_R1_"), index=None)
    R2_meth_df.to_csv(bias_meth_file.replace("_bias_", "_R2_"), index=None)

    all_dfs = []
    for r, df in enumerate([R1_meth_df, R2_meth_df]):
        for context in ["X", "Z", "H"]:
            h_count = (
                df[df.quantity.isin([context.lower(), context])]
                .groupby("position")
                .sum()
            )
            h_count_meth = df[df.quantity == context].groupby("position").sum()

            meth_frac = h_count_meth.total_cnt / h_count.total_cnt
            meth_frac = meth_frac.reset_index()

            meth_frac["read"] = f"{r+1}_{context}"

            all_dfs.append(meth_frac)

    all_dfs = pd.concat(all_dfs)
    all_dfs = all_dfs.reset_index()
    all_dfs.total_cnt = all_dfs.total_cnt * 100

    all_dfs.to_csv(bias_meth_file.replace("_R1_", "_bias_"), index=None)


def aggregate_quad_parts(indir, sample, batch, context):
    dir_split = f"{indir}/{sample}/split"
    files = os.listdir(dir_split)
    batch_pattern = f"_batch_{batch}_{context}.json"
    batch_jsons = sorted([f for f in files if batch_pattern in f])

    agg_batch_json_file = f"{dir_split}/quad_agg_{batch}_{context}.json"

    if os.path.isfile(agg_batch_json_file):
        print(agg_batch_json_file, " exists, skip")
        return

    data_agg = {}
    for p in batch_jsons:
        sub_parts_of_batch = f"{dir_split}/{p}"
        with open(sub_parts_of_batch, "r") as json_file:
            data_sub = json.load(json_file)
            # print(p, len(data_sub))
            for k in data_sub:
                if data_agg.get(k) is not None:
                    data_agg[k].extend(data_sub[k])
                else:
                    data_agg[k] = data_sub[k]

    with open(agg_batch_json_file, "w") as json_file:
        json.dump(data_agg, json_file)


def write_mtx(mtx_dir, batch, state, csr):
    if not os.path.exists(mtx_dir):
        try:
            os.makedirs(mtx_dir)
            print(f"{mtx_dir} created")
        except Exception as e:
            print(f"A {e} occurred, {mtx_dir} have already been created")
    else:
        print(f"{mtx_dir} already exists")

    csr_file = f"{mtx_dir}/b_{batch}_{state}.mtx.gz"

    with gzip.open(csr_file, "wb") as out:
        scipy.io.mmwrite(out, csr)


def fasta_index_to_windows(fasta_index, window_size):
    i = 0
    chrs = pd.read_table(fasta_index, header=None)
    chrom_len_dict = dict(zip(chrs[0], chrs[1]))
    chr_idx_dict = {}
    for chrom in chrom_len_dict:
        last_idx = chrom_len_dict[chrom] // window_size
        for chr_id in range(last_idx + 1):
            chr_idx_dict[(chrom, chr_id)] = i
            i += 1
    return chr_idx_dict


def make_count_sparse_mtx_batch_windows(
    indir, sample, batch, window, chr_idx_dict, context
):
    dir_split = f"{indir}/{sample}/split"
    agg_batch_json_file = f"{dir_split}/quad_agg_{batch}_{context}.json"
    # agg_batch_json_file = f"{dir_split}/quads_part_001_batch_{batch}.json"

    window_mtx_dir = f"{indir}/{sample}/counts_w_{window}_m{context}"
    csr_file = f"{window_mtx_dir}/b_{batch}_score.mtx.gz"

    if os.path.isfile(csr_file):
        print(csr_file, " exists, skip")
        return

    start_time = time.time()

    with open(agg_batch_json_file, "r") as json_file:
        data_sub = json.load(json_file)

    elapsed = time.time() - start_time
    print(f"opened {agg_batch_json_file} after {elapsed:.2f}s", flush=True)

    rows_idx = []
    cols_idx = []

    row_col_values_meth = []
    row_col_values_notmeth = []

    rows_idx_score = []
    cols_idx_score = []
    row_col_values_score = []

    batch_bcs = list(data_sub.keys())

    for idx, bc in enumerate(batch_bcs):
        total_bases = len(data_sub[bc])

        elapsed = time.time() - start_time
        print(f"total_bases for {bc} = {total_bases} in {elapsed:.2f}s", flush=True)

        if total_bases < 1000:
            continue

        # convert ['chr_100', 'Z'] to ['chr_100_Z'] and count dedup
        # maybe save this stats later
        # one call also do this in aggregate_quad_parts

        #dedup = {}

        #for b in data_sub[bc]:
        #    seq_counter(dedup, b)

        # convert back to ['chr', '100', 'Z'] format
        #triplets = [d.split("_") for d in dedup.keys()]
        
        dedup = set(data_sub[bc])
        triplets = [d.split("_") for d in dedup]
        triplets = pd.DataFrame(triplets)

        elapsed = time.time() - start_time
        print(
            f"deduped bases for {bc} = {triplets.shape[0]} in {elapsed:.2f}s",
            flush=True,
        )

        triplets[1] = triplets[1].astype("int")
        triplets[3] = (
            triplets[1] // window
        )  # bin base position into size equal to window

        # count Z and z per bin

        Z_cnt = triplets[triplets[2] == "Z"].groupby([0, 3]).size()
        z_cnt = triplets[triplets[2] == "z"].groupby([0, 3]).size()
        Z_cnt.name = "Z_cnt"
        z_cnt.name = "z_cnt"

        mrg = pd.merge(Z_cnt, z_cnt, how="outer", left_index=True, right_index=True)
        mrg = mrg.fillna(0)
        cell_mC = mrg.sum()

        mrg["cov"] = mrg.sum(axis=1)

        cell_mC_ratio = cell_mC.loc["Z_cnt"] / (
            cell_mC.loc["Z_cnt"] + cell_mC.loc["z_cnt"]
        )
        # mrg['meth'] = mrg.Z_cnt / mrg['cov']
        # mrg['ratio'] = mrg['meth'] / cell_mC_ratio
        mrg["diff"] = mrg.Z_cnt / mrg["cov"] - cell_mC_ratio

        diff_pos = mrg[(mrg["diff"] > 0) & (mrg["cov"] > 1)]["diff"] / (
            1 - cell_mC_ratio
        )
        diff_neg = mrg[(mrg["diff"] < 0) & (mrg["cov"] > 1)]["diff"] / cell_mC_ratio

        elapsed = time.time() - start_time
        print(f"dataframes built in {elapsed:.2f}s", flush=True)

        row_vals = (np.ones(len(mrg), dtype=int) * idx).tolist()
        col_vals = [chr_idx_dict[key] for key in mrg.index]

        rows_idx.extend(row_vals)  # cells
        cols_idx.extend(col_vals)  # genes

        row_col_values_meth.extend(mrg["Z_cnt"].tolist())
        row_col_values_notmeth.extend(mrg["z_cnt"].tolist())

        for diff in [diff_pos, diff_neg]:
            row_vals = (np.ones(len(diff), dtype=int) * idx).tolist()
            col_vals = [chr_idx_dict[key] for key in diff.index]
            rows_idx_score.extend(row_vals)  # cells
            cols_idx_score.extend(col_vals)  # genes
            row_col_values_score.extend(diff.tolist())

        elapsed = time.time() - start_time
        print(f"all matrices built in {elapsed:.2f}s", flush=True)

    shape = (len(batch_bcs), len(chr_idx_dict))

    csr = csr_matrix(
        (row_col_values_meth, (rows_idx, cols_idx)), shape=shape, dtype="float32"
    )
    csr.eliminate_zeros()
    write_mtx(window_mtx_dir, batch, "meth", csr)

    csr = csr_matrix(
        (row_col_values_notmeth, (rows_idx, cols_idx)), shape=shape, dtype="float32"
    )
    csr.eliminate_zeros()
    write_mtx(window_mtx_dir, batch, "notmeth", csr)

    csr = csr_matrix(
        (row_col_values_score, (rows_idx_score, cols_idx_score)),
        shape=shape,
        dtype="float32",
    )
    write_mtx(window_mtx_dir, batch, "score", csr)

    elapsed = time.time() - start_time
    print(f"all matrices saved in {elapsed:.2f}s", flush=True)


def make_count_sparse_mtx_batch_genes(indir, sample, batch, gencode, context):
    dir_split = f"{indir}/{sample}/split"
    agg_batch_json_file = f"{dir_split}/quad_agg_{batch}_{context}.json"

    df_gtf_genes = pd.read_csv(gencode)
    bed_genes = pybedtools.BedTool.from_dataframe(df_gtf_genes)

    gene_idx_map = dict(zip(df_gtf_genes["gene_id"], df_gtf_genes.index))

    gene_mtx_dir = f"{indir}/{sample}/counts_gene_m{context}"
    csr_file = f"{gene_mtx_dir}/b_{batch}_score.mtx.gz"

    if os.path.isfile(csr_file):
        print(csr_file, " exists, skip")
        return

    start_time = time.time()

    with open(agg_batch_json_file, "r") as json_file:
        data_sub = json.load(json_file)

    elapsed = time.time() - start_time
    print(f"opened {agg_batch_json_file} after {elapsed:.2f}s", flush=True)

    rows_idx = []
    cols_idx = []

    row_col_values_meth = []
    row_col_values_notmeth = []

    rows_idx_score = []
    cols_idx_score = []
    row_col_values_score = []

    batch_bcs = list(data_sub.keys())

    for idx, bc in enumerate(batch_bcs):
        total_bases = len(data_sub[bc])

        elapsed = time.time() - start_time
        print(f"total_bases for {bc} = {total_bases} in {elapsed:.2f}s", flush=True)

        if total_bases < 1000:
            continue

        # convert ['chr_100', 'Z'] to ['chr_100_Z'] and count dedup
        # maybe save this stats later
        # one call also do this in aggregate_quad_parts

        dedup = {}

        for b in data_sub[bc]:
            seq_counter(dedup, b)

        # convert back to ['chr', '100', 'Z'] format

        triplets = [d.split("_") for d in dedup.keys()]
        triplets = pd.DataFrame(triplets)

        elapsed = time.time() - start_time
        print(
            f"deduped bases for {bc} = {triplets.shape[0]} in {elapsed:.2f}s",
            flush=True,
        )

        triplets[3] = triplets[1]
        triplets = triplets[[0, 1, 3, 2]].copy()
        triplets = triplets[triplets[0].str.contains("chr")].copy()
        bed_bases = pybedtools.BedTool.from_dataframe(triplets)

        overlaps_bed = bed_bases.intersect(bed_genes, wb=True)
        df = overlaps_bed.to_dataframe()
        df = df[["thickEnd", "name"]].copy()
        df.columns = ["gene_id", "meth"]
        Z_cnt = df[df.meth == "Z"].groupby("gene_id").size()
        z_cnt = df[df.meth == "z"].groupby("gene_id").size()

        Z_cnt.name = "Z_cnt"
        z_cnt.name = "z_cnt"

        mrg = pd.merge(Z_cnt, z_cnt, how="outer", left_index=True, right_index=True)
        mrg = mrg.fillna(0)
        cell_mC = mrg.sum()

        mrg["cov"] = mrg.sum(axis=1)

        cell_mC_ratio = cell_mC.loc["Z_cnt"] / (
            cell_mC.loc["Z_cnt"] + cell_mC.loc["z_cnt"]
        )
        # mrg['meth'] = mrg.Z_cnt / mrg['cov']
        # mrg['ratio'] = mrg['meth'] / cell_mC_ratio
        mrg["diff"] = mrg.Z_cnt / mrg["cov"] - cell_mC_ratio

        diff_pos = mrg[(mrg["diff"] > 0) & (mrg["cov"] > 1)]["diff"] / (
            1 - cell_mC_ratio
        )
        diff_neg = mrg[(mrg["diff"] < 0) & (mrg["cov"] > 1)]["diff"] / cell_mC_ratio

        elapsed = time.time() - start_time
        print(f"dataframes built in {elapsed:.2f}s", flush=True)

        row_vals = (np.ones(len(mrg), dtype=int) * idx).tolist()

        col_vals = [gene_idx_map[key] for key in mrg.index]
        # col_vals = [chr_idx_dict[key] for key in mrg.index]

        rows_idx.extend(row_vals)  # cells
        cols_idx.extend(col_vals)  # genes

        row_col_values_meth.extend(mrg["Z_cnt"].tolist())
        row_col_values_notmeth.extend(mrg["z_cnt"].tolist())

        for diff in [diff_pos, diff_neg]:
            row_vals = (np.ones(len(diff), dtype=int) * idx).tolist()
            col_vals = [gene_idx_map[key] for key in diff.index]
            rows_idx_score.extend(row_vals)  # cells
            cols_idx_score.extend(col_vals)  # genes
            row_col_values_score.extend(diff.tolist())

        elapsed = time.time() - start_time
        print(f"all matrices built in {elapsed:.2f}s", flush=True)

    shape = (len(batch_bcs), len(gene_idx_map))

    csr = csr_matrix(
        (row_col_values_meth, (rows_idx, cols_idx)), shape=shape, dtype="float32"
    )
    csr.eliminate_zeros()
    write_mtx(gene_mtx_dir, batch, "meth", csr)

    csr = csr_matrix(
        (row_col_values_notmeth, (rows_idx, cols_idx)), shape=shape, dtype="float32"
    )
    csr.eliminate_zeros()
    write_mtx(gene_mtx_dir, batch, "notmeth", csr)

    csr = csr_matrix(
        (row_col_values_score, (rows_idx_score, cols_idx_score)),
        shape=shape,
        dtype="float32",
    )
    write_mtx(gene_mtx_dir, batch, "score", csr)

    elapsed = time.time() - start_time
    print(f"all matrices saved in {elapsed:.2f}s", flush=True)

def make_allc_tsv(indir, sample, batch, context):
    
    agg_batch_json_file = f"{indir}/{sample}/split/quad_agg_{batch}_{context}.json"
    
    print(agg_batch_json_file)
    with open(agg_batch_json_file, "r") as json_file:
        data_sub = json.load(json_file)
        
    allcools_dir = f"{indir}/{sample}/allcools"
    if not os.path.exists(allcools_dir):
        os.makedirs(allcools_dir)
        print(f"{allcools_dir} created")
    else:
        print(f"{allcools_dir} already exists")
    
    if context.startswith('CpG'):
        bed_context_col = 'CG'
    else:
        bed_context_col = 'CA'
    
    batch_bcs = list(data_sub.keys())

    for idx, bc in enumerate(batch_bcs):
        total_bases = len(data_sub[bc])

        print(bc, total_bases)
        
        allc_tsv_file = f'{allcools_dir}/{bc}_{context}.tsv'
        
        if os.path.isfile(allc_tsv_file):
            print(allc_tsv_file, " exists, skip", flush=True)
            continue
        
        dedup = set(data_sub[bc])
        triplets = [d.split("_") for d in dedup]
        triplets = pd.DataFrame(triplets)

        Z_cnt = triplets[triplets[2] == "Z"].groupby([0, 1]).size()
        z_cnt = triplets[triplets[2] == "z"].groupby([0, 1]).size()
        
        Z_cnt.name = "Z_cnt"
        z_cnt.name = "z_cnt"
        
        mrg = pd.merge(Z_cnt, z_cnt, how="outer", left_index=True, right_index=True)
        mrg = mrg.fillna(0)
        cell_mC = mrg.sum()
        mrg = mrg.reset_index()
        mrg['strand'] = '+'
        mrg['cov'] = mrg['z_cnt'] + mrg['Z_cnt']
        mrg['sig_meth'] = 1
        mrg['context'] = bed_context_col
        mrg = mrg[[0, 1, 'strand', 'context', 'Z_cnt', 'cov', 'sig_meth']]
        mrg['Z_cnt'] = mrg['Z_cnt'].astype(int)
        mrg['cov'] = mrg['cov'].astype(int)
        mrg.to_csv(allc_tsv_file,index=None,header=None,sep='\t')
        

def load_mtx(file):
    
    return scipy.io.mmread(file)


def stack_mtx_windows(indir, sample, window, chr_idx_dict, context, cores):
    with open(f"{indir}/{sample}/{sample}_whitelist_batches.json", "r") as file:
        bc_splits = json.load(file)

    all_bcs_list = []
    for b in bc_splits:
        all_bcs_list.extend(b)

    chr_idx_string = [idx[0] + "_" + str(idx[1]) for idx in chr_idx_dict.keys()]
    print("chr_idx_dict lenght = ", len(chr_idx_dict), window, context, flush=True)

    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        mtx_folder = f"{indir}/{sample}/counts_w_{window}_m{context}"

        for mtx_type in ["notmeth", "meth", "score"]:
            adata_file = f"{mtx_folder}/adata_{mtx_type}.h5ad"

            if os.path.isfile(adata_file):
                print(adata_file, " exists, skip", flush=True)
                continue

            files = os.listdir(mtx_folder)

            file_pattern = f"_{mtx_type}.mtx.gz"
            files = sorted([f"{mtx_folder}/{f}" for f in files if file_pattern in f])
            print(files[0])
            print("load all mtx", flush=True)

            all_mtx = list(executor.map(load_mtx, files))

            print("stacking", flush=True)
            merged_matrix = vstack(all_mtx)
            merged_matrix = csr_matrix(merged_matrix, dtype="float32")
            print("AnnData making", flush=True)
            adata = AnnData(merged_matrix)
            adata.obs.index = all_bcs_list
            adata.var.index = chr_idx_string
            print("AnnData writing", flush=True)
            adata.write_h5ad(adata_file, compression="gzip")

        coverage_adata_file = f"{mtx_folder}/adata_coverage.h5ad"

        if os.path.isfile(coverage_adata_file):
            print(coverage_adata_file, " exists, skip", flush=True)
            pass
        else:
            print(coverage_adata_file, " will be made", flush=True)

            adatas = []
            for mtx_type in ["notmeth", "meth"]:
                adatas.append(sc.read_h5ad(f"{mtx_folder}/adata_{mtx_type}.h5ad"))

            adatas[0].X = adatas[0].X + adatas[1].X
            adatas[0].write_h5ad(coverage_adata_file, compression="gzip")


def stack_mtx_genes(indir, sample, gene_idx_list, context, cores):
    with open(f"{indir}/{sample}/{sample}_whitelist_batches.json", "r") as file:
        bc_splits = json.load(file)

    all_bcs_list = []
    for b in bc_splits:
        all_bcs_list.extend(b)

    print("gene_idx_string length = ", len(gene_idx_list), flush=True)

    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        mtx_folder = f"{indir}/{sample}/counts_gene_m{context}"

        for mtx_type in ["notmeth", "meth", "score"]:
            adata_file = f"{mtx_folder}/adata_{mtx_type}.h5ad"

            if os.path.isfile(adata_file):
                print(adata_file, " exists, skip", flush=True)
                continue

            files = os.listdir(mtx_folder)

            file_pattern = f"_{mtx_type}.mtx.gz"
            files = sorted([f"{mtx_folder}/{f}" for f in files if file_pattern in f])
            print(files[0])
            print("load all mtx", flush=True)

            all_mtx = list(executor.map(load_mtx, files))

            print("stacking", flush=True)
            merged_matrix = vstack(all_mtx)
            merged_matrix = csr_matrix(merged_matrix, dtype="float32")
            print("AnnData making", flush=True)
            adata = AnnData(merged_matrix)
            adata.obs.index = all_bcs_list
            adata.var.index = gene_idx_list
            print("AnnData writing", flush=True)
            adata.write_h5ad(adata_file, compression="gzip")

        coverage_adata_file = f"{mtx_folder}/adata_coverage.h5ad"

        if os.path.isfile(coverage_adata_file):
            print(coverage_adata_file, " exists, skip", flush=True)
            pass
        else:
            print(coverage_adata_file, " will be made", flush=True)

            adatas = []
            for mtx_type in ["notmeth", "meth"]:
                adatas.append(sc.read_h5ad(f"{mtx_folder}/adata_{mtx_type}.h5ad"))

            adatas[0].X = adatas[0].X + adatas[1].X
            adatas[0].write_h5ad(coverage_adata_file, compression="gzip")


def Mbias_parser(Mbias_filename):
    context_dict = {}
    current_key = None
    all_lines = []
    with open(Mbias_filename, "r") as file:
        for line in file:
            line = line.strip()
            all_lines.append(line)
            if "context" in line:
                current_key = line
                context_dict[current_key] = []
            elif current_key is not None and line != "":
                if "====" not in line and "position" not in line:
                    pieces = line.split("\t")
                    if "" not in pieces:
                        context_dict[current_key].append(pieces)

        all_dfs = [
            pd.DataFrame(context_dict[key]).astype(float) for key in context_dict.keys()
        ]
        for i, df in enumerate(all_dfs):
            df["context"] = list(context_dict.keys())[i]
        result_df = pd.concat(all_dfs, ignore_index=True)
        result_df.columns = [
            "pos",
            "meth_cnt",
            "unmeth_cnt",
            "prct_meth",
            "cov",
            "context",
        ]

        result_df.to_csv(Mbias_filename.replace(".txt", ".csv"))
