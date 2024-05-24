import pysam
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import edlib
import json
import subprocess
import gzip
import scipy.io
from scipy.sparse import csr_matrix
import re
import mappy
import matplotlib.pyplot as plt
import seaborn as sns

N_read_extract = 100000

print(N_read_extract)


def seq_counter(seq_dict, seq_instance):
    if seq_dict.get(seq_instance) is None:
        seq_dict[seq_instance] = 1
    else:
        seq_dict[seq_instance] += 1


def quad_dict_store(quad_dict, quad_key, quad_items):
    if quad_dict.get(quad_key) is None:
        quad_dict[quad_key] = [quad_items]
    else:
        quad_dict[quad_key].extend([quad_items])


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

        command_R1 = f"zcat {R1} | split -a 3 -l {int(lines)} -d --additional-suffix=.fastq - {split_R1_name}"
        command_R2 = command_R1.replace("_R1", "_R2")
        command_R3 = command_R1.replace("_R1", "_R3")

        subprocess.call(f"{command_R1} & {command_R2} & {command_R3}", shell=True)


def quality_calc(seq, quals, bases_dict, quals_dict):
    for i in range(len(seq)):
        if bases_dict.get(str(i)) is None:
            bases_dict[str(i)] = {}
            seq_counter(bases_dict[str(i)], seq[i])
        else:
            seq_counter(bases_dict[str(i)], seq[i])

        if quals_dict.get(str(i)) is None:
            quals_dict[str(i)] = {}
            seq_counter(quals_dict[str(i)], quals[i])
        else:
            seq_counter(quals_dict[str(i)], quals[i])


def quality_df(quals_dict):
    quals_df = pd.DataFrame(quals_dict)
    quals_df = quals_df.T
    quals_df = quals_df.fillna(0)
    quals_df = quals_df.stack()
    quals_df = quals_df.reset_index()
    # quals_df.columns = ['base', 'quality', 'tot_count']
    # quals_df['mult'] = quals_df.quality * quals_df.tot_count
    # quals_df_grouped = quals_df.groupby('base').sum()
    quals_df.columns = ["position", "base_qual", "tot_count"]
    quals_df.position = quals_df.position.astype("int")
    # quals_df[quals_df.position.isin(np.arange(10))]
    counts_df = quals_df.groupby("position").sum()
    quals_df["position_cnt"] = quals_df.position.apply(
        lambda x: counts_df.loc[x].tot_count
    )
    quals_df["freq"] = quals_df.tot_count / quals_df.position_cnt * 100
    return quals_df


def extract_clean_fastq(indir, sample, part, limit):
    i = 0

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
        for r1, r2, r3 in tqdm(zip(R1, R2, R3)):
            i += 1

            seq1 = r1.sequence
            seq2 = r2.sequence
            seq3 = r3.sequence

            len1 = len(seq1)
            len3 = len(seq3)

            # bc = seq2 for sciMET data, matched already
            # seq_counter(bcs_dict,bc)

            if do_qc and i % 500 == 0:
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

                    R1_clean.write(f"@{r1.name}_1_{bc}\n")
                    R1_clean.write(f"{r1.sequence[11:-2]}\n")
                    R1_clean.write("+\n")
                    R1_clean.write(f"{r1.quality[11:-2]}\n")

                    R3_clean.write(f"@{r3.name}_2_{bc}\n")
                    R3_clean.write(f"{r3.sequence[2:-2]}\n")
                    R3_clean.write("+\n")
                    R3_clean.write(f"{r3.quality[2:-2]}\n")

            if i > N_read_extract and limit:
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
        return

    bcs_dict = {}

    samfile = pysam.AlignmentFile(bam, "r")

    for read in tqdm(samfile.fetch(until_eof=True)):
        i += 1

        bc = read.qname.split("_")[-1]

        seq_counter(bcs_dict, bc)

        if i > N_read_extract and limit:
            break

    with open(bcs_json, "w") as json_file:
        json.dump(bcs_dict, json_file)


def tag_bam_with_barcodes(indir, sample, part, limit):
    matching_csv = pd.read_csv(f"{indir}/{sample}/{sample}_raw_to_tenx_whitelist.csv")
    raw_to_tenx = dict(zip(matching_csv.bc, matching_csv.tenx_whitelist))

    sam_tag = f"{indir}/{sample}/split/{sample}.part_{part}_tagged.bam"
    sam = f"{indir}/{sample}/split/output_{part}/{sample}_R1.part_{part}_clean_bismark_bt2_pe.bam"

    if os.path.isfile(sam_tag):
        print(sam_tag, " exists, skip")
        return

    samfile = pysam.AlignmentFile(sam, "rb", threads=8)
    tagged_bam = pysam.AlignmentFile(sam_tag, "wb", template=samfile, threads=8)

    i = 0
    for read in tqdm(samfile.fetch(until_eof=True)):
        i += 1
        raw_barcode = read.qname.split("_")[-1]
        if raw_barcode in raw_to_tenx:
            matched_barcode = raw_to_tenx[raw_barcode]

            read.set_tag("CB", matched_barcode)
            tagged_bam.write(read)
        if i > N_read_extract and limit:
            break

    tagged_bam.close()
    samfile.close()


def aggregate_bc_dicts(indir, sample):
    dir_split = f"{indir}/{sample}/split/"
    files = os.listdir(dir_split)
    jsons = sorted([f for f in files if "_bcs.json" in f])

    agg_read_csv = f"{indir}/{sample}/{sample}_agg_cnt_raw_bcs.csv"

    if os.path.isfile(agg_read_csv):
        print(agg_read_csv, " exists, skip")
        return

    data_agg = {}

    for i in tqdm(range(len(jsons))):
        with open(f"{dir_split}{jsons[i]}", "r") as json_file:
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


def processing_matching(indir, sample, AS_min=10):
    all_AS = []
    all_pairs = []

    samfile = pysam.AlignmentFile(f"{indir}/{sample}/{sample}_matching.sam", "rb")

    for read in tqdm(samfile.fetch()):
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
    matched.to_csv(
        f"{indir}/{sample}/{sample}_matching_raw_to_tenx_passing.csv", index=None
    )


def filered_barcodes(indir, sample, read_cnt_min=10000):
    matched = pd.read_csv(f"{indir}/{sample}/{sample}_matching_raw_to_tenx_passing.csv")
    matched_sub = matched[["read_cnt", "tenx_whitelist"]]
    matched_sub_grouped = matched_sub.groupby("tenx_whitelist").sum()

    agg_bcs = matched_sub_grouped[matched_sub_grouped.read_cnt > 0]
    agg_bcs = agg_bcs.sort_values(by="read_cnt", ascending=False)
    agg_bcs["log10_read_cnt"] = np.log10(agg_bcs.read_cnt)

    plt.figure(figsize=(4, 3))
    log10_ranks = np.log10(np.arange(1, len(agg_bcs) + 1))
    log10_cnts = agg_bcs.log10_read_cnt
    plt.plot(log10_ranks, log10_cnts)  # ,label='Rank Plot of Reads')
    plt.xlabel("Log10 Ranks")
    plt.ylabel("Log10 Read Counts")
    plt.savefig(f"{indir}/{sample}/{sample}_rankplot.pdf", bbox_inches="tight")

    plt.figure(figsize=(4, 3))
    sns.histplot(log10_cnts[log10_cnts > 2], bins=100)
    plt.xlabel("Log10 Read Counts")
    plt.savefig(f"{indir}/{sample}/{sample}_histogram.pdf", bbox_inches="tight")

    agg_bcs_filtered = agg_bcs[agg_bcs.read_cnt > read_cnt_min]
    matched_filtered = matched[
        matched.tenx_whitelist.isin(agg_bcs_filtered.index)
    ].copy()
    print("saving raw_to_tenx_whitelist")
    matched_filtered.to_csv(
        f"{indir}/{sample}/{sample}_raw_to_tenx_whitelist.csv", index=None
    )

    bcs_random = agg_bcs_filtered.read_cnt.sample(frac=1)
    bcs = agg_bcs_filtered.loc[bcs_random.index]
    bcs = bcs[["read_cnt"]].copy()
    bcs = bcs.reset_index()
    bcs.columns = ["bc", "count"]
    bcs = bcs.set_index("bc")
    print("saving whitelist")
    bcs.to_csv(f"{indir}/{sample}/{sample}_whitelist.csv")


def save_quad_batch_json(indir, sample, part, context, limit):
    batch = str(1).zfill(3)
    batch_json = (
        f"{indir}/{sample}/split/quads_part_{part}_batch_{batch}_{context}.json"
    )

    if os.path.isfile(batch_json):
        print(batch_json, " exists, skip")
        return

    bcs = pd.read_csv(f"{indir}/{sample}/{sample}_whitelist.csv")
    matching_csv = pd.read_csv(f"{indir}/{sample}/{sample}_raw_to_tenx_whitelist.csv")
    raw_to_tenx = dict(zip(matching_csv.bc, matching_csv.tenx_whitelist))

    dir_split = f"{indir}/{sample}/split"

    meth_file = f"{dir_split}/output_{part}/{context}_{sample}_R1.part_{part}_clean_bismark_bt2_pe.txt.gz"

    # sub_batch_N = int(np.sqrt(len(bcs)))+1 # better load balancing for very large datasets
    sub_batch_N = 20

    bc_splits = np.array_split(bcs.bc, sub_batch_N)

    quad_dict = {}

    # for bc in bcs.bc.to_list():
    #    quad_dict[bc] = []

    i = 0

    conversion = False
    if context == "Non_CpG_context":
        conversion = True

    Non_CpG_to_CpG_dict = {"x": "z", "h": "z", "X": "Z", "H": "Z"}

    all_failed_bc = []

    with gzip.open(meth_file, "rt") as f:
        for line in tqdm(f):
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


def write_mtx(indir, sample, batch, window, context, state, csr):
    window_mtx_dir = f"{indir}/{sample}/counts_w_{window}_m{context}"

    if not os.path.exists(window_mtx_dir):
        try:
            os.makedirs(window_mtx_dir)
            print(f"{window_mtx_dir} created")
        except:
            print(f"{window_mtx_dir} already created")
    else:
        print(f"{window_mtx_dir} already exists")

    csr_file = f"{window_mtx_dir}/b_{batch}_{state}.mtx.gz"

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

    with open(agg_batch_json_file, "r") as json_file:
        data_sub = json.load(json_file)

    rows_idx = []
    cols_idx = []

    row_col_values_meth = []
    row_col_values_notmeth = []

    rows_idx_score = []
    cols_idx_score = []
    row_col_values_score = []

    batch_bcs = list(data_sub.keys())

    for idx, bc in enumerate(tqdm(batch_bcs)):
        if len(data_sub[bc]) == 0:
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

    shape = (len(batch_bcs), len(chr_idx_dict))

    csr = csr_matrix(
        (row_col_values_meth, (rows_idx, cols_idx)), shape=shape, dtype="float32"
    )
    csr.eliminate_zeros()
    write_mtx(indir, sample, batch, window, context, "meth", csr)

    csr = csr_matrix(
        (row_col_values_notmeth, (rows_idx, cols_idx)), shape=shape, dtype="float32"
    )
    csr.eliminate_zeros()
    write_mtx(indir, sample, batch, window, context, "notmeth", csr)

    csr = csr_matrix(
        (row_col_values_score, (rows_idx_score, cols_idx_score)),
        shape=shape,
        dtype="float32",
    )
    write_mtx(indir, sample, batch, window, context, "score", csr)
