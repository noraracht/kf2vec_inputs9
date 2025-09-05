# Usage: python closeast2taxid.py /Users/nora/Documents/ml_metagenomics/cami_long_reads/tax_assignment /Users/nora/Documents/ml_metagenomics/65K_TOL/assembly_summary_refseq-after20190108-sampled_500_lineage.txt


import argparse
import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
import fnmatch
from collections import Counter
from scipy.stats import poisson



# Parse inputs
parser = argparse.ArgumentParser()
parser.add_argument('input_dir')
parser.add_argument('name_to_id')
parser.add_argument('length_file')
args = parser.parse_args()

input_dir  = args.input_dir
name_to_id = args.name_to_id
length_file = args.length_file

# Mapping files
df_taxid = pd.read_csv(name_to_id, index_col=None, header=0, sep='\t')
df_length = pd.read_csv(length_file, index_col=None, header=None, sep='\t')
df_length[1] = df_length[1].astype(int)

my_map_dict_lens = pd.Series(df_length[1].values, index=df_length[0]).to_dict()
my_map_dict_sp = pd.Series(df_taxid['species'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_gen = pd.Series(df_taxid['genus'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_fam = pd.Series(df_taxid['family'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_ordr = pd.Series(df_taxid['order'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_cls = pd.Series(df_taxid['class'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_phyl = pd.Series(df_taxid['phylum'].values, index=df_taxid["tree_label"]).to_dict()
my_map_dict_spkindm = pd.Series(df_taxid['superkingdom'].values, index=df_taxid["tree_label"]).to_dict()

# Output files
out_species_file = os.path.join(input_dir, "summary_subset_v10_taxid.species")
out_log_file = os.path.join(input_dir, "summary_subset_v10.log")

# Remove old files if exist
for f in [out_species_file, out_log_file]:
    if os.path.exists(f):
        os.remove(f)

# Initialize final result containers
df_species = pd.DataFrame(columns=["SEQUENCEID", "TAXID"])
df_log = pd.DataFrame(columns=["SEQUENCEID", "NumColumns", "Rank", "AssignedTax"])

# Find all subtree files and sort numerically
subtree_files = [f for f in os.listdir(input_dir) if re.match(r"apples_input_di_mtrx_subtree_\d+\.csv", f)]
subtree_files_sorted = sorted(subtree_files, key=lambda x: int(re.search(r"subtree_(\d+)\.csv", x).group(1)))

for dist_mtrx in subtree_files_sorted:
    df = pd.read_csv(os.path.join(input_dir, dist_mtrx), index_col=0, header=0, sep='\t')

    for query_label in df.index:
        out = df.loc[query_label].sort_values(ascending=True)
        columns_sorted = out.index

        column_subset = []
        maxdelt = -1
        p = 0
        i = 0
        first = -1
        selected = 5000

        # Safe lookup for length
        clen = my_map_dict_lens.get(query_label, 1)

        for col in columns_sorted:
            a = float(out[col])
            if first == -1:
                first = a
            i += 1
            if (maxdelt != -1 and (a-p)/maxdelt>3) or (maxdelt==-1 and first!=a and a/p>1.05) or 1-poisson.cdf(a*clen/100, first*clen/100)<0.1:
                if selected == 5000:
                    selected = i
                break
            if first != a and a-p>maxdelt:
                maxdelt = a-p
            p = a

        column_subset = columns_sorted[:i-1]

        # Counters with safe get
        column_subset_sp = Counter([my_map_dict_sp.get(n1, -1) for n1 in column_subset if my_map_dict_sp.get(n1, -1) != -1])
        column_subset_gen = Counter([my_map_dict_gen.get(n1, -1) for n1 in column_subset if my_map_dict_gen.get(n1, -1) != -1])
        column_subset_fam = Counter([my_map_dict_fam.get(n1, -1) for n1 in column_subset if my_map_dict_fam.get(n1, -1) != -1])
        column_subset_ordr = Counter([my_map_dict_ordr.get(n1, -1) for n1 in column_subset if my_map_dict_ordr.get(n1, -1) != -1])
        column_subset_cls = Counter([my_map_dict_cls.get(n1, -1) for n1 in column_subset if my_map_dict_cls.get(n1, -1) != -1])
        column_subset_phyl = Counter([my_map_dict_phyl.get(n1, -1) for n1 in column_subset if my_map_dict_phyl.get(n1, -1) != -1])
        column_subset_spkindm = Counter([my_map_dict_spkindm.get(n1, -1) for n1 in column_subset if my_map_dict_spkindm.get(n1, -1) != -1])

        # Majority vote
        maj_thresh = 0.75
        assigned = False
        for rank, counter in [('sp', column_subset_sp), ('gen', column_subset_gen),
                              ('fam', column_subset_fam), ('ordr', column_subset_ordr),
                              ('cls', column_subset_cls), ('phyl', column_subset_phyl),
                              ('spkindm', column_subset_spkindm)]:
            if len(counter) > 0 and counter.most_common(1)[0][1]/sum(counter.values()) > maj_thresh:
                s = counter.most_common(1)[0][0]
                df_species = pd.concat([df_species, pd.DataFrame([[query_label, s]], columns=df_species.columns)], ignore_index=True)
                df_log = pd.concat([df_log, pd.DataFrame([[query_label, i-1, rank, s]], columns=df_log.columns)], ignore_index=True)
                assigned = True
                break

        if not assigned:
            df_species = pd.concat([df_species, pd.DataFrame([[query_label, 1]], columns=df_species.columns)], ignore_index=True)
            df_log = pd.concat([df_log, pd.DataFrame([[query_label, i-1, "root", 1]], columns=df_log.columns)], ignore_index=True)

# Write final output once
with open(out_species_file, "w") as fo_sp:
    fo_sp.write("@Version:0.9.0\n")
    fo_sp.write("@SampleID:marmgCAMI2_short_read_pooled_gold_standard_assembly\n")
    fo_sp.write("@@SEQUENCEID\tTAXID\n")
    df_species.to_csv(fo_sp, sep='\t', index=False, header=False)

df_log.to_csv(out_log_file, sep='\t', index=False)
print("All subtrees processed, results written to single file.")

