"""Describe here what this rule accomplishes."""

import pandas as pd
import numpy as np


# Settings
edger_results_labels = snakemake.params.edger_results_labels
cufflinks_results_labels = snakemake.params.cufflinks_results_labels

# input
edger_results = snakemake.input.edger_results
cufflinks_results = snakemake.input.cufflinks_results
ids_no_diff_expr = snakemake.input.ids_no_diff_expr

#output
ids_with_diff_expr = snakemake.output.ids_with_diff_expr



def load_and_filter_diff_expr_data(path,ids,comparison,program,fdr_thresh):
    """Return new dataframe that has standardized and filtered the DE input tables.

    `path` (str):
        location of input file

    `ids` (dataframe):
        with the following columns
           - tcons_id
           - xloc_id
           - gene_id_external
           - gene_id_internal

    `comparison` (str):
        describe the RNA-seq analysis run ('midgut', 'salivary gland', etc)

    `program` (str):
        one of ['edger', 'cufflinks']

    `fdr_thresh` (float):
        defining multiple testing significance threshold above which DE tests should NOT be reported
    """
    column_conversions = {'edger': {'Gene_Name': 'gene_id_external',
                                    'Gene_ID': 'xloc_id',
                                    'logFC': 'lg2_fc',
                                    'PValue': 'p',
                                    'FDR': 'fdr'},
                          'cufflinks': {'gene': 'gene_id_external',
                                        'gene_id': 'xloc_id',
                                        'log2.fold_change.': 'lg2_fc',
                                        'p_value': 'p',
                                        'q_value': 'fdr'},
                         }

    keep_columns = ["de_id", "xloc_id", "tcons_id","gene_id_external","gene_id_internal","lg2_fc","p","fdr","comparison","program"]

    de_id_program_map = {'edger': 'EDGR',
                         'cufflinks': 'CUFF',
                        }

    # Load
    df = pd.read_csv(path, sep='\t')

    # Convert Columns
    df = df.rename(columns=column_conversions[program])

    # Make missing fdr values NaN
    df['fdr'] = df.fdr.apply(lambda i: np.nan if i == '-' else i)

    # Filter for fdr
    df = df.query(""" fdr <= 0.05 """).copy()

    # Add Columns
    df['program'] = program
    df['comparison'] = comparison
    df['de_id'] = generate_de_ids(df=df,
                                  de_type=de_id_program_map[program],
                                  type_mod='|{comparison}'.format(comparison=comparison),
                                  nlen=7)

    # Join external and internal IDS
    df = pd.merge(left=df, right=ids_no_diff_expr,
                  how='left',
                  on=None, left_on=None, right_on=None,
                  left_index=False, right_index=False,
                  sort=False, suffixes=('_x', '_y'), copy=True, indicator=False).fillna('-')

    # Retain only needed columns
    df = df[keep_columns]

    # Return dataframe
    return df.copy()


def generate_de_ids(df,de_type,type_mod='',nlen=7):
    """Generate unique tracking IDs for each statistical test of diff expr."""
    template = '{de}{mod}_{{0:0{nlen}d}}'.format(de=de_type, mod=type_mod, nlen=nlen)

    return [template.format(n) for n in range(1,len(df)+1)]


ids_no_diff_expr = pd.read_csv(ids_no_diff_expr)

table_list = []

# Load EDGER DE reults
for name, path in zip(edger_results_labels, edger_results):
    df = load_and_filter_diff_expr_data(path=path,
                                        ids=ids_no_diff_expr,
                                        comparison=name,
                                        program='edger',
                                        fdr_thresh=0.05)
    table_list.append(df)

# Load CUFFLINKS DE reults
for name, path in zip(cufflinks_results_labels, cufflinks_results):
    df = load_and_filter_diff_expr_data(path=path,
                                        ids=ids_no_diff_expr,
                                        comparison=name,
                                        program='cufflinks',
                                        fdr_thresh=0.05)
    table_list.append(df)

# Concat all result files into single dataframe
combined = pd.concat(objs=table_list, axis=0)

# Write out the resulting dataframe
combined.to_csv(path_or_buf=ids_with_diff_expr,
                sep=',',
                header=True, index=False,)
