"""Output table that tracks TCONS, XLOC, GENE_ID_EXTERNAL, GENE_ID_INTERNAL.

Here, 'external' means the gene id of the species used to get the RNA-seq Tx info.
Likewise, 'internal' means the gene id from the species that you are actually interested in.
"""

import pandas as pd


# Settings

# input
cuffcmp_tracking = snakemake.input.cuffcmp_tracking
ortholog_table = snakemake.input.ortholog_table


#output
ids_no_diff_expr = snakemake.output.ids_no_diff_expr





def load_ortho_table(path):
    """Load and return a dataframe representing only the needed columns from the ortholog table."""
    df = pd.read_csv(path, sep=',')

    df = df.iloc[:,0:2]
    df.columns = ["gene_id_external","gene_id_internal"]

    return df

def load_tracking_table(path):
    """Load and return a dataframe representing only the needed columns from the cuffcompare tracking table."""
    df = pd.read_csv(path, sep='\t', names=["tcons_id","xloc_id","gene_id_tx_id","class_code","info"])

    df['gene_id_external'] = df.gene_id_tx_id.apply(lambda i: i.split('|')[0] if i != '-' else i)
    columns = ["tcons_id","xloc_id","gene_id_external"]

    return df[columns].copy()

def combine_tracking_and_orthologs(tracking, orthologs):
    """Merge and return the combined table."""
    df = pd.merge(left=tracking, right=orthologs,
                  how="outer",
                  on="gene_id_external", left_on=None, right_on=None,
                  left_index=False, right_index=False,
                  sort=False, suffixes=('_x', '_y'), copy=True, indicator=False).fillna('-')
    return df



tracking_df = load_tracking_table(path=cuffcmp_tracking)
ortho_df = load_ortho_table(path=ortholog_table)

combined_df = combine_tracking_and_orthologs(tracking=tracking_df, orthologs=ortho_df)


combined_df.to_csv(path_or_buf=ids_no_diff_expr,
                   sep=',',
                   header=True, index=False,)
