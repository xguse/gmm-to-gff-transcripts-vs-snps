#!/usr/bin/env python
"""Functions used to compare results from original SNPs and those from SNPs in 'strong linkage' with the original SNPs."""

# Imports
import pandas as pd

from munch import Munch


# Metadata
__author__ = "Gus Dunn"
__email__ = "w.gus.dunn@gmail.com"


def filter_by_distance(df, d=5000):
    """Return only rows that are within ``d`` absolute distance of target SNP."""
    return df.query(""" abs(distance) <= {d} """.format(d=d))

def set_comparison_table(a, b, a_name=None, b_name=None):
    """Return `Munch` of actual setwise comparison and the resulting count table."""
    a = set(a)
    b = set(b)

    data = [pd.Series((sorted(list(a))),name='{a}'.format(a=a_name)),
            pd.Series((sorted(list(b))),name='{b}'.format(b=b_name)),
            pd.Series((sorted(list(a & b))),name='{a} AND {b}'.format(a=a_name,b=b_name)),
            pd.Series((sorted(list(a - b))),name='{a} NOT {b}'.format(a=a_name,b=b_name)),
            pd.Series((sorted(list(b - a))),name='{b} NOT {a}'.format(a=a_name,b=b_name)),
            pd.Series((sorted(list(a | b))),name='{a} OR {b}'.format(a=a_name,b=b_name)),
            pd.Series((sorted(list(a ^ b))),name='one OR other NOT both'),]

    df = pd.DataFrame(data).T
    counts = pd.DataFrame(df.count(),columns=['count'])

    return Munch({'genes':df, 'counts':counts})


# Loading
def load_snps_near_de(path):
    """Load file with appropriate function."""
    return pd.read_csv(path)

def load_nearest_k(path):
    """Load file with appropriate function."""
    return pd.read_excel(path)

# Set up
def split_official_novel(df):
    """Return ``Munch`` of ``set`` objs representing the 'official_annotations' or 'novel_mapped_tx' for ``df``."""
    split = Munch()
    try:
        split.official_genes = set(df.query(""" feature_set_name == 'official_annotations' """).proximal_id)
        split.novel_tx = set(df.query(""" feature_set_name == 'novel_mapped_tx' """).proximal_id)
    except AttributeError:
        split.official_genes = set(df.query(""" feature_set_name == 'official_annotations' """).name)
        split.novel_tx = set(df.query(""" feature_set_name == 'novel_mapped_tx' """).name)

    return split

# # ?
# def get_coords(gname, win, source='snps_near_de'):
#     """Return a coords object suitable for use in `get_genome_image`."""
#     valid = ['snps_near_de', 'nearest_k']
#     if source not in valid:
#         raise ValueError('Valid values are {values}'.format(values=vaild))
#
#     coords = Munch()
#
#     if source == 'snps_near_de':
#         name_col = 'proximal_id'
#     elif source == 'nearest_k':
#         name_col = 'name'
#
#
#     chr, loc = data.query(""" {name_col} == '{gname}' """.format(name_col=name_col, gname=gname))[['SNP_chrom','SNP_end']].astype(tuple)
#     coords.chr =
