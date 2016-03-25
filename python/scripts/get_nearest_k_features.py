"""For each SNP file, produce a bed representing the nearest k gene or mapped transcript features and its distance from the SNP."""

import pybedtools as pbt
import pandas as pd

k_number = snakemake.params.k_number

snp_beds = snakemake.input.snp_beds
gene_model_subtracted = snakemake.input.gene_model_subtracted
gene_models = snakemake.input.gene_models


nearest_features_beds = snakemake.output.nearest_features_beds
snps_in_features = snakemake.output.snps_in_features


headers = ["SNP_chrom",
           "SNP_start",
           "SNP_end",
           "feature_set_name",
           "chrom",
           "chromStart",
           "chromEnd",
           "name",
           "score",
           "strand",
           "thickStart",
           "thickEnd",
           "itemRgb",
           "blockCount",
           "blockSizes",
           "blockStarts",
           "distance"
          ]




for snp_bed, nearest_bed, feature_hit_file in zip(snp_beds, nearest_features_beds, snps_in_features):

    snp_bed = pbt.BedTool(snp_bed)
    gene_model_subtracted_bed = pbt.BedTool(gene_model_subtracted)
    gene_models_bed = pbt.BedTool(gene_models)

    k_nearest = snp_bed.closest([gene_model_subtracted_bed.fn,
                                 gene_models_bed.fn],
                                k=k_number,
                                names=['novel_mapped_tx', 'official_annotations'],
                                D='ref',    # Include SIGNED distances from SNP based on the ref genome
                                t='all',    # Return all members of a distance "tie"
                                mdb='each', # Return `k_number` of neighboors for each `names`
                                )

    k_nearest.saveas(nearest_bed)

    nearest_df = pd.read_csv(nearest_bed, sep="\t", names=headers)
    nearest_df

    in_features = nearest_df.query(""" abs(distance) <= 0 """)
    in_features.to_excel(feature_hit_file, index=False)
