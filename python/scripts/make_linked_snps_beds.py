"""Convert vcftools outputs of linked SNPs to BED files."""

import pandas as pd

snp_file = snakemake.input.snp_file
snp_bed = snakemake.output.snp_bed


snps = pd.read_csv(snp_file, sep='\t')


bed_table = snps[["CHR2","POS2"]].copy()
bed_table = bed_table.rename(columns={"POS2":"END", "CHR2":"CHROM"}).copy()
bed_table['START'] = bed_table.END - 1
bed_table = bed_table[["CHROM","START","END"]].sort_values(by=["CHROM","START"]).copy()

bed_table.to_csv(snp_bed, header=False, index=False, sep='\t')
