"""Snakemake file."""

# See tutorial at: http://tiny.cc/snakemake_tutorial

import os

import yaml

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

from python.functions import *

ORIGINAL_CONFIG_AS_STRING = yaml.dump(config, default_flow_style=False)


#### COMMON RUN SPECIFICS ####

RUN_NAME = config["COMMON"]["RUN_NAME"]
OUT_DIR = "{base_dir}/{run_name}".format(base_dir=config["COMMON"]["OUT_DIR"], run_name=RUN_NAME)
OUT_DIR_MAIN = "{out_base}/main".format(out_base=OUT_DIR)
OUT_DIR_PREP = "{out_base}/prep".format(out_base=OUT_DIR)


# collect meta-target 'inputs'
PREP = []
MAIN = []

############ BEGIN PIPELINE RULES ############


#### SAVE_RUN_CONFIG ####
SAVE_RUN_CONFIG_OUT = OUT_DIR+"/{RUN_NAME}.yaml".format(RUN_NAME=RUN_NAME)

rule save_run_config:
    input:
    output:
        file=SAVE_RUN_CONFIG_OUT

    run:
        with open(output.file, 'w') as cnf_out:
            cnf_out.write(ORIGINAL_CONFIG_AS_STRING)

PREP.append(rules.save_run_config.output)
MAIN.append(rules.save_run_config.output)

# ------------------------- #
#### ANNOTATIONS_VIA_FASTA ####
ANNOTATIONS_VIA_FASTA = config["ANNOTATIONS_VIA_FASTA"]
ANNOTATIONS_VIA_FASTA_OUT = OUT_DIR_MAIN+"/annotations_via_fasta"
ANNOTATIONS_XLS = ANNOTATIONS_VIA_FASTA_OUT+"/annotations_via_fasta.xls"
TX_FASTA = ANNOTATIONS_VIA_FASTA["TX_FASTA"]
ORTHOLOG_TABLE = ANNOTATIONS_VIA_FASTA["ORTHOLOG_TABLE"]

# ---
rule annotations_via_fasta:
    input:
        tx_fasta=TX_FASTA,
        ortholog_table=ORTHOLOG_TABLE
    output:
        annotations_xls=ANNOTATIONS_XLS

    script:
        "python/scripts/annotations_via_fasta.py"



MAIN.append(rules.annotations_via_fasta.output)

# ------------------------- #
#### FILTER_PSL_TO_BED ####
FILTER_PSL_TO_BED = config["FILTER_PSL_TO_BED"]
PSL = FILTER_PSL_TO_BED["PSL"]

FILTER_PSL_TO_BED_OUT = OUT_DIR_MAIN+"/filter_psl_to_bed"
BED_FROM_PSL = FILTER_PSL_TO_BED_OUT+"/filtered_bed_from_psl.bed"
TX_LENGTH_VS_HITS = FILTER_PSL_TO_BED_OUT+"/tx_length_vs_hits.png"
FILTERED_TX_DATA = FILTER_PSL_TO_BED_OUT+"/filtered_tx_data.csv"

config["FILTERED_TX_DATAFRAME"] = pd.DataFrame()

# ---
rule filter_psl_to_bed:
    params:
        bed_from_psl_coverage=FILTER_PSL_TO_BED["BED_FROM_PSL_COVERAGE"],
        bed_from_psl_qsize=FILTER_PSL_TO_BED["BED_FROM_PSL_QSIZE"]
    input:
        psl=PSL,
        orthos=rules.annotations_via_fasta.output.annotations_xls
    output:
        bed_from_psl=BED_FROM_PSL,
        tx_length_vs_hits=TX_LENGTH_VS_HITS,
        filtered_tx_data=FILTERED_TX_DATA,

    script:
        "python/scripts/filter_psl_to_bed.py"



MAIN.append(rules.filter_psl_to_bed.output)

# ------------------------- #
#### SUBTRACT_GENE_MODELS ####
SUBTRACT_GENE_MODELS = config["SUBTRACT_GENE_MODELS"]
GENE_MODELS_BED = SUBTRACT_GENE_MODELS["GENE_MODELS_BED"]

SUBTRACT_GENE_MODELS_OUT = OUT_DIR_MAIN+"/subtract_gene_models"
GENE_MODEL_SUBTRACTED = SUBTRACT_GENE_MODELS_OUT+"/gene_model_subtracted.bed"

# ---
rule subtract_gene_models:
    input:
        gene_models_bed=GENE_MODELS_BED,
        bed_from_psl=rules.filter_psl_to_bed.output.bed_from_psl,

    output:
        gene_model_subtracted=GENE_MODEL_SUBTRACTED,

    script:
        "python/scripts/subtract_gene_models.py"



MAIN.append(rules.subtract_gene_models.output)

# ------------------------- #
#### MAKE_SNP_BEDS ####
MAKE_SNP_BEDS = config["MAKE_SNP_BEDS"]
SCAFFOLD_NAME_MAP = MAKE_SNP_BEDS["SCAFFOLD_NAME_MAP"]
DO_CLEANING = MAKE_SNP_BEDS["DO_CLEANING"]
SNP_FILES = MAKE_SNP_BEDS["SNP_FILES"]
P_THRESH = MAKE_SNP_BEDS["P_THRESH"]

MAKE_SNP_BEDS_OUT = OUT_DIR_MAIN+"/make_snp_beds"
SNP_BEDS = ["{path}/{basename}.bed".format(path=MAKE_SNP_BEDS_OUT, basename=os.path.splitext(os.path.basename(x))[0]) for x in SNP_FILES]

# ---
rule make_snp_beds:
    params:
        do_cleaning=DO_CLEANING,
        scaffold_name_map=SCAFFOLD_NAME_MAP,
        p_thresh=P_THRESH
    input:
        SNP_FILES
    output:
        SNP_BEDS

    script:
        "python/scripts/make_snp_beds.py"



MAIN.append(rules.make_snp_beds.output)

# ------------------------- #
#### SORT_BED_FILES ####
# ---
rule sort_bed_files:
    input:
        snp_beds=rules.make_snp_beds.output,
        gene_model_subtracted=rules.subtract_gene_models.output.gene_model_subtracted,
        gene_models_bed=rules.subtract_gene_models.input.gene_models_bed

    output:
        sorted_status=OUT_DIR_MAIN+"/sort_bed_files/sorted_status"

    shell:
        """

            for i in {input}
            do
                sort -k1,1 -k2,2n $i > $i.tmp
                mv $i.tmp $i
            done
            touch {output}
        """

MAIN.append(rules.sort_bed_files.output)

# ------------------------- #
#### GET_NEAREST_K_FEATURES ####
GET_NEAREST_K_FEATURES = config["GET_NEAREST_K_FEATURES"]
K_NUMBER = GET_NEAREST_K_FEATURES["K"]

GET_NEAREST_K_FEATURES_OUT = OUT_DIR_MAIN+'/get_nearest_k_features'
NEAREST_FEATURES_BEDS = ["{path}/{basename}.nearest.bed".format(path=GET_NEAREST_K_FEATURES_OUT, basename=os.path.splitext(os.path.basename(x))[0]) for x in SNP_BEDS]

SNPS_IN_FEATURES = ["{path}/{basename}.snps_in_features.xls".format(path=GET_NEAREST_K_FEATURES_OUT, basename=os.path.splitext(os.path.basename(x))[0]) for x in SNP_BEDS]

# ---
rule get_nearest_k_features:
    params:
        k_number=K_NUMBER
    input:
        snp_beds=rules.make_snp_beds.output,
        gene_model_subtracted=rules.subtract_gene_models.output.gene_model_subtracted,
        gene_models=rules.subtract_gene_models.input.gene_models_bed

    output:
        nearest_features_beds=NEAREST_FEATURES_BEDS,
        snps_in_features=SNPS_IN_FEATURES,

    script:
        "python/scripts/get_nearest_k_features.py"



MAIN.append(rules.get_nearest_k_features.output)

# ------------------------- #
#### MAKE_ID_TABLE_NO_DIFF_EXPR ####
MAKE_ID_TABLE_NO_DIFF_EXPR = config["MAKE_ID_TABLE_NO_DIFF_EXPR"]

CUFFCMP_TRACKING = MAKE_ID_TABLE_NO_DIFF_EXPR["CUFFCMP_TRACKING"]
ORTHOLOG_TABLE = ANNOTATIONS_VIA_FASTA["ORTHOLOG_TABLE"]

MAKE_ID_TABLE_NO_DIFF_EXPR_OUT = OUT_DIR_MAIN+'/make_id_table_no_diff_expr'
IDS_NO_DIFF_EXPR = MAKE_ID_TABLE_NO_DIFF_EXPR_OUT+'/ids_no_diff_expr.csv'

# ---
rule make_id_table_no_diff_expr:
    input:
        cuffcmp_tracking=CUFFCMP_TRACKING,
        ortholog_table=ORTHOLOG_TABLE,

    output:
        ids_no_diff_expr=IDS_NO_DIFF_EXPR,

    script:
        "python/scripts/make_id_table_no_diff_expr.py"



MAIN.append(rules.make_id_table_no_diff_expr.output)

# ------------------------- #
#### MAKE_ID_TABLE_WITH_DIFF_EXPR ####
MAKE_ID_TABLE_WITH_DIFF_EXPR = config["MAKE_ID_TABLE_WITH_DIFF_EXPR"]

IDS_NO_DIFF_EXPR = IDS_NO_DIFF_EXPR

EDGER_RESULTS_LABELS = list(MAKE_ID_TABLE_WITH_DIFF_EXPR["EDGER_RESULTS_INFO"].keys())
CUFFLINKS_RESULTS_LABELS = list(MAKE_ID_TABLE_WITH_DIFF_EXPR["CUFFLINKS_RESULTS_INFO"].keys())

EDGER_RESULTS = [MAKE_ID_TABLE_WITH_DIFF_EXPR["EDGER_RESULTS_INFO"][label] for label in EDGER_RESULTS_LABELS]
CUFFLINKS_RESULTS = [MAKE_ID_TABLE_WITH_DIFF_EXPR["CUFFLINKS_RESULTS_INFO"][label] for label in CUFFLINKS_RESULTS_LABELS]

MAKE_ID_TABLE_WITH_DIFF_EXPR_OUT = OUT_DIR_MAIN+'/make_id_table_with_diff_expr'
IDS_WITH_DIFF_EXPR = MAKE_ID_TABLE_WITH_DIFF_EXPR_OUT+'/ids_with_diff_expr.csv'



# ---
rule make_id_table_with_diff_expr:
    params:
        edger_results_labels=EDGER_RESULTS_LABELS,
        cufflinks_results_labels=CUFFLINKS_RESULTS_LABELS,
    input:
        edger_results=EDGER_RESULTS,
        cufflinks_results=CUFFLINKS_RESULTS,
        ids_no_diff_expr=IDS_NO_DIFF_EXPR,

    output:
        ids_with_diff_expr=IDS_WITH_DIFF_EXPR,

    script:
        "python/scripts/make_id_table_with_diff_expr.py"



MAIN.append(rules.make_id_table_with_diff_expr.output)

# ------------------------- #
#### SNPS_NEAR_HOMOLOGOUS_DE ####
SNPS_NEAR_HOMOLOGOUS_DE = config["SNPS_NEAR_HOMOLOGOUS_DE"]

# params
SNP_DISTANCE_FROM_GENE = SNPS_NEAR_HOMOLOGOUS_DE["SNP_DISTANCE_FROM_GENE"]
DE_TABLE_CHUNKSIZE = SNPS_NEAR_HOMOLOGOUS_DE["DE_TABLE_CHUNKSIZE"]
RUN_PARALLEL = SNPS_NEAR_HOMOLOGOUS_DE["RUN_PARALLEL"]
GENOME_BROWSER_URL = SNPS_NEAR_HOMOLOGOUS_DE["GENOME_BROWSER_URL"]


# input
NEAREST_FEATURES_BEDS = NEAREST_FEATURES_BEDS
IDS_WITH_DIFF_EXPR = IDS_WITH_DIFF_EXPR
GENE_MODELS_BED = SUBTRACT_GENE_MODELS["GENE_MODELS_BED"]

# output

SNPS_NEAR_HOMOLOGOUS_DE_OUT = OUT_DIR_MAIN+'/snps_near_homologous_de'
SNPS_NEAR_HOMOLOGOUS_DE_PATH = SNPS_NEAR_HOMOLOGOUS_DE_OUT+'/snps_near_homologous_de_distance_{distance}.csv'.format(distance=SNP_DISTANCE_FROM_GENE)



# ---
rule snps_near_homologous_de:
    params:
        snp_distance_from_gene=SNP_DISTANCE_FROM_GENE,
        de_table_chunksize=DE_TABLE_CHUNKSIZE,
        run_parallel=RUN_PARALLEL,
        genome_browser_url=GENOME_BROWSER_URL,

    input:
        nearest_features_beds=NEAREST_FEATURES_BEDS,
        ids_with_diff_expr=IDS_WITH_DIFF_EXPR,
        gene_models_bed=GENE_MODELS_BED,

    output:
        snps_near_homologous_de_path=SNPS_NEAR_HOMOLOGOUS_DE_PATH,

    shell:
        """python python/scripts/snps_near_homologous_de.py \
        {input.nearest_features_beds} \
        {input.ids_with_diff_expr} \
        {output.snps_near_homologous_de_path} \
        --distance {params.snp_distance_from_gene} \
        --chunksize {params.de_table_chunksize} \
        --parallel {params.run_parallel} \
        --url '{params.genome_browser_url}' \
        --bed {input.gene_models_bed} \
        """



MAIN.append(rules.snps_near_homologous_de.output)

# ------------------------- #


#### MAIN ####
# ---
rule main:
    input:
        MAIN
