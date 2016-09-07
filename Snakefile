"""Snakemake file."""

# See tutorial at: http://tiny.cc/snakemake_tutorial
from pathlib import Path

import os

import yaml

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

import munch

from python.functions import *

def pathify_by_key_ends(dictionary):
    """Return a dict that has had all values with keys marked as '*_PATH' or '*_DIR' converted to Path() instances."""
    for key, value in dictionary.items():
        if isinstance(value, dict):
            pathify_by_key_ends(value)
        elif key.endswith("_PATH") or key.endswith("_DIR"):
            dictionary[key] = Path(value)

    return dictionary


class MyRun(object):

    """Initialize and manage information common to the whole run."""

    def __init__(self, cfg):
        """Initialize common information for a run."""
        assert isinstance(cfg, dict)

        common = cfg["COMMON"]

        self.globals = munch.Munch()
        self.cfg = cfg
        self.name = common["RUN_NAME"]
        self.d = common["SHARED"]
        self.out_dir = Path("{base_dir}/{run_name}".format(base_dir=common["OUT_DIR"],
                                                           run_name=self.name
                                                           )
                            )
        self.log_dir = self.out_dir / "logs"

class MyRule(object):

    """Manage the initialization and deployment of rule-specific information."""

    def __init__(self, run, name):
        """Initialize logs, inputs, outputs, params, etc for a single rule."""
        assert isinstance(run, MyRun)

        self.run = run
        self.name = name.lower()
        self.log_dir = run.log_dir / self.name
        self.log = self.log_dir / "{name}.log".format(name=self.name)
        self.out_dir = run.out_dir / self.name
        self.i = munch.Munch() # inputs
        self.o = munch.Munch() # outputs
        self.p = munch.Munch() # params

        self._import_config_dict()

    def _import_config_dict(self):
        """Inport configuration values set for this rule so they are directly accessable as attributes."""
        try:
            for key, val in self.run.cfg[self.name.upper()].items():
                self.__setattr__(key, val)
            self.cfg = True
        except KeyError:
            self.cfg = False



#### COMMON RUN STUFF ####
ORIGINAL_CONFIG_AS_STRING = yaml.dump(config, default_flow_style=False)
config_ = pathify_by_key_ends(config)
config_ = munch.munchify(config)

RUN = MyRun(cfg=config_)

RUN_NAME = config["COMMON"]["RUN_NAME"]
OUT_DIR = "{base_dir}/{run_name}".format(base_dir=config["COMMON"]["OUT_DIR"], run_name=RUN_NAME)
OUT_DIR_MAIN = "{out_base}".format(out_base=OUT_DIR)
OUT_DIR_PREP = "{out_base}".format(out_base=OUT_DIR)


# collect meta-target 'inputs'
PREP = []
MAIN = []

############ BEGIN PIPELINE RULES ############


#### SAVE_RUN_CONFIG ####
SAVE_RUN_CONFIG = MyRule(run=RUN, name="SAVE_RUN_CONFIG")
SAVE_RUN_CONFIG_OUT = str(SAVE_RUN_CONFIG.out_dir)

rule save_run_config:
    input:
    output:
        file=SAVE_RUN_CONFIG_OUT

    run:
        with open(output.file, 'w') as cnf_out:
            cnf_out.write(ORIGINAL_CONFIG_AS_STRING)

PREP.append(rules.save_run_config.output)
MAIN.append(rules.save_run_config.output)


# ------------ prep ------------- #
#### MAKE_OOC ####
MAKE_OOC = MyRule(run=RUN, name="MAKE_OOC")
#MAKE_OOC = config["MAKE_OOC"]
GENOME = MAKE_OOC.GENOME
OLIGO_LEN = MAKE_OOC.OLIGO_LEN


OOC = str(MAKE_OOC.out_dir / "{OLIGO_LEN}.ooc".format(OLIGO_LEN=OLIGO_LEN))

# ---
rule make_ooc:
    params:
        oligo_len=OLIGO_LEN
    input:
        genome=GENOME,
    output:
        ooc=OOC

    shell:
        "blat {input.genome} "
        "/dev/null /dev/null "
        "-tileSize={params.oligo_len} -makeOoc={output.ooc}"

PREP.append(rules.make_ooc.output)

# ------------ prep ------------- #
#### RUN_BLAT ####
RUN_BLAT = MyRule(run=RUN, name="RUN_BLAT")
#RUN_BLAT = config["RUN_BLAT"]
TRANSCRIPTS = str(RUN_BLAT.TRANSCRIPTS)
PSL_NAME = RUN_BLAT.PSL_NAME


PSL = str(RUN_BLAT.out_dir / "{psl_name}".format(psl_name=PSL_NAME))

# ---
rule run_blat:
    params:
        oligo_len=OLIGO_LEN
    input:
        genome=GENOME,
        transcripts=TRANSCRIPTS,
        ooc=rules.make_ooc.output.ooc
    output:
        psl=PSL

    shell:
        "blat {input.genome} {input.transcripts} "
        "-q=dna -t=dna -ooc={input.ooc} "
        "{output.psl} "

PREP.append(rules.run_blat.output)

# # ------------ prep ------------- #
# #### GTF_JUST_EXONS ####
# GTF_JUST_EXONS = config["GTF_JUST_EXONS"]
# GTF_PATH = GTF_JUST_EXONS["GTF_PATH"]
#
# GTF_BASE_NAME = os.path.splitext(os.path.basename(GTF_PATH))
# GTF_EXONS = OUT_DIR_PREP+"/{base_name}.exons.gtf".format(base_name=GTF_BASE_NAME[0])
#
# # ---
# rule gtf_just_exons:
#     input:
#         gtf_path=GTF_PATH,
#     output:
#         gtf_exons=GTF_EXONS,
#
#     shell:
#         """awk '/\texon\t/' \
#         < {input.gtf_path} | \
#         sort -k 1,1 -k 4,4n > \
#         {output.gtf_exons}
#         """

# PREP.append(rules.gtf_just_exons.output)

# ------------ prep ------------- #
#### GTF_TO_BED ####
GTF_TO_BED = MyRule(run=RUN, name="GTF_TO_BED")
#GTF_TO_BED = config["GTF_TO_BED"]
GTF = str(GTF_TO_BED.GTF)

BED_BASE_NAME = os.path.splitext(os.path.basename(GTF))

BED = str(GTF_TO_BED.out_dir / "{base_name}.bed".format(base_name=BED_BASE_NAME[0]))
GTF_DB =  str(GTF_TO_BED.out_dir / "{base_name}.gtf.db".format(base_name=BED_BASE_NAME[0]))

# ---
rule gtf_to_bed:
    input:
        gtf=GTF,
    output:
        bed=BED,
        gtf_db=GTF_DB,

    script:
        "python/scripts/gtf_to_bed.py"

PREP.append(rules.gtf_to_bed.output)

# ------------ prep ------------- #
#### SORT_GENES_BED ####
SORT_GENES_BED = MyRule(run=RUN, name="SORT_GENES_BED")
GTF = str(GTF_TO_BED.GTF)

BED_SORTED = str(SORT_GENES_BED.out_dir / "{base_name}.sorted.bed".format(base_name=BED_BASE_NAME[0]))



# ---
rule sort_genes_bed:
    input:
        bed=rules.gtf_to_bed.output.bed,
    output:
        bed_sorted=BED_SORTED,

    shell:
        "sort -k 1,1 -k 2,2n "
        "{input.bed} > "
        "{output.bed_sorted}"


PREP.append(rules.sort_genes_bed.output)


# ------------ main ------------- #
#### ANNOTATIONS_VIA_FASTA ####
ANNOTATIONS_VIA_FASTA = MyRule(run=RUN, name="ANNOTATIONS_VIA_FASTA")
#ANNOTATIONS_VIA_FASTA = config["ANNOTATIONS_VIA_FASTA"]
ANNOTATIONS_VIA_FASTA_OUT = str(ANNOTATIONS_VIA_FASTA.out_dir)
ANNOTATIONS_XLS = ANNOTATIONS_VIA_FASTA_OUT+"/annotations_via_fasta.xls"
TX_FASTA = ANNOTATIONS_VIA_FASTA.TX_FASTA
ORTHOLOG_TABLE = ANNOTATIONS_VIA_FASTA.ORTHOLOG_TABLE

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

# ------------ main ------------- #
#### FILTER_PSL_TO_BED ####
FILTER_PSL_TO_BED = MyRule(run=RUN, name="FILTER_PSL_TO_BED")
#FILTER_PSL_TO_BED = config["FILTER_PSL_TO_BED"]
PSL = PSL # defined in meta-target: PREP

FILTER_PSL_TO_BED_OUT = str(FILTER_PSL_TO_BED.out_dir)
BED_FROM_PSL = FILTER_PSL_TO_BED_OUT+"/filtered_bed_from_psl.bed"
TX_LENGTH_VS_HITS = FILTER_PSL_TO_BED_OUT+"/tx_length_vs_hits.png"
FILTERED_TX_DATA = FILTER_PSL_TO_BED_OUT+"/filtered_tx_data.csv"

config["FILTERED_TX_DATAFRAME"] = pd.DataFrame()

# ---
rule filter_psl_to_bed:
    params:
        bed_from_psl_coverage=FILTER_PSL_TO_BED.BED_FROM_PSL_COVERAGE,
        bed_from_psl_qsize=FILTER_PSL_TO_BED.BED_FROM_PSL_QSIZE
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

# ------------ main ------------- #
#### SUBTRACT_GENE_MODELS ####
SUBTRACT_GENE_MODELS = MyRule(run=RUN, name="SUBTRACT_GENE_MODELS")
#SUBTRACT_GENE_MODELS = config["SUBTRACT_GENE_MODELS"]
GENE_MODELS_BED = rules.sort_genes_bed.output.bed_sorted

SUBTRACT_GENE_MODELS_OUT = str(SUBTRACT_GENE_MODELS.out_dir)
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


# ------------ main ------------- #
#### MAKE_SNP_BEDS ####
MAKE_SNP_BEDS = MyRule(run=RUN, name="MAKE_SNP_BEDS")
#MAKE_SNP_BEDS = config["MAKE_SNP_BEDS"]
SCAFFOLD_NAME_MAP = MAKE_SNP_BEDS.SCAFFOLD_NAME_MAP
DO_CLEANING = MAKE_SNP_BEDS.DO_CLEANING
SNP_FILES = [Path(snp_file) for snp_file in MAKE_SNP_BEDS.SNP_FILES]
P_THRESH = MAKE_SNP_BEDS.P_THRESH

MAKE_SNP_BEDS_OUT = str(MAKE_SNP_BEDS.out_dir)
SNP_BEDS = ["{path}/{basename}.bed".format(path=MAKE_SNP_BEDS_OUT, basename=f.stem) for f in SNP_FILES]
SNP_FILES = [str(f) for f in SNP_FILES]


# ---
rule make_snp_beds:
    params:
        do_cleaning=DO_CLEANING,
        scaffold_name_map=SCAFFOLD_NAME_MAP,
        p_thresh=P_THRESH
    input:
        snp_files=SNP_FILES
    output:
        snp_beds=SNP_BEDS

    script:
        "python/scripts/make_snp_beds.py"


MAIN.append(rules.make_snp_beds.output)

# ------------ main ------------- #
#### SORT_BED_FILES ####
SORT_BED_FILES = MyRule(run=RUN, name="SORT_BED_FILES")

# ---
rule sort_bed_files:
    input:
        snp_beds=rules.make_snp_beds.output,
        gene_model_subtracted=rules.subtract_gene_models.output.gene_model_subtracted,
        gene_models_bed=rules.subtract_gene_models.input.gene_models_bed

    output:
        sorted_status=str(SORT_BED_FILES.out_dir / "sorted_status")

    shell:
        "for i in {input}; "
        "do; "
        "    sort -k1,1 -k2,2n $i > $i.tmp; "
        "    mv $i.tmp $i; "
        "done; "
        "touch {output}"

MAIN.append(rules.sort_bed_files.output)

# ------------ main ------------- #
#### GET_NEAREST_K_FEATURES ####
GET_NEAREST_K_FEATURES = MyRule(run=RUN, name="GET_NEAREST_K_FEATURES")
#GET_NEAREST_K_FEATURES = config["GET_NEAREST_K_FEATURES"]
K_NUMBER = GET_NEAREST_K_FEATURES.K

GET_NEAREST_K_FEATURES_OUT = str(GET_NEAREST_K_FEATURES.out_dir)
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

# ------------ main ------------- #
#### MAKE_ID_TABLE_NO_DIFF_EXPR ####
MAKE_ID_TABLE_NO_DIFF_EXPR = MyRule(run=RUN, name="MAKE_ID_TABLE_NO_DIFF_EXPR")
#MAKE_ID_TABLE_NO_DIFF_EXPR = config["MAKE_ID_TABLE_NO_DIFF_EXPR"]

CUFFCMP_TRACKING = MAKE_ID_TABLE_NO_DIFF_EXPR.CUFFCMP_TRACKING
ORTHOLOG_TABLE = ANNOTATIONS_VIA_FASTA.ORTHOLOG_TABLE

MAKE_ID_TABLE_NO_DIFF_EXPR_OUT = str(MAKE_ID_TABLE_NO_DIFF_EXPR.out_dir)
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

# ------------ main ------------- #
#### MAKE_ID_TABLE_WITH_DIFF_EXPR ####
MAKE_ID_TABLE_WITH_DIFF_EXPR = MyRule(run=RUN, name="MAKE_ID_TABLE_WITH_DIFF_EXPR")
#MAKE_ID_TABLE_WITH_DIFF_EXPR = config["MAKE_ID_TABLE_WITH_DIFF_EXPR"]

IDS_NO_DIFF_EXPR = IDS_NO_DIFF_EXPR

EDGER_RESULTS_LABELS = list(MAKE_ID_TABLE_WITH_DIFF_EXPR.EDGER_RESULTS_INFO.keys())
CUFFLINKS_RESULTS_LABELS = list(MAKE_ID_TABLE_WITH_DIFF_EXPR.CUFFLINKS_RESULTS_INFO.keys())

EDGER_RESULTS = [MAKE_ID_TABLE_WITH_DIFF_EXPR.EDGER_RESULTS_INFO[label] for label in EDGER_RESULTS_LABELS]
CUFFLINKS_RESULTS = [MAKE_ID_TABLE_WITH_DIFF_EXPR.CUFFLINKS_RESULTS_INFO[label] for label in CUFFLINKS_RESULTS_LABELS]

MAKE_ID_TABLE_WITH_DIFF_EXPR_OUT = str(MAKE_ID_TABLE_WITH_DIFF_EXPR.out_dir)
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

# ------------ main ------------- #
#### SNPS_NEAR_HOMOLOGOUS_DE ####
SNPS_NEAR_HOMOLOGOUS_DE = MyRule(run=RUN, name="SNPS_NEAR_HOMOLOGOUS_DE")
#SNPS_NEAR_HOMOLOGOUS_DE = config["SNPS_NEAR_HOMOLOGOUS_DE"]

# params
SNP_DISTANCE_FROM_GENE = SNPS_NEAR_HOMOLOGOUS_DE.SNP_DISTANCE_FROM_GENE
DE_TABLE_CHUNKSIZE = SNPS_NEAR_HOMOLOGOUS_DE.DE_TABLE_CHUNKSIZE
RUN_PARALLEL = SNPS_NEAR_HOMOLOGOUS_DE.RUN_PARALLEL
GENOME_BROWSER_URL = SNPS_NEAR_HOMOLOGOUS_DE.GENOME_BROWSER_URL


# input
NEAREST_FEATURES_BEDS = NEAREST_FEATURES_BEDS
IDS_WITH_DIFF_EXPR = IDS_WITH_DIFF_EXPR
GENE_MODELS_BED = rules.sort_genes_bed.output.bed_sorted

# output

SNPS_NEAR_HOMOLOGOUS_DE_OUT = str(SNPS_NEAR_HOMOLOGOUS_DE.out_dir)
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
#### GET_GENO_R2_POS ####
GET_GENO_R2_POS = MyRule(run=RUN, name="GET_GENO_R2_POS")
VCF_POPS = [str(vcf_pop.stem) for vcf_pop in RUN.d.VCF_POPULATION_DIR.glob('*.vcf')]

# input
GET_GENO_R2_POS.i.snp_bed = str(MAKE_SNP_BEDS.out_dir / "{vcf_pop}.bed")

# output
GET_GENO_R2_POS.o.snp_pos_wldcd = str(GET_GENO_R2_POS.out_dir / "{vcf_pop}.pos")
GET_GENO_R2_POS.o.snp_pos_expanded = expand(GET_GENO_R2_POS.o.snp_pos_wldcd, vcf_pop=VCF_POPS)

# ---
rule get_geno_r2_pos:
    log:
        path=str(GET_GENO_R2_POS.log)

    input:
        snp_bed=GET_GENO_R2_POS.i.snp_bed,

    output:
        GET_GENO_R2_POS.o.snp_pos_wldcd,

    shell:
        "cut -d$'\t' -f1,3 {input.snp_bed} > {output}"
        # "&> {log.path}"

MAIN.append(GET_GENO_R2_POS.o.snp_pos_expanded)

# ------------------------- #
#### GET_LINKED_SNPS ####
GET_LINKED_SNPS = MyRule(run=RUN, name="GET_LINKED_SNPS")

# logs
GET_LINKED_SNPS.log = "{dir}/{{vcf_pop}}.log".format(dir=GET_LINKED_SNPS.log.parent)

# params
GET_LINKED_SNPS.p.ld_window_bp = GET_LINKED_SNPS.PARAMS.LD_WINDOW_BP
GET_LINKED_SNPS.p.min_r2 = GET_LINKED_SNPS.PARAMS.MIN_R2
GET_LINKED_SNPS.p.out_prefix = str(GET_LINKED_SNPS.out_dir / "{vcf_pop}")

# input
GET_LINKED_SNPS.i.vcf = str(GET_LINKED_SNPS.IN.VCF_DIR / "{vcf_pop}.vcf")
GET_LINKED_SNPS.i.snp_list = str(GET_GENO_R2_POS.out_dir / "{vcf_pop}.pos") # CHROM\tPOS

# output
GET_LINKED_SNPS.o.linked_snps_wldcd = GET_LINKED_SNPS.p.out_prefix + ".list.geno.ld"
GET_LINKED_SNPS.o.linked_snps_expanded = expand(GET_LINKED_SNPS.o.linked_snps_wldcd, vcf_pop=VCF_POPS)


# ---
rule get_linked_snps:
    log:
        path=str(GET_LINKED_SNPS.log),
    params:
        ld_window_bp=GET_LINKED_SNPS.p.ld_window_bp,
        min_r2=GET_LINKED_SNPS.p.min_r2,
        out_prefix=GET_LINKED_SNPS.p.out_prefix,

    input:
        vcf=GET_LINKED_SNPS.i.vcf,
        snp_list=GET_LINKED_SNPS.i.snp_list,

    output:
        linked_snps=GET_LINKED_SNPS.o.linked_snps_wldcd,

    shell:
        "vcftools --vcf {input.vcf} "
        # "--geno-r2 "
        "--geno-r2-positions {input.snp_list} "
        "--ld-window-bp {params.ld_window_bp} "
        "--min-r2 {params.min_r2} "
        "--out {params.out_prefix} "
        "&> {log.path}"

MAIN.append(GET_LINKED_SNPS.o.linked_snps_expanded)


# ------------ meta-targets ------------- #

#### PREP ####
# ---
rule prep:
    input:
        PREP

#### MAIN ####
# ---
rule main:
    input:
        MAIN

#### ALL ####
# ---
rule all:
    input:
        MAIN,
        PREP

#### PUT_PSL ####
psl_dir_ = "{run_dir}/run_blat".format(run_dir=str(RUN.out_dir))
# ---
rule put_psl:
    shell:
        "mkdir -p {psl_dir_}; "
        "tar --strip-components=3 -xvmf backup_data/gmm_transcripts_v_GfusI1_scaffold.dna2dna.psl.tar.gz -C {psl_dir_}; "
