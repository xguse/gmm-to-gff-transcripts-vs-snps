"""Snakemake file."""
# See tutorial at: http://tiny.cc/snakemake_tutorial
from pathlib import Path
from glob import glob
from pprint import pprint

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


VCF_POPS = [str(vcf_pop.stem) for vcf_pop in RUN.d.VCF_POPULATION_DIR.glob('*.vcf')]
SNP_SOURCES = ['original','linked']

# collect meta-target 'inputs'
PREP = []
PREP_SNPS = []
ANALYZE_SNPS = []
DEBUG = []


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
PREP_SNPS.append(rules.save_run_config.output)


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

# ------------ prep ------------- #
#### GTF_TO_BED ####
GTF_TO_BED = MyRule(run=RUN, name="GTF_TO_BED")
#GTF_TO_BED = config["GTF_TO_BED"]

# inputs
GTF_TO_BED.i.gtf = str(GTF_TO_BED.GTF_PATH)

# outputs
GTF_TO_BED.o.bed = str(GTF_TO_BED.out_dir / "{base_name}.bed".format(base_name=GTF_TO_BED.GTF_PATH.stem))
GTF_TO_BED.o.gtf_db = str(GTF_TO_BED.out_dir / "{base_name}.gtf.db".format(base_name=GTF_TO_BED.GTF_PATH.stem))



# ---
rule gtf_to_bed:
    input:
        gtf=GTF_TO_BED.i.gtf,
    output:
        bed=GTF_TO_BED.o.bed,
        gtf_db=GTF_TO_BED.o.gtf_db,

    script:
        "python/scripts/gtf_to_bed.py"

PREP.append(GTF_TO_BED.o.bed)
PREP.append(GTF_TO_BED.o.gtf_db)

# ------------ prep ------------- #
#### SORT_GENES_BED ####
SORT_GENES_BED = MyRule(run=RUN, name="SORT_GENES_BED")


SORT_GENES_BED.o.bed_sorted = str(SORT_GENES_BED.out_dir / "{base_name}.sorted.bed".format(base_name=GTF_TO_BED.GTF_PATH.stem))



# ---
rule sort_genes_bed:
    input:
        bed=GTF_TO_BED.o.bed,
    output:
        bed_sorted=SORT_GENES_BED.o.bed_sorted,

    shell:
        "sort -k 1,1 -k 2,2n "
        "{input.bed} > "
        "{output.bed_sorted}"


PREP.append(SORT_GENES_BED.o.bed_sorted)


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



PREP_SNPS.append(rules.annotations_via_fasta.output)

# ------------ main ------------- #
#### FILTER_PSL_TO_BED ####
FILTER_PSL_TO_BED = MyRule(run=RUN, name="FILTER_PSL_TO_BED")
#FILTER_PSL_TO_BED = config["FILTER_PSL_TO_BED"]
PSL = PSL # defined in meta-target: PREP


FILTER_PSL_TO_BED.o.bed_from_psl = str(FILTER_PSL_TO_BED.out_dir)+"/filtered_bed_from_psl.bed"
FILTER_PSL_TO_BED.o.tx_length_vs_hits = str(FILTER_PSL_TO_BED.out_dir)+"/tx_length_vs_hits.png"
FILTER_PSL_TO_BED.o.filtered_tx_data = str(FILTER_PSL_TO_BED.out_dir)+"/filtered_tx_data.csv"

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
        bed_from_psl=FILTER_PSL_TO_BED.o.bed_from_psl,
        tx_length_vs_hits=FILTER_PSL_TO_BED.o.tx_length_vs_hits,
        filtered_tx_data=FILTER_PSL_TO_BED.o.filtered_tx_data,

    script:
        "python/scripts/filter_psl_to_bed.py"



PREP_SNPS.append(FILTER_PSL_TO_BED.o.bed_from_psl)
PREP_SNPS.append(FILTER_PSL_TO_BED.o.tx_length_vs_hits)
PREP_SNPS.append(FILTER_PSL_TO_BED.o.filtered_tx_data)

# ------------ main ------------- #
#### SUBTRACT_GENE_MODELS ####
SUBTRACT_GENE_MODELS = MyRule(run=RUN, name="SUBTRACT_GENE_MODELS")
#SUBTRACT_GENE_MODELS = config["SUBTRACT_GENE_MODELS"]


SUBTRACT_GENE_MODELS.o.gene_model_subtracted = str(SUBTRACT_GENE_MODELS.out_dir)+"/gene_model_subtracted.bed"

# ---
rule subtract_gene_models:
    input:
        gene_models_bed=SORT_GENES_BED.o.bed_sorted,
        bed_from_psl=FILTER_PSL_TO_BED.o.bed_from_psl,

    output:
        gene_model_subtracted=SUBTRACT_GENE_MODELS.o.gene_model_subtracted,

    script:
        "python/scripts/subtract_gene_models.py"

PREP_SNPS.append(SUBTRACT_GENE_MODELS.o.gene_model_subtracted)


# ------------ main ------------- #
#### MAKE_SNP_BEDS ####
MAKE_SNP_BEDS = MyRule(run=RUN, name="MAKE_SNP_BEDS")
#MAKE_SNP_BEDS = config["MAKE_SNP_BEDS"]

# inputs
MAKE_SNP_BEDS.i.snp_files_wldcd = '{}'.format(MAKE_SNP_BEDS.ORIG_SNP_DIR) + '/{vcf_pop}.txt'
MAKE_SNP_BEDS.i.snp_files_expanded = expand(MAKE_SNP_BEDS.i.snp_files_wldcd,vcf_pop=VCF_POPS)

# outputs
# SNP_BEDS = ["{path}/{basename}.bed".format(path=MAKE_SNP_BEDS_OUT, basename=f.stem) for f in SNP_FILES]
MAKE_SNP_BEDS.o.snp_beds_wldcd = str(MAKE_SNP_BEDS.out_dir / '{vcf_pop}.original.bed')
MAKE_SNP_BEDS.o.snp_beds_expanded = expand(MAKE_SNP_BEDS.o.snp_beds_wldcd,vcf_pop=VCF_POPS)

# print(MAKE_SNP_BEDS.o.snp_beds_expanded)

# ---
rule make_snp_beds:
    params:
        do_cleaning=MAKE_SNP_BEDS.DO_CLEANING,
        scaffold_name_map=MAKE_SNP_BEDS.SCAFFOLD_NAME_MAP,
        p_thresh=MAKE_SNP_BEDS.P_THRESH
    input:
        snp_files=MAKE_SNP_BEDS.i.snp_files_wldcd
    output:
        snp_beds=MAKE_SNP_BEDS.o.snp_beds_wldcd

    script:
        "python/scripts/make_snp_beds.py"

PREP_SNPS.append(MAKE_SNP_BEDS.o.snp_beds_expanded)

# ------------ main ------------- #
#### SORT_BED_FILES ####
SORT_BED_FILES = MyRule(run=RUN, name="SORT_BED_FILES")

SORT_BED_FILES.o.sorted_status = str(SORT_BED_FILES.out_dir / "sorted_status")
# ---
rule sort_bed_files:
    input:
        snp_beds=MAKE_SNP_BEDS.o.snp_beds_expanded,
        gene_model_subtracted=SUBTRACT_GENE_MODELS.o.gene_model_subtracted,
        gene_models_bed=SORT_GENES_BED.o.bed_sorted

    output:
        sorted_status=SORT_BED_FILES.o.sorted_status

    shell:
        "for i in {input}; "
        "do; "
        "    sort -k1,1 -k2,2n $i > $i.tmp; "
        "    mv $i.tmp $i; "
        "done; "
        "touch {output}"

PREP_SNPS.append(SORT_BED_FILES.o.sorted_status)



# ------------------------- #
#### GET_GENO_R2_POS ####
GET_GENO_R2_POS = MyRule(run=RUN, name="GET_GENO_R2_POS")

# input
GET_GENO_R2_POS.i.snp_bed = str(MAKE_SNP_BEDS.out_dir / "{vcf_pop}.original.bed")

# output
GET_GENO_R2_POS.o.snp_pos_wldcd = str(GET_GENO_R2_POS.out_dir / "{vcf_pop}.original.pos")
GET_GENO_R2_POS.o.snp_pos_expanded = expand(GET_GENO_R2_POS.o.snp_pos_wldcd, vcf_pop=VCF_POPS)

# ---
rule get_geno_r2_pos:
    log:
        path=str(GET_GENO_R2_POS.log)

    input:
        snp_bed=GET_GENO_R2_POS.i.snp_bed,
        sorted_sentinel=SORT_BED_FILES.o.sorted_status

    output:
        GET_GENO_R2_POS.o.snp_pos_wldcd,

    shell:
        "cut -d$'\t' -f1,3 {input.snp_bed} > {output}"
        # "&> {log.path}"

PREP_SNPS.append(GET_GENO_R2_POS.o.snp_pos_expanded)

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
GET_LINKED_SNPS.i.snp_list = str(GET_GENO_R2_POS.out_dir / "{vcf_pop}.original.pos") # CHROM\tPOS

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

# DEBUG.append(GET_LINKED_SNPS.o.linked_snps_expanded)
PREP_SNPS.append(GET_LINKED_SNPS.o.linked_snps_expanded)

# ------------------------- #
#### MAKE_LINKED_SNPS_BEDS ####
MAKE_LINKED_SNPS_BEDS = MyRule(run=RUN, name="MAKE_LINKED_SNPS_BEDS")

# input
MAKE_LINKED_SNPS_BEDS.i.linked_snps_wldcd = GET_LINKED_SNPS.o.linked_snps_wldcd

# output
MAKE_LINKED_SNPS_BEDS.o.bed_wldcd = str(MAKE_LINKED_SNPS_BEDS.out_dir / "{vcf_pop}.linked.bed")
MAKE_LINKED_SNPS_BEDS.o.bed_expanded = expand(MAKE_LINKED_SNPS_BEDS.o.bed_wldcd, vcf_pop=VCF_POPS)

# ---
rule make_linked_snps_beds:
    log:
        path=str(MAKE_LINKED_SNPS_BEDS.log)

    input:
        snp_file=MAKE_LINKED_SNPS_BEDS.i.linked_snps_wldcd,

    output:
        snp_bed=MAKE_LINKED_SNPS_BEDS.o.bed_wldcd,

    script:
        "python/scripts/make_linked_snps_beds.py"

DEBUG.append(MAKE_LINKED_SNPS_BEDS.o.bed_expanded)
PREP_SNPS.append(MAKE_LINKED_SNPS_BEDS.o.bed_expanded)


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



ANALYZE_SNPS.append(rules.make_id_table_no_diff_expr.output)

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



ANALYZE_SNPS.append(rules.make_id_table_with_diff_expr.output)




# ------------ main ------------- #
#### GET_NEAREST_K_FEATURES ####
GET_NEAREST_K_FEATURES = MyRule(run=RUN, name="GET_NEAREST_K_FEATURES")


# outputs
GET_NEAREST_K_FEATURES.o.nearest_features_bed_wldcd = str(GET_NEAREST_K_FEATURES.out_dir / "{vcf_pop}.{snp_source}.nearest.bed")
GET_NEAREST_K_FEATURES.o.nearest_features_bed_expanded = expand(GET_NEAREST_K_FEATURES.o.nearest_features_bed_wldcd, vcf_pop=VCF_POPS, snp_source=SNP_SOURCES)
GET_NEAREST_K_FEATURES.o.snps_in_features_wldcd = str(GET_NEAREST_K_FEATURES.out_dir / "{vcf_pop}.{snp_source}.snps_in_features.xls")
GET_NEAREST_K_FEATURES.o.snps_in_features_expanded = expand(GET_NEAREST_K_FEATURES.o.snps_in_features_wldcd, vcf_pop=VCF_POPS, snp_source=SNP_SOURCES)


# ---
rule get_nearest_k_features:
    params:
        k_number=GET_NEAREST_K_FEATURES.PARAMS.K
    input:
        snp_beds=[str(path) for path in list(RUN.out_dir.glob('**/*.linked.bed'))+list(RUN.out_dir.glob('**/*.original.bed'))],
        gene_model_subtracted=SUBTRACT_GENE_MODELS.o.gene_model_subtracted,
        gene_models=SORT_GENES_BED.o.bed_sorted,
        sorted_sentinel=SORT_BED_FILES.o.sorted_status,
        linked_snp_beds=MAKE_LINKED_SNPS_BEDS.o.bed_expanded,

    output:
        nearest_features_beds=GET_NEAREST_K_FEATURES.o.nearest_features_bed_expanded,
        snps_in_features=GET_NEAREST_K_FEATURES.o.snps_in_features_expanded,

    script:
        "python/scripts/get_nearest_k_features.py"

ANALYZE_SNPS.append(GET_NEAREST_K_FEATURES.o.nearest_features_bed_expanded)
ANALYZE_SNPS.append(GET_NEAREST_K_FEATURES.o.snps_in_features_expanded)

# DEBUG.append(GET_NEAREST_K_FEATURES.o.nearest_features_bed_expanded)
# DEBUG.append(GET_NEAREST_K_FEATURES.o.snps_in_features_expanded)


# ------------ main ------------- #
#### SNPS_NEAR_HOMOLOGOUS_DE ####
SNPS_NEAR_HOMOLOGOUS_DE = MyRule(run=RUN, name="SNPS_NEAR_HOMOLOGOUS_DE")
#SNPS_NEAR_HOMOLOGOUS_DE = config["SNPS_NEAR_HOMOLOGOUS_DE"]

# params
SNPS_NEAR_HOMOLOGOUS_DE.p.snp_distance_from_gene = SNPS_NEAR_HOMOLOGOUS_DE.SNP_DISTANCE_FROM_GENE
SNPS_NEAR_HOMOLOGOUS_DE.p.de_table_chunksize = SNPS_NEAR_HOMOLOGOUS_DE.DE_TABLE_CHUNKSIZE
SNPS_NEAR_HOMOLOGOUS_DE.p.run_parallel = SNPS_NEAR_HOMOLOGOUS_DE.RUN_PARALLEL
SNPS_NEAR_HOMOLOGOUS_DE.p.genome_browser_url = SNPS_NEAR_HOMOLOGOUS_DE.GENOME_BROWSER_URL

# input
## set vcf_pop manually but leave kind as wildcard
nearest_features_bed_wldcd_mod_template = GET_NEAREST_K_FEATURES.o.nearest_features_bed_wldcd.replace('{snp_source}','{{snp_source}}')

SNPS_NEAR_HOMOLOGOUS_DE.i.nearest_features_bed_list_wldcds = [nearest_features_bed_wldcd_mod_template.format(vcf_pop=pop) for pop in VCF_POPS]

# output
SNPS_NEAR_HOMOLOGOUS_DE.o.snps_near_homologous_de_path_wldcd = str(SNPS_NEAR_HOMOLOGOUS_DE.out_dir / 'snps_near_homologous_de_distance_{distance}.{{snp_source}}.csv'.format(distance=SNPS_NEAR_HOMOLOGOUS_DE.p.snp_distance_from_gene))

SNPS_NEAR_HOMOLOGOUS_DE.o.snps_near_homologous_de_path_expanded = expand(SNPS_NEAR_HOMOLOGOUS_DE.o.snps_near_homologous_de_path_wldcd,
                                                                         snp_source=SNP_SOURCES)

# ---
rule snps_near_homologous_de:
    params:
        snp_distance_from_gene=SNPS_NEAR_HOMOLOGOUS_DE.p.snp_distance_from_gene,
        de_table_chunksize=SNPS_NEAR_HOMOLOGOUS_DE.p.de_table_chunksize,
        run_parallel=SNPS_NEAR_HOMOLOGOUS_DE.p.run_parallel,
        genome_browser_url=SNPS_NEAR_HOMOLOGOUS_DE.p.genome_browser_url,

    input:
        nearest_features_beds=SNPS_NEAR_HOMOLOGOUS_DE.i.nearest_features_bed_list_wldcds,
        ids_with_diff_expr=IDS_WITH_DIFF_EXPR,
        gene_models_bed=SORT_GENES_BED.o.bed_sorted,

    output:
        snps_near_homologous_de_path=SNPS_NEAR_HOMOLOGOUS_DE.o.snps_near_homologous_de_path_wldcd,

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

ANALYZE_SNPS.append(SNPS_NEAR_HOMOLOGOUS_DE.o.snps_near_homologous_de_path_expanded)










# ------------ meta-targets ------------- #

#### DEBUG ####
# ---
rule debug:
    input:
        DEBUG

#### PREP ####
# ---
rule prep:
    input:
        PREP

#### PREP_SNPS ####
# ---
rule prep_snps:
    input:
        PREP_SNPS

#### ANALYZE_SNPS ####
# ---
rule analyze_snps:
    input:
        ANALYZE_SNPS


#### PUT_PSL ####
psl_dir_ = "{run_dir}/run_blat".format(run_dir=str(RUN.out_dir))
# ---
rule put_psl:
    shell:
        "mkdir -p {psl_dir_}; "
        "tar --strip-components=3 -xvmf backup_data/gmm_transcripts_v_GfusI1_scaffold.dna2dna.psl.tar.gz -C {psl_dir_}; "
        "touch {psl_dir_}/gmm_transcripts_v_GfusI1_scaffold.dna2dna.psl"

# #### SAFE_CLEAN ####
#
# # ---
# rule safe_clean:
#     shell:
#         "rm -rf; "
#         "tar --strip-components=3 -xvmf backup_data/gmm_transcripts_v_GfusI1_scaffold.dna2dna.psl.tar.gz -C {psl_dir_}; "
#         "touch {psl_dir_}/gmm_transcripts_v_GfusI1_scaffold.dna2dna.psl"
#
