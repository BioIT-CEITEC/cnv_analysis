import os
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]


#### Reference info processing
if config["lib_ROI"] == "rna":
    config["lib_ROI"] = "wgs"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

config = BR.load_organism()

#### FOLDERS
reference_directory = config["reference_dir"]

# ####################################
# # VARIALBES FROM CONFIG
used_SV_callers = []
used_CNV_callers = []
if config["use_gatk_cnv"]:
    used_CNV_callers.append("gatk_cnv")
if config["use_cnvkit"]:
    used_CNV_callers.append("cnvkit")
if config["use_jabCoNtool"]:
    used_CNV_callers.append("jabCoNtool")
if config["use_control_freec"]:
    used_CNV_callers.append("control_freec")
if config["use_manta"]:
    used_SV_callers.append("manta")
if config["use_gridss"]:
    used_SV_callers.append("gridss")


wildcard_constraints:
     tumor_normal = "tumor|normal|sample",


####################################
# SEPARATE RULES
include: "rules/cnvkit.smk"
include: "rules/gatk_cnv.smk"
include: "rules/jabCoNtool.smk"
include: "rules/control_freec.smk"
include: "rules/manta.smk"
include: "rules/delly.smk"
# include: "rules/gridss.smk"
include: "rules/svdb.smk"
include: "rules/variant_postprocessing.smk"
include: "rules/common_prep.smk"
# include: "rules/vep.smk"



####################################
# RULE ALL
def all_inputs(wildcards):
    input_dict = {}
    input_dict["final_report"] = "reports/final_SV_report.html",
    if config["create_cohort_data"] == True:
        input_dict["cohort_data_update_tag"] = "cohort_data/cohort_data_updated"
    return input_dict

rule all:
    input: unpack(all_inputs)
