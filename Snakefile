import os
import pandas as pd

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]


##### Config processing #####
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

#### Reference info processing

#### Setting reference from lib_ROI
if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]

#### Setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"))
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

#### FOLDERS
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

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
if config["use_cnmops"]:
    used_SV_callers.append("cnmops")
if config["use_panel_cnmops"]:
    used_SV_callers.append("panel_cnmops")
if config["use_exome_depth"]:
    used_SV_callers.append("exome_depth")


wildcard_constraints:
     tumor_normal = "tumor|normal|sample",


####################################
# SEPARATE RULES
include: "rules/cnvkit.smk"
include: "rules/gatk_cnv.smk"
include: "rules/jabCoNtool.smk"
include: "rules/control_freec.smk"
include: "rules/variant_postprocessing.smk"
include: "rules/common_prep.smk"




####################################
# RULE ALL
def all_inputs(wildcards):
    input_dict = {}
    input_dict["final_report"] = "CNV_varcalls/final_CNV_visualization.html",
    if config["create_cohort_data"] == True:
        input_dict["cohort_data_update_tag"] = "cohort_data/cohort_data_updated"
    return input_dict

rule all:
    input: unpack(all_inputs)
