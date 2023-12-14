#############################################################
# wrapper for rule: cnv_computation
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: cnv_computation \n##\n")
f.close()

shell.executable("/bin/bash")

# if hasattr(snakemake.input, "snp_bed"):
#     panel_snps_filename = snakemake.input.snp_bed
#     snp_AF_files = "snp_AF " + " ".join(snakemake.input.snp_AF)
#     if snakemake.params.calling_type == "tumor_normal":
#         normal_snp_AF_files = "normal_snp_AF " + " ".join(snakemake.input.normal_snp_AF)
# else:
#     panel_snps_filename = "no_use_snps"
#     snp_AF_files = ""
#     normal_snp_AF_files = ""


if hasattr(snakemake.input, "cohort_data"):
    previous_cohort_data = snakemake.input.cohort_data
else:
    previous_cohort_data = "no_previous_cohort_data"

if hasattr(snakemake.input, "snp_bed"):
    panel_snps_filename = snakemake.input.snp_bed
else:
    panel_snps_filename = "no_use_snps"

if hasattr(snakemake.input, "GC_profile_file"):
    GC_normalization_file = snakemake.input.GC_profile_file
else:
    GC_normalization_file = "no_GC_norm"

if hasattr(snakemake.input, "cytoband_file"):
    cytoband_file = snakemake.input.cytoband_file
else:
    cytoband_file = "no_cytoband"

if snakemake.params.lib_ROI == "wgs":
    library_type = "wgs"
else:
    library_type = "panel"

default_norm_cov_text = snakemake.params.default_norm_cov

if snakemake.params.calling_type == "tumor_normal":
    norm_cov_sample_params = " norm_cov " + " ".join(snakemake.input.normal_sample_cov)
else:
    norm_cov_sample_params = ""

    command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/jabConTool_main.R" \
                    + " " + snakemake.output.all_res_prob_tab \
                    + " " + snakemake.input.region_bed\
                    + " " + panel_snps_filename \
                    + " " + snakemake.params.calling_type \
                    + " " + library_type \
                    + " " + GC_normalization_file \
                    + " " + cytoband_file \
                    + " " + previous_cohort_data \
                    + " " + str(snakemake.params.jabCoNtool_predict_TL) \
                    + " " + str(snakemake.params.max_CNV_occurance_in_cohort) \
                    + " '" + ';'.join(str(element) for element in snakemake.params.default_norm_cov) + "' " \
                    + " cov " + " ".join(snakemake.input.sample_cov)\
                    + norm_cov_sample_params \
                    + " >> " + log_filename + " 2>&1"


f = open(log_filename + "_Rargs", 'w')
f.write(" ".join(command.split(" ")[2:-3]) + "\n")
f.close()

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()
shell(command)

# command = "touch " + snakemake.output.all_res_prob_tab


