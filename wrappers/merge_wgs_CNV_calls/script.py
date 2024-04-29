#############################################################
# wrapper for rule: merge_CNVs
#############################################################
import os
import sys
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: merge_CNVs \n##\n")
f.close()

shell.executable("/bin/bash")

CNV_call_file_list = []

if hasattr(snakemake.input, "gatk_cnv_variants"):
    CNV_call_file_list.append("gatk_cnv")
    CNV_call_file_list.extend(snakemake.input.gatk_cnv_variants)
    CNV_call_file_list.append("gatk_cnv_end")

if hasattr(snakemake.input, "jabCoNtool_variants"):
    CNV_call_file_list.append("jabCoNtool")
    CNV_call_file_list.append(snakemake.input.jabCoNtool_variants)
    CNV_call_file_list.append("jabCoNtool_end")

if hasattr(snakemake.input, "control_freec_variants"):
    CNV_call_file_list.append("control_freec")
    CNV_call_file_list.extend(snakemake.input.control_freec_variants)
    CNV_call_file_list.append("control_freec_end")

if len(CNV_call_file_list) > 0:
    command = ("Rscript "+os.path.abspath(os.path.dirname(__file__))+"/process_and_format_CNV.R" \
              + " " + snakemake.output.all_vars_tsv \
              + " " + str(snakemake.params.overlap) \
              + " " + str(snakemake.input.region_bed) \
              + " " + str(snakemake.params.library_type) \
              + " '" + str(snakemake.params.direct_sample_CNV_compare) + "'" \
              + " " + " ".join(CNV_call_file_list) \
              + " >> " + log_filename + " 2>&1")

    f = open(log_filename + "_Rargs", 'w')
    f.write(" ".join(command.split(" ")[2:-3]) + "\n")
    f.close()

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
    f.close()

    shell(command)
else:
    print >> sys.stderr, "No callers set for a sample."
    sys.exit(1)


