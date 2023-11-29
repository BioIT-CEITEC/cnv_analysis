######################################
# wrapper for rule: create_cohort_data
######################################
import os
import subprocess
from snakemake.shell import shell
import re

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: create_cohort_data \n##\n")
f.close()

shell.executable("/bin/bash")

if hasattr(snakemake.input, "previous_cohort_data"):
    command = "tar -xzf " + snakemake.input.previous_cohort_data + " -C cohort_data/ "
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
else:
    command = "mkdir -p cohort_data/cohort_data/ "
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

if hasattr(snakemake.input, "cnvkit_normal_coverage_inputs"):
    command = "mkdir -p cohort_data/cohort_data/cnvkit "
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    f = open(log_filename, 'at')
    for filename in snakemake.input.cnvkit_normal_coverage_inputs:
        out_filename = filename.replace("CNV_varcalls/","")
        out_filename = out_filename.replace("/cnvkit/", "_")
        command = "cp " + filename + " cohort_data/cohort_data/cnvkit/" + out_filename
        f.write("## COMMAND: " + command + "\n")
        shell(command)
    f.close()


if hasattr(snakemake.input, "jabCoNtool_cohort_info"):
    command = "mkdir -p cohort_data/cohort_data/jabCoNtool "
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "cp " + snakemake.input.jabCoNtool_cohort_info + " cohort_data/cohort_data/jabCoNtool/cohort_info_tab.tsv"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

command = "cd cohort_data/; tar -czvf cohort_cnv_info.tar.gz cohort_data >> ../" + log_filename + " 2>&1; cd ../ "
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

command = "touch " + snakemake.output.update_finished_checkfile
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)
