######################################
# wrapper for rule: per_sample_coverage_counting
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: per_sample_coverage_counting \n##\n")
f.close()

command = "bedtools coverage -sorted " + \
          " -a " + snakemake.input.region_bed + \
          " -b " + str(snakemake.input.bam) + \
          " -g " + snakemake.input.ref_dict + \
          " > " + snakemake.output.cov_tab

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)
