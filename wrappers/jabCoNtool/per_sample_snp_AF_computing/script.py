######################################
# wrapper for rule: per_sample_snp_AF_computing
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: per_sample_snp_AF_computing \n##\n")
f.close()

command = "alleleCounter " + \
            " -r " + snakemake.input.ref + \
            " -l " + snakemake.input.snp_tsv + \
            " -b  " + str(snakemake.input.bam) + \
            " -o " + snakemake.output.snp_tab + \
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)