######################################
# wrapper for rule: germline_strelka
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: gridss \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("gridss --version 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: manta "+version+"\n")
f.close()



# #not implemented --targeted
# if snakemake.params.library_scope == "wgs":
#     scope = ""
# else:
#     scope = " --exome --callRegions " + snakemake.input.regions_gz + " "



command = "ln -fst " + os.path.dirname(snakemake.input.ref) + " " + " ".join(snakemake.input.bwa_index)
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# gridss --reference <reference.fa> --output <output.vcf.gz> [--threads n] [--jar gridss.jar] [--workingdir <directory>] [--jvmheap 30g] [--blacklist <exclude_list.bed>] [--steps All|PreProcess|Assemble|Call] [--configuration gridss.properties] [--maxcoverage 50000] [--labels input1,input2,...] input1.bam [input2.bam [...]]
# [--blacklist <exclude_list.bed>] [--maxcoverage 50000] input1.bam [input2.bam [...]]

command = "gridss --jvmheap 30g" + \
          " --reference " + snakemake.input.ref + \
          " --output " + snakemake.output.vcf + \
          " --threads " + str(snakemake.threads) + \
          " --workingdir " + os.path.dirname(snakemake.output.vcf) + \
          " " + " ".join(snakemake.input.bams) + \
          " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


command = "rm " + " ".join([os.path.dirname(snakemake.input.ref) + "/" + os.path.basename(bwa_index_file) for bwa_index_file in snakemake.input.bwa_index])

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)



