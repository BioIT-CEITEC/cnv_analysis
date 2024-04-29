######################################
# wrapper for rule: germline_strelka
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: somatic manta \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("configManta.py --version 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: manta "+version+"\n")
f.close()

#not implemented --targeted
if snakemake.params.library_scope == "wgs":
    scope = ""
else:
    scope = " --exome --callRegions " + snakemake.input.regions_gz + " "


shell("rm -fR " + snakemake.params.dir)

if snakemake.params.calling_type == "tumor_normal":
    command = "configManta.py" + \
              scope + \
              " --normalBam " + snakemake.input.normal + \
              " --tumorBam " + snakemake.input.tumor + \
              " --referenceFasta " + snakemake.input.ref + \
              " --runDir " + snakemake.params.dir + \
              " >> " + log_filename + " 2>&1"
else:
    command = "configManta.py" + \
              scope + \
              " --tumorBam " + snakemake.input.tumor + \
              " --referenceFasta " + snakemake.input.ref + \
              " --runDir " + snakemake.params.dir + \
              " >> " + log_filename + " 2>&1"



f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = snakemake.params.dir + "/runWorkflow.py -m local -j " + str(snakemake.threads) + " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# remove workspace - lots of files messing copying speeds
shell("rm -fR " + snakemake.params.dir + "/workspace")


command = "gunzip -c " + snakemake.params.manta_sv_vcf + " > " + snakemake.output.vcf

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)



# command = "gatk MergeVcfs " + \
#         " -I "+ snakemake.params.dir + "/results/variants/somatic.indels.vcf.gz "+\
#         " -I "+ snakemake.params.dir + "/results/variants/somatic.snvs.vcf.gz "+\
#         " -R "+ snakemake.input.ref +\
#         " -D " + snakemake.input.dict +\
#         " -O " + snakemake.output.vcf +\
#         " >> " + log_filename + " 2>&1"
#
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)


# command =  "cp " + snakemake.params.dir + "/results/variants/somatic.indels.vcf.gz " + str(snakemake.output.vcf) + ".gz"
#
# f = open(log_filename, 'a+')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
#
# command =  "gunzip " + str(snakemake.output.vcf) + ".gz"
#
# f = open(log_filename, 'a+')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
