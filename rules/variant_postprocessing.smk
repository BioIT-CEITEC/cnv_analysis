def process_and_format_CNV_inputs(wildcards):
    if config["calling_type"] == "tumor_normal":
        sample_name_list=sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist()
    else:
        sample_name_list=sample_tab.sample_name

    input_dict = {}
    if config["use_gatk_cnv"]:
        input_dict["gatk_cnv_variants"] = expand("structural_varcalls/{sample_name}/gatk_cnv/CNV_varcalls.vcf",sample_name = sample_name_list)
    if config["use_cnvkit"]:
        input_dict["cnvkit_variants"] = expand("structural_varcalls/{sample_name}/cnvkit/CNV_calls.cns",sample_name = sample_name_list)
    if config["use_jabCoNtool"]:
        input_dict["jabCoNtool_variants"] = "structural_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    if config["use_control_freec"]:
        input_dict["control_freec_variants"] = expand("structural_varcalls/{sample_name}/control_freec/CNV_varcalls.tsv",sample_name = sample_name_list)
    if config["lib_ROI"] == "wgs":
        input_dict["region_bed"] = expand("structural_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0]
    else:
        input_dict["region_bed"] = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]

    return input_dict

rule process_and_format_CNV:
    input:
        unpack(process_and_format_CNV_inputs)
    output:
        all_vars_tsv= "final_CNV_results/CNV_variants.tsv",
    params:
        lib_ROI = config["lib_ROI"],
        overlap = 0.6
    log:
        "logs/process_and_format_CNV.log",
    threads: 8
    conda:  "../wrappers/process_and_format_CNV/env.yaml"
    script: "../wrappers/process_and_format_CNV/script.py"



rule final_alignment_report:
    input:  all_vars_tsv= "final_CNV_results/CNV_variants.tsv",
    output: html = "reports/final_SV_report.html"
    # params: sample_name = sample_tab.sample_name,
    #         config = "./config.json"
    # conda: "../wrappers/final_alignment_report/env.yaml"
    # script: "../wrappers/final_alignment_report/script.Rmd"
    shell:
        "mkdir -p reports; touch {output.html}"

def create_cohort_data_inputs(wildcards):
    input_dict = {}
    if config["use_cnvkit"]:
        if config["calling_type"] == "tumor_normal":
            input_dict["cnvkit_normal_coverage_inputs"] = set(expand("structural_varcalls/{sample_name}/cnvkit/normal.{tag}targetcoverage.cnn",sample_name=sample_tab.loc[sample_tab.tumor_normal == "normal", "donor"].tolist(),tag=["", "anti"]))
        else:
            if len(sample_tab.index) > 4:
                input_dict["cnvkit_normal_coverage_inputs"] = set(expand("structural_varcalls/{sample_name}/cnvkit/tumor.{tag}targetcoverage.cnn",sample_name=sample_tab.sample_name.tolist(),tag=["", "anti"]))

    if config["use_jabCoNtool"]:
        input_dict["jabCoNtool_cohort_info"] = "structural_varcalls/all_samples/jabCoNtool/cohort_info_tab.tsv"

    if config["use_cohort_data"]:
        input_dict["previous_cohort_data"] = "cohort_data/cohort_cnv_info.tar.gz"

    return input_dict

rule create_cohort_data:
    input:  unpack(create_cohort_data_inputs)
    output: update_finished_checkfile = "cohort_data/cohort_data_updated"
    log: "logs/create_cohort_data.log"
    conda:  "../wrappers/create_cohort_data/env.yaml"
    script: "../wrappers/create_cohort_data/script.py"

# vcf_cnvkit = lambda wildcards: expand("cnv_sv/cnvkit_vcf/{sample_name}.vcf",sample_name=sample_tab.sample_name),
# segment_regions = lambda wildcards: expand("cnv_sv/cnvkit_batch/{sample_name}.cnr",sample_name=sample_tab.sample_name),
# vcf_gatk = lambda wildcards: expand("cnv_sv/gatk_cnv_vcf/{sample_name}.vcf",sample_name=sample_tab.sample_name),
# manta_som_sv_vcf = lambda wildcards: expand("cnv_sv/manta/{sample_name}/results/variants/somaticSV.vcf.gz",donor=sample_tab.sample_name),
# pindel_vcf = lambda wildcards: expand("cnv_sv/pindel/{sample_name}.vcf",donor=sample_tab.sample_name),
# svdb = lambda wildcards: expand("cnv_sv/svdb_query/{sample_name}.svdb_query.vcf",sample_name=sample_tab.sample_name)

# def final_report_inputs(wildcards):
#     input = {}
#     if config["organism"] == "homo_sapiens":
#         tag = "svdb_query"
#     else:
#         tag = "merged"
#
#     if config["calling_type"] == "tumor_normal":
#         if len(used_SV_callers) > 0:
#             input['svdb'] = expand("final_SV_calls/{sample_name}.{tag}.vcf",tag = tag,sample_name=sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist())
#     else:
#         if len(used_SV_callers) > 0:
#             input['svdb'] = expand("final_SV_calls/{sample_name}.{tag}.vcf",tag = tag,sample_name=sample_tab.sample_name)
#     return input
#
# rule final_report:
#     input:  unpack(final_report_inputs)
#     output: "report/final_report.html"
#     shell:
#         "touch {output}"



# def process_and_format_annot_variants_inputs(wildcards):
#     varcall_type = []
#     if len(used_CNV_callers) > 0:
#         varcall_type.append("CNV")
#     if len(used_SV_callers) > 0:
#         varcall_type.append("SV")
#
#     if config["calling_type"] == "tumor_normal":
#         sample_name_list=sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist()
#     else:
#         sample_name_list=sample_tab.sample_name
#
#     return expand("structural_varcalls/{sample_name}.{varcall_type}.final_variants.tsv",varcall_type = varcall_type,sample_name=sample_name_list)
#
#
# rule process_and_format_annot_variants:
#     input:  var_tabs = process_and_format_annot_variants_inputs,
#     output: all_vars_tsv = "final_variant_table.tsv",
#     log:    "logs/postprocess_and_format_annot_variants.log"
#     threads: 1
#     resources:
#         mem_mb=8000
#     # params: reference = config["reference"],
#     #         min_variant_frequency = str(config["min_variant_frequency"]),
#     #         format = config["format"],
#     #         anno_gtf = expand("{ref_dir}/annot/{ref_name}.gtf",ref_dir = reference_directory,ref_name = config["reference"]),
#     #         create_cohort_data = config["create_cohort_data"],
#     #         batch_name = config["entity_name"],
#     #         ref_dir= reference_directory,
#     #         organism=config["organism"],
#     #         mut_load_output_filename= "mutation_loads.xlsx",
#     # conda:  "../wrappers/process_and_format_annot_variants/env.yaml"
#     # script: "../wrappers/process_and_format_annot_variants/script.py"
#     shell:
#         "touch {output}"
