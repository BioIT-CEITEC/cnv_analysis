def merge_wgs_CNV_calls_inputs(wildcards):
    input_dict = {}
    if config["use_gatk_cnv"]:
        input_dict["gatk_cnv_variants"] = expand("CNV_varcalls/{sample_name}/gatk_cnv/CNV_varcalls.vcf",sample_name = sample_tab.sample_name)
    if config["use_jabCoNtool"]:
        input_dict["jabCoNtool_variants"] = "CNV_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    if config["use_control_freec"]:
        input_dict["control_freec_variants"] = expand("CNV_varcalls/{sample_name}/control_freec/CNV_varcalls.tsv",sample_name = sample_tab.sample_name)
    input_dict["region_bed"] = expand("CNV_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0]
    return input_dict

def wgs_CNV_formating_and_visualisation_inputs(wildcards):
    input_dict = {}
    input_dict["all_vars_tsv"] = "CNV_varcalls/all_samples/all_merged_CNV.tsv",
    if config["use_jabCoNtool"]:
        input_dict["jabCoNtool_variants"] = "CNV_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    input_dict["region_bed"] = expand("CNV_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0]
    return input_dict

def merge_trg_CNV_calls_inputs(wildcards):
    input_dict = {}
    if config["use_cnvkit"]:
        input_dict["cnvkit_variants"] = expand("CNV_varcalls/{sample_name}/cnvkit/CNV_calls.cns",sample_name = sample_tab.sample_name)
    if config["use_jabCoNtool"]:
        input_dict["jabCoNtool_variants"] = "CNV_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv"
    input_dict["region_bed"] = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]

    return input_dict

if config["lib_ROI"] == "wgs":
    #MERGING VARIANT
    rule merge_wgs_CNV_calls:
        input:
            unpack(merge_wgs_CNV_calls_inputs)
        output:
            all_vars_tsv= "CNV_varcalls/all_samples/all_merged_CNV.tsv",
        params:
            lib_ROI = config["lib_ROI"],
            overlap = 0.6
        log:
            "logs/process_and_format_CNV.log",
        threads: 8
        conda:  "../wrappers/process_and_format_CNV/env.yaml"
        script: "../wrappers/process_and_format_CNV/script.py"

    #MERGING VARIANT
    rule wgs_CNV_formating_and_visualisation:
        input:
            unpack(wgs_CNV_formating_and_visualisation_inputs)
        output:
            final_res = "CNV_varcalls/final_CNV_visualization.html",
        params:
            lib_ROI = config["lib_ROI"],
            overlap = 0.6
        log:
            "logs/process_and_format_CNV.log",
        threads: 8
        conda:  "../wrappers/process_and_format_CNV/env.yaml"
        script: "../wrappers/process_and_format_CNV/script.py"

else:
    #MERGING VARIANT
    rule merge_trg_CNV_calls:
        input:
            unpack(merge_trg_CNV_calls_inputs)
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

    #MERGING VARIANT
    rule trg_CNV_formating_and_visualisation:
        input:
            unpack(merge_trg_CNV_calls_inputs)
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



# rule final_alignment_report:
#     input:  all_vars_tsv= "final_CNV_results/CNV_variants.tsv",
#     output: html = "CNV_varcalls/final_CNV_report.html"
#     # params: sample_name = sample_tab.sample_name,
#     #         config = "./config.json"
#     # conda: "../wrappers/final_alignment_report/env.yaml"
#     # script: "../wrappers/final_alignment_report/script.Rmd"
#     shell:
#         "mkdir -p reports; touch {output.html}"

#here add all info from previous var calls
def create_cohort_data_inputs(wildcards):
    input_dict = {}
    if config["use_cnvkit"]:
        if config["calling_type"] == "tumor_normal":
            input_dict["cnvkit_normal_coverage_inputs"] = set(expand("CNV_varcalls/{sample_name}/cnvkit/normal.{tag}targetcoverage.cnn",sample_name=sample_tab.loc[sample_tab.tumor_normal == "normal", "donor"].tolist(),tag=["", "anti"]))
        else:
            if len(sample_tab.index) > 4:
                input_dict["cnvkit_normal_coverage_inputs"] = set(expand("CNV_varcalls/{sample_name}/cnvkit/tumor.{tag}targetcoverage.cnn",sample_name=sample_tab.sample_name.tolist(),tag=["", "anti"]))

    if config["use_jabCoNtool"]:
        input_dict["jabCoNtool_cohort_info"] = "CNV_varcalls/all_samples/jabCoNtool/cohort_info_tab.tsv"

    if config["use_cohort_data"]:
        input_dict["previous_cohort_data"] = "cohort_data/cohort_cnv_info.tar.gz"

    return input_dict

rule create_cohort_data:
    input:  unpack(create_cohort_data_inputs)
    output: update_finished_checkfile = "cohort_data/cohort_data_updated"
    log: "logs/create_cohort_data.log"
    conda:  "../wrappers/create_cohort_data/env.yaml"
    script: "../wrappers/create_cohort_data/script.py"
