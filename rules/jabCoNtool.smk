def get_region_bed_input(wildcards):
    if config["lib_ROI"] != "wgs":
        return expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
    else:
        return expand("CNV_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0]

rule jabCoNtool_per_sample_coverage:
    input:  bam = "mapped/{sample_name}.bam",
            region_bed = get_region_bed_input,
            ref_dict= expand("{ref_dir}/seq/{ref_name}.fa.fai",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output: cov_tab = "CNV_varcalls/{sample_name}/jabCoNtool/region_coverage.tsv",
    log:    "logs/{sample_name}/jabCoNtool/get_coverage.log"
    threads: 8
    resources: mem=10
    conda:  "../wrappers/jabCoNtool/per_sample_coverage_computing/env.yaml"
    script: "../wrappers/jabCoNtool/per_sample_coverage_computing/script.py"

rule jabCoNtool_per_sample_snp_AF:
    input:  bam = "mapped/{sample_name}.bam",
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            snp_tsv = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.tsv",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: snp_tab = "CNV_varcalls/{sample_name}/jabCoNtool/{tumor_normal}.snpAF.tsv",
    log:    "logs/{sample_name}/jabCoNtool/{tumor_normal}_get_snpAF.log"
    threads: 8
    resources: mem=10
    conda:  "../wrappers/jabCoNtool/per_sample_snp_AF_computing/env.yaml"
    script: "../wrappers/jabCoNtool/per_sample_snp_AF_computing/script.py"


def jabCoNtool_cnv_computation_inputs(wildcards):
    input_dict = {}
    input_dict["sample_cov"] = set(expand("CNV_varcalls/{sample_name}/jabCoNtool/sample.region_coverage.tsv",sample_name=sample_tab.sample_name.tolist()))
    if config["jabCoNtool_use_snps"] == True:
        input_dict["snp_AF"] = set(expand("CNV_varcalls/{sample_name}/jabCoNtool/sample.snpAF.tsv",sample_name=sample_tab.sample_name.tolist()))
    if config["use_cohort_data"] == True:
        input_dict["cohort_data"] = "cohort_data/cohort_data/jabCoNtool/cohort_info_tab.tsv"
    if config["lib_ROI"] == "wgs":
        input_dict["region_bed"] = expand("CNV_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0]
        if config["jabCoNtool_normalize_to_GC"] == True:
            input_dict["GC_profile_file"] = expand("CNV_varcalls/all_samples/GC_profile_{window_size}.cnp",window_size=config["wgs_bin_size"])[0]
        if config["jabCoNtool_remove_centromeres"] == True:
            input_dict["cytoband_file"] = expand("{ref_dir}/other/cytoband/{ref}.cytoband.tsv",ref_dir=reference_directory,ref = config["reference"])[0]
    else:
        input_dict["region_bed"] = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
    if config["jabCoNtool_use_snps"] == True:
        input_dict["snp_bed"] = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.tsv",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
    return input_dict


rule jabCoNtool_cnv_computation:
    input: unpack(jabCoNtool_cnv_computation_inputs)
    output: all_res_prob_tab="CNV_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv",
            cohort_info_tab="CNV_varcalls/all_samples/jabCoNtool/cohort_info_tab.tsv"
    params: jabCoNtool_predict_TL = config["jabCoNtool_predict_TL"],
            calling_type = config["calling_type"],
            lib_ROI= config["lib_ROI"],
            max_CNV_occurance_in_cohort = config["max_CNV_occurance_in_cohort"]
    log:    "logs/all_samples/jabCoNtool/cnv_computation.log",
    threads: workflow.cores
    conda:  "../wrappers/jabCoNtool/cnv_computation/env.yaml"
    script: "../wrappers/jabCoNtool/cnv_computation/script.py"


# rule jabCoNtool_get_per_sample_res:
#     input:  all_res_prob_tab="CNV_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv"
#     output: CNV_res="CNV_varcalls/{sample_name}/jabCoNtool/CNV_varcalls.tsv",
#     log:    "logs/{sample_name}/jabCoNtool/get_per_sample_res.log",
#     conda:  "../wrappers/jabCoNtool/get_per_sample_res/env.yaml"
#     script: "../wrappers/jabCoNtool/get_per_sample_res/script.py"
