
rule get_binned_bed_from_dict:
    input:  ref_dict = expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output: bed = expand("CNV_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0],
    log:    "logs/get_binned_bed_from_dict.log"
    threads: 1
    resources: mem=1
    params: window_size=config["wgs_bin_size"]
    conda:  "../wrappers/common_prep/get_binned_bed_from_dict/env.yaml"
    script: "../wrappers/common_prep/get_binned_bed_from_dict/script.py"

rule get_binned_gc_content:
    input:  region_bed = expand("CNV_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0],
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output: gc_content_file = expand("CNV_varcalls/all_samples/GC_profile_{window_size}.cnp",window_size=config["wgs_bin_size"])[0],
    log:    "logs/get_binned_gc_content.log"
    threads: 1
    resources: mem=20
    params: window_size=config["wgs_bin_size"]
    conda:  "../wrappers/common_prep/get_binned_gc_content/env.yaml"
    script: "../wrappers/common_prep/get_binned_gc_content/script.py"

rule unzip_cohort_data:
    input:  cohort_gz = "cohort_data/cohort_cnv_info.tar.gz"
    output: jabContooldata = "cohort_data/cohort_data/jabCoNtool/cohort_info_tab.tsv",
    shell:  "cd cohort_data/; tar -xzf ../{input.cohort_gz}; cd ../"
