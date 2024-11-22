process PREPROCESSING {

    input:
    path reference
    path reference_index

    output:
    path("preprocessed/"), emit: preprocessed

    script:
    """
    mkdir -p preprocessed

    Rscript get_binned_bed_from_dict.R \
        ${params.assembly}.fa.fai \
        preprocessed/binned_genome_${params.wgs_bin_size}.bed \
        ${params.wgs_bin_size}

    bedtools nuc -fi ${params.assembly}.fa \
        -bed ppreprocessed/binned_genome_${params.wgs_bin_size}.bed \
        > GC_profile_${params.wgs_bin_size}.cnp.tmp

    Rscript get_binned_gc_content.R \
        GC_profile_${params.wgs_bin_size}.cnp.tmp \
        preprocessed/GC_profile_${params.wgs_bin_size}.cnp
    """

    stub:
    """
    mkdir -p preprocessed
    touch preprocessed/binned_genome_${params.wgs_bin_size}.bed
    touch preprocessed/GC_profile_${params.wgs_bin_size}.cnp
    """
}