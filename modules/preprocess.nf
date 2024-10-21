process PREPROCESSING {

    input:
    path(reference)
    path(reference_index)

    output:
    tuple path("binned_genome_${params.wgs_bin_size}.bed"), path("GC_profile_${params.wgs_bin_size}.cnp"), emit: preprocessed

    script:

    """
    Rscript get_binned_bed_from_dict.R \
        ${params.assembly}.fa.fai \
        binned_genome_${params.wgs_bin_size}.bed \
        ${params.wgs_bin_size}

    bedtools nuc -fi ${params.assembly}.fa \
        -bed binned_genome_${params.wgs_bin_size}.bed \
        > GC_profile_${params.wgs_bin_size}.cnp.tmp

    Rscript get_binned_gc_content.R \
        GC_profile_${params.wgs_bin_size}.cnp.tmp \
        GC_profile_${params.wgs_bin_size}.cnp
    """
}