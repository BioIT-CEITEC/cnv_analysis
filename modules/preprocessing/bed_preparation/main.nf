process PREPROCESSING {

    input:
    path reference_fasta
    path reference_index

    output:
    path("preprocessed/binned_genome_${params.wgs_bin_size}.bed"), emit: binned_genome
    path("preprocessed/GC_profile_${params.wgs_bin_size}.cnp"), emit: gc_profile

    script:
    """
    mkdir -p preprocessed

    Rscript get_binned_bed_from_dict.R \
        ${reference_index} \
        preprocessed/binned_genome_${params.wgs_bin_size}.bed \
        ${params.wgs_bin_size}

    bedtools nuc -fi ${reference_fasta} \
        -bed preprocessed/binned_genome_${params.wgs_bin_size}.bed \
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