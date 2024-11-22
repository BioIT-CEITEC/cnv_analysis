process PREPARE_REGIONS_CNVKIT {

    publishDir "structural_varcalls/all_samples/cnvkit", mode: 'copy'

    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null)

    input:
    ch_reference_fasta

    output:
    tuple path("*.target.bed"), path("*.antitarget.bed"), emit: cnvkit_beds

    script:

        """
        cnvkit.py access $params.organism_fasta -o reference_bed.bed
        cnvkit.py autobin $bams -t $regions -g reference_bed.bed
        mv ${params.lib_ROI}.target.bed target.bed
        mv ${params.lib_ROI}.antitarget.bed antitarget.bed
        """
}