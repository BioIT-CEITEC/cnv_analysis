process PREPARE_REGIONS_CNVKIT {

    publishDir "structural_varcalls/all_samples/cnvkit", mode: 'copy'

    //TODO: add proper conda environment in env.yaml
    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null) 

    input:
    path all_bams // channel of list of all the bam files
    path reference_fasta // channel to the reference fasta file
    path lib_ROI // channel to the regions of interest bed file

    output:
    tuple path("target.bed"), path("antitarget.bed"), emit: cnvkit_beds

    script:

        """
        cnvkit.py access ${reference_fasta} -o reference_bed.bed
        cnvkit.py autobin ${all_bams} -t ${lib_ROI} -g reference_bed.bed
        mv ${lib_ROI}.target.bed target.bed
        mv ${lib_ROI}.antitarget.bed antitarget.bed
        """

    stub:
        """
        touch target.bed
        touch antitarget.bed
        """
}