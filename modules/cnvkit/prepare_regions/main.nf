process PREPARE_REGIONS_CNVKIT {

    publishDir "structural_varcalls/all_samples/cnvkit", mode: 'copy'

    //TODO: add proper conda environment in env.yaml
    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null) 

    input:
    path reference_fasta // channel to the reference fasta file
    path lib_ROI // channel to the regions of interest bed file

    output:
    path("results/target.bed"), emit: target_regions
    path("results/antitarget.bed"), emit: antitarget_regions

    script:

        """
        mkdir -p results
        cnvkit.py access ${reference_fasta} -o reference_bed.bed
        cnvkit.py autobin mapped/*.bam -t ${lib_ROI} -g reference_bed.bed
        mv ${lib_ROI}.target.bed results/target.bed
        mv ${lib_ROI}.antitarget.bed results/antitarget.bed
        """

    stub:
        """
        mkdir -p results
        touch target.bed
        touch antitarget.bed
        """
}