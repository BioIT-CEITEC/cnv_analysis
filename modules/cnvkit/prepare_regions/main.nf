process PREPARE_REGIONS_CNVKIT {

    conda "${moduleDir}/../env.yaml"

    input:
    path reference_fasta // channel to the reference fasta file
    path lib_ROI // channel to the regions of interest bed file
    path bam_files // channel to the the list of bam files for the autobinning process

    output:
    tuple path("results/target.bed"), path("results/antitarget.bed"), emit: prepared_regions

    script:

      def panel = lib_ROI.toString().replace('.bed', '') 

        """
        mkdir -p results
        cnvkit.py access ${reference_fasta} -o reference_bed.bed
        cnvkit.py autobin *.bam -t ${lib_ROI} -g reference_bed.bed
        mv ${panel}.target.bed results/target.bed
        mv ${panel}.antitarget.bed results/antitarget.bed
        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
}