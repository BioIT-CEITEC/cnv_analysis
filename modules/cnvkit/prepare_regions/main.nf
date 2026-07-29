process PREPARE_REGIONS_CNVKIT {

    conda "${moduleDir}/../env.yaml"

    input:
    path reference_fasta
    path lib_ROI

    output:
    tuple path("results/target.bed"), path("results/antitarget.bed"), emit: prepared_regions

    script:
        """
        mkdir -p results
        cnvkit.py access ${reference_fasta} -o access.bed
        cnvkit.py target ${lib_ROI} --split -o results/target.bed
        cnvkit.py antitarget ${lib_ROI} -g access.bed --avg-size 20000 -o antitarget_raw.bed

        awk '\$1 ~ /^(chr)?([0-9]{1,2}|X|Y|MT|M)\$/' antitarget_raw.bed > results/antitarget.bed
        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
}