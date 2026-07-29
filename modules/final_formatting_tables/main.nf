process FINAL_FORMATTING_TABLES {

    publishDir "structural_varcalls/", mode: 'copy'
    conda "${moduleDir}/env.yaml"
    cache false

    input:
    path(cnv_varcalls)
    path(dna_panel)
    path(gtf_file)
    path(coverage_files)

    output:
    path("cohort_results/all_samples_merged.tsv"),    emit: all_merged
    path("cohort_results/all_samples_annotated.tsv"), emit: all_annotated
    path("cohort_results/all_samples_combined.tsv"),  emit: all_combined
    path("cohort_results/all_samples_smoothed.tsv"),  emit: all_smoothed
    path("cohort_results/final_report.html"),         emit: html_report
    path("cohort_results/coverage/*.tsv"), optional: true, emit: coverage_out

    script:
    """
    Rscript ${projectDir}/bin/combine_final_tables.R ./ cohort_results ${params.smooth_gap_bp} ${dna_panel}

    mkdir -p cohort_results/coverage
    for f in *.region_coverage.tsv; do
        [ -f "\$f" ] && cp "\$f" cohort_results/coverage/
    done
    """

    stub:
    """
    mkdir -p cohort_results/coverage
    touch cohort_results/all_samples_merged.tsv
    touch cohort_results/all_samples_annotated.tsv
    touch cohort_results/all_samples_combined.tsv
    touch cohort_results/all_samples_smoothed.tsv
    touch cohort_results/final_report.html
    touch cohort_results/coverage/stub.region_coverage.tsv
    """

}