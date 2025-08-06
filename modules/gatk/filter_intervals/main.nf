process FILTER_INTERVALS_GATK {

    conda "${moduleDir}/../env.yaml"

    input:
    path interval_list
    path annotated_intervals
    path lib_ROI
    path read_counts

    output:
    path("*.gc.filtered.interval_list"), emit: qc_filtered_regions_gatk

    script:

      def panel = lib_ROI.toString().replace('.bed', '') 
      
        """
        gatk FilterIntervals \\
          -L ${interval_list} \\
          --annotated-intervals ${annotated_intervals} \\
          -I ${read_counts.join(" -I ")} \\
          -imr OVERLAPPING_ONLY \\
          -O ${panel}.gc.filtered.interval_list

        
        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
        
}