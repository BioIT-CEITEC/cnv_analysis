class Utils {

    public static set_data_tags(conf) {
        conf.read_pair_tags = conf.is_paired ? ['_R1', '_R2'] : ['']
        conf.read_pair_qc_tags = conf.is_paired ? ['R1', 'R2'] : ['SE']
        conf.paired_tags = conf.is_paired ? 'PE' : 'SE'
        conf.read_pair_dmtex_tags = conf.is_paired ? ['_R1', '_R2'] : ['_R1']
        return conf
    }

    public static List<Map> loadSample(Map conf) {
        List<Map> samplesList = []
        conf.samples.each { key, value ->
            Map<String, Object> row = value
            row['index'] = key
            samplesList << row
        }
        conf.new_samples = samplesList
        return samplesList
    }

    static def load_lib_ROI(Map conf) {
        conf.panel = conf.lib_ROI != "wgs" ? conf.lib_ROI.split('_')[0..-2].join('_') : "wgs"
        return conf
    }

    public static List<Map> readCsvFile(String filePath) {
        def csvData = []
        def lines = new File(filePath).readLines()
        def headers = lines[0].split(',')
        lines.drop(1).each { line ->
            def values = line.split(',')
            def row = [:]
            headers.eachWithIndex { header, idx ->
                row[header] = values[idx]
            }
            csvData << row
        }
        return csvData
    }

    public static load_organism(conf) {

        //def organism_tab = readCsvFile("${conf.globalResources}/reference_info/organism_tab.tsv")
        def organism_tab = readCsvFile("/Users/alejandromedaglia/Documents/nf_core_test/tutorial/data.csv")

        def organism_data = organism_tab.find { it.assembly == conf.assembly }

        conf.new_release = !(conf.containsKey('release')) || conf.release == "UNK_UNK" ? organism_data.release : conf.release.split("_")[-1]

        conf.kegg_code = organism_data.kegg_code
        conf.reference_dir = "${conf.globalResources}/references/${conf.organism}/${conf.assembly}"
        conf.organism_fasta = "${conf.reference_dir}/seq/${conf.assembly}.fa"
        conf.organism_ucsc = "${conf.reference_dir}/seq/${conf.assembly}.fa.fai.ucsc"
        conf.organism_gtf = "${conf.reference_dir}/annot/${conf.new_release}/${conf.assembly}.gtf"
        conf.organism_gtf_cellranger = "${conf.reference_dir}/annot/${conf.new_release}/${conf.assembly}_cellranger.gtf"
        conf.organism_cds_fasta = "${conf.reference_dir}/annot/${conf.new_release}/${conf.assembly}.cds.fa"
        conf.organism_cdna_fasta = "${conf.reference_dir}/annot/${conf.new_release}/${conf.assembly}.cdna.fa"
        conf.organism_star = "${conf.reference_dir}/tool_data/STAR/${conf.new_release}/SAindex"
        conf.organism_star_solo = "${conf.reference_dir}/tool_data/STAR/${conf.new_release}/STAR_cellranger/SAindex"
        conf.organism_rsem = "${conf.reference_dir}/tool_data/RSEM/${conf.new_release}/${conf.assembly}.idx.fa"
        conf.organism_salmon = "${conf.reference_dir}/tool_data/Salmon/${conf.new_release}"
        conf.organism_salmon_gentrome = "${conf.reference_dir}/tool_data/Salmon/${conf.new_release}/Salmon_decoy/gentrome.fa"
        conf.organism_kallisto = "${conf.reference_dir}/tool_data/Kallisto/${conf.new_release}/Kallisto"
        conf.organism_picard_bed12 = "${conf.reference_dir}/annot/${conf.new_release}/Picard/${conf.assembly}.bed12"
        conf.organism_picard_refFlat = "${conf.reference_dir}/annot/${conf.new_release}/Picard/${conf.assembly}.refFlat"
        conf.organism_ncbi_general = "${conf.reference_dir}/seq/BOWTIE2_fastq_screen/${conf.assembly}.ncbi.fna"
        conf.organism_ncbi_gff = "${conf.reference_dir}/seq/BOWTIE2_fastq_screen/${conf.assembly}.ncbi.gff"
        conf.organism_ncbi_rRNA = "${conf.reference_dir}/seq/BOWTIE2_fastq_screen/${conf.assembly}.ncbi.rRNA.fasta"
        conf.organism_ncbi_tRNA = "${conf.reference_dir}/seq/BOWTIE2_fastq_screen/${conf.assembly}.ncbi.tRNA.fasta"
        conf.organism_bwa = "${conf.reference_dir}/tool_data/BWA/${conf.assembly}.bwt"
        conf.organism_bowtie2 = "${conf.reference_dir}/tool_data/Bowtie2/${conf.new_release}/${conf.assembly}.1.bt2"
        conf.organism_vep_dir = "${conf.reference_dir}/annot/${conf.new_release}/vep"
        conf.organism_chr_sizes = "${conf.reference_dir}/seq/${conf.assembly}.chrom.sizes"
        conf.organism_dict = "${conf.reference_dir}/seq/${conf.assembly}.dict"
        conf.organism_snp_bed = "${conf.reference_dir}/seq/${conf.assembly}.snp.bed"
        conf.organism_custom_DB_folder = "${conf.reference_dir}/others/custom_variant_annot_DBs"
        conf.organism_cadd_db_snvs = "${conf.reference_dir}/others/CADD_scores_DB/whole_genome_SNVs.tsv.gz"
        conf.organism_cadd_db_indels = "${conf.reference_dir}/others/CADD_scores_DB/gnomad.genomes.r3.0.indel.tsv.gz"
        conf.organism_dbsnp = "${conf.reference_dir}/others/dbSNP/common_all.vcf.gz"
        conf.organism_cytoband = "${conf.reference_dir}/others/cytoband/${conf.assembly}.cytoband.tsv"
        conf.organism_svdb = "${conf.reference_dir}/others/svdb/gnomad_v2.1_sv.sites.vcf"
        conf.organism_transcriptome = "${conf.reference_dir}/tool_data/Cellranger/refdata-gex-${conf.assembly}"
        conf.organism_gmap_folder = "${conf.reference_dir}/tool_data/GMAP/${conf.new_release}"
        conf.organism_splicesites = "${conf.reference_dir}/tool_data/GMAP/${conf.new_release}/${conf.assembly}.splicesites"
        conf.organism_introns = "${conf.reference_dir}/tool_data/GMAP/${conf.new_release}/${conf.assembly}.introns"
        conf.organism_map_splice = "${conf.reference_dir}/tool_data/GMAP/${conf.new_release}/${conf.assembly}.maps/${conf.assembly}.splicesites.iit"
        conf.organism_map_intron = "${conf.reference_dir}/tool_data/GMAP/${conf.new_release}/${conf.assembly}.maps/${conf.assembly}.introns.iit"

        if (conf.containsKey('lib_ROI')) {
            load_lib_ROI(conf)
            conf.organism_dna_panel = "${conf.reference_dir}/others/DNA_ROI/${conf.panel}/${conf.panel}.bed"
            conf.organism_snps_panel = "${conf.reference_dir}/others/snp/${conf.panel}/${conf.panel}_snps.bed"
            conf.organism_interval_list = "${conf.reference_dir}/others/DNA_ROI/${conf.panel}/${conf.panel}.interval_list"
        }
    }
}