import nextflow.Nextflow
import nextflow.splitter.SplitterEx

class Utils {

public static parseInputVC(inputSh, projectDir, log) {

    def samplesList = inputSh.collect { entry ->
        def meta = [
            sample_name : entry.sample_name
        ]

        def bamPath = entry.bam_path ?: "${projectDir}/mapped/${meta.sample_name}.bam"
        def baiPath = entry.bai_path ?: "${projectDir}/mapped/${meta.sample_name}.bam.bai"

        return [meta, bamPath, baiPath]
    }

    return [samples: samplesList]
}

    public static List<Map> loadSample(Map conf) {
        List<Map> samplesList = []
        conf.samples.each { key, value ->
            def row = [
                sample_name : value.sample_name,
                index       : key,
                bam_path    : value.bam_path,
                bai_path    : value.bai_path
            ]
            samplesList << row
        }
        conf.new_samples = samplesList
        return samplesList
    }

    public static set_data_tags(conf) {
        conf.read_pair_tags = conf.is_paired ? ['_R1','_R2'] : ['']
        conf.read_pair_qc_tags = conf.is_paired ? ['R1','R2'] : ['SE']
        conf.paired_tags = conf.is_paired ? 'PE' : 'SE'
        conf.read_pair_dmtex_tags = conf.is_paired ? ['_R1','_R2'] : ['_R1']
        return conf
    }

    static def load_lib_ROI(Map conf) {
        conf.panel = conf.lib_ROI != "wgs" ? conf.lib_ROI.split('_')[0..-2].join('_') : "wgs"
        return conf
    }

    public static List<Map> readCsvFile(String filePath) {
        def csvData = []
        def lines = new File(filePath).readLines()
        def headers = lines[0].split('\t')
        lines.drop(1).each { line ->
            def values = line.split('\t')
            def row = [:]
            headers.eachWithIndex { header, idx ->
                row[header] = values[idx]
            }
            csvData << row
        }
        return csvData
    }

    static def mixAndCollectVarcalls(ch_current, ch_new) {
        return ch_current
            .mix(ch_new)
            .groupTuple()
            .map { tuple ->
                def meta  = tuple[0]
                def files = tuple[1..-1].flatten()
                [meta, files]
            }
    }

}
