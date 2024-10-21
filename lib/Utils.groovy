import nextflow.Nextflow
import nextflow.splitter.SplitterEx

class Utils {

public static parseInput_VC(input_sh, normal_tumor, log) {

    def input_sh_path = Utils.getFilePath(input_sh)
    def inputs = nextflow.splitter.SplitterEx.splitCsv(input_sh_path, [header: true])
        .groupBy { it['patient_id'] }
        .collect { patient_id, entries ->

            def meta = [id: patient_id]
            def normal = entries.find { it.donor == 'normal' }?.sample_name
            def tumor = entries.find { it.donor == 'tumor' }?.sample_name

            def patientInputs
            if (normal_tumor) {
                if (normal && tumor) {
                    patientInputs = [[meta.id], [("mapped/${tumor}.bam"), ("mapped/${tumor}.bam.bai"), ("mapped/${normal}.bam"), ("mapped/${normal}.bam.bai")]]
                } else {
                    log.warn "Missing normal-tumor pair for patient_id: ${patient_id}"
                    patientInputs = [null, null]
                }
            } else {
                patientInputs = [[meta.id], [("mapped/${tumor}.bam"), ("mapped/${tumor}.bam.bai"), [], []]]
            }
            return patientInputs
        }
}

    static public getFilePath(path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

}