workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    Channel
        .fromPath(file('spreadsheet.csv'))
        .splitCsv(header: true, sep: ',', strip: true)
        .toList()
        .flatMap { process_sample_names(it) }
        .set { ch_input_data }
    emit:
    samples
}

def process_sample_names(List input_data) {
    input_data
        .groupBy { row -> row.patient_id }
        .collect { patient_id, entries -> 
            def meta = [:]
            meta.id = patient_id
            def normal = entries.find { it.donor == 'normal' }?.sample_name
            def tumor = entries.find { it.donor == 'tumor' }?.sample_name
            def samples = []

            if (params.normal_tumor) {
                if (normal && tumor) {
                    samples = [[id:meta.id], [file("mapped/" + normal), file("mapped/" + tumor)]]
                } else {
                    log.warn "Missing normal-tumor pair for patient_id: $patient_id"
                    return [ [ id:meta.id ], null ] 
                }
            } else {
                samples = [ meta, [ file("mapped/" + tumor) ] ]
            }

            return samples 
        }
}