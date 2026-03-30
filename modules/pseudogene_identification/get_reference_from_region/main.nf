process GETREFERENCE_FROM_REGION {

    tag "${region}"
    publishDir "pseudogene/test", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple path(reference), val(region)
    path gene_bed
    path pseudogene_bed


    output:
    tuple val(region), path("${reference.baseName}_${region}.fasta"), emit: pseudogene_reference

    script:
    def chrom = region.split(':')[0]
    def start = region.split(':')[1].split('-')[0]
    def end = region.split(':')[1].split('-')[1]

    def strand_1 = GetStrand(gene_bed, chrom, start, end)
    def strand_2 = GetStrand(pseudogene_bed, chrom, start, end)
    def strand = strand_1 != '-1' ? strand_1 : (strand_2 != '-1' ? strand_2 : '-')
    """
    python ${projectDir}/bin/get_reference.py --region ${region} --reference ${reference} --output ${reference.baseName}_${region}.fasta --strand ${strand}
    """
}

/**
 * Retrieves the strandedness for a given genomic region from a BED file.
 *
 * @param bedFile Path to the BED file (gene_bed or pseudogene_bed).
 * @param chrom   Chromosome name (e.g., chr7 or 7).
 * @param start   Start coordinate (0-based).
 * @param end     End coordinate (1-based).
 * @return        Strandedness as a string ('+' or '-'). Defaults to '-' if not specified.
 */
def GetStrand(bedFile, chrom, start, end) {
    def strand = '-1'
    chrom = chrom.toString()
    start = start.toInteger()
    end = end.toInteger()
    def bed = new File(bedFile.toString())
    if (!bed.exists()) {
        return strand
    }
    bed.eachLine { line ->
        if (line.startsWith('#') || line.trim().isEmpty()) {
            return
        }
        def fields = line.split('\t')
        if (fields.size() < 3) {
            // BED requires at least 4 fields: chrom, start, end, name
            return -1
        }
        def bedChrom = fields[0].toString()
        def bedStart = fields[1].toInteger()
        def bedEnd = fields[2].toInteger()

        if (bedChrom == chrom && start == bedStart && end == bedEnd) {
            if (fields.size() >= 6) {
                def bedStrand = fields[5].trim().toString()
                if (bedStrand == '+' || bedStrand == '-') {
                    strand = bedStrand
                }
            }
        }
    }
    return strand
}