process CNVKIT_CALL {

    tag "${meta.id}"
    
    input:
    tuple val(meta), path()
}