process CHLOE {
    tag "$meta.id"
    label 'process_low'

    // Chloe (Small lab, UWA) is not on bioconda. This is willmundlab/chloe:0.1.17 with
    // permissions fixed for non-root users and procps added (see containers/chloe).
    container "ghcr.io/johannesnicolaus/nanoplasty-chloe:0.1.17"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.chloe.gff"), emit: gff
    tuple val(meta), path("${prefix}.chloe.fa") , emit: fasta, optional: true
    tuple val("${task.process}"), val('chloe'), val('0.1.17'), emit: versions_chloe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Writable depot for Julia's compile cache, then the packages bundled in the image
    export JULIA_DEPOT_PATH=\$PWD/.julia:/root/.julia

    mkdir -p input
    cp -L ${fasta} input/${prefix}.fa
    # --no-transform keeps the input orientation, so coordinates match the assembly
    julia --project=/usr/local/bin/chloe /usr/local/bin/chloe/chloe.jl annotate \\
        --no-transform \\
        -o . \\
        ${args} \\
        input/${prefix}.fa
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '##gff-version 3' > ${prefix}.chloe.gff
    """
}
