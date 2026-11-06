process CIRCOS_PLOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pycirclize_python:731acb4cc6a99f3c' :
        'community.wave.seqera.io/library/pycirclize_python:ad50b4c46d2a3f7e' }"

    input:
    tuple val(meta), path(fasta), path(gff), path(depth_bed), path(junction)

    output:
    tuple val(meta), path("${prefix}.circos.png"), emit: png
    tuple val(meta), path("${prefix}.circos.svg"), emit: svg
    tuple val(meta), path("${prefix}.summary.tsv"), emit: summary
    tuple val("${task.process}"), val('pycirclize'), eval("python -c 'import pycirclize; print(pycirclize.__version__)'"), emit: versions_pycirclize, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD/.mpl
    plot_organelle.py \\
        --fasta ${fasta} \\
        ${gff ? "--gff ${gff}" : ''} \\
        ${depth_bed ? "--depth ${depth_bed}" : ''} \\
        ${junction ? "--junction ${junction}" : ''} \\
        --sample ${meta.id} \\
        --type ${meta.type} \\
        ${meta.genome_range ? "--genome-range ${meta.genome_range}" : ''} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.circos.png ${prefix}.circos.svg
    printf 'Sample\\n${meta.id}\\n' > ${prefix}.summary.tsv
    """
}
