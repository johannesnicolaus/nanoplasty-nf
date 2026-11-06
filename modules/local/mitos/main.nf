process MITOS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::mitos=2.1.10"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path refdb   // MITOS2 reference: a directory such as refseq89m, or its .tar.bz2 archive

    output:
    tuple val(meta), path("${prefix}")           , emit: results
    tuple val(meta), path("${prefix}/result.gff"), emit: gff
    tuple val("${task.process}"), val('mitos'), val('2.1.10'), emit: versions_mitos, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p refdir
    if [ -d ${refdb} ]; then
        ln -s \$(readlink -f ${refdb}) refdir/
    else
        tar -xjf ${refdb} -C refdir
    fi
    REF=\$(ls refdir | head -n 1)

    mkdir -p ${prefix}
    runmitos \\
        -i ${fasta} \\
        -o ${prefix} \\
        -r \$REF \\
        -R \$PWD/refdir \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    echo '##gff-version 3' > ${prefix}/result.gff
    """
}
