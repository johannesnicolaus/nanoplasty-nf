process ORGANELLE_GFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pycirclize_python:731acb4cc6a99f3c' :
        'community.wave.seqera.io/library/pycirclize_python:ad50b4c46d2a3f7e' }"

    input:
    // mode 'normalise': annotations = [ annotator_gff ]
    // mode 'homology' : annotations = [ miniprot_gff, blastn_tsv, trnascan_gff ]
    tuple val(meta), val(mode), val(source), path(annotations, stageAs: "in?/*")

    output:
    tuple val(meta), path("${prefix}.gff3"), emit: gff
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def files = annotations instanceof List ? annotations : [annotations]
    def cmd = mode == 'homology'
        ? "homology --miniprot ${files[0]} --blast ${files[1]} --trnascan ${files[2]}"
        : "normalise --gff ${files[0]} --source ${source}"
    """
    organelle_gff.py ${cmd} ${args} --out ${prefix}.gff3
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '##gff-version 3' > ${prefix}.gff3
    """
}
