process WRAP_ASSEMBLY {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87660348a4d032cc08d6ac8101ff9d4d27164300d39c9472043bb5001f9ba801/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.25--704fb80d75745bd4' }"

    input:
    tuple val(meta), path(fasta)
    val wrap_length

    output:
    tuple val(meta), path("${prefix}.wrapped.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.fasta.fai")    , emit: fai
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Append the first <wrap_length> bases of each contig to its end so that
    // reads crossing the arbitrary start/end of a circular genome map in one piece.
    """
    cp -L ${fasta} ${prefix}.fasta 2>/dev/null || true
    samtools faidx ${prefix}.fasta
    awk -v W=${wrap_length} '
        /^>/ { if (name) print_seq(); name = \$1; seq = ""; next }
        { seq = seq \$0 }
        END { if (name) print_seq() }
        function print_seq(   w, s) {
            w = (W < length(seq)) ? W : length(seq)
            s = seq substr(seq, 1, w)
            print name
            for (i = 1; i <= length(s); i += 80) print substr(s, i, 80)
        }
    ' ${prefix}.fasta > ${prefix}.wrapped.fasta
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.wrapped.fasta ${prefix}.fasta.fai
    """
}
