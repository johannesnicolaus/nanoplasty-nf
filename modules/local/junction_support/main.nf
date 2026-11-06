process JUNCTION_SUPPORT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87660348a4d032cc08d6ac8101ff9d4d27164300d39c9472043bb5001f9ba801/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.25--704fb80d75745bd4' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fai)
    val flank

    output:
    tuple val(meta), path("${prefix}.junction.tsv"), emit: tsv
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // A read supports circularity when its primary alignment on the wrapped
    // reference covers <flank> bp on both sides of the original contig end.
    """
    printf 'contig\\tlength\\tspanning_reads\\n' > ${prefix}.junction.tsv
    samtools view -F 0x904 ${bam} \\
    | awk -v F=${flank} -v FAI=${fai} '
        BEGIN { while ((getline line < FAI) > 0) { split(line, a, "\\t"); len[a[1]] = a[2]; n[a[1]] = 0 } }
        {
            span = 0; cig = \$6
            while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
                op = substr(cig, RLENGTH, 1); num = substr(cig, 1, RLENGTH - 1) + 0
                if (op ~ /[MDN=X]/) span += num
                cig = substr(cig, RLENGTH + 1)
            }
            L = len[\$3]; f = (F * 2 < L) ? F : int(L / 4)
            if (\$4 <= L - f && \$4 + span - 1 >= L + f) n[\$3]++
        }
        END { for (c in len) printf "%s\\t%d\\t%d\\n", c, len[c], n[c] }
    ' >> ${prefix}.junction.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'contig\\tlength\\tspanning_reads\\n' > ${prefix}.junction.tsv
    """
}
