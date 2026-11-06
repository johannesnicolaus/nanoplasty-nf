process NANOPLASTY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Image built by the NANOPlasty repository's GitHub Action (ghcr.io/johannesnicolaus/nanoplasty).
    // With --nanoplasty_script the container is dropped and the script runs with conda
    // (-profile conda) or with the host's tools.
    container "${ params.nanoplasty_script ? '' : (params.nanoplasty_container ?: '') }"

    input:
    tuple val(meta), path(reads), path(seed), path(chloro_reference)
    path nanoplasty_script   // optional: [] to use NANOPlasty0.2.pl from PATH (the container)
    path plant_mito_db       // optional: directory holding a 'plant_mito' BLAST database

    output:
    tuple val(meta), path("${prefix}.fasta")                    , emit: assembly
    tuple val(meta), path("out/assembly_*option_*.fasta")       , emit: options , optional: true
    tuple val(meta), path("out/log_*.txt")                      , emit: log
    tuple val(meta), path("out/quality_scores_*.txt")           , emit: quality , optional: true
    tuple val(meta), path("out/structural_variation_*.vcf")     , emit: sv      , optional: true
    tuple val(meta), path("out/Assembled_reads_*.fasta")        , emit: reads   , optional: true
    tuple val(meta), path("config_nanoplasty.txt")              , emit: config
    tuple val("${task.process}"), val('NANOPlasty'), val(VERSION), emit: versions_nanoplasty, topic: versions
    tuple val("${task.process}"), val('blastn'), eval("blastn -version 2>&1 | sed '1!d; s/^.*blastn: //; s/ .*\$//'"), emit: versions_blastn, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    VERSION = '0.2'
    prefix  = task.ext.prefix ?: "${meta.id}"
    def cmd        = nanoplasty_script ? "perl ${nanoplasty_script}" : 'NANOPlasty0.2.pl'
    def is_ont     = meta.platform == 'ont'
    def np_reads   = is_ont ? reads : ''
    def pb_reads   = is_ont ? '' : reads
    // nf-schema turns empty samplesheet cells into [], so test for a real value
    def min_len    = meta.min_read_length instanceof Integer ? meta.min_read_length : ''
    def np_min     = is_ont ? min_len : ''
    def pb_min     = is_ont ? '' : min_len
    def pm_db      = plant_mito_db ? "\$PWD/${plant_mito_db}/plant_mito" : ''
    def chloro_ref = chloro_reference ? "\$PWD/${chloro_reference}" : ''
    def seed_file  = seed ? "\$PWD/${seed}" : ''
    // NANOPlasty parses the config line by line with "<Key><spaces>= <value>";
    // keys and section headers must match its parser exactly.
    """
    cat > config_nanoplasty.txt <<EOF
    Project:
    -----------------------
    Project name          = ${prefix}
    Type                  = ${meta.type}
    Genome Range          = ${meta.genome_range ?: ''}
    Save assembled reads  = ${params.save_assembled_reads ? 'yes' : ''}
    Seed Input            = ${seed_file}
    Extend seed directly  = 
    Threads               = ${task.cpus}
    Chloroplast sequence  = ${chloro_ref}
    Plant mito DB         = ${pm_db}
    Output path           = \$PWD/out/

    Nanopore reads:
    -----------------------
    Nanopore reads        = ${np_reads ? "\$PWD/${np_reads}" : ''}
    Local DB and NP reads = 
    R10                   = ${is_ont && meta.r10 ? 'yes' : ''}
    Min read length NP    = ${np_min}

    PacBio reads:
    -----------------------
    PacBio reads          = ${pb_reads ? "\$PWD/${pb_reads}" : ''}
    Local DB and PB reads = 
    Min read length PB    = ${pb_min}
    EOF
    sed -i 's/^    //' config_nanoplasty.txt

    ${cmd} -c config_nanoplasty.txt

    # NANOPlasty appends the type to the project name (e.g. <sample>_chloro), so match by glob.
    # When it resolves a complete circle in more than one configuration (e.g. the two
    # orientations of a plastid's single-copy region), the 'option' files hold those circles
    # while 'final' holds the separate contigs, so prefer option_1, then final, then raw.
    option=\$(ls out/assembly_*_option_1.fasta 2>/dev/null | head -n 1)
    final=\$(ls out/assembly_final_*.fasta 2>/dev/null | head -n 1)
    raw=\$(ls out/assembly_raw_*.fasta 2>/dev/null | head -n 1)
    if [ -s "\$option" ]; then
        cp "\$option" ${prefix}.fasta
    elif [ -s "\$final" ]; then
        cp "\$final" ${prefix}.fasta
    elif [ -s "\$raw" ]; then
        cp "\$raw" ${prefix}.fasta
    else
        echo "NANOPlasty produced no assembly for ${prefix}; tail of its log:" >&2
        tail -n 40 out/log_*.txt >&2 || true
        exit 1
    fi

    # Keep the read databases and temporary files out of the published results
    rm -rf out/TMP out/BLAST_DB_* out/tmp_sequences_*
    """

    stub:
    VERSION = '0.2'
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p out
    touch config_nanoplasty.txt out/log_${prefix}.txt
    printf '>${prefix}\\nACGTACGTACGTACGTACGT\\n' > ${prefix}.fasta
    """
}
