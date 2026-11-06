//
// Long-read QC: NanoPlot and seqkit stats on the raw reads, optional Chopper filtering
//
include { NANOPLOT     } from '../../../modules/nf-core/nanoplot/main'
include { SEQKIT_STATS } from '../../../modules/nf-core/seqkit/stats/main'
include { CHOPPER      } from '../../../modules/nf-core/chopper/main'

workflow READ_QC {

    take:
    ch_reads     // channel: [ val(meta), path(reads) ]
    filter_reads // boolean

    main:
    def ch_multiqc_files = channel.empty()

    NANOPLOT ( ch_reads )
    SEQKIT_STATS ( ch_reads )
    ch_multiqc_files = ch_multiqc_files
        .mix( NANOPLOT.out.txt.map { _meta, txt -> txt } )
        .mix( SEQKIT_STATS.out.stats.map { _meta, tsv -> tsv } )

    def ch_out = ch_reads
    if (filter_reads) {
        CHOPPER ( ch_reads, [] )
        ch_out = CHOPPER.out.fastq
    }

    emit:
    reads         = ch_out           // channel: [ val(meta), path(reads) ]
    multiqc_files = ch_multiqc_files // channel: path
}
