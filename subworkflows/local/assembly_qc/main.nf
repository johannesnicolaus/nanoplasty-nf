//
// Assembly QC: map the reads back to the (wrapped) circular assembly, then
// coverage, mapping stats, reads spanning the start/end, and QUAST stats.
//
include { WRAP_ASSEMBLY    } from '../../../modules/local/wrap_assembly/main'
include { JUNCTION_SUPPORT } from '../../../modules/local/junction_support/main'
include { MINIMAP2_ALIGN   } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_STATS   } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH         } from '../../../modules/nf-core/mosdepth/main'
include { QUAST            } from '../../../modules/nf-core/quast/main'

workflow ASSEMBLY_QC {

    take:
    ch_assembly // channel: [ val(meta), path(fasta) ]
    ch_reads    // channel: [ val(meta), path(reads) ]
    wrap_length // integer
    flank       // integer

    main:
    WRAP_ASSEMBLY ( ch_assembly, wrap_length )

    def ch_align = ch_reads
        .join( WRAP_ASSEMBLY.out.fasta )
        .multiMap { meta, reads, ref ->
            reads: [ meta, reads ]
            ref:   [ meta, ref ]
        }
    MINIMAP2_ALIGN ( ch_align.reads, ch_align.ref, true, 'bai', false, false )

    def ch_bam = MINIMAP2_ALIGN.out.bam.join( MINIMAP2_ALIGN.out.index )

    SAMTOOLS_STATS ( ch_bam, [ [:], [], [] ] )
    MOSDEPTH ( ch_bam.map { meta, bam, bai -> [ meta, bam, bai, [] ] }, [ [:], [] ], false )
    JUNCTION_SUPPORT ( ch_bam.join( WRAP_ASSEMBLY.out.fai ), flank )
    QUAST ( ch_assembly, [ [:], [] ], [ [:], [] ] )

    def ch_multiqc_files = channel.empty()
        .mix( SAMTOOLS_STATS.out.stats.map { _meta, f -> f } )
        .mix( MOSDEPTH.out.global_txt.map { _meta, f -> f } )
        .mix( MOSDEPTH.out.summary_txt.map { _meta, f -> f } )
        .mix( QUAST.out.tsv.map { _meta, f -> f } )

    emit:
    bam           = ch_bam                      // channel: [ val(meta), path(bam), path(bai) ]
    per_base      = MOSDEPTH.out.per_base_bed   // channel: [ val(meta), path(bed.gz) ]
    junction      = JUNCTION_SUPPORT.out.tsv    // channel: [ val(meta), path(tsv) ]
    multiqc_files = ch_multiqc_files            // channel: path
}
