//
// Organelle annotation, chosen by assembly type:
//   mito        -> MITOS2
//   chloro      -> Chloe (default) or the homology annotator
//   plant_mito  -> homology annotator (miniprot + BLASTN rRNAs + tRNAscan-SE)
// Every branch ends in a gene-level GFF3 carrying a gene_class attribute.
//
include { MITOS                                   } from '../../../modules/local/mitos/main'
include { CHLOE                                   } from '../../../modules/local/chloe/main'
include { ORGANELLE_GFF                           } from '../../../modules/local/organelle_gff/main'
include { MINIPROT_ALIGN                          } from '../../../modules/nf-core/miniprot/align/main'
include { TRNASCANSE                              } from '../../../modules/nf-core/trnascanse/main'
include { BLAST_MAKEBLASTDB                       } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                            } from '../../../modules/nf-core/blast/blastn/main'

workflow ANNOTATION {

    take:
    ch_assembly      // channel: [ val(meta), path(fasta) ]
    mitos_refdb      // path: MITOS2 reference (directory or .tar.bz2)
    chloro_annotator // string: 'chloe' or 'homology'
    homology_refs    // map: [ chloro: [proteins, rrna], plant_mito: [proteins, rrna] ]

    main:
    def ch_gff_in = channel.empty()

    def ch_by_tool = ch_assembly.branch { meta, _fasta ->
        mitos:    meta.type == 'mito'
        chloe:    meta.type == 'chloro' && chloro_annotator == 'chloe'
        homology: true
    }

    //
    // Animal mitochondria: MITOS2
    //
    MITOS ( ch_by_tool.mitos, mitos_refdb )
    ch_gff_in = ch_gff_in.mix( MITOS.out.gff.map { meta, gff -> [ meta, 'normalise', 'MITOS2', gff ] } )

    //
    // Angiosperm plastids: Chloe
    //
    CHLOE ( ch_by_tool.chloe )
    ch_gff_in = ch_gff_in.mix( CHLOE.out.gff.map { meta, gff -> [ meta, 'normalise', 'Chloe', gff ] } )

    //
    // Plant mitochondria (and plastids with --chloro_annotator homology):
    // protein, rRNA and tRNA evidence merged into best non-overlapping hits
    //
    def ch_refs = channel.of(
        [ [ id: 'chloro'     ], file(homology_refs.chloro[0],     checkIfExists: true), file(homology_refs.chloro[1],     checkIfExists: true) ],
        [ [ id: 'plant_mito' ], file(homology_refs.plant_mito[0], checkIfExists: true), file(homology_refs.plant_mito[1], checkIfExists: true) ],
    )
    BLAST_MAKEBLASTDB ( ch_refs.map { meta, _prot, rrna -> [ meta, rrna ] }, [] )

    // Pair every homology sample with the reference set for its organelle
    def ch_homology = ch_by_tool.homology
        .map { meta, fasta -> [ meta.type, meta, fasta ] }
        .combine( ch_refs.map { meta, prot, _rrna -> [ meta.id, prot ] }, by: 0 )
        .combine( BLAST_MAKEBLASTDB.out.db.map { meta, db -> [ meta.id, db ] }, by: 0 )
        .multiMap { _type, meta, fasta, prot, db ->
            fasta: [ meta, fasta ]
            prot:  [ meta, prot ]
            db:    [ meta, db ]
        }

    MINIPROT_ALIGN ( ch_homology.prot, ch_homology.fasta )
    BLAST_BLASTN ( ch_homology.fasta, ch_homology.db, [], [], false )
    TRNASCANSE ( ch_homology.fasta, false, false, true, false )

    ch_gff_in = ch_gff_in.mix(
        MINIPROT_ALIGN.out.gff
            .join( BLAST_BLASTN.out.txt )
            .join( TRNASCANSE.out.gff )
            .map { meta, prot_gff, blast, trna_gff -> [ meta, 'homology', 'homology', [ prot_gff, blast, trna_gff ] ] }
    )

    ORGANELLE_GFF ( ch_gff_in )

    emit:
    gff        = ORGANELLE_GFF.out.gff // channel: [ val(meta), path(gff3) ]
    mitos      = MITOS.out.results     // channel: [ val(meta), path(dir) ]
    chloe      = CHLOE.out.gff         // channel: [ val(meta), path(gff) ]
}
