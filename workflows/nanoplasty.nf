/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NANOPLASTY as NANOPLASTY_ASSEMBLE } from '../modules/local/nanoplasty/main'
include { CIRCOS_PLOT            } from '../modules/local/circos_plot/main'
include { READ_QC                } from '../subworkflows/local/read_qc'
include { ASSEMBLY_QC            } from '../subworkflows/local/assembly_qc'
include { ANNOTATION             } from '../subworkflows/local/annotation'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanoplasty_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOPLASTY {

    take:
    ch_samplesheet // channel: [ val(meta), path(reads), path(seed), path(chloro_reference), path(assembly) ]
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: read QC (and optional filtering)
    //
    READ_QC (
        ch_samplesheet.map { meta, reads, _seed, _chloro, _asm -> [ meta, reads ] },
        params.filter_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix( READ_QC.out.multiqc_files )

    //
    // MODULE: NANOPlasty assembly, unless the samplesheet already gives an assembly
    //
    def ch_inputs = ch_samplesheet
        .map { meta, _reads, seed, chloro, asm -> [ meta, seed, chloro, asm ] }
        .join( READ_QC.out.reads )
        .branch { meta, seed, chloro, asm, reads ->
            given:    asm
                return [ meta, asm ]
            assemble: true
                return [ meta, reads, seed, chloro ]
        }

    NANOPLASTY_ASSEMBLE (
        ch_inputs.assemble,
        params.nanoplasty_script ? file(params.nanoplasty_script, checkIfExists: true) : [],
        params.plant_mito_db     ? file(params.plant_mito_db,     checkIfExists: true) : []
    )
    def ch_assembly = NANOPLASTY_ASSEMBLE.out.assembly.mix( ch_inputs.given )

    //
    // SUBWORKFLOW: assembly QC
    //
    ASSEMBLY_QC (
        ch_assembly,
        READ_QC.out.reads,
        params.wrap_length,
        params.junction_flank
    )
    ch_multiqc_files = ch_multiqc_files.mix( ASSEMBLY_QC.out.multiqc_files )

    //
    // SUBWORKFLOW: annotation
    //
    def ch_gff = channel.empty()
    if (!params.skip_annotation) {
        ANNOTATION (
            ch_assembly,
            file(params.mitos_refdb),
            params.chloro_annotator,
            [
                chloro:     [
                    params.chloro_proteins     ?: "${projectDir}/assets/refs/chloro_proteins.faa",
                    params.chloro_rrnas        ?: "${projectDir}/assets/refs/chloro_rrna.fna",
                ],
                plant_mito: [
                    params.plant_mito_proteins ?: "${projectDir}/assets/refs/plant_mito_proteins.faa",
                    params.plant_mito_rrnas    ?: "${projectDir}/assets/refs/plant_mito_rrna.fna",
                ],
            ]
        )
        ch_gff = ANNOTATION.out.gff
    }

    //
    // MODULE: circular plot and per-sample summary
    //
    def ch_plot_in = ch_assembly
        .join( ch_gff, remainder: true )
        .join( ASSEMBLY_QC.out.per_base )
        .join( ASSEMBLY_QC.out.junction )
        .map { meta, fasta, gff, depth, junction -> [ meta, fasta, gff ?: [], depth, junction ] }
    CIRCOS_PLOT ( ch_plot_in )

    def ch_summary = CIRCOS_PLOT.out.summary
        .map { _meta, tsv -> tsv }
        .collectFile(
            name: 'organelle_summary_mqc.tsv',
            keepHeader: true,
            skip: 1,
            sort: true,
            storeDir: "${outdir}/summary"
        )
    ch_multiqc_files = ch_multiqc_files.mix( ch_summary )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'nanoplasty_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'nanoplasty'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
