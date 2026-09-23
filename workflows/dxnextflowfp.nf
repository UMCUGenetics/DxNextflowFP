/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWAMEM2_MEM } from '../modules/nf-core/bwamem2/mem/main'
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from '../modules/nf-core/samtools/merge/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'


include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { BAM_VARIANTCALLING_INTERVALS } from '../subworkflows/UMCUGenetics/bam_variantcalling_intervals/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dxnextflowfp_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DXNEXTFLOWFP {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    ch_genome_fasta = Channel.fromPath("${params.genome_fasta}").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_genome_fasta_index = Channel.fromPath("${params.genome_fasta}.fai").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_genome_dict = Channel.fromPath("${params.genome_dict}").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_bwa_index = Channel.fromPath("${params.bwa_index}*").map{ file -> [file.getSimpleName(), file] }.groupTuple().collect()
    ch_dbsnp = Channel.fromPath("${params.dbsnp}").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_dbsnp_index = Channel.fromPath("${params.dbsnp}.tbi").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_intervals = Channel.fromPath("${params.intervals}").map{ file -> [file.getSimpleName(), file] }.collect()

    ///
    /// Workflow
    ///
    FASTQC(ch_samplesheet)
    TRIMGALORE(ch_samplesheet)
    BWAMEM2_MEM(TRIMGALORE.out.reads, ch_bwa_index, ch_genome_fasta, true)
    BWAMEM2_MEM.out.bam
        .map{ meta, bam -> [ meta - meta.subMap('rg_id', 'flowcell'), bam ] }
        .groupTuple().branch{
            single: it[1].size() == 1
            multiple: it[1].size() > 1
            }
        .set{ bams }
 
    // If there are no samples to merge, skip the process
    SAMTOOLS_MERGE(bams.multiple, ch_genome_fasta.join(ch_genome_fasta_index), "bai")
    prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)

    SAMTOOLS_INDEX(prepared_bam)

    ch_bam_bai = prepared_bam.join(SAMTOOLS_INDEX.out.index)

    //UMI dedup
    BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS(ch_bam_bai, true, false)

    BAM_VARIANTCALLING_INTERVALS(
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bam,
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.index,
        ch_genome_fasta,
        ch_genome_fasta_index,
        ch_genome_dict,
        ch_intervals,
        ch_dbsnp,
        ch_dbsnp_index
    )

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

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
            name:  'dxnextflowfp_software_'  + 'mqc_'  + 'versions.yml',
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
                [id: 'dxnextflowfp'],
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
