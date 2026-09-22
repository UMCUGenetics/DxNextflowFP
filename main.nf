#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowFP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowFP
----------------------------------------------------------------------------------------


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { extractFastqPairFromDir } from './modules/local/utils/fastq.nf'

include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'

include { FASTQC } from './modules/nf-core/fastqc/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main'
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'

// Subworkflows nf-core/UMCUGenetics
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from './subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { BAM_VARIANTCALLING_INTERVALS } from './subworkflows/UMCUGenetics/bam_variantcalling_intervals/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Create reference file channels, add meta values
    ch_genome_fasta = Channel.fromPath("${params.genome_fasta}").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_genome_fasta_index = Channel.fromPath("${params.genome_fasta}.fai").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_genome_dict = Channel.fromPath("${params.genome_dict}").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_bwa_index = Channel.fromPath("${params.bwa_index}*").map{ file -> [file.getSimpleName(), file] }.groupTuple().collect()
    ch_dbsnp = Channel.fromPath("${params.dbsnp}").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_dbsnp_index = Channel.fromPath("${params.dbsnp}.tbi").map{ file -> [file.getSimpleName(), file] }.collect()
    ch_intervals = Channel.fromPath("${params.intervals}").map{ file -> [file.getSimpleName(), file] }.collect()


    // Input channel
    ch_fastq = extractFastqPairFromDir(params.input, params.outdir)
    def outdir = params.outdir

    // Trim FASTQs
    TRIMGALORE(ch_fastq)

    // Mapping
    BWAMEM2_MEM(TRIMGALORE.out.reads, ch_bwa_index, ch_genome_fasta, true)

    // Merge multiple lane samples and index
    BWAMEM2_MEM.out.bam
        .map{ meta, bam -> [ meta - meta.subMap('rg_id', 'flowcell'), bam ] }
        .groupTuple().branch{
            single: it[1].size() == 1
            multiple: it[1].size() > 1
            }
        .set{ bams }
 
    // If there are no samples to merge, skip the process
    SAMTOOLS_MERGE(bams.multiple, ch_genome_fasta.join(ch_genome_fasta_index))
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


    // QC
    FASTQC(ch_fastq)


    // MultiQC
    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
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
            name:  'Fingerprint_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )



    ch_multiqc_logo = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()


    
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    
    ch_multiqc_files = ch_multiqc_files
        .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.stats.collect{it[1]}.ifEmpty([]))
        .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.flagstat.collect{it[1]}.ifEmpty([]))
        .mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))
        .mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))
        .mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    
    def multiqc_config = "${projectDir}/assets/multiqc_config.yml"
    def multiqc_logo = ""


    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'fp'],
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
workflow.onComplete {
    def analysis_id = params.outdir.split('/')[-1]
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: analysis_id,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "FP Workflow Successful: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html, attach: "${params.outdir}/QC/multiqc_report.html")
    } else {
        def subject = "FP Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}*/
