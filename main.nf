#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowFP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowFP
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; } from 'plugin/nf-validation'
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { extractFastqPairFromDir } from './modules/local/utils/fastq.nf'

include { FASTQC } from './modules/nf-core/fastqc/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'
include { GATK4_HAPLOTYPECALLER } from './modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOTYPEGVCFS } from './modules/nf-core/gatk4/genotypegvcfs/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_UMITOOLS } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'
include { UMITOOLS_DEDUP } from './modules/nf-core/umitools/dedup/main' 
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'

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
    ch_intervals = Channel.fromPath("${params.intervals}").collect()


    // Input channel
    ch_fastq = extractFastqPairFromDir(params.input, params.outdir)

    // Trim FASTQs
    TRIMGALORE(ch_fastq)

    // Mapping
    BWAMEM2_MEM(TRIMGALORE.out.reads, ch_bwa_index, true)

    // Merge multiple lane samples and index
    BWAMEM2_MEM.out.bam
        .map{ meta, bam -> [ meta - meta.subMap('rg_id', 'flowcell'), bam ] }
        .groupTuple().branch{
            single: it[1].size() == 1
            multiple: it[1].size() > 1
            }
        .set{ bams }
 
    // If there are no samples to merge, skip the process
    SAMTOOLS_MERGE(bams.multiple, ch_genome_fasta, ch_genome_fasta_index)
    prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)
    SAMTOOLS_INDEX(prepared_bam)

    ch_bam_bai = prepared_bam.join(SAMTOOLS_INDEX.out.bai)

    //UMI dedup
    UMITOOLS_DEDUP(ch_bam_bai, true)
    SAMTOOLS_INDEX_UMITOOLS(UMITOOLS_DEDUP.out.bam)
    ch_bam_bai = UMITOOLS_DEDUP.out.bam.join(SAMTOOLS_INDEX_UMITOOLS.out.bai) 

    // Variant calling
    GATK4_HAPLOTYPECALLER(
        ch_bam_bai.combine(ch_intervals).map{ meta, bam, bai, intervals -> [meta, bam, bai, intervals, [] ] },
        ch_genome_fasta, ch_genome_fasta_index, ch_genome_dict, ch_dbsnp, ch_dbsnp_index
    )

    GATK4_GENOTYPEGVCFS(
        GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi).combine(ch_intervals).map{
            meta, vcf, tbi , intervals -> [meta, vcf, tbi, intervals, []]
        },
        ch_genome_fasta.map{ meta, file -> [file] },
        ch_genome_fasta_index.map{ meta, file -> [file] },
        ch_genome_dict.map{ meta, file -> [file] },
        ch_dbsnp.map{ meta, file -> [file] },
        ch_dbsnp_index.map{ meta, file -> [file] }
    )


    // QC
    FASTQC(ch_fastq)

    // Softare versions
    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_UMITOOLS.out.versions)
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))


    // MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        Channel.empty().toList(),
        Channel.empty().toList()
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
}
