include { GATK4_GENOTYPEGVCFS   } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'

workflow BAM_FP {
    take:
    ch_bam        // channel: [ val(meta), [ bam ] ]
    ch_bai        // channel: [ val(meta), [ bai ] ]
    ch_ref        // channel: [ val(meta), [ ref.fasta ] ]
    ch_ref_index  // channel: [ val(meta), [ ref.fasta.fai ] ]
    ch_ref_dict   // channel: [ val(meta), [ ref.dict ] ]
    ch_intervals  // channel: [ val(meta), [ intervals ] ]
    ch_db_snp     // channel: [ val(meta), [ db_snp.vcf ] ]
    ch_db_snp_tbi // channel: [ val(meta), [ db_snp.vcf.tbi ] ]

    main:
    // Create channel for haplotypecaller with meta, bam, bai, intervals and empty optional. Discard meta from intervals.
    ch_bam
      .join(ch_bai)
      .combine(ch_intervals)
      .map{meta, bam, bai, meta2, intervals -> 
        return tuple([meta, bam, bai, intervals, [] ])
      }
      .set{ ch_bam_bai_intervals }

    GATK4_HAPLOTYPECALLER(
        ch_bam_bai_intervals,
        ch_ref,
        ch_ref_index,
        ch_ref_dict,
        ch_db_snp,
        ch_db_snp_tbi
    )
    // Create channel for Genotypegvcfs by using meta, vcf, index, intervals and empty optional.
    GATK4_HAPLOTYPECALLER.out.vcf
      .join(GATK4_HAPLOTYPECALLER.out.tbi)
      .combine(ch_intervals)
      .map{meta, vcf, tbi, meta2, intervals ->
        return tuple([meta, vcf, tbi, intervals, [] ])
      }
      .set{ ch_vcf_tbi_intervals }

    GATK4_GENOTYPEGVCFS(
      ch_vcf_tbi_intervals,
      ch_ref,
      ch_ref_index,
      ch_ref_dict,
      ch_db_snp,
      ch_db_snp_tbi
    )

    emit:
    fp_vcf        = GATK4_GENOTYPEGVCFS.out.vcf           // channel: [ val(meta), vcf.gz ]
    fp_vcf_tbi    = GATK4_GENOTYPEGVCFS.out.tbi           // channel: [ val(meta), vcf.gz.tbi ]

}
