include { GATK4_GENOTYPEGVCFS   } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'

workflow BAM_VARIANTCALLING_INTERVALS {
    take:
    ch_bam        // channel: [ meta, bam ]
    ch_bai        // channel: [ meta, bai ]
    ch_ref        // channel: [ meta, ref.fasta ]
    ch_ref_index  // channel: [ meta, ref.fasta.fai ]
    ch_ref_dict   // channel: [ meta, ref.dict ]
    ch_intervals  // channel: [ meta, intervals ]
    ch_db_snp     // channel: [ meta, db_snp.vcf ]
    ch_db_snp_tbi // channel: [ meta, db_snp.vcf.tbi ]

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
    vcf    = GATK4_GENOTYPEGVCFS.out.vcf           // channel: [ meta, vcf.gz ]
    tbi    = GATK4_GENOTYPEGVCFS.out.tbi           // channel: [ meta, vcf.gz.tbi ]

}
