//
// VARSCAN somatic variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_MPILEUP as MPILEUP_NORMAL               } from '../../../modules/nf-core/samtools/mpileup/main'
include { SAMTOOLS_MPILEUP as MPILEUP_TUMOR                } from '../../../modules/nf-core/samtools/mpileup/main'
include { VARSCAN_SOMATIC                                  } from '../../../modules/nf-core/varscan/somatic/main'
include { VARSCAN_PROCESSSOMATIC as VARSCAN_PROCESS_SNV    } from '../../../modules/nf-core/varscan/processsomatic/main'
include { VARSCAN_PROCESSSOMATIC as VARSCAN_PROCESS_INDEL  } from '../../../modules/nf-core/varscan/processsomatic/main'
include { VARSCAN_FPFILTER as VARSCAN_FILTER_SNV           } from '../../../modules/nf-core/varscan/fpfilter/main'
include { VARSCAN_FPFILTER as VARSCAN_FILTER_INDEL         } from '../../../modules/nf-core/varscan/fpfilter/main'
include { BAMREADCOUNT                                     } from '../../../modules/nf-core/bamreadcount/main'

workflow BAM_VARIANT_CALLING_SOMATIC_VARSCAN {
    take:
    cram          // channel: [mandatory] [ meta, cram1, crai1, cram2, crai2 ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    _fasta_fai     // channel: [mandatory] [ meta, fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed ] or [ ] if no intervals

    main:
    versions = Channel.empty()

    // run samtools mpileup for normal and tumor samples
    ch_bam_for_mpileup_normal = cram
        .map{ meta, normal_cram, _normal_crai, _tumor_cram, _tumor_crai -> [ meta, normal_cram, intervals ] }
    ch_fasta_for_mpileup_normal = fasta.map{ _meta, fa -> [ fa ] }
    MPILEUP_NORMAL(ch_bam_for_mpileup_normal, ch_fasta_for_mpileup_normal)

    ch_bam_for_mpileup_tumor = cram
        .map{ meta, _normal_cram, _normal_crai, tumor_cram, _tumor_crai -> [ meta, tumor_cram, intervals ] }
    ch_fasta_for_mpileup_tumor = fasta.map{ _meta, fa -> [ fa ] }
    MPILEUP_TUMOR(ch_bam_for_mpileup_tumor, ch_fasta_for_mpileup_tumor)

    // run varscan somatic variant calling
    ch_mpileup_for_varscan_somatic = MPILEUP_NORMAL.out.mpileup
        .combine(MPILEUP_TUMOR.out.mpileup, by: 0)
    VARSCAN_SOMATIC(ch_mpileup_for_varscan_somatic)

    VARSCAN_PROCESS_SNV(VARSCAN_SOMATIC.out.vcf_snvs)
    VARSCAN_PROCESS_INDEL(VARSCAN_SOMATIC.out.vcf_indels)

    ch_bam_for_bamreadcount = cram.map{ meta, _normal_cram, _normal_crai, tumor_cram, tumor_crai -> [ meta, tumor_cram, tumor_crai ] }
    ch_fasta_for_bamreadcount = fasta.map{ _meta, fa -> [ fa ] }
    ch_bed_for_bamreadcount = intervals.map{ interval -> [ interval ] }
    BAMREADCOUNT(ch_bam_for_bamreadcount, ch_fasta_for_bamreadcount, ch_bed_for_bamreadcount)

    // filter variants
    ch_process_snv_to_filter = VARSCAN_PROCESS_SNV.out.somatic_hc_vcf
        .combine(BAMREADCOUNT.out.rc, by: 0)
    ch_process_indel_to_filter = VARSCAN_PROCESS_INDEL.out.somatic_hc_vcf
        .combine(BAMREADCOUNT.out.rc, by: 0)

    VARSCAN_FILTER_SNV(ch_process_snv_to_filter)
    VARSCAN_FILTER_INDEL(ch_process_indel_to_filter)

    snv_vcf   = VARSCAN_FILTER_SNV.out.pass_vcf
    indel_vcf = VARSCAN_FILTER_INDEL.out.pass_vcf

    snv_vcf = VARSCAN_SOMATIC.out.vcf_snvs
    indel_vcf = VARSCAN_SOMATIC.out.vcf_indels

    // add variantcaller to meta map
    vcf = Channel.empty().mix(snv_vcf, indel_vcf).map{ meta, vcf -> [ meta + [ variantcaller:'varscan' ], vcf ] }

    versions = versions.mix(MPILEUP_NORMAL.out.versions)
    versions = versions.mix(MPILEUP_TUMOR.out.versions)
    versions = versions.mix(VARSCAN_SOMATIC.out.versions)
    versions = versions.mix(VARSCAN_PROCESS_SNV.out.versions)
    versions = versions.mix(VARSCAN_PROCESS_INDEL.out.versions)
    versions = versions.mix(VARSCAN_FILTER_SNV.out.versions)
    versions = versions.mix(VARSCAN_FILTER_INDEL.out.versions)
    versions = versions.mix(BAMREADCOUNT.out.versions)

    emit:
    vcf

    versions
}
