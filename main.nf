#!/usr/bin/env nextflow

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { DOWNSAMPLE } from './modules/downsample.nf'
include { BBMERGE } from './modules/bbmerge.nf'
include { UNICYCLER } from './modules/unicycler.nf'
include { CIRCULARIZE } from './modules/circularize.nf'
include { BWA } from './modules/bwa.nf'
include { BAM_READ_COUNT } from './modules/bam_read_count.nf'
include { PLANNOTATE } from './modules/plannotate.nf'
include { FIX_START } from './modules/fix_start.nf'
include { SUMMARIZE } from './modules/summarize.nf'
include { ANALYZE_BAM_READ_COUNT } from './modules/analyze_bam_read_count.nf'
include { PLASMIDMAP } from './modules/plasmidMap.nf'

workflow {
  
    validateParameters()
    log.info paramsSummaryLog(workflow)

    fq_ch = Channel
          .fromFilePairs( params.input +  "/*_R{1,2}_001.fastq.gz",  flat: true )
           .map { it -> tuple( it[0], it[1], it[2])}
           
    downsample_out = DOWNSAMPLE( fq_ch )

    bbmerge_out = BBMERGE(downsample_out)

    unicycler_gfa = UNICYCLER(bbmerge_out.LR_ch)
    
    circular_out = CIRCULARIZE(unicycler_gfa)
    
    fix_start_out = FIX_START(circular_out.assembled_fa)
    
    align_in = bbmerge_out.LR_ch.join(fix_start_out)
   
    bam_out = BWA(align_in)
    
    metrics_out = bam_out.metrics.collect().mix(circular_out.circle_csv.collect()).collect()
     
    summary_output = SUMMARIZE(metrics_out)
    
    bam_read_count_report = BAM_READ_COUNT(bam_out.bam.join(fix_start_out))
    
    gbk = PLANNOTATE(fix_start_out)

    PLASMIDMAP(gbk.collect())
    
    ANALYZE_BAM_READ_COUNT(bam_read_count_report.collect())
}
