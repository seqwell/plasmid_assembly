process BWA {
    tag "${pair_id}"
  
    input:
    tuple val(pair_id),  path(fq),  path(fa)

    output:
     tuple val(pair_id), path("*bam"), emit: bam
     path ("*.csv"), emit: metrics
     tuple val(pair_id), path("*.position.read.coverage.txt"), emit: cov_ch

    script:
    """
    if [ -s $fa ]; then
        bwa index $fa
        bwa mem $fa $fq | samtools view -bh -F2048 - | samtools sort > ${pair_id}.bam
        samtools depth -a ${pair_id}.bam > ${pair_id}.depth.csv
        samtools view -c ${pair_id}.bam >${pair_id}.count.csv
        samtools mpileup -B -q 0 -Q 0 -f $fa ${pair_id}.bam | awk '{print \$2, \$3, \$4}' > ${pair_id}.position.read.coverage.txt
    else
        touch ${pair_id}.depth.csv
        touch ${pair_id}.count.csv
        touch ${pair_id}.bam
        touch ${pair_id}.position.read.coverage.txt
    fi
    """
}
