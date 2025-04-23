
process DOWNSAMPLE {
  
     tag "${pair_id}"

     input:
     tuple val(pair_id),  path(read1), path(read2) 

     output:
     tuple val(pair_id),  path('*_R1_001.fastq.gz'), path('*_R2_001.fastq.gz') 

     """
     if [ $params.downsample -gt 0 ]; then
     seqtk sample -s 14 ${read1} $params.downsample | gzip > ${pair_id}.${params.downsample}_R1_001.fastq.gz
     seqtk sample -s 14 ${read2} $params.downsample | gzip > ${pair_id}.${params.downsample}_R2_001.fastq.gz
     else
     ln -s ${read1} ${pair_id}_full_R1_001.fastq.gz
     ln -s ${read2} ${pair_id}_full_R2_001.fastq.gz
     fi
     """
}
