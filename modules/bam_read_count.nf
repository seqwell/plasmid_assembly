process BAM_READ_COUNT {
  tag "${pair_id}"
  
  input:
  tuple val(pair_id), path(bam), path(fa) 
  
  output:
  path("*.txt")  
  
  """
  bam-readcount -f $fa -w 0 $bam > ${pair_id}.bam.readcount.txt

  """
  
}
