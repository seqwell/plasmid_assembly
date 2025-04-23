process FIX_START {
  tag "${pair_id}"
  
  input:
  tuple val( pair_id), path(fa) 
  
  output:
  tuple val( pair_id), path("${pair_id}.fixstart.fasta")  
  
  """
  circlator fixstart $fa  ${pair_id}.fixstart

  """
}
