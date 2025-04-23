process ANALYZE_BAM_READ_COUNT {
  
   
  publishDir path: "${params.output}/per_base_data", pattern: '*.csv', mode: 'copy'
  
  input:
  path(read_count) 
  
  output:
  path("*.csv") 
  
  """
  
  ls | grep bam.readcount.txt > read_count_file_list
  
  cat read_count_file_list | sed 's|.bam.readcount.txt||' > file_name
  
  paste file_name read_count_file_list > info
  
  while read pair_id read_count_file; do
  
  analyze_bam_read_count.r \$pair_id \$read_count_file
  
  done < info
  
  
  """

  
}

