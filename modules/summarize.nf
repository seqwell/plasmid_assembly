process SUMMARIZE {
  publishDir path: "${params.output}/", pattern: '*.csv', mode: 'copy'
  publishDir path: "${params.output}/COVERAGE", pattern: '*.png', mode: 'copy'

  input:
  path(metrics) 

  output:
  path("*") 
  
	"""
	SummarizeAssembly.py 
	"""
}
