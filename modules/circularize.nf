process CIRCULARIZE {
  
    tag "${pair_id}"
    publishDir path: "${params.output}/GFA", pattern: '*.final.gfa', mode: 'copy'
    publishDir path: "${params.output}/FASTA", pattern: '*.final.fasta', mode: 'copy'

    input:
    tuple val(pair_id), path(gfa)

    output:
     tuple val(pair_id), path("*.final.gfa"), path("*.final.fasta"), emit: circle_out
     path('*.info.csv'), emit: circle_csv
     tuple val(pair_id),  path("*.final.fasta"), emit: assembled_fa
     tuple val(pair_id), path("*.final.fasta"), emit: assembled_fa2

    script:
    """
    if [ -s $gfa ]; then
        Graph.py
    else
        touch ${pair_id}.final.gfa
        touch ${pair_id}.final.fasta
        touch ${pair_id}.info.csv
    fi
    """
}
