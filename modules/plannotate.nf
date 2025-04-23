process PLANNOTATE {
    tag "${pair_id}"
    
    publishDir path: "${params.output}/GBK", pattern: '*.gbk', mode: 'copy'
 
    input:
    tuple val(pair_id), path(fa)

    output:
    path("*.gbk") 

    script:
    """
    mkdir ${pair_id}
    plannotate batch -i $fa --output ${pair_id}
    cp ${pair_id}/${pair_id}*pLann.gbk  ${pair_id}.gbk
    """
}
