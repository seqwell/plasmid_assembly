process PLASMIDMAP {
    
   
    publishDir path: "${params.output}/plasmidMap", pattern: '*.png', mode: 'copy'

    input:
    path(gbk_file) 

    output:
    path("*.png")

    script:
    """
    plasmidMap.r
    """
}
