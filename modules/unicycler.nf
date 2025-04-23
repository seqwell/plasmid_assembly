process UNICYCLER {
    tag "${pair_id}"
  
    input:
    tuple val(pair_id),  path(norm_fq)


    output:
    tuple val(pair_id),  path('*.gfa') 

    script:
    """
    if [ -s $norm_fq ]; then
     
     python3 /opt/unicycler/unicycler-runner.py  \
       -s ${norm_fq} \
       -o ${pair_id} \
       --kmers 27,47,63,77,89,99,107,115,121,127 \
       --spades_path /opt/conda/envs/fastq-assembly/bin/spades.py \
       --tblastn_path /opt/conda/envs/fastq-assembly/bin/tblastn \
       --makeblastdb_path /opt/conda/envs/fastq-assembly/bin/makeblastdb \
       || touch ${pair_id}.gfa

     fi

     if [ -f ${pair_id}/assembly.gfa ]; then
        mv ${pair_id}/assembly.gfa ${pair_id}.gfa
     else
        touch ${pair_id}.gfa
     fi

    """
}
