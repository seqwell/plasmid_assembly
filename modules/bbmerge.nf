process BBMERGE {
  
    tag "${pair_id}"
    
    publishDir path: "${params.output}/norm_fastq", pattern: '*.norm.merged.fq', mode: 'copy'
    input:
    tuple val(pair_id),  path(read1), path(read2) 

    output:
    tuple val(pair_id), path("*.norm.merged.fq"), emit: LR_ch
   


    """
    bbduk.sh -Xmx2G \
    in=$read1 \
    in2=$read2 \
    out=${pair_id}.clean.fastq.gz \
    ref=/opt/conda/envs/fastq-assembly/opt/bbmap-38.90-3/resources/nextera.fa.gz \
    k=21 mink=15 ktrim=n hdist=2 tpe tbo \
    overwrite=true

    bbmerge-auto.sh -Xmx2G \
    in=${pair_id}.clean.fastq.gz \
    out=${pair_id}.merged.fq \
    extend2=20 iterations=10 \
    k=60 ecc=true \
    ecctadpole=true \
    reassemble=true \
    rem=true \
    merge=true \
    strict=true

    bbnorm.sh -Xmx2G \
    in=${pair_id}.merged.fq \
    target=100 \
    maxdepth=150 \
    fixspikes=t \
    overwrite=true \
    out=${pair_id}.norm.merged.fq

    

    """
 
}

