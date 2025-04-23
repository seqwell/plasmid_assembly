#!/usr/bin/bash



input='./inputs'
output='./outputs/assembly_out'
downsample=0


nextflow run  main.nf \
-profile docker \
-work-dir './outputs/work/' \
--input $input \
--output $output \
--downsample $downsample \
-resume  -bg 
        
        
