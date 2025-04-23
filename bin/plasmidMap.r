#!/usr/bin/env Rscript

install.packages("devtools")
devtools::install_github("bradyajohnston/plasmapr")

library(plasmapR)
library(dplyr)
library(ggplot2)
library(purrr)


#https://bradyajohnston.github.io/plasmapR/

get_pair_id = function(x){
  
  file_base_name = basename(x)
  pair_id = gsub(".gbk", "", file_base_name)
  return(pair_id)
  
}


file_path_list = list.files( pattern = "gbk")
pair_id_list = purrr::map_chr( file_path_list, get_pair_id)

print(file_path_list)
print(pair_id_list)


get_plasmidMap = function(file_path, pair_id){
  
  plasmid_data <- read_gb(file_path)
  
  plot_plasmid(plasmid_data,name = pair_id) +
    theme_bw() +  
    theme(panel.background = element_rect(fill = "white", color = "black")) 
  
  ggsave( paste0( pair_id, ".plasmidMap.png"), plot = last_plot(), width = 10, height = 10, dpi = 300)
}


purrr::walk2(file_path_list, pair_id_list, get_plasmidMap )




