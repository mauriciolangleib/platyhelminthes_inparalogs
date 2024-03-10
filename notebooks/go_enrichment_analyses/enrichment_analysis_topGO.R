# loading libraries
library(tidyverse)
library(magrittr)
library(glue)
library(furrr)
library(topGO)

# creating directory to allocate goslim for organisms
if(!dir.exists('../results/goslim_annot_for_organisms')){
  system('mkdir ../results/goslim_annot_for_organisms')
  }

# creating directory to allocate gene universes
if(!dir.exists('../results/gene_universes_for_organisms')){
  system('mkdir ../results/gene_universes_for_organisms')
  }

if(!dir.exists('../results/goslim_enrichment_for_organisms')){
  system('mkdir ../results/goslim_enrichment_for_organisms')
  }

# creating table with files to be used
files_structure.tibble = tibble(genes_of_interest = list.files('../../topGO_nuevo', pattern = 'genesInparalogos.*')) %>%
dplyr::mutate(., organism_code = genes_of_interest %>%
                                  str_split(., '\\.') %>%
                                  purrr::map_chr(2),
                 output_file = glue('../results/enrichment_results/{organism_code}.tsv'))

# running analyses with topGO under 'weight01' algorithm
files_structure.tibble %>%
  dplyr::transmute(., enrichment_run = purrr::pmap(list(genes_of_interest, organism_code, output_file), ~{
    # setting variables
    genes_of_interest = ..1
    organism_code = ..2
    output_file = ..3
    
    # loading inparalogs
    readLines(glue('../../topGO_nuevo/{genes_of_interest}')) -> genes.inparalogos
    # subsetting total GOslims for the organism
    readr::read_tsv('../results/PROTEOMA_TOTAL_goslim.annot.tsv', col_names = TRUE) %>%
      dplyr::filter(str_detect(query_name, organism_code)) %>%
    # saving total GOslims in the organism
      readr::write_tsv(., glue('../results/goslim_annot_for_organisms/{organism_code}.go_slim_annot.tsv'), col_names = FALSE)
    # loading table with GOSlim annotation for the organism
    readMappings(glue('../results/goslim_annot_for_organisms/{organism_code}.go_slim_annot.tsv')) -> anotacionGO
    
    # reading gene universe in the organism
    readr::read_tsv('../results/PROTEOMA_TOTAL_goslim.annot.tsv', col_names = TRUE) %>%
      dplyr::filter(str_detect(query_name, organism_code)) %>%
      .$query_name %>%
      unique() %>%
      writeLines(text = ., glue('../results/gene_universes_for_organisms/{organism_code}.gene_universe.txt'))
    readLines(glue('../results/gene_universes_for_organisms/{organism_code}.gene_universe.txt')) -> totalGenes
    geneList <- factor(as.integer(totalGenes %in% genes.inparalogos))
    names(geneList) <- totalGenes
    
    # running for ontologies
    ontologies = c('BP', 'MF', 'CC') 
    for(i in seq_along(ontologies)){
      # generate topGOdata object
      sampleGOdata = new("topGOdata", ontology=ontologies[i], allGenes=geneList, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = anotacionGO)
      
      # performing analysis
      resultFisher = runTest(sampleGOdata, algorithm = "weight01", statistic = "fisher")
      
      # obtaining data frame
      GOtermsdeInteres.data.frame = GenTable(sampleGOdata, 
                                             weightFisher = resultFisher, 
                                             orderBy = "weightFisher", 
                                             ranksOf = "weightFisher", 
                                             topNodes = sum(score(resultFisher) < .05))
      # saving result
      readr::write_csv(GOtermsdeInteres.data.frame %>% as_tibble(), 
                  col_names=FALSE, 
                  path=glue('../results/goslim_enrichment_for_organisms/{organism_code}_goslims.{ontologies[i]}.csv'))
      }
    
    
                                               })
               )
