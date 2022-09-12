# loading functions
source('../R_scripts/functions.R')

# loading furrr library
library(tidyverse)
library(furrr)
library(doParallel)

# loading table of GO terms and creating GO slim annotation if it doesnt exist
eggnog_cols = c('query_name',     
              'seed_eggNOG_ortholog',
              'seed_ortholog_evalue',    
              'seed_ortholog_score',     
              'predicted_gene_name',     
              'GO_terms',        
              'KEGG_pathways',   
              'Annotation_tax_scope',    
              'OGs',     
              'bestOG|evalue|score',     
              'COG cat', 
              'eggNOG annot')

eggnog_annot.tibble = read_tsv('../data/PROTEOMA_TOTAL.diamondCOGs.emapper.annotations',
          comment = '#',
          col_names = eggnog_cols)

if(!file.exists('../results/goterm2goslim_annotation.tsv')){
  # seteo la parte del computo en paralelo
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  # charging eggnog annotation
  go_slims.table = eggnog_annot.tibble %>%
  # selecting gene ontology terms and query names
  dplyr::select(GO_terms) %>%
  # applying map_to_slim function
  tidyr::separate_rows(., 'GO_terms', sep = ',') %>%
  unique() %>%
  #group_split(GO_terms) %>%
  list() %>%
  purrr::map_dfr(., ~{
  # separating rows by GO terms
    tablaGOs = .x 
    # applying map_to_slim
    tablaGOslims <- foreach(i=1:length(as.character(tablaGOs$GO_terms)),
	                .combine=rbind,
	                .packages=c('magrittr', 'reshape2', 'dplyr', 'stringr')
	                 ) %dopar% {
      
     map_to_slim = function(GOterm, go_obo, slim_obo){
      	linea <- paste("map_to_slim.py --term='", GOterm, "' --slim_out='direct' ", go_obo, " ", slim_obo, sep = "")
      	system(linea, intern = TRUE) %>%
      	grep(., pattern = "GO:", value = TRUE )-> resultado
      	return(resultado)
      }
      
	   # obtengo el codigo de la especie en base al nombre de archivo, y filtro en la tabla de metadata
	   map_to_slim(as.character(tablaGOs$GO_terms[i]),
				  go_obo="../data/mapeoAGOSlims/go.obo",
	                          slim_obo="../data/mapeoAGOSlims/goslim_generic.obo") -> go_slims

	   if (length(go_slims) != 0){
		   dplyr::mutate(tablaGOs[i,],
		 		go_slim = go_slims) -> temporal
	   }

	   if (length(go_slims) == 0){
	           dplyr::mutate(tablaGOs[i,],
	                        go_slim = "NA") -> temporal
	   }

	   temporal #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
	}
  
  })
  
  # stop cluster
  stopCluster(cl)   

  # saving table
  write_tsv(go_slims.table, '../results/goterm2goslim_annotation.tsv')
  }

# creating GOSlim annotation table 
goterm2goslim.tibble = readr::read_tsv('../results/goterm2goslim_annotation.tsv', col_names = T)
if(!file.exists('../results/PROTEOMA_TOTAL_goslim.annot.tsv')){
	## merging
	goslim_annotations.tibble = eggnog_annot.tibble %>%
		dplyr::select(query_name, GO_terms) %>%
		tidyr::separate_rows(., 'GO_terms', sep = ',') %>%
		dplyr::left_join(x = ., 
				 y = goterm2goslim.tibble %>%
				 	tidyr::separate_rows(., 'GO_terms', sep = ','),
				 by = c('GO_terms' = 'GO_terms')) %>%
		dplyr::select(query_name, go_slim) %>%
		group_split(query_name) %>%
		purrr::map_dfr(., ~{tibble(query_name = unique(.x$query_name),
				       GOSlim = paste(unique(.x$go_slim), collapse = ','))})
	# saving
	readr::write_tsv(goslim_annotations.tibble, '../results/PROTEOMA_TOTAL_goslim.annot.tsv')
}
