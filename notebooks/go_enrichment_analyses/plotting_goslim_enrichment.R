# loading libraries
library(tidyverse)
library(magrittr)
library(reshape2)
library(glue)
#library(corrplot)

# loading functions
source('functions.R')

# loading species code table
tabla_codigos = readr::read_tsv('../../../TablaCodigos.tsv', col_names = TRUE)

# loading results
#goslim_table = tibble(files = list.files('../results/goslim_enrichment_for_organisms', pattern = '.csv$', full.names = TRUE)) %>%
#  dplyr::mutate(ontology = files %>% str_split(., '/') %>% purrr::map_chr(4) %>%
#                                      str_split(., '\\.') %>% purrr::map_chr(2),
#               organism_code = files %>% str_split(., '/') %>% purrr::map_chr(4) %>%
#                                      str_split(., '\\.') %>% purrr::map_chr(1) %>%
#                                      str_replace_all(., '_goslims', '')) %>%
#  dplyr::left_join(x = .,
#                  y = tabla_codigos,
#                  by = c('organism_code' = 'code')) %>%
#  dplyr::mutate(., charging = purrr::pmap(list(files, ontology, organism_code, Species, class), ~{
#    # defining variables
#    files = ..1
#    ontology = ..2
#    organism_code = ..3
#    Species = ..4
#    class = ..5
#    
#    # parsing enriched GOs
#    readr::read_csv(files, 
#                    col_names = c('GOSlim', 'Description', 'Background count', 'Count', 'Expected count', 'p_value'),
#                    col_types = c(col_character(), col_character(), col_double(), col_double(), col_double(), col_character())
#                   ) %>%
#    dplyr::mutate(p_value = p_value %>% str_replace_all(., '< ', '') %>% as.numeric())
#  })
#                  ) %>%
#  unnest() %>%
#dplyr::select(ontology, Species, GOSlim, Description, p_value) %>%
#dplyr::mutate(p_value = -log10(p_value)) %>%
#tidyr::pivot_wider(., id_cols = c('ontology', 'GOSlim', 'Description'), names_from = 'Species', values_from = 'p_value')

# loading results
read_tsv('../data/tablaConGoSlims.porEspecie.global.tsv', col_types = cols(col_character(),
                                                                   col_character(),
                                                                   col_character(), 
                                                                   col_character(),
                                                                   col_character(),
                                                                   col_character(), 
                                                                   col_character(),  
                                                                   col_character()
                                                                                )) %>%
  dplyr::rename(go_term = V1, go_slim = go_slim,
                Species = species,
                go_slim_description = TERM.y,
                go_slim_ontology = ONTOLOGY.x,
                go_term_description = TERM.x,
                go_term_ontology = ONTOLOGY.y) %>%
  mutate(Species = Species %>% purrr::map_chr(., ~{
                                                    .x %>% str_split(., '_') %>% unlist() -> texto
                                                    dos = texto[2] 
                                                    texto[1] %>% str_to_upper() %>% stringr::str_trunc(., width = 2, 'right', '.') -> uno
                                                    glue('{uno} {dos}')
                                                  }),
         go_term_description = glue('{go_term} - {go_term_description}')
       ) %>%
  #filter(., go_term_ontology == 'MF') %>%
  dplyr::mutate(Species = Species %>% str_replace_all(., 'H. taeniformis', 'H. taeniaeformis') %>% str_replace_all(., 'S. matthei', 'S. mattheei')) %>%
  dplyr::select(Species, go_slim, go_slim_description, go_slim_ontology) %>% 
  dcast(go_slim + go_slim_description + go_slim_ontology ~ Species, fill = 0) %>%
  as_tibble() %>%
  dplyr::rename(Description = 'go_slim_description', ontology = 'go_slim_ontology', GOSlim = 'go_slim') %>%
  mutate_if(., is.double, ~ if_else(. == 0, 0, 1)) %>%
  dplyr::mutate(ontology = ontology %>% str_replace_all(., pattern = 'BP', replacement = 'Biological process') %>%
                                        str_replace_all(., pattern = 'MF', replacement = 'Molecular function') %>%
                                        str_replace_all(., pattern = 'CC', replacement = 'Cellular component')) %>%
  dplyr::mutate(ontology = ontology %>% as.factor(), Description = Description %>% as.factor()) %>%
  dplyr::filter(!Description %in% c('biological_process','molecular_function','cellular_component'))  -> goslim_table
  
# creating matrix
goslim_matrix = goslim_table %>% dplyr::select(-c(ontology,GOSlim)) %>% column_to_rownames('Description') %>% as.matrix
rownames(goslim_matrix) = goslim_table %>% dplyr::mutate(goslim_description = glue('{GOSlim}-{Description}')) %>% .$goslim_description

# creating directory
if(!dir.exists('../results/plots')){
  system('mkdir ../results/plots')
  }

# saving table
#readr::write_csv(goslim_table, '../results/plots/goslim_pvalues.csv', col_names = TRUE)

# loading table
#goslim_enrichment.tibble = readr::read_csv('../results/plots/goslim_pvalues.csv', col_names = TRUE)

# plotting
# creating annotation variables
annot.color.col <- list('Class'=c(Cestoda = "dodgerblue3", 
                                  Trematoda = 'coral3',
                                  Monogenea = 'olivedrab4',
                                  Turbellaria = 'darkorange'),
                        'Ontology' = c(`Biological process` = 'chartreuse4',
                                       `Molecular function` = 'darkorchid3', 
                                       `Cellular component` = 'deepskyblue2'))

orden_especies = c("S. mediterranea", "M. lignano", 
                   "P. xenopodis", "G. salaris",
                   "S. rodhaini", "S. mattheei", "S. margrebowiei", "S. mansoni", "S. japonicum", "S. haematobium", 
                   "S. curassoni", "T. regenti", "F. hepatica", "E. caproni", "O. viverrini", "C. sinensis",
                   "S. erinaceieuropaei", "S. solidus", "D. latum", "H. taeniaeformis",
                   "T. saginata", "T. asiatica", "T. solium", "E. canadensis", "E. multilocularis", "E. granulosus",
                   "M. corti", "H. nana", "H. microstoma", "H. diminuta")

TablaCodigos = tibble(Species = c("D. latum", "E. canadensis", "E. granulosus",
                                  "E. multilocularis", "H. taeniaeformis", "H. diminuta",
                                  "H. microstoma", "H. nana", "M. corti",
                                  "S. solidus", "S. erinaceieuropaei", "T. asiatica",
                                  "T. saginata", "T. solium", "C. sinensis", "E. caproni",
                                  "F. hepatica", "O. viverrini", "S. curassoni", "S. haematobium",
                                  "S. japonicum", "S. mansoni", "S. margrebowiei",
                                  "S. mattheei", "S. rodhaini", "T. regenti",
                                  "G. salaris", "P. xenopodis", "M. lignano", "S. mediterranea"),
                      code = c("dpy", "ecd", "Ecg", "ecm", "htf", "hyd", "hym", "hyn",
                               "msc", "scs", "spe", "taa", "tas", "tae", "cln", "ecp",
                               "Fsc", "opv", "scc", "sch", "scj", "scm", "scw", "sct",
                               "scr", "trc", "gys", "ppx", "Msl", "stm"),
                      class = c(rep("Cestoda", 14), rep("Trematoda", 12), rep("Monogenea", 2), rep("Turbellaria", 2))
)

epa = data.frame(Class = factor(x = TablaCodigos$class))
rownames(epa) = TablaCodigos$Species

epa$Class %<>%forcats::fct_relevel(., c('Cestoda', 'Trematoda', 'Monogenea', 'Turbellaria'))

opa = data.frame(`Ontology` = goslim_table$ontology %>%as.factor())

rownames(opa) =  goslim_table %>% 
  dplyr::mutate(goslim_description = glue('{GOSlim}-{Description}')) %>% 
  .$goslim_description

orden_filas = opa %>%
  rownames_to_column('GOslim description') 

orden_filas = bind_rows(orden_filas %>% dplyr::filter(Ontology == 'Biological process'),
                        orden_filas %>% dplyr::filter(Ontology == 'Molecular function'),
                        orden_filas %>% dplyr::filter(Ontology == 'Cellular component'))

orden_filas = orden_filas$`GOslim description`


#hecho = hecho.casi
#rownames(hecho) = rownames(hecho.casi)
#hecho = hecho[rownames(opa), rownames(epa)] %>% rownames_to_column('momentaneo') %>% filter(!str_detect(momentaneo, 'NA_') & momentaneo != 'NA') %>% column_to_rownames('momentaneo')
  

library(pheatmap)

pheatmap(mat = goslim_matrix[orden_filas,rev(orden_especies)], 
         #annCol = epa, 
         annotation_col = epa,
         width = (19*1.3), 
         height = (27*1.3), 
         cellheight = 15,
         cellwidth = 24,
         filename = '../results/plots/GOSlim_significance_presenciaAusencia.pdf',
         #annRow = opa, 
         annotation_row = opa, 
         annotation_names_row = F, annotation_names_col = F,
         border_color = 'black', 
         #color = c('#ff2d00', '#fa461e', '#f55832', '#ef6744', '#e87555', '#df8166', '#d68d77', '#cb9889', '#bea29a'), 
	 color = c('white', '#cb9889'),
         na_col = 'white',
         #labCol = 1:30,
         #             color = 'heat',
         fontsize = 15, #annLegend = FALSE,
         #Colv = NA,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         angle_col = c("90"),
         annotation_colors = annot.color.col, 
         #display_numbers = T, number_color = 'white',fontsize_number = 2,
         #show_rownames = F, show_colnames = F,
         #annotation_colors = annot.color.col,
         #annColors = annot.color.col,
         legend = T)

