# analisis pero con los GO terms
library(readr)
library(dplyr)
library(magrittr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(glue)
library(tibble)
#library(viridis)
#library(KEGGprofile)
#library(KEGGREST)
library(ggplot2)
#library(pathfindR)
#library(clusterProfiler)
library(stringr)


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

read_tsv('../data/tablaConGoSlims.porEspecie.global.tsv', col_types = cols(col_character(),
                                                                   col_character(),
                                                                   col_character(), 
                                                                   col_character(),
                                                                   col_character(),
                                                                   col_character(), 
                                                                   col_character(),  
                                                                   col_character()
                                                                                )) %>%
  dplyr::rename(go_term = go_slim, go_slim = V1,
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
  rowwise() %>%
  dplyr::mutate(Species = case_when(Species == 'H. taeniformis' ~ 'H. taeniaeformis', Species == 'S. matthei' ~ 'S. mattheei', TRUE ~ Species)) %>%
  #filter(., go_term_ontology == 'MF') %>%
  dplyr::select(Species, go_slim_description, go_slim_ontology) %>% 
  dcast(go_slim_description + go_slim_ontology ~ Species, fill = 0) %>%
  as_tibble() -> data.desordenada

  data.desordenada$go_slim_ontology %<>% as.factor
  data.desordenada$go_slim_description %<>% as.factor
  
  data.desordenada %>%
    mutate_if(., is.character, ~ if_else(. == '0', 0, 1)) -> data.desordenada
  
  data.desordenada$go_slim_ontology %<>% as.character
  data.desordenada$go_slim_description %<>% as.character

  data.desordenada$go_slim_ontology %<>% str_replace_all(., pattern = 'BP', replacement = 'Biological process')
  data.desordenada$go_slim_ontology %<>% str_replace_all(., pattern = 'MF', replacement = 'Molecular function')
  data.desordenada$go_slim_ontology %<>% str_replace_all(., pattern = 'CC', replacement = 'Cellular component')
  
  opa = data.frame(`Ontology` = factor(x = data.desordenada$go_slim_ontology))
  rownames(opa) = data.desordenada$go_slim_description
  rownames_to_column(opa, 'je') %>%
    arrange(., Ontology) %>%
    column_to_rownames(., 'je') -> opa
  
  
  library(forcats)
  
  epa = data.frame(Class = factor(x = TablaCodigos$class))
  rownames(epa) = TablaCodigos$Species
  
  epa$Class %<>%forcats::fct_relevel(., c('Cestoda', 'Trematoda', 'Monogenea', 'Turbellaria'))
  #epa = data.frame(Bicho = TablaCodigos$Species, Class = factor(x = TablaCodigos$class))
  

  
  library(RColorBrewer)
  
  data.desordenada %>%
    #filter(., !descripcion == 'Metabolic pathways') %>%
    dplyr::select(-go_slim_ontology) %>%
    #mutate_if(., is.double, ~ if_else(. > 0, 1, 0)) %>% ------> Esta parte la saque para experimentar
    filter(., !go_slim_description %in% c('biological_process','molecular_function','cellular_component')) %>% 
    column_to_rownames('go_slim_description') -> hecho.casi
  
  #hecho.casi %>%
  #  mutate_if(., is.numeric, ~if_else(. > 0, 1, 0)) -> hecho
  hecho = hecho.casi
  rownames(hecho) = rownames(hecho.casi)
  hecho = hecho[rownames(opa), rownames(epa)] %>% rownames_to_column('momentaneo') %>% filter(!str_detect(momentaneo, 'NA_') & momentaneo != 'NA' & !is.na(`D. latum`)) %>% column_to_rownames('momentaneo')

  annot.color.col <- list('Class'=c(Cestoda = "dodgerblue3", 
                                    Trematoda = 'coral3',
                                    Monogenea = 'olivedrab4',
                                    Turbellaria = 'darkorange'),
                          'Ontology' = c(`Biological process` = 'chartreuse4',
                                                 `Molecular function` = 'darkorchid3', 
                                                 `Cellular component` = 'deepskyblue2'))
  
  categories <- data.frame(Sex = factor(sample(c("1", "2"),size = 10,replace = T),labels = c("Male", "Female")), Stage= factor(sample(c('Patient10','Patient9','Patient5'),size = 10,replace = T), labels = c('I','II','III')))
  
  orden_especies = c("S. mediterranea", "M. lignano", 
                     "P. xenopodis", "G. salaris",
                     "S. rodhaini", "S. mattheei", "S. margrebowiei", "S. mansoni", "S. japonicum", "S. haematobium", 
                     "S. curassoni", "T. regenti", "F. hepatica", "E. caproni", "O. viverrini", "C. sinensis",
                     "S. erinaceieuropaei", "S. solidus", "D. latum", "H. taeniaeformis",
                     "T. saginata", "T. asiatica", "T. solium", "E. canadensis", "E. multilocularis", "E. granulosus",
                     "M. corti", "H. nana", "H. microstoma", "H. diminuta")
  
  pheatmap(mat = hecho[,rev(orden_especies)], 
           #annCol = epa, 
           annotation_col = epa,
           width = (15*1), 
           height = (20*1), 
           filename = '../results/GO_Terms_presenciaAusencia.pdf',
           #annRow = opa, 
           annotation_row = opa, 
           annotation_names_row = F, annotation_names_col = F,
           border_color = 'black', #LOSAQUEPARAEXPERIMENTARcolor = c('white', 'tomato1'), 
           #labCol = 1:30,
           #             color = 'heat',
           fontsize = 10, #annLegend = FALSE,
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
      
