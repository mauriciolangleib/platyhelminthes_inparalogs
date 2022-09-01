# reading inparalogs list
read_tsv('../../data_inparalogos/multigenes_genes_enlistados.inparalogos.tsv') -> multigenes_genes_enlistados.inparalogos
read_tsv('../../data_inparalogos/indexados.tsv') -> indexados.tibble

#la vieja #read_tsv('../../../indexado_final_inparalogos.tibble.tsv') %>% 
#la vieja #  tidyr::separate_rows(., Genes, sep = ',') %>% 
#la vieja #  dplyr::mutate(., Genes = Genes %>% str_replace_all(., ' *', '')) -> index_monophyletic_groups.tibble

read_tsv('../../../index_grupos_inparalogos/results/indexado_final_inparalogos.tibble.tsv') %>% 
  tidyr::separate_rows(., Genes, sep = ',') %>% 
  dplyr::mutate(., Genes = Genes %>% str_replace_all(., ' *', '')) -> index_monophyletic_groups.tibble

# loading expression data
rnaquant_finished.tibble = readr::read_tsv('../results/expression_tables/rnaquant_tpms.tsv', col_names = TRUE)

# loading metadata
samples_metadata.tibble = list.files('../data/run_metadata', recursive = T, pattern = '.tsv$', full.names = T) %>%
  as.list() %>%
  purrr::map_dfr(., ~{
    # getting Study
    study_tag = .x %>% str_split('/') %>% purrr::map_chr(5) %>% str_split('\\.') %>% purrr::map_chr(1)
    readr::read_tsv(.x, col_names = T) %>% dplyr::mutate(Study = study_tag) %>% dplyr::select(Study, Sample, Run, Condition) %>% unique()}
                )

# merging with expression data
multigenes_genes_enlistados.inparalogos %>%
  dplyr::full_join(x = ., y = indexados.tibble, by = c('Gene' = 'codigo')) %>%
  #dplyr::mutate(WBPS = WBPS %>% str_replace_all(., '$', '\\.1')) %>%
  dplyr::select(Gene, archivo, WBPS) %>%
  dplyr::right_join(x = ., y = rnaquant_finished.tibble, by = c('WBPS' = 'GeneID')) %>%
  dplyr::filter(!is.na(species) & !is.na(archivo)) %>%
  dplyr::left_join(x = ., y = index_monophyletic_groups.tibble, by = c('WBPS' = 'Genes')) -> rnaquant.tibble

rnaquant.tibble %<>% 
  group_by(monophyletic_group_code, WBPS, species, Condition, `Multimapping reads`) %>%
  dplyr::mutate(ymax = max(TPM), ymin = min(TPM)) %>%
  ungroup()

# creating dictionary run/study & run/sample
samples.dict = samples_metadata.tibble %>% dplyr::select(Study, Sample, Run) %>% unique() %>% .$Sample
names(samples.dict) = samples_metadata.tibble %>% dplyr::select(Study, Sample, Run) %>% unique() %>% .$Run

study.dict = samples_metadata.tibble %>% dplyr::select(Study, Sample, Run) %>% unique() %>% .$Study
names(study.dict) = samples_metadata.tibble %>% dplyr::select(Study, Sample, Run) %>% unique() %>% .$Run

# saving the employed table, with some modifications
rnaquant.tibble %>%
  dplyr::filter(`Multimapping reads` == 'Not allowed') %>%
  dplyr::select(-c(`Multimapping reads`, species, Gene, archivo, ymax, ymin)) %>%
  dplyr::rename(`Gene code` = 'WBPS',
                `Inparalogs Group Code` = 'monophyletic_group_code',
                `Run` = sample,
                `TPM per run` = 'TPM',
                `Median TPM among replicates` = 'TPM median') %>%
  dplyr::mutate(Sample = Run %>% tidytidbits::lookup_chr(dict = samples.dict, default = NA),
                Study = Run %>% tidytidbits::lookup_chr(dict = study.dict, default = NA)) -> rnaquant.tosave.tibble

rnaquant.tosave.tibble$Condition %<>% as.factor()
rnaquant.tosave.tibble$Condition %<>% forcats::fct_recode(., `Scolex-Neck` = 'Adult, scolex',
                                                     `Mid-section` = 'Adult, mid-section',
                                                     `Posterior` = 'Adult, posterior',
                                                     `Adult (whole)` = 'Adult, whole',
                                                     `Larvae (mid-metamorphosis)` = 'mid-metamorphosis larvae',
                                                     Miracidia = 'miracidia',
                                                    `First stage sporocyst` = 'first stage sporocyst',
                                                    `Second stage sporocyst` = 'second stage sporocyst',
                                                     Cercariae = 'Cercariae',
                                                    `Schistosomulae (3 hr)` = '3 hour schistosomule',
                                                    `Schistosomulae (24 hr)` = '24 hour schistosomule',
                                                    `Juvenile` = 'juvenile males and females NMRI',
                                                    `Adult` = 'adult males and females NMRI')  

rnaquant.tosave.tibble %>% readr::write_tsv(., '../results/expression_tables/rnaquant_tpms.toplot.tsv')

# creating directories
list('../results/expression_plots', '../results/expression_plots/smansoni', '../results/expression_plots/hmicrostoma') %>%
purrr::map(., ~{
  if (!dir.exists(.x)) {
   system(glue('mkdir {.x}'))
   }
})

# plotting
rnaquant.tibble %>%
  group_split(species) %>%
  purrr::map(., ~{
  # defining variable for .x
  expr_table = .x %>% dplyr::rename(`Dev. stage` = 'Condition')
    
  # creating pdf
    if(unique(expr_table$species) == 'hmicrostoma'){
      # creating a directory to allocate plots
      if(!dir.exists('hmicrostoma')){
        system('mkdir hmicrostoma')
        }
      
      # adjusting table and splitting in two sets: developmental stages and spatial data
      expr_table %>% dplyr::filter(str_detect(`Dev. stage`, 'mid-section|posterior|scolex')) -> expr_table.spatial
      expr_table %>% dplyr::filter(!str_detect(`Dev. stage`, 'mid-section|posterior|scolex')) -> expr_table.developmental
      
      expr_table.spatial$`Dev. stage` %<>% as.factor()

      expr_table.spatial$`Dev. stage` %<>% forcats::fct_recode(., `Scolex-Neck` = 'Adult, scolex',
                                                     `Mid-section` = 'Adult, mid-section',
                                                     `Posterior` = 'Adult, posterior')  
      
      expr_table.spatial$`Dev. stage` %<>% forcats::fct_relevel(., 'Scolex-Neck', 'Mid-section', 'Posterior')
      
      #expr_table.spatial %<>% tidyr::complete(., WBPS, archivo, `Dev. stage`, fill = list(TPM = 0, `TPM median` = 0)) %>% filter(!is.na(`Dev. stage`))

      
      expr_table.developmental$`Dev. stage` %<>% as.factor()

      expr_table.developmental$`Dev. stage` %<>% forcats::fct_recode(., `Adult (whole)` = 'Adult, whole',
                                                     `Larvae (mid-metamorphosis)` = 'mid-metamorphosis larvae')  

      expr_table.developmental$`Dev. stage` %<>% forcats::fct_relevel(., 'Larvae (mid-metamorphosis)', 'Adult (whole)')
      
      #expr_table.developmental %<>% tidyr::complete(., WBPS, archivo, `Dev. stage`, fill = list(TPM = 0, `TPM median` = 0)) %>% filter(!is.na(`Dev. stage`))

      
      } else if (unique(expr_table$species) == 'smansoni') {
      # creating a directory to allocate plots
      if(!dir.exists('smansoni')){
        system('mkdir smansoni')
        }
      
      # adjusting table
      expr_table$`Dev. stage` %<>% as.factor()
      
      expr_table$`Dev. stage` %<>% forcats::fct_recode(., Miracidia = 'miracidia',
                                        `First stage sporocyst` = 'first stage sporocyst',
                                        `Second stage sporocyst` = 'second stage sporocyst',
                                         Cercariae = 'Cercariae',
                                         `Schistosomulae (3 hr)` = '3 hour schistosomule',
                                         `Schistosomulae (24 hr)` = '24 hour schistosomule',
                                     `Juvenile` = 'juvenile males and females NMRI',
                                     `Adult` = 'adult males and females NMRI')  
  
      expr_table$`Dev. stage` %<>% forcats::fct_relevel(., 'Miracidia', 'First stage sporocyst', 'Second stage sporocyst', 'Cercariae',
                                                      'Schistosomulae (3 hr)', 'Schistosomulae (24 hr)',
                                                      'Juvenile', 'Adult')
      }
                 
  # plotting
  if(unique(.x$species) == 'smansoni'){
    expr_table %>%
      dplyr::filter(`Multimapping reads` == 'Not allowed') %>%
      group_split(monophyletic_group_code) %>%
      purrr::map(., ~{
        # plotting
        .x %>%
          ggplot(data = ., mapping = aes(x = `Dev. stage`, y = `TPM median`, color = WBPS, group = WBPS)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                legend.title = element_text()) +
          geom_line() +
          geom_point(mapping = aes(y = TPM, fill = WBPS, color = WBPS, group = WBPS)) +
          xlab('Developmental stage') + ylab('Expression level (TPM)') +
          geom_hline(yintercept=5, linetype="dashed", color = "indianred4", size=0.5) +
          labs(fill = 'Gene code (WormBase ParaSite)', color = 'Gene code (WormBase ParaSite)', group = 'Gene code (WormBase ParaSite)') +
          theme(legend.title = element_text(size = 7),
                legend.text = element_text(size = 6)) +
          #geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1) +
          ggtitle(glue('Monophyletic group: {unique(.x$monophyletic_group_code)}')) -> plot
        
        # printing
        ggsave(filename = glue('../results/expression_plots/smansoni/{unique(.x$monophyletic_group_code)}.developmental.pdf'), plot, device = 'pdf')
      })
                   
    } else if(unique(.x$species) == 'hmicrostoma'){
    expr_table.spatial %>%
      dplyr::filter(`Multimapping reads` == 'Not allowed') %>%
      group_split(monophyletic_group_code) %>%
      purrr::map(., ~{
        # plotting
        .x %>%
          ggplot(data = ., mapping = aes(x = `Dev. stage`, y = `TPM median`, color = WBPS, group = WBPS)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                legend.title = element_text()) +
          geom_line() +
          geom_point(mapping = aes(y = TPM, fill = WBPS, color = WBPS, group = WBPS)) +
          xlab('Spatial section') + ylab('Expression level (TPM)') +
          geom_hline(yintercept=5, linetype="dashed", color = "indianred4", size=0.5) +
          labs(fill = 'Gene code (WormBase ParaSite)', color = 'Gene code (WormBase ParaSite)', group = 'Gene code (WormBase ParaSite)') +
          theme(legend.title = element_text(size = 7),
                legend.text = element_text(size = 6)) +
          #geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1) +
          ggtitle(glue('Monophyletic group: {unique(.x$monophyletic_group_code)}')) -> plot.spat
        
        # printing
        ggsave(filename = glue('../results/expression_plots/hmicrostoma/{unique(.x$monophyletic_group_code)}.spatial.pdf'), plot.spat, device = 'pdf')
      })

    expr_table.developmental %>%
      dplyr::filter(`Multimapping reads` == 'Not allowed') %>%
      group_split(monophyletic_group_code) %>%
      purrr::map(., ~{
        # plotting
        .x %>%
          ggplot(data = ., mapping = aes(x = `Dev. stage`, y = `TPM median`, color = WBPS, group = WBPS)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                legend.title = element_text()) +
          geom_line() +
          geom_point(mapping = aes(y = TPM, fill = WBPS, color = WBPS, group = WBPS)) +
          xlab('Developmental stage') + ylab('Expression level (TPM)') +
          geom_hline(yintercept=5, linetype="dashed", color = "indianred4", size=0.5) +
          labs(fill = 'Gene code (WormBase ParaSite)', color = 'Gene code (WormBase ParaSite)', group = 'Gene code (WormBase ParaSite)') +
          theme(legend.title = element_text(size = 7),
                legend.text = element_text(size = 6)) +
          #geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1) +
          ggtitle(glue('Monophyletic group: {unique(.x$monophyletic_group_code)}')) -> plot.dev
        
        # printing
        ggsave(filename = glue('../results/expression_plots/hmicrostoma/{unique(.x$monophyletic_group_code)}.developmental.pdf'), plot.dev, device = 'pdf')
        })
    }        
  })

# creating unified plots and embeding them
system('pdfunite ../results/expression_plots/smansoni/*developmental* ../results/expression_plots/smansoni_developmental.pdf')
system('pdfunite ../results/expression_plots/hmicrostoma/*spatial* ../results/expression_plots/hmicrostoma_spatial.pdf')
system('pdfunite ../results/expression_plots/hmicrostoma/*developmental* ../results/expression_plots/hmicrostoma_developmental.pdf')
