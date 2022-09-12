library(tidyverse)
library(magrittr)
#library(plotly)
library(glue)

# creating supplementary figure 7
# loading data
dir_tabla_correlaciones = '../../RNASeq/vPaper/results/molecular_evolution_correlation_plots'
tabla_correlaciones = readr::read_tsv(glue('{dir_tabla_correlaciones}/data_supplementary_table_5.tsv'), col_names = TRUE)

# plotting average cophenetic synonymous distance vs Pearson's r (supp. fig 7a)
tabla_plot_a = tabla_correlaciones %>%
  dplyr::select(`Inparalogs Group Code`, `Gen X`, `Gen Y`, `Organism`, `Pearson r`, `Cophenetic dS`) %>%
  unique() %>%
  dplyr::mutate(Organism = Organism %>% str_replace_all(., 'S.mansoni', 'S. mansoni')) %>%
  dplyr::filter(`Organism` == 'S. mansoni') %>%
  dplyr::filter(!is.na(`Pearson r`)) %>%
  dplyr::mutate(`Pearson r` = case_when(between(`Pearson r`, -1, -0.8) ~ '-1 .. -0.8', 
                                            between(`Pearson r`, -0.8, -0.6) ~ '-0.8 .. -0.6', 
                                            between(`Pearson r`, -0.6, -0.4) ~ '-0.6 .. -0.4', 
                                            between(`Pearson r`, -0.4, -0.2) ~ '-0.4 .. -0.2', 
                                            between(`Pearson r`, -0.2, 0.0) ~ '-0.2 .. 0.0', 
                                            between(`Pearson r`, 0.0, 0.2) ~ '0.0 .. 0.2', 
                                            between(`Pearson r`, 0.2, 0.4) ~ '0.2 .. 0.4', 
                                            between(`Pearson r`, 0.4, 0.6) ~ '0.4 .. 0.6', 
                                            between(`Pearson r`, 0.6, 0.8) ~ '0.6 .. 0.8', 
                                            between(`Pearson r`, 0.8, 1.0) ~ '0.8 .. 1.0' 
                                            )
                )

tabla_plot_a$`Pearson r` %<>% as.factor
tabla_plot_a$`Pearson r` %<>% forcats::fct_relevel('-1 .. -0.8', '-0.8 .. -0.6', '-0.6 .. -0.4', '-0.4 .. -0.2', '-0.2 .. 0.0', '0.0 .. 0.2', '0.2 .. 0.4', '0.4 .. 0.6', '0.6 .. 0.8', '0.8 .. 1.0')

tabla_plot_b = tabla_correlaciones %>%
  dplyr::select(`Inparalogs Group Code`, `Gen X`, `Gen Y`, `Organism`, `Spearman r`, `Cophenetic dS`) %>%
  unique() %>%
  dplyr::mutate(Organism = Organism %>% str_replace_all(., 'S.mansoni', 'S. mansoni')) %>%
  dplyr::filter(`Organism` == 'S. mansoni') %>%
  dplyr::filter(!is.na(`Spearman r`)) %>%
  dplyr::mutate(`Spearman r` = case_when(between(`Spearman r`, -1, -0.8) ~ '-1 .. -0.8', 
                                        between(`Spearman r`, -0.8, -0.6) ~ '-0.8 .. -0.6', 
                                        between(`Spearman r`, -0.6, -0.4) ~ '-0.6 .. -0.4', 
                                        between(`Spearman r`, -0.4, -0.2) ~ '-0.4 .. -0.2', 
                                        between(`Spearman r`, -0.2, 0.0) ~ '-0.2 .. 0.0', 
                                        between(`Spearman r`, 0.0, 0.2) ~ '0.0 .. 0.2', 
                                        between(`Spearman r`, 0.2, 0.4) ~ '0.2 .. 0.4', 
                                        between(`Spearman r`, 0.4, 0.6) ~ '0.4 .. 0.6', 
                                        between(`Spearman r`, 0.6, 0.8) ~ '0.6 .. 0.8', 
                                        between(`Spearman r`, 0.8, 1.0) ~ '0.8 .. 1.0' 
  )
  )

tabla_plot_b$`Spearman r` %<>% as.factor
tabla_plot_b$`Spearman r` %<>% forcats::fct_relevel('-1 .. -0.8', '-0.8 .. -0.6', '-0.6 .. -0.4', '-0.4 .. -0.2', '-0.2 .. 0.0', '0.0 .. 0.2', '0.2 .. 0.4', '0.4 .. 0.6', '0.6 .. 0.8', '0.8 .. 1.0')


# plotting
plot_a  = tabla_plot_a %>%
    dplyr::group_by(`Pearson r`) %>%
    dplyr::mutate(mean_cophenetic_ds = mean(`Cophenetic dS`), sd_cophenetic_ds = sd(`Cophenetic dS`),
                  lower = mean_cophenetic_ds-sd_cophenetic_ds, upper = mean_cophenetic_ds+sd_cophenetic_ds) %>%
    dplyr::select(`Pearson r`, `mean_cophenetic_ds`, `sd_cophenetic_ds`, lower, upper) %>%
    unique() %>%
  ggplot(data = ., mapping = aes(x = `Pearson r`, y = mean_cophenetic_ds)) +
  geom_col() +
  ylab('Average cophenetic synonymous (dS) distance') +
  xlab("Pearson's r") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)

plot_b = tabla_plot_b %>%
  dplyr::group_by(`Spearman r`) %>%
  dplyr::mutate(mean_cophenetic_ds = mean(`Cophenetic dS`), sd_cophenetic_ds = sd(`Cophenetic dS`),
                lower = mean_cophenetic_ds-sd_cophenetic_ds, upper = mean_cophenetic_ds+sd_cophenetic_ds) %>%
  dplyr::select(`Spearman r`, `mean_cophenetic_ds`, `sd_cophenetic_ds`, lower, upper) %>%
  unique() %>%
  ggplot(data = ., mapping = aes(x = `Spearman r`, y = mean_cophenetic_ds)) +
  geom_col() +
  ylab('Average cophenetic synonymous (dS) distance') +
  xlab("Spearman's r") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)

# creating mixed plot and saving
grafico_final = cowplot::plot_grid(plot_a, plot_b,
                                   nrow = 1,
                                   labels = c('A', 'B'))

# lo exporto
grafico_final %>% ggsave(filename = '../results/supplementary_figure_7a.pdf', plot = ., device = 'pdf', width=12, height=6)


# getting manhattan distance and number of inparalogs per inparalog group (supp. figure 7b)
# getting number of inparalogs per group and creating a dictionary
index_inparalogs = readr::read_tsv('../../index_grupos_inparalogos/results/indexado_final_inparalogos.tibble.tsv', 
                                   col_names = TRUE) %>%
  tidyr::separate_rows(data = ., 'Genes', sep = ', ') %>%
  group_by(monophyletic_group_code) %>%
  dplyr::mutate(`Number of inparalogs` = n()) %>%
  dplyr::select(monophyletic_group_code, `Number of inparalogs`) %>%
  unique()

index_inparalogs_number.dict = index_inparalogs$`Number of inparalogs`
names(index_inparalogs_number.dict) = index_inparalogs$monophyletic_group_code

manhattan_distance.table = tabla_correlaciones %>%
  dplyr::select(`Inparalogs Group Code`, `Gen X`, `Gen Y`, `Manhattan distance`, `Type of Data`) %>%
  unique() %>%
  # filtering NA values
  dplyr::filter(!is.na(`Manhattan distance`)) %>%
  # counting number of inparalogs in each group
  dplyr::mutate(`Number of Inparalogs` = `Inparalogs Group Code` %>% tidytidbits::lookup_dbl(dict = index_inparalogs_number.dict)) %>%
  dplyr::mutate(`Number of Inparalogs` = case_when(`Number of Inparalogs` == 2 ~ 'Two genes',
                                                   `Number of Inparalogs` == 3 ~ 'Three genes',
                                                   `Number of Inparalogs` == 4 ~ 'Four genes',
                                                   `Number of Inparalogs` > 4 ~ '> Four genes')
                ) %>%
  group_by(`Inparalogs Group Code`) %>%
  dplyr::mutate(`Manhattan distance` = mean(`Manhattan distance`))
# converting variable of number of inparalogs into factors and ordering them
manhattan_distance.table$`Number of Inparalogs` %<>% as.factor()
manhattan_distance.table$`Number of Inparalogs` %<>% forcats::fct_relevel('Two genes', 'Three genes', 'Four genes', '> Four genes')

# plotting
plot_7b = manhattan_distance.table %>%
  dplyr::filter(str_detect(`Inparalogs Group Code`, 'SCM')) %>%
  dplyr::select(-c(`Gen X`, `Gen Y`)) %>%
  unique() %>%
  group_by(`Number of Inparalogs`) %>%
  dplyr::mutate(ymax = max(`Manhattan distance`), ymin = min(`Manhattan distance`), mean = mean(`Manhattan distance`)) %>%
  ungroup() %>%
  dplyr::arrange(`Number of Inparalogs`, `Manhattan distance`) %>%
  dplyr::mutate(`Order` = 1:nrow(.)) %>%
ggplot(data = ., mapping = aes(y = `Manhattan distance`, x = Order, fill = `Number of Inparalogs`)) +
  geom_col() +
  ylab('Manhattan distance') +
  xlab('') +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
plot_7b %>% ggsave(filename = '../../Graficos_dNdS/results/supplementary_figure_7_manhattan.pdf', plot = ., device = 'pdf', width=12, height=6)
