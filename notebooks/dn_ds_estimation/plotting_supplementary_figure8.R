# creating supplementary figure 8
library(ggridges)
library(cowplot)
library(tidyverse)
library(glue)
library(magrittr)
library(tidytidbits)

# loading dN/dS table
tabla_dnds = readr::read_tsv('../../analisis_evolucion_molecular.emprolijando/results/tabla_final.formato_xlsx.tsv') %>% dplyr::filter(dS > 0.1)

# loading taxa table
tabla_codigos = readr::read_tsv('../../TablaCodigos.tsv')

# getting total number of inparalogous group in each species
inparalogous_groups_per_species = readr::read_tsv('../../index_grupos_inparalogos/results/indexado_final_inparalogos.tibble.tsv', col_names = T) %>% group_by(Species) %>% summarise(count = n())
## creating dictionary for this values
inparalogous_groups_per_species.dict = inparalogous_groups_per_species$count
names(inparalogous_groups_per_species.dict) = inparalogous_groups_per_species$Species

# organizing data according to its class
## creating dictionaries
species.dict = tabla_codigos$Species
names(species.dict) = tabla_codigos$code

classes.dict = tabla_codigos$class
names(classes.dict) = tabla_codigos$code

## mutating to add class, modifiy name of species 
## ordering species to be coherent with the manuscript

#tabla_dnds %<>%
#  dplyr::mutate(organismo = group_inparalogos %>% stringr::str_split('_') %>% purrr::map_chr(2, ~.x) %>% stringr::str_to_lower()) %>%
#  tidyr::unnest(cols = c(organismo)) %>%
#  dplyr::mutate(organismo = case_when(organismo %in% c('ecg', 'fsc', 'msl') ~ stringr::str_to_title(organismo), TRUE ~ organismo)) %>%
#  dplyr::mutate(Class = organismo %>% tidytidbits::lookup_chr(., dict = classes.dict) %>% as.factor(),
#                organismo = organismo %>% tidytidbits::lookup_chr(., dict = species.dict) %>% as.factor()) %>%
#  dplyr::rename(Species = 'organismo')

tabla_dnds$Species %<>% forcats::fct_relevel(., "S. mediterranea", "M. lignano", 
                                             "P. xenopodis", "G. salaris",
                                             "S. rodhaini", "S. mattheei", "S. margrebowiei", "S. mansoni", "S. japonicum", "S. haematobium", 
                                             "S. curassoni", "T. regenti", "F. hepatica", "E. caproni", "O. viverrini", "C. sinensis",
                                             "S. erinaceieuropaei", "S. solidus", "D. latum", "H. taeniaeformis",
                                             "T. saginata", "T. asiatica", "T. solium", "E. canadensis", "E. multilocularis", "E. granulosus",
                                             "M. corti", "H. nana", "H. microstoma", "H. diminuta")

# plotting with ggridges
c(rep("darkorange",2),
  rep("olivedrab4",2),
  rep("coral3",12),
  rep("dodgerblue3", 14)) -> coloresEspecies

final_plot = tabla_dnds %>%
  # binning data
  #dplyr::mutate(dN_dS_mlc2 = case_when(dN_dS_mlc2 == 0 ~ '0 - 0.3', 
  #                                     between(dN_dS_mlc2, 0, 0.3) ~ '0 - 0.3', 
  #                                     between(dN_dS_mlc2, 0.3, 0.6) ~ '0.3 - 0.6', 
  #                                     between(dN_dS_mlc2, 0.6, 0.9) ~ '0.6 - 0.9' , 
  #                                     between(dN_dS_mlc2, 0-9, 1.1) ~ '0.9 - 1.1' , 
  #                                     dN_dS_mlc2 >= 1.1 ~ '>1.1')) %>% 
  dplyr::mutate(`Global dN/dS` = case_when(between(`Global dN/dS`, 0, 1.0) ~ `Global dN/dS`, 
                                       `Global dN/dS` > 1.0 ~ 1.1)) %>% 
  ggplot(data = ., mapping = aes(x = `Global dN/dS`)) +
  #geom_col() +
  #xlim(0, 1.100001) +
  #geom_tile() +
  #ggridges::geom_density_ridges(mapping = aes(x = dN_dS_mlc2, y = Species), 
  #                              stat = 'density', scale = 0.9) +
  theme(axis.text.y = element_text(hjust = 1.0, face = c("italic"), colour = coloresEspecies)) +
  ggridges::geom_density_ridges(mapping = aes(x = `Global dN/dS`, y = Species), 
                                stat = 'binline', scale = 0.75, binwidth = 0.25) +
  #ggplot2::scale_x_binned(breaks = c(0, 0.3)) +
  #geom_histogram(breaks = c(0, 0.3, 0.6, 0.9, 1.1), aes(y = stat(density))) +
  labs(title = "Supplementary Figure 8: dN/dS distribution among species") +
  xlab('dN/dS') +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c('0 - 0.25', '0.25 - 0.50', '0.50 - 0.75', '0.75 - 1.0', '> 1.0'))
  #facet_wrap(~Class+Species, dir = 'v')

  # saving plot to results folder
  ggplot2::ggsave(filename = '../results/supplementary_figure8.pdf', final_plot, height = 12, width = 12)

# making Supplementary Figure 8a. Number of inparalogous groups evolving with codon positions with dN/dS > 1
## re-loading dN/dS table
tabla_dnds = readr::read_tsv('../../analisis_evolucion_molecular.emprolijando/results/tabla_final.formato_xlsx.tsv') 

tabla_dnds$Species %<>% forcats::fct_relevel(., "S. mediterranea", "M. lignano", 
                                             "P. xenopodis", "G. salaris",
                                             "S. rodhaini", "S. mattheei", "S. margrebowiei", "S. mansoni", "S. japonicum", "S. haematobium", 
                                             "S. curassoni", "T. regenti", "F. hepatica", "E. caproni", "O. viverrini", "C. sinensis",
                                             "S. erinaceieuropaei", "S. solidus", "D. latum", "H. taeniaeformis",
                                             "T. saginata", "T. asiatica", "T. solium", "E. canadensis", "E. multilocularis", "E. granulosus",
                                             "M. corti", "H. nana", "H. microstoma", "H. diminuta")

plot_8a1.tibble = tabla_dnds %>% 
  dplyr::mutate(Type = case_when((dS <= 0.1) ~ 'dS <= 0.1',
				 (`Significative (LRT > 6)` == '***') & (dS > 0.1) ~ 'LRT > 6',
				 is.na(`Significative (LRT > 6)`) ~ 'Missing data',
				 (`Significative (LRT > 6)` == '-') & (dS > 0.1) ~ 'LRT < 6',
				)
	       ) %>%
  group_by(Species, Type) %>%
  summarise(count = n()) %>% 
  #dplyr::filter(`Significative (LRT > 6)` == '***' & hay_sitio_BEB == 1) %>%
  ungroup() %>%
  #tidyr::complete(Species) %>%
  #dplyr::mutate(`Significative (LRT > 6)` = `Significative (LRT > 6)` %>% replace_na(0),
  #	        count = count %>% replace_na(0)) %>%
  dplyr::mutate(total_count = Species %>% tidytidbits::lookup_dbl(., dict = inparalogous_groups_per_species.dict),
	       	 `Percentage of inparalogs groups` = count*100/total_count) 

plot_8a1.tibble$`Type` %<>% as.factor()
plot_8a1.tibble$`Type` %<>% forcats::fct_relevel(., 'Missing data', 'dS <= 0.1', 'LRT < 6', 'LRT > 6')

plot_8a1 = plot_8a1.tibble %>%
    ggplot2::ggplot(data = ., mapping = aes(x = reorder(Species, desc(Species)), y = `Percentage of inparalogs groups`, fill = Type)) +
    theme(axis.text.x = element_text(hjust = 1.0, face = c("italic"), colour = rev(coloresEspecies), angle = 45),
          plot.title = element_text(hjust = 0.5)) +
    #theme(axis.text.x = element_blank(),
    #      plot.title = element_text(hjust = 0.5)) +
    geom_col() +
    xlab('') +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    ylab('Groups of inparalogs (%)') +
    ggplot2::ggtitle('Distribution of dN/dS calculations for inparalogous groups')

plot_8a2 = tabla_dnds %>%
	dplyr::filter(dS > 0.1) %>%
	dplyr::select(Species, `# DB`, `# sites (BEB pp > 0.95)`, `Significative (LRT > 6)`) %>%
	unique() %>%
	dplyr::filter(`Significative (LRT > 6)` == '***') %>%
	ungroup() %>%
	tidyr::complete(Species) %>%
	ggplot(data = ., mapping = aes(x = Species, y = `# sites (BEB pp > 0.95)`)) +
	#geom_violin() +
	geom_boxplot() +
	theme(axis.text.y = element_text(hjust = 1.0, face = c("italic"), colour = coloresEspecies),
          plot.title = element_text(hjust = 0.5)) +
	ylab('Number of codons evolving at dN/dS > 1') +
	xlab('Species') +
	coord_flip() #+ 
	#scale_y_log10()

plot_8a = cowplot::plot_grid(plot_8a1, plot_8a2,
                                   nrow = 2,
                           labels = c('', ''))
ggplot2::ggsave(filename = '../results/supplementary_figure8a.pdf', plot_8a, height = 10, width = 10)
