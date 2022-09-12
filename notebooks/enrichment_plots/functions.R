armaPlotEnrichment = function (multigenes_genes_enlistados, TablaCodigos, coloresEspecies, categoria, cat_1, cat_1_no, label_interes, label_plot, color_1, color_2, resultado_enriquecimiento) {
  multigenes_genes_enlistados %>%
	group_by(organismo) %>%
	count(homologia, !! categoria) %>%
	tidyr::spread(., key = !! categoria, value = n) %>%
	mutate(p_secr = (!! cat_1)/( (!! cat_1) + (!! cat_1_no) ), p_non_secr = (1 - p_secr)) %>%
	tidyr::gather('p_secr', 'p_non_secr', key = "Cellular", value = "value") %>%
	dplyr::select(-c( !! cat_1_no , !! cat_1)) %>%
	rename(Species = organismo, variable = homologia) -> conteos_fisher

	str_replace_all(conteos_fisher$Cellular, "^p_secr$", label_interes) %>%
	str_replace_all(., "^p_non_secr$", "Other products") -> conteos_fisher$Cellular

	as.factor(conteos_fisher$Species) -> conteos_fisher$Species
	as.factor(conteos_fisher$variable) -> conteos_fisher$variable
	as.factor(conteos_fisher$Cellular) -> conteos_fisher$Cellular

	conteos_fisher %>%
		merge(y = ., by.y = "Species",
			x = TablaCodigos, by.x = "code") %>%
	dplyr::arrange(., desc(class)) %>%
	dplyr::select(Species, variable, Cellular, value, class) -> conteos_fisher


	conteos_fisher$Species <- factor(conteos_fisher$Species, levels = rev(TablaCodigos$Species))

	# Defino colores
#	c(rep("darkorange",2),
#	  rep("olivedrab4",2),
#	  rep("coral3",12),
#	  rep("dodgerblue3", 14)) -> coloresEspecies

	conteos_fisher %>% as_tibble() -> conteos_fisher

	library(forcats)

	conteos_fisher$variable %>%
		fct_recode(., Inparalogs = "inparalogo",
			      'Other homologs' = "homologo_no_inparalogo") -> conteos_fisher$variable

  label_enrichment <- paste0(strwrap(label_plot, 113), sep="", collapse="\n")
  
  resultado_enriquecimiento %>% 
  rename(code = "organismo") %>%
  left_join(TablaCodigos, by = "code") %>%
  dplyr::select(Species, significativo) %>%
  filter(significativo == "TRUE") %>%
  mutate(., variable = "Inparalogs", Cellular = "Significant", value = 1) -> annotacion
	
  annotacion$Species %<>% as.factor()
  annotacion$variable %<>% as.factor()
  annotacion$Cellular %<>% as.factor()
	
# ordeno las especies
conteos_fisher$Species %<>% as.factor()
conteos_fisher$Species %<>% forcats::fct_relevel(., "S. mediterranea", "M. lignano", 
	                                            "P. xenopodis", "G. salaris",
	                                            "S. rodhaini", "S. mattheei", "S. margrebowiei", "S. mansoni", "S. japonicum", "S. haematobium", 
	                                             "S. curassoni", "T. regenti", "F. hepatica", "E. caproni", "O. viverrini", "C. sinensis",
	                                            "S. erinaceieuropaei", "S. solidus", "D. latum", "H. taeniaeformis",
	                                             "T. saginata", "T. asiatica", "T. solium", "E. canadensis", "E. multilocularis", "E. granulosus",
	                                             "M. corti", "H. nana", "H. microstoma", "H. diminuta")
	
cols = c(color_1,color_2)
	ggplot(data = conteos_fisher, mapping = aes(x = Species, y = value, color = Cellular)) +
	  ggplot2::geom_col(aes(y=value, fill = Cellular), position = "stack", width = 0.75, color = rep("transparent",nrow(conteos_fisher)))  +
	  coord_flip(clip = 'off') +
	  #scale_x_discrete(expand = c(1,1)) +
	  scale_alpha(range = c(0.1, 1)) +
	  theme(panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.background = element_rect(fill = "gray95"),
	#        plot.margin = unit(c(2, .8, 2, .8), "cm"),
	        panel.spacing = unit(0.25, "lines"),
	        axis.line = element_line(colour = "white"),
	        axis.title.x=element_blank(),
	        axis.text.x=element_blank(),
	        axis.ticks.x=element_blank(),
	        axis.text.y = element_text(hjust = 1.0, face = c("italic"), colour = coloresEspecies,family = "serif", size = 14),
	        strip.text = element_text(family = "serif", size = 14, colour = "gray20",margin = margin(t = 5,r=5,b=5,l=5)),
	        #strip.background = element_rect(size = 5),
	        legend.text = element_text(family = "serif", size = 14),
	        axis.text.y.left = element_text(family = "serif", size = 14),
          plot.caption = element_text(vjust = 0, size = 14, hjust = 0),
	        legend.title = element_blank(),
	        axis.ticks.y = element_blank(),
	        axis.line.x = element_blank(),
	        axis.line.y = element_blank()) +
	scale_colour_manual(values = cols, aesthetics = c("colour", "fill")) +
	  facet_wrap(~variable) + labs(x="", y="") +
	  geom_text(data = annotacion, label = "*", 
			inherit.aes = TRUE, nudge_y = 0.05, colour = "black", size = 6,
			label.padding = unit(0.15, "lines")) %>% #+
#    labs(caption = label_enrichment) %>%

    return(.)
}
