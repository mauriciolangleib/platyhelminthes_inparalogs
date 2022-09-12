# cargo funciones
source('functions.R')

label_enrichment_secr_platys = "Figure 1. Enrichment analysis for secreted products among inparalogs for the phylum Platyhelminthes. Fishers exact test was performed in order to assess if secreted products were enriched among inparalogs. Bars indicate the proportion of secreted (orange) and non-secreted (red) products for inparalogs (right panel) and other homologs (left panel). Genes belonging to singleton families were excluded from the analysis (see Materials and Methods). Colours indicate the class of the organism: Cestoda, blue; Trematoda, red; Monogenea, green; Turbellaria; orange.  FDR method was performed to correct for multiple hypothesis testing. (*; FDR < 0.05.)."

# levanto tabla correcta de codigos
TablaCodigos = readr::read_tsv('../../TablaCodigos.tsv', col_names = TRUE)

# Defino colores
        c(rep("darkorange",2),
          rep("olivedrab4",2),
          rep("coral3",12),
          rep("dodgerblue3", 14)) -> coloresEspecies.platy

armaPlotEnrichment(multigenes_genes_enlistados = multigenes_genes_enlistados, 
			TablaCodigos = TablaCodigos,
			coloresEspecies = coloresEspecies.platy,
			categoria = quo(secrecion), 
			cat_1 = quo(secretado), 
			cat_1_no = quo(no_secretado), 
			label_interes = "Secreted products", 
			label_plot = label_enrichment_secr_platys, 
			resultado_enriquecimiento = resultado_enriquecimiento.sec,
			color_2 = "darkorange", 
			color_1 = "indianred4") -> plot_plat_secr
plot_plat_secr %>%
ggsave(.,
	device = "pdf",
	width = 13,
	height = 12,
	filename = "../results/platys_secr_enrichment_widthmalo.pdf",
	dpi = 1000)

plot_plat_secr

resultado_enriquecimiento.sec %>%
	as_tibble() %>%
	merge(y = ., by.y = "organismo", x = TablaCodigos, by.x = "code") %>%
	arrange(., desc(class)) %>%
	dplyr::select(Species, class, p.values, FDR, significativo) %>%
	rename(Class = class, Significance = significativo) %>%
	knitr::kable(., caption = "Table Y. Epa...") %>%
	kable_styling(bootstrap_options = c("striped"))

#```

## Tabla de proteinas transmembrana
#```{r tabla_transmem, echo = FALSE, fig.width = 10, fig.height = 8}
# saco conteos de fisher para grafico, transmem

label_enrichment_transmem_platys = "Figure 2. Enrichment analysis for transmembrane products among inparalogs for the phylum Platyhelminthes. Fishers exact test was performed in order to assess if transmembrane products were enriched among inparalogs. Bars indicate the proportion of transmembrane (orange) and non-transmembrane (red) products for inparalogs (right panel) and other homologs (left panel). Genes belonging to singleton families were excluded from the analysis (see Materials and Methods). Colours indicate the class of the organism: Cestoda, blue; Trematoda, red; Monogenea, green; Turbellaria; orange.  FDR method was performed to correct for multiple hypothesis testing. (*; FDR < 0.05.)."

armaPlotEnrichment(multigenes_genes_enlistados = multigenes_genes_enlistados,
			TablaCodigos = TablaCodigos,
			coloresEspecies = coloresEspecies.platy,
                        categoria = quo(transmembrana),
                        cat_1 = quo(transmembrana),
                        cat_1_no = quo(no_transmembrana),
                        label_interes = "Transmembrane products",
                        label_plot = label_enrichment_transmem_platys,
			resultado_enriquecimiento = resultado_enriquecimiento.transmem,
                        color_1 = "deepskyblue4",
                        color_2 = "deepskyblue3") -> plot_plat.transmem
plot_plat.transmem  %>%
ggsave(.,
	device = "pdf",
	width = 13,
	height = 12,
	filename = "../results/platys_transmem_enrichment_widthmalo.pdf",
	dpi = 1000)

plot_plat.transmem

resultado_enriquecimiento.transmem %>%
	as_tibble() %>%
	merge(y = ., by.y = "organismo", x = TablaCodigos, by.x = "code") %>%
	arrange(., desc(class)) %>%
	dplyr::select(Species, class, p.values, FDR, significativo) %>%
	rename(Class = class, Significance = significativo) %>%
	knitr::kable(., caption = "Table Z. Epaa") %>%
	kable_styling(., bootstrap_options = c("striped"))

load("../data/analisisInparalogos.moluscos.RData")

TablaCodigos.setMoluscos = tibble(Species = c("Aplysia californica", "Crassostrea gigas",
						"Lottia gigantea", "Biomphalaria glabrata",
						"Octopus bimaculoides", "Crassostrea virginica",
						"Mizuhopecten yessoensis", "Pomacea canaliculata",
						"Helobdella robusta"),
					code = c("AAA","AAB","AAC","AAD","AAE","AAF","AAG", "AAJ", "AAM"),
					class = c(rep("Mollusca", 8), "Annelidae"))

# al cargar la imagen multigenes_genes_enlistados pasa a ser el de moluscos

# armo colores
c(rep("olivedrab4",1),
	rep("darkorange",8)) -> coloresEspecies.molusc

armaPlotEnrichment(multigenes_genes_enlistados = multigenes_genes_enlistados, 
			coloresEspecies = coloresEspecies.molusc,
			TablaCodigos = TablaCodigos.setMoluscos,
			categoria = quo(secrecion), 
			cat_1 = quo(secretado), 
			cat_1_no = quo(no_secretado), 
			label_interes = "Secreted products", 
			label_plot = label_enrichment_secr_platys, 
			resultado_enriquecimiento = resultado_enriquecimiento.sec.moluscos,
			color_2 = "darkorange", 
			color_1 = "indianred4") -> plot_plat_molusc
plot_plat_molusc

plot_plat_molusc %>%
	ggsave(.,
	device = "pdf",
	width = 13,
	height = 4,
	filename = "../results/molusc_secr_enrichment.pdf",
	dpi = 1000)

#```

## Analisis de enriquecimiento (proteinas transmembrana)
#```{r tabla_transmem_molusc, echo = FALSE, fig.width = 10, fig.height = 8}

armaPlotEnrichment(multigenes_genes_enlistados = multigenes_genes_enlistados,
			TablaCodigos = TablaCodigos.setMoluscos,
			coloresEspecies = coloresEspecies.molusc,
                        categoria = quo(transmembrana),
                        cat_1 = quo(transmembrana),
                        cat_1_no = quo(no_transmembrana),
                        label_interes = "Transmembrane products",
                        label_plot = "pa",
			resultado_enriquecimiento = resultado_enriquecimiento.transmem.moluscos,
                        color_1 = "deepskyblue4",
                        color_2 = "deepskyblue3") -> molusc.transmem
molusc.transmem

molusc.transmem %>%
	ggsave(.,
	device = "pdf",
	width = 13,
	height = 4,
	filename = "../results/mollusc_transmem_enrichment.pdf",
	dpi = 1000)
