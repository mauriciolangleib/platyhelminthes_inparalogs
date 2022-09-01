# plots involving dN and dS vs cophenetic distances, and related analyses

# loading tables regarding inparalogs and inparalog groups
read_tsv('../../data_inparalogos/multigenes_genes_enlistados.inparalogos.tsv') -> multigenes_genes_enlistados.inparalogos
read_tsv('../../data_inparalogos/indexados.tsv') -> indexados.tibble
#la vieja #read_tsv('../../../indexado_final_inparalogos.tibble.tsv') %>% 
#la vieja #  tidyr::separate_rows(., Genes, sep = ',') %>% 
#la vieja #  dplyr::mutate(., Genes = Genes %>% str_replace_all(., ' *', '')) -> index_monophyletic_groups.tibble

read_tsv('../../../index_grupos_inparalogos/results/indexado_final_inparalogos.tibble.tsv') %>% 
  tidyr::separate_rows(., Genes, sep = ',') %>% 
  dplyr::mutate(., Genes = Genes %>% str_replace_all(., ' *', '')) -> index_monophyletic_groups.tibble
  
# defining a fue useful functions
source('../R_scripts/functions_dNdS_corr.R')

# listing MLC files
dir_resultados = '../../../analisis_evolucion_molecular/resueltos'

dicc_indexados = indexados.tibble$WBPS
names(dicc_indexados) = indexados.tibble$codigo

index.tibble = readr::read_tsv('../data/interpro_index.tsv', col_names = TRUE)

index.tibble %<>% 
	# se descarta dato referente al codigo interno que se usa en este trabajo
	dplyr::select(-code) %>% unique() %>% 
	# se descartan filas para las que no se tienen datos
	dplyr::filter(!is.na(info)) %>% 
	# se agrupa segun gen
	group_by(Gene) %>% 
	# se colapsa la informacion que se tiene en una unica fila
	dplyr::mutate(info = paste(info, collapse = '\n')) %>%
	# se desagrupa
	ungroup() %>%
	# saco las filas repetidas, quedandome solo con una
	unique()

list.files(dir_resultados, pattern = "^out$", full.names = T, recursive = T) %>%
tibble(archivo_mlc = .) %>%
dplyr::filter(str_detect(archivo_mlc, 'codeml_dir') & str_detect(archivo_mlc, 'M1')) %>%
dplyr::mutate(grupo_inparalogos = archivo_mlc %>% str_replace_all(., ".codeml_dir.*$", ""),
          		organismo = archivo_mlc %>% str_split(., "/") %>% map_chr(6))  -> mlc1s.tibble

list.files(dir_resultados, pattern = "^out$", full.names = T, recursive = T) %>%
tibble(archivo_mlc = .) %>%
dplyr::filter(str_detect(archivo_mlc, 'codeml_dir') & str_detect(archivo_mlc, 'M2')) %>%
dplyr::mutate(grupo_inparalogos = archivo_mlc %>% str_replace_all(., ".codeml_dir.*$", ""),
          		organismo = archivo_mlc %>% str_split(., "/") %>% map_chr(6))  -> mlc2s.tibble

    
merge(x = mlc1s.tibble, by.x = "grupo_inparalogos",
	y = mlc2s.tibble, by.y = "grupo_inparalogos") %>%
	as_tibble() %>%
	dplyr::select(-organismo.y) %>%
	dplyr::rename(archivo_mlc1 = archivo_mlc.x,
			archivo_mlc2 = archivo_mlc.y,
			organismo = organismo.x) -> mlcs_corr.tibble
      
mlcs_corr.tibble %>%
	mutate(., grados_de_libertad = pmap_dbl(list(archivo_mlc2,archivo_mlc1), ~ try(lector_np(..1) - lector_np(..2))),
		  lnL_1 = pmap_dbl(list(archivo_mlc1), ~try(lector_lnL(..1))),
		  lnL_2 = pmap_dbl(list(archivo_mlc2), ~try(lector_lnL(..1))),
		  sitios_evolucion_BEB = (pmap(list(archivo_mlc2), ~ parser_sitios_BEB(..1)) ),
      sitios_evolucion_BEB_totales = pmap_dbl(list(sitios_evolucion_BEB), ~ {filter(.x, sitio == 1) %>% nrow()}) ,
      sitios_evolucion_BEB_AltaProb = pmap_dbl(list(sitios_evolucion_BEB), ~ {filter(.x, sitiosAltaProb == 1) %>% nrow()}) ,
		  N_mlc2 = (pmap_dbl(list(archivo_mlc2), ~{parser_dNs(..1) %>% .$`dN` %>% .[1] }) ),
      S_mlc2 = (pmap_dbl(list(archivo_mlc2), ~{parser_dNs(..1) %>% .$`dS` %>% .[1] }) ),
      dN_dS_mlc2 = (pmap_dbl(list(archivo_mlc2), ~{parser_dNs(..1) %>% .$`dN/dS` %>% .[1] }) )
		) %>%
    dplyr::select(-sitios_evolucion_BEB) -> datos_codeml.tibble

# saving result
if (!dir.exists('../results/molecular_evolution_correlation_plots')) {
  system('mkdir ../results/molecular_evolution_correlation_plots')
}

write_tsv(datos_codeml.tibble, '../results/molecular_evolution_correlation_plots/datos_codeml.tibble.tsv')

purrr::pmap_dfr(list(datos_codeml.tibble$lnL_1, datos_codeml.tibble$lnL_2, datos_codeml.tibble$grados_de_libertad), ~LRT(lnL_1 = ..1, lnL_2 = ..2, grados_libertad =..3)) -> resultados_LRT.tibble
write_tsv(resultados_LRT.tibble, '../results/molecular_evolution_correlation_plots/resultados_LRT.tibble.tsv')

bind_cols(datos_codeml.tibble, resultados_LRT.tibble) -> tabla_final.tibble
write_tsv(tabla_final.tibble, '../results/molecular_evolution_correlation_plots/tabla_final.tibble.tsv')

# loading rna-seq data and changing some colnames in order to be coherent with defined functions
rnaquant_finished.tibble = readr::read_tsv('../results/expression_tables/rnaquant_tpms.tsv', col_names = TRUE) %>%
	dplyr::filter(`Multimapping reads` == 'Not allowed') %>%
	dplyr::rename(gene_id = 'GeneID',
		      dev_stage = 'Condition',
		      Study = 'sample', 
		      tpm_median = 'TPM median') %>%
	dplyr::select(-c(TPM)) %>%
	dplyr::mutate(species = case_when(species == 'hmicrostoma' ~ 'hym', species == 'smansoni' ~ 'scm'))

# filtro para considerar los datos de S. mansoni e H. microstoma
tabla_final.tibble %>% 
  dplyr::filter(str_detect(string = organismo, pattern = 'scm|hym')) %>%
# tomo en cuenta un analisis por grupo monofiletico
  group_split(grupo_inparalogos) %>%
  purrr::map(., ~{
    grupo = .x
    
    # cargo las secuencias del grupo monofiletico
    inparalogos_grupo = ape::read.tree(grupo$grupo_inparalogos) %>% .$tip.label
    inparalogos_grupo.WBPS = indexados.tibble %>% filter(codigo %in% inparalogos_grupo) %>% .$WBPS
    nombre_grupo = grupo$grupo_inparalogos %>% str_split(., '/') %>% unlist() %>% .[length(.)]
    org = grupo$organismo
    
    # cargo el arbol de este grupo de inparalogos y calculo la distancia cofenetica de dS
    list.files(glue('../../../analisis_evolucion_molecular/resueltos/{org}/{nombre_grupo}.codeml_dir'), pattern = 'out', full.names = T, recursive = T) %>% 
    .[str_detect(string = ., pattern='M0~')] %>%
    treeio::read.codeml_mlc(.) -> data_codeml_mlc
    
    cophenetic_dS = calculador_cophenetic_dS(data_codeml_mlc)
    
    # saco los datos de RNA-Seq y veo la correlacion    
    
    # opero sobre los datos
  ## en caso de que el organismo sea H. microstoma
    if (unique(grupo$organismo) == 'hym') {
      resultados_RNASeq.TPMs.filtrado.hymenolepis.tsv = rnaquant_finished.tibble %>% dplyr::filter(species == 'hym') %>% dplyr::select(-species)
      # arreglo los factores 
      resultados_RNASeq.TPMs.filtrado.hymenolepis.tsv$dev_stage %<>% forcats::fct_recode(., `Scolex-Neck` = 'Adult, scolex',
                                                     `Mid-section` = 'Adult, mid-section',
                                                     `Posterior` = 'Adult, posterior',
						     `Adult (whole)` = 'Adult, whole',
                                                     `Larvae (mid-metamorphosis)` = 'mid-metamorphosis larvae')  
  
  resultados_RNASeq.TPMs.filtrado.hymenolepis.tsv$dev_stage %<>% forcats::fct_relevel(., 'Scolex-Neck', 'Mid-section', 'Posterior', 'Larvae (mid-metamorphosis)', 'Adult (whole)')
    
      resultados_RNASeq.TPMs.org =  resultados_RNASeq.TPMs.filtrado.hymenolepis.tsv
# clasifico en base al tipo de datos
      resultados_RNASeq.TPMs.org %<>% 
      	dplyr::mutate(tipo_datos = case_when(dev_stage == 'Mid-section' | dev_stage == 'Posterior' | dev_stage == 'Scolex-Neck' ~ 'spatial',
					     dev_stage == 'Adult (whole)' | dev_stage == 'Larvae (mid-metamorphosis)' ~ 'developmental'))
					     
     # corro el analisis sobre el set 'spatial'
     resultados_RNASeq.TPMs.org %>%
     	dplyr::filter(., tipo_datos == 'spatial') %>%
	dplyr::select(-c(tipo_datos, `Multimapping reads`)) %>%
	analisis_expresion(rna_seq_data=., cophenetic_dS=cophenetic_dS, inparalogos_grupo.WBPS = inparalogos_grupo.WBPS, org = org) -> resultado_spatial
	
     # corro el analisis sobre el set 'developmental'
     resultados_RNASeq.TPMs.org %>%
     	dplyr::filter(., tipo_datos == 'developmental') %>%
	dplyr::select(-tipo_datos) %>%
	analisis_expresion(rna_seq_data=., cophenetic_dS=cophenetic_dS, inparalogos_grupo.WBPS = inparalogos_grupo.WBPS, org = org) -> resultado_developmental
     # si existe tabla en estos resultados les pongo un tag
     if(!is.null(resultado_spatial)){
     	resultado_spatial %<>%
		dplyr::mutate(tipo_datos = 'spatial')
     }
     
     if(!is.null(resultado_developmental)){
     	resultado_developmental %<>%
		dplyr::mutate(tipo_datos = 'developmental')
     }
     
     # uno los resultados
     rbind(resultado_spatial, resultado_developmental) -> resultado_combinado
     # lo imprimo a pantalla
     return(resultado_combinado)
    }
    
## en caso de que el organismo sea S. mansoni
    if (unique(grupo$organismo) == 'scm') {
      resultados_RNASeq.TPMs.filtrado.tsv = rnaquant_finished.tibble %>% dplyr::filter(species == 'scm') %>% dplyr::select(-species)
      # arreglo los factores de las etapas de desarrollo
        resultados_RNASeq.TPMs.filtrado.tsv$dev_stage %<>% forcats::fct_recode(., Miracidia = 'miracidia',
					`First stage sporocyst` = 'first stage sporocyst',				       
                                        `Second stage sporocyst` = 'second stage sporocyst',
                                         Cercariae = 'Cercariae',
                                         `Schistosomulae (3 hr)` = '3 hour schistosomule',
                                         `Schistosomulae (24 hr)` = '24 hour schistosomule',
                                     `Juvenile` = 'juvenile males and females NMRI',
                                     `Adult` = 'adult males and females NMRI')  
  
  resultados_RNASeq.TPMs.filtrado.tsv$dev_stage %<>% forcats::fct_relevel(., 'Miracidia', 'First stage sporocyst', 'Second stage sporocyst', 'Cercariae',
                                                      'Schistosomulae (3 hr)', 'Schistosomulae (24 hr)',
                                                      'Juvenile', 'Adult')
						      
      resultados_RNASeq.TPMs.org = resultados_RNASeq.TPMs.filtrado.tsv
      # corro el analisis      
      analisis_expresion(rna_seq_data=resultados_RNASeq.TPMs.org, cophenetic_dS=cophenetic_dS, inparalogos_grupo.WBPS = inparalogos_grupo.WBPS, org = org) -> resultado_mansoni
      
      # si existe, le pongo un tag a los datos
      if(!is.null(resultado_mansoni)){
     	resultado_mansoni %<>%
		dplyr::mutate(tipo_datos = 'developmental')
     }
     
     return(resultado_mansoni)
    }
    
    
  }) %>%
  # limpio la tabla
  purrr::reduce(.,bind_rows) %>%
  dplyr::filter(!is.na(x) | !is.na(y)) -> correlaciones_spearman_cofeneticdS

# saving the table
# creating dictionary gene/inparalog group
read_tsv('../../../indexado_final_inparalogos.tibble.tsv') %>% 
  tidyr::separate_rows(., Genes, sep = ',') %>% 
  dplyr::mutate(., Genes = Genes %>% str_replace_all(., ' *', '')) -> index_monophyletic_groups.tibble

inparalog_groups.dict = index_monophyletic_groups.tibble$monophyletic_group_code
names(inparalog_groups.dict) = index_monophyletic_groups.tibble$Genes

correlaciones_spearman_cofeneticdS %>%
	dplyr::mutate(`Inparalogs Group Code` = x %>% tidytidbits::lookup_chr(dict = inparalog_groups.dict, default = NA),
		      inpagroup.y = y %>% tidytidbits::lookup_chr(dict = inparalog_groups.dict, default = NA),
		      organismo = case_when(organismo == 'hym' ~ 'H. microstoma', organismo == 'scm' ~ 'S.mansoni'),
		      correlacion = correlacion %>% str_to_title(),
		      tipo_datos = tipo_datos %>% str_to_title()) -> correlaciones_spearman_cofeneticdS.tosave

# testing that each 2uple belongs to the sample inparalog group -> it is ok
correlaciones_spearman_cofeneticdS.tosave %>%
	dplyr::filter(`Inparalogs Group Code` != inpagroup.y)

# saving the table
correlaciones_spearman_cofeneticdS.tosave %<>%
	dplyr::left_join(x = ., y = (index.tibble %>% dplyr::select(Gene, info)), by = c('x' = 'Gene')) %>%
	dplyr::rename(info_x = 'info') %>%
	dplyr::left_join(x = ., y = (index.tibble %>% dplyr::select(Gene, info)), by = c('y' = 'Gene')) %>%
	dplyr::rename(info_y = 'info') %>%
	dplyr::rename(`Gen X` = 'x', `Gen Y` = 'y', `Organism` = 'organismo', `Cophenetic dS` = 'cophenetic dS', `Type of Data` = 'tipo_datos',
		      `InterPro Annotation Gen X` = 'info_x', `InterPro Annotation Gen Y` = 'info_y') %>%
	dplyr::select(-inpagroup.y) %>%
	tidyr::pivot_wider(data = ., id_cols = c(`Gen X`, `Gen Y`, `Organism`, `Cophenetic dS`, `Type of Data`, `Inparalogs Group Code`, `InterPro Annotation Gen X`, `InterPro Annotation Gen Y`), 
			    names_from = `correlacion`, values_from = `r`) %>%
	dplyr::rename(`Pearson r` = 'Pearson', `Spearman r` = 'Spearman', `Manhattan distance` = 'Manhattan') %>%
	dplyr::relocate(., `Inparalogs Group Code`, `Gen X`, `Gen Y`, `Organism`, `Pearson r`, `Spearman r`, `Manhattan distance`, `Cophenetic dS`, `Type of Data`, 
		    `InterPro Annotation Gen X`, `InterPro Annotation Gen Y`)

correlaciones_spearman_cofeneticdS.tosave %>% readr::write_tsv(., '../results/molecular_evolution_correlation_plots/data_supplementary_table_5.tsv')
