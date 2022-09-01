ordenador_tablas = function(tabla) {
             # reordeno los genes, para tener criterio coherente con la tabla de dS cofenetico
             for (i in seq_along(tabla$x)) {
                if (tabla$x[i] < tabla$y[i]) {
                tabla$y[i] -> y.tmp
                tabla$x[i] -> x.tmp
                tabla$y[i] = x.tmp
                tabla$x[i] = y.tmp
                    }
              }
           return(tabla)
     }
     
calculador_cophenetic_dS = function(data_codeml_mlc) {
    # obtengo arbol de este archivo mlc
        data_codeml_mlc %>% treeio::get.tree() -> arbol_codeml
        # obtengo distancia cofenetica, tomando como largo de rama el dS calculado por PAML
        arbol_codeml %>%
          treeio::as.treedata() %>%
            as_tibble() %>%
            dplyr::full_join(x = ., y = (data_codeml_mlc %>% get.data()), by = c('node')) %>%
            dplyr::select(parent, node, dS, label) %>%
            dplyr::rename(branch.length = 'dS') %>%
            ape::as.phylo() %>%
            ape::cophenetic.phylo(x = .) -> cophenetic_dS
        # modifico un poco la tabla de distancia cofeneticas
        cophenetic_dS %<>%
        corrr::as_cordf() %>%
        corrr::stretch(x = ., remove.dups=T, na.rm=T) %>%
        dplyr::rename(`cophenetic dS` = 'r')
	
	cophenetic_dS %<>%
		dplyr::mutate(x = tidytidbits::lookup_chr(dict = dicc_indexados, x, default = NA),
			      y = tidytidbits::lookup_chr(dict = dicc_indexados, y, default = NA)) %>%
		ordenador_tablas(tabla = .)

        # devuelvo el resultado
        return(cophenetic_dS)
}

calculador_cophenetic_dN = function(data_codeml_mlc) {
    # obtengo arbol de este archivo mlc
        data_codeml_mlc %>% treeio::get.tree() -> arbol_codeml
        # obtengo distancia cofenetica, tomando como largo de rama el dS calculado por PAML
        arbol_codeml %>%
          treeio::as.treedata() %>%
            as_tibble() %>%
            dplyr::full_join(x = ., y = (data_codeml_mlc %>% get.data()), by = c('node')) %>%
            dplyr::select(parent, node, dN, label) %>%
            dplyr::rename(branch.length = 'dN') %>%
            ape::as.phylo() %>%
            ape::cophenetic.phylo(x = .) -> cophenetic_dN
        # modifico un poco la tabla de distancia cofeneticas
        cophenetic_dN %<>%
        corrr::as_cordf() %>%
        corrr::stretch(x = ., remove.dups=T, na.rm=T) %>%
        dplyr::rename(`cophenetic dN` = 'r')
        
	cophenetic_dN %<>%
		dplyr::mutate(x = tidytidbits::lookup_chr(dict = dicc_indexados, x, default = NA),
			      y = tidytidbits::lookup_chr(dict = dicc_indexados, y, default = NA)) %>%
		ordenador_tablas(tabla = .)
	
        # devuelvo el resultado
        return(cophenetic_dN)
}

parserW = function(mlc){
	system(glue("grep 'w:' {mlc}"), intern = TRUE) %>%
		str_split(" ") %>%
		unlist() %>%
		parse_character() %>%
		.[!is.na(.)] %>%
	return(.)
}

parser_dNs = function(mlc){
	paste("sed -n '/dN & dS for each branch/,/Naive Empirical Bayes.*$/p' ",
			mlc,
			" | grep -v 'Naive' | grep -v '^$' | grep -v 'each branch' | grep -v 'Time'",
			sep = "") %>%
	system(., intern = TRUE) %>%
	read_delim(delim = ' ',
				col_types = cols(col_character(), col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number()),
				col_names = TRUE, trim_ws = TRUE) %>%
	return(.)
}

parser_sitios_BEB = function(mlc){
	paste("sed -n '/Bayes Empirical Bayes (BEB) analysis/,/The grid (see ternary graph for p0-p1).*$/p' ",
                        mlc,
                        " | grep -v 'Bayes' | grep -v '^$' | grep -v 'Positively' | grep -v 'amino' | grep -v 'The grid'",
                        sep = "") %>%
	system(., intern = TRUE) %>%
	.[-1] %>%
	str_squish() %>%
	str_replace_all(., "\\*","") -> tmp
	if (length(tmp) >= 1) {
	tmp %>%
	#read_delim(., ' ', col_types = cols(col_double(), col_character(),
	#				col_double(), col_character(), col_character(), col_character()),
	#			col_names = c('Site', 'Amino acid', 'Pr(w>1)', 'post', 'masmenos', 'tmp'),
	#		skip = 1) %>%
	str_split(., ' ') %>%
	purrr::map_dfr(., ~{
				enframe(.x) %>%
				dplyr::select(-name) -> tmp
				tibble(`Site` = tmp$value[1],
					`Amino acid` = tmp$value[2],
					`Pr(w>1)` = tmp$value[3],
					`post` = tmp$value[4],
					`masmenos` = tmp$value[5],
					`tmp` = tmp$value[6])
				}
			) %>%
	mutate(post_mean = glue('{post}+-{tmp}'),
			sitio = 1,
			sitiosAltaProb = case_when(`Pr(w>1)` >= 0.95 ~ 1,
							`Pr(w>1)` < 0.95 ~ 0)
		) %>%
	dplyr::select(-c(post,masmenos,tmp)) -> resultado
	}
	if (length(tmp) == 0){
		tibble("Site" = NA, "Amino acid" = NA, "Pr(w>1)" = NA, "post_mean" = NA, "sitio" = 0, "sitiosAltaProb" = 0) -> resultado
	}
	return(resultado)
}

# Defino funcion que lee valor de lnL
lector_lnL = function (archivo_mlc){
	glue("grep 'lnL' {archivo_mlc}") %>%
	system(., intern = TRUE) %>%
		str_split(., pattern = "\\): ") %>%
		unlist(use.names = F) %>%
		.[2] %>%
		str_replace(string = ., pattern = "^ *", replacement="") %>%
		str_split(., pattern = " ") %>%
		unlist(use.names = F) %>%
		.[1] %>%
		as.numeric(.) %>%
	return(.)
}
# Defino funcion que lee el numero de parametros empleado
lector_np = function (archivo_mlc) {
	glue("grep 'lnL' {archivo_mlc}") %>%
		system(., intern = TRUE) %>%
		strsplit(., split = "):") %>%
		unlist() %>%
		.[1] %>%
		strsplit(., split = "np:") %>%
		unlist() %>%
		.[2] %>%
		as.numeric() %>%
	return(.)
}

# Defino funcion que haga LRT, que defina si el valor es significativo y devuelva dicho valor
LRT = function(lnL_1, lnL_2, grados_libertad){
	# Obtenemos valor de 2*delta(lnL_2 - lnL_1)
	LR = 2*(lnL_2 - lnL_1)
	# Obtenemos valor de Chi2 al grado de libertad dado
	chi = qchisq(p = .95, df = grados_libertad)
	# Se evalua si el valor es significativo
	tibble(Likelihood_ratio = LR,
		Grados_de_libetad = grados_libertad,
		chi_2 = chi) %>%
		mutate(Significativo = case_when((LR > chi) ~ 1,
						 (LR < chi) ~ 0)) %>%
	return(.)
}

completa_genes_expresion = function(datos_expresion, inparalogos_grupo.WBPS) {
	faltantes = setdiff(inparalogos_grupo.WBPS, unique(datos_expresion$gene_id))
		    if (length(faltantes > 0)) { 
			    datos_expresion %<>%
		            	tibble::add_row(., gene_id = faltantes) 
	            } 
	# devuelvo la tabla, haya sido modificada o no
	return(datos_expresion)
}


## calculador para correlaciones de Spearman/Pearson
calculador_correlaciones = function(datos_expresion.corr, metodo) {
	datos_expresion.corr %>%
	     	 # llevo a un formato amigable para la correlacion con el paquete corrr
	     	 tidyr::pivot_longer(data = ., cols = (colnames(.) %>% .[-c(1)]), names_to = 'dev_stage', values_to = 'tpm_median') %>%
	         tidyr::pivot_wider(data = ., id_cols = c('gene_id', 'dev_stage'), names_from = 'gene_id', values_from = 'tpm_median') %>%
	         # saco las columnas categoricas
	      	 dplyr::select(-c(dev_stage)) %>%
	         # llevo la etapa de desarrollo nombre de fila
	         #column_to_rownames(., var = 'dev_stage') %>%
	         # saco correlacion de Spearman
	         corrr::correlate(., method = metodo, use = 'everything') %>%
	         corrr::shave(.) %>%
	         corrr::stretch(., na.rm=T, remove.dups=T) -> correlacion
		 
		 # devuelvo al usuario
		 return(correlacion)
}

## funcion para calcular distancia de Manhattan u otras
calcula_distancias = function(datos_expresion.paraDistancias, metodo) {
	datos_expresion.paraDistancias %<>%
		 	tidyr::pivot_longer(data = ., cols = (colnames(.) %>% .[-c(1)]), names_to = 'dev_stage', values_to = 'tpm_median') %>%
		 	tidyr::pivot_wider(data = ., id_cols = c('gene_id', 'dev_stage'), names_from = 'dev_stage', values_from = 'tpm_median') 
	
		# armo un diccionario
		diccionario.genes = datos_expresion.paraDistancias$gene_id
		names(diccionario.genes) = 1:length(datos_expresion.paraDistancias$gene_id)
	
		datos_expresion.paraDistancias %>%
			# paso los genes a nombre de fila
			column_to_rownames(., var = 'gene_id') %>%
			# calculo la correlacion por metodo establecido
			stats::dist(x = ., method = metodo) %>%
			broom::tidy(x = .) -> correlacion.distancia
		# modifico un poco	
		correlacion.distancia %<>% 
			dplyr::rename(x = 'item1', y = 'item2') %>%
			dplyr::rename(r = 'distance') %>%
			dplyr::filter(!is.na(r))
			
		correlacion.distancia$x %<>% as.character()
		correlacion.distancia$y %<>% as.character()
		# devuelvo la tabla
		return(correlacion.distancia)
}

## defino funcion de analisis de RNA-Seq con datos de cofenetica de dS
analisis_expresion = function(rna_seq_data, cophenetic_dS, inparalogos_grupo.WBPS, org){ 
    resultados_RNASeq.TPMs.org = rna_seq_data 
    resultados_RNASeq.TPMs.org %>% 
	        # filtro datos de RNA-Seq para quedarme con los inparalogos
	        dplyr::filter(gene_id %in% inparalogos_grupo.WBPS) %>%
		# ordeno en base a factores
		dplyr::arrange(dev_stage) %>%
		# saco los datos de estudio
		dplyr::select(-Study) %>%
		unique() %>%
	        # pivoteo la tabla a formato ancho
	        tidyr::pivot_wider(data = ., id_cols = c('dev_stage', 'gene_id'), names_from = 'dev_stage', values_from = 'tpm_median') -> datos_expresion
    
    # opero si alguno de los genes tiene datos de expresion
    if(nrow(datos_expresion) > 0){ 
	# agrego los gene_id que no estan ahi
	datos_expresion %<>% completa_genes_expresion(., inparalogos_grupo.WBPS)    
      	
	# SACO CORRELACION DE SPEARMAN
	datos_expresion %>%
	 # modifico valores de expresion menores a uno, reemplazandolos por cero
	 dplyr::mutate_at(., colnames(.)[-c(1)], ~if_else(. < 1, 0, .)) %>%
	 dplyr::mutate_at(., colnames(.)[-c(1)], ~if_else(is.na(.), 0, .)) -> datos_expresion.corr
	 
	 correlacion.spearman = calculador_correlaciones(datos_expresion.corr, metodo = 'spearman')
      	 
	 # SACO CORRELACION DE PEARSON
	 correlacion.pearson = calculador_correlaciones(datos_expresion.corr, metodo = 'pearson')
	 
	 # SACO DISTANCIA DE MANHATTAN
	 datos_expresion %>%
	 # saco los tags que refieren al estudio, dato no menor ya que pasara a considerarse el conjunto de valores para cada gen como perteneciente a una misma entidad
	 #dplyr::select(-Study) %>%
	 # no incluido de momento # modifico valores de expresion menores a uno, reemplazandolos por NA
	 # no incluido de momento # dplyr::mutate_at(., colnames(datos_expresion)[-c(1)], ~if_else(. < 1, na_if(.,.), .)) 
	 # modifico, transformando los NA en cero, ya que stats::dist se comporta raro con los NA. Ademas, esta asuncion no es arriesgada, en efecto esos genes no se detectaron para el estudio hecho
	 dplyr::mutate_at(., colnames(.)[-c(1)], ~if_else(is.na(.), 0, .)) -> datos_expresion.paraManhattan
	 
	 correlacion.manhattan = calcula_distancias(datos_expresion.paraDistancias = datos_expresion.paraManhattan, metodo = 'manhattan')
	
	 # UNO LOS DATOS
	 # ordeno los datos
	 correlacion.spearman %<>% ordenador_tablas()
	 correlacion.pearson %<>% ordenador_tablas()
	 correlacion.manhattan %<>% ordenador_tablas()
	 
         # uno las tablas con la distancia cofenetica y las devuelvo con un tag
	 if (nrow(correlacion.spearman) > 0) { 
         dplyr::full_join(x = correlacion.spearman, y = cophenetic_dS, by = c('x' = 'x', 'y' = 'y')) %>%
         # agrego un tag con el organismo
         dplyr::mutate(organismo = org, correlacion = 'spearman') -> spearman.tibble } else {
	 spearman.tibble = tibble(x = NA, y = NA, organismo = org, correlacion = 'spearman')
	 }
	 
	  if (nrow(correlacion.pearson) > 0) { 
         dplyr::full_join(x = correlacion.pearson, y = cophenetic_dS, by = c('x' = 'x', 'y' = 'y')) %>%
         # agrego un tag con el organismo
         dplyr::mutate(organismo = org, correlacion = 'pearson') -> pearson.tibble } else {
	 pearson.tibble = tibble(x = NA, y = NA, organismo = org, correlacion = 'pearson')
	 }
	 
	  if (nrow(correlacion.manhattan) > 0) { 
         dplyr::full_join(x = correlacion.manhattan, y = cophenetic_dS, by = c('x' = 'x', 'y' = 'y')) %>%
         # agrego un tag con el organismo
         dplyr::mutate(organismo = org, correlacion = 'manhattan') -> manhattan.tibble } else {
	 manhattan.tibble = tibble(x = NA, y = NA, organismo = org, correlacion = 'manhattan')
	 }
	 
         # devuelvo la union
	 return(bind_rows(spearman.tibble, pearson.tibble, manhattan.tibble))
	 
         } 
    }
