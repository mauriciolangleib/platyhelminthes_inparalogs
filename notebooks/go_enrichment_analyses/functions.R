# Defino funcion que corre el script map_to_slim.py
map_to_slim = function(GOterm, go_obo, slim_obo){
	linea <- paste("map_to_slim.py --term='", GOterm, "' --slim_out='direct' ", go_obo, " ", slim_obo, sep = "")
	system(linea, intern = TRUE) %>%
		grep(., pattern = "GO:", value = TRUE )-> resultado
	return(resultado)
}

# defining parserGOs
	parserGOs = function(tabla_gos, especie, genero, clase){
		# leo la tabla
		read.table(tabla_gos, sep = '\t', header = F) %>%
		# agrego la metadata puesta
		select(., V1, V3) %>%
		mutate(.,
			species = especie,
			genre = genero,
			class = clase,
			phylum = "platyhelminthes") %>%
		rename(., p_value = V3) %>%
		# devuelvo el data.frame
		return(.)
	}
