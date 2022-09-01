correlaciones_spearman_cofeneticdS %>%
# filtro los que no tiene o bien distancia o bien Manhattan distance
  dplyr::filter(!is.na(r) & !is.na(`cophenetic dS`) & correlacion == 'manhattan') %>%
# modifico labels
dplyr::mutate(organismo = case_when(organismo == 'hym' ~ 'H. microstoma', organismo == 'scm' ~ 'S. mansoni')) %>%
dplyr::mutate(tipo_datos = str_to_title(tipo_datos)) %>%
dplyr::rename(`Type of Data` = 'tipo_datos') %>%
# uno datos de anotacion de InterPro que dan en WBPS
dplyr::left_join(x = ., y = (index.tibble %>% dplyr::select(Gene, info)), by = c('x' = 'Gene')) %>%
dplyr::rename(info_x = 'info') %>%
dplyr::left_join(x = ., y = (index.tibble %>% dplyr::select(Gene, info)), by = c('y' = 'Gene')) %>%
dplyr::rename(info_y = 'info') %>%
  datatable(., extensions = 'Buttons', filter = 'top', options = list(
    pageLength = 15,
    autoWidth = TRUE,
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
  ))
