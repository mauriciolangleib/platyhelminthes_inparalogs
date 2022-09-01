correlaciones_spearman_cofeneticdS %>%
# filtro los que no tiene o bien distancia o bien r-pearson
  dplyr::filter(!is.na(r) & !is.na(`cophenetic dS`) & correlacion == 'pearson') %>%
# modifico labels
dplyr::mutate(organismo = case_when(organismo == 'hym' ~ 'H. microstoma', organismo == 'scm' ~ 'S. mansoni')) %>%
dplyr::mutate(tipo_datos = str_to_title(tipo_datos)) %>%
dplyr::rename(`Type of Data` = 'tipo_datos') %>%
ggplot(data = ., mapping = aes(x = `cophenetic dS`, y = r, color = `Type of Data`, fill = `Type of Data`)) +
geom_point() +
facet_wrap(~organismo, dir = 'v') +
labs(title = "Correlation between Pearson's r coefficient and cophenetic dS distance") +
xlab('Cophenetic dS distance') +
ylab("Pearson's r coefficient") -> pearson_vs_cofeneticDs.plot
#pearson_vs_cofeneticDs.plot
plotly::ggplotly(pearson_vs_cofeneticDs.plot)
