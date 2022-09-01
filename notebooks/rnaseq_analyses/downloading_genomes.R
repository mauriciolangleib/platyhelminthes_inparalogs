# si no existe creo el directorio donde va a ser bajada la data
if(!dir.exists('../data/genomes')){
  system('mkdir ../data/genomes')
  }

# me bajo los proteomas de S. mansoni y H. microstoma de las versiones WBPSv9.0-14.0
tibble(links = c('ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/hymenolepis_microstoma/PRJEB124/hymenolepis_microstoma.PRJEB124.WBPS11.genomic_softmasked.fa.gz',
                 'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/hymenolepis_microstoma/PRJEB124/hymenolepis_microstoma.PRJEB124.WBPS11.canonical_geneset.gtf.gz',
                 'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS11.genomic_softmasked.fa.gz',
                 'ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS11/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS11.canonical_geneset.gtf.gz'),
       feature = c('hymenolepis_microstoma.PRJEB124.WBPS11.genomic_softmasked.fa.gz',
                   'hymenolepis_microstoma.PRJEB124.WBPS11.canonical_geneset.gtf.gz',
                   'schistosoma_mansoni.PRJEA36577.WBPS11.genomic_softmasked.fa.gz',
                   'schistosoma_mansoni.PRJEA36577.WBPS11.canonical_geneset.gtf.gz')
                ) %>%
  dplyr::mutate(WBPSv = rep(c('v11.0'), 4),
                organismo = c(rep('hmicrostoma', 2), rep('smansoni', 2)),
                org_dir = glue('../data/genomes/{organismo}'),
                version_dir = glue('{org_dir}/{WBPSv}')) %>%
# bajo la data
  dplyr::transmute(., run = pmap(list(org_dir, version_dir, feature, links), ~{
                                      # si no existen, creo los directorios
                                      c(..1, ..2) %>%
                                       as.list() %>%
                                       purrr::map(., ~{
                                        if (!dir.exists(.x)){
                                          system(glue('mkdir {.x}'))
                                          }
                                       })
    
                                      # si no esta bajado el archivo de proteinas, lo bajo
                                      if(!file.exists(glue('{..2}/{..3}'))){
                                        # ingreso a la carpeta en cuestion
                                        dir_trabajo = getwd()
                                        
                                        # voy al directorio donde se va a agarrar la version
                                        setwd(..2)
                                        
                                        # bajo el archivo
                                        system(glue('wget {..4}'))
                                        
                                        # vuelvo al dir donde estaba trabajando
                                        setwd(dir_trabajo)
                                         }
                                      })
                  )
