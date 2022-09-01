# listing wanted files
# creating tibble containing employed reads
tribble(~project_id, ~species,
        'ERP000427', 'schistosoma_mansoni',
        'ERP000427', 'schistosoma_mansoni',
        'ERP017466', 'schistosoma_mansoni',
        'ERP002113', 'schistosoma_mansoni',
        'ERP004459', 'hymenolepis_microstoma') -> rnaseq_readfiles.tibble

# getting connection from SRAdb
sqlfile <- '/home/mauricio/SRAmetadb.sqlite'
if(!file.exists('/home/mauricio/SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)

# getting data regarding employed runs
runs_metadata.tibble = list.files('../data/run_metadata', recursive = TRUE, pattern = '.tsv', full.names = T) %>%
        as.list() %>%
        purrr::map_dfr(., ~{readr::read_tsv(.x, col_names = TRUE)})

# getting sample SRA accesion numbers
rnaseq_readfiles.tibble %<>%
        dplyr::mutate(., getting_sra_Accessions = pmap(list(project_id), ~{
                                                                       sraConvert(in_acc = .x, 'sra', sra_con)
                                                                        })
                        ) %>%
        tidyr::unnest() %>%
        # creating tibble containing file and directory architectures employed to allocate the data
        dplyr::mutate(species_dir = glue('../data/fastq_files/{species}'),
                      sra_dir = glue('{species_dir}/{project_id}')) %>%
        # filtering to get only runs belonging to the selected set of experiments
        dplyr::filter(run %in% runs_metadata.tibble$Run)

# creating flag with FastQ.gz filename
rnaseq_readfiles.tibble %<>%
  dplyr::mutate(FastQ_gz.filename = glue('{run}.fastq.gz'))

# listing already downloaded files
downloaded_fastqgz = tibble(fastq_gz.fulldir.downloaded = list.files('../data/fastq_files/schistosoma_mansoni.mal.backup', full.names = T, recursive = T, pattern = 'fastq.gz$'),
                            fastq_gz.downloaded = list.files('../data/fastq_files/schistosoma_mansoni.mal.backup', recursive = T, pattern = 'fastq.gz$') %>% str_split('/') %>% purrr::map_chr(2)
                            )

# merging tables by FastQ.gz filename
dplyr::left_join(x = rnaseq_readfiles.tibble,
                 y = downloaded_fastqgz,
                 by = c('FastQ_gz.filename' = 'fastq_gz.downloaded')) %>%
 dplyr::mutate(dest_file = glue('{sra_dir}/{FastQ_gz.filename}')) %>% 
 dplyr::filter(!is.na(fastq_gz.fulldir.downloaded)) %>% 
 # moving those downloaded files that were wanted
 dplyr::transmute(., moving_files = pmap(list(fastq_gz.fulldir.downloaded, dest_file), ~{
                                                    system(glue('mv {..1} {..2}'))
                                                    })
                  )
