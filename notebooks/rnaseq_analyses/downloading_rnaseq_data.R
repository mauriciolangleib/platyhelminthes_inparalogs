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

# downloading data
rnaseq_readfiles.tibble %>%
  dplyr::transmute(., download_run = purrr::pmap(list(species_dir, sra_dir, run), ~{
                                                    # creating directories if they doesnt exist
                                                    list('../data', '../data/fastq_files', ..1, ..2) %>%
                                                        purrr::map(., ~{
                                                                if(!dir.exists(.x)){
                                                                        system(glue('mkdir {.x}'))
                                                                        }
                                                        })
               
                                                    # getting samples with parallel-fastq-dump
                                                    final_file_1 = glue('{..2}/{..3}_1.fastq.gz')
          
                                                    if(!file.exists(final_file_1)){
                                                        #system(glue('parallel-fastq-dump -t 16 -s {..3} -O {..2}'))
                                                         #system(glue('fastq-dump --split-3 --gzip -O {..2} {..3}'))
							 #system(glue('fastq-dump --split-3 --gzip -O {..2} -A {..3}'))
							  system(glue('fastq-dump --split-3 --gzip -O . -A {..3}'))  
						
						    # moving
						    system(glue('mv {..3}_*.fastq.gz {..2}'))
							
                                                    # compressing data with gunzip
                                                        #system(glue('gzip {..2}/{..3}.fastq'))
                                                            
                                                    # removing standar FastQ file
                                                    if(file.exists(glue('{..2}/{..3}.fastq'))){
                                                        system(glue('rm {..2}/{..3}.fastq'))
                                                            }
                                                        }
          
                                                      })
                    )
