# listing reads
list.files('../data/fastq_files', full.names = T, recursive = T, pattern = '.fastq.gz$') %>%
# creating a file architecture for trimmed reads
tibble(fastq_file = .) %>%
# creating directories
  dplyr::mutate(trimmed_out = fastq_file %>% str_replace_all(., 'fastq_files', 'trimmed_reads') %>% str_replace_all(., '.fastq.gz', '.trimmed.fastq.gz'),
                unpaired_out = trimmed_out %>% str_replace_all(., '_1', '_1_U') %>% str_replace_all(., '_2', '_2_U'),
                species = trimmed_out %>% str_split('/') %>% purrr::map_chr(4),
                SRA_tag = trimmed_out %>% str_split('/') %>% purrr::map_chr(5),
                `Sample tag` = fastq_file %>% str_split('/') %>% purrr::map_chr(6) %>% str_replace_all(., '_.*$', ''),
                species_dir = glue('../data/trimmed_reads/{species}'),
                SRA_dir = glue('{species_dir}/{SRA_tag}'),
                run_tag = fastq_file %>% str_replace_all(., '_..fastq.gz$', ''),
                genome_index = case_when(species == 'hymenolepis_microstoma' ~ '../data/genomes/indexes/hmicrostoma',
                                         species == 'schistosoma_mansoni' ~ '../data/genomes/indexes/smansoni'),
                STAR_dir = case_when(species == 'hymenolepis_microstoma' ~ '../results/STAR/hmicrostoma',
                                     species == 'schistosoma_mansoni' ~ '../results/STAR/smansoni'),
                sample_aln.dir = glue('{STAR_dir}/{SRA_tag}'),
                sample_prefix = glue('{sample_aln.dir}/{`Sample tag`}_STAR_aln'),
		STAR_dir.def = glue('{STAR_dir}_default'),
		sample_aln.dir.def = glue('{STAR_dir.def}/{SRA_tag}'),
		sample_prefix.def = glue('{sample_aln.dir.def}/{`Sample tag`}_STAR_aln')
                ) -> file_architecture.run

# pivoting table
dplyr::left_join(x = file_architecture.run %>% 
                      dplyr::filter(str_detect(fastq_file, '_1')) %>% 
                      dplyr::rename(trimmed_out_1 = 'trimmed_out', unpaired_out_1 = 'unpaired_out', fastq_file_1 = 'fastq_file'),
                 y = file_architecture.run %>% 
                      dplyr::filter(str_detect(fastq_file, '_2')) %>% 
                      dplyr::rename(trimmed_out_2 = 'trimmed_out', unpaired_out_2 = 'unpaired_out', fastq_file_2 = 'fastq_file') %>%
                      dplyr::select(-c(species, SRA_tag, species_dir, SRA_dir, `Sample tag`, genome_index, STAR_dir, sample_aln.dir, sample_prefix, STAR_dir.def, sample_aln.dir.def, sample_prefix.def)),
                 by = c('run_tag' = 'run_tag')) -> file_architecture.run

adapter.dir = '../data'
adapter.file = glue('{adapter.dir}/TruSeq3-PE.fa')

# creating STAR indexes
list('../data/genomes/indexes', '../data/genomes/indexes/hmicrostoma', '../data/genomes/indexes/smansoni', '../results', '../results/STAR') %>%
  purrr::map(., ~{
    if(!dir.exists('')){
      system(glue('mkdir {.x}'))
      }
  })

hmicrostoma_genome.dir = '../data/genomes/hmicrostoma/v11.0'
hmicrostoma_genome.file = 'hymenolepis_microstoma.PRJEB124.WBPS11.genomic_softmasked.fa'
hmicrostoma_gtf.file = 'hymenolepis_microstoma.PRJEB124.WBPS11.canonical_geneset.gtf'

smansoni_genome.dir = '../data/genomes/smansoni/v11.0'
smansoni_genome.file = 'schistosoma_mansoni.PRJEA36577.WBPS11.genomic_softmasked.fa'
smansoni_gtf.file = 'schistosoma_mansoni.PRJEA36577.WBPS11.canonical_geneset.gtf'

if(!file.exists('../data/genomes/indexes/hmicrostoma/SAindex')){
  system(glue('STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ../data/genomes/indexes/hmicrostoma/ --genomeFastaFiles {hmicrostoma_genome.dir}/{hmicrostoma_genome.file} --sjdbGTFfile {hmicrostoma_genome.dir}/{hmicrostoma_gtf.file} --sjdbOverhang 100'))
}

if(!file.exists('../data/genomes/indexes/smansoni/SAindex')){
  system(glue('STAR --runThreadN 16  --runMode genomeGenerate --genomeDir ../data/genomes/indexes/smansoni/ --genomeFastaFiles {smansoni_genome.dir}/{smansoni_genome.file} --sjdbGTFfile {smansoni_genome.dir}/{smansoni_gtf.file} --sjdbOverhang 100'))
}

# running trimming
file_architecture.run %>%
  dplyr::transmute(trimmomatic_run = pmap(list(species_dir, SRA_dir, fastq_file_1, fastq_file_2, trimmed_out_1, trimmed_out_2, unpaired_out_1, unpaired_out_2, genome_index, STAR_dir, sample_aln.dir, sample_prefix, STAR_dir.def, sample_aln.dir.def, sample_prefix.def), ~{
                                          # creating directories
                                          list('../data/trimmed_reads', ..1, ..2, ..10, ..11, ..13, ..14) %>%
                                            purrr::map(., ~{
                                              if(!dir.exists(.x)){
                                                system(glue('mkdir {.x}'))
                                                }
                                               })
    
                                          # trimming reads
                                          if (!file.exists(..5)) {
                                            system(glue('trimmomatic PE {..3} {..4} {..5} {..7} {..6} {..8} -threads 16 ILLUMINACLIP:{adapter.file}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36'))
                                            }
                   
                                          # running STAR 
                                          if (!file.exists(glue('{..12}Aligned.sortedByCoord.out.bam'))){
                                            # without multimappers
                                            system(glue('STAR --genomeDir {..9} --runThreadN 16 --readFilesIn {..5} {..6} --outFileNamePrefix {..12} --readFilesCommand zcat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard'))
                                            }
					  # running STAR with default parameters (no multimapping)
                                          if (!file.exists(glue('{..15}Aligned.sortedByCoord.out.bam'))){
                                            # without multimappers
                                            system(glue('STAR --genomeDir {..9} --runThreadN 16 --readFilesIn {..5} {..6} --outFileNamePrefix {..15} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard'))
                                            }	  
					})
                  )
