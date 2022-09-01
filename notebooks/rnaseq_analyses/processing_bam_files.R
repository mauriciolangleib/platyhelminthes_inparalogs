# processing mapped data into a unique table of feature counts

# creating variables regarding employed genomes
hmicrostoma_genome.dir = '../data/genomes/hmicrostoma/v11.0'
hmicrostoma_genome.file = 'hymenolepis_microstoma.PRJEB124.WBPS11.genomic_softmasked.fa'
hmicrostoma_gtf.file = 'hymenolepis_microstoma.PRJEB124.WBPS11.canonical_geneset.gtf'

smansoni_genome.dir = '../data/genomes/smansoni/v11.0'
smansoni_genome.file = 'schistosoma_mansoni.PRJEA36577.WBPS11.genomic_softmasked.fa'
smansoni_gtf.file = 'schistosoma_mansoni.PRJEA36577.WBPS11.canonical_geneset.gtf'

# listing bam files
list.files('../results/STAR', pattern = 'bam', full.names = T, recursive = T) %>%
# creating file architecture
tibble(bam_file = .) %>% 
  dplyr::mutate(species = bam_file %>% str_split('/') %>% purrr::map_chr(4),
                SRA_tag = bam_file %>% str_split('/') %>% purrr::map_chr(5),
                Sample = bam_file %>% str_split('/') %>% purrr::map_chr(6) %>% str_replace_all(., '_.*$', ''),
                expression_file.org_dir = glue('../results/expression_tables/{species}'),
                expression_file.dir = glue('{expression_file.org_dir}/{SRA_tag}'),
                `Multimapping reads` = case_when(str_detect(bam_file, '_default/') ~ 'Allowed (N <= 10)',
                                                 !str_detect(bam_file, '_default/') ~ 'Not allowed'),
                expression_file.out = glue('{expression_file.dir}/{species}_{Sample}.STAR.tsv'),
                gtf_file = case_when(str_detect(species, 'hmicrostoma') ~ glue('{hmicrostoma_genome.dir}/{hmicrostoma_gtf.file}'),
                                     str_detect(species, 'smansoni') ~ glue('{smansoni_genome.dir}/{smansoni_gtf.file}'))
                ) -> file_architecture.tibble

# creating directory allocating results
if(!dir.exists('../results/expression_tables')){
  system(glue('mkdir ../results/expression_tables'))
  }

# running Rsubread
file_architecture.tibble %>%
  dplyr::transmute(., rsubread_run = pmap(list(expression_file.org_dir, expression_file.dir, bam_file, gtf_file, expression_file.out, Sample, species, `Multimapping reads`), ~try({
# creating directories
list(..1, ..2) %>%
  purrr::map(., ~{
  if(!dir.exists(.x)){
    system(glue('mkdir {.x}'))
    }
  })

# reading BAM files with featureCounts
readed_bam = Rsubread::featureCounts(..3, annot.ext = as.character(..4), countMultiMappingReads = FALSE, isPairedEnd = TRUE, isGTFAnnotationFile = TRUE, GTF.attrType = 'gene_id', nthreads = 16) 

# creating tibble
rnaseq_counts = readed_bam$annotation %>% as_tibble() %>% dplyr::mutate(counts = readed_bam$counts, species = ..7, sample = ..6)
rnaseq_counts %<>% dplyr::rename_with(., .fn = ~{'counts'}, starts_with('counts')) %>% dplyr::mutate(`Multimapping reads` = ..8)

rnaseq_counts    
                                                }))
                  ) -> rnaquant.tibble

# hay que chequear que parece que algunos runs no salieron bien!!!!.
# por temas de memoria

rnaquant.tibble %<>%
   group_split(rsubread_run) %>% 
   purrr::keep(., ~{!'rsubread_run' %in% colnames(unnest(.x))}) %>% 
   purrr::map_dfr(., ~{.x %>% unnest()})


# pivoting, obtaining TPMs and TPM medians across replicates and saving files
## getting list of studies belonging to same condition
samples_metadata.tibble = list.files('../data/run_metadata', recursive = T, pattern = '.tsv$', full.names = T) %>%
  as.list() %>%
  purrr::map_dfr(., ~{readr::read_tsv(.x, col_names = T) %>% dplyr::select(Run, Condition) %>% unique()})

# calculating TPM from counts
## obtaining transcript lengths
transcript_lengths.tibble = list('../data/genomes/smansoni/v11.0/schistosoma_mansoni.PRJEA36577.WBPS11.annotations.gff3',
     '../data/genomes/hmicrostoma/v11.0/hymenolepis_microstoma.PRJEB124.WBPS11.annotations.gff3') %>%
  purrr::map_dfr(., ~{
  rtracklayer::readGFF(filepath= .x, version = 3) %>% as_tibble()
  }) %>%
  # selecting names, seq. start & end, type, ID and name
  dplyr::select(type, start, end, ID, Name) %>%
  # subsetting for transcripts
  dplyr::filter(type == 'gene') %>%
  # estimating transcript lengths
  dplyr::mutate(., length = abs(end-start)) %>%
  # selecting transcript Name and length
  dplyr::select(Name, length)

# defining function to calculate TPMs from counts and transcript length
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

rnaquant.tibble %>%
  dplyr::select(GeneID, counts, species, sample, `Multimapping reads`) %>%
  dplyr::left_join(x = ., y = samples_metadata.tibble, by = c('sample' = 'Run')) %>%
  dplyr::left_join(x = ., y = transcript_lengths.tibble, by = c('GeneID' = 'Name')) %>%
  group_by(species, sample, `Multimapping reads`) %>%
  dplyr::mutate(rate = counts/length,
                TPM =  rate*1e6/sum(rate)) %>%
  ungroup() %>%
  group_by(., `GeneID`, `species`, `Condition`, `Multimapping reads`) %>%
  dplyr::summarise(TPM = TPM, sample = sample, `TPM median` = median(TPM), `Multimapping reads` = `Multimapping reads`) -> rnaquant_finished.tibble

# saving table
rnaquant_finished.tibble %>%
      readr::write_tsv(., '../results/expression_tables/rnaquant_tpms.tsv', col_names = T)
