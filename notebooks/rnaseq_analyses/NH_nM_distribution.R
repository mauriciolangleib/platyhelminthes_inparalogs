# plots analyzing distribtution of multimapping reads vs number of allowed mismatches

# loading BAM file with nM and NH tags
mapping_processed.list = Rsamtools::scanBam(file = '../results/STAR/smansoni_default/ERR022872_STAR_alnAligned.sortedByCoord.out.bam', param = ScanBamParam(tag=c("nM", "NH")))

# creating tibble
mapping_processed.tibble = tibble(nM = mapping_processed.list[[1]]$tag$nM,
                                  NH = mapping_processed.list[[1]]$tag$NH,
                                  run = 'ERR022872_STAR_alnAligned.sortedByCoord.out.bam')

# classifying based on thresholds
mapping_processed.tibble %>% 
  #group_by(nM, NH) %>% 
  #summarise(count = n()) %>%
  ggplot(data = ., mapping = aes(x = nM, y = NH)) + 
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") -> nM_nH.plot

## plotting
nM_nH.plot %>%
  ggsave(filename = 'nM_vs_NH.pdf', plot = ., device = 'pdf')

# trying another way
mapping_processed.tibble %>% 
  group_by(nM, NH) %>% 
  summarise(count = n()) %>%
  #dplyr::filter(NH == 0) %>%
  ungroup() %>%
  ggplot(data = ., mapping = aes(x = nM, y = count, group = NH, color = NH, fill = NH)) + 
  geom_point() + 
  geom_line() -> nM_vs_NH.lines

nM_vs_NH.lines %>%
  ggsave(filename = 'nM_vs_NH.points_and_lines.pdf', plot = ., device = 'pdf')

# geom density
mapping_processed.tibble %>% 
  ggplot(data = ., mapping = aes(x = NH)) + 
    geom_density(alpha = .1) +
    facet_wrap(~nM) -> nM_vs_NH.density

nM_vs_NH.density %>%
  ggsave(filename = 'nM_vs_NH.density.pdf', plot = ., device = 'pdf')
