library(dplyr)
library(reshape2)
library(readr)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(tidyr)
library(baseline)
library(smoother)
library(dbscan)
library(ggpmisc)
library(rtracklayer)
library(ggbeeswarm)
library(Rtsne)
library(pracma)
devtools::load_all('breaktools')

main = function()
{
  repeatmasker_df = repeatmasker_read("data/mm10/annotation/ucsc_repeatmasker.tsv")
  samples_df = readr::read_tsv("data/tlx_samples.tsv")
  tlx_df = tlx_read_many(samples_df)


  baits_df = tlx_identify_baits(tlx_df, breaksite_size=19)
  libsize_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_group_i) %>%
    dplyr::summarize(sample_size=sum(!tlx_control), control_size=sum(tlx_control))


  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_mark_bait_chromosome(tlx_df)
  tlx_df = tlx_mark_bait_junctions(tlx_df, 1.5e6)
  tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)
  tlx_df = tlx_df %>%
    dplyr::filter(!tlx_is_bait_junction) %>%
    dplyr::select(-Seq) %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::ungroup()
  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction), keep.extra.columns=T, ignore.strand=T)

}