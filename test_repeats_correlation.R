setwd("~/Workspace/Vivien")
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
library(Gviz)
devtools::load_all('breaktools')

main = function()
{
  repeatmasker_df = repeatmasker_read("genomes/mm10/annotation/ucsc_repeatmasker.tsv")
  repeatmasker_df %>%
    dplyr::arrange(repeatmasker_chrom, repeatmasker_start, repeatmasker_strand) %>%
    dplyr::mutate(repeatmasker_score=1, repeatmasker_name=paste0(repeatmasker_family, "--", repeatmasker_class, "--", repeatmasker_name)) %>%
    dplyr::select(repeatmasker_chrom, repeatmasker_start, repeatmasker_end, repeatmasker_name, repeatmasker_score, repeatmasker_strand) %>%
    readr::write_tsv(file="genomes/mm10/annotation/ucsc_repeatmasker.bed", col_names=F)
  x = repeatmasker_df %>%
    dplyr::filter(repeatmasker_family=="centr") %>%
    dplyr::arrange(repeatmasker_chrom, repeatmasker_start, repeatmasker_strand) %>%
    dplyr::mutate(repeatmasker_score=1, repeatmasker_name=paste0(repeatmasker_family, "--", repeatmasker_class, "--", repeatmasker_name)) %>%
    dplyr::select(repeatmasker_chrom, repeatmasker_start, repeatmasker_end, repeatmasker_name, repeatmasker_score, repeatmasker_strand) %>%
    readr::write_tsv(file="genomes/mm10/annotation/ucsc_repeatmasker_centr.bed", col_names=F)

  samples_df = readr::read_tsv("data/tlx_samples.tsv") %>%
    dplyr::mutate(group_short=dplyr::case_when(
      group=="Perental cell (NXP010)"~"NXP010",
      group=="Ctnna2 allelic deletion (47/5)" ~ "fndr",
      group=="Ctnna2 allelic+promoter deletion (18/4)" ~ "fndr-prom1",
      group=="Ctnna2 allelic+promoter deletion (38/3)" ~ "fndr-prom2"
    ))

  tlx_df = tlx_read_many(samples_df)
  tlx_df = tlx_mark_rand_chromosomes(tlx_df)
  tlx_df = tlx_mark_bait_chromosome(tlx_df)
  tlx_df = tlx_mark_bait_junctions(tlx_df, 1.5e6)
  tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)
  tlx_fdf = tlx_df %>%
    dplyr::filter(!tlx_is_bait_junction & !tlx_is_rand_chrom) %>%
    dplyr::select(-Seq) %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::ungroup()

  repeatmasker_df = repeatmasker_df %>%
    dplyr::mutate(repeatmasker_group=paste(repeatmasker_df$repeatmasker_class, "/", repeatmasker_df$repeatmasker_family, "/")) %>%
    dplyr::group_by(repeatmasker_group) %>%
    dplyr::mutate(repeatmasker_group_n=n()) %>%
    dplyr::ungroup()

  repeatmasker_width_df = data.frame()
  unique_groups = unique(repeatmasker_df$repeatmasker_group)
  for(rm in unique_groups) {
    writeLines(paste(which(rm==unique_groups), "/", length(unique_groups)))
    rm_df = repeatmasker_df %>% dplyr::filter(repeatmasker_group==rm)
    rm_ranges = GenomicRanges::makeGRangesFromDataFrame(rm_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end))
    rm_reduced_ranges = GenomicRanges::reduce(rm_ranges)
    repeatmasker_width_df = rbind(repeatmasker_width_df, data.frame(repeat_group=rm, repeats_count=nrow(rm_df), total_width=sum(width(rm_reduced_ranges))))
  }

  pdf("report/breaks_explained_by_repeats_sorted2.pdf", width=12, height=8)
  for(extend in c(0, 10, 1000, 10000)) {
    gc()
    repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(repeatmasker_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end), keep.extra.columns=T)
    tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_fdf %>% dplyr::mutate(seqnames=Rname, start=Junction-extend, end=Junction+extend, Strand="*"), keep.extra.columns=T)
    tlx2repeatmasker_df = as.data.frame(mergeByOverlaps(tlx_ranges, repeatmasker_ranges))

    x = tlx2repeatmasker_df %>%
      dplyr::filter(repeatmasker_group_n>=100) %>%
      dplyr::inner_join(samples_df, by=c("tlx_sample"="sample")) %>%
      dplyr::mutate(group_long=group) %>%
      dplyr::distinct(group_long, tlx_sample, repeatmasker_group, Junction, .keep_all=T) %>%
      dplyr::group_by(group_long, tlx_sample, repeatmasker_group) %>%
      dplyr::summarize(breaks_explained=n()) %>%
      dplyr::inner_join(repeatmasker_width_df, by=c("repeatmasker_group"="repeat_group")) %>%
      dplyr::inner_join(tlx_df %>% dplyr::group_by(tlx_sample) %>% dplyr::summarize(breaks_total=n()), by=c("tlx_sample")) %>%
      dplyr::mutate(breaks_explained=breaks_explained/breaks_total) %>%
      dplyr::mutate(facet=paste(group_long, tlx_sample)) %>%
      # dplyr::filter(breaks_explained>0.001) %>%
      dplyr::arrange(breaks_explained/total_width) %>%
      dplyr::mutate(i=1:n())
    x_reduced = x %>%
      dplyr::group_by(repeatmasker_group) %>%
      dplyr::arrange(dplyr::desc(breaks_explained)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    p = ggplot(x, aes(x=i, y=breaks_explained/total_width, color=repeatmasker_group)) +
      geom_point() +
      ggtitle(paste0("Extending breaks symetrically: ", extend)) +
      ggrepel::geom_text_repel(aes(label=repeatmasker_group), data=x %>% dplyr::group_by(repeatmasker_group) %>% dplyr::arrange(dplyr::desc(breaks_explained)) %>% dplyr::slice(1)) +
      labs(x="Sorted sample/repeat_type pairs", y="breaks explained / total repeats width") +
      guides(color="none")
    p = ggplot(x, aes(x=breaks_explained*100, y=total_width, color=repeatmasker_group)) +
      geom_point() +
      ggtitle(paste0("Extending breaks symetrically: ", extend)) +
      ggrepel::geom_text_repel(aes(label=repeatmasker_group), data=x_reduced) +
      guides(color="none") +
      labs(x="breaks explained, %", y="total width of repeats")
    print(p)
  }
  dev.off()
}