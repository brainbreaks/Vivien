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

gviz_tlxdf_track = function(tlx_df, type=c("density", "breaks"), group=c("all", "group", "treatment", "sample", "group+treatment", "sample+treatment"), split=c("strand", "none", "treatment"), genome_info, chrom, range, normalize=T, show_legend=F, cex=1)
{
  tlx_df_plot = tlx_df %>%
    dplyr::mutate(Strand=ifelse(Strand=="-1", "-", "+")) %>%
    dplyr::mutate(
      track_split=dplyr::case_when(
        split[1]=="strand" ~ Strand,
        split[1]=="none" ~ "All",
        split[1]=="treatment" ~ ifelse(tlx_control, "Control", "Treatment"),
        T ~ "WHAAAAT?"
      )) %>%
    dplyr::mutate(track_split_sign=-(1 - 2*(track_split %in% c("+", "Treatment", "All")))) %>%
    dplyr::mutate(track_name=dplyr::case_when(
      group[1]=="all" ~ "All samples",
      group[1]=="group" ~ tlx_group,
      group[1]=="treatment" ~ ifelse(tlx_control, "Control", "Treatment"),
      group[1]=="sample" ~ tlx_sample,
      group[1]=="group+treatment" ~ paste0(tlx_group, "(", ifelse(tlx_control, "Control", "Treatment"), ")"),
      group[1]=="sample+treatment" ~ paste0(tlx_sample, "(", ifelse(tlx_control, "Control", "Treatment"), ")"),
    ))

  libsizes_df = tlx_df_plot %>%
    dplyr::group_by(track_name) %>%
    dplyr::summarize(libsize=n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(libfactor=ifelse(normalize, max(libsize)/libsize, 1))
  tlx_df_plot = tlx_df_plot %>%
    dplyr::filter(Rname==chrom & (dplyr::between(Rstart, range[1], range[2]) | dplyr::between(Rend, range[1], range[2]))) %>%
    dplyr::inner_join(libsizes_df, by="track_name")

  plot_ranges = list()
  split_factors_df = list()
  ylim = c(0,0)

  print("Calculating graphs data")
  for(s in unique(tlx_df_plot$track_name)) {
    tlx_df_s = tlx_df_plot %>% dplyr::filter(track_name==s)

    if(type=="density") {
      density_df = tlx_df_s %>%
        dplyr::group_by(libfactor, seqnames, track_split, track_split_sign) %>%
        dplyr::do((function(z) {
          zz<<-z
          d = density(z$start, from=range[1], to=range[2], bw=1e4, n=1000)
          d.bin = d$x[2L] - d$x[1L]
          d.count = zoo::rollsum(d$y, 2)*d.bin*nrow(z)/2
          data.frame(start=d$x[-length(d$x)], end=d$x[-1], score=z$track_split_sign[1] * d.count * z$libfactor[1])
        })(.)) %>%
        dplyr::ungroup() %>%
        reshape2::dcast(seqnames + start + end ~ track_split, value.var="score", fill=0)
      ymax = pmax(ylim[2], max(abs(density_df[,-(1:3), drop=F])))
      ylim = c(-ymax, ymax)
      split_factors_df[[s]] = colnames(density_df)[-(1:3)]
      plot_ranges[[s]] = GenomicRanges::makeGRangesFromDataFrame(density_df, seqinfo=genome_info, keep.extra.columns=T)
      #
      # track_tlx = Gviz::DataTrack(plot_ranges[[s]], chromosome=chrom, ylim=ylim, type="hist", window=-1, windowSize=10000, groups=split_factors_df[[s]])
      # Gviz::plotTracks(trackList=track_tlx)
    } else {
      if(split!="none") {
        split_factors_df[[s]] = rep(tlx_df_plot %>% dplyr::arrange(track_split_sign) %>% dplyr::distinct(track_split) %>% .$track_split, each=2)
        carpet_df = tlx_df_s %>%
          dplyr::group_by(track_name, seqnames, start, end) %>%
          dplyr::summarise(
            minus1=ifelse(any(track_split_sign<0), -1, NA_real_),
            minus2=ifelse(any(track_split_sign<0), 0, NA_real_),
            plus1=ifelse(any(track_split_sign>0), 0, NA_real_),
            plus2=ifelse(any(track_split_sign>0), 1, NA_real_),) %>%
          dplyr::ungroup()
        ylim = c(-1,1)
      } else {
        split_factors_df[[s]] = tlx_df_plot %>% dplyr::arrange(track_split_sign) %>% dplyr::distinct(track_split) %>% .$track_split
        carpet_df = tlx_df_s %>%
          dplyr::group_by(track_name, seqnames, start, end) %>%
          dplyr::summarise(any1=0, any2=1) %>%
          dplyr::ungroup()
        ylim = c(0,1)
      }
      plot_ranges[[s]] = GenomicRanges::makeGRangesFromDataFrame(carpet_df, seqinfo=genome_info, keep.extra.columns=T)
    }
  }

  print("Generatig graphs")

  track_tlx = c()
  for(s in names(plot_ranges)) {
    ss<<-s
    print(s)
    if(type=="density") {
      track_tlx = c(track_tlx, Gviz::DataTrack(plot_ranges[[s]], chromosome=chrom, start=range[1], end=range[2], name=s, ylim=ylim, type="hist", window=-1, windowSize=10000, groups=split_factors_df[[s]], legend=show_legend, cex.legend=cex*1.5, cex.title=cex*1))
    } else {
      track_tlx = c(track_tlx, Gviz::DataTrack(plot_ranges[[s]], chromosome=chrom, start=range[1], end=range[2], name=s, howAxis=F, type="h", groups=split_factors_df[[s]], legend=show_legend, cex.legend=cex*1.5, cex.title=cex*1))
    }
  }

  track_tlx
}

main = function()
{
  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info = with(readr::read_tsv("genomes/mm10/annotation/mm10.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
           GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm10", length(seqnames))))
  genome_info = genome_info[paste0("chr", c(1:19, "X", "Y"))]

  repeatmasker_df = repeatmasker_read("genomes/mm10/annotation/ucsc_repeatmasker.tsv")
  samples_df = readr::read_tsv("data/tlx_samples.tsv")
  tlx_df = tlx_read_many(samples_df)

    # libsize_df = tlx_df %>%
    # dplyr::group_by(tlx_group, tlx_group_i) %>%
    # dplyr::summarize(sample_size=sum(!tlx_control), control_size=sum(tlx_control))
  # baits_df = tlx_identify_baits(tlx_df, breaksite_size=19)


  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_mark_bait_chromosome(tlx_df)
  tlx_df = tlx_mark_bait_junctions(tlx_df, 1.5e6)
  tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)
  tlx_df = tlx_df %>%
    dplyr::filter(!tlx_is_bait_junction) %>%
    dplyr::select(-Seq) %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::ungroup()


  chrom = "chr6"
  range = c(76.8e6, 78.2e6)
  track_axis = Gviz::GenomeAxisTrack(chromosome=chrom, cex=5)
  # track_ideogram = Gviz::IdeogramTrack(genome="mm10", chromosome=chrom, cex=8)
  group_names=c("Perental cell (NXP010)"="founder", "Ctnna2 allelic deletion (47/5)"="-allele\n47/5", "Ctnna2 allelic+promoter deletion (18/4)"="-allele\n-prom1\n18/4", "Ctnna2 allelic+promoter deletion (38/3)"="-allele\n-prom2\n38/3")
  track_tlx = gviz_tlxdf_track(
    tlx_df %>% dplyr::mutate(tlx_sample=paste0(group_names[tlx_group], "\n", tlx_sample, "\n", ifelse(tlx_control, "DMSO", "APH"))),
    chrom=chrom,
    range=range,
    normalize=T,
    genome_info=genome_info,
    type="density",
    group="sample",
    split="strand",
    show_legend=F,
    cex=1)
  gviz_params = list(
        trackList=c(track_tlx, track_axis),
        from=range[1], to=range[2],
        fontcolor.title="black", background.title="transparent")

  pdf("report/break_density_ctnna2.pdf", height=20, width=10)
  do.call(Gviz::plotTracks, gviz_params)
  dev.off()
}