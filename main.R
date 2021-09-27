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

  genes_df = readr::read_tsv("genomes/mm10/annotation/refGene.bed", col_names=F) %>%
    setNames(c("gene_chrom", "gene_start", "gene_end", "gene_name", "gene_score", "gene_strand")) %>%
    dplyr::filter(!grepl("_rev", gene_name)) %>%
    dplyr::mutate(gene_length=gene_end-gene_start) %>%
    dplyr::filter(gene_length>=200e3 & gene_chrom=="chr6") %>%
    dplyr::arrange(dplyr::desc(gene_length)) %>%
    data.frame()

  repeatmasker_df = repeatmasker_read("genomes/mm10/annotation/ucsc_repeatmasker.tsv")
  # repeatmasker_df %>%
  #   dplyr::arrange(repeatmasker_chrom, repeatmasker_start, repeatmasker_strand) %>%
  #   dplyr::mutate(repeatmasker_score=1, repeatmasker_name=paste0(repeatmasker_family, "--", repeatmasker_class, "--", repeatmasker_name)) %>%
  #   dplyr::select(repeatmasker_chrom, repeatmasker_start, repeatmasker_end, repeatmasker_name, repeatmasker_score, repeatmasker_strand) %>%
  #   readr::write_tsv(file="genomes/mm10/annotation/ucsc_repeatmasker.bed", col_names=F)

  samples_df = readr::read_tsv("data/tlx_samples.tsv") %>%
    dplyr::mutate(group_short=dplyr::case_when(
      group=="Perental cell (NXP010)"~"NXP010",
      group=="Ctnna2 allelic deletion (47/5)" ~ "fndr",
      group=="Ctnna2 allelic+promoter deletion (18/4)" ~ "fndr-prom1",
      group=="Ctnna2 allelic+promoter deletion (38/3)" ~ "fndr-prom2"
    )) %>%
    dplyr::mutate(treatment=ifelse(control, "DMSO", "APH"))

  copynumber_df = readr::read_tsv("data/copynumber.tsv")

  tlx_df = tlx_read_many(samples_df)
  tlx_df = tlx_mark_dust(tlx_df)
  tlx_df = tlx_mark_rand_chromosomes(tlx_df)
  tlx_df = tlx_mark_bait_chromosome(tlx_df)
  tlx_df = tlx_mark_bait_junctions(tlx_df, 1.5e6)
  tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)
  tlx_df = tlx_df %>% dplyr::inner_join(copynumber_df, by="Qname")

  #
  # Export data for IGV
  #
  # tlxcov_export_df = tlx_coverage(tlx_df, group="sample", extsize=5e2, exttype="symmetrical") %>%
  #   dplyr::filter(!tlx_control)
  # for(s in unique(tlxcov_export_df$tlx_sample)) {
  #   smpl = samples_df %>% dplyr::filter(sample==s)
  #   tlxcov_export_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov_export_df %>% dplyr::filter(tlx_sample==smpl$sample) %>% dplyr::select(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end, score=tlxcov_pileup), keep.extra.columns=T)
  #   tlxcov_file_basename = paste0(smpl$group_short, " - ", ifelse(smpl$control, "DMSO", "APH"), "_", smpl$group_i, " - ", smpl$sample, ".bedgraph")
  #   tlxcov_file = paste0("data/pileup5e2/", tlxcov_file_basename)
  #   rtracklayer::export.bedGraph(tlxcov_export_ranges, tlxcov_file)
  # }



  #
  # Analize variability using different normalization strategy
  #
  roi_df = readr::read_tsv("data/roi.tsv")
  roi_ranges = GenomicRanges::makeGRangesFromDataFrame(roi_df %>% dplyr::mutate(seqnames=roi_chrom, start=roi_start, end=roi_end), keep.extra.columns=T)
  libsizes_df = tlx_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::summarize(
      libsize_unnormalized=100,
      libsize_all=n(),
      libsize_nobait=sum(!tlx_is_bait_junction),
      libsize_norepeats=sum(is.na(tlx_repeatmasker_class)),
      libsize_nobait_norepeats=sum(!tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_baitchr_nobait=sum(tlx_is_bait_chromosome & !tlx_is_bait_junction),
      libsize_baitchr_norepeats=sum(tlx_is_bait_chromosome & is.na(tlx_repeatmasker_class)),
      libsize_baitchr_nobait_norepeats=sum(tlx_is_bait_chromosome & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_lowcomplex=sum(!tlx_has_dust),
      libsize_lowcomplex_nobait=sum(!tlx_has_dust & !tlx_is_bait_junction),
      libsize_lowcomplex_norepeats=sum(!tlx_has_dust & is.na(tlx_repeatmasker_class)),
      libsize_lowcomplex_nobait_norepeats=sum(!tlx_has_dust & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_lowcomplex_baitchr_nobait=sum(!tlx_has_dust & tlx_is_bait_chromosome & !tlx_is_bait_junction),
      libsize_lowcomplex_baitchr_norepeats=sum(!tlx_has_dust & tlx_is_bait_chromosome & is.na(tlx_repeatmasker_class)),
      libsize_lowcomplex_baitchr_nobait_norepeats=sum(!tlx_has_dust & tlx_is_bait_chromosome & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_cv=sum(copynumber_count<=2),
      libsize_cv_nobait=sum(copynumber_count<=2 & !tlx_is_bait_junction),
      libsize_cv_norepeats=sum(copynumber_count<=2 & is.na(tlx_repeatmasker_class)),
      libsize_cv_nobait_norepeats=sum(copynumber_count<=2 & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_cv_baitchr_nobait=sum(copynumber_count<=2 & tlx_is_bait_chromosome & !tlx_is_bait_junction),
      libsize_cv_baitchr_norepeats=sum(copynumber_count<=2 & tlx_is_bait_chromosome & is.na(tlx_repeatmasker_class)),
      libsize_cv_baitchr_nobait_norepeats=sum(copynumber_count<=2 & tlx_is_bait_chromosome & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_cv_lowcomplex=sum(copynumber_count<=2 & !tlx_has_dust),
      libsize_cv_lowcomplex_nobait=sum(copynumber_count<=2 & !tlx_has_dust & !tlx_is_bait_junction),
      libsize_cv_lowcomplex_norepeats=sum(copynumber_count<=2 & !tlx_has_dust & is.na(tlx_repeatmasker_class)),
      libsize_cv_lowcomplex_nobait_norepeats=sum(copynumber_count<=2 & !tlx_has_dust & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_cv_lowcomplex_baitchr_nobait=sum(copynumber_count<=2 & !tlx_has_dust & tlx_is_bait_chromosome & !tlx_is_bait_junction),
      libsize_cv_lowcomplex_baitchr_norepeats=sum(copynumber_count<=2 & !tlx_has_dust & tlx_is_bait_chromosome & is.na(tlx_repeatmasker_class)),
      libsize_cv_lowcomplex_baitchr_nobait_norepeats=sum(copynumber_count<=2 & !tlx_has_dust & tlx_is_bait_chromosome & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class))
    ) %>%
    reshape2::melt(id.vars="tlx_sample", variable.name="libsize_var", value.name="libsize_val") %>%
    dplyr::mutate(libsize_var=gsub("libsize_", "", as.character(libsize_var))) %>%
    dplyr::ungroup()


  #
  # Evaluate variability
  #
  tlx2roi_norm_df = as.data.frame(IRanges::mergeByOverlaps(roi_ranges, tlx_ranges)) %>%
    dplyr::inner_join(libsizes_df, by="tlx_sample") %>%
    dplyr::group_by(tlx_sample, roi_gene, roi_chrom, roi_start, roi_end, roi_peak, libsize_var, libsize_factor) %>%
    dplyr::summarize(breaks=dplyr::case_when(
      libsize_var=="unnormalized"~n(),
      libsize_var=="all"~n(),
      libsize_var=="nobait"~sum(!tlx_is_bait_junction),
      libsize_var=="norepeats"~sum(is.na(tlx_repeatmasker_class)),
      libsize_var=="nobait_norepeats"~sum(!tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_var=="baitchr_nobait"~sum(tlx_is_bait_chromosome & !tlx_is_bait_junction),
      libsize_var=="baitchr_norepeats"~sum(tlx_is_bait_chromosome & is.na(tlx_repeatmasker_class)),
      libsize_var=="baitchr_nobait_norepeats"~sum(tlx_is_bait_chromosome & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_var=="lowcomplex"~sum(!tlx_has_dust),
      libsize_var=="lowcomplex_nobait"~sum(!tlx_has_dust & !tlx_is_bait_junction),
      libsize_var=="lowcomplex_norepeats"~sum(!tlx_has_dust & is.na(tlx_repeatmasker_class)),
      libsize_var=="lowcomplex_nobait_norepeats"~sum(!tlx_has_dust & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_var=="lowcomplex_baitchr_nobait"~sum(!tlx_has_dust & tlx_is_bait_chromosome & !tlx_is_bait_junction),
      libsize_var=="lowcomplex_baitchr_norepeats"~sum(!tlx_has_dust & tlx_is_bait_chromosome & is.na(tlx_repeatmasker_class)),
      libsize_var=="lowcomplex_baitchr_nobait_norepeats"~sum(!tlx_has_dust & tlx_is_bait_chromosome & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_var=="cv"~sum(copynumber_count<=2),
      libsize_var=="cv_nobait"~sum(copynumber_count<=2 & !tlx_is_bait_junction),
      libsize_var=="cv_norepeats"~sum(copynumber_count<=2 & is.na(tlx_repeatmasker_class)),
      libsize_var=="cv_nobait_norepeats"~sum(copynumber_count<=2 & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_var=="cv_baitchr_nobait"~sum(copynumber_count<=2 & tlx_is_bait_chromosome & !tlx_is_bait_junction),
      libsize_var=="cv_baitchr_norepeats"~sum(copynumber_count<=2 & tlx_is_bait_chromosome & is.na(tlx_repeatmasker_class)),
      libsize_var=="cv_baitchr_nobait_norepeats"~sum(copynumber_count<=2 & tlx_is_bait_chromosome & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_var=="cv_lowcomplex"~sum(copynumber_count<=2 & !tlx_has_dust),
      libsize_var=="cv_lowcomplex_nobait"~sum(copynumber_count<=2 & !tlx_has_dust & !tlx_is_bait_junction),
      libsize_var=="cv_lowcomplex_norepeats"~sum(copynumber_count<=2 & !tlx_has_dust & is.na(tlx_repeatmasker_class)),
      libsize_var=="cv_lowcomplex_nobait_norepeats"~sum(copynumber_count<=2 & !tlx_has_dust & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class)),
      libsize_var=="cv_lowcomplex_baitchr_nobait"~sum(copynumber_count<=2 & !tlx_has_dust & tlx_is_bait_chromosome & !tlx_is_bait_junction),
      libsize_var=="cv_lowcomplex_baitchr_norepeats"~sum(copynumber_count<=2 & !tlx_has_dust & tlx_is_bait_chromosome & is.na(tlx_repeatmasker_class)),
      libsize_var=="cv_lowcomplex_baitchr_nobait_norepeats"~sum(copynumber_count<=2 & !tlx_has_dust & tlx_is_bait_chromosome & !tlx_is_bait_junction & is.na(tlx_repeatmasker_class))
    )) %>%
    dplyr::mutate(breaks_norm=breaks*libsize_factor) %>%
    dplyr::inner_join(samples_df, by=c("tlx_sample"="sample")) %>%
    dplyr::group_by(group, roi_gene, roi_chrom, roi_start, roi_end) %>%
    dplyr::mutate(breaks_norm_rel=breaks_norm/max(breaks_norm)) %>%
    dplyr::ungroup()

  # tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction), ignore.strand=T, keep.extra.columns=T)
  # tlx2roi_df = as.data.frame(IRanges::mergeByOverlaps(roi_ranges, tlx_ranges)) %>%
  #   dplyr::group_by(tlx_sample, roi_gene, roi_chrom, roi_start, roi_end, roi_peak) %>%
  #   dplyr::summarize(breaks=n())
  # tlx2roi_norm_df = tlx2roi_df %>%
  #   dplyr::inner_join(libsizes_df, by="tlx_sample") %>%
  #   dplyr::mutate(breaks_norm=breaks*libsize_factor) %>%
  #   dplyr::inner_join(samples_df, by=c("tlx_sample"="sample")) %>%
  #   dplyr::group_by(group, roi_gene, roi_chrom, roi_start, roi_end) %>%
  #   dplyr::mutate(breaks_norm_rel=breaks_norm/max(breaks_norm)) %>%
  #   dplyr::ungroup()

  #
  # Genes vs normalization strategies heatmap
  #
  x = tlx2roi_norm_df %>%
    dplyr::filter(!control & (roi_gene!="Ctnna2" | grepl("47|NXP010", group)) & grepl("Ccser1|Ctnna2|Grid2|Magi1|Sox5", roi_gene)) %>%
    dplyr::mutate(roi_gene=ifelse(roi_gene=="Ctnna2", paste(roi_gene, "(47/5 & NXP010)"), roi_gene)) %>%
    dplyr::group_by(roi_gene, libsize_var, group) %>%
    dplyr::summarize(sd=sd(breaks_norm_rel), mean=mean(breaks_norm_rel))
  pheatmap::pheatmap(x %>% dplyr::mutate(cov=sd/mean) %>% reshape2::dcast(libsize_var ~ roi_gene, value.var="cov", fun.aggregate=mean) %>% tibble::column_to_rownames("libsize_var"))

  #
  # Variability analysis
  #
  ggplot(x, aes(x=reorder(libsize_var, sd/mean, median), y=sd/mean, group=libsize_var)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(aes(fill=roi_gene), shape=21, size=2, position=position_jitter()) +
    coord_flip()

  #
  # Samples vs genes heatmap
  #
  x = tlx2roi_norm_df %>%
    dplyr::filter(libsize_var=="all" & !control & grepl("Ccser1|Ctnna2|Grid2|Magi1|Sox5", roi_gene)) %>%
    dplyr::mutate(tlx_sample=paste(group_short, tlx_sample), roi_gene=ifelse(roi_gene=="Ctnna2", paste(roi_gene, "(47/5 & NXP010)"), roi_gene)) %>%
    dplyr::group_by(roi_gene) %>%
    dplyr::mutate(breaks_norm=breaks_norm/sum(breaks_norm)) %>%
    reshape2::dcast(tlx_sample ~ roi_gene, value.var="breaks_norm") %>%
    tibble::column_to_rownames("tlx_sample")
  pheatmap::pheatmap(x, cluster_rows=F)

  # Now we can analize the data knowing which normalization to use
  tlx_fdf = tlx_df %>%
    dplyr::filter(!tlx_is_bait_junction & !tlx_is_rand_chrom & copynumber_count<=2) %>%
    dplyr::select(-Seq) %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::ungroup()
  libsizes_fdf = tlx_fdf %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::summarize(libsize_size=n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(libsize_factor=max(libsize_size)/libsize_size)
  tlx_franges = GenomicRanges::makeGRangesFromDataFrame(tlx_fdf %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction), ignore.strand=T, keep.extra.columns=T)
  tlx2roi_df = as.data.frame(IRanges::mergeByOverlaps(roi_ranges, tlx_franges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(tlx_sample, roi_gene) %>%
    dplyr::summarize(breaks_count=n()) %>%
    dplyr::inner_join(libsizes_fdf, by="tlx_sample") %>%
    dplyr::mutate(breaks_norm_count=breaks_count*libsize_factor) %>%
    dplyr::select(-libsize_factor) %>%
    dplyr::inner_join(samples_df, by=c("tlx_sample"="sample")) %>%
    dplyr::group_by(group, roi_gene) %>%
    dplyr::mutate(breaks_norm=breaks_norm_count/max(breaks_norm_count)) %>%
    dplyr::ungroup()

  pdf("report/difference_in_breaks.pdf", width=11, height=8)
  ggplot(tlx2roi_df, aes(x=group_short, y=breaks_norm, fill=treatment)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(aes(fill=treatment), shape=21, size=2, position=position_jitterdodge(jitter.width=0.2)) +
    coord_flip() +
    facet_wrap(~roi_gene, scales="free")
  dev.off()

  #
  # Export breaks with copy number >= 50
  #
  tlxcov_export_df = tlx_coverage(tlx_df %>% dplyr::filter(copynumber_count>=50 & !tlx_control), group="sample", extsize=5e2, exttype="symmetrical")
  for(s in unique(tlxcov_export_df$tlx_sample)) {
    smpl = samples_df %>% dplyr::filter(sample==s)
    tlxcov_export_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov_export_df %>% dplyr::filter(tlx_sample==smpl$sample) %>% dplyr::select(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end, score=tlxcov_pileup), keep.extra.columns=T)
    tlxcov_file_basename = paste0(smpl$group_short, " - ", ifelse(smpl$control, "DMSO", "APH"), "_", smpl$group_i, " - ", smpl$sample, ".bedgraph")
    tlxcov_file = paste0("data/cv50_pileup5e2/C_", tlxcov_file_basename)
    rtracklayer::export.bedGraph(tlxcov_export_ranges, tlxcov_file)
  }

  #
  # HERE we plot density plots (DOESN'T WORK YET)
  #
  tlxcov_fdf = tlxcov_df %>%
    dplyr::filter(tlxcov_chrom==chrom & (dplyr::between(tlxcov_start, range[1], range[2]) | dplyr::between(tlxcov_end, range[1], range[2]))) %>%
    dplyr::inner_join(libsizes_df %>% dplyr::filter(libsize_factor=="cv"), by="tlx_sample") %>%
    dplyr::mutate(tlxcov_pileup_norm=tlxcov_pileup*libsize_factor)
  ggplot(tlxcov_fdf) +
    geom_step(aes(x=tlxcov_start, y=tlxcov_pileup_norm, color=tlx_sample)) +
    facet_wrap(~libsize_var, scales="free")

  repeatmasker_df %>%
    dplyr::mutate(m=paste(repeatmasker_class, repeatmasker_family, repeatmasker_name, sep="__"), score=1) %>%
    dplyr::select(repeatmasker_chrom, repeatmasker_start, repeatmasker_end, m, score, repeatmasker_strand) %>%
    dplyr::arrange(gsub("chr", "", repeatmasker_chrom), repeatmasker_start) %>%
    readr::write_tsv("data/repeatmasker.bed", col_names=F)

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

  pdf("report/break_density_ctnna2.pdf", height=24, width=10)
  do.call(Gviz::plotTracks, gviz_params)
  dev.off()
}