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

gviz_tlxdf_track = function(tlx_df, group=c("all", "sample", "group"), genome_info, chrom, range, name="", produce_density=F, split=c("strand", "none", "treatment"), show_legend=F, cex=1)
{
  tlx_df_plot = tlx_df %>% dplyr::mutate(strand=ifelse(Strand=="-1", "-", "+"))

  if(!split_strands) {
    tlx_df_plot = tlx_df_plot %>% dplyr::mutate(strand=".")
  }

  track_tlx = c()
  for(s in unique(tlx_df_plot$tlx_sample)) {
    tlx_df_s = tlx_df_plot %>% dplyr::filter(tlx_sample==s)

    if(produce_density) {
      density_df = tlx_df_s %>%
        dplyr::group_by(seqnames, strand) %>%
        dplyr::do((function(z) {
          zz<<-z
          d = density(z$start, from=range[1], to=range[2], bw=1e4, n=1000)
          d.bin = d$x[2L] - d$x[1L]
          d.count = zoo::rollsum(d$y, 2)*d.bin*nrow(z)/2

          data.frame(start=d$x[-length(d$x)], end=d$x[-1], score=ifelse(z$strand[1]=="-", -1, 1) * d.count)
        })(.)) %>%
        dplyr::ungroup()

      if(split_strands) {
        groups_factor = c("+", "-")
        density_df = density_df %>% reshape2::dcast(seqnames + start + end ~ strand, value.var="score")
      } else {
        groups_factor = NULL
        density_df = density_df %>% dplyr::select(seqnames, start, end, score)
      }

      density_ranges = GenomicRanges::makeGRangesFromDataFrame(density_df, seqinfo=genome_info, keep.extra.columns=T)
      return(Gviz::DataTrack(density_ranges, name=name, chromosome=chrom, type="hist", window=-1, windowSize=10000, groups=groups_factor, legend=show_legend, cex.legend=cex*1.5, cex.title=cex*1))
    } else {
      groups_factor = NULL
      if(split_strands) {
        groups_factor = rep(c("+", "-"), each=2)
        tlx_df_s = tlx_df_s %>%
          dplyr::group_by(seqnames, start, end) %>%
          dplyr::summarise(plus1=ifelse(any(strand=="+"), 0, NA_real_), plus2=ifelse(any(strand=="+"), 1, NA_real_), minus1=ifelse(any(strand=="-"), -1, NA_real_),  minus2=ifelse(any(strand=="-"), 0, NA_real_)) %>%
          dplyr::ungroup()
      } else {
        tlx_df_s = tlx_df_s %>%
          dplyr::group_by(seqnames, start, end) %>%
          dplyr::summarise(any1=0, any2=1) %>%
          dplyr::ungroup()
      }

      bed_sum_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df_s, seqinfo=genome_info, keep.extra.columns=T)
      track_tlx = c(track_tlx, Gviz::DataTrack(bed_sum_ranges, name=name, chromosome=chrom, showAxis=F, type="h", groups=groups_factor, legend=show_legend, cex.legend=cex*1.5, cex.title=cex*1))
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
  track_tlx = c()
  for(s in unique(tlx_df$tlx_sample)) {
    tlx_sdf = tlx_df %>% dplyr::filter(tlx_sample==s)
    track_tlx = c(track_tlx, gviz_tlxdf_track(tlx_sdf, chrom=chrom, range=range, genome_info=genome_info, name=s, produce_density=T, split_strands=T, show_legend=F, cex=1))
  }
  gviz_params = list(
        trackList=c(track_tlx, track_axis),
        from=range[1], to=range[2],
        fontcolor.title="black", background.title="transparent")
  do.call(Gviz::plotTracks, gviz_params)

  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction), keep.extra.columns=T, ignore.strand=T)
  tlxcov_df = tlx_coverage(tlx_df, group="sample", extsize=1e4, exttype="symmetrical")

}