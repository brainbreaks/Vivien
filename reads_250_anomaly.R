
library(dplyr)
library(ggplot2)
devtools::load_all('breaktools')

main = function()
{
  repeatmasker_df = repeatmasker_read("genomes/mm10/annotation/ucsc_repeatmasker.tsv")
  samples_df = readr::read_tsv("data/tlx_samples.tsv") %>%
    dplyr::mutate(group_short=dplyr::case_when(
      group=="Perental cell (NXP010)"~"NXP010",
      group=="Ctnna2 allelic deletion (47/5)" ~ "fndr",
      group=="Ctnna2 allelic+promoter deletion (18/4)" ~ "fndr-prom1",
      group=="Ctnna2 allelic+promoter deletion (38/3)" ~ "fndr-prom2"
    ))

  tlx_df = tlx_read_many(samples_df)
  tlx_df = tlx_mark_dust(tlx_df)
  tlx_df = tlx_mark_rand_chromosomes(tlx_df)
  tlx_df = tlx_mark_bait_chromosome(tlx_df)
  tlx_df = tlx_mark_bait_junctions(tlx_df, 1.5e6)
  tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)

  tlx_fdf = tlx_df %>%
    dplyr::filter(!tlx_control & tlx_is_bait_chromosome & !tlx_is_bait_junction) %>%
    dplyr::mutate(Seq_length_grouped=cut(Seq_length, breaks=seq(10, 500, by=10)))

  ggplot(tlx_fdf) +
    geom_density(aes(x=Seq_length, color=tlx_sample))
  tlx_fdf.sum = tlx_fdf %>%
    dplyr::group_by(tlx_sample, Seq_length_grouped) %>%
    dplyr::summarize(dust_length=mean(tlx_dust_prop, na.rm=T), has_dust_prop=mean(tlx_has_dust, na.rm=T)/Seq_length)
  ggplot(tlx_fdf.sum) +
    geom_boxplot(aes(x=Seq_length_grouped, y=dust_length)) +
    labs(y="Proportion of junctions identified by dustmasker (low complexity)", y="") +
    theme_grey(base_size = 16) +
    coord_flip()

  tlx_f2df = tlx_fdf %>%
    dplyr::filter(Seq_length_grouped %in% c("(250,260]", "(290,300]")) %>%
    dplyr::group_by(Seq_length_grouped) %>%
    dplyr::mutate(Seq_length_grouped_n=n()) %>%
    tidyr::separate_rows(tlx_repeatmasker_class, sep=", ") %>%
    dplyr::group_by(tlx_repeatmasker_class) %>%
    dplyr::mutate(tlx_repeatmasker_class_n=n()) %>%
    dplyr::group_by(Seq_length_grouped, tlx_repeatmasker_class, Seq_length_grouped_n, tlx_repeatmasker_class_n) %>%
    dplyr::summarize(n=n()) %>%
    dplyr::group_by(tlx_repeatmasker_class) %>%
    dplyr::mutate(count=n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(count==2) %>%
    dplyr::group_by(tlx_repeatmasker_class) %>%
    dplyr::do((function(z) {
      zz<<-z
      z.matrix = matrix(c(
        z$n[z$Seq_length_grouped=="(250,260]"], z$Seq_length_grouped_n[z$Seq_length_grouped=="(250,260]"]-z$n[z$Seq_length_grouped=="(250,260]"],
        z$n[z$Seq_length_grouped=="(290,300]"], z$Seq_length_grouped_n[z$Seq_length_grouped=="(290,300]"]-z$n[z$Seq_length_grouped=="(290,300]"]), ncol=2, byrow=T)
      z.test = fisher.test(z.matrix)
      data.frame(class=z$tlx_repeatmasker_class[1], pvalue=z.test$p.value, odds=z.test$estimate)
    })(.))

  ggplot(tlx_f2df) +
    geom_bar(aes(x=reorder(class, odds), y=log2(odds), fill=-log10(pvalue)), stat="identity") +
    theme_grey(base_size = 16) +
    coord_flip()

  tlx_fdf = tlx_df %>%
    dplyr::mutate(tlx_dust_length=ifelse(is.na(tlx_dust_length), 0, tlx_dust_length)) %>%
    dplyr::filter(!is.na(tlx_dust_length)) %>%
    dplyr::mutate(tlx_dust_length_group=cut(tlx_dust_length, breaks=seq(-1, 1250, by=10))) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::mutate(sample_n=n()) %>%
    dplyr::group_by(tlx_sample, tlx_dust_length_group) %>%
    dplyr::summarize(breaks=n(), breaks_prop=breaks/sample_n[1]) %>%
    dplyr::filter(breaks>100) %>%
    data.frame()

  y = tlx_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::mutate(sample_n=n()) %>%
    dplyr::select(-dplyr::matches("copynumber_count")) %>%
    dplyr::inner_join(copynumber_df, by="Qname") %>%
    dplyr::mutate(copynumber_count_group=cut(copynumber_count, breaks=c(0,1,5,10,50,100,200,300,350, 400, 450, 500, 1000, 2000))) %>%
    dplyr::group_by(tlx_sample, copynumber_count_group) %>%
    dplyr::summarize(breaks=n(), breaks_prop=breaks/sample_n[1]) %>%
    dplyr::group_by(copynumber_count_group) %>%
    dplyr::mutate(breaks_mean=mean(breaks))

  ggplot(y) +
    geom_boxplot(aes(x=copynumber_count_group, y=breaks_prop, fill=breaks_mean))+
    labs(y="Proportion of breaks in sample", x="Number of sequence copies in\n reference genome") +
    coord_flip() +
    scale_y_continuous(breaks=seq(0, 1, by=0.1)) +
    theme_grey(base_size=18)

  y = tlx_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::mutate(sample_n=n()) %>%
    tidyr::separate_rows(tlx_repeatmasker_class, sep=", ") %>%
    dplyr::inner_join(x, by="Qname") %>%
    dplyr::mutate(alignment_n_group=cut(alignment_n, breaks=c(0,1,5,10,50,100,200,300,350, 400, 450, 500, 1000))) %>%
    dplyr::group_by(tlx_repeatmasker_class, tlx_sample, alignment_n_group) %>%
    dplyr::summarize(breaks=n(), breaks_prop=breaks/sample_n[1]) %>%
    dplyr::group_by(alignment_n_group) %>%
    dplyr::mutate(breaks_mean=mean(breaks))


  ggplot(y) +
    geom_boxplot(aes(x=alignment_n_group, y=log10(breaks_prop), fill=breaks_mean))+
    labs(y="Proportion of breaks in sample", x="Number of allignment") +
    facet_wrap(~tlx_repeatmasker_class) +
    coord_flip() +
    theme_grey(base_size=18)


  seqs_df = tlx_fdf %>%
    dplyr::filter(Seq_length_grouped %in% c("(250,260]", "(290,300]")) %>%
    tidyr::separate_rows(tlx_repeatmasker_class, sep=", ") %>%
    dplyr::filter(tlx_has_dust & !is.na(tlx_repeatmasker_class) & tlx_repeatmasker_class=="Simple_repeat") %>%
    dplyr::distinct(Qname, .keep_all=T) %>%
    tidyr::separate_rows(tlx_dust_coordinates, sep="; ") %>%
    tidyr::separate(tlx_dust_coordinates, into=c("dust_start", "dust_end")) %>%
    dplyr::mutate(dust_start=as.numeric(dust_start), dust_end=as.numeric(dust_end)) %>%
    dplyr::mutate(SeqDust=substr(Seq, dust_start, dust_end), dust_len=dust_end-dust_start) %>%
    dplyr::filter(dust_len>150)

  ggplot(seqs_df) +
    geom_density(aes(x=dust_len, color=Seq_length_grouped))
}