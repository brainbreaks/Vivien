library(dplyr)
library(ggplot2)
library(readr)
devtools::load_all('breaktools')

main = function() {
  tlx_df = tlx_read_many(samples_df) %>%
    dplyr::mutate(QSeq=substr(Seq, Qstart, Qend))

  writeLines(with(tlx_df, paste0(">", Qname, "\n", QSeq)), con="tmp/qseq.fasta")
  qsec_allignments = readr::read_tsv("tmp/qseq", col_names=F, skip=68)

  qsec_allignments_df = readr::read_tsv("tmp/qseq", col_names=F, skip=68) %>%
    dplyr::group_by(X1) %>%
    dplyr::summarize(n=n())

  qsec_allignments_df = qsec_allignments_df %>%
    setNames(c("Qname", "copynumber_count")) %>%
    readr::write_tsv("")
}