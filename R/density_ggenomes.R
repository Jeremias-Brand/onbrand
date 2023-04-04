#' Calculate feature density for genomic data
#'
#' This function takes the output of `gggenomes::read_fai` and `gggenomes::read_bed` as input,
#' computes the coverage and feature density, and returns a dataframe with the calculated density.
#'
#' @param seqs A dataframe obtained from reading a FASTA index file (FAI) with `gggenomes::read_fai`,
#'   filtered to include only the desired sequences.
#' @param data A dataframe obtained from reading a BED file with `gggenomes::read_bed`, filtered to
#'   include only the desired sequences.
#' @param binwidth A numeric value specifying the width of the bins used for calculating feature
#'   density. Defaults to 1e5.
#' @param varname A character string specifying the name of the density variable in the output
#'   dataframe. Defaults to "feat_density".
#' @return A dataframe containing the calculated feature density.
#' @importFrom dplyr filter rename
#' @importFrom GenomicRanges makeGRangesFromDataFrame coverage tileGenome binnedAverage
#' @examples
#' # Assuming you have the required input files, "genome.fa.fai" and "bed_file.bed":
#' seqs <- gggenomes::read_fai("genome.fa.fai") %>%
#'   filter(grepl("Chr", seq_id))
#' data <- gggenomes::read_bed("bed_file.bed") %>%
#'   filter(grepl("Chr", seq_id))
#'
#' # Calculate feature density with default parameters:
#' density_df <- density_ggenomes(seqs, data)
density_ggenomes <- function(seqs, data, binwidth=1e5, varname = "feat_density") {
  # working with gggenome imported files
  # seqs:  gggenomes::read_fai("genome.fa.fai") %>%
  # filter(grepl("Chr", seq_id))
  # data: gggenomes::read_bed(bed) %>% filter(grepl("Chr", seq_id))
  glen <- seqs$length
  names(glen) <- seqs$seq_id

  satGR <- data %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seq_id",
                                            seqinfo = glen )
  cov <- coverage(satGR)
  bins <- tileGenome(glen, tilewidth=binwidth, cut.last.tile.in.chrom=TRUE)
  # this object is useful in the GenomicRanges universe but not good for gggenomes
  den <- binnedAverage(bins, cov, varname = varname)
  # We need a few tweaks
  df <- as.data.frame(den) %>%
    dplyr::rename(seq_id = seqnames)
  # only specific strand symbols are allowed.
  df$strand <- "."

  return(df)
}
