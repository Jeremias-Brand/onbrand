#' Create a combined dataframe of chromosomes from two input files
#'
#' This function reads two tab-separated input files containing chromosome names and sizes,
#' filters and arranges them based on the specified prefixes, and creates a combined dataframe
#' with chromosome names, start, and end positions.
#'
#' @param chrom_first A character string specifying the path to the first input file containing
#'   chromosome names and sizes.
#' @param chrom_second A character string specifying the path to the second input file containing
#'   chromosome names and sizes.
#' @param chr_first A character string specifying the prefix for filtering and arranging chromosomes
#'   from the first input file. Defaults to "chr".
#' @param chr_second A character string specifying the prefix for filtering and arranging chromosomes
#'   from the second input file. Defaults to the value of `chr_first`.
#' @return A dataframe containing combined chromosomes with their names (chr), start positions (start),
#'   and end positions (end).
#' @importFrom dplyr filter arrange mutate select
#' @examples
#' # Example input files:
#' # chrom_first.txt
#' # chr1  200
#' # chr2  300
#' # chr3  150
#' # chrom_second.txt
#' # chr4  400
#' # chr5  100
#' # chr6  200
#'
#' # Example function call:
#' makechromosomes("chrom_first.txt", "chrom_second.txt")
makechromsomes <- function(chrom_first, chrom_second, chr_first = "chr", chr_second = chr_first) {
  h1 <- read.table(chrom_first)
  names(h1) <- c("chrom", "size")
  h1 <- filter(h1, grepl(chr_first, chrom)) %>%
    arrange(chrom)
  h1_xlim <- matrix(c(rep(0, nrow(h1)), h1$size), ncol=2)

  h2 <- read.table(chrom_second)
  names(h2) <- c("chrom", "size")
  h2 <- filter(h2, grepl(chr_second, chrom)) %>%
    arrange(desc(chrom))
  h2_xlim <- matrix(c(rep(0, nrow(h2)), h2$size), ncol=2)

  genome <- rbind(h1,h2)
  genome_xlim <- matrix(c(rep(0, nrow(genome)), genome$size), ncol=2)

  chromosomes <- genome %>%
    mutate(chr = chrom, start = 0, end = size) %>%
    select(chr, start, end)
  return(chromosomes)
}
