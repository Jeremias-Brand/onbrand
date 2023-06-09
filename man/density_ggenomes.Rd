% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/density_ggenomes.R
\name{density_ggenomes}
\alias{density_ggenomes}
\title{Calculate feature density for genomic data}
\usage{
density_ggenomes(seqs, data, binwidth = 1e+05, varname = "feat_density")
}
\arguments{
\item{seqs}{A dataframe obtained from reading a FASTA index file (FAI) with `gggenomes::read_fai`,
filtered to include only the desired sequences.}

\item{data}{A dataframe obtained from reading a BED file with `gggenomes::read_bed`, filtered to
include only the desired sequences.}

\item{binwidth}{A numeric value specifying the width of the bins used for calculating feature
density. Defaults to 1e5.}

\item{varname}{A character string specifying the name of the density variable in the output
dataframe. Defaults to "feat_density".}
}
\value{
A dataframe containing the calculated feature density.
}
\description{
This function takes the output of `gggenomes::read_fai` and `gggenomes::read_bed` as input,
computes the coverage and feature density, and returns a dataframe with the calculated density.
}
\examples{
# Assuming you have the required input files, "genome.fa.fai" and "bed_file.bed":
seqs <- gggenomes::read_fai("genome.fa.fai") \%>\%
  filter(grepl("Chr", seq_id))
data <- gggenomes::read_bed("bed_file.bed") \%>\%
  filter(grepl("Chr", seq_id))

# Calculate feature density with default parameters:
density_df <- density_ggenomes(seqs, data)
}
