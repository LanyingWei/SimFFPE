\name{calcPhredScoreProfile}
\alias{calcPhredScoreProfile}
\title{
Estimate Phred score profile for FFPE read simulation
}
\description{
Calculate Phred score profile from the entire BAM file or reads in subsampled 
regions.
}
\usage{
calcPhredScoreProfile(bamFilePath, mapqFilter = 0, maxFileSize = 1, 
targetRegions = NULL, subsampleRatio = NA, subsampleRegionLength = 1e+05, 
disableSubsampling = FALSE, threads = 1)
}
\arguments{
  \item{bamFilePath}{
BAM file to be processed.
}
  \item{mapqFilter}{
Filter for mapping quality. Reads with mapping quality below this value will be
excluded from calculation.
}
  \item{maxFileSize}{
The maximum file size (in GB) that allows processing of the entire BAM file. If 
disableSubsampling is set to false, BAM file larger than this size will be 
subsampled for calculation. 
}
  \item{targetRegions}{
A DataFrame or GenomicRanges object representing target regions for calculation. 
Use it for targeted sequencing / WES data, or when you need to manually select 
subsampled regions (set disableSubsampling to true in this case). If it is a 
DataFrame, the first column should be the chromosome, the second the start 
position and the third the end position. Please use one-based coordinate systems 
(the first base should be marked with 1 but not 0).
}
  \item{subsampleRatio}{
Subsample ratio. Together with subsampleRegionLength to determine subsampled 
regions. When subsampleRatio is not given, it will be assigned 
the value of maxFileSize divided by the input BAM file size. Range: 0 to 1.
}
  \item{subsampleRegionLength}{
Length of each subsampled region. Unit: base pair (bp).
}
  \item{disableSubsampling}{
Force to use the entire BAM file for calculation when set to true.
}
  \item{threads}{
Number of threads used. Multi-threading can speed up the process.
}
}
\details{
Calculate positional Phred score profile from the entire BAM file or reads in 
subsampled regions. A Phred score profile will be returned, which can then 
be used in read simulation.
}
\value{
A matrix will be returned. Each row of the matrix represents a position in the 
read (from begin to end), and each column the Phred quality score of 
base-calling error probabilities. The value in the matrix represents the 
positional Phred score proportion.
}
\author{
Lanying Wei <lanying.wei@uni-muenster.de>
}

\seealso{
\code{\link{SimFFPE}}, \code{\link{readSimFFPE}},
\code{\link{targetReadSimFFPE}}
}
\examples{
bamFilePath <- system.file("extdata", "example.bam", package = "SimFFPE")
regionPath <- system.file("extdata", "regionsBam.txt", package = "SimFFPE")
regions <- read.table(regionPath)
PhredScoreProfile <- calcPhredScoreProfile(bamFilePath, targetRegions = regions)
}
