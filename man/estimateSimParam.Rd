\name{estimateSimParam}
\alias{estimateSimParam}
\title{
Estimate Phred score profile and other basic parameters for FFPE read simulation
}
\description{
Estimate Phred score profile, read length, median / MAD of insert length and
some FFPE-artifact related statistics from the whole BAM file or reads in
subsampled regions.
}
\usage{
estimateSimParam(bamFilePath, mapqFilter = 0, maxFileSize = 1,
targetRegions = NULL, subsampleRatio = NA, subsampleRegionLength = 1e+05,
disableSubsampling = FALSE, threads = 1, calcPhredProfile = TRUE)
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
Maximum allowed file size to process on the whole BAM file (unit:GB). BAM file
larger than this value will be subsampled for parameter estimation if
disableSubsampling is set to be false.
}
  \item{targetRegions}{
Target regions for calculation. Use it for targeted sequencing / WES data, or
whenever manual selected subsample regions are desired (set
disableSubsampling to true in this case).
}
  \item{subsampleRatio}{
Ratio for subsampling. Together with subsampleRegionLength to determine
the subsampled regions. When subsampleRatio is not given, it will be
determined as the maxFileSize divided by the input BAM file size.
}
  \item{subsampleRegionLength}{
Length of each subsampled region. Unit: base pair (bp).
}
  \item{disableSubsampling}{
Force to use the whole BAM file for calculation when set to true.
}
  \item{threads}{
Number of threads used. Multi-threading can speed up the process.
}
  \item{calcPhredProfile}{
Calculate and return positional Phred score profile when set to true.
Phred score profile is needed for the read simulation.
}
}
\details{
Estimate Phred score profile, read length, median / median absolute deviation
of insert length and some FFPE-artifact related statistics from the whole
BAM file or reads in the subsampled regions. Phred score profile will be
returned if \code{calcPhredProfile} is set to true. This function can be used to
help determine some parameters used in FFPE read simulation. It can also be
applied to the simulated bam file / paired fresh frozen sample of FFPE sample
for comparison of certain statistics.
}
\value{
A data frame object will be returned if calcPhredProfile is set to true.
Each row represents a position in the read (from begin to end), and each column
the Phred quality score of base-calling error probabilities. The value in the
data frame represents the positional Phred score proportion.
}

\author{
Lanying Wei <lanying.wei@uni-muenster.de>
}

\seealso{
\code{\link{SimFFPE}}, \code{\link{readSimFFPE}},
\code{\link{targetReadSimFFPE}}
}
\examples{
bamFilePath1 <- system.file("extdata", "example1.bam", package = "SimFFPE")
PhredScoreProfile <- estimateSimParam(bamFilePath1)

bamFilePath2 <- system.file("extdata", "example2.bam", package = "SimFFPE")
regionPath <- system.file("extdata", "regionsBam.txt", package = "SimFFPE")
regions <- fread(regionPath)
estimateSimParam(bamFilePath2, calcPhredProfile = FALSE,
                 targetRegions = regions)
}
