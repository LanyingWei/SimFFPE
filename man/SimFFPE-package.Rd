\name{SimFFPE-package}
\alias{SimFFPE-package}
\alias{SimFFPE}
\docType{package}
\title{
\packageTitle{SimFFPE}
}
\description{
\packageDescription{SimFFPE}
}
\details{

\packageDESCRIPTION{SimFFPE}
The NGS (Next-Generation Sequencing) reads from FFPE
(Formalin-Fixed Paraffin-Embedded) samples contain numerous artificial chimeric
reads. These reads are derived from the combination of two single-stranded DNA 
(ss-DNA) fragments with short reverse complementary sequences. The combined 
ss-DNA may come from adjacent or distant regions. This package simulates these 
artifacts as well as normal reads for FFPE samples. The simulation can cover 
whole genome, or several chromosomes, or large regions, or whole exome, or 
targeted regions. It also supports enzymatic / random fragmentation and 
paired-end / single-end sequencing simulations. Fine-tuning can be performed 
for desired simulation results, and multi-threading can help reduce the runtime. 
Please check the package vignette for the guidance of fine-tuning.
\packageIndices{SimFFPE}
There are three available functions for NGS read simulation of FFPE samples:

1. \code{\link{calcPhredScoreProfile}}: Calculate positional Phred score profile 
from BAM file for read simulation.

2. \code{\link{readSimFFPE}}: Simulate noisy NGS reads of FFPE samples on whole
genome, or several chromosomes, or large regions.

3. \code{\link{targetReadSimFFPE}}: Simulate noisy NGS reads of FFPE samples in
exonic / targeted regions.
}
\author{
\packageAuthor{SimFFPE}

Maintainer: \packageMaintainer{SimFFPE}
}

\keyword{ package }
\seealso{
\code{\link{calcPhredScoreProfile}}, \code{\link{readSimFFPE}},
\code{\link{targetReadSimFFPE}}
}
\examples{

PhredScoreProfilePath <- system.file("extdata", "PhredScoreProfile2.rds",
                                     package = "SimFFPE")
PhredScoreProfile <- as.matrix(read.table(PhredScoreProfilePath, skip = 1))
colnames(PhredScoreProfile) <- read.table(PhredScoreProfilePath, 
                                          nrows = 1, 
                                          colClasses = "character")

referencePath <- system.file("extdata", "example.fasta", package = "SimFFPE")
reference <- readDNAStringSet(referencePath)

## Simulate reads of the first three sequences of the reference genome

sourceSeq <- reference[1:3]
outFile1 <- paste0(tempdir(), "/sim1")
readSimFFPE(sourceSeq, referencePath, PhredScoreProfile, outFile1, 
            coverage = 80, enzymeCut = TRUE, threads = 4)

## Simulate reads for targeted regions

bamFilePath <- system.file("extdata", "example.bam", package = "SimFFPE")
regionPath <- system.file("extdata", "regionsBam.txt", package = "SimFFPE")
regions <- read.table(regionPath)
PhredScoreProfile <- calcPhredScoreProfile(bamFilePath, targetRegions = regions)

regionPath <- system.file("extdata", "regionsSim.txt", package = "SimFFPE")
targetRegions <- read.table(regionPath)

outFile <- paste0(tempdir(), "/sim2")
targetReadSimFFPE(referencePath, PhredScoreProfile, targetRegions, outFile,
                  coverage = 120, readLen = 100, meanInsertLen = 150, 
                  sdInsertLen = 40, enzymeCut = FALSE)
}
