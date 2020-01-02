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
fragments with short reverse complementary sequences from nearby or distant
regions. This package simulates these artifacts as well as normal reads for
FFPE samples. The simulation can cover whole genome, or several chromosomes, or
large regions, or whole exome, or targeted regions. It also supports
enzymatic / random fragmentation and paired-end / single-end sequencing
simulations. Fine-tuning can be performed for desired simulation results, and
multi-threading can help reduce the runtime. Please check the package vignette
for the guidance of fine-tuning.

\packageIndices{SimFFPE}
There are three available functions for NGS read simulation of FFPE samples:

1. \code{\link{estimateSimParam}}: Estimate Phred score profile and other basic
parameters for FFPE read simulation.

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
\code{\link{estimateSimParam}}, \code{\link{readSimFFPE}},
\code{\link{targetReadSimFFPE}}
}
\examples{

PhredScoreProfilePath <- system.file("extdata", "PhredScoreProfile2.txt",
                                     package = "SimFFPE")
PhredScoreProfile <- as.matrix(read.table(PhredScoreProfilePath, skip = 1))
colnames(PhredScoreProfile) <- read.table(PhredScoreProfilePath, 
                                          nrows = 1, colClasses = "character")

referencePath <- system.file("extdata", "example.fasta", package = "SimFFPE")
reference <- readDNAStringSet(referencePath)

## Simulate reads of the first three sequences of the reference genome

sourceSeq <- reference[1:3]
outFile1 <- paste0(tempdir(), "/sim1")
readSimFFPE(sourceSeq, reference, PhredScoreProfile, outFile1, overWrite = TRUE,
            coverage = 80, readLen = 150, meanInsertLen = 250, sdInsertLen = 80,
            enzymeCut = TRUE, threads = 2)


## Simulate reads for targeted regions

bamFilePath <- system.file("extdata", "example.bam", package = "SimFFPE")
regionPath <- system.file("extdata", "regionsBam.txt", package = "SimFFPE")
regions <- read.table(regionPath)
PhredScoreProfile <- estimateSimParam(bamFilePath, calcPhredProfile = TRUE,
                                      targetRegions = regions)

outFile <- paste0(tempdir(), "/sim2")
targetReadSimFFPE(reference, PhredScoreProfile, targetRegions, outFile,
                  coverage = 120, readLen = 100, meanInsertLen = 150, 
                  sdInsertLen = 40, enzymeCut = FALSE)
}
