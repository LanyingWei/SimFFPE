\name{targetReadSimFFPE}
\alias{targetReadSimFFPE}
\title{
Simulate noisy NGS reads of FFPE samples in exonic / targeted regions
}
\description{
NGS read from FFPE samples contain numerous artificial chimeric reads. This
function simulates these artifacts as well as normal reads for FFPE samples
within defined regions.
}
\usage{
targetReadSimFFPE(reference, PhredScoreProfile, targetRegions, outFile, 
coverage, readLen, meanInsertLen, sdInsertLen, enzymeCut = FALSE,
covRatioFFPE = 0.08, localMatchRatio = 0.1, padding = NA, minGap = NA,
matchPadding = 5000, windowLen = 10000, matchWinLen = 10000, meanSeedLen = 5,
sdSeedLen = 3, revChimericProb = 0.8, spikeWidth = 1500, mutationRate = 0.005,
noiseRate = 0.0015, highNoiseRate = 0.08, highNoiseProb = 0.015,
pairedEnd = TRUE, prefix = "SimFFPE", threads = 1, localChimeric = TRUE,
distantChimeric = TRUE, normalReads = TRUE, overWrite = FALSE)
}
\arguments{
  \item{reference}{
A DNAStringSet object of the whole reference genome.
}
  \item{PhredScoreProfile}{
A matrix representing the positional Phred score proportion. Each row of the 
matrix represents a position in the read (from begin to end), and each column 
the Phred quality score of base-calling error probabilities. The profile can be 
calculated from BAM file using the \code{\link{estimateSimParam}} function.
}
  \item{targetRegions}{
A DataFrame or GenomicRanges object representing the exonic / targeted regions 
to simulate. If it is a DataFrame, the first column should be the chromosome, 
the second the start position and the third the end position. Please use 
one-based coordinate systems (the first base should be marked with 1 but not 0).
}
  \item{outFile}{
Output file path for the FASTQ file with simulated reads. Please include the
name of the output file, e.g. "/tmp/sim1".
}
  \item{coverage}{
Coverage of the simulated file.
}
  \item{readLen}{
Read length of the simulated reads.
}
  \item{meanInsertLen}{
Mean insert length for simulation (normally distributed).
}
  \item{sdInsertLen}{
Standard deviation of the insert length for simulation (normally distributed).
}
  \item{enzymeCut}{
Simulate enzymatic fragmentation if it is set to true, otherwise simulate random
fragmentation.
}
  \item{covRatioFFPE}{
Proportion of reads generated from artifact chimeric sequences (the number of
reads from chimeric sequences divided by the number of all reads including
chimeric reads and normal reads). Range: 0 to 1.
}
  \item{localMatchRatio}{
Proportion of chimeric reads that are formed by the combination of two nearby
DNA fragments (the number of chimeric reads of nearby DNA combination divided by
the number of all chimeric reads including chimeric reads derived from nearby or
distant DNA combination). Range: 0 to 1.
}
  \item{padding}{
Length of padding of the input target regions. The padding length will be
added to both sides of the target regions. If this value is not given,
the value of input meanInsertLen divided by two will be used.
Range: natural numbers. Unit: base pair (bp).
}
  \item{minGap}{
Minimal allowed length of gap between target regions. Regions with a gap smaller
than this value will be merged. If this value is not given, the value of input
readLen will be used. Range: natural numbers. Unit: base pair (bp).
}
  \item{matchPadding}{
Length of padding to the modified target regions for match searching in nearby
regions. The padding length will be added to both sides of the modified regions
(after padding, merging and dividing). The larger the matchPadding, the higher
the probability of successfully finding a short reverse complementary match, but
the greater the allowed maximum distance between two nearby DNA fragments for
combination and the longer the runtime. Suggested range: 2500-10000.
Unit: base pair (bp).
}
  \item{windowLen}{
The local window length used in the simulation of the combination of nearby DNA
fragments. To simulate nearby DNA fragments combination, if the length of a
certain target region is larger than the windowLen, this region will then be
divided into smaller regions of windowLen, and the short reverse
complementary sequences will be searched within the same region (plus match
padding). The larger the window, the higher the probability of successfully
finding a short reverse complementary match, but the greater the allowed maximum
distance between two nearby DNA fragments for combination and the longer the
runtime. Suggested range: 5000-20000. Unit: base pair (bp).
}
  \item{matchWinLen}{
The target window length used in the simulation of the combination of distant
DNA fragments. To simulate distant DNA fragments combination, the reference DNA
sequences are divided into small regions with same window length, and the short
reverse complementary sequence will be searched in a random region. The larger
the window, the higher the probability of successfully finding a short reverse
complementary match, but the longer the runtime. Suggested range: 5000-20000.
Unit: base pair (bp).
}
  \item{meanSeedLen}{
Mean length of short seed sequences used to search for reverse complementary
matches in chimeric sequences simulation (normally distributed length, with a
built-in minimum of 2 and maximum of 20). The longer it is, the lower the
probability of successfully finding a short reverse complementary match.
Suggested range: 4-8. Unit: base pair (bp).
}
  \item{sdSeedLen}{
Standard deviation of short seed sequences used to search for reverse
complementary matches in chimeric sequences simulation (normally distributed
length, with a built-in minimum of 2 and maximum of 20). Suggested range: 2-5.
}

  \item{revChimericProb}{
The probability of local-chimeric sequences to be generated through the
combination of a DNA fragment and its duplicate in reverse direction (the number
of local-chimeric sequences generated by the matches on the same strand in
reverse direction divided by the number of all local-chimeric sequences
generated by the matches on the same or complementary strand). For paired end
sequencing simulation, the larger the probability, the greater the proportion of
improper pair-orientated reads with LL / RR pair orientation, and the smaller
with RL pair orientation. For distant-chimeric sequences, this probability is
with a built-in value of 0.5. Range: 0 to 1. Suggested range: 0.7-0.9 based on
the real FFPE samples to simulate.
}

  \item{spikeWidth}{
The width of chimeric read spike used in the simulation of the combination of
distant DNA fragments. In real FFPE samples, the chimeric reads formed by
distant DNA combination are unevenly distributed along the chromosome. Some
regions are enriched in these distant-chimeric reads and some are not. The
length of these regions are of similar scale and therefore a defined width is
used for simulation. Suggested range: 1500-2000. Unit: base pair (bp).
}
  \item{mutationRate}{
Mutation rate for each base in chimeric sequences. In the whole chimeric
sequence formation process, some biological-level errors might occur and lead to
mutations. The probability is set to be even for all four basic type of
nucleotides. Range: 0-0.75.
}
  \item{noiseRate}{
Noise rate for each base in reads. This is used for sequencing-level errors.
The probability is set to be even for all four basic type of nucleotides.
Range: 0-0.75.
}
  \item{highNoiseRate}{
A second noise rate for each base in reads. This is used to simulate the
sequencing-level errors. The probability is set to be even for all four basic
type of nucleotides. In some real sequencing data, some reads are much more
noisy than others and this parameter can be used to simulate the noises on these
reads. Range: 0-0.75.
}
  \item{highNoiseProb}{
Probability of reads to be simulated with highNoiseRate other than noiseRate.
Range: 0-1.
}
  \item{pairedEnd}{
Simulate paired end sequencing when set to true.
}
  \item{prefix}{
Prefix for read names. When reads from different runs of simulation have to be
merged, please make sure that they have different prefixes.
}
  \item{threads}{
Number of threads used. Multi-threading can speed up the process.
}
  \item{localChimeric}{
Generate reads from local-chimeric sequences if it is set to true. If it is set
to false, skip this process.
}
  \item{distantChimeric}{
Generate reads from distant-chimeric sequences if it is set to true. If it is
set to false, skip this process.
}
  \item{normalReads}{
Generate reads from normal sequences if it is set to true. If it is set to
false, skip this process.
}
  \item{overWrite}{
Overwrite the file if file with same output path exists and it is set to true.
If file with same output path exists and it is set to false, reads will be
appended to the existing file.
}
}
\details{
The NGS (Next-Generation Sequencing) reads from FFPE
(Formalin-Fixed Paraffin-Embedded) samples contain numerous artificial chimeric
reads. These reads are derived from the combination of two single-stranded DNA
fragments with short reverse complementary sequences. This function simulates
these artifacts as well as normal reads for FFPE samples within defined regions.
The DNA fragments to be combined can be derived from nearby regions or
distant regions. In the output fastq file these reads are distinguished by
"localChimeric", "distantChimeric" and "Normal" in their names. The parameter
PhredScoreProfile can be derived by the function \code{\link{estimateSimParam}}.
To simulate whole genome sequencing (WGS) or to simulate reads on several large
regions / full chromosomes, please use the function \code{\link{readSimFFPE}}.
}
\value{
NULL. Reads will be written to the output FASTQ file.
}
\author{
Lanying Wei <lanying.wei@uni-muenster.de>
}
\note{
When fine-tuning is needed, simulate reads from part of the regions instead of
all the target regions to save the runtime. Please check the package vignette
for the guidance of fine-tuning.
}
\seealso{
\code{\link{SimFFPE}}, \code{\link{estimateSimParam}}, \code{\link{readSimFFPE}}
}
\examples{
PhredScoreProfilePath <- system.file("extdata", "PhredScoreProfile1.txt",
                                     package = "SimFFPE")
PhredScoreProfile <- as.matrix(read.table(PhredScoreProfilePath, skip = 1))
colnames(PhredScoreProfile) <- read.table(PhredScoreProfilePath, 
                                          nrows = 1, colClasses = "character")

referencePath <- system.file("extdata", "example.fasta", package = "SimFFPE")
reference <- readDNAStringSet(referencePath)

regionPath <- system.file("extdata", "regionsSim.txt", package = "SimFFPE")
targetRegions <- read.table(regionPath)

outFile <- paste0(tempdir(), "/sim3")
targetReadSimFFPE(reference, PhredScoreProfile, targetRegions, outFile,
                  coverage = 120, readLen = 100, meanInsertLen = 150, 
                  sdInsertLen = 40, enzymeCut = FALSE)

}
