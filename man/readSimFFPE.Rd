\name{readSimFFPE}
\alias{readSimFFPE}
\title{
Simulate noisy NGS reads of FFPE samples for whole genome / several chromosomes
/ large regions
}
\description{
NGS data from FFPE samples contain numerous artificial chimeric reads. These 
chimeric reads are formed through the combination of two single-stranded DNA 
(ss-DNA). This function simulates these artificial reads as well as normal reads 
for FFPE samples on whole genome, or several chromosomes, or large regions.
}
\usage{
readSimFFPE(sourceSeq, referencePath, PhredScoreProfile, outFile, coverage, 
readLen = 150, meanInsertLen = 250, sdInsertLen = 80, enzymeCut = FALSE, 
chimericRatio = 0.08, localMatchRatio = 0.1, windowLen = 10000, 
matchWinLen = 10000, meanLogSeedLen = 1.7, sdLogSeedLen = 0.4, 
seedPassRate = 0.78, sdTargetDist = 120, sameStrandProb = 0.5, 
spikeWidth = 1500, betaShape1 = 0.5, betaShape2 = 0.5, 
sameTarRegionProb = 0, chimMutRate = 0.005, noiseRate = 0.0015, 
highNoiseRate = 0.08, highNoiseProb = 0.015, pairedEnd = TRUE, 
prefix = "SimFFPE", threads = 1, localChimeric = TRUE, 
distantChimeric = TRUE, normalReads = TRUE, overWrite = FALSE)
}
\arguments{
  \item{sourceSeq}{
A DNAStringSet object of DNA sequences used for simulation. It can cover the
entire reference genome or selected chromosomes or chromosome regions.
}
  \item{referencePath}{
Path to the reference genome.
}
  \item{PhredScoreProfile}{
A matrix representing the positional Phred score proportion. Each row of the 
matrix represents a position in the read (from begin to end), and each column 
the Phred quality score of base-calling error probabilities. The profile can be 
calculated from BAM file using the \code{\link{calcPhredScoreProfile}} function.
}
  \item{outFile}{
Output file path for the FASTQ file with simulated reads. Please include the
name of the output file without extension, e.g. "/tmp/sim1".
}
  \item{coverage}{
Coverage of the simulation.
}
  \item{readLen}{
Read length of the simulation.
}
  \item{meanInsertLen}{
Mean insert length for the simulation (normally distributed).
}
  \item{sdInsertLen}{
Standard deviation of the insert length for simulation (normally distributed).
}
  \item{enzymeCut}{
Simulate enzymatic fragmentation if it is set to true, otherwise simulate random
fragmentation.
}
  \item{chimericRatio}{
Proportion of artificial chimeric fragments (chimeric fragments / chimeric or 
normal fragments). Range: 0 to 1.
}
  \item{localMatchRatio}{
Proportion of adjacent ss-DNA combination (adjacent ss-DNA combination / 
adjacent or distant ss-DNA combination). Range: 0 to 1.
}
  \item{windowLen}{
The window length used in adjacent ss-DNA combination simulation. To simulate 
adjacent ss-DNA combinations, input DNA sequences are divided into small windows 
of equal size, and short complementary pairs are searched within the same window
. Suggested range: 5000-20000. Unit: base pair (bp).
}
  \item{matchWinLen}{
The target window length used in distant ss-DNA simulation. To simulate distant 
ss-DNA combinations, the target sequences are searched in a random window. 
Suggested range: 5000-20000. Unit: base pair (bp).
}
  \item{meanLogSeedLen}{
Mean of log scaled seed length (bp). Seeds are used to search for complementary 
targets. The mapping of seed and target links two ss-DNA together, yielding 
artificial chimeric fragments. The seed length follows a log-normal distribution
. See \code{\link[stats]{rlnorm}} for more details.
}
  \item{sdLogSeedLen}{
Standard deviation of log scaled seed length (bp).
}
  \item{seedPassRate}{
Proportion of seeds successfully forming chimeric fragments. Adjust this value 
when the percentage of chimeric reads in the output file is different from 
the parameter "chimericRatio".
}
  \item{sdTargetDist}{
Standard deviation of the normal distribution (mean = 0) used to simulate target 
selection probability. In adjacent ss-DNA combinations, when there are multiple 
targets for a seed, one target will be selected for combination. Target 
selection probability is simulated using the distance between seed and target. 
The smaller the distance, the larger the probability.
}
  \item{sameStrandProb}{
Probability of seed and target from the same DNA strand (same strand 
ss-DNA combination / same or complementary strand ss-DNA combination). Only 
valid for adjacent ss-DNA combination. For paired end sequencing, the larger the 
probability, the greater the proportion of improperly paired reads with LL / RR 
pair orientation, and the smaller with RL pair orientation. Range: 0 to 1.
}
  \item{spikeWidth}{
The width of chimeric read spike used to simulate distant ss-DNA combinations. 
In real FFPE samples, the chimeric reads formed by distant DNA combination are 
unevenly distributed along the chromosome. Some regions are enriched in these 
reads while some others are scarce. The length of these regions are of similar 
scale; therefore, a defined width is used for simulation. Suggested range: 
1500-2000. Unit: base pair (bp).
}
  \item{betaShape1}{
Shape parameter a of beta distribution used to model the unevenly distributed
distant ss-DNA combinations. The number of seeds in each "spike" follows a "U" 
shaped beta distribution. Use this parameter to adjust the shape of the curve. 
See \code{\link[stats]{rbeta}} for more details. Range: 0-1.
}
  \item{betaShape2}{
Shape parameter b of beta distribution used to model the unevenly distributed
distant ss-DNA combinations. The number of seeds in each "spike" follows a "U" 
shaped beta distribution. Use this parameter to adjust the shape of the curve. 
See \code{\link[stats]{rbeta}} for more details. Range: 0-1.
}
  \item{sameTarRegionProb}{
Probability of neighboring seeds to search targets in same random region for 
distant ss-DNA combination simulation. The larger the value, the more the false 
positive translocation variants.
}
  \item{chimMutRate}{
Mutation rate for each base in chimeric fragments. In the chimeric fragment 
formation process, biological-level errors might occur and lead to mutations on 
these artificial fragments. For all four basic types of nucleotides, the 
substitution probability is set equal. Range: 0-0.75.
}
  \item{noiseRate}{
Noise rate for each base in reads. This is used for sequencing-level errors.
The probability is set equal for all four basic types of nucleotides. Range: 
0-0.75.
}
  \item{highNoiseRate}{
A second noise rate for each base in reads. In some real sequencing data, 
some reads are much more noisy than others. This parameter can be used for this 
situation. Range: 0-0.75.
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
Generate reads from adjacent ss-DNA combinations if it is set to true. If it is 
set to false, skip this process.
}
  \item{distantChimeric}{
Generate reads from distant ss-DNA combinations if it is set to true. If it is
set to false, skip this process.
}
  \item{normalReads}{
Generate reads from normal fragments if it is set to true. If it is set to
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
(ss-DNA) fragments with short reverse complementary sequences. This function 
simulates these artificial reads as well as normal reads for FFPE samples 
on whole genome / several chromosomes / large regions. The combined ss-DNA may 
come from adjacent or distant regions. In the output fastq file these reads are 
distinguished by prefixes "localChimeric", "distantChimeric" and "Normal" in
their names. The parameter PhredScoreProfile can be calculated by the function
\code{\link{calcPhredScoreProfile}}. To simulate whole exome sequencing (WES) or
targeted sequencing, please use the function \code{\link{targetReadSimFFPE}}.
}
\value{
NULL. Reads will be written to the output FASTQ file.
}
\author{
Lanying Wei <lanying.wei@uni-muenster.de>
}
\note{
When fine-tuning is needed, simulate reads from certain areas / chromosomes
instead of the entire genome to save the runtime. Please check the package
vignette for the guidance of fine-tuning.
}
\seealso{
\code{\link{SimFFPE}}, \code{\link{calcPhredScoreProfile}},
\code{\link{targetReadSimFFPE}}
}
\examples{
PhredScoreProfilePath <- system.file("extdata", "PhredScoreProfile2.txt",
                                      package = "SimFFPE")
PhredScoreProfile <- as.matrix(read.table(PhredScoreProfilePath, skip = 1))
colnames(PhredScoreProfile) <- read.table(PhredScoreProfilePath, 
                                          nrows = 1, 
                                          colClasses = "character")

referencePath <- system.file("extdata", "example.fasta", package = "SimFFPE")
reference <- readDNAStringSet(referencePath)

## Simulate reads of the first three sequences of reference genome

sourceSeq <- reference[1:3]
outFile1 <- paste0(tempdir(), "/sim1")
readSimFFPE(sourceSeq, referencePath, PhredScoreProfile, outFile1,
            enzymeCut = FALSE, coverage=80, threads = 4)


## Simulate reads of defined regions on the first two sequences of reference
## genome

sourceSeq2 <- DNAStringSet(lapply(reference[1:2], function(x) x[1:10000]))
outFile2 <- paste0(tempdir(), "/sim2")
readSimFFPE(sourceSeq2, referencePath, PhredScoreProfile, outFile2,
            coverage = 80, enzymeCut = TRUE, threads = 1)


## Simulate reads of defined regions on the second and the third sequence of 
## reference genome and merge them with existing reads (a different prefix is 
## needed)

sourceSeq3 <- DNAStringSet(lapply(reference[2:3], function(x) x[1:10000]))
readSimFFPE(sourceSeq3, referencePath, PhredScoreProfile, outFile2,
            prefix = "simFFPE2", coverage = 80, enzymeCut = TRUE, 
            threads = 1, overWrite = FALSE)
}

