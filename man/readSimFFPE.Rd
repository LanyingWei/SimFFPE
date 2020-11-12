\name{readSimFFPE}
\alias{readSimFFPE}
\title{
Simulate normal and artifact chimeric reads in NGS data of FFPE samples for 
whole genome / several chromosomes / large regions
}
\description{
NGS data from FFPE samples contain numerous artifact chimeric reads. These 
chimeric reads are formed through the combination of two single-stranded DNA 
(ss-DNA) with short reverse complementary regions (SRCR). This function 
simulates these artifact chimeric reads as well as normal reads for FFPE samples 
on the whole genome, or several chromosomes, or large regions. To simplify the 
simulation, the genome is divided into small windows, and SRCRs are found within
the same window (adjacent ss-DNA combination) or between different windows 
(distant ss-DNA simulation).
}
\usage{
readSimFFPE(sourceSeq, referencePath, PhredScoreProfile, outFile, coverage, 
readLen=150, meanInsertLen=250, sdInsertLen=80, enzymeCut=FALSE, 
chimericProp=0.1, sameChrProp=0.43, windowLen=5000, adjChimProp=0.63, 
sameStrandProp=0.65, meanLogSRCRLen=1.8, sdLogSRCRLen=0.55, maxSRCRLen=32, 
meanLogSRCRDist=4.7, sdLogSRCRDist=0.35, distWinLen=5000, spikeWidth = 1500, 
betaShape1=0.5, betaShape2=0.5, sameTarRegionProb=0, adjFactor=1.65, 
distFactor=1.65, chimMutRate=0.003, noiseRate=0.0015, highNoiseRate=0.08, 
highNoiseProp=0.01, pairedEnd=TRUE, prefix="SimFFPE", threads=1, 
adjChimeric=TRUE, distChimeric=TRUE, normalReads=TRUE, overWrite=FALSE)
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
  \item{chimericProp}{
Proportion of artifact chimeric fragments (chimeric fragments / chimeric or 
normal fragments). Range: 0 to 1.
}
  \item{sameChrProp}{
Proportion of artifact chimeric fragments that are derived from the combination 
of two ss-DNA coming from the same chromosome. Range: 0 to 1.
}
  \item{windowLen}{
The window length used in adjacent ss-DNA combination simulation. To simulate 
adjacent ss-DNA combinations, input DNA sequences are divided into small windows 
of equal size, and short reverse complementary regions are searched within the 
same window to form artifact chimeric fragments. Unit: base pair (bp).
}
  \item{adjChimProp}{
Proportion of adjacent ss-DNA combinations among same chromosomal ss-DNA 
combinations. Range: 0 to 1.
}
  \item{sameStrandProp}{
Proportion of same-strand ss-DNA combinations among adjacent ss-DNA combinations
. For paired end sequencing, the larger the proportion, the greater the 
proportion of improperly paired reads with LL / RR pair orientation, and the 
smaller the proportion with RL pair orientation. Range: 0 to 1.
}
  \item{meanLogSRCRLen}{
Mean of log scaled length of the short reverse complementary regions (SRCR) in 
artifact chimeric fragments. SRCRs links two ss-DNA together, yielding artifact 
chimeric fragments. The length of SRCR follows a log-normal distribution. See 
\code{\link[stats]{rlnorm}} for more details. Unit: base pair (bp).
}
  \item{sdLogSRCRLen}{
Standard deviation of log scaled length of the short reverse complementary 
regions.
}
  \item{maxSRCRLen}{
Maximum length of the short reverse complementary regions. Unit: base pair (bp).
}
  \item{meanLogSRCRDist}{
Mean of log scaled original genomic distance of the short reverse complementary 
regions(SRCR) in artifact chimeric fragments. SRCRs links two ss-DNA together, 
yielding artifact chimeric fragments. The distance of SRCR is the original 
genomic distance between the two short reverse complementary segments, which 
follows a log-normal distribution in simulation. For log-normal distribution, 
see \code{\link[stats]{rlnorm}} for more details. Unit: base pair (bp).
}
  \item{sdLogSRCRDist}{
Standard deviation of log scaled original genomic distance of the short reverse 
complementary regions(SRCR) in artifact chimeric fragments.
}
  \item{distWinLen}{
The window length used in distant ss-DNA simulation. To simulate distant 
ss-DNA combinations, the short reverse complementary regions(SRCR) are searched 
between different windows. Unit: base pair (bp).
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
Probability of two distant ss-DNA combination events coming from the same two 
different windows.
}
  \item{adjFactor}{
Increase this value if the number of simulated adjacent chimeric reads is 
smaller than expected (sameChrProp * adjChimProp), decrease if the opposite is 
true.
}
  \item{distFactor}{
Increase this value if the number of simulated distant chimeric reads is smaller 
than expected, decrease if the opposite is true.
}
  \item{chimMutRate}{
Mutation rate for each base in chimeric fragments. In the chimeric fragment 
formation process, biological-level errors might occur and lead to mutations on 
these artifact fragments. For all four basic types of nucleotides, the 
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
  \item{highNoiseProp}{
Proportion of reads to be simulated with highNoiseRate other than noiseRate.
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
  \item{adjChimeric}{
Generate reads from adjacent ss-DNA combinations if it is set to true. If it is 
set to false, skip this process.
}
  \item{distChimeric}{
Generate reads from distant ss-DNA combinations if it is set to true. If it is
set to false, skip this process.
}
  \item{normalReads}{
Generate reads from normal fragments if it is set to true. If it is set to
false, skip this process.
}
  \item{overWrite}{
Overwrite the file if file with the same output path exists and it is set to 
true. If file with same output path exists and it is set to false, reads will be
appended to the existing file.
}
}
\details{
The NGS (Next-Generation Sequencing) reads from FFPE (Formalin-Fixed 
Paraffin-Embedded) samples contain numerous artifact chimeric reads (ACRS), 
which can lead to false positive structural variant calls. These ACRs are 
derived from the combination of two single-stranded DNA (ss-DNA) fragments with 
short reverse complementary regions (SRCR). This function simulates these 
artifact chimeric reads as well as normal reads for FFPE samples on the whole 
genome / several chromosomes / large regions. To simplify the simulation, the 
genome is divided into small windows, and SRCRs are found within the same 
window (adjacent ss-DNA combination) or between different windows (distant 
ss-DNA simulation). For adjacent ss-DNA combination events, the original genomic 
distance between and strands of two combined SRCRs are also simulated based on 
real data. In the output fastq file, reads are distinguished by prefixes 
"adjChimeric", "distChimeric" and "Normal" in their names. The parameter 
PhredScoreProfile can be calculated by the function 
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
instead of the entire genome to save the run-time. Please check the package
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
colnames(PhredScoreProfile)  <- 
    strsplit(readLines(PhredScoreProfilePath)[1], "\t")[[1]]

referencePath <- system.file("extdata", "example.fasta", package = "SimFFPE")
reference <- readDNAStringSet(referencePath)

## Simulate reads of the first three sequences of reference genome

sourceSeq <- reference[1:3]
outFile1 <- paste0(tempdir(), "/sim1")
readSimFFPE(sourceSeq, referencePath, PhredScoreProfile, outFile1,
            enzymeCut = FALSE, coverage=80, threads = 2)


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
