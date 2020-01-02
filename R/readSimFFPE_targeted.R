targetReadSimFFPE <- function(reference, PhredScoreProfile,
                              targetRegions, outFile, coverage,
                              readLen, meanInsertLen, sdInsertLen,
                              enzymeCut=FALSE, covRatioFFPE=0.1,
                              localMatchRatio=0.1,
                              padding=NA, minGap=NA, matchPadding=5000,
                              windowLen=10000, matchWinLen=10000,
                              meanSeedLen=5, sdSeedLen=3,
                              revChimericProb=0.8, spikeWidth=1500,
                              mutationRate=0.005, noiseRate=0.0015,
                              highNoiseRate=0.08, highNoiseProb=0.015,
                              pairedEnd=TRUE, prefix="SimFFPE", threads=1,
                              localChimeric=TRUE, distantChimeric=TRUE,
                              normalReads=TRUE, overWrite=FALSE) {
    
    if (nrow(PhredScoreProfile) != readLen) {
        stop("The number of rows in PhredScoreProfile should be the same as", 
             " the input readLen")
    }
    
    if (overWrite) {
        overWriteFastq(outFile = outFile, pairedEnd = pairedEnd)
    }

    covFFPE <- coverage * covRatioFFPE / 0.85
    covNorm <- coverage - covFFPE

    if(!pairedEnd) {
        enzymeCut <- FALSE
    }

    reference <- reference[width(reference) >= matchWinLen]
    names(reference) <- unlist(lapply(names(reference),
                                      function(x) strsplit(x, " ")[[1]][1]))
    
    if (is(targetRegions, "GenomicRanges")) {
        targetRegions <- data.frame(seqnames(targetRegions), 
                                    start(targetRegions), 
                                    end(targetRegions))
    }
    colnames(targetRegions) <- c("chr", "start", "end")
    targetRegions <- targetRegions[targetRegions$chr %in% names(reference), ]
    targetRegions <- targetRegions[targetRegions$end >= targetRegions$start, ]

    if (distantChimeric) {
        allChr <- names(reference)
        totalLen <- do.call("sum", lapply(reference, length))
        chrProb <- unlist(lapply(reference, length)) / totalLen
    }

    if (is.na(padding)) {
        padding <- ceiling(meanInsertLen / 2)
    }
    if (is.na(minGap)) {
        minGap <- readLen
    }

    nReadsLocal <- 0
    nReadsDistant <- 0
    nReadsNorm <- 0

    cl <- makeCluster(threads)
    registerDoParallel(cl)
    message("Using ", threads, " thread(s) for simulation.")
    for (chr in allChr[allChr %in% unique(targetRegions$chr)]) {
        chrName <- chr
        targetRegion <- targetRegions[targetRegions$chr == chrName, ]

        if (nrow(targetRegion) > 0) {
            fullSeq <- reference[[chr]]
            message("Simulating FFPE reads on the chromosome ", chr, "...")

            targetRegion$start <- unlist(lapply(targetRegion$start,
                                                function(x)
                                                    max(0, x - padding)))
            chrLen <- length(fullSeq)
            targetRegion$end <- unlist(lapply(targetRegion$end,
                                              function(x)
                                                  min(chrLen, x + padding)))

            targetRegion <- targetRegion[order(targetRegion$start),]

            if (nrow(targetRegion) > 1) {
                starts <- targetRegion$start
                starts <- starts[2:length(starts)]
                ends <- targetRegion$end
                ends <- ends[seq_len(length(ends) - 1)]
                ifMerge <- (starts - ends) <= minGap
                allStartPos <- c(targetRegion$start[1], starts[!ifMerge])
                allEndPos <- c(ends[!ifMerge],
                               targetRegion$end[nrow(targetRegion)])
                keep <- (allEndPos - allStartPos + 1) >= readLen
                allStartPos <- allStartPos[keep]
                allEndPos <- allEndPos[keep]
                startPoses <- unlist(mapply(function(x, y) seq(from = x,
                                                               to = y,
                                                               by = windowLen),
                                            allStartPos, allEndPos))
                endPoses <- mapply(function(x, y) seq(from = x,
                                                      to = y,
                                                      by = windowLen) +
                                       windowLen - 1,
                                   allStartPos, allEndPos)
                for (i in seq_along(endPoses)) {
                    endPoses[[i]][endPoses[[i]] > allEndPos[i]] <- allEndPos[i]
                }
                endPoses <- unlist(endPoses)
                endPoses[endPoses > length(fullSeq)] <- length(fullSeq)
                allStartPos <- unlist(allStartPos)
                allEndPos <- unlist(allEndPos)
            } else {
                startPoses <- targetRegion$start
                endPoses <- targetRegion$end
                allStartPos <- targetRegion$start
                allEndPos <- targetRegion$end
            }

            if (localChimeric) {
                if (pairedEnd) {
                    nSeeds <- round((endPoses - startPoses + 1) *
                                        covFFPE * localMatchRatio / readLen / 2)
                } else {
                    nSeeds <- round((endPoses - startPoses + 1) *
                                        covFFPE * localMatchRatio / readLen)
                }
                targetStarts <- startPoses - matchPadding
                targetEnds <- endPoses + matchPadding
                remove <- targetStarts < 0 | targetEnds > length(fullSeq)
                targetStarts[remove] <- startPoses[remove]
                targetEnds[remove] <- endPoses[remove]
                paddings <- rep(matchPadding, length(targetStarts))
                paddings[remove] <- 0
                startPos <- NULL
                endPos <- NULL
                targetStart <- NULL
                targetEnd <- NULL
                nSeed <- NULL
                paddingLen <- NULL
                simReads <- foreach(startPos = startPoses,
                                    endPos = endPoses,
                                    targetStart = targetStarts,
                                    targetEnd = targetEnds,
                                    nSeed = nSeeds,
                                    paddingLen = paddings,
                                    .combine = "c",
                                    .inorder = TRUE,
                                    .verbose = FALSE,
                                    .errorhandling    = "remove",
                                    .packages = c("Biostrings"),
                                    .export = c(
                                        "targetRegionalChimericReads",
                                        "targetRegionalChimericSeqs",
                                        "generateEnzymicCutSeq",
                                        "addRandomMutation",
                                        "generateReads",
                                        "addNoise",
                                        "rtruncnorm",
                                        "round")
                ) %dopar% {
                    targetRegionalChimericReads(
                        fullSeq = fullSeq,
                        startPos = startPos,
                        endPos = endPos,
                        targetStart = targetStart,
                        targetEnd = targetEnd,
                        padding = paddingLen,
                        nSeed = nSeed,
                        meanSeedLen = meanSeedLen,
                        sdSeedLen = sdSeedLen,
                        meanInsertLen = meanInsertLen,
                        sdInsertLen = sdInsertLen,
                        enzymeCut = enzymeCut,
                        readLen = readLen,
                        mutationRate = mutationRate,
                        noiseRate = noiseRate,
                        highNoiseRate = highNoiseRate,
                        highNoiseProb = highNoiseProb,
                        pairedEnd = pairedEnd,
                        revChimericProb = revChimericProb)
                }
                readsToFastq(simReads = simReads,
                             PhredScoreProfile = PhredScoreProfile,
                             prefix = prefix,
                             prefix2 = "localChimeric",
                             chr = chr,
                             pairedEnd = pairedEnd,
                             outFile = outFile,
                             threads = threads)
                message("Generated ", length(simReads), 
                        " local chimeric reads ",
                        "on chromosome ", chr, ".")
                nReadsLocal <- nReadsLocal + length(simReads)
            }

            if (distantChimeric) {
                if (pairedEnd) {
                    nSeeds <- round((allEndPos - allStartPos + 1) * covFFPE *
                                        (1 - localMatchRatio) / readLen / 2)
                } else {
                    nSeeds <- round((allEndPos - allStartPos + 1) * covFFPE *
                                        (1 - localMatchRatio) / readLen)
                }

                simReads <- foreach(startPos = allStartPos,
                                    windowLen = allEndPos - allStartPos + 1,
                                    nSeed = nSeeds,
                                    .combine = "c",
                                    .inorder = TRUE,
                                    .verbose = FALSE,
                                    .packages = "Biostrings",
                                    .errorhandling = "remove",
                                    .export = c("generateDistantChimericReads",
                                                "generateDistantChimericSeqs",
                                                "generateReads",
                                                "addRandomMutation",
                                                "addNoise",
                                                "generateTargetSeqs",
                                                "rbeta", "round",
                                                "rtruncnorm")
                ) %dopar% {
                    generateDistantChimericReads(
                        fullSeq = fullSeq,
                        reference = reference,
                        startPos = startPos,
                        windowLen = windowLen,
                        matchWinLen = matchWinLen,
                        nSeed = nSeed,
                        meanSeedLen = meanSeedLen,
                        sdSeedLen = sdSeedLen,
                        meanInsertLen = meanInsertLen,
                        sdInsertLen = sdInsertLen,
                        readLen = readLen,
                        allChr = allChr,
                        chrProb = chrProb,
                        mutationRate = mutationRate,
                        noiseRate = noiseRate,
                        highNoiseRate = highNoiseRate,
                        highNoiseProb = highNoiseProb,
                        pairedEnd = pairedEnd,
                        spikeWidth = spikeWidth,
                        revChimericProb = 0.5)
                }
                readsToFastq(simReads = simReads,
                             PhredScoreProfile = PhredScoreProfile,
                             prefix = prefix,
                             prefix2 = "distantChimeric",
                             chr = chr,
                             pairedEnd = pairedEnd,
                             outFile = outFile,
                             threads = threads)

                message("Generated ", length(simReads), 
                        " distant chimeric reads ",
                        "on chromosome ", chr, ".")
                nReadsDistant <- nReadsDistant + length(simReads)
            }

            if (normalReads) {
                if (pairedEnd) {
                    nSeeds <- round((allEndPos - allStartPos + 1) * covNorm /
                                        readLen / 2)
                } else {
                    nSeeds <- round((allEndPos - allStartPos + 1) * covNorm
                                    / readLen)
                }
                simReads <- foreach(start = allStartPos,
                                    end = allEndPos,
                                    nSeed = nSeeds,
                                    .combine = "c",
                                    .inorder = TRUE,
                                    .verbose = FALSE,
                                    .errorhandling  = "remove",
                                    .packages = c("Biostrings"),
                                    .export = c("addNoise",
                                                "rtruncnorm",
                                                "generateNormalReads")
                ) %dopar% {
                    generateNormalReads(
                        fullSeq = fullSeq[start:end],
                        nSeq = nSeed,
                        meanInsertLen = meanInsertLen,
                        sdInsertLen = sdInsertLen,
                        readLen = readLen,
                        noiseRate = noiseRate,
                        highNoiseRate = highNoiseRate,
                        highNoiseProb = highNoiseProb,
                        pairedEnd = pairedEnd)
                }

                readsToFastq(simReads = simReads,
                             PhredScoreProfile = PhredScoreProfile,
                             prefix = prefix,
                             prefix2 = "Normal",
                             chr = chr,
                             pairedEnd = pairedEnd,
                             outFile = outFile,
                             threads = threads)

                message("Generated ", length(simReads), " normal reads ",
                        "on chromosome ", chr, ".")
                nReadsNorm <- nReadsNorm + length(simReads)
            }
        }
    }
    stopCluster(cl)
    message("In totoal ", nReadsLocal, " local chimeric reads, ",
            nReadsDistant, " distant chimeric reads, ",
            nReadsNorm, " normal reads were generated.",
            "\nAlltogether ", nReadsLocal + nReadsDistant + nReadsNorm,
            " reads were generated.", "\nSimulation done.")
}
