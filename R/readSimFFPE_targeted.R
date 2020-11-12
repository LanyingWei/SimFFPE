targetReadSimFFPE <- function(referencePath, PhredScoreProfile,
                              targetRegions, outFile, coverage,
                              readLen=150, meanInsertLen=250, sdInsertLen=80,
                              enzymeCut=FALSE, padding=50, minGap=5, 
                              chimericProp=0.1, sameChrProp=0.43, 
                              windowLen=5000, adjChimProp=0.63, 
                              sameStrandProp=0.65, meanLogSRCRLen=1.8, 
                              sdLogSRCRLen=0.55, maxSRCRLen=32,
                              meanLogSRCRDist=4.7, sdLogSRCRDist=0.35, 
                              distWinLen=5000, spikeWidth=1500, 
                              betaShape1=0.5, betaShape2=0.5, 
                              sameTarRegionProb=0, adjFactor = 1.3, 
                              distFactor = 1.3, chimMutRate=0.003, 
                              noiseRate=0.0015, highNoiseRate=0.08, 
                              highNoiseProp=0.01, pairedEnd=TRUE, 
                              prefix="SimFFPE", threads=1,
                              adjChimeric=TRUE, distChimeric=TRUE,
                              normalReads=TRUE, overWrite=FALSE) {
    
    if(missing(referencePath)) stop("referencePath is required")
    if(missing(PhredScoreProfile)) stop("PhredScoreProfile is required")
    if(missing(targetRegions)) stop("targetRegions is required")
    if(missing(outFile)) stop("outFile is required")
    if(missing(coverage)) stop("coverage is required")
    
    if(max(betaShape1, betaShape2) > 1) {
        stop ("betaShape1 and betaShape2 should not be greater than 1")
    }
    if(min(betaShape1, betaShape2) <= 0) {
        stop ("betaShape1 and betaShape2 should be greater than 0")
    }
    
    if (nrow(PhredScoreProfile) != readLen) {
        stop("The number of rows in PhredScoreProfile should be the same as", 
             " the input readLen")
    }
    
    if (overWrite) {
        overWriteFastq(outFile = outFile, pairedEnd = pairedEnd)
    }
    
    covFFPE <- coverage * chimericProp
    covNorm <- coverage * (1 - chimericProp)
    adjChimProp <- sameChrProp * adjChimProp
    sameStrandProp <- sameStrandProp * 1.05
    
    if(!pairedEnd) {
        enzymeCut <- FALSE
    }
    
    reference <- readDNAStringSet(referencePath)
    reference <- reference[width(reference) >= distWinLen]
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
    
    
    nReadsLocal <- 0
    nReadsDistant <- 0
    nReadsNorm <- 0
    
    cl <- makeCluster(threads)
    registerDoParallel(cl)
    message("Using ", threads, " thread(s) for simulation.")
    for (chr in unique(targetRegions$chr)) {
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
            
            if (adjChimeric) {
                
                if (pairedEnd) {
                    nSeeds <- round((endPoses - startPoses + 1)  * adjFactor * 
                                        covFFPE * adjChimProp / readLen / 2)
                } else {
                    nSeeds <- round((endPoses - startPoses + 1) * adjFactor * 
                                        covFFPE * adjChimProp / readLen)
                }
                matchPaddings <- 
                    round((windowLen - (endPoses - startPoses + 1)) / 2)
                matchPaddings[matchPaddings < 0] <- 0
                targetStarts <- startPoses - matchPaddings
                targetEnds <- endPoses + matchPaddings
                remove <- targetStarts < 0 | targetEnds > length(fullSeq)
                targetStarts[remove] <- startPoses[remove]
                targetEnds[remove] <- endPoses[remove]
                matchPaddings[remove] <- 0
                startPos <- NULL
                endPos <- NULL
                targetStart <- NULL
                targetEnd <- NULL
                nSeed <- NULL
                paddingLen <- NULL
                simReads <- foreach(targetStart = targetStarts,
                                    targetEnd = targetEnds,
                                    nSeed = nSeeds,
                                    paddingLen = matchPaddings,
                                    .combine = "c",
                                    .inorder = TRUE,
                                    .verbose = FALSE,
                                    .errorhandling    = "remove",
                                    .packages = c("Biostrings"),
                                    .export = c(
                                        "targetRegionalChimericReads",
                                        "generateRegionalChimericSeqs",
                                        "generateEnzymicCutSeq",
                                        "addRandomMutation",
                                        "generateReads",
                                        "addNoise",
                                        "rtruncnorm",
                                        "round")
                ) %dopar% {
                    targetRegionalChimericReads(
                        targetSeq = fullSeq[targetStart:targetEnd],
                        targetPadding = paddingLen,
                        nSeed = nSeed,
                        meanLogSRCRLen = meanLogSRCRLen,
                        sdLogSRCRLen = sdLogSRCRLen,
                        maxSRCRLen = maxSRCRLen,
                        meanLogSRCRDist = meanLogSRCRDist,
                        sdLogSRCRDist = sdLogSRCRDist,
                        meanInsertLen = meanInsertLen,
                        sdInsertLen = sdInsertLen,
                        enzymeCut = enzymeCut,
                        readLen = readLen,
                        chimMutRate = chimMutRate,
                        noiseRate = noiseRate,
                        highNoiseRate = highNoiseRate,
                        highNoiseProp = highNoiseProp,
                        pairedEnd = pairedEnd,
                        sameStrandProp = sameStrandProp)
                }
                readsToFastq(simReads = simReads,
                             PhredScoreProfile = PhredScoreProfile,
                             prefix = prefix,
                             prefix2 = "adjChimeric",
                             chr = chr,
                             pairedEnd = pairedEnd,
                             outFile = outFile,
                             threads = threads)
                message("Generated ", length(simReads), 
                        " adjacent chimeric reads ",
                        "on chromosome ", chr, ".")
                nReadsLocal <- nReadsLocal + length(simReads)
                simReads <- NULL
            }
            
            if (distChimeric) {
                
                allChr <- names(reference)
                totalLen <- do.call("sum", lapply(reference, length))
                chrProb <- unlist(lapply(reference, length)) / totalLen
                tmpProb <- (sameChrProp - adjChimProp)/(1 - adjChimProp)
                chrProb[chr] <- 0
                chrProb <- chrProb / sum(chrProb) * (1 - tmpProb)
                chrProb[chr] <- tmpProb
                
                if (pairedEnd) {
                    allSeeds <- round((allEndPos - allStartPos + 1) * 
                                          covFFPE * distFactor / 
                                          readLen / 2)
                } else {
                    allSeeds <- round((allEndPos - allStartPos + 1) *
                                          covFFPE * distFactor / readLen)
                }
                
                if (length(allStartPos) > 1000) {
                    nTimes <- ceiling(length(allStartPos) / 1000 / threads)
                } else {
                    nTimes <- 1
                }
                
                nReadsDistantChr <- 0
                firstID <- 1
                
                for (j in seq_len(nTimes)) {
                    startPos <- split(allStartPos, 
                                      ceiling(seq_along(allStartPos) / 
                                                  1000 / threads))[[j]]
                    endPos <- split(allEndPos, 
                                    ceiling(seq_along(allEndPos) / 
                                                1000 / threads))[[j]]
                    
                    nSeeds <- split(allSeeds, 
                                    ceiling(seq_along(allSeeds) / 
                                                1000 / threads))[[j]]
                    tmp <- ceiling(length(startPos) / threads)
                    startPos <- split(startPos, 
                                      ceiling(seq_along(startPos) / tmp))
                    endPos <- split(endPos, 
                                    ceiling(seq_along(endPos) / tmp))
                    nSeeds <- split(nSeeds, 
                                    ceiling(seq_along(nSeeds) / tmp))
                    
                    simReads <- 
                        foreach(i = seq_along(startPos),
                                .combine = "c",
                                .inorder = TRUE,
                                .verbose = FALSE,
                                .packages = "Biostrings",
                                .errorhandling = "remove",
                                .export = c("targetDistantChimericReads",
                                            "generateTarDistChimericReads",
                                            "generateDistantChimericSeqs",
                                            "generateReads",
                                            "addRandomMutation",
                                            "addNoise",
                                            "generateTargetSeqs",
                                            "rbeta", "round",
                                            "rtruncnorm")
                        ) %dopar% {
                            targetDistantChimericReads(
                                fullSeq = fullSeq,
                                referencePath = referencePath,
                                startPoses = startPos[[i]],
                                endPoses = endPos[[i]],
                                distWinLen = distWinLen,
                                nSeeds = nSeeds[[i]],
                                meanLogSRCRLen = meanLogSRCRLen,
                                sdLogSRCRLen = sdLogSRCRLen,
                                maxSRCRLen = maxSRCRLen,
                                meanInsertLen = meanInsertLen,
                                sdInsertLen = sdInsertLen,
                                readLen = readLen,
                                allChr = allChr,
                                chrProb = chrProb,
                                chimMutRate = chimMutRate,
                                noiseRate = noiseRate,
                                highNoiseRate = highNoiseRate,
                                highNoiseProp = highNoiseProp,
                                pairedEnd = pairedEnd,
                                spikeWidth = spikeWidth,
                                betaShape1 = betaShape1,
                                betaShape2 = betaShape2,
                                sameTarRegionProb = sameTarRegionProb,
                                sameStrandProp = 0.5)
                        }
                    readsToFastq(simReads = simReads,
                                 PhredScoreProfile = PhredScoreProfile,
                                 prefix = prefix,
                                 prefix2 = "distChimeric",
                                 chr = chr,
                                 firstID = firstID,
                                 pairedEnd = pairedEnd,
                                 outFile = outFile,
                                 threads = threads)
                    if (pairedEnd) {
                        firstID = firstID + length(simReads) / 2
                    } else {
                        firstID = firstID + length(simReads)
                    }
                    nReadsDistantChr <- nReadsDistantChr + length(simReads)
                    simReads <- NULL
                }
                message("Generated ", nReadsDistantChr,
                        " distant chimeric reads ",
                        "on chromosome ", chr, ".")
                nReadsDistant <- nReadsDistant + nReadsDistantChr
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
                        highNoiseProp = highNoiseProp,
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
                simReads <- NULL
            }
        }
    }
    stopCluster(cl)
    message("In totoal ", nReadsLocal, " adjacent chimeric reads, ",
            nReadsDistant, " distant chimeric reads, ",
            nReadsNorm, " normal reads were generated.",
            "\nAlltogether ", nReadsLocal + nReadsDistant + nReadsNorm,
            " reads were generated. ", nReadsLocal + nReadsDistant, 
            " reads (", round((nReadsLocal + nReadsDistant) / 
                                  (nReadsLocal + nReadsDistant + nReadsNorm), 
                              4) * 100, 
            "%, adjusted by chimericProp) are artifact chimeric reads.",
            "\nOf all artifact chimeric reads, ", nReadsLocal, " reads (",
            round(nReadsLocal / (nReadsLocal + nReadsDistant), 4) * 100,
            "%, adjusted by sameChrProp * adjChimProp)",
            " are adjacent chimeric reads.",
            "\nSimulation done.")
}

