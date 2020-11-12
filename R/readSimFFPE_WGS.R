readSimFFPE <- function(sourceSeq, referencePath, PhredScoreProfile, outFile,
                        coverage, readLen=150, meanInsertLen=250, 
                        sdInsertLen=80, enzymeCut=FALSE, chimericProp=0.1, 
                        sameChrProp=0.43, windowLen=5000, adjChimProp=0.63, 
                        sameStrandProp=0.65, meanLogSRCRLen=1.8, 
                        sdLogSRCRLen=0.55, maxSRCRLen=32,
                        meanLogSRCRDist=4.7, sdLogSRCRDist=0.35, 
                        distWinLen=5000, spikeWidth = 1500, 
                        betaShape1=0.5, betaShape2=0.5, sameTarRegionProb=0,
                        adjFactor=1.65, distFactor=1.65, 
                        chimMutRate=0.003, noiseRate=0.0015,
                        highNoiseRate=0.08, highNoiseProp=0.01,
                        pairedEnd=TRUE, prefix="SimFFPE", threads=1,
                        adjChimeric=TRUE, distChimeric=TRUE,
                        normalReads=TRUE, overWrite=FALSE) {
    
    if(missing(sourceSeq)) stop("sourceSeq is required")
    if(missing(referencePath)) stop("referencePath is required")
    if(missing(PhredScoreProfile)) stop("PhredScoreProfile is required")
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
    
    reference <- readDNAStringSet(referencePath)
    reference <- reference[width(reference) >= distWinLen]
    
    names(reference) <- unlist(lapply(names(reference),
                                      function(x) strsplit(x, " ")[[1]][1]))
    names(sourceSeq) <- unlist(lapply(names(sourceSeq),
                                      function(x) strsplit(x, " ")[[1]][1]))
    nReadsLocal <- 0
    nReadsDistant <- 0
    nReadsNorm <- 0
    
    cl <- makeCluster(threads)
    registerDoParallel(cl)
    message("Using ", threads, " thread(s) for simulation.")
    
    for (i in seq_along(sourceSeq)) {
        chr <- names(sourceSeq)[i]
        message("Simulating FFPE reads on the chromosome ", chr, "...")
        fullSeq <- sourceSeq[[i]]
        
        if (adjChimeric) {
            startPoses <- seq(1, length(fullSeq), by=windowLen)
            endPoses <- startPoses + windowLen - 1
            endPoses[length(endPoses)] <- length(fullSeq)
            if(pairedEnd) {
                nSeedsAll <- round(covFFPE * adjChimProp  * adjFactor * 
                                       length(fullSeq) / readLen / 2)
            } else {
                enzymeCut <- FALSE
                nSeedsAll <- round(covFFPE * adjChimProp  * adjFactor * 
                                       length(fullSeq) / readLen)
            }
            
            nSeedsRegional <- floor(nSeedsAll / length(startPoses))
            nSeeds <- rep(nSeedsRegional, length(startPoses))
            restSeeds <- nSeedsAll - sum(nSeeds)
            nSeeds[length(startPoses)] <- 
                round(nSeedsRegional / windowLen * 
                          (length(fullSeq) - startPoses[length(startPoses)]+1))
            tmp <- sample(seq_along(startPoses), restSeeds)
            nSeeds[tmp] <- nSeeds[tmp] + 1

            startPos <- NULL
            endPos <- NULL
            nSeed <- NULL

            simReads <-
                foreach(startPos = startPoses,
                        endPos = endPoses,
                        nSeed = nSeeds,
                        .combine = "c",
                        .inorder = TRUE,
                        .verbose = FALSE,
                        .errorhandling    = "remove",
                        .packages = c("Biostrings"),
                        .export = c("generateRegionalChimericReads",
                                    "generateRegionalChimericSeqs",
                                    "generateEnzymicCutSeq",
                                    "addRandomMutation",
                                    "addNoise",
                                    "generateReads",
                                    "rtruncnorm")
                ) %dopar% {
                    generateRegionalChimericReads(
                        seq = fullSeq[startPos:endPos],
                        nSeed = nSeed,
                        meanLogSRCRLen = meanLogSRCRLen,
                        sdLogSRCRLen = sdLogSRCRLen,
                        maxSRCRLen=maxSRCRLen,
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
            message("Generated ", length(simReads), " adjacent chimeric reads ",
                    "on chromosome ", chr, ".")
            nReadsLocal <- nReadsLocal + length(simReads)
            simReads <- NULL
        }
        
        if (distChimeric) {
            
            allChr <- names(reference)
            totalLen <- do.call("sum", lapply(reference, length))
            chrProb <- unlist(lapply(reference, length)) / totalLen
            tmpProb <- (sameChrProp - adjChimProp) / (1 - adjChimProp)
            chrProb[chr] <- 0
            chrProb <- chrProb / sum(chrProb) * (1 - tmpProb)
            chrProb[chr] <- tmpProb
            
            if (pairedEnd) {
                nSeedDistant <- round(covFFPE * (1 - adjChimProp) * distFactor *
                                          length(fullSeq) / readLen / 2)
            } else {
                nSeedDistant <- round(covFFPE * (1 - adjChimProp) * distFactor *
                                          length(fullSeq) / readLen)
            }
            
            
            if (nSeedDistant > 1e+4) {
                distWindow <- round(length(fullSeq) / (nSeedDistant / 1e+4))
                nTimes <- ceiling(length(fullSeq) / distWindow / threads)
                nSeedDistant <- 
                    round(nSeedDistant *  distWindow / length(fullSeq))
            } else {
                distWindow <- length(fullSeq)
                nTimes <- 1
            }
            
            allStartPos <- seq(1, length(fullSeq), by = distWindow)
            nReadsDistantChr <- 0
            firstID <- 1
            for (j in seq_len(nTimes)) {
                tmp <- min(j * threads, length(allStartPos))
                startPos <- allStartPos[((j - 1) * threads + 1):tmp]
                simReads <-
                    foreach(startPos = startPos,
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
                            referencePath = referencePath,
                            startPos = startPos,
                            windowLen = distWindow,
                            distWinLen = distWinLen,
                            nSeed = nSeedDistant,
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
                            sameStrandProp= 0.5)
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
                nSeqNorm <- round(covNorm * length(fullSeq) / readLen / 2)
            } else {
                nSeqNorm <- round(covNorm * length(fullSeq) / readLen)
            }
            if (nSeqNorm > 5e+6) {
                nSeqPerTime <- c(rep(5e+6, floor(nSeqNorm / 5e+6)), 
                                 nSeqNorm %% 5e+6)
                nTimes <- ceiling(nSeqNorm / 5e+6)
            } else {
                nSeqPerTime <- nSeqNorm
                nTimes <- 1
            }
            nReadsNormChr <- 0
            firstID <- 1
            for (j in seq_len(nTimes)) {
                
                nSeqPerCluster <- round(nSeqPerTime[j] / threads)
                
                simReads <- foreach(i = seq_len(threads),
                                    .combine = "c",
                                    .inorder = TRUE,
                                    .verbose = FALSE,
                                    #.errorhandling  = "remove",
                                    .packages = c("Biostrings"),
                                    .export = c("addNoise",
                                                "rtruncnorm",
                                                "generateNormalReads")
                ) %dopar% {
                    generateNormalReads(
                        fullSeq = fullSeq,
                        nSeq = nSeqPerCluster,
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
                             firstID = firstID,
                             pairedEnd = pairedEnd,
                             outFile = outFile,
                             threads = threads)
                if (pairedEnd) {
                    firstID = firstID + length(simReads) / 2
                } else {
                    firstID = firstID + length(simReads)
                }
                nReadsNormChr <- nReadsNormChr + length(simReads)
                simReads <- NULL
            }
            message("Generated ", nReadsNormChr, " normal reads ",
                    "on chromosome ", chr, ".")
            nReadsNorm <- nReadsNorm + nReadsNormChr
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
