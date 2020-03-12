readSimFFPE <- function(sourceSeq, referencePath, PhredScoreProfile, outFile,
                        coverage, readLen=150, meanInsertLen=250, 
                        sdInsertLen=80, enzymeCut=FALSE, chimericRatio=0.08, 
                        localMatchRatio=0.1, windowLen=10000, matchWinLen=10000,
                        meanLogSeedLen=1.7, sdLogSeedLen=0.4, seedPassRate=0.78, 
                        sdTargetDist=120, sameStrandProb=0.5, spikeWidth = 1500, 
                        betaShape1=0.5, betaShape2=0.5, sameTarRegionProb=0,
                        chimMutRate=0.005, noiseRate=0.0015,
                        highNoiseRate=0.08, highNoiseProb=0.015,
                        pairedEnd=TRUE, prefix="SimFFPE", threads=1,
                        localChimeric=TRUE, distantChimeric=TRUE,
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
    
    covFFPE <- coverage * chimericRatio / seedPassRate
    covNorm <- coverage * (1 - chimericRatio)
    
    if(pairedEnd) {
        nSeedRegional <-
            round(covFFPE * localMatchRatio * windowLen / readLen / 2)
    } else {
        enzymeCut <- FALSE
        nSeedRegional <- round(covFFPE * localMatchRatio * windowLen / readLen)
    }
    
    reference <- readDNAStringSet(referencePath)
    reference <- reference[width(reference) >= matchWinLen]
    
    names(reference) <- unlist(lapply(names(reference),
                                      function(x) strsplit(x, " ")[[1]][1]))
    names(sourceSeq) <- unlist(lapply(names(sourceSeq),
                                      function(x) strsplit(x, " ")[[1]][1]))
    
    if (distantChimeric) {
        allChr <- names(reference)
        totalLen <- do.call("sum", lapply(reference, length))
        chrProb <- unlist(lapply(reference, length)) / totalLen
    }
    reference <- NULL
    
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
        if (localChimeric) {
            startPos <- NULL
            simReads <-
                foreach(startPos = seq(1, length(fullSeq),
                                       by=windowLen),
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
                        fullSeq = fullSeq,
                        startPos = startPos,
                        windowLen = windowLen,
                        nSeed = nSeedRegional,
                        meanLogSeedLen = meanLogSeedLen,
                        sdLogSeedLen = sdLogSeedLen,
                        sdTargetDist = sdTargetDist,
                        meanInsertLen = meanInsertLen,
                        sdInsertLen = sdInsertLen,
                        enzymeCut = enzymeCut,
                        readLen = readLen,
                        chimMutRate = chimMutRate,
                        noiseRate = noiseRate,
                        highNoiseRate = highNoiseRate,
                        highNoiseProb = highNoiseProb,
                        pairedEnd = pairedEnd,
                        sameStrandProb = sameStrandProb)
                }
            
            readsToFastq(simReads = simReads,
                         PhredScoreProfile = PhredScoreProfile,
                         prefix = prefix,
                         prefix2 = "localChimeric",
                         chr = chr,
                         pairedEnd = pairedEnd,
                         outFile = outFile,
                         threads = threads)
            message("Generated ", length(simReads), " local chimeric reads ",
                    "on chromosome ", chr, ".")
            nReadsLocal <- nReadsLocal + length(simReads)
            simReads <- NULL
        }
        
        if (distantChimeric) {
            if (pairedEnd) {
                nSeedDistant <- round(covFFPE * (1 - localMatchRatio) *
                                          length(fullSeq) / readLen / 2)
            } else {
                nSeedDistant <- round(covFFPE * (1 - localMatchRatio) *
                                          length(fullSeq) / readLen)
            }
            
            if (nSeedDistant > 1e+6) {
                distWindow <- round(length(fullSeq) / (nSeedDistant / 1e+6))
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
                            matchWinLen = matchWinLen,
                            nSeed = nSeedDistant,
                            meanLogSeedLen = meanLogSeedLen,
                            sdLogSeedLen = sdLogSeedLen,
                            meanInsertLen = meanInsertLen,
                            sdInsertLen = sdInsertLen,
                            readLen = readLen,
                            allChr = allChr,
                            chrProb = chrProb,
                            chimMutRate = chimMutRate,
                            noiseRate = noiseRate,
                            highNoiseRate = highNoiseRate,
                            highNoiseProb = highNoiseProb,
                            pairedEnd = pairedEnd,
                            spikeWidth = spikeWidth,
                            betaShape1 = betaShape1,
                            betaShape2 = betaShape2,
                            sameTarRegionProb = sameTarRegionProb,
                            sameStrandProb= 0.5)
                    }
                
                readsToFastq(simReads = simReads,
                             PhredScoreProfile = PhredScoreProfile,
                             prefix = prefix,
                             prefix2 = "distantChimeric",
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
                        highNoiseProb = highNoiseProb,
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
    
    message("In totoal ", nReadsLocal, " local chimeric reads, ",
            nReadsDistant, " distant chimeric reads, ",
            nReadsNorm, " normal reads were generated.",
            "\nAlltogether ", nReadsLocal + nReadsDistant + nReadsNorm,
            " reads were generated.", "\nSimulation done.")
}
