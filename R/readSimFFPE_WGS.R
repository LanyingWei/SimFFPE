readSimFFPE <- function(sourceSeq, reference, PhredScoreProfile,
                        coverage, readLen, meanInsertLen, sdInsertLen, outFile,
                        enzymeCut=FALSE, covRatioFFPE=0.05, localMatchRatio=0.1,
                        windowLen=10000, matchWinLen=10000,
                        meanSeedLen=5, sdSeedLen=3,
                        revChimericProb=0.8, spikeWidth = 1500,
                        mutationRate=0.005, noiseRate=0.0015,
                        highNoiseRate=0.08, highNoiseProb=0.015,
                        pairedEnd=TRUE, prefix="SimFFPE", threads=1,
                        localChimeric=TRUE, distantChimeric=TRUE,
                        normalReads=TRUE, overWrite=FALSE) {
    if (overWrite) {
        overWriteFastq(outFile = outFile, pairedEnd = pairedEnd)
    }

    covFFPE <- coverage * covRatioFFPE / 0.85
    covNorm <- coverage - covFFPE

    if(pairedEnd) {
        nSeedRegional <-
            round(covFFPE * localMatchRatio * windowLen / readLen / 2)
    } else {
        enzymeCut <- FALSE
        nSeedRegional <- round(covFFPE * localMatchRatio * windowLen / readLen)
    }

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

    nReadsLocal <- 0
    nReadsDistant <- 0
    nReadsNorm <- 0

    cl <- makeCluster(threads)
    registerDoParallel(cl)
    message(paste("Using", threads, "thread(s) for simulation."))

    for (i in seq_along(sourceSeq)) {
        chr <- names(sourceSeq)[i]
        message(paste0("Simulating FFPE reads on the chromosome ", chr, "..."))
        fullSeq <- sourceSeq[[i]]
        if (localChimeric) {
            startPos <- NULL
            simReads <-
                foreach(startPos = seq(1, length(fullSeq),
                                       by=windowLen),
                        .combine = "c",
                        .inorder = TRUE,
                        .verbose = FALSE,
                        #.errorhandling    = "remove",
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
            message(paste0("Generated ", length(simReads),
                           " local chimeric reads ",
                           "on chromosome ", chr, "."))
            nReadsLocal <- nReadsLocal + length(simReads)
        }

        if (distantChimeric) {
            if (pairedEnd) {
                nSeedDistant <- round(covFFPE * (1 - localMatchRatio) *
                                          length(fullSeq)  / readLen / 2)
            } else {
                nSeedDistant <- round(covFFPE * (1 - localMatchRatio) *
                                          length(fullSeq)  / readLen)
            }
            nSeedDistant <- round(nSeedDistant / threads)
            distWindow <- ceiling(length(fullSeq) / threads)
            simReads <-
                foreach(startPos = seq(1, length(fullSeq),
                                       by = distWindow),
                        .combine = "c",
                        .inorder = TRUE,
                        .verbose = FALSE,
                        .packages = "Biostrings",
                        #.errorhandling = "remove",
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
                        windowLen = distWindow,
                        matchWinLen = matchWinLen,
                        nSeed = nSeedDistant,
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

            message(paste0("Generated ", length(simReads),
                           " distant chimeric reads ",
                           "on chromosome ", chr, "."))
            nReadsDistant <- nReadsDistant + length(simReads)
        }

        if (normalReads) {
            if (pairedEnd) {
                nSeqNorm <- round(covNorm * length(fullSeq) / readLen / 2)
            } else {
                nSeqNorm <- round(covNorm * length(fullSeq) / readLen)
            }
            nSeqPerCluster <- round(nSeqNorm / threads)
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
                         pairedEnd = pairedEnd,
                         outFile = outFile,
                         threads = threads)

            message(paste0("Generated ", length(simReads), " normal reads ",
                           "on chromosome ", chr, "."))
            nReadsNorm <- nReadsNorm + length(simReads)
        }
    }

    stopCluster(cl)

    message(paste("In totoal", nReadsLocal, "local chimeric reads,",
                  nReadsDistant, "distant chimeric reads",
                  nReadsNorm, "normal reads were generated.",
                  "Alltogether", nReadsLocal + nReadsDistant + nReadsNorm,
                  "reads were generated."))
    message("Simulation done.")
}
