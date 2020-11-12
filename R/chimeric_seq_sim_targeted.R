targetRegionalChimericReads <- function(targetSeq, targetPadding,
                                        nSeed, meanLogSRCRLen, sdLogSRCRLen,
                                        maxSRCRLen, meanLogSRCRDist, 
                                        sdLogSRCRDist, meanInsertLen, 
                                        sdInsertLen, enzymeCut, readLen, 
                                        chimMutRate, noiseRate,
                                        highNoiseRate, highNoiseProp,
                                        pairedEnd, sameStrandProp) {

    if (length(targetSeq) >= readLen & nSeed > 1) {
        chimericSeqs <-
            generateRegionalChimericSeqs(seq = targetSeq,
                                         nSeed = nSeed,
                                         meanLogSRCRLen = meanLogSRCRLen,
                                         sdLogSRCRLen = sdLogSRCRLen,
                                         maxSRCRLen = maxSRCRLen,
                                         meanLogSRCRDist = meanLogSRCRDist,
                                         sdLogSRCRDist = sdLogSRCRDist,
                                         meanInsertLen = meanInsertLen,
                                         sdInsertLen = sdInsertLen,
                                         readLen = readLen,
                                         enzymeCut = enzymeCut,
                                         sameStrandProp = sameStrandProp, 
                                         targetPadding = targetPadding)

        if(length(chimericSeqs) != 0) {
            chimericSeqs <- addRandomMutation(chimericSeqs, chimMutRate)
            simReads <- generateReads(seqs = chimericSeqs,
                                      readLen = readLen,
                                      noiseRate = noiseRate,
                                      highNoiseRate = highNoiseRate,
                                      highNoiseProp = highNoiseProp,
                                      pairedEnd = pairedEnd)
            return(simReads)
        }
    }
}

targetDistantChimericReads <- function(fullSeq, referencePath, startPoses,
                                       endPoses, distWinLen, nSeeds,
                                       meanLogSRCRLen, sdLogSRCRLen, maxSRCRLen, 
                                       meanInsertLen, sdInsertLen, readLen, 
                                       allChr, chrProb, chimMutRate, 
                                       noiseRate, highNoiseRate, highNoiseProp,
                                       pairedEnd, spikeWidth, betaShape1, 
                                       betaShape2, sameTarRegionProb, 
                                       sameStrandProp) {
    reference <- readDNAStringSet(referencePath)
    reference <- reference[width(reference) >= distWinLen]
    
    names(reference) <- unlist(lapply(names(reference),
                                      function(x) strsplit(x, " ")[[1]][1]))
    res <- NULL
    for (i in seq_along(startPoses)) {
        startPos <- startPoses[i]
        endPos <- endPoses[i]
        nSeed <- nSeeds[i]
        simulatedReads <- 
            generateTarDistChimericReads(seq = fullSeq[startPos:endPos], 
                                         reference =  reference, 
                                         distWinLen = distWinLen, 
                                         nSeed = nSeed, 
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
        res <- c(res, simulatedReads)
    }
    return(res)
}

generateTarDistChimericReads <- function(seq, reference, distWinLen, nSeed, 
                                         meanLogSRCRLen, sdLogSRCRLen, 
                                         maxSRCRLen, meanInsertLen, sdInsertLen, 
                                         readLen, allChr, chrProb, chimMutRate,
                                         noiseRate, highNoiseRate,
                                         highNoiseProp, pairedEnd, spikeWidth,
                                         betaShape1, betaShape2, 
                                         sameTarRegionProb, sameStrandProp) {
    
    if (length(seq) >= readLen && nSeed > 0) {
        chimericSeqs <-
            generateDistantChimericSeqs(seq = seq,
                                        reference = reference,
                                        allChr = allChr,
                                        chrProb = chrProb,
                                        distWinLen = distWinLen,
                                        nSeed = nSeed,
                                        meanLogSRCRLen = meanLogSRCRLen,
                                        sdLogSRCRLen = sdLogSRCRLen,
                                        maxSRCRLen = maxSRCRLen,
                                        meanInsertLen = meanInsertLen,
                                        sdInsertLen = sdInsertLen,
                                        readLen = readLen,
                                        spikeWidth = spikeWidth,
                                        betaShape1 = betaShape1, 
                                        betaShape2 = betaShape2,
                                        sameStrandProp = sameStrandProp,
                                        sameTarRegionProb = sameTarRegionProb)
        
        if(length(chimericSeqs) != 0) {
            chimericSeqs <-
                addRandomMutation(seqs = chimericSeqs, chimMutRate)
            simReads <- generateReads(seqs = chimericSeqs,
                                      readLen = readLen,
                                      noiseRate = noiseRate,
                                      highNoiseRate = highNoiseRate,
                                      highNoiseProp = highNoiseProp,
                                      pairedEnd = pairedEnd)
            return(simReads)
        }
    }
}

