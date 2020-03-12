generateRegionalChimericSeqs <- function(seq, nSeed, meanLogSeedLen, 
                                         sdLogSeedLen, meanInsertLen, 
                                         sdInsertLen, sdTargetDist,
                                         sameStrandProb, enzymeCut, readLen) {
    nSeed <- round(nSeed / 2)
    seqs <- list(seq, reverseComplement(seq))
    if (enzymeCut) {
        chimericSeqs <- NULL
    } else {
        chimericSeqs <- DNAStringSet()
    }
    doubleMin <- .Machine$double.xmin
    
    for (seq in seqs) {
        seedLen <- round(rlnorm(nSeed, meanlog = meanLogSeedLen, 
                                sdlog = sdLogSeedLen))
        
        start <- sample(seq_len(length(seq) - max(seedLen) + 1),
                        nSeed, replace = TRUE)
        revSeq <- reverse(seq)
        compSeq <- complement(seq)

        seeds <- extractAt(x = compSeq,
                           at = IRanges(start = start, width = seedLen))

        start <- start[!grepl('N', seeds)]
        seeds <- seeds[!grepl('N', seeds)]

        if (length(seeds) != 0) {

            startList <- list()
            revMatch <- matchPDict(seeds, revSeq)
            startList[['revSeq']] <-  start(revMatch)
            startList[['revSeqRev']] <- length(seq) - end(revMatch) + 1
            dist <- startList[['revSeqRev']] - start
            select <- lapply(dist, function(x)
                if (length(x) > 1) {
                    sample(x, 1, 
                           prob = dnorm(x, mean = 0, 
                                        sd = sdTargetDist) + doubleMin)
                } else x)
            revSeqStarts <- mapply(function(x, y) x[y][1],
                                   startList[['revSeq']], which(dist == select))
            startList[['compSeq']] <- start(matchPDict(seeds, compSeq))
            dist <- startList[['compSeq']] - start
            select <- lapply(dist, function(x)
                if (length(x) > 2) {
                    sample(x[x!=0], 1, 
                           prob = dnorm(x[x!=0], mean = 0, 
                                        sd = sdTargetDist) + doubleMin)
                } else x[x!=0])
            compSeqStarts <- mapply(function(x, y) x[y][1],
                                    startList[['compSeq']], 
                                    which(dist == select))
            whichRev <- runif(length(seeds)) <= sameStrandProb
            whichSeq <- NULL
            whichSeq[whichRev] <- "rev"
            whichSeq[!whichRev] <- "comp"
            seqStarts <- NULL
            seqStarts[whichRev] <- revSeqStarts[whichRev]
            seqStarts[!whichRev] <- compSeqStarts[!whichRev]
            noSeqStarts <- is.na(seqStarts)
            seqStarts[noSeqStarts] <- revSeqStarts[noSeqStarts]
            whichSeq[noSeqStarts & !(is.na(seqStarts))] <- "rev"
            noSeqStarts <- is.na(seqStarts)
            seqStarts[noSeqStarts] <-  compSeqStarts[noSeqStarts]
            whichSeq[noSeqStarts & !(is.na(seqStarts))] <- "comp"
            whichSeq <- whichSeq[!is.na(seqStarts)]
            seeds <- seeds[!is.na(seqStarts)]
            start <- start[!is.na(seqStarts)]
            seqStarts <- seqStarts[!is.na(seqStarts)]

            fragLens <- as.integer(
                rtruncnorm(n = length(seeds),
                           a = readLen,
                           b = Inf,
                           mean = meanInsertLen,
                           sd = sdInsertLen))
        }
        if (enzymeCut) {
            for (i in seq_along(seeds)) {
                if (whichSeq[i] == "rev") {
                    chimericSeq <-
                        generateEnzymicCutSeq(sourceSeq = seq,
                                              targetSeq = revSeq,
                                              sourceSeqMapStart = start[i],
                                              targetSeqMapStart = seqStarts[i],
                                              mapLen = width(seeds)[i],
                                              targetStrand = 'rev',
                                              fragLen = fragLens[i])
                } else if (whichSeq[i] == "comp") {
                    chimericSeq <-
                        generateEnzymicCutSeq(sourceSeq = seq,
                                              targetSeq = compSeq,
                                              sourceSeqMapStart = start[i],
                                              targetSeqMapStart = seqStarts[i],
                                              mapLen = width(seeds)[i],
                                              targetStrand = 'comp',
                                              fragLen = fragLens[i])
                }
                chimericSeqs <- c(chimericSeqs, chimericSeq)
            }

        } else if (length(seeds) > 0) {
            mapLens <- width(seeds)
            sourceFragLens <- round(runif(length(seeds)) * (fragLens - mapLens))
            targetFragLens <- fragLens - mapLens - sourceFragLens

            sourceFragStarts <- start - sourceFragLens
            sourceFragEnds <- start + mapLens - 1
            targetFragStarts <- seqStarts + mapLens
            targetFragEnds <- seqStarts + mapLens + targetFragLens - 1
            keep <- sourceFragStarts > 0 & sourceFragEnds <= length(seq) &
                targetFragStarts > 0 & targetFragEnds <= length(seq)
            revKeep <- whichSeq == "rev" & keep
            compKeep <- whichSeq == "comp" & keep
            if (any(revKeep)) {
                sourceFrags <-
                    extractAt(seq, IRanges(sourceFragStarts[revKeep],
                                           sourceFragEnds[revKeep]))
                targetFrags <-
                    complement(extractAt(revSeq,
                                         IRanges(targetFragStarts[revKeep],
                                                 targetFragEnds[revKeep])))
                chimericSeqs <- c(chimericSeqs, xscat(sourceFrags, targetFrags))
            }
            if (any(compKeep)) {
                sourceFrags <-
                    extractAt(seq, IRanges(sourceFragStarts[compKeep],
                                           sourceFragEnds[compKeep]))
                targetFrags <-
                    complement(extractAt(compSeq,
                                         IRanges(targetFragStarts[compKeep],
                                                 targetFragEnds[compKeep])))
                chimericSeqs <- c(chimericSeqs, xscat(sourceFrags, targetFrags))
            }
        }
    }
    if (enzymeCut) {
        chimericSeqs <- DNAStringSet(chimericSeqs)
    }
    chimericSeqs <- chimericSeqs[width(chimericSeqs) >= readLen]
    return (chimericSeqs)
}



generateRegionalChimericReads <- function(fullSeq, startPos, windowLen,
                                          nSeed, meanLogSeedLen, sdLogSeedLen,
                                          sdTargetDist, meanInsertLen, 
                                          sdInsertLen, enzymeCut, readLen, 
                                          chimMutRate, noiseRate, highNoiseRate, 
                                          highNoiseProb, pairedEnd, 
                                          sameStrandProb) {

    if ((startPos + windowLen - 1) <= length(fullSeq)) {
        tempSeq <- fullSeq[startPos:(startPos + windowLen - 1)]
    } else {
        tempSeq <- fullSeq[startPos:length(fullSeq)]
        nSeed <- round(length(tempSeq) / windowLen * nSeed)
    }
    if (length(tempSeq) >= readLen && nSeed > 1) {
        chimericSeqs <-
            generateRegionalChimericSeqs(seq = tempSeq,
                                         nSeed = nSeed,
                                         meanLogSeedLen = meanLogSeedLen,
                                         sdLogSeedLen = sdLogSeedLen,
                                         sdTargetDist = sdTargetDist,
                                         meanInsertLen = meanInsertLen,
                                         sdInsertLen = sdInsertLen,
                                         readLen = readLen,
                                         enzymeCut = enzymeCut,
                                         sameStrandProb = sameStrandProb)

        if(length(chimericSeqs) != 0) {
            chimericSeqs <- addRandomMutation(chimericSeqs, chimMutRate)
            simReads <- generateReads(seqs = chimericSeqs,
                                      readLen = readLen,
                                      noiseRate = noiseRate,
                                      highNoiseRate = highNoiseRate,
                                      highNoiseProb = highNoiseProb,
                                      pairedEnd = pairedEnd)
            return(simReads)
        }
    }
}


generateDistantChimericSeqs <- function(seq, reference, nSeed, 
                                        meanLogSeedLen, sdLogSeedLen, 
                                        meanInsertLen, sdInsertLen,
                                        allChr, chrProb, matchWinLen, readLen,
                                        sameStrandProb, spikeWidth,  betaShape1, 
                                        betaShape2, sameTarRegionProb) {
    chimericSeqs <- DNAStringSet()
    seqs <- list(seq, reverseComplement(seq))
    spikeStarts <- seq(1, length(seq), by = spikeWidth)
    rbetas <- rbeta(length(spikeStarts), betaShape1, betaShape2)
    nSeed <- nSeed / 2 / (betaShape1 / (betaShape1 + betaShape2)) 
    
    for (j in seq_along(seqs)) {
        seq <- seqs[[j]]
        if (length(seq) <= spikeWidth) {
            start <- sample(seq_along(seq), nSeed * rbetas[1], replace = TRUE)
        } else {
            start <- NULL
            nSpikeSeed <- nSeed / length(spikeStarts)
            for (i in seq_len(length(spikeStarts) - 1)) {
                start <- c(start, sample(spikeStarts[i]:spikeStarts[i + 1],
                                         round(rbetas[i] * nSpikeSeed),
                                         replace = TRUE))
            }
            start <- c(start,
                       sample(spikeStarts[i + 1]:(length(seq)),
                              round(rbetas[i + 1] * nSpikeSeed *
                                        (length(seq) - spikeStarts[i + 1] + 1) /
                                        spikeWidth),
                              replace = TRUE))
        }

        if (j == 2) {
            start <- length(seq) - start + 1
        }
        
        start <- sort(start)
        
        seedLen <- round(rlnorm(length(start), meanlog = meanLogSeedLen, 
                                sdlog = sdLogSeedLen))
        
        end <- start + seedLen - 1
        start <- start[end <= length(seq)]
        end <- end[end <= length(seq)]

        seeds <- extractAt(x = complement(seq),
                           at = IRanges(start = start, end = end))

        start <- start[!grepl('N', seeds)]
        end <- end[!grepl('N', seeds)]
        seeds <- seeds[!grepl('N', seeds)]

        if (length(seeds) != 0) {
            targetSeqs <- 
                generateTargetSeqs(reference, nSeqs = length(seeds),
                                   allChr = allChr, chrProb = chrProb,
                                   matchWinLen = matchWinLen, 
                                   sameTarRegionProb = sameTarRegionProb)
            if (j == 2) {
                targetSeqs <- reverseComplement(targetSeqs)
            }
            revSeqs <- reverse(targetSeqs)
            compSeqs <- complement(targetSeqs)

            revSeqStarts <- mapply(function(x, y) {
                tmp <- start(matchPattern(x, y))
                if (length(tmp) != 0) sample(tmp, 1) else NA
                }, seeds, revSeqs)
            compSeqStarts <- mapply(function(x, y) {
                tmp <- start(matchPattern(x, y))
                if (length(tmp) != 0) sample(tmp, 1) else NA
                }, seeds, compSeqs)
            
            whichRev <- runif(length(seeds)) <= sameStrandProb
            whichSeq <- NULL
            whichSeq[whichRev] <- "rev"
            whichSeq[!whichRev] <- "comp"
            seqStarts <- NULL
            seqStarts[whichRev] <- revSeqStarts[whichRev]
            seqStarts[!whichRev] <- compSeqStarts[!whichRev]
            noSeqStarts <- is.na(seqStarts)
            seqStarts[noSeqStarts] <-  revSeqStarts[noSeqStarts]
            whichSeq[noSeqStarts & !(is.na(seqStarts))] <- "rev"
            noSeqStarts <- is.na(seqStarts)
            seqStarts[noSeqStarts] <-  compSeqStarts[noSeqStarts]
            whichSeq[noSeqStarts & !(is.na(seqStarts))] <- "comp"

            keep <- !is.na(seqStarts)
            seeds <- seeds[keep]
            whichSeq <- whichSeq[keep]
            start <- start[keep]
            revSeqs <- revSeqs[keep]
            compSeqs <- compSeqs[keep]
            seqStarts <- seqStarts[keep]

            fragLens <- as.integer(rtruncnorm(n = length(seeds),
                                              a = readLen,
                                              b = Inf,
                                              mean = meanInsertLen,
                                              sd = sdInsertLen))
            mapLens <- width(seeds)
            sourceFragLens <- round(runif(length(seeds)) *
                                        (fragLens - mapLens))
            targetFragLens <- fragLens - mapLens - sourceFragLens

            sourceFragStarts <- start - sourceFragLens
            sourceFragEnds <- start + mapLens - 1
            targetFragStarts <- seqStarts + mapLens
            targetFragEnds <- seqStarts + mapLens + targetFragLens - 1
            keep <- sourceFragStarts > 0 & sourceFragEnds <= length(seq) &
                targetFragStarts > 0 & targetFragEnds <= width(revSeqs)
            revKeep <- whichSeq == "rev" & keep
            compKeep <- whichSeq == "comp" & keep
            if (any(revKeep)) {
                sourceFrags <- extractAt(seq, IRanges(sourceFragStarts[revKeep],
                                                      sourceFragEnds[revKeep]))
                ranges <- IRangesList(mapply(function(x, y) IRanges(x, y),
                                             targetFragStarts[revKeep],
                                             targetFragEnds[revKeep]))
                targetFrags <- complement(unlist(extractAt(revSeqs[revKeep],
                                                           ranges)))
                chimericSeqs <- c(chimericSeqs, xscat(sourceFrags, targetFrags))
            }
            if (any(compKeep)) {
                sourceFrags <- extractAt(seq,
                                         IRanges(sourceFragStarts[compKeep],
                                                 sourceFragEnds[compKeep]))
                ranges <- IRangesList(mapply(function(x, y) IRanges(x, y),
                                             targetFragStarts[compKeep],
                                             targetFragEnds[compKeep]))
                targetFrags <- complement(unlist(extractAt(compSeqs[compKeep],
                                                           ranges)))
                chimericSeqs <- c(chimericSeqs, xscat(sourceFrags, targetFrags))
            }
        }
    }
    chimericSeqs <- chimericSeqs[width(chimericSeqs) >= readLen]
    return(chimericSeqs)
}

generateDistantChimericReads <- function(fullSeq, referencePath, startPos,
                                         windowLen, matchWinLen,
                                         nSeed, meanLogSeedLen, sdLogSeedLen,
                                         meanInsertLen, sdInsertLen, readLen,
                                         allChr, chrProb, chimMutRate,
                                         noiseRate, highNoiseRate,
                                         highNoiseProb, pairedEnd, spikeWidth,
                                         betaShape1, betaShape2, 
                                         sameTarRegionProb, sameStrandProb) {
    reference <- readDNAStringSet(referencePath)
    reference <- reference[width(reference) >= matchWinLen]
    
    names(reference) <- unlist(lapply(names(reference),
                                      function(x) strsplit(x, " ")[[1]][1]))
    
    if ((startPos + windowLen - 1) <= length(fullSeq)) {
        tempSeq <- fullSeq[startPos:(startPos + windowLen - 1)]
    } else {
        tempSeq <- fullSeq[startPos:length(fullSeq)]
        nSeed <- round(length(tempSeq) / windowLen * nSeed)
    }

    if (length(tempSeq) >= readLen && nSeed > 0) {
        chimericSeqs <-
            generateDistantChimericSeqs(seq = tempSeq,
                                        reference = reference,
                                        allChr = allChr,
                                        chrProb = chrProb,
                                        matchWinLen = matchWinLen,
                                        nSeed = nSeed,
                                        meanLogSeedLen = meanLogSeedLen,
                                        sdLogSeedLen = sdLogSeedLen,
                                        meanInsertLen = meanInsertLen,
                                        sdInsertLen = sdInsertLen,
                                        readLen = readLen,
                                        spikeWidth = spikeWidth,
                                        betaShape1 = betaShape1, 
                                        betaShape2 = betaShape2,
                                        sameStrandProb = sameStrandProb,
                                        sameTarRegionProb = sameTarRegionProb)

        if(length(chimericSeqs) != 0) {
            chimericSeqs <-
                addRandomMutation(seqs = chimericSeqs, chimMutRate)
            simReads <- generateReads(seqs = chimericSeqs,
                                      readLen = readLen,
                                      noiseRate = noiseRate,
                                      highNoiseRate = highNoiseRate,
                                      highNoiseProb = highNoiseProb,
                                      pairedEnd = pairedEnd)
            return(simReads)
        }
    }
}

generateTargetSeqs <- function(reference, nSeqs, allChr, chrProb, 
                               matchWinLen, sameTarRegionProb) {
    targetSeqsAll <- DNAStringSet()
    while (length(targetSeqsAll) < nSeqs) {
        targetChrs <- sample(allChr,
                             prob = chrProb,
                             size = nSeqs - length(targetSeqsAll),
                             replace = TRUE)
        targetChrs <- table(targetChrs)
        for (chr in names(targetChrs)) {
            startPoses <-
                sample(seq_len(length(reference[[chr]]) - matchWinLen + 1),
                       targetChrs[chr], replace = TRUE)
            endPoses <- startPoses + matchWinLen - 1
            targetSeqs <- extractAt(x = reference[[chr]],
                                    at = IRanges(start = startPoses,
                                                 end = endPoses))
            targetSeqs <- targetSeqs[!grepl("N", targetSeqs)]
            targetSeqsAll <- c(targetSeqsAll, targetSeqs)
        }
    }
    targetSeqsAll <- targetSeqsAll[sample(seq_along(targetSeqsAll))]
    if (sameTarRegionProb > 0) {
        repTarget <- which(runif(nSeqs) <= sameTarRegionProb)
        targetSeqsAll[(repTarget + 1)] <- targetSeqsAll[(repTarget)]
        targetSeqsAll <- targetSeqsAll[seq_len(nSeqs)]
    }
    return(targetSeqsAll)
}

