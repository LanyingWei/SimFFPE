generateRegionalChimericSeqs <- function(seq, nSeed, meanLogSRCRLen, 
                                         sdLogSRCRLen, maxSRCRLen, 
                                         meanInsertLen, sdInsertLen, 
                                         meanLogSRCRDist, sdLogSRCRDist, 
                                         sameStrandProp, enzymeCut, readLen,
                                         targetPadding = 0) {
    nSeed <- round(nSeed / 2)
    seqs <- list(seq, reverseComplement(seq))
    if (enzymeCut) {
        chimericSeqs <- NULL
    } else {
        chimericSeqs <- DNAStringSet()
    }
    doubleMin <- .Machine$double.xmin
    
    for (seq in seqs) {
        SRCRLen <- round(rlnorm(nSeed, meanlog = meanLogSRCRLen, 
                                sdlog = sdLogSRCRLen))
        SRCRLen <- SRCRLen[SRCRLen > 0 & SRCRLen <= maxSRCRLen]
        while(length(SRCRLen) < nSeed) {
            SRCRLen <- c(SRCRLen, round(rlnorm(nSeed, meanlog = meanLogSRCRLen, 
                                               sdlog = sdLogSRCRLen)))
            SRCRLen <- SRCRLen[SRCRLen > 0 & SRCRLen <= maxSRCRLen]
        }
        SRCRLen <- SRCRLen[seq_len(nSeed)]
        start <- sample(1:(length(seq) - targetPadding - max(SRCRLen) + 1),
                        nSeed, replace = TRUE)
        revSeq <- reverse(seq)
        compSeq <- complement(seq)

        seeds <- extractAt(x = compSeq,
                           at = IRanges(start = start, width = SRCRLen))

        start <- start[!grepl('N', seeds)]
        seeds <- seeds[!grepl('N', seeds)]

        if (length(seeds) != 0) {

            startList <- list()
            revMatch <- matchPDict(seeds, revSeq)
            startList[['revSeq']] <-  start(revMatch)
            startList[['revSeqRev']] <- length(revSeq) - end(revMatch) + 1
            dist <- abs(startList[['revSeqRev']] - start)
            select <- lapply(dist, function(x)
                if (length(x) > 1) {
                    sample(x, 1, 
                           prob = dlnorm(x, meanlog = meanLogSRCRDist, 
                                         sdlog = sdLogSRCRDist) + doubleMin)
                } else x)
            revSeqStarts <- mapply(function(x, y) x[y][1],
                                   startList[['revSeq']], which(dist == select))
            startList[['compSeq']] <- start(matchPDict(seeds, compSeq))
            dist <- abs(startList[['compSeq']] - start)
            select <- lapply(dist, function(x)
                if (length(x) > 2) {
                    sample(x[x!=0], 1, 
                           prob = dlnorm(x[x!=0], meanlog = meanLogSRCRDist, 
                                         sdlog = sdLogSRCRDist) + doubleMin)
                } else x[x!=0])
            compSeqStarts <- mapply(function(x, y) x[y][1],
                                    startList[['compSeq']], 
                                    which(dist == select))
            whichRev <- runif(length(seeds)) <= sameStrandProp
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
            sourceFragEnds <- start - 1
            targetFragStarts <- seqStarts
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



generateRegionalChimericReads <- function(seq, nSeed, meanLogSRCRLen, 
                                          sdLogSRCRLen, maxSRCRLen, 
                                          meanLogSRCRDist, sdLogSRCRDist, 
                                          meanInsertLen, sdInsertLen, enzymeCut, 
                                          readLen, chimMutRate, noiseRate, 
                                          highNoiseRate, highNoiseProp, 
                                          pairedEnd, sameStrandProp) {

    if (length(seq) >= readLen && nSeed > 1) {
        chimericSeqs <-
            generateRegionalChimericSeqs(seq = seq,
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
                                         sameStrandProp = sameStrandProp)

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


generateDistantChimericSeqs <- function(seq, reference, nSeed, 
                                        meanLogSRCRLen, sdLogSRCRLen, 
                                        maxSRCRLen, meanInsertLen, sdInsertLen,
                                        allChr, chrProb, distWinLen, readLen,
                                        sameStrandProp, spikeWidth,  betaShape1, 
                                        betaShape2, sameTarRegionProb) {
    chimericSeqs <- DNAStringSet()
    seqs <- list(seq, reverseComplement(seq))
    spikeStarts <- seq(1, length(seq), by = spikeWidth)
    rbetas <- rbeta(length(spikeStarts), betaShape1, betaShape2)
    nSeed <- nSeed / 2 / (betaShape1 / (betaShape1 + betaShape2)) 
    
    for (j in seq_along(seqs)) {
        seq <- seqs[[j]]
        if (length(seq) <= spikeWidth) {
            start <- sample(seq_along(seq), round(nSeed * rbetas[1]), 
                            replace = TRUE)
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
        
        SRCRLen <- round(rlnorm(length(start), meanlog = meanLogSRCRLen, 
                                sdlog = sdLogSRCRLen))
        SRCRLen <- SRCRLen[SRCRLen > 0 & SRCRLen <= maxSRCRLen]
        while(length(SRCRLen) < length(start)) {
            SRCRLen <- c(SRCRLen, round(rlnorm(length(start), 
                                               meanlog = meanLogSRCRLen, 
                                               sdlog = sdLogSRCRLen)))
            SRCRLen <- SRCRLen[SRCRLen > 0 & SRCRLen <= maxSRCRLen]
        }
        SRCRLen <- SRCRLen[seq_len(length(start))]
        
        end <- start + SRCRLen - 1
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
                                   distWinLen = distWinLen, 
                                   sameTarRegionProb = sameTarRegionProb)
            if (j == 2) {
                targetSeqs <- reverseComplement(targetSeqs)
            }
            revSeqs <- reverse(targetSeqs)
            compSeqs <- complement(targetSeqs)
            
            revSeqStarts <- mapply(function(x, y) {
                tmp <- start(matchPattern(x, y))
                if (length(tmp) != 0) tmp[1] else NA
            }, seeds, revSeqs)
            compSeqStarts <- mapply(function(x, y) {
                tmp <- start(matchPattern(x, y))
                if (length(tmp) != 0) tmp[1] else NA
            }, seeds, compSeqs)
            
            whichRev <- runif(length(seeds)) <= sameStrandProp
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
            sourceFragEnds <- start - 1
            targetFragStarts <- seqStarts
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
                                         windowLen, distWinLen, nSeed, 
                                         meanLogSRCRLen, sdLogSRCRLen, 
                                         maxSRCRLen, meanInsertLen, sdInsertLen, 
                                         readLen, allChr, chrProb, chimMutRate,
                                         noiseRate, highNoiseRate,
                                         highNoiseProp, pairedEnd, spikeWidth,
                                         betaShape1, betaShape2, 
                                         sameTarRegionProb, sameStrandProp) {
    reference <- readDNAStringSet(referencePath)
    reference <- reference[width(reference) >= distWinLen]
    
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

generateTargetSeqs <- function(reference, nSeqs, allChr, chrProb, 
                               distWinLen, sameTarRegionProb) {
    targetSeqsAll <- DNAStringSet()
    while (length(targetSeqsAll) < nSeqs) {
        targetChrs <- sample(allChr,
                             prob = chrProb,
                             size = nSeqs - length(targetSeqsAll),
                             replace = TRUE)
        targetChrs <- table(targetChrs)
        for (chr in names(targetChrs)) {
            startPoses <-
                sample(seq_len(length(reference[[chr]]) - distWinLen + 1),
                       targetChrs[chr], replace = TRUE)
            endPoses <- startPoses + distWinLen - 1
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

