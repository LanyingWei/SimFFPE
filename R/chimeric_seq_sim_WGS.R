generateRegionalChimericSeqs <- function(seq, nSeed, meanSeedLen, sdSeedLen,
                                         meanInsertLen, sdInsertLen,
                                         revChimericProb, enzymeCut, readLen,
                                         minSeedLen=2, maxSeedLen=20) {
    nSeed <- round(nSeed / 2)
    seqs <- list(seq, reverseComplement(seq))
    if (enzymeCut) {
        chimericSeqs <- NULL
    } else {
        chimericSeqs <- DNAStringSet()
    }
    for (seq in seqs) {
        seedLen <- as.integer(
            rtruncnorm(nSeed, minSeedLen, maxSeedLen, meanSeedLen, sdSeedLen))
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
            nearestStart <- which.min(abs(startList[['revSeqRev']] - start))
            revSeqStarts <- mapply(function(x, y) x[y],
                                   startList[['revSeq']], nearestStart)
            startList[['compSeq']] <- start(matchPDict(seeds, compSeq))

            notSameStart <- which((startList[['compSeq']] - start)!=0)
            compSeqStarts <- (mapply(function(x, y) x[which.min(y)],
                                     startList[['compSeq']], notSameStart))
            is.na(compSeqStarts) <- lengths(compSeqStarts) == 0
            compSeqStarts <- unlist(compSeqStarts)
            whichRev <- runif(length(seeds)) <= revChimericProb
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

        } else if (length(seeds) >0) {
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
                                          nSeed, meanSeedLen, sdSeedLen,
                                          meanInsertLen, sdInsertLen,enzymeCut,
                                          readLen, mutationRate, noiseRate,
                                          highNoiseRate, highNoiseProb,
                                          pairedEnd, revChimericProb) {

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
                                         meanSeedLen = meanSeedLen,
                                         sdSeedLen = sdSeedLen,
                                         meanInsertLen = meanInsertLen,
                                         sdInsertLen = sdInsertLen,
                                         readLen = readLen,
                                         enzymeCut = enzymeCut,
                                         revChimericProb = revChimericProb)

        if(length(chimericSeqs) != 0) {
            chimericSeqs <- addRandomMutation(chimericSeqs, mutationRate)
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


generateDistantChimericSeqs <- function(seq, reference, nSeed, meanSeedLen,
                                        sdSeedLen, meanInsertLen, sdInsertLen,
                                        allChr, chrProb, matchWinLen, readLen,
                                        revChimericProb, spikeWidth,
                                        minSeedLen=2, maxSeedLen=20) {
    chimericSeqs <- DNAStringSet()
    seqs <- list(seq, reverseComplement(seq))
    spikeStarts <- seq(1, length(seq), by = spikeWidth)
    rbetas <- rbeta(length(spikeStarts), 0.5, 0.5)

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

        seedLen <- as.integer(rtruncnorm(length(start), minSeedLen, maxSeedLen,
                                         meanSeedLen, sdSeedLen))
        end <- start + seedLen - 1
        start <- start[end <= length(seq)]
        end <- end[end <= length(seq)]

        seeds <- extractAt(x = complement(seq),
                           at = IRanges(start = start, end = end))

        start <- start[!grepl('N', seeds)]
        end <- end[!grepl('N', seeds)]
        seeds <- seeds[!grepl('N', seeds)]

        if (length(seeds) != 0) {

            targetSeqs <- generateTargetSeqs(reference, nSeqs = length(seeds),
                                             allChr = allChr, chrProb = chrProb,
                                             matchWinLen = matchWinLen)

            if (j == 2) {
                targetSeqs <- reverseComplement(targetSeqs)
            }
            revSeqs <- reverse(targetSeqs)
            compSeqs <- complement(targetSeqs)

            revSeqStarts <- mapply(function(x, y) start(matchPattern(x, y))[1],
                                   seeds, revSeqs)
            compSeqStarts <- mapply(function(x, y) start(matchPattern(x, y))[1],
                                    seeds, compSeqs)
            whichRev <- runif(length(seeds)) <= revChimericProb
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

generateDistantChimericReads <- function(fullSeq, reference, startPos,
                                         windowLen, matchWinLen,
                                         nSeed, meanSeedLen, sdSeedLen,
                                         meanInsertLen, sdInsertLen, readLen,
                                         allChr, chrProb, mutationRate,
                                         noiseRate, highNoiseRate,
                                         highNoiseProb, pairedEnd, spikeWidth,
                                         revChimericProb) {

    if ((startPos + windowLen - 1) <= length(fullSeq)) {
        tempSeq <- fullSeq[startPos:(startPos + windowLen - 1)]
    } else {
        tempSeq <- fullSeq[startPos:length(fullSeq)]
        nSeed <- round(length(tempSeq) / windowLen * nSeed)
    }

    if (length(tempSeq) >= readLen && nSeed >0) {
        chimericSeqs <-
            generateDistantChimericSeqs(seq = tempSeq,
                                        reference = reference,
                                        allChr = allChr,
                                        chrProb = chrProb,
                                        matchWinLen = matchWinLen,
                                        nSeed = nSeed,
                                        meanSeedLen = meanSeedLen,
                                        sdSeedLen = sdSeedLen,
                                        meanInsertLen = meanInsertLen,
                                        sdInsertLen = sdInsertLen,
                                        readLen = readLen,
                                        spikeWidth = spikeWidth,
                                        revChimericProb = revChimericProb)

        if(length(chimericSeqs)!=0) {
            chimericSeqs <-
                addRandomMutation(seqs = chimericSeqs, mutationRate)
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

generateTargetSeqs <- function(reference, nSeqs, allChr, chrProb, matchWinLen) {
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
    targetSeqsAll <- targetSeqsAll[sample(seq_along(targetSeqsAll)), ]
    return(targetSeqsAll)
}

