targetRegionalChimericSeqs <- function(seq, targetSeq, padding,
                                       nSeed, meanSeedLen, sdSeedLen,
                                       meanInsertLen, sdInsertLen,
                                       readLen, revChimericProb, enzymeCut,
                                       minSeedLen=2, maxSeedLen=20) {
    nSeed <- round(nSeed / 2)
    seqs <- list(seq, reverseComplement(seq))
    if (enzymeCut) {
        chimericSeqs <- NULL
    } else {
        chimericSeqs <- DNAStringSet()
    }
    for (seq in seqs) {
        start <- sample(seq_along(seq), nSeed, replace = TRUE)
        seedLen <- as.integer(
            rtruncnorm(nSeed, minSeedLen, maxSeedLen, meanSeedLen, sdSeedLen))
        end  <- start + seedLen - 1
        start <- start[end <= length(seq)]
        seedLen <- seedLen[end <= length(seq)]

        seeds <- extractAt(x = complement(seq),
                           at = IRanges(start = start, width = seedLen))

        start <- start[!grepl('N', seeds)]
        seeds <- seeds[!grepl('N', seeds)]

        revSeq <- reverse(targetSeq)
        compSeq <- complement(targetSeq)

        if (length(seeds) != 0) {

            startList <- list()
            revMatch <- matchPDict(seeds, revSeq)
            startList[['revSeq']] <-  start(revMatch)
            startList[['revSeqRev']] <-
                length(seq) - end(revMatch) + padding + 1
            nearestStart <- which.min(abs(startList[['revSeqRev']] - start))
            revSeqStarts <- mapply(function(x, y) x[y],
                                   startList[['revSeq']], nearestStart)
            startList[['compSeq']] <- start(matchPDict(seeds, compSeq))
            notSameStart <-
                which((startList[['compSeq']] - padding - start) != 0)
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

        } else if (length(seeds) > 0) {
            mapLens <- width(seeds)
            sourceFragLens <- round(runif(length(seeds)) * (fragLens - mapLens))
            targetFragLens <- fragLens - mapLens - sourceFragLens

            sourceFragStarts <- start - sourceFragLens
            sourceFragEnds <- start + mapLens - 1
            targetFragStarts <- seqStarts + mapLens
            targetFragEnds <- seqStarts + mapLens + targetFragLens - 1
            keep <- sourceFragStarts > 0 & sourceFragEnds <= length(seq) &
                targetFragStarts > 0 & targetFragEnds <= length(targetSeq)
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
    return(chimericSeqs)
}

targetRegionalChimericReads <- function(fullSeq, startPos, endPos,
                                        targetStart, targetEnd, padding,
                                        nSeed, meanSeedLen, sdSeedLen,
                                        meanInsertLen, sdInsertLen,
                                        enzymeCut, readLen, targetPadding,
                                        mutationRate, noiseRate,
                                        highNoiseRate, highNoiseProb,
                                        pairedEnd, revChimericProb) {

    tempseq <- fullSeq[startPos:endPos]
    targetSeq <- fullSeq[targetStart:targetEnd]

    if (length(tempseq) >= readLen & nSeed > 1) {
        chimericSeqs <-
            targetRegionalChimericSeqs(seq = tempseq,
                                       targetSeq = targetSeq,
                                       padding = padding,
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
