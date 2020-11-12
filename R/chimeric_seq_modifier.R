generateEnzymicCutSeq <- function(sourceSeq, targetSeq, sourceSeqMapStart,
                                  targetSeqMapStart, mapLen, fragLen,
                                  targetStrand=NA) {
    if (targetStrand == 'comp' &&
        sourceSeqMapStart - targetSeqMapStart >= 0.8 * fragLen &&
        sourceSeqMapStart - targetSeqMapStart <= 1.2 * fragLen) {
        sourceFragLen <-
            round(runif(1) * (sourceSeqMapStart - targetSeqMapStart))
        targetFragLen <-
            sourceSeqMapStart - targetSeqMapStart - mapLen - sourceFragLen
    } else if (targetStrand == 'rev' &&
               abs(length(targetSeq) - targetSeqMapStart +
                   2 - sourceSeqMapStart) <= fragLen) {
        startDif <-
            length(targetSeq) - targetSeqMapStart +  2 -
            sourceSeqMapStart - mapLen
        fixLen <- abs(startDif)
        restLen <- fragLen - fixLen - mapLen
        shortLen <- round(restLen / 2)
        longLen <- fixLen + shortLen
        if (startDif > 0) {
            sourceFragLen <- shortLen
            targetFragLen <- longLen
        } else {
            sourceFragLen <- longLen
            targetFragLen <- shortLen
        }
    } else {
        sourceFragLen <- round(runif(1) * (fragLen - mapLen))
        targetFragLen <- fragLen - mapLen - sourceFragLen
    }
    
    sourceFragStart <- sourceSeqMapStart - sourceFragLen
    sourceFragEnd <- sourceSeqMapStart + mapLen - 1
    targetFragStart <- targetSeqMapStart + mapLen
    targetFragEnd <- targetSeqMapStart + mapLen + targetFragLen - 1
    
    chimericSeq <- NULL
    if (sourceFragStart > 0 && sourceFragEnd <= length(sourceSeq) &&
        targetFragStart > 0 && targetFragEnd <= length(targetSeq)) {
        sourceFrag <- sourceSeq[sourceFragStart:sourceFragEnd]
        targetFrag <- complement(targetSeq[targetFragStart:targetFragEnd])
        chimericSeq <- xscat(sourceFrag, targetFrag)
    }
    return (chimericSeq)
}

addRandomMutation <- function(seqs, mutationRate) {
    mutationRate <- mutationRate / 0.75
    width <- width(seqs)
    maxWidth <- max(width)
    seqs <- padAndClip(seqs,
                       IRanges(1, maxWidth),
                       Rpadding.letter = ".",
                       Lpadding.letter = ".")
    noiseIndex <- matrix(nrow = length(seqs), ncol = maxWidth)
    noises <- character(length = length(seqs))
    for (i in seq_along(seqs)) {
        noiseIndex[i, ] <- runif(maxWidth) <= mutationRate
        noises[i] <- paste0(sample(c("A", "C", "G", "T"),
                                   size = sum(noiseIndex[i, ]),
                                   replace = TRUE),
                            collapse = "")
    }

    seqs <- replaceLetterAt(x = seqs,
                            at = noiseIndex,
                            letter =  noises)

    seqs <- padAndClip(seqs, IRanges(1, width = width),
                       Lpadding.letter = ".", Rpadding.letter = ".")
    return(seqs)
}
