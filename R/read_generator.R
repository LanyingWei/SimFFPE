generateReads <- function(seqs, readLen, pairedEnd, noiseRate,
                          highNoiseRate, highNoiseProb) {
    noiseRate <- noiseRate / 0.75
    highNoiseRate <- highNoiseRate / 0.75
    revComp <- reverseComplement(seqs)
    nSeq <- length(seqs)
    if (pairedEnd) {
        ranges <- IRangesList(rep(list(IRanges(start = 1, width = readLen)),
                                  nSeq))
        reads1 <- unlist(extractAt(x = seqs, at = ranges))
        reads2 <- unlist(extractAt(x = revComp, at = ranges))
        removedReads <- union(which(grepl("N", reads1)),
                              which(grepl("N", reads2)))
        if (length(removedReads) > 0) {
            reads <- c(reads1[-removedReads], reads2[-removedReads])
        } else {
            reads <- c(reads1, reads2)
        }
        nReads <- length(reads)
        if (nReads > 0) {
            noiseRates <- rep(noiseRate, length(reads))
            noiseRates[sample(seq_len(nReads),
                              size = nReads * highNoiseProb,
                              replace = TRUE)] <- highNoiseRate
            reads <- addNoise(reads, noiseRates=noiseRates, readLen = readLen)
            reads <- as.character(reads)
            mid <- nReads / 2
            mid2 <- round(mid / 2)
            reads1 <- c(rbind(reads[seq_len(mid2)],
                              reads[(mid + 1):(mid + mid2)]))
            reads2 <- c(rbind(reads[(mid + mid2 + 1):nReads],
                              reads[(mid2 + 1):mid]))
            reads <- c(reads1, reads2)
            return(reads)
        }

    } else {
        mid <- round(nSeq / 2)
        ranges <- IRangesList(rep(list(IRanges(start = 1, width = readLen)),
                                  mid))
        reads1 <- unlist(extractAt(x = seqs, at = ranges))
        reads2 <- unlist(extractAt(x = revComp, at = ranges))
        reads <- c(reads1, reads2)

        reads <- reads[!grepl("N", reads)]
        nReads <- length(reads)
        if (nReads > 0) {
            noiseRates <- rep(noiseRate, length(reads))
            noiseRates[sample(seq_len(nReads),
                              size = nReads * highNoiseProb,
                              replace = TRUE)] <- highNoiseRate
            reads <- addNoise(reads, noiseRates=noiseRates, readLen = readLen)
            reads <- as.character(reads)
            return(reads)
        }
    }
}

addNoise <- function (reads, noiseRates, readLen) {
    noiseIndex <- matrix(nrow = length(reads), ncol = readLen)
    noises <- character(length = length(reads))
    for (i in seq_along(reads)) {
        noiseIndex[i, ] <- runif(readLen) <= noiseRates[i]
        noises[i] <- paste0(sample(c("A", "C", "G", "T"),
                                   size = sum(noiseIndex[i, ]),
                                   replace = TRUE), collapse = "")
    }


    reads <- replaceLetterAt(x = reads,
                             at = noiseIndex,
                             letter =  noises)

    return(reads)
}


generateNormalReads <- function(fullSeq, nSeq, meanInsertLen, sdInsertLen,
                                readLen, noiseRate, highNoiseRate,
                                highNoiseProb, pairedEnd) {


    starts <- sample(seq_len(length(fullSeq) - meanInsertLen + 1),
                     size = nSeq, replace = TRUE)
    starts <- starts[starts > 0]
    insertLens <- as.integer(rtruncnorm(n = length(starts), a = readLen,
                                        mean = meanInsertLen,
                                        sd = sdInsertLen))
    ends <- starts + insertLens - 1
    starts <- starts[ends < length(fullSeq)]
    ends <- ends[ends < length(fullSeq)]

    revComp <- reverseComplement(fullSeq)
    noiseRate <- noiseRate / 0.75
    highNoiseRate <- highNoiseRate / 0.75
    fullLen <- length(fullSeq)

    if (pairedEnd && length(starts) > 0) {
        reads1 <- extractAt(x = fullSeq,
                            at = IRanges(start = starts, width = readLen))
        reads2 <- extractAt(x = revComp,
                            at = IRanges(start = fullLen - ends + 1,
                                         width = readLen))

        removedReads <- union(which(grepl("N", reads1)),
                              which(grepl("N", reads2)))
        if (length(removedReads) > 0) {
            reads <- c(reads1[-removedReads], reads2[-removedReads])
        } else {
            reads <- c(reads1, reads2)
        }

        nReads <- length(reads)
        if (nReads > 0) {
            noiseRates <- rep(noiseRate, length(reads))
            noiseRates[sample(seq_len(nReads),
                              size = nReads * highNoiseProb,
                              replace = TRUE)] <- highNoiseRate
            reads <- addNoise(reads, noiseRates=noiseRates, readLen = readLen)
            reads <- as.character(reads)
            mid <- nReads / 2
            mid2 <- round(mid / 2)
            reads1 <- c(rbind(reads[seq_len(mid2)],
                              reads[(mid + 1):(mid + mid2)]))
            reads2 <- c(rbind(reads[(mid + mid2 + 1):nReads],
                              reads[(mid2 + 1):mid]))
            reads <- c(reads1, reads2)
            return(reads)
        }

    } else if (!pairedEnd && length(starts) > 0){
        mid <- round(nSeq / 2)
        reads1 <- extractAt(x = fullSeq,
                            at = IRanges(start = starts[seq_len(mid)],
                                         width = readLen))
        reads2 <-
            extractAt(x = revComp,
                      at = IRanges(start = fullLen - ends[(mid+1):nSeq] + 1,
                                   width = readLen))
        reads <- c(reads1, reads2)

        reads <- reads[!which(grepl("N", reads))]
        nReads <- length(reads)
        if (nReads > 0) {
            noiseRates <- rep(noiseRate, length(reads))
            noiseRates[sample(seq_len(nReads),
                              size = nReads * highNoiseProb,
                              replace = TRUE)] <- highNoiseRate
            reads <- addNoise(reads, noiseRates=noiseRates, readLen = readLen)
            reads <- as.character(reads)
            return(reads)
        }
    }
}
