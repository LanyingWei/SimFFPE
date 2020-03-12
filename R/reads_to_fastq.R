readsToFastq <- function(simReads, PhredScoreProfile, prefix, prefix2,
                         chr, outFile, pairedEnd, threads, firstID=1) {
    if (length(simReads) > 0) {
        PhredScores <- generatePhredScores(
            nReads = length(simReads),
            PhredScoreProfile = PhredScoreProfile,
            threads = threads,
            existedCluster = TRUE)

        readNames <- generateReadNames(nReads = length(simReads),
                                       prefix = paste0(
                                           prefix, '-',
                                           prefix2, '-',
                                           chr),
                                       firstID = firstID,
                                       pairedEnd = pairedEnd)

        writeReadsToFastq(simReads = simReads, PhredScores = PhredScores,
                          readNames = readNames, outFile = outFile,
                          overWrite = FALSE, pairedEnd = pairedEnd)
    }
}

generateReadNames <- function(nReads, prefix='SimFFPE',
                              firstID=1, pairedEnd=TRUE) {
    if(nReads == 0 ) return()
    if (pairedEnd) {
        name1 <- paste0('@', prefix, '-',
                        firstID:(nReads / 2 + firstID - 1), '/1')
        name2 <- paste0('@', prefix, '-',
                        firstID:(nReads / 2 + firstID - 1), '/2')
        names <- c(rbind(name1,name2))
    } else{
        names <- paste0('@', prefix, '-', firstID:(nReads + firstID - 1))
    }
    return(names)
}


generatePhredScores <- function(nReads, PhredScoreProfile,
                                threads=1, existedCluster=FALSE) {
    if (!existedCluster) {
        cl <- makeCluster(threads)
        registerDoParallel(cl)
    }

    PhredScores <- NULL
    n <- ceiling(nReads/threads)

    PhredScores <-
        foreach(i = seq_len(threads),
                .combine='c',
                .inorder = TRUE,
                .verbose = FALSE,
                .export = c('generatePhredScore')) %dopar% {
                    generatePhredScore(nReads = n,
                                       PhredScoreProfile = PhredScoreProfile)
                }

    if (!existedCluster) {
        stopCluster(cl)
    }
    PhredScores <- PhredScores[seq_len(nReads)]
    return(PhredScores)
}


generatePhredScore <- function(nReads, PhredScoreProfile) {
    PhredScores <- data.frame((matrix(NA, nrow = nReads,
                                      ncol = nrow(PhredScoreProfile))))
    PhredScoreLetters <- colnames(PhredScoreProfile)
    for (i in seq_len(nrow(PhredScoreProfile))) {
        message(paste("Generating the", i, "row of Phred Scores"))
        positionalPhredScores <- sample(PhredScoreLetters,
                                        size = nReads,
                                        replace = TRUE,
                                        prob = PhredScoreProfile[i, ])
        PhredScores[, i] <- positionalPhredScores
    }
    PhredScores <- do.call('paste0', PhredScores)
    return(PhredScores)
}

overWriteFastq <- function (outFile, pairedEnd) {
    if (pairedEnd) {
        fastq1 <- file(paste0(outFile, '_1.fq'), 'w')
        fastq2 <- file(paste0(outFile, '_2.fq'), 'w')
        close(fastq1)
        close(fastq2)
    } else {
        fastq <- file(paste0(outFile, '.fq'), 'w')
        close(fastq)
    }
}

writeReadsToFastq <- function(simReads, PhredScores, readNames, outFile,
                              pairedEnd, overWrite=FALSE) {
    if (overWrite) {
        overWriteFastq(outFile = outFile, pairedEnd = pairedEnd)
    }

    if (pairedEnd) {
        fastq1 <- file(paste0(outFile, '_1.fq'), 'a')
        fastq2 <- file(paste0(outFile, '_2.fq'), 'a')
        index1 <- seq(from = 1, to = length(simReads) - 1, by = 2)
        index2 <- seq(from = 2, to = length(simReads), by = 2)
        readInfo <- c(rbind(readNames[index1], simReads[index1],
                            rep("+", length(index1)), PhredScores[index1]))
        writeLines(readInfo, fastq1)
        close(fastq1)
        readInfo <- c(rbind(readNames[index2], simReads[index2],
                            rep("+", length(index2)), PhredScores[index2]))
        writeLines(readInfo, fastq2)
        close(fastq2)
    } else {
        fastq <- file(paste0(outFile, '.fq'), 'a')
        readInfo <- c(rbind(readNames, simReads,
                            rep('+', length(simReads)), PhredScores))
        writeLines(readInfo, fastq)
        close(fastq)
    }
}
