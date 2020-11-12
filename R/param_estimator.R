countPhredScores <- function(qual) {
    qual <- as.character(qual)
    qual <- lapply(qual, function(x) unlist(strsplit(x, "")))
    qual <- lapply(qual, function(x) {names(x) <- seq(1, length(x)); return(x)})
    qual <- as.data.frame(do.call("bind_rows", qual))
    qualProb <- apply(qual, 2, table)
    qualProb <- as.matrix(do.call("bind_rows", qualProb))
    qualProb <- qualProb[, order(phred2ASCIIOffset(colnames(qualProb)))]
    qualProb[is.na(qualProb)] <- 0
    return(qualProb)
}


findReads <- function(bamFile,
                      tag = character(0),
                      what=scanBamWhat(),
                      which=IRangesList(),
                      isSupplementary=NA,
                      reverseComplement=NA,
                      mapqFilter=0) {

    f <- scanBamFlag(isDuplicate = FALSE,
                     isSupplementaryAlignment = isSupplementary)

    p <- ScanBamParam(flag = f,
                      tag = tag,
                      simpleCigar = FALSE,
                      reverseComplement = reverseComplement,
                      what = what,
                      which = which,
                      mapqFilter = mapqFilter)

    bam <- scanBam(bamFile, param = p)[[1]]
    return(bam)
}


calcPhredScoreProfile <- function(bamFilePath, mapqFilter=0, maxFileSize=1,
                                  targetRegions=NULL, subsampleRatio=NA,
                                  subsampleRegionLength=1e+5,
                                  disableSubsampling=FALSE, threads=1) {
    
    if(missing(bamFilePath)) stop("bamFilePath is required")
    
    bamFile <- BamFile(bamFilePath)
    seqInfo <- scanBamHeader(bamFile)$targets
    
    colnames(targetRegions) <- c("chr", "start", "end")
    targetRegions <- targetRegions[targetRegions$chr %in% names(seqInfo), ]
    
    if (length(targetRegions) != 0 & is(targetRegions, "GenomicRanges")) {
        targetRegions <- data.frame(seqnames(targetRegions), 
                                    start(targetRegions), 
                                    end(targetRegions))
    }
    
    if (file.info(bamFilePath)$size / 1e+9 <= maxFileSize |
        disableSubsampling) {
        if (length(targetRegions) == 0) {
            regions <- GRanges(names(seqInfo), IRanges(1, seqInfo))
            regions <- 
                unlist(slidingWindows(regions, width=50000, step = 50000))
        } else {
            colnames(targetRegions) <- c("chr", "start", "end")
            regions <- GRanges(targetRegions$chr, IRanges(targetRegions$start,
                                                          targetRegions$end))
        }
        if (file.info(bamFilePath)$size / 1e+9 <= maxFileSize) {
            message("BAM filesize smaller than ",
                    maxFileSize,
                    " gigabytes, calculate Phred quality score profile ",
                    "based on all reads.")
        } else {
            message("Subsampling disabled, ",
                    "calculate Phred quality score profile based on all reads.")
        }
    } else {
        if (is.na(subsampleRatio)) {
            subsampleRatio <-
                maxFileSize * 1e+9 / file.info(bamFilePath)$size
        }
        message("BAM filesize greater than ",
                maxFileSize,
                " gigabytes, calculation based on subsampled reads ",
                "(subsample ratio: ",
                subsampleRatio, ").")
        if (length(targetRegions) == 0) {
            nRegions <- seqInfo * subsampleRatio / subsampleRegionLength
            seqInfo <- seqInfo[nRegions >= 1]
            nRegions <- ceiling(nRegions[nRegions >= 1])
            if (length(seqInfo) == 0) {
                stop("subsampleRegionLength is too large ",
                     "with this subsampleRatio.")
            }
            subsampleRegions <- GRanges()
            for (i in seq_along(seqInfo)) {
                startPos <- sample(seq(1, seqInfo[i], subsampleRegionLength),
                                   nRegions[i], replace = FALSE)
                subsampleRegion <-
                    GRanges(seqnames = names(seqInfo[i]),
                            ranges = IRanges(startPos,
                                             startPos +
                                                 subsampleRegionLength - 1))

                suppressWarnings(subsampleRegions <-
                                     c(subsampleRegions, subsampleRegion))
            }
            regions <- sort(subsampleRegions)
        } else {
            regions <- GRanges(targetRegions$chr, IRanges(targetRegions$start,
                                                          targetRegions$end))
            regions <- regions[sample(seq_along(regions),
                                      round(length(regions) * subsampleRatio))]
        }
    }
    
    message("Get reads from BAM file in ", length(regions), " regions.", 
            "\nUsing ", threads, " thread(s).")
    cl <- makeCluster(threads)
    registerDoParallel(cl)

    readInfo <- foreach(i = seq_along(regions),
                        .combine = "c",
                        .inorder = TRUE,
                        .verbose = FALSE,
                        .errorhandling = "remove",
                        .packages = c("Biostrings", "dplyr",
                                      "GenomicRanges", "Rsamtools"),
                        .export = c("findReads", "countPhredScores")
    ) %dopar% {
            temp <- findReads(bamFile,
                              what = "qual",
                              which = regions[i],
                              reverseComplement=TRUE,
                              mapqFilter = mapqFilter)
            if (length(temp$qual) > 1) {
                temp$qual <- countPhredScores(temp$qual)
                return(temp)
            }
    }
    stopCluster(cl)
    quals <- unname(readInfo[names(readInfo) == "qual"])
    if (length(quals) == 0) {
        stop("No reads found.")
    }
    readLen <- max(unlist(lapply(quals, nrow)))
    PhredScores <- unique(unlist(lapply(quals, colnames)))
    PhredScores <- PhredScores[order(phred2ASCIIOffset(PhredScores))]
    for (i in seq_along(quals)) {
        qual <- quals[[i]]
        PhredScoreDiff <- setdiff(PhredScores, colnames(qual))
        if (length(PhredScoreDiff) > 0 ) {
            temp = matrix(0, nrow = nrow(qual), ncol = length(PhredScoreDiff))
            colnames(temp) <- PhredScoreDiff
            quals[[i]] <- cbind(qual, temp)
            quals[[i]] <- quals[[i]][, PhredScores]
        }
        rowDiff <- readLen - nrow(qual)
        if (rowDiff > 0) {
            temp = matrix(0, nrow = rowDiff, ncol = length(PhredScores))
            quals[[i]] <- rbind(quals[[i]], temp)
        }
    }
    qualProb <- Reduce("+", quals)
    qualProb <- qualProb / rowSums(qualProb)
    return(qualProb)
}

