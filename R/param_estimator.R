calcPhredScoreProfile <- function(qual) {
    qual <- as.character(qual)
    qual <- lapply(qual, function(x) unlist(strsplit(x, "")))
    qual <- lapply(qual, function(x) {names(x) <- seq(1, length(x)); return(x)})
    qual <- as.data.frame(do.call("bind_rows", qual))
    qualProb <- mapply("/", apply(qual, 2, table), apply(qual, 2, length),
                       SIMPLIFY = FALSE)
    qualProb <- as.data.frame(do.call("bind_rows", qualProb))
    qualProb <- qualProb[, order(phred2ASCIIOffset(names(qualProb)))]
    qualProb[is.na(qualProb)] <- 0
    message('Calculation done.')
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
                     isSupplementaryAlignment = isSupplementary,)

    p <- ScanBamParam(flag = f,
                      tag = tag,
                      simpleCigar = FALSE,
                      reverseComplement = reverseComplement,
                      what = what,
                      which = which,
                      mapqFilter = mapqFilter)

    bam <- scanBam(bamFile, param = p)[[1]]
    message(paste("Found", length(bam), "reads."))
    return(bam)
}


estimateSimParam <- function(bamFilePath, mapqFilter=0, maxFileSize=1,
                             targetRegions=NULL, subsampleRatio=NA,
                             subsampleRegionLength=1e+5,
                             disableSubsampling=FALSE, threads=1,
                             calcPhredProfile=TRUE) {

    bamFile <- BamFile(bamFilePath)
    seqInfo <- scanBamHeader(bamFile)$targets
    if (file.info(bamFilePath)$size / 1e+9 <= maxFileSize |
        disableSubsampling) {
        if (length(targetRegions) == 0) {
            regions <- GRanges(names(seqInfo), IRanges(1, seqInfo))
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
        message(paste0("BAM filesize greater than ",
                       maxFileSize,
                       " gigabytes, calculation based on subsampled reads ",
                       "(subsample ratio: ",
                       subsampleRatio, ")."))
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
            colnames(targetRegions) <- c("chr", "start", "end")
            regions <- GRanges(targetRegions$chr, IRanges(targetRegions$start,
                                                          targetRegions$end))
            regions <- regions[sample(seq_along(regions),
                                      round(length(regions) * subsampleRatio))]
        }
    }
    message(paste("Subsample BAM file in", length(regions), "regions."))
    message(paste("Using", threads, "thread(s)."))
    cl <- makeCluster(threads)
    registerDoParallel(cl)

    if (calcPhredProfile) {
        what = c("qual","isize", "flag", "qwidth")
    } else {
        what = c("isize", "flag", "qwidth")
    }

    readInfo <- foreach(i = seq_along(regions),
                        .combine = "c",
                        .inorder = TRUE,
                        .verbose = FALSE,
                        .packages = c("Biostrings", "dplyr",
                                      "GenomicRanges", "Rsamtools"),
                        .export = c("findReads")
    ) %dopar% {
        findReads(bamFile,
                  what = what,
                  tag=c("NM"),
                  which = regions[i],
                  reverseComplement=TRUE,
                  mapqFilter = mapqFilter)
    }
    stopCluster(cl)
    insertSize <- unlist(unname(readInfo[names(readInfo) == "isize"]))
    if (length(insertSize) == 0) {
        stop("No reads found.")
    }
    message(paste("Calculating based on",
                  length(insertSize), "reads..."))
    readLen <- max(unlist(unname(readInfo[names(readInfo) == "qwidth"])),
                   na.rm = TRUE)
    message("Read length: ", readLen)
    medianIns <- median(insertSize[insertSize >= readLen], na.rm = TRUE)
    MADIns <- mad(insertSize[insertSize >= readLen], na.rm = TRUE)
    message("Estimated median insert size (for insertsize >= read length): ",
            medianIns)
    message("Estimated median absolute deviation of ",
            "insert size (for insertsize >= read length): ", MADIns)
    flag <- unlist(unname(readInfo[names(readInfo) == "flag"]))
    isSupAlign <- vapply(flag, function(x) as.integer(intToBits(x))[12],
                         numeric(1))
    isProperPaired <- vapply(flag, function(x) as.integer(intToBits(x))[2],
                             numeric(1))
    supAlignProp <- sum(isSupAlign)/length(isSupAlign)
    message("Estimated proportion of supplementary alignment: ", supAlignProp)
    IPPProp <- sum(!isProperPaired) / length(isProperPaired)
    message("Estimated proportion of read NOT mapped in proper pair: ", IPPProp)
    message("Estimated proportion of supplementary alignments in reads ",
            "NOT mapped in proper pair: ",
            (sum(!isProperPaired & isSupAlign)) / sum(!isProperPaired))
    NMTable <- table(unname(unlist(readInfo[names(readInfo) == "tag"])))
    NMTable <- NMTable / sum(NMTable)
    message("Proportion of reads with N edit distance:")
    message(paste0(capture.output(NMTable), collapse = "\n"))

    if (calcPhredProfile) {
        qual <- do.call("c", unname(readInfo[names(readInfo) == "qual"]))
        qualProb <- calcPhredScoreProfile(qual)
        return(qualProb)
    }
}

