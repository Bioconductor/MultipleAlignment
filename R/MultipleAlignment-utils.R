### =========================================================================
### MultipleAlignment utilities
### -------------------------------------------------------------------------
###

setMethod("consensusMatrix","MultipleAlignment",
    function(x, as.prob=FALSE, baseOnly=FALSE)
    {
        strings <- unmasked(x)
        if (maskednrow(x) > 0)
            strings <- strings[- as.integer(rowmask(x))]
        m <- callGeneric(strings, as.prob=as.prob, baseOnly=baseOnly,
                         width=ncol(x))
        if (maskedncol(x) > 0)
            m[, as.integer(colmask(x))] <- NA
        m
    }
)

setMethod("consensusString","MultipleAlignment",
    function(x, ambiguityMap, threshold, codes)
    {
        consensus <- rep.int("#", ncol(x))
        if (ncol(x) > 0) {
            cmat <- consensusMatrix(x, width=ncol(x))
            ngaps <- cmat["-",]
            cmat <- cmat[codes, , drop=FALSE]
            colsums <- colSums(cmat, na.rm=TRUE)
            nzsum <- which(colsums > 1e-6)
            if (length(nzsum) > 0) {
                colsums[- nzsum] <- 1  # to avoid division by 0
                gaplocs <- which(as.logical(ngaps > colsums))
                nzsum <- setdiff(nzsum, gaplocs)
                cmat <- cmat / rep(colsums, each=nrow(cmat))
                consensus[gaplocs] <- "-"
                consensus[nzsum] <-
                  safeExplode(consensusString(cmat[, nzsum, drop=FALSE],
                                              ambiguityMap=ambiguityMap,
                                              threshold=threshold))
            }
            consensus <- paste(consensus, collapse="")
        }
        consensus
    }
)

setMethod("consensusString","DNAMultipleAlignment",
    function(x, ambiguityMap=IUPAC_CODE_MAP, threshold=0.25)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold,
                       codes=names(IUPAC_CODE_MAP))
    }
)

setMethod("consensusString","RNAMultipleAlignment",
    function(x,
            ambiguityMap=
            structure(as.character(RNAStringSet(DNAStringSet(IUPAC_CODE_MAP))),
                      names=
                      as.character(RNAStringSet(DNAStringSet(names(IUPAC_CODE_MAP))))),
            threshold=0.25)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold,
                       codes=
                       as.character(RNAStringSet(DNAStringSet(names(IUPAC_CODE_MAP)))))
    }
)

setMethod("consensusString","AAMultipleAlignment",
    function(x, ambiguityMap="?", threshold=0.5)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold,
                       codes=names(AMINO_ACID_CODE))
    }
)

setGeneric("consensusViews", signature="x",
    function(x, ...) standardGeneric("consensusViews")
)

setMethod("consensusViews","MultipleAlignment",
    function(x, ambiguityMap, threshold)
    {
        cmask <- colmask(x)
        if (length(cmask) > 0)
            colmask(x) <- NULL
        consensus <-
          consensusString(x, ambiguityMap=ambiguityMap, threshold=threshold)
        if (length(consensus) == 0)
            consensus <- ""
        Views(BString(consensus), gaps(cmask, start=1, end=ncol(x)))
    }
)

setMethod("consensusViews","DNAMultipleAlignment",
    function(x, ambiguityMap=IUPAC_CODE_MAP, threshold=0.25)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold)
    }
)

setMethod("consensusViews","RNAMultipleAlignment",
    function(x,
             ambiguityMap=as.character(RNAStringSet(DNAStringSet(IUPAC_CODE_MAP))),
             threshold=0.25)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold)
    }
)

setMethod("consensusViews","AAMultipleAlignment",
    function(x, ambiguityMap="?", threshold=0.5)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold)
    }
)

setMethod("alphabetFrequency","MultipleAlignment",
    function(x, as.prob=FALSE, collapse=FALSE)
    {
        if (collapse) {
            strings <- as(x, sprintf("%sStringSet", seqtype(x)))
            m <- callGeneric(strings, as.prob=as.prob, collapse=TRUE)
        } else {
            rmask <- rowmask(x)
            if (length(rmask) > 0)
                rowmask(x) <- NULL
            strings <- as(x, sprintf("%sStringSet", seqtype(x)))
            m <- callGeneric(strings, as.prob=as.prob)
            if (length(rmask) > 0)
                m[as.integer(rmask),] <- NA
        }
        m
    }
)

