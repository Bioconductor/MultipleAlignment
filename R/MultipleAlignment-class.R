### =========================================================================
### MultipleAlignment objects
### -------------------------------------------------------------------------
###

## virtual class for multiple alignments
setClass("MultipleAlignment",
    representation("VIRTUAL",
                   unmasked="XStringSet",
                   rowmask="NormalIRanges",
                   colmask="NormalIRanges")
)

## concrete classes for multiple alignments
setClass("DNAMultipleAlignment",
    contains="MultipleAlignment",
    representation=representation(unmasked="DNAStringSet")
)
setClass("RNAMultipleAlignment",
    contains="MultipleAlignment",
    representation=representation(unmasked="RNAStringSet")
)
setClass("AAMultipleAlignment",
    contains="MultipleAlignment",
    representation=representation(unmasked="AAStringSet")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.MultipleAlignment.unmasked <- function(object)
{
    if (length(unique(nchar(object@unmasked))) > 1)
        "alignments have an unequal number of characters"
    else
        NULL
}

.valid.MultipleAlignment.rowmask <- function(object)
{
    rng <- range(object@rowmask)
    if (length(rng) > 0 &&
        (start(rng) < 1 || end(rng) > length(object@unmasked)))
        "'rowmask' contains values outside of [1, nrow]"
    else
        NULL
}

.valid.MultipleAlignment.colmask <- function(object)
{
    rng <- range(object@colmask)
    if (length(rng) > 0 &&
        (start(rng) < 1 || end(rng) > nchar(object@unmasked)[1L]))
        "'colmask' contains values outside of [1, ncol]"
    else
        NULL
}

setValidity("MultipleAlignment",
    function(object)
    {
        problems <-
          c(.valid.MultipleAlignment.unmasked(object),
            .valid.MultipleAlignment.rowmask(object),
            .valid.MultipleAlignment.colmask(object))
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###
### The *MultipleAlignment class definitions were moved from Biostrings to
### the MultipleAlignment package in BioC 3.24.

update_MultipleAlignment_object <- function(object, ..., verbose=FALSE)
{
    object_class <- class(object)
    class_package_attr <- attr(object_class, "package")
    if (class_package_attr != "MultipleAlignment") {
        if (verbose)
            message("[updateObject] class package attribute ",
                    "on ", object_class, " object ",
                    "is \"", class_package_attr, "\".\n",
                    "[updateObject] Setting it to \"MultipleAlignment\" ... ",
                    appendLF=FALSE)
        class(object) <- class(get(object_class, mode="function")())
        if (verbose)
            message("OK")
    }
    callNextMethod()
}

setMethod("updateObject", "DNAMultipleAlignment",
          update_MultipleAlignment_object)
setMethod("updateObject", "RNAMultipleAlignment",
          update_MultipleAlignment_object)
setMethod("updateObject", "AAMultipleAlignment",
          update_MultipleAlignment_object)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("unmasked", "MultipleAlignment", function(x) x@unmasked)

setMethod("rownames","MultipleAlignment", function(x) names(x@unmasked))
setReplaceMethod("rownames", "MultipleAlignment",
    function(x,value)
    {
        names(x@unmasked) <- value
        x
    }
)

setGeneric("rowmask", signature="x", function(x) standardGeneric("rowmask"))
setMethod("rowmask", "MultipleAlignment", function(x) x@rowmask)

.setMask <- function(mask, append, invert, length, value){
    if (!isTRUEorFALSE(invert))
        stop("'invert' must be TRUE or FALSE")
    if (invert)
        value <- gaps(value, start=1L, end=length)
    append <- match.arg(append, choices=c("union", "replace", "intersect"))
    value <- switch(append,
                  "union"=union(mask, value),
                  "replace"=value,
                  "intersect"=intersect(mask, value))
    as(value, "NormalIRanges")
}
setGeneric("rowmask<-", signature=c("x", "value"),
    function(x, append="union", invert=FALSE, value)
           standardGeneric("rowmask<-")
)

setReplaceMethod("rowmask", signature(x="MultipleAlignment", value="NULL"),
    function(x, append="replace", invert=FALSE, value)
        callGeneric(x, append=append, invert=invert, value=new("NormalIRanges"))
)
setReplaceMethod("rowmask",
    signature(x="MultipleAlignment", value="ANY"),
    function(x, append="union", invert=FALSE, value)
    {
        if (!is(value, "IRanges"))
            value <- as(value, "IRanges")
        value <- .setMask(mask=rowmask(x), append=append, invert=invert,
                        length=dim(x)[1], value=value)
        initialize(x, rowmask = value)
    }
)

setGeneric("colmask", signature="x", function(x) standardGeneric("colmask"))
setMethod("colmask", "MultipleAlignment", function(x) x@colmask)

setGeneric("colmask<-", signature=c("x", "value"),
    function(x, append="union", invert=FALSE, value)
           standardGeneric("colmask<-")
)
setReplaceMethod("colmask", signature(x="MultipleAlignment", value="NULL"),
    function(x, append="replace", invert=FALSE, value)
        callGeneric(x, append=append, invert=invert,
                    value=new("NormalIRanges"))
)
setReplaceMethod("colmask",
    signature(x="MultipleAlignment", value="ANY"),
    function(x, append="union", invert=FALSE, value)
    {
        if (!is(value, "IRanges"))
            value <- as(value, "IRanges")
        value <- .setMask(mask=colmask(x), append=append, invert=invert,
                        length=dim(x)[2], value=value)
        initialize(x, colmask = value)
    }
)

setMethod("maskMotif", signature(x="MultipleAlignment", motif="ANY"),
    function(x, motif, min.block.width=1, ...)
    {
        cmask <- colmask(x)
        if (length(colmask(x)) > 0)
            colmask(x) <- NULL
        string <- consensusString(x)
        if (length(string) == 1) {
            string <- gsub("#", "+", string)
            string <- as(string, xsbaseclass(unmasked(x)))
            maskedString <-
                callGeneric(string, motif, min.block.width=min.block.width, ...)
            newmask <- nir_list(masks(maskedString))[[1L]]
            colmask(x) <- union(newmask, cmask)
        }
        x
    }
)

setGeneric("maskGaps", signature="x",
    function(x, ...) standardGeneric("maskGaps")
)
setMethod("maskGaps", "MultipleAlignment",
    function(x, min.fraction=0.5, min.block.width=3)
    {
        if (!isSingleNumber(min.fraction) || min.fraction < 0 ||
            min.fraction > 1)
            stop("'min.fraction' must be a number in [0, 1]")
        if (!isSingleNumber(min.block.width) || min.block.width <= 0)
            stop("'min.block.width' must be a positive integer")
        if (!is.integer(min.block.width))
            min.block.width <- as.integer(min.block.width)

        cmask <- colmask(x)
        if (length(colmask(x)) > 0)
            colmask(x) <- NULL
        m <- consensusMatrix(x)
        newmask <- (m["-",] / colSums(m)) >= min.fraction
        newmask[is.na(newmask)] <- FALSE
        newmask <- as(newmask, "NormalIRanges")
        newmask <- newmask[width(newmask) >= min.block.width]
        colmask(x) <- union(newmask, cmask)
        x
    }
)

setMethod("nrow","MultipleAlignment", function(x) length(x@unmasked))
setMethod("ncol","MultipleAlignment",
    function(x) {if (nrow(x) == 0L) 0L else nchar(x@unmasked[[1L]])}
)
setMethod("dim","MultipleAlignment", function(x) c(nrow(x), ncol(x)))

setGeneric("maskednrow", signature="x",
    function(x) standardGeneric("maskednrow")
)
setMethod("maskednrow", "MultipleAlignment",
    function(x) sum(width(rowmask(x)))
)

setGeneric("maskedncol", signature="x",
    function(x) standardGeneric("maskedncol")
)
setMethod("maskedncol", "MultipleAlignment",
    function(x) sum(width(colmask(x)))
)

setGeneric("maskeddim", signature="x",
    function(x) standardGeneric("maskeddim")
)
setMethod("maskeddim", "MultipleAlignment",
    function(x) c(maskednrow(x), c(maskedncol(x)))
)

setMethod("maskedratio", "MultipleAlignment",
    function(x) maskeddim(x) / dim(x)
)

setMethod("nchar", "MultipleAlignment",
    function(x) ncol(x) - maskedncol(x)
)

setMethod("seqtype", "MultipleAlignment",
    function(x) seqtype(unmasked(x))
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

.new_MultipleAlignment <- function(seqtype, x, start, end, width, use.names,
                                   rowmask, colmask)
{
    ans_class <- paste0(seqtype, "MultipleAlignment")
    unmasked <- get(paste0(seqtype, "StringSet"))(x=x,
                                                  start=start,
                                                  end=end,
                                                  width=width,
                                                  use.names=use.names)
    if (is.null(rowmask))
        rowmask <- new("NormalIRanges")
    if (is.null(colmask))
        colmask <- new("NormalIRanges")
    new(ans_class, unmasked=unmasked, rowmask=rowmask, colmask=colmask)
}

DNAMultipleAlignment <- function(x=character(),
                                 start=NA, end=NA, width=NA, use.names=TRUE,
                                 rowmask=NULL, colmask=NULL)
{
    .new_MultipleAlignment("DNA", x, start, end, width, use.names,
                           rowmask, colmask)
}

RNAMultipleAlignment <- function(x=character(),
                                 start=NA, end=NA, width=NA, use.names=TRUE,
                                 rowmask=NULL, colmask=NULL)
{
    .new_MultipleAlignment("RNA", x, start, end, width, use.names,
                           rowmask, colmask)
}

AAMultipleAlignment <- function(x=character(),
                                start=NA, end=NA, width=NA, use.names=TRUE,
                                rowmask=NULL, colmask=NULL)
{
    .new_MultipleAlignment("AA", x, start, end, width, use.names,
                           rowmask, colmask)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show detailed pager view of unmasked sequences.
###

setMethod("detail", "MultipleAlignment",
    function(object, invertColMask=FALSE, hideMaskedCols=TRUE, ...)
    {
        ## We don't want a permanent file for this
        FH <- tempfile(pattern = "tmpFile", tmpdir = tempdir())
        ## Then write out a temp file with the correct stuff in it
        write.MultAlign(object, FH, invertColMask=invertColMask,
                        showRowNames=TRUE,
                        hideMaskedCols=hideMaskedCols)
        ## use file.show() to display
        file.show(FH)
    }
)

## TODO: explore if we can make a textConnection() to save time writing to disk.
## page() allows this, but seems to be actually slower than what we do here,
## also, it looks considerably less nice than when it comes from a file...


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("MultipleAlignment", "DNAStringSet",
    function(from) DNAStringSet(as.character(from))
)
setAs("MultipleAlignment", "RNAStringSet",
    function(from) RNAStringSet(as.character(from))
)
setAs("MultipleAlignment", "AAStringSet",
    function(from) AAStringSet(as.character(from))
)
setAs("MultipleAlignment", "BStringSet",
    function(from) BStringSet(as.character(from))
)

setAs("character", "DNAMultipleAlignment",
      function(from) DNAMultipleAlignment(from)
)
setAs("character", "RNAMultipleAlignment",
      function(from) RNAMultipleAlignment(from)
)
setAs("character", "AAMultipleAlignment",
      function(from) AAMultipleAlignment(from)
)

setMethod("as.character", "MultipleAlignment",
    function(x, use.names=TRUE)
        apply(as.matrix(x, use.names=use.names), 1, paste, collapse="")
)

setMethod("as.matrix", "MultipleAlignment",
    function(x, use.names=TRUE)
    {
        if (nrow(x) == 0) {
            m <- matrix(character(0), nrow=0, ncol=0)
        } else {
            m <- as.matrix(unmasked(x), use.names=use.names)
            if (maskednrow(x) > 0)
                m <- m[- as.integer(rowmask(x)), , drop=FALSE]
            if (maskedncol(x) > 0)
                m <- m[, - as.integer(colmask(x)), drop=FALSE]
        }
        m
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

xsbaseclass <- Biostrings:::xsbaseclass
.namesW <- Biostrings:::.namesW
toSeqSnippet <- Biostrings:::toSeqSnippet
add_colors <- Biostrings:::add_colors
add_colors.default <- Biostrings:::add_colors.default
add_colors.DNA <- Biostrings:::add_colors.DNA
add_colors.RNA <- Biostrings:::add_colors.RNA
add_colors.AA <- Biostrings:::add_colors.AA

.MultipleAlignment.show_frame_header <-
function (iW, with.names)
{
    cat(format("", width = iW + 1), sep="")
    if (with.names) {
        cat(format(" aln", width = getOption("width") - iW - .namesW - 1),
            format("names", width=.namesW, justify="left"), sep="")
    } else {
        cat(" aln")
    }
    cat("\n")
}

.MultipleAlignment.show_frame_line <-
function (x, i, iW)
{
    snippetWidth <- getOption("width") - 2L - iW
    if (!is.null(names(x)))
        snippetWidth <- snippetWidth - .namesW - 1L
    snippet <- toSeqSnippet(x[[i]], snippetWidth)
    if (!is.null(names(x))) {
        snippet_class <- class(snippet)
        snippet <- format(snippet, width=snippetWidth)
        class(snippet) <- snippet_class
    }
    cat(format(paste("[", i, "]", sep=""), width=iW, justify="right"), " ",
        add_colors(snippet), sep="")
    if (!is.null(names(x))) {
        snippet_name <- names(x)[i]
        if (is.na(snippet_name))
            snippet_name <- "<NA>"
        else if (nchar(snippet_name) > .namesW)
            snippet_name <- paste0(substr(snippet_name, 1L, .namesW - 3L),
                                   #compact_ellipsis)
                                   "...")
        cat(" ", snippet_name, sep="")
    }
    cat("\n")
}

.MultipleAlignment.show_frame <-
function (x, half_nrow=9L)
{
    lx <- length(x)
    iW <- nchar(as.character(lx)) + 2
    .MultipleAlignment.show_frame_header(iW, !is.null(names(x)))
    if (lx <= 2 * half_nrow + 1) {
        for (i in seq_len(lx))
            .MultipleAlignment.show_frame_line(x, i, iW)
    } else {
        for (i in seq_len(half_nrow))
            .MultipleAlignment.show_frame_line(x, i, iW)
        cat(format("...", width=iW, justify="right"), "...\n")
        for (i in (lx - half_nrow + 1L):lx)
            .MultipleAlignment.show_frame_line(x, i, iW)
    }
}

setMethod("show", "MultipleAlignment",
    function(object)
    {
        nr <- nrow(object)
        nc <- ncol(object)
        cat(class(object), " with ", nr, ifelse(nr == 1, " row and ",
            " rows and "), nc, ifelse(nc == 1, " column\n", " columns\n"),
            sep = "")
        if (nr != 0) {
            strings <- unmasked(object)
            mdim <- maskeddim(object)
            if (sum(mdim) > 0) {
                if (mdim[1] > 0) {
                    strings <- BStringSet(strings)
                    maskStrings <-
                      rep(BStringSet(paste(rep.int("#", ncol(object)),
                                           collapse="")), mdim[1])
                    i <- as.integer(rowmask(object))
                    if (!is.null(rownames(object)))
                        names(maskStrings) <- rownames(object)[i]
                    strings[i] <- maskStrings
                }
                if (mdim[2] > 0) {
                    strings <- as.matrix(strings)
                    strings[, as.integer(colmask(object))] <- "#"
                    strings <- BStringSet(apply(strings, 1, paste, collapse=""))
                }
            }
            .MultipleAlignment.show_frame(strings)
        }
    }
)

