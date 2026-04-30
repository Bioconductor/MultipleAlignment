### =========================================================================
### Read/write MultipleAlignment objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read function.
###

.read_XStringSet <- Biostrings:::.read_XStringSet

## markupPattern specifies which lines to skip
.read.MultipleAlignment.splitRows <-
function(rows, markupPattern)
{
    markupLines <- grep(markupPattern, rows, perl=TRUE)
    alnLines <- gaps(as(markupLines, "IRanges"), start=1, end=length(rows))
    nseq <- unique(width(alnLines))
    if (length(nseq) != 1)
        stop("missing alignment rows")
    rows <- extractROWS(rows, alnLines)
    spaces <- regexpr("\\s+", rows)
    ids <- substr(rows, 1L, spaces - 1L)
    nsplits <- length(rows) %/% nseq
    if (!identical(ids, rep.int(head(ids, nseq), nsplits)))
        stop("alignment rows out of order")
    alns <- substr(rows, spaces + attr(spaces, "match.length"), nchar(rows))
    structure(do.call(paste,
                      c(split(alns, rep(seq_len(nsplits), each=nseq)), sep="")),
              names = head(ids, nseq))
}

.read.Stockholm <-
function(filepath)
{
    rows <- scan(filepath, what = "", sep = "\n", strip.white = TRUE,
                   quiet = TRUE, blank.lines.skip = FALSE)
    if (length(rows) < 3 ||
        !identical(grep("^# STOCKHOLM", rows[1L]), 1L))
        stop("invalid Stockholm file")
    chartr(".", "-",
           .read.MultipleAlignment.splitRows(rows, "(^\\s*|^#.*|^//\\s*)$"))
}

.read.ClustalAln <-
function(filepath)
{
    rows <- scan(filepath, what = "", sep = "\n", strip.white = TRUE,
                   quiet = TRUE, blank.lines.skip = FALSE)
    if (length(rows) < 3 ||
        !identical(grep("^CLUSTAL", rows[1L]), 1L) ||
        !identical(rows[2:3], c("","")))
        stop("invalid Clustal aln file")
    rows <- tail(rows, -3)
    rows <- sub("^(\\S+\\s+\\S+)\\s*\\d*$", "\\1", rows)
    .read.MultipleAlignment.splitRows(rows, "^(\\s|\\*|:|\\.)*$")
}

## In order to recycle .read.MultipleAlignment.splitRows().
## I need to have the names on each row.
.read.PhylipAln <-
function(filepath, maskGen=FALSE)
{
    rows <- scan(filepath, what = "", sep = "\n", strip.white = TRUE,
                 quiet = TRUE, blank.lines.skip = FALSE)
    if (length(rows) < 1 ||
        !identical(grep("^\\d+?\\s\\d+?", rows[1L]), 1L))
        stop("invalid Phylip file")
    ##(mask+num rows + blank line)
    nameLength <- as.numeric(sub("(\\d+).*$","\\1", rows[1])) +1
    rows <- tail(rows, -1)
    names <- character()
    names[nameLength] <- "" ## an empty string is ALWAYS the last "name"
    offset <- 0L
    for(i in seq_len(length(rows))){
        if(i<nameLength){
            rows[i] <- sub("(^\\S+)\\s+(\\S+)", "\\1\\|\\2", rows[i])
            rows[i] <- gsub("\\s", "", rows[i])
            rows[i] <- sub("\\|", " ", rows[i])
            names[i] <- sub("(\\S+).*$","\\1",rows[i])
        } else {
            rows[i] <- gsub("\\s", "", rows[i])
            rows[i] <- paste(names[i %% nameLength], rows[i])
        }
    }
    rows <- c(" ",rows)
    if(maskGen==FALSE){ ## filter out the Mask values OR blank lines
        .read.MultipleAlignment.splitRows(rows, "^(Mask|\\s)")
    } else {## only retrieve the Mask values
        if(length(grep("^(?!Mask)",rows, perl=TRUE))==length(rows)){
            return(as(IRanges(),"NormalIRanges"))
        } else {
            msk <- .read.MultipleAlignment.splitRows(rows, "^(?!Mask)")
            ## THEN cast them to be a NormalIRanges object.
            splt <- strsplit(msk,"") ## split up all chars
            names(splt) <- NULL ## drop the name
            splt <- unlist(splt) ## THEN unlist
            lsplt <- as.logical(as.numeric(splt)) ## NOW you can get a logical
            return(gaps(as(lsplt,"NormalIRanges"))) ## gaps() inverts mask
        }
    }
}

.checkFormat <- function(filepath, format){
    if (missing(format)) {
        ext <- tolower(sub(".*\\.([^.]*)$", "\\1", filepath))
        format <- switch(ext, "sto" = "stockholm", "aln" = "clustal", "fasta")
    } else {
        format <- match.arg(tolower(format), c("fasta", "stockholm", "clustal",
                                               "phylip"))
    }
    format
}

.read.MultipleAlignment <-
function(filepath, format, seqtype, ...)
{
    format <- .checkFormat(filepath, format)
    switch(format,
           "stockholm" = .read.Stockholm(filepath),
           "clustal" = .read.ClustalAln(filepath),
           "phylip" = .read.PhylipAln(filepath),
           .read_XStringSet(filepath, format,
                            nrec=-1L, skip=0L, seek.first.rec=FALSE,
                            use.names=TRUE, seqtype=seqtype))
    ##TODO: BUGs with stockholm??
}

.read.MultipleMask <-
function(filepath, format)
{
    format <- .checkFormat(filepath, format)
    switch(format,
           "stockholm" = as(IRanges(),"NormalIRanges"),
           "clustal" = as(IRanges(),"NormalIRanges"),
           "phylip" = .read.PhylipAln(filepath, maskGen=TRUE),
           as(IRanges(),"NormalIRanges"))
}


readDNAMultipleAlignment <-
function(filepath, format)
{
    DNAMultipleAlignment(.read.MultipleAlignment(filepath, format, "DNA"),
                         rowmask=as(IRanges(),"NormalIRanges"),
                         colmask=.read.MultipleMask(filepath,format))
}

readRNAMultipleAlignment <-
function(filepath, format)
{
    RNAMultipleAlignment(.read.MultipleAlignment(filepath, format, "RNA"),
                         rowmask=as(IRanges(),"NormalIRanges"),
                         colmask=.read.MultipleMask(filepath,format))
}

readAAMultipleAlignment <-
function(filepath, format)
{
    AAMultipleAlignment(.read.MultipleAlignment(filepath, format, "AA"),
                        rowmask=as(IRanges(),"NormalIRanges"),
                        colmask=.read.MultipleMask(filepath,format))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Write functions.
###

## helper to chop up strings into pieces.
.strChop <- function(x, chopsize=10)
{
    chunks <- breakInChunks(nchar(x), chunksize=chopsize)
    vapply(seq_along(chunks),
           function(i) substr(x, start=start(chunks)[i], stop=end(chunks)[i]),
           character(1))
}

## We just have to just insert line spaces
.insertSpaces <- function(str){
  str <- .strChop(str)
  paste(str, collapse=" ")
}

write.MultAlign <- function(x, filepath, invertColMask, showRowNames,
                            hideMaskedCols){
    if(!is(x, "MultipleAlignment"))
        stop("'x' must be a MultipleAlignment object or derivative")

    ## 1st, we need to capture the colmask as a vector that can be included
    msk <- colmask(x)
    dims <- dim(x)
    if(invertColMask==FALSE){
        msk <- gaps(msk, start=1, end=dims[2])
    }
    ##If we are hiding the masked cols, then we don't care about the mask
    if(hideMaskedCols){
        hasMask <- FALSE
    } else {## If we show masked cols, drop mask before as.character()
        colmask(x) <- NULL
        if(length(msk) > 0){ hasMask<-TRUE } else { hasMask <- FALSE }
    }
    if(hasMask){ dims[1] <- dims[1]+1 }
    ## Massage to character vector
    ch <- as.character(x)
    ch <- unlist(lapply(ch, .insertSpaces))
    ## Convert mask to string format
    if(hasMask){
        mskInd <- as.integer(msk) ## index that should be masked
        mskCh <- paste(as.character(replace(rep(1,dim(x)[2]), mskInd, 0)),
                        collapse="")
        mskCh <- .insertSpaces(mskCh)
    }
    ## Split up the output into lines, but grouped into a list object
    ch <- setNames(lapply(ch, .strChop, chopsize=55), ch)
    ## Again consider mask, split, name & cat on (if needed)
    if(hasMask){
        mskCh <- .strChop(mskCh, chopsize=55)
        ch <- c(list(Mask = mskCh), ch)
    }
    ## 1) precalculate the max length of the names and then
    maxLen <- max(nchar(names(ch)))
    ## 2) make a string of that many spaces into a row
    stockSpc <- paste(rep(" ", maxLen), collapse="")
    ## 3) will need to buffer all names() to be that length.
    bufferSpacing<-function(name){
        spc <- paste(rep(" ", maxLen - nchar(name)), collapse="")
        paste(name, spc, sep="")
    }
    ## 4) append set of blank strings to list for "gapped-rows" later
    ch <- c(ch, list(rep("",length(ch[[1]]))))
    ## 5) Iterate so that all the rows are interleaved together
    output <- character(length(ch[[1]])*length(ch))
    for(i in seq_len(length(ch[[1]]))){
        for(j in seq_len(length(ch))){
            if(i==1) {
            output[j] <- paste(unlist(lapply(names(ch[j]), bufferSpacing)),
                               "   ",ch[[j]][i],sep="")
            } else {
                if(showRowNames){
                  output[(length(ch)*(i-1)) + j] <-
                    paste(unlist(lapply(names(ch[j]), bufferSpacing)),
                            "   ",ch[[j]][i],sep="")
                } else {
                    output[(length(ch)*(i-1)) + j] <- paste(stockSpc, "   ",
                                                            ch[[j]][i], sep="")
                }
            }
        }
    }
    ## drop trailing spaces
    output <-  gsub("\\s+$","", output)
    ## remove the extra end line
    output <- head(output, n=-1L)
    ## finally attach the dims
    if (hasMask) {
        ##Honestly not sure if I need a "W" here or what it means?
        output <- c(paste("", paste(c(dims, ""), collapse=" "), collapse=" "),
                    output)
    } else {
        output <- c(paste("", paste(dims, collapse=" "), collapse=" "),
                    output)
    }
    writeLines(output, filepath)
}

write.phylip <- function(x, filepath){
    write.MultAlign(x, filepath, invertColMask=TRUE, showRowNames=FALSE,
                       hideMaskedCols=FALSE)
}

