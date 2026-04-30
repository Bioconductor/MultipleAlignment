strings_DNAMultipleAlignment <- function()
{
    c(string1 = "AAGGTCTCCA-GCCTGCCCTTCAGTGTGGAGGCGCTCATG--TCGGACA",
      string2 = "AAGGTCTCCA-GCCTGCCCTTCAGCGTGGAGGCGCTCATG--TCCGACA",
      string3 = "CATTTATATATGGTCCCCCTCCCCCCAAGAAACACACATAGTTTTGACA")
}

make_DNAMultipleAlignment <- function()
{
    DNAMultipleAlignment(strings_DNAMultipleAlignment())
}

test_that("DNAMultipleAlignment construction", {

    ## --- Empty DNAMultipleAlignment object ---

    malign <- DNAMultipleAlignment()
    expect_true(validObject(malign, test=TRUE, complete=TRUE))
    expect_identical(as.character(unmasked(malign)), as.character(DNAStringSet()))
    expect_identical(rownames(malign), NULL)
    expect_identical(rowmask(malign), new("NormalIRanges"))
    expect_identical(colmask(malign), new("NormalIRanges"))
    expect_identical(as.character(maskMotif(malign, "GC")), character())
    expect_identical(as.character(maskGaps(malign)), character())
    expect_identical(nrow(malign), 0L)
    expect_identical(ncol(malign), 0L)
    expect_identical(dim(malign), c(0L, 0L))
    expect_identical(maskednrow(malign), 0L)
    expect_identical(maskedncol(malign), 0L)
    expect_identical(maskeddim(malign), c(0L, 0L))
    expect_identical(maskedratio(malign), c(0L/0L, 0L/0L))
    expect_identical(nchar(malign), 0L)
    expect_identical(seqtype(malign), "DNA")
    expect_identical(as.character(malign), character(0))

    ## TODO: Move these tests to test-MultipleAlignment-utils.R
    expect_identical(consensusMatrix(malign),
                     matrix(integer(), nrow=length(DNA_ALPHABET), ncol=0,
                            dimnames=list(DNA_ALPHABET, NULL)))
    expect_identical(consensusString(malign), character())
    expect_identical(as.character(consensusViews(malign)), character())
    expect_identical(alphabetFrequency(malign),
                     matrix(integer(), nrow=0, ncol=length(DNA_ALPHABET),
                            dimnames=list(NULL, DNA_ALPHABET)))
    expect_identical(alphabetFrequency(malign, collapse=TRUE),
                     structure(integer(length(DNA_ALPHABET)), names=DNA_ALPHABET))

    ## --- Unnamed DNAMultipleAlignment object ---

    malign <- make_DNAMultipleAlignment()
    rownames(malign) <- NULL
    expect_true(validObject(malign, test=TRUE))
    expect_identical(as.character(unmasked(malign)),
                     unname(strings_DNAMultipleAlignment()))
    expect_identical(rownames(malign), NULL)
    expect_identical(rowmask(malign), new("NormalIRanges"))
    expect_identical(colmask(malign), new("NormalIRanges"))
    expect_identical(as.character(maskMotif(malign, "GC", fixed=FALSE)),
                     c("AAGGTCTCCA-TCCTTCAGTGGAGTCATG--TCGGACA",
                       "AAGGTCTCCA-TCCTTCAGTGGAGTCATG--TCCGACA",
                       "CATTTATATATCCCTCCCCAAGAAACATAGTTTTGACA"))
    expect_identical(as.character(maskGaps(malign)),
                     unname(strings_DNAMultipleAlignment()))
    expect_identical(nrow(malign), length(strings_DNAMultipleAlignment()))
    expect_identical(ncol(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    expect_identical(dim(malign),
                     c(length(strings_DNAMultipleAlignment()),
                       nchar(strings_DNAMultipleAlignment())[[1]]))
    expect_identical(maskednrow(malign), 0L)
    expect_identical(maskedncol(malign), 0L)
    expect_identical(maskeddim(malign), c(0L, 0L))
    expect_identical(maskedratio(malign), c(0, 0))
    expect_identical(nchar(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    expect_identical(seqtype(malign), "DNA")
    expect_identical(as.character(malign), unname(strings_DNAMultipleAlignment()))

    ## TODO: Move these tests to test-MultipleAlignment-utils.R
    expect_identical(consensusMatrix(malign)[1:4, 1:4],
                     rbind(A=c(2L,3L,0L,0L),
                           C=c(1L,0L,0L,0L),
                           G=c(0L,0L,2L,2L),
                           T=c(0L,0L,1L,1L)))
    expect_identical(consensusString(malign),
                     "MAKKTMTMYA-GSYYSCCCTYCMSYSWRGARRCRCWCATR--TYBGACA")
    expect_identical(as.character(consensusViews(malign)),
                     "MAKKTMTMYA-GSYYSCCCTYCMSYSWRGARRCRCWCATR--TYBGACA")
    expect_identical(alphabetFrequency(malign)[,1:4],
                     cbind(A=c(8L,8L,15L),
                           C=c(14L,16L,16L),
                           G=c(14L,13L,5L),
                           T=c(10L,9L,13L)))
    expect_identical(alphabetFrequency(malign, collapse=TRUE)[1:4],
                     c(A=31L, C=46L, G=32L, T=32L))

    ## --- Named DNAMultipleAlignment object ---

    malign <- make_DNAMultipleAlignment()
    expect_true(validObject(malign, test=TRUE))
    expect_identical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    expect_identical(rownames(malign), names(strings_DNAMultipleAlignment()))
    expect_identical(rowmask(malign), new("NormalIRanges"))
    expect_identical(colmask(malign), new("NormalIRanges"))
    expect_identical(as.character(maskMotif(malign, "GC", fixed=FALSE)),
                     c(string1="AAGGTCTCCA-TCCTTCAGTGGAGTCATG--TCGGACA",
                       string2="AAGGTCTCCA-TCCTTCAGTGGAGTCATG--TCCGACA",
                       string3="CATTTATATATCCCTCCCCAAGAAACATAGTTTTGACA"))
    expect_identical(as.character(maskGaps(malign)),
                     strings_DNAMultipleAlignment())
    expect_identical(nrow(malign), length(strings_DNAMultipleAlignment()))
    expect_identical(ncol(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    expect_identical(dim(malign),
                     c(length(strings_DNAMultipleAlignment()),
                       nchar(strings_DNAMultipleAlignment())[[1]]))
    expect_identical(maskednrow(malign), 0L)
    expect_identical(maskedncol(malign), 0L)
    expect_identical(maskeddim(malign), c(0L, 0L))
    expect_identical(maskedratio(malign), c(0, 0))
    expect_identical(nchar(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    expect_identical(seqtype(malign), "DNA")
    expect_identical(as.character(malign), strings_DNAMultipleAlignment())

    ## TODO: Move these tests to test-MultipleAlignment-utils.R
    expect_identical(consensusMatrix(malign)[1:4, 1:4],
                     rbind(A=c(2L,3L,0L,0L),
                           C=c(1L,0L,0L,0L),
                           G=c(0L,0L,2L,2L),
                           T=c(0L,0L,1L,1L)))
    expect_identical(consensusString(malign),
                     "MAKKTMTMYA-GSYYSCCCTYCMSYSWRGARRCRCWCATR--TYBGACA")
    expect_identical(as.character(consensusViews(malign)),
                     "MAKKTMTMYA-GSYYSCCCTYCMSYSWRGARRCRCWCATR--TYBGACA")
    expect_identical(alphabetFrequency(malign)[,1:4],
                     cbind(A=c(8L,8L,15L),
                           C=c(14L,16L,16L),
                           G=c(14L,13L,5L),
                           T=c(10L,9L,13L)))
    expect_identical(alphabetFrequency(malign, collapse=TRUE)[1:4],
                     c(A=31L, C=46L, G=32L, T=32L))
})

### The tests below are the original RUnit-style tests that were copied as-is
### from Biostrings. They've been broken for years but unfortunately the
### original authors/maintainers of the MultipleAlignment code have left the
### project.
### They need to be FIXED and migrated to testthat.

### FIXME!!
BROKEN_test_DNAMultipleAlignment_mask_some <- function()
{
    malign <- make_DNAMultipleAlignment()
    rowmask(malign) <- IRanges(2,2)
    colmask(malign) <- IRanges(c(1,21,43), c(10,35,49))
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    checkIdentical(rownames(malign), names(strings_DNAMultipleAlignment()))
    checkIdentical(rowmask(malign), asNormalIRanges(IRanges(2,2)))
    checkIdentical(colmask(malign), asNormalIRanges(IRanges(c(1,21,43), c(10,35,49))))
    checkIdentical(as.character(maskMotif(malign, "GC", fixed=FALSE)),
                   c(string1="-TCCTTCATG--", string3="TCCCTACATAGT"))
    checkIdentical(as.character(maskGaps(malign, min.block.width=1)),
                   c(string1="GCCTGCCCTTCATG", string3="GGTCCCCCTACATA"))
    checkIdentical(nrow(malign), length(strings_DNAMultipleAlignment()))
    checkIdentical(ncol(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    checkIdentical(dim(malign),
                   c(length(strings_DNAMultipleAlignment()),
                     nchar(strings_DNAMultipleAlignment())[[1]]))
    checkIdentical(maskednrow(malign), 1L)
    checkIdentical(maskedncol(malign), 32L)
    checkIdentical(maskeddim(malign), c(1L, 32L))
    checkIdentical(maskedratio(malign), c(1/3, 32/49))
    checkIdentical(nchar(malign), 17L)
    checkIdentical(seqtype(malign), "DNA")
    checkIdentical(as.character(malign),
                   c(string1="-GCCTGCCCTTCATG--", string3="TGGTCCCCCTACATAGT"))
    checkIdentical(consensusMatrix(malign)[1:4, 1:4],
                   rbind(A=rep(NA_integer_,4),
                         C=rep(NA_integer_,4),
                         G=rep(NA_integer_,4),
                         T=rep(NA_integer_,4)))
    checkIdentical(consensusString(malign),
                   "##########TGSYYSCCCT###############WCATRGT#######")
    checkIdentical(as.character(consensusViews(malign)),
                   c("TGSYYSCCCT", "WCATRGT"))
    checkIdentical(alphabetFrequency(malign)[,1:4],
                   cbind(A=c(1L,NA,3L),
                         C=c(6L,NA,6L),
                         G=c(3L,NA,3L),
                         T=c(4L,NA,5L)))
    checkIdentical(alphabetFrequency(malign, collapse=TRUE)[1:4],
                   c(A=4L, C=12L, G=6L, T=9L))
}

### FIXME!!
BROKEN_test_DNAMultipleAlignment_mask_all_rows <- function()
{
    malign <- make_DNAMultipleAlignment()
    rowmask(malign) <- IRanges(1,3)
    colmask(malign) <- IRanges(c(1,21,43), c(10,35,49))
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    checkIdentical(rownames(malign), names(strings_DNAMultipleAlignment()))
    checkIdentical(rowmask(malign), asNormalIRanges(IRanges(1,3)))
    checkIdentical(colmask(malign), asNormalIRanges(IRanges(c(1,21,43), c(10,35,49))))
    checkIdentical(as.character(maskMotif(malign, "GC", fixed=FALSE)), character())
    checkIdentical(as.character(maskGaps(malign)), character())
    checkIdentical(nrow(malign), length(strings_DNAMultipleAlignment()))
    checkIdentical(ncol(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    checkIdentical(dim(malign),
                   c(length(strings_DNAMultipleAlignment()),
                     nchar(strings_DNAMultipleAlignment())[[1]]))
    checkIdentical(maskednrow(malign), 3L)
    checkIdentical(maskedncol(malign), 32L)
    checkIdentical(maskeddim(malign), c(3L, 32L))
    checkIdentical(maskedratio(malign), c(3/3, 32/49))
    checkIdentical(nchar(malign), 17L)
    checkIdentical(seqtype(malign), "DNA")
    checkIdentical(as.character(malign), character(0))
    checkIdentical(consensusMatrix(malign)[1:4, 1:4],
                   rbind(A=rep(NA_integer_,4),
                         C=rep(NA_integer_,4),
                         G=rep(NA_integer_,4),
                         T=rep(NA_integer_,4)))
    checkIdentical(consensusString(malign),
                   "#################################################")
    checkIdentical(alphabetFrequency(malign)[,1:4],
                   cbind(A=rep(NA_integer_, 3),
                         C=rep(NA_integer_, 3),
                         G=rep(NA_integer_, 3),
                         T=rep(NA_integer_, 3)))
    checkIdentical(alphabetFrequency(malign, collapse=TRUE)[1:4],
                   c(A=0L, C=0L, G=0L, T=0L))
}

