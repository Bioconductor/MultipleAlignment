\name{MultipleAlignment-IO}

\alias{MultipleAlignment-IO}
\alias{MultipleAlignment_IO}

% Read functions:
\alias{readDNAMultipleAlignment}
\alias{readRNAMultipleAlignment}
\alias{readAAMultipleAlignment}

% Write functions:
\alias{write.phylip}

\title{Read/write MultipleAlignment objects}

\description{
   Functions to read/write \link{MultipleAlignment} objects.
}

\usage{
## Read functions:
readDNAMultipleAlignment(filepath, format)
readRNAMultipleAlignment(filepath, format)
readAAMultipleAlignment(filepath, format)

## Write funtions:
write.phylip(x, filepath)
}

\arguments{
  \item{x}{
    A \link{MultipleAlignment} object
  }
  \item{filepath}{
    A character vector (of arbitrary length when reading, of length 1
    when writing) containing the paths to the files to read or write.
    Note that special values like \code{""} or \code{"|cmd"} (typically
    supported by other I/O functions in R) are not supported here. Also
    \code{filepath} cannot be a connection.
  }
  \item{format}{
    Either \code{"fasta"} (the default), \code{"stockholm"}, \code{"phylip"}, or
    \code{"clustal"}.
  }
}

\author{P. Aboyoun and M. Carlson}

\seealso{
  \link{MultipleAlignment}
}

\examples{
example(MultipleAlignment)  # make MultipleAlignment object 'autoMasked'
autoMasked

## write out the alignment object (with current masks) to Phylip format
write.phylip(autoMasked, filepath=tempfile())
}

\keyword{utilities}
\keyword{manip}
