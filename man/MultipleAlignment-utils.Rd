\name{MultipleAlignment-utils}

\alias{MultipleAlignment-utils}
\alias{MultipleAlignment_utils}

\alias{consensusMatrix,MultipleAlignment-method}
\alias{consensusString,MultipleAlignment-method}
\alias{consensusString,DNAMultipleAlignment-method}
\alias{consensusString,RNAMultipleAlignment-method}
\alias{consensusString,AAMultipleAlignment-method}
\alias{consensusViews}
\alias{consensusViews,MultipleAlignment-method}
\alias{consensusViews,DNAMultipleAlignment-method}
\alias{consensusViews,RNAMultipleAlignment-method}
\alias{consensusViews,AAMultipleAlignment-method}
\alias{alphabetFrequency,MultipleAlignment-method}


\title{MultipleAlignment utilities}

\description{
   A small set of convenient utilities function to operate on a
   \link{MultipleAlignment} object.
}

\details{
  In the code snippets below, \code{x} is a \link{MultipleAlignment} object.

  \describe{
    \item{\code{consensusMatrix(x, as.prob, baseOnly)}:}{
      Creates an integer matrix containing the column frequencies of
      the underlying alphabet with masked columns being represented
      with \code{NA} values. If \code{as.prob} is \code{TRUE}, then
      probabilities are reported, otherwise counts are reported (the
      default). If \code{baseOnly} is \code{TRUE}, then the non-base
      letters are collapsed into an \code{"other"} category.
    }
    \item{\code{consensusString(x, ...)}:}{
      Creates a consensus string for \code{x} with the symbol \code{"#"}
      representing a masked column. See \code{\link{consensusString}}
      for details on the arguments.
    }
    \item{\code{consensusViews(x, ...)}:}{
      Similar to the \code{consensusString} method. It returns a 
      \link{XStringViews} on the consensus string containing subsequence
      contigs of non-masked columns. Unlike the \code{consensusString}
      method, the masked columns in the underlying string contain a
      consensus value rather than the \code{"#"} symbol.
    }
    \item{\code{alphabetFrequency(x, as.prob, collapse)}:}{
      Creates an integer matrix containing the row frequencies of
      the underlying alphabet. If \code{as.prob} is \code{TRUE}, then
      probabilities are reported, otherwise counts are reported (the
      default). If \code{collapse} is \code{TRUE}, then returns the
      overall frequency instead of the frequency by row.
    }
  }
}

\value{
  See description of each method above for what they return.
}

\author{P. Aboyoun and M. Carlson}

\seealso{
  \link{MultipleAlignment}
}

\examples{
example(MultipleAlignment)  # make MultipleAlignment object 'autoMasked'
autoMasked

## calculate frequencies
alphabetFrequency(autoMasked)
consensusMatrix(autoMasked, baseOnly=TRUE)[, 84:90]

## get consensus values
consensusString(autoMasked)
consensusViews(autoMasked)
}

\keyword{utilities}
\keyword{manip}
