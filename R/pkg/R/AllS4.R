##' Returns the name of the aligner used
##' @param x A \code{\linkS4class{SeqStore}} object
##' @return The name of the aligner used, eg. \code{bwa}, \code{bowtie}, etc.
setGeneric("aligner", function(x, ...) standardGeneric("aligner"))


##' Query \code{\linkS4class{SeqStore}} object to see if data is from a paired end run
##'
##' @exportMethod isPaired
##'
##' @param x A \code{\linkS4class{SeqStore}} object
##' @return A logical indicating if data is paired.
setGeneric("isPaired", function(x, ...) standardGeneric("isPaired"))

setGeneric("header", function(x, ...) standardGeneric("header"))

setGeneric("swapStrand", function(x, ...) standardGeneric("swapStrand"))
