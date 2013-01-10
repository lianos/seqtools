##' Calculates the TPM for the number of reads \code{x} in a sample with a total
##' number of reads \code{y}
setGeneric("TPM", signature=c('x', 'normer'),
function(x, normer, ...) {
  standardGeneric("TPM")
})

setMethod("TPM", c(x="numeric", normer="numeric"),
function(x, normer, ...) {
  (1e6 / normer) * x
})

##' Calculates the number of read counts needed to get a TPM value of \code{x}
##' from a sample with \code{y}-many reads.
setGeneric("iTPM", signature=c('x', 'normer'),
function(x, normer, ...) {
  standardGeneric("iTPM")
})

setMethod("iTPM", c(x="numeric", normer='numeric'),
function(x, normer, ...) {
  (normer * x) / 1e6
})
