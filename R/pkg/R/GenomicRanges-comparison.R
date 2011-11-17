## setMethod("order", "GenomicRanges",
##     function(..., ignore.strand=FALSE, na.last=TRUE, decreasing=FALSE) {
##       if (!isTRUEorFALSE(decreasing)) {
##         stop("'decreasing' must be TRUE or FALSE")
##       }
##       args <- list(...)
##       if (length(args) != 1L) {
##         stop("ordering more than one GRanges is undefined")
##       }
##       x <- args[[1L]]
##       seqnames.order <- seqlevels(x)

##       ordered <- integer(length(x))
##       iranges <- ranges(x)
##       sofar <- 0L

##       ## Logical operations over the vectors are noticably faster than the
##       ## same over Rle's returned from seqnames() and strand()
##       all.seqnames <- as.vector(seqnames(x))
##       is.seqname <- lapply(seqnames.order, function(seqname) {
##         all.seqnames == seqname
##       })

##       if (ignore.strand) {
##         for (which.seq in is.seqname) {
##           take <- which(which.seq)
##           n <- length(take)
##           if (n > 0L) {
##             o <- order(iranges[take], decreasing=decreasing)
##             ordered[seq_along(o) + sofar] <- take[o]
##             sofar <- sofar + n
##           }
##         }
##       } else {
##         all.strands <- as.vector(strand(x))
##         is.strand <- lapply(levels(strand()), function(.strand) {
##           all.strands == .strand
##         })

##         for (which.seq in is.seqname) {
##           for (which.strand in is.strand) {
##             take <- which(which.seq & which.strand)
##             n <- length(take)
##             if (n > 0L) {
##               o <- order(iranges[take], decreasing=decreasing)
##               ordered[seq_along(o) + sofar] <- take[o]
##               sofar <- sofar + n
##             }
##           }
##         }
##       }

##       ordered
##     })

## setMethod("sort", "GenomicRanges",
##     function(x, decreasing=FALSE, ...) {
##       x[order(x, decreasing=decreasing)]
##     })

