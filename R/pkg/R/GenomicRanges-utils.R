## Thanks to Tim Triche, Jr. 
setMethod("$", "GenomicRanges", function(x, name) {
  elementMetadata(x)[, name]
})

setMethod("$<-", "GenomicRanges", function(x, name, value) {
  elementMetadata(x)[[name]] <- value
  x
})

setMethod("[[", "GenomicRanges", function(x, i, j, ...) {
  switch(i, strand=strand(x), seqnames=seqnames(x), ranges=ranges(x),
         values(x)[[i, j, ...]])
})
