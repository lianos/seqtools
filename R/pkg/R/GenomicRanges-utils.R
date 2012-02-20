## Thanks to Tim Triche, Jr. 
setMethod("$", "GRanges", function(x, name) {
  elementMetadata(x)[, name]
})

setMethod("$<-", "GRanges", function(x, name, value) {
  elementMetadata(x)[[name]] <- value
  x
})
