load("~/projects/chip-seq/features.two.chunks.figure.RData")

toDF <- function(x){
  if(is.data.frame(x)){
    return(as.data.frame(x))
  }
  if(is.list(x)){
    return(lapply(x, toDF))
  }
  str(x)
  stop("do not know how to convert object to data.frame")
}

chip.seq <- list()
for(data.name in names(features.two.chunks.figure)){
  chip.seq[[data.name]] <- toDF(features.two.chunks.figure[[data.name]])
}

save(chip.seq, file="../data/chip.seq.RData", compress="xz")
