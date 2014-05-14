load("~/projects/scaffold-viz/insert.list.coverage.RData")

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

scaffolds <- list()
for(data.name in names(insert.list.coverage)){
  scaffolds[[data.name]] <- toDF(insert.list.coverage[[data.name]])
}

save(scaffolds, file="../data/scaffolds.RData", compress="xz")
