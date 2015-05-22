require(data.table)

## adapted from https://github.com/tdhock/PeakSegJoint/blob/d0def4825863f248a85ca08db10e1960814c093d/exec/Step4v-viz-train-errors.R

argv <- "~/exampleData/PeakSegJoint-chunks/1"

chunk.dir <- argv
chunks.dir <- dirname(chunk.dir)
trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
problems.RData <- file.path(chunk.dir, "problems.RData")
regions.RData <- file.path(chunk.dir, "regions.RData")

robjs <- load(regions.RData)
regions$region.i <- 1:nrow(regions)
chrom <- paste(regions$chrom[1])

## Filter counts by limits (area around regions).
limits <- regions[, .(chromStart=min(chromStart),
                      chromEnd=max(chromEnd))]
limits[, bases := chromEnd - chromStart ]
limits[, expand := as.integer(bases/10) ]
limits[, min := chromStart - expand]
limits[, max := chromEnd + expand]
lim.vec <- with(limits, c(min, max))/1e3

mobjs <- load(trained.model.RData)
pobjs <- load(problems.RData)

## Some fake data to demonstrate how to sub-sample a coverage
## data.frame using approx.
fake.segs <- data.frame(chromStart=c(0, 1, 3),
                   chromEnd=c(1, 3, 10),
                   coverage=1:3)
x <- 1:10
y <- with(fake.segs, {
  approx(chromStart+1, coverage, x,
         method="constant", rule=2)
})$y
fake.points <- data.frame(x,y)
## ggplot()+
##   geom_segment(aes(chromStart+0.5, coverage,
##                    xend=chromEnd+0.5, yend=coverage),
##                data=fake.segs)+
##   geom_point(aes(x, y), data=fake.points)

## Read and subsample count data to plot width.
width.pixels <- 1500
counts.RData.vec <- Sys.glob(file.path(chunk.dir, "*", "*.RData"))
counts.by.sample <- list()
for(counts.RData.path in counts.RData.vec){
  objs <- load(counts.RData.path)
  sample.id <- sub(".RData$", "", basename(counts.RData.path))
  a.data <- with(counts, {
    approx(chromStart+1, count,
           as.integer(seq(limits$min, limits$max, l=width.pixels)),
           method="constant", rule=2)
  })
  counts.by.sample[[sample.id]] <- with(a.data, {
    data.table(sample.id,
               base=x,
               count=y)
  })
}
some.counts <- do.call(rbind, counts.by.sample)

step2.problems.by.res <- list()
prob.labels.by.res <- list()
modelSelection.by.problem <- list()
peaks.by.problem <- list()
first.selection.list <-
  list(bases.per.problem=train.errors.picked$bases.per.problem)
regions.by.problem <- list()
for(res.str in names(step2.data.list)){
  res.data <- step2.data.list[[res.str]]
  step2.problems.by.res[[res.str]] <- 
    data.table(sample.id="problems", res.data$problems)
  bases.vec <- with(res.data$problems, problemEnd-problemStart)
  bases.per.problem <- as.integer(res.str)
  prob.labels.by.res[[res.str]] <-
    data.table(sample.id="problems",
               bases.per.problem,
               mean.bases=as.integer(mean(bases.vec)),
               problems=nrow(res.data$problems))
  for(problem.i in 1:nrow(res.data$problems)){
    problem <- res.data$problems[problem.i, ]
    problem.name <- paste(problem$problem.name)
    problem.dot <- paste0(gsub("[:-]", ".", problem.name), "peaks")
    error <- step2.error.list[[problem.name]]
    model <- step2.model.list[[problem.name]]
    if(!is.null(model)){
      if(is.data.frame(model$peaks)){
        peaks.by.problem[[problem.dot]] <- data.table(problem, model$peaks)
      }
      if(is.null(error$peaks)){
        first.selection.list[[problem.dot]] <- 0
      }else{
        ms <- 
        first.selection.list[[problem.dot]] <- error$peaks$peaks[1]
      }
      
      ms <- if(is.list(error$problem$error.regions)){
        regions.by.peaks <- list()
        for(peaks.str in names(error$problem$error.regions)){
          regions.df <- error$problem$error.regions[[peaks.str]]
          regions.df$peaks <- as.integer(peaks.str)
          regions.by.peaks[[peaks.str]] <-
            data.table(problem, regions.df)
        }
        regions.by.problem[[problem.dot]] <- do.call(rbind, regions.by.peaks)
        error$problem$modelSelection
      }else{
        data.frame(model$modelSelection,
                   errors=NA)
      }      
      modelSelection.by.problem[[problem.dot]] <-
        data.table(problem, ms)
    }
  }
}
step2.problems <- do.call(rbind, step2.problems.by.res)
prob.labels <- do.call(rbind, prob.labels.by.res)
prob.labels$problem.i <- max(prob.labels$problems)
prob.labels$chromStart <- limits$chromStart

PSJ <-
  list(problem.labels=prob.labels,
       problems=step2.problems,
       filled.regions=regions,
       coverage=some.counts,
       first=first.selection.list,
       modelSelection.by.problem=modelSelection.by.problem,
       regions.by.problem=regions.by.problem,
       peaks.by.problem=peaks.by.problem,
       error.total.chunk=res.error,
       error.total.all=train.errors)

convert <- function(L.or.DF){
  cat(length(L.or.DF), class(L.or.DF), "\n", sep=" ")
  if(is.data.frame(L.or.DF)){
    data.frame(L.or.DF)
  }else if(is.list(L.or.DF)){
    lapply(L.or.DF, convert)
  }else L.or.DF
}

PSJ <- convert(PSJ)

save(PSJ, file="../data/PSJ.RData")

