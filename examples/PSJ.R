library(animint)

load("../data/PSJ.RData")

res.error <- PSJ$error.total.chunk

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

## prob.regions are the black segments that show which regions are
## mapped to which segmentation problems.
all.regions <- do.call(rbind, PSJ$regions.by.problem)
prob.regions.names <-
  c("bases.per.problem", "problem.i", "problem.name",
    "chromStart", "chromEnd")
prob.regions <- unique(data.frame(all.regions)[, prob.regions.names])
prob.regions$sample.id <- "problems"

all.modelSelection <- do.call(rbind, PSJ$modelSelection.by.problem)
modelSelection.errors <-
  all.modelSelection[!is.na(all.modelSelection$errors), ]
penalty.range <-
  with(all.modelSelection, c(min(max.log.lambda), max(min.log.lambda)))
penalty.mid <- mean(penalty.range)

coverage.counts <- table(PSJ$coverage$sample.id)
facet.rows <- length(coverage.counts)+1
dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval-1)/2 + 1
res.error$min.bases.per.problem <- res.error$bases.per.problem/dval2
res.error$max.bases.per.problem <- res.error$bases.per.problem*dval2

modelSelection.labels <- unique(with(all.modelSelection, {
  data.frame(problem.name=problem.name,
             bases.per.problem=bases.per.problem,
             problemStart=problemStart,
             problemEnd=problemEnd,
             min.log.lambda=penalty.mid,
             peaks=max(peaks)+0.5)
}))

cat("constructing data viz\n")
print(system.time({
  viz <-
    list(coverage=ggplot()+
           geom_segment(aes(chromStart/1e3, problem.i,
                            xend=chromEnd/1e3, yend=problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name),
                        data=prob.regions)+
           ggtitle("select problem")+
           geom_text(aes(chromStart/1e3, problem.i,
                         showSelected=bases.per.problem,
                         label=sprintf("%d problems mean size %.1f kb",
                           problems, mean.bases/1e3)),
                     data=PSJ$problem.labels,
                     hjust=0)+
           geom_segment(aes(problemStart/1e3, problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name,
                            xend=problemEnd/1e3, yend=problem.i),
                        size=5,
                        data=PSJ$problems)+
           scale_y_continuous("aligned read coverage",
                              breaks=function(limits){
                                floor(limits[2])
                              })+
           scale_linetype_manual("error type",
                                 limits=c("correct", 
                                   "false negative",
                                   "false positive"
                                          ),
                                 values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
           scale_x_continuous(paste("position on chr11",
                                    "(kilo bases = kb)"))+
           coord_cartesian(xlim=c(118167.406, 118238.833))+
           geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                             fill=annotation),
                         alpha=0.5,
                         color="grey",
                         data=PSJ$filled.regions)+
           scale_fill_manual(values=ann.colors)+
           theme_bw()+
           theme_animint(width=1500, height=facet.rows*100)+
           theme(panel.margin=grid::unit(0, "cm"))+
           facet_grid(sample.id ~ ., labeller=function(var, val){
             sub("McGill0", "", sub(" ", "\n", val))
           }, scales="free")+
           geom_line(aes(base/1e3, count),
                     data=PSJ$coverage,
                     color="grey50"),

         resError=ggplot()+
           ggtitle("select problem size")+
           ylab("minimum percent incorrect regions")+
           geom_tallrect(aes(xmin=min.bases.per.problem,
                             xmax=max.bases.per.problem,
                             clickSelects=bases.per.problem),
                         alpha=0.5,
                         data=res.error)+
           scale_x_log10()+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(res.error, chunks="this"))+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(PSJ$error.total.all, chunks="all")),

         modelSelection=ggplot()+
           geom_segment(aes(min.log.lambda, peaks,
                            xend=max.log.lambda, yend=peaks,
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(all.modelSelection, what="peaks"),
                        size=5)+
           geom_text(aes(min.log.lambda, peaks,
                         showSelected=problem.name,
                         showSelected2=bases.per.problem,
                         label=sprintf("%.1f kb in problem %s",
                           (problemEnd-problemStart)/1e3, problem.name)),
                     data=data.frame(modelSelection.labels, what="peaks"))+
           geom_segment(aes(min.log.lambda, as.integer(errors),
                            xend=max.log.lambda, yend=as.integer(errors),
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(modelSelection.errors, what="errors"),
                        size=5)+
           ggtitle("select number of samples with 1 peak")+
           ylab("")+
           facet_grid(what ~ ., scales="free"),
         
         title="PeakSegJoint model of 4 H3K4me3 ChIP-seq samples",

         first=PSJ$first)

  ## For every problem there is a selector (called problem.dot) for the
  ## number of peaks in that problem. So in this for loop we add a few
  ## layers with aes_string(clickSelects=problem.dot) or
  ## aes_string(showSelected=problem.dot) to the coverage and
  ## modelSelection plots.
  for(problem.dot in names(PSJ$modelSelection.by.problem)){
    regions.dt <- PSJ$regions.by.problem[[problem.dot]]
    regions.dt[[problem.dot]] <- regions.dt$peaks
    if(!is.null(regions.dt)){
      viz$coverage <- viz$coverage+
        geom_tallrect(aes_string(xmin="chromStart/1e3",
                                 xmax="chromEnd/1e3",
                                 linetype="status",
                                 showSelected=problem.dot,
                                 showSelected2="bases.per.problem"),
                      data=data.frame(regions.dt),
                      fill=NA,
                      color="black")
    }
    if(problem.dot %in% names(PSJ$peaks.by.problem)){
      peaks <- PSJ$peaks.by.problem[[problem.dot]]
      peaks[[problem.dot]] <- peaks$peaks
      prob.peaks.names <-
        c("bases.per.problem", "problem.i", "problem.name",
          "chromStart", "chromEnd", problem.dot)
      prob.peaks <- unique(data.frame(peaks)[, prob.peaks.names])
      prob.peaks$sample.id <- "problems"
      viz$coverage <- viz$coverage +
        geom_segment(aes_string("chromStart/1e3", "0",
                                xend="chromEnd/1e3", yend="0",
                                clickSelects="problem.name",
                                showSelected=problem.dot,
                                showSelected2="bases.per.problem"),
                     data=peaks, size=7, color="deepskyblue")+
        geom_segment(aes_string("chromStart/1e3", "problem.i",
                                xend="chromEnd/1e3", yend="problem.i",
                                clickSelects="problem.name",
                                showSelected=problem.dot,
                                showSelected2="bases.per.problem"),
                     data=prob.peaks, size=7, color="deepskyblue")
    }
    modelSelection.dt <- PSJ$modelSelection.by.problem[[problem.dot]]
    modelSelection.dt[[problem.dot]] <- modelSelection.dt$peaks
    viz$modelSelection <- viz$modelSelection+
      geom_tallrect(aes_string(xmin="min.log.lambda", 
                               xmax="max.log.lambda", 
                               clickSelects=problem.dot,
                               showSelected="problem.name",
                               showSelected2="bases.per.problem"),
                    data=modelSelection.dt, alpha=0.5)
  }
}))

cat("compiling data viz\n")
print(system.time({
  animint2dir(viz, out.dir="PSJ")
}))

