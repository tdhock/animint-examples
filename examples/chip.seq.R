requireGitHub::requireGitHub(
  "tdhock/animint@9151224296018fc0876a57acce05cf59d02e7b9e",
  "tdhock/ggplot2@98cefe4d653ce8f214177b66dc030c2f3c725ffb")

load("../data/chip.seq.RData")

## Some constants for the labels.
ann.colors <- c("noDifference"='#f6f4bf',
                "difference"="#ff7d7d")
fp.fn.map <- c(false.positive="noDifference",
               false.negative="difference")
fp.fn.colors <- ann.colors[fp.fn.map]
names(fp.fn.colors) <- names(fp.fn.map)
total.vars <- 168
set.linetypes <- c(validation="dotted", train="solid")
## nonzero pair errors, for dots on the select sample plot.
pair.errors <-
  rbind(data.frame(chip.seq$fp.pairs, error.type="false.positive"),
        data.frame(chip.seq$fn.pairs, error.type="false.negative"))
nonzero.pairs <- subset(pair.errors, errors > 0)
nonzero.pairs$error.type <-
  factor(nonzero.pairs$error.type, c("false.positive", "false.negative"))
## false positive/negative rates for the error/model complexity plot.
fp.fn <-
  rbind(data.frame(chip.seq$fp.set, error.type="false.positive"),
        data.frame(chip.seq$fn.set, error.type="false.negative"))
fp.fn$error.type <-
  factor(fp.fn$error.type, c("false.positive", "false.negative"))
fp.fn.nonzero <- subset(fp.fn, errors > 0)
unaligned <-
  list(chroms=ggplot()+
       theme_animint(width=250, height=220)+
       geom_segment(aes(0, chr.int, xend=bases/1e6, yend=chr.int),
                    data=chip.seq$chroms, color="grey")+
       geom_text(aes(0, chr.int, label=paste0("chr", chr)),
                 data=chip.seq$set.info, hjust=1)+
       ggtitle("select annotated set")+
       guides(color="none")+
       xlim(-25, 250)+
       xlab("position on chromosome (mega base pairs)")+
       theme(axis.line.y=element_blank(),
             axis.text.y=element_blank(),
             axis.title.y=element_blank(),
             axis.ticks.y=element_blank())+
       geom_point(aes((chromEnd+chromStart)/2/1e6, chr.int, 
                      clickSelects=set.name),
                  data=chip.seq$set.info, size=5),
       
       samples=ggplot()+
       scale_x_continuous("false positive/negative rate (percent)",
                          breaks=c(0, 50, 100), limits=c(-250, 100))+
       theme_animint(width=300, height=220)+
       scale_fill_manual(values=fp.fn.colors)+
       ggtitle("select samples")+
       ylab("sample")+
       theme(axis.line.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank()),
       
       error=ggplot()+
       theme_animint(height=220)+
       ggtitle("select model complexity and set")+
       xlab("model complexity -log(lambda)")+
       ylab("incorrect/total annotations (percent)")+
       geom_vline(aes(xintercept=complexity,
                      clickSelects=complexity.i),
                  data=chip.seq$models, size=15, alpha=1/2)+
       scale_color_manual(values=fp.fn.colors)+
       geom_line(aes(complexity, percent.error,
                     linetype=set.name, group=set.name,
                     clickSelects=set.name),
                 data=chip.seq$error, size=5, alpha=3/4)+
       scale_linetype_manual(values=set.linetypes)+
       geom_line(aes(complexity, percent/2, 
                     showSelected=set.name, group=error.type),
                 data=fp.fn, size=4.5, color="grey")+
       geom_line(aes(complexity, percent/2, color=error.type, size=error.type,
                     showSelected=set.name, group=error.type),
                 data=fp.fn)+
       scale_size_manual(values=c(false.negative=1.5, false.positive=3.5))+
       geom_text(aes(6, 40, 
                     label=sprintf("%d/%d=%.1f%% false positive bases",
                       errors, bases, percent),
                     showSelected=complexity.i,
                     showSelected4=set.name),
                 data=chip.seq$fp.set)+
       geom_text(aes(6, 33, 
                     label=sprintf("%d/%d=%.1f%% false negative bases",
                       errors, bases, percent),
                     showSelected=complexity.i,
                     showSelected4=set.name),
                 data=chip.seq$fn.set)+
       geom_text(aes(6, 26,
                     label=sprintf("%d/%d nonzero coefficients",
                       variables, total.vars),
                     showSelected=complexity.i),
                 data=chip.seq$nonzero),

       roc=ggplot()+
       theme_animint(width=200, height=220)+
       scale_linetype_manual(values=set.linetypes)+
       guides(size="none", linetype="none")+
       geom_path(aes(FPR, TPR, group=set.name,
                     linetype=set.name, size=set.name,
                     showSelected=complexity.i, clickSelects=set.name),
                 data=chip.seq$roc.curves)+
       geom_point(aes(FPR, TPR, showSelected=complexity.i,
                      clickSelects=set.name),
                  data=chip.seq$roc.points, color="violet", size=5)+
       scale_size_manual(values=c(train=3, validation=5))+
       ggtitle("ROC curves")+
       xlab("False positive rate")+
       ylab("True positive rate"),
       
       probability=ggplot()+
       theme_animint(width=1300, height=300)+
       ggtitle("learned difference function")+
       theme(axis.line.x=element_blank(),
             axis.text.x=element_blank(),
             axis.title.x=element_blank(),
             axis.ticks.x=element_blank())+
       geom_rect(aes(xmin=min.norm, xmax=max.norm,
                     ymin=0, ymax=1,
                     fill=annotation,
                     showSelected=set.name),
                     data=chip.seq$regions)+
       geom_text(aes(0.5, -0.05, label=paste0(width.bp, " bases on chr", chr),
                     showSelected=set.name),
                 data=chip.seq$bases)+
       geom_text(aes(0, -0.05, label=first.base,
                     showSelected=set.name),
                 data=chip.seq$bases, hjust=0)+
       geom_text(aes(1, -0.05, label=last.base,
                     showSelected=set.name),
                 data=chip.seq$bases, hjust=1)+
       scale_fill_manual(values=ann.colors)+
       geom_text(aes(mid.norm, 1.05,
                     label=sprintf("%d/%d=%.1f%% false positive bases",
                       errors, bases, percent),
                     showSelected=complexity.i,
                     showSelected2=sample1,
                     showSelected3=sample2,
                     showSelected4=set.name),
                 data=chip.seq$fp.pairs)+
       geom_text(aes(mid.norm, 1.05,
                     label=sprintf("%d/%d=%.1f%% false negative bases",
                       errors, bases, percent),
                     showSelected=complexity.i,
                     showSelected2=sample1,
                     showSelected3=sample2,
                     showSelected4=set.name),
                 data=chip.seq$fn.pairs)+
       geom_hline(yintercept=1/2, color="grey")+
       geom_ribbon(aes(mid.norm, ymin=min.prob, ymax=max.prob,
                       showSelected=sample1,
                       showSelected2=sample2,
                       showSelected3=complexity.i,
                       showSelected4=set.name),
                   chunk_vars=c("sample1","sample2","complexity.i","set.name"),
                   data=chip.seq$probability, color="blue")+
       ylab("probability of difference"),
       
       signal=ggplot()+
       theme_animint(width=1300, height=300)+
       ggtitle("ChIP-seq signal pair")+
       theme(axis.line.x=element_blank(),
             axis.text.x=element_blank(),
             axis.title.x=element_blank(),
             axis.ticks.x=element_blank())+
       ylab("<---- one signal ----- another signal ->")+
       scale_fill_manual(values=ann.colors)+
       geom_rect(aes(xmin=min.norm, xmax=max.norm,
                     ymin=-1, ymax=1,
                     fill=annotation,
                     showSelected=set.name),
                     data=chip.seq$regions)+
       geom_text(aes(0, 1, label=sprintf("%s %s max=%.1f",
                             cell.type, sample1, max),
                     showSelected=set.name,
                     showSelected2=sample1),
                 data=chip.seq$signal.max$sample1, hjust=0)+
       geom_rect(aes(xmin=normStart, xmax=normEnd,
                     ymin=0, ymax=signal.norm,
                     showSelected=set.name,
                     showSelected2=sample1),
                 data=chip.seq$signal.segments$sample1)+
       geom_text(aes(0, -1, label=sprintf("%s %s max=%.1f",
                             cell.type, sample2, max),
                     showSelected=set.name,
                     showSelected2=sample2),
                 data=chip.seq$signal.max$sample2, hjust=0)+
       geom_rect(aes(xmin=normStart, xmax=normEnd,
                     ymin=-signal.norm, ymax=0,
                     showSelected=set.name,
                     showSelected2=sample2),
                 data=chip.seq$signal.segments$sample2)+
       geom_hline(yintercept=0, color="white"),
       
       duration=list(complexity.i=2000))
## WORKS: key aesthetic in ggplot2 and key function in D3 ensures
## constancy in the select samples geom_points.

## TODO: when I click to change the model complexity.i selection, I
## expect the select samples geom_points, ROC curves, and learned
## difference function blue geom_ribbon to have a smooth transition
## duration over 2 seconds. Smooth transitions make sense here because
## nearby values of complexity.i represent similar models, thus I
## expect some continuity. However when I select another sample, I do
## not expect the blue geom_ribbon to have a smooth transition
## duration, because the functions/data in different pairs of profiles
## are not related. It should cut from the first data to the next, as
## the signal pairs do.
for(selector.name in names(chip.seq$samples)){
  sample.df <- chip.seq$samples[[selector.name]]
  sample.df$x <- -5
  y.fact <- if(selector.name=="sample1")-1 else 1
  other.name <- if(selector.name=="sample1")"sample2" else "sample1"
  sample.df$y <- sample.df$y * y.fact
  sample.id <- sample.df[[selector.name]]
  sample.df$label <- paste(sample.df$cell.type, sample.id)
  rownames(sample.df) <- sample.id
  nonzero.names <- as.character(nonzero.pairs[[selector.name]])
  nonzero.pairs$y <- sample.df[nonzero.names, "y"]
  sample.df$xmin <- 0
  sample.df$xmax <- 100
  sample.df$ymin <- sample.df$y-1/2
  sample.df$ymax <- sample.df$y+1/2
  text.df <- sample.df
  text.df$y <- text.df$y-1/2
  nonzero.pairs$key <- with(nonzero.pairs, {
    paste(sample1, sample2, error.type)
  })
  unaligned$samples <- unaligned$samples+
    geom_text(aes_string(clickSelects=selector.name, x="x", y="y",
                         showSelected="set.name", label="label"),
              data=text.df, hjust=1)+
    geom_rect(aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                         showSelected="set.name",
                         clickSelects=selector.name),
              data=sample.df, alpha=1/2)+
    geom_point(aes_string(clickSelects=other.name, x="percent", y="y",
                          key="key",
                          showSelected="set.name", showSelected2="complexity.i",
                          fill="error.type"),
               data=nonzero.pairs, color="black", size=3, alpha=0.55)
}
gg2animint(unaligned, "chip-seq-unaligned")

## TODO: the probability plot and the signal pair plot are combined
## into a single facetted plot below, which enforces that they have
## aligned x axes. However, the animint compiler/renderer do not yet
## handle facets.
chip.seq$bases$facet <- "probability"
fp.fn.pairs <- with(chip.seq, {
  rbind(data.frame(fp.pairs, error.type="false positive"),
        data.frame(fn.pairs, error.type="false negative"))
})
fp.fn.pairs$facet <- "probability"
chip.seq$probability$facet <- "probability"
fp.fn.set <- with(chip.seq, {
  rbind(data.frame(fp.set, error.type="false positive", y=40),
        data.frame(fn.set, error.type="false negative", y=33))
})
fp.fn.set$label <- with(fp.fn.set, {
  sprintf("%d/%d=%.1f%% %s bases",
          errors, bases, percent, error.type)
})
facet.data <- chip.seq
for(data.name in c("signal.max", "signal.segments")){
  for(selector.name in c("sample1", "sample2")){
    facet.data[[data.name]][[selector.name]]$facet <- "signal"
  }
}
aligned <-
  list(chroms=ggplot()+
       theme_animint(width=250, height=220)+
       geom_segment(aes(0, chr.int, xend=bases/1e6, yend=chr.int),
                    data=chip.seq$chroms, color="grey")+
       geom_text(aes(0, chr.int, label=paste0("chr", chr)),
                 data=chip.seq$set.info, hjust=1)+
       ggtitle("select annotated set")+
       guides(color="none")+
       xlim(-25, 250)+
       xlab("position on chromosome (mega base pairs)")+
       theme(axis.line.y=element_blank(),
             axis.text.y=element_blank(),
             axis.title.y=element_blank(),
             axis.ticks.y=element_blank())+
       geom_point(aes((chromEnd+chromStart)/2/1e6, chr.int, 
                      clickSelects=set.name),
                  data=chip.seq$set.info, size=5),
       
       samples=ggplot()+
       scale_x_continuous("false positive/negative rate (percent)",
                          breaks=c(0, 50, 100), limits=c(-250, 100))+
       theme_animint(width=300, height=220)+
       scale_fill_manual(values=fp.fn.colors)+
       ggtitle("select samples")+
       ylab("sample")+
       theme(axis.line.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank()),
       
       error=ggplot()+
       theme_animint(height=220)+
       ggtitle("select model complexity and set")+
       xlab("model complexity -log(lambda)")+
       ylab("incorrect/total annotations (percent)")+
       geom_vline(aes(xintercept=complexity,
                      clickSelects=complexity.i),
                  data=chip.seq$models, size=15, alpha=1/2)+
       scale_color_manual(values=fp.fn.colors)+
       geom_line(aes(complexity, percent.error,
                     linetype=set.name, group=set.name,
                     clickSelects=set.name),
                 data=chip.seq$error, size=5, alpha=3/4)+
       scale_linetype_manual(values=set.linetypes)+
       geom_line(aes(complexity, percent/2, 
                     showSelected=set.name, group=error.type),
                 data=fp.fn, size=4.5, color="grey")+
       geom_line(aes(complexity, percent/2, color=error.type, size=error.type,
                     showSelected=set.name, group=error.type),
                 data=fp.fn)+
       scale_size_manual(values=c(false.negative=1.5, false.positive=3.5))+
       geom_text(aes(6, y, label=label,
                     showSelected=complexity.i,
                     showSelected4=set.name),
                 data=fp.fn.set)+
       geom_text(aes(6, 26,
                     label=sprintf("%d/%d nonzero coefficients",
                       variables, total.vars),
                     showSelected=complexity.i),
                 data=chip.seq$nonzero),

       roc=ggplot()+
       theme_animint(width=200, height=220)+
       scale_linetype_manual(values=set.linetypes)+
       guides(size="none", linetype="none")+
       geom_path(aes(FPR, TPR, group=set.name,
                     linetype=set.name, size=set.name,
                     showSelected=complexity.i, clickSelects=set.name),
                 data=chip.seq$roc.curves)+
       geom_point(aes(FPR, TPR, showSelected=complexity.i,
                      clickSelects=set.name),
                  data=chip.seq$roc.points, color="violet", size=5)+
       scale_size_manual(values=c(train=3, validation=5))+
       ggtitle("ROC curves")+
       xlab("False positive rate")+
       ylab("True positive rate"),
       
       probSignal=ggplot()+
       theme_animint(width=1300, height=400)+
       ggtitle("ChIP-seq signal pair and learned difference function")+
       theme(axis.line.x=element_blank(),
             axis.text.x=element_blank(),
             axis.title.x=element_blank(),
             axis.ticks.x=element_blank())+
       geom_rect(aes(xmin=min.norm, xmax=max.norm,
                     ymin=ymin, ymax=1,
                     fill=annotation,
                     showSelected=set.name),
                     data=rbind(data.frame(chip.seq$regions,
                       ymin=0, facet="probability"),
                       data.frame(chip.seq$regions,
                                  ymin=-1, facet="signal")))+
       geom_text(aes(0.5, -0.05, label=paste0(width.bp, " bases on chr", chr),
                     showSelected=set.name),
                 data=chip.seq$bases)+
       geom_text(aes(0, -0.05, label=first.base,
                     showSelected=set.name),
                 data=chip.seq$bases, hjust=0)+
       geom_text(aes(1, -0.05, label=last.base,
                     showSelected=set.name),
                 data=chip.seq$bases, hjust=1)+
       scale_fill_manual(values=ann.colors)+
       geom_text(aes(mid.norm, 1.05,
                     label=sprintf("%d/%d=%.1f%% %s bases",
                       errors, bases, percent, error.type),
                     showSelected=complexity.i,
                     showSelected2=sample1,
                     showSelected3=sample2,
                     showSelected4=set.name),
                 data=fp.fn.pairs)+
       ylab("probability of difference")+
       theme_animint(width=1300, height=300)+
       geom_text(aes(0, 1, label=sprintf("%s %s max=%.1f",
                             cell.type, sample1, max),
                     showSelected=set.name,
                     showSelected2=sample1),
                 data=chip.seq$signal.max$sample1, hjust=0)+
       geom_rect(aes(xmin=normStart, xmax=normEnd,
                     ymin=0, ymax=signal.norm,
                     showSelected=set.name,
                     showSelected2=sample1),
                 data=chip.seq$signal.segments$sample1)+
       geom_text(aes(0, -1, label=sprintf("%s %s max=%.1f",
                             cell.type, sample2, max),
                     showSelected=set.name,
                     showSelected2=sample2),
                 data=chip.seq$signal.max$sample2, hjust=0)+
       geom_rect(aes(xmin=normStart, xmax=normEnd,
                     ymin=-signal.norm, ymax=0,
                     showSelected=set.name,
                     showSelected2=sample2),
                 data=chip.seq$signal.segments$sample2)+
       geom_hline(aes(yintercept=y),
                  data=data.frame(y=0, facet="signal"),
                  color="white")+
       geom_hline(aes(yintercept=y),
                  data=data.frame(y=1/2, facet="probability"),
                  color="grey")+
       facet_grid(facet ~ ., scales="free_y")+
       theme_bw()+theme(panel.margin=grid::unit(0, "cm"))+
       geom_ribbon(aes(mid.norm, ymin=min.prob, ymax=max.prob,
                       showSelected=sample1,
                       showSelected2=sample2,
                       showSelected3=complexity.i,
                       showSelected4=set.name),
                   chunk_vars=c("sample1","sample2","complexity.i","set.name"),
                   data=chip.seq$probability, color="blue"),
       
       duration=list(complexity.i=2000))
for(selector.name in names(chip.seq$samples)){
  sample.df <- chip.seq$samples[[selector.name]]
  sample.df$x <- -5
  y.fact <- if(selector.name=="sample1")-1 else 1
  other.name <- if(selector.name=="sample1")"sample2" else "sample1"
  sample.df$y <- sample.df$y * y.fact
  sample.id <- sample.df[[selector.name]]
  sample.df$label <- paste(sample.df$cell.type, sample.id)
  rownames(sample.df) <- sample.id
  nonzero.names <- as.character(nonzero.pairs[[selector.name]])
  nonzero.pairs$y <- sample.df[nonzero.names, "y"]
  sample.df$xmin <- 0
  sample.df$xmax <- 100
  sample.df$ymin <- sample.df$y-1/2
  sample.df$ymax <- sample.df$y+1/2
  text.df <- sample.df
  text.df$y <- text.df$y-1/2
  nonzero.pairs$key <- with(nonzero.pairs, {
    paste(sample1, sample2, error.type)
  })
  aligned$samples <- aligned$samples+
    geom_text(aes_string(clickSelects=selector.name, x="x", y="y",
                         showSelected="set.name", label="label"),
              data=text.df, hjust=1)+
    geom_rect(aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                         showSelected="set.name",
                         clickSelects=selector.name),
              data=sample.df, alpha=1/2)+
    geom_point(aes_string(clickSelects=other.name, x="percent", y="y",
                          key="key",
                          showSelected="set.name", showSelected2="complexity.i",
                          fill="error.type"),
               data=nonzero.pairs, color="black", size=3, alpha=0.55)
}
gg2animint(aligned, "chip-seq-aligned")

## TODO: theme_animint(updateAxes=c("x", "y")) Update the x and y axes
## after the selection changes, instead of drawing geom_text labels
## for the min/max values.
prob.df <-
  rbind(head(chip.seq$probability, 2000),
        tail(chip.seq$probability, 2000))
prob.df <- chip.seq$probability
fac <- function(x)factor(x, c("sample1", "probability", "sample2"))
prob.df$facet <- fac("probability")
updateAxes <-
  list(chroms=ggplot()+
       theme_animint(width=250, height=220)+
       geom_segment(aes(0, chr.int, xend=bases/1e6, yend=chr.int),
                    data=chip.seq$chroms, color="grey")+
       geom_text(aes(0, chr.int, label=paste0("chr", chr)),
                 data=chip.seq$set.info, hjust=1)+
       ggtitle("select annotated set")+
       guides(color="none")+
       xlim(-25, 250)+
       xlab("position on chromosome (mega base pairs)")+
       theme(axis.line.y=element_blank(),
             axis.text.y=element_blank(),
             axis.title.y=element_blank(),
             axis.ticks.y=element_blank())+
       geom_point(aes((chromEnd+chromStart)/2/1e6, chr.int, 
                      clickSelects=set.name),
                  data=chip.seq$set.info, size=5),
       
       samples=ggplot()+
       scale_x_continuous("false positive/negative rate (percent)",
                          breaks=c(0, 50, 100), limits=c(-250, 100))+
       theme_animint(width=300, height=220)+
       scale_fill_manual(values=fp.fn.colors)+
       ggtitle("select samples")+
       ylab("sample")+
       theme(axis.line.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank()),
       
       error=ggplot()+
       theme_animint(height=220)+
       ggtitle("select model complexity and set")+
       xlab("model complexity -log(lambda)")+
       ylab("incorrect/total annotations (percent)")+
       geom_vline(aes(xintercept=complexity,
                      clickSelects=complexity.i),
                  data=chip.seq$models, size=15, alpha=1/2)+
       scale_color_manual(values=fp.fn.colors)+
       geom_line(aes(complexity, percent.error,
                     linetype=set.name, group=set.name,
                     clickSelects=set.name),
                 data=chip.seq$error, size=5, alpha=3/4)+
       scale_linetype_manual(values=set.linetypes)+
       geom_line(aes(complexity, percent/2, 
                     showSelected=set.name, group=error.type),
                 data=fp.fn, size=4.5, color="grey")+
       geom_line(aes(complexity, percent/2, color=error.type, size=error.type,
                     showSelected=set.name, group=error.type),
                 data=fp.fn)+
       scale_size_manual(values=c(false.negative=1.5, false.positive=3.5))+
       geom_text(aes(6, y, label=label,
                     showSelected=complexity.i,
                     showSelected4=set.name),
                 data=fp.fn.set)+
       geom_text(aes(6, 26,
                     label=sprintf("%d/%d nonzero coefficients",
                       variables, total.vars),
                     showSelected=complexity.i),
                 data=chip.seq$nonzero),

       roc=ggplot()+
       theme_animint(width=200, height=220)+
       scale_linetype_manual(values=set.linetypes)+
       guides(size="none", linetype="none")+
       geom_path(aes(FPR, TPR, group=set.name,
                     linetype=set.name, size=set.name,
                     showSelected=complexity.i, clickSelects=set.name),
                 data=chip.seq$roc.curves)+
       geom_point(aes(FPR, TPR, showSelected=complexity.i,
                      clickSelects=set.name),
                  data=chip.seq$roc.points, color="violet", size=5)+
       scale_size_manual(values=c(train=3, validation=5))+
       ggtitle("ROC curves")+
       xlab("False positive rate")+
       ylab("True positive rate"),
       
       probSignal=ggplot()+
       theme_animint(width=1300, height=400, updateAxes=c("x", "y"))+
       ggtitle("ChIP-seq signal pair and learned difference function")+
       xlab("position on chromosome (base pairs)")+
       geom_tallrect(aes(xmin=chromStart, xmax=chromEnd,
                         fill=annotation,
                         showSelected=set.name),
                 data=chip.seq$regions)+
       scale_fill_manual(values=ann.colors)+
       ylab("probability of difference")+
       theme_animint(width=1300, height=300)+
       geom_hline(aes(yintercept=y),
                  data=data.frame(y=1/2, facet=fac("probability")),
                  color="grey")+
       ## NOTE: free_y means that the 3 panels going from top to
       ## bottom have different axes.
       theme_bw()+theme(panel.margin=grid::unit(0, "cm"))+
       facet_grid(facet ~ ., scales="free_y")+
       facet_grid(facet ~ set.name, scales="free")+
       geom_ribbon(aes(mid.base, ymin=min.prob, ymax=max.prob,
                       showSelected=sample1,
                       showSelected2=sample2,
                       showSelected3=complexity.i,
                       showSelected4=set.name),
                   chunk_vars=c("sample1","sample2","complexity.i","set.name"),
                   data=prob.df, color="blue"),
       
       duration=list(complexity.i=2000))
sig.df <- chip.seq$signal.segments[[1]]
sets <- split(sig.df, sig.df$set.name)
for(selector.name in names(chip.seq$samples)){
  ## probSignal plot geom_rect and geom_text.
  max.df <- chip.seq$signal.max[[selector.name]]
  max.df$base.min <- NA
  for(set.name in names(sets)){
    is.set <- max.df$set.name == set.name
    max.df[is.set, "base.min"] <- sets[[set.name]]$chromStart[1]
  }
  max.df$label <- paste(max.df$cell.type, max.df[[selector.name]])
  max.df$facet <- fac(selector.name)
  seg.df <- chip.seq$signal.segments[[selector.name]]
  seg.df$facet <- fac(selector.name)
  updateAxes$probSignal <- updateAxes$probSignal+
    geom_text(aes_string(x="base.min", y="max", label="label",
                  showSelected="set.name",
                  showSelected2=selector.name),
              data=max.df, hjust=0)+
    geom_rect(aes_string(xmin="chromStart", xmax="chromEnd",
                         ymin=0, ymax="signal",
                         showSelected="set.name",
                         showSelected2=selector.name),
              data=seg.df)
  ## Select samples plot.
  sample.df <- chip.seq$samples[[selector.name]]
  sample.df$x <- -5
  y.fact <- if(selector.name=="sample1")-1 else 1
  other.name <- if(selector.name=="sample1")"sample2" else "sample1"
  sample.df$y <- sample.df$y * y.fact
  sample.id <- sample.df[[selector.name]]
  sample.df$label <- paste(sample.df$cell.type, sample.id)
  rownames(sample.df) <- sample.id
  nonzero.names <- as.character(nonzero.pairs[[selector.name]])
  nonzero.pairs$y <- sample.df[nonzero.names, "y"]
  sample.df$xmin <- 0
  sample.df$xmax <- 100
  sample.df$ymin <- sample.df$y-1/2
  sample.df$ymax <- sample.df$y+1/2
  text.df <- sample.df
  text.df$y <- text.df$y-1/2
  nonzero.pairs$key <- with(nonzero.pairs, {
    paste(sample1, sample2, error.type)
  })
  updateAxes$samples <- updateAxes$samples+
    geom_text(aes_string(clickSelects=selector.name, x="x", y="y",
                         showSelected="set.name", label="label"),
              data=text.df, hjust=1)+
    geom_rect(aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                         showSelected="set.name",
                         clickSelects=selector.name),
              data=sample.df, alpha=1/2)+
    geom_point(aes_string(clickSelects=other.name, x="percent", y="y",
                          key="key",
                          showSelected="set.name", showSelected2="complexity.i",
                          fill="error.type"),
               data=nonzero.pairs, color="black", size=3, alpha=0.55)
}
print(updateAxes)
gg2animint(updateAxes, "chip-seq-updateAxes")

## WORKS: The probability plot and the signal pair plot are combined
## into a single non-facetted plot below, with the probability
## function drawn on top of the signals.
combined <-
  list(chroms=ggplot()+
       theme_animint(width=250, height=220)+
       geom_segment(aes(0, chr.int, xend=bases/1e6, yend=chr.int),
                    data=chip.seq$chroms, color="grey")+
       geom_text(aes(0, chr.int, label=paste0("chr", chr)),
                 data=chip.seq$set.info, hjust=1)+
       ggtitle("select annotated set")+
       guides(color="none")+
       xlim(-25, 250)+
       xlab("position on chromosome (mega base pairs)")+
       theme(axis.line.y=element_blank(),
             axis.text.y=element_blank(),
             axis.title.y=element_blank(),
             axis.ticks.y=element_blank())+
       geom_point(aes((chromEnd+chromStart)/2/1e6, chr.int, 
                      clickSelects=set.name),
                  data=chip.seq$set.info, size=5),
       
       samples=ggplot()+
       scale_x_continuous("false positive/negative rate (percent)",
                          breaks=c(0, 50, 100), limits=c(-250, 100))+
       theme_animint(width=300, height=220)+
       scale_fill_manual(values=fp.fn.colors)+
       ggtitle("select samples")+
       ylab("sample")+
       theme(axis.line.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank()),
       
       error=ggplot()+
       theme_animint(height=220)+
       ggtitle("select model complexity and set")+
       xlab("model complexity -log(lambda)")+
       ylab("incorrect/total annotations (percent)")+
       geom_vline(aes(xintercept=complexity,
                      clickSelects=complexity.i),
                  data=chip.seq$models, size=15, alpha=1/2)+
       scale_color_manual(values=fp.fn.colors)+
       geom_line(aes(complexity, percent.error,
                     linetype=set.name, group=set.name,
                     clickSelects=set.name),
                 data=chip.seq$error, size=5, alpha=3/4)+
       scale_linetype_manual(values=set.linetypes)+
       geom_line(aes(complexity, percent/2, 
                     showSelected=set.name, group=error.type),
                 data=fp.fn, size=4.5, color="grey")+
       geom_line(aes(complexity, percent/2, color=error.type, size=error.type,
                     showSelected=set.name, group=error.type),
                 data=fp.fn)+
       scale_size_manual(values=c(false.negative=1.5, false.positive=3.5))+
       geom_text(aes(6, y, label=label,
                     showSelected=complexity.i,
                     showSelected4=set.name),
                 data=fp.fn.set)+
       geom_text(aes(6, 26,
                     label=sprintf("%d/%d nonzero coefficients",
                       variables, total.vars),
                     showSelected=complexity.i),
                 data=chip.seq$nonzero),

       roc=ggplot()+
       theme_animint(width=200, height=220)+
       scale_linetype_manual(values=set.linetypes)+
       guides(size="none", linetype="none")+
       geom_path(aes(FPR, TPR, group=set.name,
                     linetype=set.name, size=set.name,
                     showSelected=complexity.i, clickSelects=set.name),
                 data=chip.seq$roc.curves)+
       geom_point(aes(FPR, TPR, showSelected=complexity.i,
                      clickSelects=set.name),
                  data=chip.seq$roc.points, color="violet", size=5)+
       scale_size_manual(values=c(train=3, validation=5))+
       ggtitle("ROC curves")+
       xlab("False positive rate")+
       ylab("True positive rate"),
       
       probSignal=ggplot()+
       theme_animint(width=1300, height=400)+
       ggtitle("ChIP-seq signal pair and learned difference function")+
       theme(axis.line.x=element_blank(),
             axis.text.x=element_blank(),
             axis.title.x=element_blank(),
             axis.ticks.x=element_blank())+
       geom_rect(aes(xmin=min.norm, xmax=max.norm,
                     ymin=0, ymax=1,
                     fill=annotation,
                     showSelected=set.name),
                     data=chip.seq$regions)+
       geom_text(aes(0.5, -0.15, label=paste0(width.bp, " bases on chr", chr),
                     showSelected=set.name),
                 data=chip.seq$bases)+
       geom_text(aes(0, -0.15, label=first.base,
                     showSelected=set.name),
                 data=chip.seq$bases, hjust=0)+
       geom_text(aes(1, -0.15, label=last.base,
                     showSelected=set.name),
                 data=chip.seq$bases, hjust=1)+
       scale_fill_manual(values=ann.colors)+
       geom_text(aes(mid.norm, 1.10,
                     label=sprintf("%d/%d=%.1f%% %s bases",
                       errors, bases, percent, error.type),
                     showSelected=complexity.i,
                     showSelected2=sample1,
                     showSelected3=sample2,
                     showSelected4=set.name),
                 data=fp.fn.pairs)+
       scale_y_continuous("probability of difference",
                          breaks=c(0, 1/2, 1))+
       theme_animint(width=1300, height=300)+
       geom_text(aes(0, 1, label=sprintf("%s %s max=%.1f",
                             cell.type, sample1, max),
                     showSelected=set.name,
                     showSelected2=sample1),
                 data=chip.seq$signal.max$sample1, hjust=0)+
       geom_rect(aes(xmin=normStart, xmax=normEnd,
                     ymin=1/2, ymax=(1+signal.norm)/2,
                     showSelected=set.name,
                     showSelected2=sample1),
                 data=chip.seq$signal.segments$sample1)+
       geom_text(aes(0, -0.06, label=sprintf("%s %s max=%.1f",
                             cell.type, sample2, max),
                     showSelected=set.name,
                     showSelected2=sample2),
                 data=chip.seq$signal.max$sample2, hjust=0)+
       geom_rect(aes(xmin=normStart, xmax=normEnd,
                     ymin=(1-signal.norm)/2, ymax=1/2,
                     showSelected=set.name,
                     showSelected2=sample2),
                 data=chip.seq$signal.segments$sample2)+
       geom_hline(yintercept=1/2, color="white")+
       geom_ribbon(aes(mid.norm, ymin=min.prob, ymax=max.prob,
                       showSelected=sample1,
                       showSelected2=sample2,
                       showSelected3=complexity.i,
                       showSelected4=set.name),
                   chunk_vars=c("sample1","sample2","complexity.i","set.name"),
                   data=chip.seq$probability, color="blue"),
       
       duration=list(complexity.i=2000))
for(selector.name in names(chip.seq$samples)){
  sample.df <- chip.seq$samples[[selector.name]]
  sample.df$x <- -5
  y.fact <- if(selector.name=="sample1")-1 else 1
  other.name <- if(selector.name=="sample1")"sample2" else "sample1"
  sample.df$y <- sample.df$y * y.fact
  sample.id <- sample.df[[selector.name]]
  sample.df$label <- paste(sample.df$cell.type, sample.id)
  rownames(sample.df) <- sample.id
  nonzero.names <- as.character(nonzero.pairs[[selector.name]])
  nonzero.pairs$y <- sample.df[nonzero.names, "y"]
  sample.df$xmin <- 0
  sample.df$xmax <- 100
  sample.df$ymin <- sample.df$y-1/2
  sample.df$ymax <- sample.df$y+1/2
  text.df <- sample.df
  text.df$y <- text.df$y-1/2
  nonzero.pairs$key <- with(nonzero.pairs, {
    paste(sample1, sample2, error.type)
  })
  combined$samples <- combined$samples+
    geom_text(aes_string(clickSelects=selector.name, x="x", y="y",
                         showSelected="set.name", label="label"),
              data=text.df, hjust=1)+
    geom_rect(aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                         showSelected="set.name",
                         clickSelects=selector.name),
              data=sample.df, alpha=1/2)+
    geom_point(aes_string(clickSelects=other.name, x="percent", y="y",
                          key="key",
                          showSelected="set.name", showSelected2="complexity.i",
                          fill="error.type"),
               data=nonzero.pairs, color="black", size=3, alpha=0.55)
}
gg2animint(combined, "chip-seq-combined")
