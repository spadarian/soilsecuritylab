#' Methods for FuzzyClusterGroup class
#' 
#' Common methods implemented for the FuzzyClusterGroup class.
#' 
#' @param x FuzzyClusterGroup object
#' @param ... arguments to be passed to subsequent methods, e.g., plot.FuzzyCluster for the plot method 
#' @method plot FuzzyClusterGroup
#' @name FuzzyClusterGroup-methods

setMethod('plot','FuzzyClusterGroup', function(x,...){
  plot_ <- lapply(x@clusters,function(z){
    c(x=z@classes,PICP=z@pred_int$PICP,MPI=z@pred_int$MPI)
  })
  plot_ <- do.call(rbind,plot_)
  if(dim(plot_)[2] != 3) stop('Nothing to plot. Try running prediction interval for this FuzzyClusterGroup')
  plot_ <- as.data.frame(plot_)
  .conf=x@confidence
  opt <- plot_[abs(plot_$PICP-.conf*100)<=0.1,]
  if (nrow(opt) == 0) opt <- plot_[plot_$PICP == plot_[which.min(abs(plot_$PICP-.conf*100)),]$PICP,]
  opt <- opt[which.min(opt$MPI),]
  par(mar=c(5,4,4,5)+.1)
  plot(plot_$x,plot_$PICP,type='o',pch=15,xlab='Cluster',ylab='PICP (%)',lab=c(length(plot_$x),5,7),bg='red',...)
  if (any(rownames(installed.packages()) %in% 'plotrix')) plotrix::draw.circle(opt$x,opt$PICP,c(0.2,0.3),lty=3) else warning('Install package plotrix to plot optimum number of clusters')
  abline(h=.conf*100,lty=2)
  par(new=T)
  plot(plot_$x,plot_$MPI,type="o",pch=17,col="grey50",xaxt="n",yaxt="n",xlab="",ylab="",bg='blue',...)
  if (any(rownames(installed.packages()) %in% 'plotrix')) plotrix::draw.circle(opt$x,opt$MPI,c(0.2,0.3),lty=3)
  abline(v=opt$x,lty=3)
  axis(4)
  mtext("MPI",side=4,line=3)
  par(xpd=NA)
  legend('top',col=c("black","grey50",'black'),lty=c(1,1,2),pch=c(15,17,NA),legend=c("PICP","MPI",expression(1-alpha)),ncol=3,bty='n',inset=-0.12)
  par(xpd=F)
})

#' @method length FuzzyClusterGroup
#' @name FuzzyClusterGroup-methods
setMethod('length','FuzzyClusterGroup',function(x) length(x@clusters))


setMethod('[','FuzzyClusterGroup',function(x,i) {
  x@clusters <- x@clusters[i]
  x
})
