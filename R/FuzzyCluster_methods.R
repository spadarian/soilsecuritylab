#' Methods for FuzzyCluster class
#' 
#' Common methods implemented for the FuzzyCluster class.
#' 
#' @param x,object FuzzyCluster object
#' @param .plot numeric (optional). If just one plot is required it can be specified here.
#' @param ... arguments to be passed to subsequent methods, e.g., plot.FuzzyCluster for the plot method 
#' @method plot FuzzyCluster
#' @examples
#' ### Create FuzzyClusterGroup
#' fkm <- fuzzyRun(cars,1.5,2:4,3)
#' 
#' ### Functions to be used on part of the FuzzyClusterGroup
#' plot(fkm@@clusters$`2`,.plot=1)
#' 
#' summary(fkm@@clusters$`2`,latex=TRUE)
#' @name FuzzyCluster-methods

setMethod('plot','FuzzyCluster',function(x,.plot=NULL,...){
  if (1 %in% .plot | is.null(.plot)) {
    simb <- attributes(x@U)$hard_clust
    simb[simb == 'Eg'] <- 0
    simb <- as.numeric(simb)+1
    colours_ <- apply(x@U[,-ncol(x@U)],1,function(y){
      tmp_ <- c(which.max(y),max(y))
      rgb(colorRamp(c('black',tmp_[1]+1))(tmp_[2]),maxColorValue=255)
    })
    if (dim(x@data)[2] == 1) plot(x@data,col=colours_,pch=simb,main=paste0(x@classes,' classes'),ylab=colnames(x@data)) else plot(as.data.frame(x@data),col=colours_,pch=simb,main=paste0(x@classes,' classes'))
    if (interactive()) par(ask=T)
  }
  if (2 %in% .plot | is.null(.plot)) {
    if (dim(x@data)[2] == 1) {
      plot(x@data,col=rgb(colorRamp(c('grey','black'))(x@U[,'Eg']),maxColorValue=255),pch=16,main='Extragrade class\nmembership',ylab=colnames(x@data))
    } else {
      layout(matrix(1:2,ncol=2), width = c(20,1),height = c(1,1))
      plot(as.data.frame(x@data),col=rgb(colorRamp(c('grey','black'))(x@U[,'Eg']),maxColorValue=255),pch=16,main='Extragrade class\nmembership')
      
      legend_image <- as.raster(matrix(rgb(colorRamp(c('grey','black'))(seq(1,0,length.out = 10)),maxColorValue=255), ncol=1))
      mar_orig <- par()$mar
      par(mar=c(2,0,2,1.5))
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      par(xpd=NA)
      text(x=1, y = seq(0.2,0.8,l=3), labels = seq(0,1,l=3))
      rasterImage(legend_image, -2, 0.2, -1,0.8)
      par(mar=mar_orig)
      par(xpd=F)
      layout(matrix(1))
    }
    if (interactive()) par(ask=T)
  }
  if (3 %in% .plot | is.null(.plot)) {
    if (dim(x@data)[2] == 1) {
      plot(x@data,col=rgb(colorRamp(c('grey','black'))(attributes(x@U)$CI),maxColorValue=255),pch=16,main='Confusion Index',ylab=colnames(x@data))
    } else {
      layout(matrix(1:2,ncol=2), width = c(20,1),height = c(1,1))
      plot(as.data.frame(x@data),col=rgb(colorRamp(c('grey','black'))(attributes(x@U)$CI),maxColorValue=255),pch=16,main='Confusion Index')
      
      legend_image <- as.raster(matrix(rgb(colorRamp(c('grey','black'))(seq(1,0,length.out = 10)),maxColorValue=255), ncol=1))
      mar_orig <- par()$mar
      par(mar=c(2,0,2,1.5))
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      par(xpd=NA)
      text(x=1, y = seq(0.2,0.8,l=3), labels = seq(0,1,l=3))
      rasterImage(legend_image, -2, 0.2, -1,0.8)
      par(mar=mar_orig)
      par(xpd=F)
      layout(matrix(1))
    }
    if (interactive()) par(ask=T)
  }
  if (4 %in% .plot | is.null(.plot)) {
    if(length(x@pred_int) > 0){
      tmp_ <- data.frame(pred=x@pred_int$predicted,cluster=factor(attributes(x@U)$hard_clust))
      boxplot(pred~cluster,tmp_,xlab='Class',ylab='Predicted')
      if (interactive()) par(ask=T)
    }
  }
  if (5 %in% .plot | is.null(.plot)) {
    if(length(x@pred_int) > 0){
      tmp_ <- data.frame(res=x@pred_int$observed-x@pred_int$predicted,cluster=factor(attributes(x@U)$hard_clust))
      boxplot(res~cluster,tmp_,xlab='Class',ylab='Residuals')
      abline(h=0,lty=2)
      if (interactive()) par(ask=T)
    }
  }
  if (interactive()) par(ask=F)
})

#' @param extragrade_label string to use as label for the extragrade class. Default set to \code{'Eg' }.
#' @param decimals number of decimals to display.
#' @param latex if \code{FALSE} (default) the summary is displayed as a regular R output. If \code{TRUE}, the output is formated as a latex table (using booktabs). 
#' @param caption if \code{latex} is set to \code{TRUE}, is the sting used as a caption for the table.
#' @param label if \code{latex} is set to \code{TRUE}, is the sting used as a label for the table.
#' @method summary FuzzyCluster
#' @name FuzzyCluster-methods

setMethod('summary','FuzzyCluster',function(object,extragrade_label='Eg',decimals=2,latex=F,caption=NULL,label=NULL){
  res_mean <- tapply(object@pred_int$predicted-object@pred_int$observed,attributes(object@U)$hard_clust,mean)[c(seq(object@classes),'Eg')]
  table <- round(cbind(rbind(object@centroids,matrix(NA,nrow=1,ncol=ncol(object@centroids))),object@pred_int$PICl,res_mean,object@pred_int$PICu),2)
  rownames(table)[nrow(table)] <- extragrade_label
  colnames(table)[(ncol(table)-2):ncol(table)] <- c('PI_L','Mean','PI_U')
  table.dim <- dim(table)
  table_W <- object@W
  if (latex) {
    if (any(rownames(installed.packages()) %in% 'xtable')) {
      colnames(table)[(ncol(table)-2):ncol(table)] <- c('$PI^L$','Mean','$PI^U$')
      table <- xtable::xtable(table,align=c('l',rep('c',table.dim[2])))
      table_capt <- capture.output(cat(print(table, booktabs=TRUE, caption.placement='top', sanitize.rownames.function=function(x){x},
                                             sanitize.colnames.function=function(x){x},NA.string='--',print.results = FALSE)))
      index<- 5 + ifelse(is.null(label),0,1) + ifelse(is.null(caption),0,1)
      table_ <- c(table_capt[1:index],'\\toprule',paste0('Cluster & \\multicolumn{',table.dim[2]-3,'}{c}{Centroids} & \\multicolumn{3}{c}{Cluster residuals} \\\\'),
                  paste0('\\cmidrule(rl){2-',table.dim[2]-2,'}\\cmidrule(rl){',table.dim[2]-1,'-',table.dim[2]+1,'}'),
                  table_capt[(index+2):length(table_capt)])
      table_ <- gsub(paste0(c('-0.',rep('0',decimals)),collapse=''),'0.00',table_,perl=T)
      
      table_1 <- table_[5:(5+table.dim[1]+7)]
      table_1 <- c('\\subfloat[ Class centroids and prediction intervals (PI)]{',paste0('\\label{',label,'1}'),table_1,'}\\\\')
      table_2 <- capture.output(cat(print(xtable::xtable(table_W,align=c('l',rep('c',dim(table_W)[2]))),booktabs=TRUE,sanitize.rownames.function=function(x){x},
                                          sanitize.colnames.function=function(x){x},NA.string='--',print.results = FALSE)))[5:(5+dim(table_W)[1]+5)]
      table_2 <- c('\\subfloat[Variance-covariance matrix]{',paste0('\\label{',label,'2}'),table_2,'}')
      table_final <- c('\\begin{table}[ht]',paste0('\\caption{',caption,'}'),paste0('\\label{',label,'}'),'\\centering',table_1,table_2,'\\end{table}')
      cat(table_final,sep='\n')
    } else {
      warning('The package xtable is not installed')
      table
    }
  } else table
})

#' @method vcov FuzzyCluster
#' @name FuzzyCluster-methods

setMethod('vcov','FuzzyCluster',function(object){
  X <- scale(object@data, scale = F)
  CX <- (1/(nrow(object@data)-1))*(t(X)%*%X)
  CX
})

#' @param newdata a \code{data.frame} object in which to look for variables with which to predict.
#' @method predict FuzzyCluster
#' @name FuzzyCluster-methods

setMethod('predict','FuzzyCluster',function(object,newdata){
  
  obs_var <- colnames(object@data)
  if (class(try(newdata[,obs_var],T)) == 'try-error') {
    stop(paste0('New data does not include all the clustered variables. Check column names or include all necessary variables.\n  Required variables: ',paste0(obs_var,collapse=', ')))
  }
  
  newdata <- data.frame(newdata[,obs_var])
  
  newdata <- na.exclude(newdata)
  
  if (object@distance == 'Euclidean') {
    dist.m <- apply(object@centroids,1,function(x){
      sqrt(rowSums((as.matrix(newdata) - matrix(rep(x,dim(newdata)[1]),dim(newdata)[1],byrow=T))^2))
    })
  } else dist.m <- sqrt(.mahaldist(as.matrix(newdata),object@centroids,object@W))
  m1_ <- dist.m^(-2/(object@phi-1))
  m2_ <- dist.m^(-2)
  pond <- (1-object@alpha)/object@alpha
  s2 <- (pond*rowSums(m2_))^(-1/(object@phi-1))
  t1 <- rowSums(m1_)
  t2 <- matrix(t1,(length(t1)),dim(dist.m)[2])+matrix(s2,(length(s2)),dim(dist.m)[2])
  U_new <- round(m1_/t2,15)
  U_e_new <- rep(1,dim(U_new)[1])-rowSums(U_new)
  U_final <- cbind(U_new,Eg=U_e_new)
  tmp <- sapply(1:ncol(U_final), function(col){
    eval(parse(text=paste0('U_final[,',col,'] > ',paste0('U_final[,',(1:ncol(U_final))[-col],']'),collapse=' & ')))
  })
  class <- character(nrow(tmp))
  for(i in 1:ncol(tmp)){
    class[which(tmp[,i])] <- colnames(U_final)[i]
  }
  message(paste0('Number of extragrades detected: ',length(which(class == 'Eg'))))
  
  PICl <- c(U_final%*%object@pred_int$PICl)
  PICu <- c(U_final%*%object@pred_int$PICu)
  
  predicted <- data.frame(PICl=PICl,PICu=PICu)
  attributes(predicted)$hard_clust <- class
  attributes(predicted)$U <- U_final
  
  predicted
})