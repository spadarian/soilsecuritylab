#' Prediction interval estimation
#'
#' This function estimates prediction intervals for a FuzzyCluster object. It shouldn't be neccesary to use it by itself because it is internally used by \code{\link{fuzzyRun}}.
#'
#' @param clust_data a \code{FuzzyCluster} object.
#' @param obs numerical vector (optional) corresponding to the observed values of a modeled variable. It is used in conjunction with the \code{pred} argument to estimimate the prediction interval for each \code{classes}.
#' @param pred numerical vector (optional) corresponding to the predicted values of a modeled variable. It is used in conjunction with the \code{obs} argument to estimimate the prediction interval for each \code{classes}.
#' @param .conf confidence level used when estimimating the prediction interval. Default value: 0.95.
#' @return \code{list}
#' @export

prediction_interval <- function(clust_data,obs=NULL,pred=NULL,.conf=0.95) {
  
  if(is.null(obs) | is.null(pred)) stop('Nothing to calculate. Provide all the inputs.')
  res <- obs-pred
  quant <- data.frame(res,class=attr(clust_data@U,'hard_clust'))
  alpha <- 1 - .conf
  quant <- tapply(quant$res,quant$class,function(x) quantile(x,c(alpha/2,1-(alpha/2)),na.rm=T))
  quant <- do.call(rbind,quant)
  quant <- data.frame(class=rownames(quant),quant)
  
  PICl <- as.matrix(quant[,2])
  PICu <- as.matrix(quant[,3])
  
  PICl[nrow(quant)] <- PICl[nrow(quant)]*2
  PICu[nrow(quant)] <- PICu[nrow(quant)]*2
  
  clust_data@U <- clust_data@U[,levels(quant$class)]
  
  PIl <- clust_data@U%*%PICl
  PIu <- clust_data@U%*%PICu
  
  PLl <- c(pred) + c(PIl)
  PLu <- c(pred) + c(PIu)
  
  PICP_matrix <- cbind(PLl,obs,PLu)
  
  PICP_count <- table(PICP_matrix[,1] <= PICP_matrix[,2] & PICP_matrix[,2] <= PICP_matrix[,3])
  
  PICP<- unname(100*PICP_count['TRUE']/sum(PICP_count))
  MPI <- sum(PLu-PLl,na.rm=T)/length(na.exclude(PLu))
  
  list(observed=obs,predicted=pred,PICl=PICl, PICu=PICu, PIl=PIl, PIu=PIu, PLl=PLl, PLu=PLu, PICP=PICP, MPI=MPI)
}