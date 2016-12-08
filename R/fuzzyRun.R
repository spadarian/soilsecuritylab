#' Fuzzy k-means with extragrade sequence
#'
#' This function generates an FuzzyClusterGroup for a sequence of cluster numbers.
#'
#' @param data a \code{data.frame} with numerical data. Each column corresponds to a dimension of a \code{n}-dimensional clustering domain.
#' @param phi a numerical value representing the degree of fuzziness or overlap of the generated clusters.
#' @param classes numerical (integer) vector with the number of cluster to generate.
#' @param disttype numerical value representing the type of distance to calculate for the clustering. Possible values are: 1 (Euclidean), 2 (Diagonal), and 3 (default, Mahalanobis).
#' @param exp_eg the expected fraction of extragrades. If not provided, it is assumed to be dependant on the number of clusters (\code{1/(1+nclusters)}).
#' @param obs numerical vector (optional) corresponding to the observed values of a modeled variable. It is used in conjunction with the \code{pred} argument to estimimate the \code{\link{prediction_interval}} for each \code{classes}.
#' @param pred numerical vector (optional) corresponding to the predicted values of a modeled variable. It is used in conjunction with the \code{obs} argument to estimimate the \code{\link{prediction_interval}} for each \code{classes}.
#' @param .conf confidence level used when estimimating the \code{\link{prediction_interval}}. Default value: 0.95.
#' @param alpha_values numerical vector (optional), corresponding to optimised alpha values. This function estimates optimum alpha values, so this argumet is mostly used at developing.
#' @return \code{FuzzyClusterGroup}
#' @references De Gruijter, J. and McBratney, A. 1988. A modified fuzzy k-means method for predictive classification. In: 1st Conference of the International Federation of Classification Societies. Ed. by H. Bock. Elsevier Science, Amsterdam: pp. 97-104.
#' @examples
#' fkm <- fuzzyRun(cars,1.5,2:6,3)
#' @importFrom DEoptim DEoptim
#' @export

fuzzyRun <- function(data,phi,classes,disttype=3,exp_eg=NULL,obs=NULL,pred=NULL,.conf=0.95,alpha_values=NULL){

  if(!(disttype %in% 1:3 )) stop('disttype should have a value of 1, 2, or 3.')
  dist_name <- c('Euclidean','Diagonal','Mahalanobis')[disttype]
  if(any(classes < 1 )) stop('Numbers of clusters should be greater than 1.')

  lincon <- function (x,y) {
    n <- length(x)
    mx <- mean(x)
    my <- mean(y)
    sx <- sum((x-mx)^2)/n
    sy <- sum((y-my)^2)/n
    sxy <- sum((x-mx)*(y-my))/n
    2*sxy/(sx+sy+(mx-my)^2)
  }

  classes <- as.integer(classes)

  ini <- proc.time()[3]
  if (is.null(alpha_values)) {
    # require(DEoptim)
    message('Initialising alpha optimisation:')
    assign('min_alpha',0.0000001,envir=soilsecuritylab.env)
    result <- lapply(classes,function(x){
      message(paste0('Optimising alpha for ',x,' clusters'))
      if (is.null(exp_eg)) exp_eg <- 1/(1+x)
      #### VTR=0.0005; NP=20; itermax=20  using optimised blas
      fkm__ <- function(...){
        tmp <- .fuzzy_extragrades_C(...)
        abs(tmp[[1]])
      }
      tmp_ <- try(DEoptim::DEoptim(fkm__,get('min_alpha',envir=soilsecuritylab.env),1,control=list(parallelType=1,VTR=0.0005,trace=F,itermax=20,NP=20,parVar=c('.mahaldist','.fuzzy_extragrades_C'),packages=c('Rcpp','soilsecuritylab')),data=as.matrix(data),phi=phi,nclass=x,disttype=disttype,maxiter=300,toldif=0.001,exp_eg=exp_eg,optim=T), silent=TRUE)
      if (class(tmp_) == 'try-error'){
        tmp_ <- DEoptim::DEoptim(fkm__,get('min_alpha',envir=soilsecuritylab.env),1,control=list(parallelType=0,VTR=0.0005,trace=F,itermax=20,NP=20,parVar=c('.mahaldist','.fuzzy_extragrades_C'),packages=c('Rcpp','soilsecuritylab')),data=as.matrix(data),phi=phi,nclass=x,disttype=disttype,maxiter=300,toldif=0.001,exp_eg=exp_eg,optim=T)
      }
      assign('min_alpha',tmp_$optim$bestmem,envir=soilsecuritylab.env)
      print(paste('*** Optim alpha:',tmp_$optim$bestmem))
      tmp_$nclass <- x
      tmp_
    })
    rm(min_alpha,envir=soilsecuritylab.env)
    message('Generating clusters:')
    result_ <- lapply(result,function(x){
      message(paste0(x$nclass,' clusters using alpha: ',x$optim$bestmem))
      if (is.null(exp_eg)) exp_eg <- 1/(1+x$nclass)
      cluster <- .fuzzy_extragrades_C(x$optim$bestmem,as.matrix(data),phi,x$nclass,disttype,300,0.001,exp_eg=exp_eg,optim=F)
      cluster <- new('FuzzyCluster',data=as.matrix(data),U=cluster$U,W=cluster$W,centroids=cluster$centroids,phi=phi,classes=as.integer(x$nclass),distance=dist_name,alpha=unname(x$optim$bestmem),`Ue_mean-Ue_req`=abs(cluster$`Ue_mean-Ue_req`),iterations=cluster$iterations)
      if (!is.null(obs) & !is.null(pred)) {
        pred_int <- prediction_interval(cluster,obs,pred,.conf=.conf)
        cluster@pred_int <- pred_int
      } else .conf <- NULL
      cluster
    })
  } else {
    if (length(alpha_values) != length(classes)) stop('An alpha value should be supplied for every class') else {
      result <- matrix(c(classes,alpha_values),ncol=2)
    }
    message('Generating clusters:')
    result_ <- lapply(result,1,function(x){
      message(paste0(x[1],' clusters using alpha: ',x[2]))
      cluster <- .fuzzy_extragrades_C(x[2],as.matrix(data),phi,as.integer(x[1]),disttype,300,0.001,exp_eg=exp_eg,optim=F)
      cluster <- new('FuzzyCluster',data=data,U=cluster$U,W=cluster$W,centroids=cluster$centroids,phi=phi,classes=as.integer(x$nclass),distance=dist_name,alpha=unname(x$optim$bestmem),`Ue_mean-Ue_req`=abs(cluster$`Ue_mean-Ue_req`),iterations=cluster$iterations)
      if (!is.null(obs) & !is.null(pred)) {
        pred_int <- prediction_interval(cluster,obs,pred,.conf=.conf)
        cluster@pred_int <- pred_int
      } else .conf <- NULL
      cluster
    })
  }
  end <- proc.time()[3]
  message(paste0('Elapsed time: ',round((end-ini)/60,2),' minutes'))

  names(result_) <- as.character(classes)

  new('FuzzyClusterGroup',clusters=result_,confidence=.conf)
}
