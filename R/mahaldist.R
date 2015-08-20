.mahaldist <- function (x,y=NULL,w=NULL) {
  if (is.null(w)) w <- solve(cov(x))
  if (is.null(y)) {
    apply(x,1,function(point){
      tmp_ <- matrix(rep(point,dim(x)[1]),dim(x)[1],byrow=T)
      Dc <- tmp_ - x
      Re(rowSums((Dc%*%w)*Conj(Dc)))
    })
  } else {
    apply(y,1,function(centroid){
      tmp_ <- matrix(rep(centroid,dim(x)[1]),dim(x)[1],byrow=T)
      Dc = tmp_ - x
      Re(rowSums((Dc%*%w)*Conj(Dc)))
    })
  }
}