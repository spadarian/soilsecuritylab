.fuzzy_extragrades <- function(alpha, data, phi, nclass, disttype, maxiter = 300, toldif = 0.001,exp_eg=1/(1+nclass),.optim=F) {
  data <- as.matrix(data)
  Ue_req <- exp_eg
  ndata <- dim(data)[1]
  ndim <- dim(data)[2]
  centroids <- matrix(0,nclass,ndim)
  dist <- matrix(0,ndata,nclass)

  U <- matrix(rnorm(ndata*nclass,0.5,0.01),ndata,nclass)
  U <- U/rowSums(U)

  if (disttype == 1){
    W <- diag(ndim)
    dist_name <- 'Euclidean'
  }
  if (disttype == 2){
    W <- diag(ndim)*cov(data)
    dist_name <- 'Diagonal'
  }
  if (disttype == 3){
    W <- solve(cov(data))
    dist_name <- 'Mahalanobis'
  }

  obj <- 0

  Ue <- round(rep(1,ndata),15)-round(rowSums(U),15)
  uphi <- U^phi
  uephi <- Ue^phi
  a1 <- (1-alpha)/alpha

  ###### Initialice ######

  c1 <- t(uphi)%*%data
  t1 <- colSums(uphi)
  t1 <- matrix(t1,(length(t1)),ndim)
  centroids <- c1/t1

  ##### Calculate distances to centroids #####

  # library(plyr)
  if (disttype == 1) { # Euclidian

    dist <- apply(centroids,1,function(x){
      sqrt(rowSums((data - matrix(rep(x,ndata),ndata,byrow=T))^2))
    })

  } else { # Diagonal and Mahalanobis
    dist <- sqrt(.mahaldist(data,centroids,W))
  }

  for (i in 1:maxiter) {
    # Calculate centroids
    ufi <- uphi-(matrix(a1,ndata,nclass)*dist^(-4)*matrix(uephi,length(uephi),nclass))
    c1 <- t(ufi)%*%data
    t1 <- colSums(ufi)
    t1 <- matrix(t1,(length(t1)),ndim)
    centroids <- c1/t1

    # Calculate distances to centroids
    if (disttype == 1) { # Euclidean
      dist <- apply(centroids,1,function(x){
        sqrt(rowSums((data - matrix(rep(x,ndata),ndata,byrow=T))^2))
      })
    } else { # Diagonal and Mahalanobis
      dist <- sqrt(.mahaldist(data,centroids,W))
    }

    # Save previous iteration
    U_old <- U
    obj_old <- obj

    # Calculate new membership matrix
    tmp <- dist^(-2/(phi-1))
    tm2 <- dist^(-2)
    s2 <- (a1*rowSums(tm2))^(-1/(phi-1))

    t1 <- rowSums(tmp)
    t2 <- matrix(t1,(length(t1)),nclass)+matrix(s2,(length(s2)),nclass)
    U <- tmp/t2
    Ue <- round(rep(1,ndata),15)-round(rowSums(U),15)
    uphi <- U^phi
    uephi <- Ue^phi

    # Calculate obj function
    o1 <- (dist^2)*uphi
    d2 <- dist^(-2)
    o2 <- uephi*rowSums(d2)
    obj <- alpha*sum(o1)+(1-alpha)*sum(o2)

    # Check for convergence
    dif <- obj_old-obj
    difU <- sqrt((U-U_old)*(U-U_old))
    Udif <- sum(difU)
    if (dif < toldif & Udif < toldif) {break}
  }

  # Update centroids

  ufi <- uphi-(matrix(a1,ndata,nclass)*dist^(-4)*matrix(uephi,length(uephi),nclass))
  c1 <- t(ufi)%*%data
  t1 <- colSums(ufi)
  t1 <- matrix(t1,(length(t1)),ndim)
  centroids <- c1/t1
  U_end <- cbind(U,Ue)
  colnames(U_end) <- c(1:nclass,'Eg')
  attributes(U_end)$hard_clust <- colnames(U_end)[apply(U_end,1,which.max)]
  Ue_mean <- mean(Ue)

  #   dist_centroids <- mahaldist(centroids)

  CI <- plyr::aaply(U_end,1,function(x){
    tmp <- sort(x,decreasing=T)
    tmp[2]/tmp[1]
  })

  attributes(U_end)$CI <- CI

  result_ <- new('FuzzyCluster',data=data,U=U_end,W=W,centroids=centroids,phi=phi,classes=as.integer(nclass),distance=dist_name,alpha=unname(alpha),`Ue_mean-Ue_req`=abs(Ue_mean-Ue_req),iterations=i)

  if (.optim) return(abs(Ue_mean-Ue_req)) else return(result_)
}
