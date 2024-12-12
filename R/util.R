## This function is based on functions in the limma package and
## corpcor package (is.positive.definite)... check the condition
## number of a matrix based on the ratio of max/min eigenvalue.  Note
## that the definition
## tol=max(dim(x))*max(sqrt(abs(e)))*.Machine$double.eps is exactly
## compatible with the conventions used in "Octave" or "Matlab".  Note
## that for weighted regression you simply use x*L which conducts
## row-wise multiplication (i.e. diag(L)%*%X not necessary). Note also
## that crossprod(X) is significantly faster than t(X)%*%X (matrix is
## symmetric so only use lower triangle).

is.fullrank <- function(x)
{
  e <- eigen(crossprod(as.matrix(x)), symmetric = TRUE, only.values = TRUE)$values
  e[1] > 0 && abs(e[length(e)]/e[1]) > max(dim(as.matrix(x)))*max(sqrt(abs(e)))*.Machine$double.eps
}


## statistical functions

RSQfunc <- function(y,y.pred,weights=NULL) {
  if(!is.null(weights)) {
    y <- y*sqrt(weights)
    y.pred <- y.pred*sqrt(weights)
  }
  y.mean <- mean(y)
  return((sum((y-y.mean)*(y.pred-y.mean))^2)/(sum((y-y.mean)^2)*sum((y.pred-y.mean)^2)))
}


## numerically sound matrix positive definite square root

sqrtm2 <- function(x) {
  
  x.eig <- eigen(x, symmetric = TRUE)
  lambda <- apply(as.matrix(x.eig$values), 1, function(x) sqrt(max(x, 0)))
  return(x.eig$vectors%*%diag(lambda)%*%t(x.eig$vectors))
  
}

## avoid division by zero

NZD <- function(a) {
  
  ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
  
}

## Function that determines the dimension of the multivariate basis
## without precomputing it.

dimbs <- function(basis="additive",degree=NULL,segments=NULL) {
  
  ## This function computes the dimension of the glp basis without the
  ## memory overhead associated with computing the glp basis itself
  ## (thanks to Zhenghua Nie)
  
  two.dimen<- function(d1,d2,nd1,pd12){
    if(d2 ==1) {
      ret <- list()
      ret$d12 <- pd12
      ret$nd1 <- nd1
      return(ret)
    }
    d12 <- d2
    if(d1-d2>0){
      for(i in 1:(d1-d2)){
        d12 <- d12+d2*nd1[i]
      }}
    if(d2>1){
      for(i in 2:d2){
        d12 <- d12 + (i*nd1[d1-i+1])
      }
    }
    d12 <- d12 + nd1[d1]   ## The maximum number
    
    nd2 <- nd1  ## Calculate nd2
    if(d1>1){
      for(j in 1:(d1-1)) {
        nd2[j] <- 0
        for(i in j:max(0,j-d2+1)) {
          if(i > 0) {
            nd2[j] <- nd2[j] + nd1[i]                  
          }
          else {
            nd2[j] <- nd2[j] + 1  ## nd1[0] always 1
          }
        }
      }
    }
    if(d2>1) {
      nd2[d1] <- nd1[d1]
      for(i in (d1-d2+1):(d1-1)) nd2[d1] <- nd2[d1]+nd1[i]
    }
    else {
      nd2[d1] <- nd1[d1]
    }
    ret <- list()
    ret$d12 <- d12
    ret$nd1 <- nd2 
    
    return(ret)
  }
  
  ## Some basic error checking
  
  if(basis!="additive" & basis!="glp" & basis!="tensor") stop(" Error: basis must be either additive, glp, or tensor")
  
  K <- cbind(degree,segments)
  
  ncol.bs <- 0

  if(basis=="additive") {
    if(any(K[,1] > 0))
      ncol.bs <- sum(rowSums(K[K[,1]!=0,,drop=FALSE])-1)
  }
  if(basis=="glp") {
    dimen <- rowSums(K[K[,1]!=0,,drop=FALSE])-1
    dimen <- dimen[dimen>0] ## Delete elements which are equal to 0.
    dimen <- sort(dimen,decreasing=TRUE) ## Sort the array to save memory when doing the computation.
    k <-length(dimen)
    if(k==0) {
      ncol.bs <- 0
    } else {
      nd1 <- rep(1,dimen[1])   ## At the beginning,  we have one for [1, 2, 3, ..., dimen[1]]
      nd1[dimen[1]] <- 0       ## nd1 represents the frequency for every element of [1, 2, 3, ..., dimen[1]]
      ncol.bs <- dimen[1]
      if(k>1) {
        for(i in 2:k) {
          dim.rt <- two.dimen(dimen[1],dimen[i],nd1,ncol.bs)
          nd1 <- dim.rt$nd1
          ncol.bs <- dim.rt$d12
        }
        ncol.bs <- dim.rt$d12+k-1
      }
    }
  }
  if(basis=="tensor") {
    if(any(K[,1] > 0))
      ncol.bs <- prod(rowSums(K[K[,1]!=0,,drop=FALSE]))
  }

  return(ncol.bs)
  
}

