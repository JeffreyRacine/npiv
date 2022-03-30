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

## S3 related functions

explodePipe <- function(formula){
  tf <- as.character(formula)  
  tf <- tf[length(tf)]

  eval(parse(text=paste("c(",
               ifelse(length(as.character(formula)) == 3,
                      'strsplit(as.character(formula)[2]," *[+] *"),',""),
               'strsplit(strsplit(tf," *[|] *")[[1]]," *[+] *"))')))
}
