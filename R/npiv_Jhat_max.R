## Function for determining upper limit of grid of J values for the
## nonparametric IV estimation of Chen, Christensen and Kankanala
## (2021). Notation-wise, I try to follow their notation closely, and
## append .x and .w where needed for clarity (described further
## below).

npiv_Jhat_max <- function(X,
                          W,
                          J.x.degree=3,
                          K.w.degree=3,
                          K.w.smooth=1,
                          knots=c("uniform","quantiles"),
                          basis=c("tensor","additive","glp"),
                          check.is.fullrank=FALSE,
                          chol.pivot=FALSE) {
  
  ## Match variable arguments to ensure they are valid
  
  basis <- match.arg(basis)
  knots <- match.arg(knots)
  
  ## Conduct some basic error checking to test for valid input
  
  if(missing(X)) stop(" must provide X")
  if(missing(W)) stop(" must provide W")
  if(K.w.degree < 0) stop("K.w.degree must be a non-negative integer")
  if(J.x.degree < 0) stop("J.x.degree must be a non-negative integer")
  
  ## If specified, check that passed objects are of full rank
  
  if(check.is.fullrank) {
    if(!is.fullrank(X)) stop("X is not of full column rank")
    if(!is.fullrank(W)) stop("W is not of full column rank")
  }
  
  ## Generate set of J K combinations given input
  ## J.x.degree and K.w.degree to search over up to a maximum resolution
  ## level of log(n, base = (2 * dim(X))), which will have a singular
  ## denominator matrix because 
  
  L.max <- max(floor(log(NROW(X), base = 2 * NCOL(X))), 3)
  J.x.segments.set <- (2^(1:L.max))-1
  J.w.segments.set <- (2^(1:L.max+K.w.smooth))-1

  ## In what follows we loop over _rows_  (makes for easy
  ##  parallelization if needed)
  
  test.val <- array(NA, dim = L.max)
  
  ## Temporary indication of where we are in the process
  
  print("Determining grid")
  
  for(ii in 1:L.max) {
    
    ## Temporary indication of where we are in the process
    
    print(paste("Row ",ii," of ",L.max))
    
    if((ii <= 2) || ((ii > 2) & (test.val[ii-2] <= 10*sqrt(NROW(X))))){
      
      J.x.segments <- J.x.segments.set[ii]
      J.w.segments <- J.w.segments.set[ii]
      
      ## Generate basis functions for W for J
      
      if(K.w.degree < J.x.degree) stop("K.w.degree must be >= J.x.degree")
      
      ## Estimate measure of ill-posedness if IV model; if regression set to 1 by default
      
      if(all(X == W)){
        
        s.hat.J <- 1
        
      } else {
        
        ## Generate basis functions for W for J
        
        if(K.w.degree==0) {
          B.w.J <- matrix(1,NROW(W),1)
        } else {
          B.w.J <- prod.spline(x=W,
                               K=cbind(rep(K.w.degree,NCOL(W)),rep(J.w.segments,NCOL(W))),
                               knots=knots,
                               basis=basis)
          
          if(basis!="tensor") {
            B.w.J <- cbind(1,B.w.J)
          }
        }
        
        ## Generate basis functions for X for J
        
        if(J.x.degree==0) {
          Psi.x.J <- matrix(1,NROW(X),1)
        } else {
          Psi.x.J <- prod.spline(x=X,
                                 K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                 knots=knots,
                                 basis=basis)
          
          if(basis!="tensor") {
            Psi.x.J <- cbind(1,Psi.x.J)
          }
        }
        
        ## Compute \hat{s}_J
        
        G.w.J.inv <- chol2inv(chol(t(B.w.J)%*%B.w.J,pivot=chol.pivot))
        G.x.J.inv <- chol2inv(chol(t(Psi.x.J)%*%Psi.x.J,pivot=chol.pivot))
        S.J <- t(Psi.x.J)%*%B.w.J
        tmp <- sqrtm(G.x.J.inv)%*%S.J%*%sqrtm(G.w.J.inv)
        
        s.hat.J <- min(svd(tmp)$d)
        
      }
      
      ## Compute test value
      
      J.x.dim <- (J.x.segments + J.x.degree)^NCOL(X)
      
      test.val[ii] <- J.x.dim*sqrt(log(J.x.dim))*max((0.1*log(NROW(X)))^4,1/s.hat.J)
      
    } else {
      
      test.val[ii] <- test.val[ii - 1]
      
    }
    
  }
  
  ## Find appropriate value
  
  L.hat.max <- which((test.val[-length(test.val)] <= 10*sqrt(NROW(X))) & (10*sqrt(NROW(X)) < test.val[-1]))[1]
  
  ## Return largest feasible grid in case empty
  
  if(is.na(L.hat.max)){
    L.hat.max <- L.max
  }
  
  ## Return values now for use in npiv_J
  
  J.x.segments.set <- J.x.segments.set[1:L.hat.max]
  J.w.segments.set <- J.w.segments.set[1:L.hat.max]
  
  J.hat.max <- (max(J.x.segments.set) + J.x.degree)^NCOL(X)
  
  alpha.hat <- min(0.5, 1/J.hat.max)
  
  ## Return a list with various objects that might be of interest to
  ## the user
  
  return(list(J.x.segments.set=J.x.segments.set,
              J.w.segments.set=J.w.segments.set,
              J.hat.max=J.hat.max,
              alpha.hat=alpha.hat))
  
}
