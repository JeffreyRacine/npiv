## Function for determining optimal spline complexity for the
## nonparametric IV estimation of Chen, Christensen and Kankanala
## (2021). Notation-wise, I try to follow their notation closely, and
## append .x and .w where needed for clarity (described further
## below).

npivJ <- function(Y,
                  X,
                  W,
                  X.eval=NULL,
                  K.w.degree=3,
                  K.w.segments=1,
                  J.x.degree=3,
                  J.x.segments=1,
                  J.x.segments.set=c(2,4,8,16,32)-1,
                  knots=c("uniform","quantiles"),
                  basis=c("tensor","additive","glp"),
                  check.is.fullrank=FALSE,
                  chol.pivot=FALSE,
                  lambda=sqrt(.Machine$double.eps)) {

    ## Match variable arguments to ensure they are valid

    basis <- match.arg(basis)
    knots <- match.arg(knots)

    ## Conduct some basic error checking to test for valid input

    if(missing(Y)) stop(" must provide Y")
    if(missing(X)) stop(" must provide X")
    if(missing(W)) stop(" must provide W")
    if(K.w.degree < 0) stop("K.w.degree must be a non-negative integer")
    if(J.x.degree < 0) stop("J.x.degree must be a non-negative integer")
    if(K.w.segments <= 0) stop("K.w.segments must be a positive integer")
    if(J.x.segments <= 0) stop("J.x.segments must be a positive integer")

    ## If specified, check that passed objects are of full rank

    if(check.is.fullrank) {
        if(!is.fullrank(Y)) stop("Y is not of full column rank")
        if(!is.fullrank(X)) stop("X is not of full column rank")
        if(!is.fullrank(W)) stop("W is not of full column rank")
    }

    ## See function npiv() in R code npiv.R for descriptions...

    ## Generate set of J1 J2 combinations given input
    ## J.x.segments.set. expand.grid() creates all possible
    ## combinations, then we select those with J2 > J1 (remove
    ## symmetric computations and those with J2=J1)

    J1.J2 <- expand.grid(J.x.segments.set,J.x.segments.set)
    J1.J2 <- subset(J1.J2,Var2 > Var1)
    
    J.x.J1.segments <- J1.J2[1,1]
    J.x.J2.segments <- J1.J2[1,2]

    ## Tim precomputes and reuses basis functions... these are
    ## computationally efficient so might not save much time, morover
    ## if parallelized could be a waste and even add overhead

    ## This line needs to be reworked/removed once I hear back from
    ## Tim about setting K for B during the search for J for Psi

    if(K.w.degree+K.w.segments < J.x.degree+J.x.segments) stop("K.w.degree+K.w.segments must be >= J.x.degree+J.x.segments")

    ## Generate basis functions for W - held constant it appears (XXX
    ## should use larger K.segments as J.segments increases??? Tim?)

    if(K.w.degree==0) {
        B.w <- matrix(1,NROW(W),1)
    } else {
        B.w <- prod.spline(x=W,
                           K=cbind(rep(K.w.degree,NCOL(W)),rep(K.w.segments,NCOL(W))),
                           knots=knots,
                           basis=basis)

        if(basis!="tensor") B.w <- cbind(1,B.w)
    }

    ## Generate basis functions for X for J1 and J2

    if(J.x.degree==0) {
        Psi.x.J1.eval <- Psi.x.J1 <- matrix(1,NROW(X),1)
        Psi.x.J2.eval <- Psi.x.J2 <- matrix(1,NROW(X),1)
    } else {
        Psi.x.J1.eval <- Psi.x.J1 <- prod.spline(x=X,
                            K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.J1.segments,NCOL(X))),
                            knots=knots,
                            basis=basis)
        Psi.x.J2.eval <- Psi.x.J2 <- prod.spline(x=X,
                            K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.J2.segments,NCOL(X))),
                            knots=knots,
                            basis=basis)

        if(!is.null(X.eval)) {
            Psi.x.J1.eval <- prod.spline(x=X,
                                      xeval=X.eval,
                                      K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.J1.segments,NCOL(X))),
                                      knots=knots,
                                      basis=basis)
            Psi.x.J2.eval <- prod.spline(x=X,
                                      xeval=X.eval,
                                      K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.J2.segments,NCOL(X))),
                                      knots=knots,
                                      basis=basis)

        }
        
        if(basis!="tensor") {
          Psi.x.J1 <- cbind(1,Psi.x.J1)
          Psi.x.J2 <- cbind(1,Psi.x.J2)
          Psi.x.J1.eval <- cbind(1,Psi.x.J1.eval)
          Psi.x.J2.eval <- cbind(1,Psi.x.J2.eval)
        }
      
    }

    ## Code re-use/storage where possible... generate all objects
    ## needed to compute the t-stat vector

    B.w.TB.w.inv <- chol2inv(chol(t(B.w)%*%B.w,pivot=chol.pivot))

    Psi.xJ1TB.wB.wTB.w.invB.w <- t(Psi.x.J1)%*%B.w%*%B.w.TB.w.inv%*%t(B.w)
    tmp.J1 <- chol2inv(chol(Psi.xJ1TB.wB.wTB.w.invB.w%*%Psi.x.J1+diag(lambda,NCOL(Psi.x.J1)),pivot=chol.pivot))%*%Psi.xJ1TB.wB.wTB.w.invB.w
    beta.J1 <- tmp.J1%*%Y
    
    U.J1 <- Y-Psi.x.J1%*%beta.J1
    err.J1 <- Psi.x.J1.eval%*%tmp.J1%*%U.J1

    Psi.xJ2TB.wB.wTB.w.invB.w <- t(Psi.x.J2)%*%B.w%*%B.w.TB.w.inv%*%t(B.w)
    tmp.J2 <- chol2inv(chol(Psi.xJ2TB.wB.wTB.w.invB.w%*%Psi.x.J2+diag(lambda,NCOL(Psi.x.J2)),pivot=chol.pivot))%*%Psi.xJ2TB.wB.wTB.w.invB.w
    beta.J2 <- tmp.J2%*%Y

    U.J2 <- Y-Psi.x.J2%*%beta.J2
    err.J2 <- Psi.x.J2.eval%*%tmp.J2%*%U.J2

    ## Compute asymptotic variances and covariances for the IV
    ## functions - these will be memory intensive for large n as they
    ## require taking the diagonal of an n x n matrix (could be a way
    ## around but simply note that as it stands the computation of the
    ## standard errors by brute force needs to be addressed)

    CJ1 <- t(Psi.x.J1)%*%B.w
    B.wUJ1 <- B.w*as.numeric(U.J1)
    rho <- CJ1%*%B.w.TB.w.inv%*%t(B.wUJ1)%*%(B.wUJ1)%*%B.w.TB.w.inv%*%t(CJ1)
    D.J1.inv <- chol2inv(chol(CJ1%*%B.w.TB.w.inv%*%t(CJ1),pivot=chol.pivot))
    D.J1.inv.rho.D.J1.inv <- D.J1.inv%*%rho%*%D.J1.inv

    asy.var.J1 <- diag(Psi.x.J1.eval%*%D.J1.inv.rho.D.J1.inv%*%t(Psi.x.J1.eval))

    CJ2 <- t(Psi.x.J2)%*%B.w
    B.wUJ2 <- B.w*as.numeric(U.J2)
    rho <- CJ2%*%B.w.TB.w.inv%*%t(B.wUJ2)%*%(B.wUJ2)%*%B.w.TB.w.inv%*%t(CJ2)
    D.J2.inv <- chol2inv(chol(CJ2%*%B.w.TB.w.inv%*%t(CJ2),pivot=chol.pivot))
    D.J2.inv.rho.D.J2.inv <- D.J2.inv%*%rho%*%D.J2.inv
    
    asy.var.J2 <- diag(Psi.x.J2.eval%*%D.J2.inv.rho.D.J2.inv%*%t(Psi.x.J2.eval))

    ## Compute the covariance

    asy.cov.J1.J2 <- diag(Psi.x.J1.eval%*%D.J1.inv%*%CJ1%*%B.w.TB.w.inv%*%t(B.wUJ1)%*%(B.wUJ2)%*%B.w.TB.w.inv%*%t(CJ2)%*%t(D.J2.inv)%*%t(Psi.x.J2.eval))

    ## Compute the denominator of the t-stat (the numerator is the
    ## difference between err.J1 and err.J2)

    asy.se <- sqrt(asy.var.J1+asy.var.J2-2*asy.cov.J1.J2)

    ## The t-stat vector - we take the sup (max) of this to determine
    ## the optimal value of J (segments/knots of the Psi.x basis)

    Z <- (err.J1-err.J2)/asy.se

    ## Return a list with various objects that might be of interest to
    ## the user

    return(list(Z=Z))

}
