## Nonparametric IV estimation per Chen, Christensen and Kankanala
## (2021), notation-wise, I try to follow their notation closely, and
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
                 J.x.segments.set=c(2,4,8,16,32,64,128),
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

    ## Tim precomputes and reuses basis functions... thet are
    ## computationally efficient so might not save much time, morover
    ## if parallelized could be a waste and even add overhead

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
        Psi.x.J1 <- prod.spline(x=X,
                            K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.J1.segments,NCOL(X))),
                            knots=knots,
                            basis=basis)
        Psi.x.J2 <- prod.spline(x=X,
                            K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.J2.segments,NCOL(X))),
                            knots=knots,
                            basis=basis)

        if(basis!="tensor") Psi.x.J1 <- cbind(1,Psi.x.J1)
        if(basis!="tensor") Psi.x.J2 <- cbind(1,Psi.x.J2)
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

           if(basis!="tensor") Psi.x.J1.eval <- cbind(1,Psi.x.J1.eval)
           if(basis!="tensor") Psi.x.J2.eval <- cbind(1,Psi.x.J2.eval)
        }
    }

    BTB.inv <- chol2inv(chol(t(B.w)%*%B.w,pivot=chol.pivot))

    Psi.xJ1TB.wB.wTB.w.invB.w <- t(Psi.x.J1)%*%B.w%*%BTB.inv%*%t(B.w)
    beta.J1 <- chol2inv(chol(Psi.xJ1TB.wB.wTB.w.invB.w%*%Psi.x.J1+diag(lambda,NCOL(Psi.x.J1)),pivot=chol.pivot))%*%Psi.xJ1TB.wB.wTB.w.invB.w%*%Y
    U.J1 <- Y-Psi.x.J1%*%beta.J1

    Psi.xJ2TB.wB.wTB.w.invB.w <- t(Psi.x.J2)%*%B.w%*%BTB.inv%*%t(B.w)
    beta.J2 <- chol2inv(chol(Psi.xJ2TB.wB.wTB.w.invB.w%*%Psi.x.J2+diag(lambda,NCOL(Psi.x.J2)),pivot=chol.pivot))%*%Psi.xJ2TB.wB.wTB.w.invB.w%*%Y
    U.J2 <- Y-Psi.x.J2%*%beta.J2

    ## Compute asymptotic variances and covariances for the IV
    ## functions - these will be memory intensive for large n as they
    ## require taking the diagonal of an n x n matrix (could be a way
    ## around but simply note that as it stands the computation of the
    ## standard errors by brute force needs to be addressed)

    CJ1 <- (t(Psi.x.J1)%*%B.w)
    B.wUJ1 <- B.w*as.numeric(U.J1)
    rho <- CJ1%*%BTB.inv%*%t(B.wUJ1)%*%(B.wUJ1)%*%BTB.inv%*%t(CJ1)
    D.J1.inv <- chol2inv(chol(CJ1%*%BTB.inv%*%t(CJ1),pivot=chol.pivot))
    D.J1.inv.rho.D.J1.inv <- D.J1.inv%*%rho%*%D.J1.inv
    asy.var.J1 <- diag(Psi.x.J1.eval%*%D.J1.inv.rho.D.J1.inv%*%t(Psi.x.J1.eval))

    CJ2 <- (t(Psi.x.J2)%*%B.w)
    B.wUJ2 <- B.w*as.numeric(U.J2)
    rho <- CJ2%*%BTB.inv%*%t(B.wUJ2)%*%(B.wUJ2)%*%BTB.inv%*%t(CJ2)
    D.J2.inv <- chol2inv(chol(CJ2%*%BTB.inv%*%t(CJ2),pivot=chol.pivot))
    D.J2.inv.rho.D.J2.inv <- D.J2.inv%*%rho%*%D.J2.inv
    asy.var.J2 <- diag(Psi.x.J2.eval%*%D.J2.inv.rho.D.J2.inv%*%t(Psi.x.J2.eval))

    ## Compute the covariance

    cov.J1.J2 <- diag(Psi.x.J1.eval%*%D.J1.inv%*%CJ1%*%BTB.inv%*%t(B.wUJ1)%*%(B.wUJ2)%*%BTB.inv%*%t(CJ2)%*%D.J2.inv%*%Psi.x.J2.eval)

    asy.se <- sqrt(asy.var.J1+asy.var.J2-2*cov.J1.J2)

    ## Return a list with various objects that might be of interest to
    ## the user

    return(list(h=h,
                h.asy.se=asy.se,
                K.w.degree=K.w.degree,
                K.w.segments=K.w.segments,
                J.x.degree=J.x.degree,
                J.x.segments=J.x.segments,
                beta=beta,
                B.w=B.w,
                Psi.x=Psi.x,
                residuals.sample=Y-Psi.x%*%beta,
                residuals.eval=Y-Psi.x.eval%*%beta)

}

