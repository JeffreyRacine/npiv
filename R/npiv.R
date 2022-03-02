## Nonparametric IV estimation per Chen, Christensen and Kankanala
## (2021), notation-wise, I try to follow their notation closely, and
## append .x and .w where needed for clarity (described further
## below).

npiv <- function(Y,
                 X,
                 W,
                 X.eval=NULL,
                 deriv.index=1,
                 deriv.order=1,
                 K.w.degree=3,
                 K.w.segments=1,
                 J.x.degree=3,
                 J.x.segments=1,
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

    ## Per Chen, Christensen and Kankanala (2021), notation-wise, Psi
    ## is the bases matrix for X, B for W. Note that I append .x to
    ## Psi and .w to B for clarity. Per Chen, Christensen and
    ## Kankanala (2021), notation-wise, Psi is of dimension J and B of
    ## dimension K. For clarity I adopt K.w.degree/.segments
    ## (corresponding to B hence "K") and J.x.degree/.segments
    ## (corresponding to Psi hence "J"). Note K >= J is necessary
    ## (number of columns of B.w must be greater than Psi.x, i.e.,
    ## there must be at least as many "instruments" as "endogenous"
    ## predictors).

    if(K.w.degree+K.w.segments < J.x.degree+J.x.segments) stop("K.w.degree+K.w.segments must be >= J.x.degree+J.x.segments")

    ## We allow for degree 0 (constant functions), and also allow for
    ## additive and generalized polynomial bases (when these are used
    ## we append a vector of ones, i.e., a "constant" term)

    ## Generate basis functions for W

    if(K.w.degree==0) {
        B.w <- matrix(1,NROW(W),1)
    } else {
        B.w <- prod.spline(x=W,
                           K=cbind(rep(K.w.degree,NCOL(W)),rep(K.w.segments,NCOL(W))),
                           knots=knots,
                           basis=basis)

        if(basis!="tensor") B.w <- cbind(1,B.w)
    }

    ## Generate basis functions for X and its derivative function

    if(J.x.degree==0) {
        Psi.x <- matrix(1,NROW(X),1)
        Psi.x.deriv <- matrix(0,NROW(X),1)
    } else {
        Psi.x <- prod.spline(x=X,
                            K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                            knots=knots,
                            basis=basis)

        Psi.x.deriv <- prod.spline(x=X,
                                  K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                  knots=knots,
                                  basis=basis,
                                  deriv.index=deriv.index,
                                  deriv=deriv.order)

        if(basis!="tensor") Psi.x <- cbind(1,Psi.x)
        if(basis!="tensor") Psi.x.deriv <- cbind(0,Psi.x.deriv)
    }

    ## Generate the NPIV coefficient vector using Choleski
    ## decomposition (computationally efficient). We first compute an
    ## object that is reused twice to avoid unnecessary computation
    ## (Psi.xTB.wB.wTB.w.invB.w, defined in Equation (3))

    Psi.xTB.wB.wTB.w.invB.w <- t(Psi.x)%*%B.w%*%chol2inv(chol(t(B.w)%*%B.w,pivot=chol.pivot))%*%t(B.w)
    beta <- chol2inv(chol(Psi.xTB.wB.wTB.w.invB.w%*%Psi.x+diag(lambda,NCOL(Psi.x)),pivot=chol.pivot))%*%Psi.xTB.wB.wTB.w.invB.w%*%Y

    ## Compute the IV function and its derivative. If evaluation data
    ## for X is provided, use it.

    if(J.x.degree==0) {
        Psi.x.eval <- Psi.x <- matrix(1,NROW(X),1)
        Psi.x.deriv.eval <- Psi.x.deriv <- matrix(0,NROW(X),1)
    } else {
        Psi.x.eval <- Psi.x <- prod.spline(x=X,
                                           K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                           knots=knots,
                                           basis=basis)

        Psi.x.deriv.eval <- Psi.x.deriv <- prod.spline(x=X,
                                                       K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                                       knots=knots,
                                                       basis=basis,
                                                       deriv.index=deriv.index,
                                                       deriv=deriv.order)
        if(!is.null(X.eval)) {
            Psi.x.eval <- prod.spline(x=X,
                                      xeval=X.eval,
                                      K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                      knots=knots,
                                      basis=basis)

           Psi.x.deriv.eval <- prod.spline(x=X,
                                           xeval=X.eval,
                                           K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                           knots=knots,
                                           basis=basis,
                                           deriv.index=deriv.index,
                                           deriv=deriv.order)
        }
        if(basis!="tensor") Psi.x <- cbind(1,Psi.x)
        if(basis!="tensor") Psi.x.eval <- cbind(1,Psi.x.eval)
        if(basis!="tensor") Psi.x.deriv <- cbind(0,Psi.x.deriv)
        if(basis!="tensor") Psi.x.deriv.eval <- cbind(0,Psi.x.deriv.eval)
    }

    ## h and h.deriv are the IV function and its derivatives,
    ## respectively.

    h <- Psi.x.eval%*%beta
    h.deriv <- Psi.x.deriv.eval%*%beta

    ## Compute asymptotic standard errors for the IV function and its
    ## derivatives (Chen and Pouzo 2012?) This will change.

    U.hat <- Y-Psi.x%*%beta
    C <- (t(Psi.x)%*%B.w)
    PP.inv <- chol2inv(chol(t(B.w)%*%B.w,pivot=chol.pivot))
    P.U <- B.w*as.numeric(U.hat)
    mho <- C%*%PP.inv%*%t(P.U)%*%(P.U)%*%PP.inv%*%t(C)
    D.inv <- chol2inv(chol(C%*%PP.inv%*%t(C),pivot=chol.pivot))
    D.inv.mho.D.inv <- D.inv%*%mho%*%D.inv
    asy.se <- sqrt(diag(Psi.x.eval%*%D.inv.mho.D.inv%*%t(Psi.x.eval)))
    asy.se.deriv <- sqrt(diag(Psi.x.deriv.eval%*%D.inv.mho.D.inv%*%t(Psi.x.deriv.eval)))

    ## Return a list with various objects that might be of interest to
    ## the user

    return(list(h=h,
                h.deriv=h.deriv,
                h.asy.se=asy.se,
                h.deriv.asy.se=asy.se.deriv,
                deriv.index=deriv.index,
                deriv.order=deriv.order,
                K.w.degree=K.w.degree,
                K.w.segments=K.w.segments,
                J.x.degree=J.x.degree,
                J.x.segments=J.x.segments,
                beta=beta,
                B.w=B.w,
                Psi.x=Psi.x,
                Psi.x.deriv=Psi.x.deriv))

}
