npiv <- function(Y,
                 X,
                 W,
                 X.eval=NULL,
                 deriv.index=1,
                 deriv.order=1,
                 W.degree=3,
                 W.segments=1,
                 X.degree=3,
                 X.segments=1,
                 knots=c("uniform","quantiles"),
                 basis=c("tensor","additive","glp"),
                 check.is.fullrank=FALSE,
                 chol.pivot=FALSE,
                 lambda=sqrt(.Machine$double.eps)) {

    ## Match arguments

    basis <- match.arg(basis)
    knots <- match.arg(knots)

    ## Basic error checking for invalid input

    if(missing(Y)) stop(" must provide Y")
    if(missing(X)) stop(" must provide X")
    if(missing(W)) stop(" must provide W")
    if(W.degree < 0) stop("W.degree must be a non-negative integer")
    if(X.degree < 0) stop("X.degree must be a non-negative integer")
    if(W.segments <= 0) stop("W.segments must be a positive integer")
    if(X.segments <= 0) stop("X.segments must be a positive integer")

    if(W.degree+W.segments < X.degree+X.segments) stop("W.degree+W.segments must be >= X.degree+X.segments")

    if(check.is.fullrank) {
        if(!is.fullrank(Y)) stop("Y is not of full column rank")
        if(!is.fullrank(X)) stop("X is not of full column rank")
        if(!is.fullrank(W)) stop("W is not of full column rank")
    }

    ## Compute the bases for W, X, and the requested derivative

    if(W.degree==0) {
        P.w <- matrix(1,NROW(W),1)
    } else {
        P.w <- prod.spline(x=W,
                           K=cbind(rep(W.degree,NCOL(W)),rep(W.segments,NCOL(W))),
                           knots=knots,
                           basis=basis)
        if(basis!="tensor") P.w <- cbind(1,P.w)
    }

    if(X.degree==0) {
        Q.x <- matrix(1,NROW(X),1)
        Q.x.deriv <- matrix(0,NROW(X),1)
    } else {
        Q.x <- prod.spline(x=X,
                            K=cbind(rep(X.degree,NCOL(X)),rep(X.segments,NCOL(X))),
                            knots=knots,
                            basis=basis)
        Q.x.deriv <- prod.spline(x=X,
                                  K=cbind(rep(X.degree,NCOL(X)),rep(X.segments,NCOL(X))),
                                  knots=knots,
                                  basis=basis,
                                  deriv.index=deriv.index,
                                  deriv=deriv.order)
        if(basis!="tensor") Q.x <- cbind(1,Q.x)
        if(basis!="tensor") Q.x.deriv <- cbind(0,Q.x.deriv)
    }

    Q.xTP.wP.wTP.w.invP.w <- t(Q.x)%*%P.w%*%chol2inv(chol(t(P.w)%*%P.w,pivot=chol.pivot))%*%t(P.w)
    beta <- chol2inv(chol(Q.xTP.wP.wTP.w.invP.w%*%Q.x+diag(lambda,NCOL(Q.x)),pivot=chol.pivot))%*%Q.xTP.wP.wTP.w.invP.w%*%Y

    ## Compute the IV function and its derivative. If evaluation data
    ## for X is provided, use it.

    if(X.degree==0) {
        Q.x.eval <- Q.x <- matrix(1,NROW(X),1)
        Q.x.deriv.eval <- Q.x.deriv <- matrix(0,NROW(X),1)
    } else {
        Q.x.eval <- Q.x <- prod.spline(x=X,
                                         K=cbind(rep(X.degree,NCOL(X)),rep(X.segments,NCOL(X))),
                                         knots=knots,
                                         basis=basis)
        Q.x.deriv.eval <- Q.x.deriv <- prod.spline(x=X,
                                                     K=cbind(rep(X.degree,NCOL(X)),rep(X.segments,NCOL(X))),
                                                     knots=knots,
                                                     basis=basis,
                                                     deriv.index=deriv.index,
                                                     deriv=deriv.order)
        if(!is.null(X.eval)) {
            Q.x.eval <- prod.spline(x=X,
                                     xeval=X.eval,
                                     K=cbind(rep(X.degree,NCOL(X)),rep(X.segments,NCOL(X))),
                                     knots=knots,
                                     basis=basis)
            Q.x.deriv.eval <- prod.spline(x=X,
                                          xeval=X.eval,
                                          K=cbind(rep(X.degree,NCOL(X)),rep(X.segments,NCOL(X))),
                                          knots=knots,
                                          basis=basis,
                                          deriv.index=deriv.index,
                                          deriv=deriv.order)
        }
        if(basis!="tensor") Q.x <- cbind(1,Q.x)
        if(basis!="tensor") Q.x.eval <- cbind(1,Q.x.eval)
        if(basis!="tensor") Q.x.deriv <- cbind(0,Q.x.deriv)
        if(basis!="tensor") Q.x.deriv.eval <- cbind(0,Q.x.deriv.eval)
    }

    h <- Q.x.eval%*%beta
    h.deriv <- Q.x.deriv.eval%*%beta

    ## Compute asymptotic standard errors for the IV function and its
    ## derivatives (Chen and Pouzo 2012?)

    U.hat <- Y-Q.x%*%beta
    C <- (t(Q.x)%*%P.w)
    PP.inv <- chol2inv(chol(t(P.w)%*%P.w,pivot=chol.pivot))
    P.U <- P.w*as.numeric(U.hat)
    mho <- C%*%PP.inv%*%t(P.U)%*%(P.U)%*%PP.inv%*%t(C)
    D.inv <- chol2inv(chol(C%*%PP.inv%*%t(C),pivot=chol.pivot))
    D.inv.mho.D.inv <- D.inv%*%mho%*%D.inv
    asy.se <- sqrt(diag(Q.x.eval%*%D.inv.mho.D.inv%*%t(Q.x.eval)))
    asy.se.deriv <- sqrt(diag(Q.x.deriv.eval%*%D.inv.mho.D.inv%*%t(Q.x.deriv.eval)))

    ## Return a list with various objects that might be of interest to
    ## the user

    return(list(h=h,
                h.deriv=h.deriv,
                h.asy.se=asy.se,
                h.deriv.asy.se=asy.se.deriv,
                deriv.index=deriv.index,
                deriv.order=deriv.order,
                W.degree=W.degree,
                X.degree=X.degree,
                W.segments=W.segments,
                X.segments=X.segments,
                beta=beta,
                P.w=P.w,
                Q.x=Q.x,
                Q.x.deriv=Q.x.deriv))

}
