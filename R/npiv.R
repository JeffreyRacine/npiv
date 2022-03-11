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

    ## Generate basis functions for X and its derivative function, and
    ## if passed evaluation data for X, basis functions for those as
    ## well

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
        if(basis!="tensor") {
          Psi.x <- cbind(1,Psi.x)
          Psi.x.eval <- cbind(1,Psi.x.eval)
          Psi.x.deriv <- cbind(0,Psi.x.deriv)
          Psi.x.deriv.eval <- cbind(0,Psi.x.deriv.eval)
        }
    }

    ## Generate the NPIV coefficient vector using Choleski
    ## decomposition (computationally efficient) and call this
    ## "beta". We first compute an object that is reused twice to
    ## avoid unnecessary computation (Psi.xTB.wB.wTB.w.invB.w, defined
    ## in Equation (3)). Note that we use Psi.x and B.x naming
    ## conventions for clarity, and the "T" and "inv" notation
    ## connotes "Transpose" and "Inverse", respectively. Intermediate
    ## results that are reused in some form are stored temporarily.

    Psi.xTB.wB.wTB.w.invB.w <- t(Psi.x)%*%B.w%*%chol2inv(chol(t(B.w)%*%B.w,pivot=chol.pivot))%*%t(B.w)
    beta <- chol2inv(chol(Psi.xTB.wB.wTB.w.invB.w%*%Psi.x+diag(lambda,NCOL(Psi.x)),pivot=chol.pivot))%*%Psi.xTB.wB.wTB.w.invB.w%*%Y
    ## Compute Hurvich, Siminoff & Tsai's corrected AIC criterion.
    trH <- dim(Psi.x)[2]
    aic.penalty <- (1+trH/length(Y))/(1-(trH+2)/length(Y))
    aic.c <- ifelse(aic.penalty > 0,
                    log(mean((Y-Psi.x%*%beta)^2)) + aic.penalty,
                    .Machine$double.xmax)

    ## Compute the IV function and its derivative, denoted h and
    ## h.deriv, respectively, computed on the evaluation data, if
    ## provided, otherwise Psi.x.eval is simply Psi.x by default.

    h <- Psi.x.eval%*%beta
    h.deriv <- Psi.x.deriv.eval%*%beta

    ## Compute asymptotic standard errors for the IV function and its
    ## derivatives (formulae of Chen and Pouzo 2012)

    U.hat <- Y-Psi.x%*%beta
    C <- t(Psi.x)%*%B.w
    B.w.TB.w.inv <- chol2inv(chol(t(B.w)%*%B.w,pivot=chol.pivot))
    P.U <- B.w*as.numeric(U.hat)
    rho <- C%*%B.w.TB.w.inv%*%t(P.U)%*%(P.U)%*%B.w.TB.w.inv%*%t(C)
    D.inv <- chol2inv(chol(C%*%B.w.TB.w.inv%*%t(C),pivot=chol.pivot))
    D.inv.rho.D.inv <- D.inv%*%rho%*%D.inv

    ## These are the n x n memory hogs... must be a more efficient
    ## method for their computation, leave for now

    asy.se <- sqrt(diag(Psi.x.eval%*%D.inv.rho.D.inv%*%t(Psi.x.eval)))
    asy.se.deriv <- sqrt(diag(Psi.x.deriv.eval%*%D.inv.rho.D.inv%*%t(Psi.x.deriv.eval)))

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
                Psi.x.deriv=Psi.x.deriv,
                residuals.sample=Y-Psi.x%*%beta,
                residuals.eval=Y-Psi.x.eval%*%beta,
                AIC.c=aic.c,
                R2=RSQfunc(Y, Psi.x%*%beta)))

}

## Stripped down version that computes the bare minimum necessary to
## compute the AIC.c function of Hurvich et al (1998). It returns the
## AIC criterion and estimated IV function h only for the sample
## realizations (the latter for, say, Monte Carlo simulations
## etc.). You can pass any combination of spline degrees and knots
## (#knots=#segments+1).

npivaic <- function(Y,
                    X,
                    W,
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

    ## Generate basis functions for X

    if(J.x.degree==0) {
        Psi.x <- matrix(1,NROW(X),1)
    } else {
        Psi.x <- prod.spline(x=X,
                            K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                            knots=knots,
                            basis=basis)

        if(basis!="tensor") Psi.x <- cbind(1,Psi.x)
    }

    ## Generate the NPIV coefficient vector using Choleski
    ## decomposition (computationally efficient) and call this
    ## "beta". We first compute an object that is reused twice to
    ## avoid unnecessary computation (Psi.xTB.wB.wTB.w.invB.w, defined
    ## in Equation (3)). Note that we use Psi.x and B.x naming
    ## conventions for clarity, and the "T" and "inv" notation
    ## connotes "Transpose" and "Inverse", respectively. Intermediate
    ## results that are reused in some form are stored temporarily.

    Psi.xTB.wB.wTB.w.invB.w <- t(Psi.x)%*%B.w%*%chol2inv(chol(t(B.w)%*%B.w,pivot=chol.pivot))%*%t(B.w)
    h <- Psi.x%*%(chol2inv(chol(Psi.xTB.wB.wTB.w.invB.w%*%Psi.x+diag(lambda,NCOL(Psi.x)),pivot=chol.pivot))%*%Psi.xTB.wB.wTB.w.invB.w%*%Y)
    ## Compute Hurvich, Siminoff & Tsai's corrected AIC criterion.
    trH <- dim(Psi.x)[2]
    aic.penalty <- (1+trH/length(Y))/(1-(trH+2)/length(Y))
    aic.c <- ifelse(aic.penalty > 0,
                    log(mean((Y-h)^2)) + aic.penalty,
                    .Machine$double.xmax)

    return(list(h=h,
                K.w.degree=K.w.degree,
                K.w.segments=K.w.segments,
                J.x.degree=J.x.degree,
                J.x.segments=J.x.segments,
                AIC.c=aic.c))

}
