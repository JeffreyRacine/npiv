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

    ## Some ideas for potential computational efficiency.
    ## Note we can also compute beta via lm.fit
    ## beta <- coef(lm.fit(fitted(lm.fit(B.w,Psi.x)),Y))

    ## To render as matrix equivalent to the above wrap in as.matrix
    ## and as.numeric
    ## b2sls <- as.matrix(as.numeric(coef(lm.fit(fitted(lm.fit(B.w,Psi.x)),Y))))

    ## We can get the covariance matrix via a call to lm() (the -1
    ## removes as added intercept which we don't want)
    ## vcov(lm(Y~fitted(lm.fit(B.w,Psi.x))-1))
    ## The question I need to answer is whether this renders the code
    ## more efficient or can be leveraged to get the correct standard
    ## errors, perhaps with some coaxing needed...

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

    ## These are the n x n memory hogs... the workaround is to provide
    ## X.eval with a reasonable number of grid points, say 100, then
    ## one can handle huge datasets without encountering this issue -
    ## this is noted in "Details" in npiv.Rd.

    asy.se <- sqrt(diag(Psi.x.eval%*%D.inv.rho.D.inv%*%t(Psi.x.eval)))
    asy.se.deriv <- sqrt(diag(Psi.x.deriv.eval%*%D.inv.rho.D.inv%*%t(Psi.x.deriv.eval)))

    ## Return a list with various objects that might be of interest to
    ## the user

    return(list(fitted=h,
                residuals=Y-Psi.x%*%beta,
                deriv=h.deriv,
                asy.se=asy.se,
                deriv.asy.se=asy.se.deriv,
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
                AIC.c=aic.c))

}

fitted.npiv <- function(object, ...){
   object$fitted
}

residuals.npiv <- function(object, ...) {
   object$residuals
}

## Function for determining optimal spline complexity for the
## nonparametric IV estimation of Chen, Christensen and Kankanala
## (2021). Notation-wise, I try to follow their notation closely, and
## append .x and .w where needed for clarity (described further
## below).

npivJ <- function(Y,
                  X,
                  W,
                  X.eval=NULL,
                  J.x.degree=3,
                  K.w.degree=3,
                  J.x.segments.set=c(2,4,8,16,32,64)-1,
                  K.w.segments.set=c(2,4,8,16,32,64)-1,
                  knots=c("uniform","quantiles"),
                  basis=c("tensor","additive","glp"),
                  eval.num=50,
                  boot.num=99,
                  alpha=0.05,
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
    if(alpha <=0 || alpha >=1) stop("alpha must lie in (0,1)")

    ## If specified, check that passed objects are of full rank

    if(check.is.fullrank) {
        if(!is.fullrank(Y)) stop("Y is not of full column rank")
        if(!is.fullrank(X)) stop("X is not of full column rank")
        if(!is.fullrank(W)) stop("W is not of full column rank")
    }

    ## Generate X.eval if not provided - pretty naive if multivariate
    ## X exists but easy on memory - user can pass a better object in
    ## such instances, so this is mainly for the univariate case (need
    ## to cast X as a matrix in case a vector is passed).

    if(is.null(X.eval)) {
        X.eval <- matrix(NA,eval.num,NCOL(X))
        for(nc in 1:NCOL(X)) X.eval[,nc] <- seq(min(as.matrix(X)[,nc]),max(as.matrix(X)[,nc]),length=eval.num)
    }

    ## See function npiv() in R code npiv.R for descriptions...

    ## Generate set of J1 J2 combinations given input
    ## J.x.segments.set. expand.grid() creates all possible
    ## combinations, then we select those with J2 > J1 (remove
    ## symmetric computations and those with J2=J1)

    J1.J2.x <- expand.grid(J.x.segments.set,J.x.segments.set)
    J1.J2.x <- subset(J1.J2.x,Var2 > Var1)

    J1.J2.w <- apply(J1.J2.x, c(1,2), function(u) K.w.segments.set[which(J.x.segments.set == u)])

    ## In what follows we loop over _rows_ of J1.J2 (makes for easy
    ##  parallelization if needed)

    Z.sup <- numeric()
    Z.sup.boot <- matrix(NA,boot.num,NROW(J1.J2.x))

    ## Temporary indication of where we are in the process

    pb <- progress_bar$new(format = " complexity determination [:bar] :percent eta: :eta",
                           clear = TRUE,
                           width= 60,
                           total = NROW(J1.J2.x))

    for(ii in 1:NROW(J1.J2.x)) {

        pb$tick()

        ## Temporary indication of where we are in the process

#        print(paste("Row ",ii," of ",NROW(J1.J2.x)))

        J.x.J1.segments <- J1.J2.x[ii,1]
        J.x.J2.segments <- J1.J2.x[ii,2]

        K.w.J1.segments <- J1.J2.w[ii,1]
        K.w.J2.segments <- J1.J2.w[ii,2]

        ## Tim precomputes and reuses basis functions... these are
        ## computationally efficient so might not save much time, morover
        ## if parallelized could be a waste and even add overhead

        ## Segments are set deterministically during search so makes
        ## sense to ensure degrees are set appropriately

        if(K.w.degree < J.x.degree) stop("K.w.degree must be >= J.x.degree")

        if(K.w.degree==0) {
            B.w.J1 <- B.w.J2 <- matrix(1,NROW(W),1)
        } else {
            B.w.J1 <- prod.spline(x=W,
                               K=cbind(rep(K.w.degree,NCOL(W)),rep(K.w.J1.segments,NCOL(W))),
                               knots=knots,
                               basis=basis)

            B.w.J2 <- prod.spline(x=W,
                               K=cbind(rep(K.w.degree,NCOL(W)),rep(K.w.J2.segments,NCOL(W))),
                               knots=knots,
                               basis=basis)

            if(basis!="tensor") {
                B.w.J1 <- cbind(1,B.w.J1)
                B.w.J2 <- cbind(1,B.w.J2)
            }
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

        B.w.J1.TB.w.J1.inv <- chol2inv(chol(t(B.w.J1)%*%B.w.J1,pivot=chol.pivot))
        B.w.J2.TB.w.J2.inv <- chol2inv(chol(t(B.w.J2)%*%B.w.J2,pivot=chol.pivot))

        Psi.xJ1TB.wB.wTB.w.invB.w <- t(Psi.x.J1)%*%B.w.J1%*%B.w.J1.TB.w.J1.inv%*%t(B.w.J1)
        tmp.J1 <- chol2inv(chol(Psi.xJ1TB.wB.wTB.w.invB.w%*%Psi.x.J1+diag(lambda,NCOL(Psi.x.J1)),pivot=chol.pivot))%*%Psi.xJ1TB.wB.wTB.w.invB.w
        beta.J1 <- tmp.J1%*%Y

        U.J1 <- Y-Psi.x.J1%*%beta.J1
        err.J1 <- Psi.x.J1.eval%*%tmp.J1%*%U.J1
        hhat.J1 <- Psi.x.J1.eval%*%beta.J1

        Psi.xJ2TB.wB.wTB.w.invB.w <- t(Psi.x.J2)%*%B.w.J2%*%B.w.J2.TB.w.J2.inv%*%t(B.w.J2)
        tmp.J2 <- chol2inv(chol(Psi.xJ2TB.wB.wTB.w.invB.w%*%Psi.x.J2+diag(lambda,NCOL(Psi.x.J2)),pivot=chol.pivot))%*%Psi.xJ2TB.wB.wTB.w.invB.w
        beta.J2 <- tmp.J2%*%Y

        U.J2 <- Y-Psi.x.J2%*%beta.J2
        err.J2 <- Psi.x.J2.eval%*%tmp.J2%*%U.J2
        hhat.J2 <- Psi.x.J2.eval%*%beta.J2

        ## Compute asymptotic variances and covariances for the IV
        ## functions - these will be memory intensive for large n as
        ## they require taking the diagonal of an n.eval x n.eval
        ## matrix (could be a way around but simply note that as it
        ## stands the computation of the standard errors by brute
        ## force needs to be addressed). Here n.eval is the number of
        ## rows in X.eval, and X.eval _must_ be supplied since X is
        ## allowed to be multivariate.

        CJ1 <- t(Psi.x.J1)%*%B.w.J1
        B.wUJ1 <- B.w.J1*as.numeric(U.J1)
        rho <- CJ1%*%B.w.J1.TB.w.J1.inv%*%t(B.wUJ1)%*%(B.wUJ1)%*%B.w.J1.TB.w.J1.inv%*%t(CJ1)
        D.J1.inv <- chol2inv(chol(CJ1%*%B.w.J1.TB.w.J1.inv%*%t(CJ1),pivot=chol.pivot))
        D.J1.inv.rho.D.J1.inv <- D.J1.inv%*%rho%*%D.J1.inv

        asy.var.J1 <- diag(Psi.x.J1.eval%*%D.J1.inv.rho.D.J1.inv%*%t(Psi.x.J1.eval))

        CJ2 <- t(Psi.x.J2)%*%B.w.J2
        B.wUJ2 <- B.w.J2*as.numeric(U.J2)
        rho <- CJ2%*%B.w.J2.TB.w.J2.inv%*%t(B.wUJ2)%*%(B.wUJ2)%*%B.w.J2.TB.w.J2.inv%*%t(CJ2)
        D.J2.inv <- chol2inv(chol(CJ2%*%B.w.J2.TB.w.J2.inv%*%t(CJ2),pivot=chol.pivot))
        D.J2.inv.rho.D.J2.inv <- D.J2.inv%*%rho%*%D.J2.inv

        asy.var.J2 <- diag(Psi.x.J2.eval%*%D.J2.inv.rho.D.J2.inv%*%t(Psi.x.J2.eval))

        ## Compute the covariance

        asy.cov.J1.J2 <- diag(Psi.x.J1.eval%*%D.J1.inv%*%CJ1%*%B.w.J1.TB.w.J1.inv%*%t(B.wUJ1)%*%(B.wUJ2)%*%B.w.J2.TB.w.J2.inv%*%t(CJ2)%*%t(D.J2.inv)%*%t(Psi.x.J2.eval))

        ## Compute the denominator of the t-stat (the numerator is the
        ## difference between err.J1 and err.J2)

        asy.se <- sqrt(asy.var.J1+asy.var.J2-2*asy.cov.J1.J2)

        ## The t-stat vector - we take the sup (max) of this to determine
        ## the optimal value of J (segments/knots of the Psi.x basis)

        Z.sup[ii] <- max((hhat.J1-hhat.J2)/asy.se)

        ## Bootstrap the sup t-stat, store in matrix Z.sup.boot, 1
        ## column per J1/J2 combination

        pbb <- progress_bar$new(format = " bootstrapping [:bar] :percent eta: :eta",
                               clear = TRUE,
                               width= 60,
                               total = boot.num)

        for(b in 1:boot.num) {
            pbb$tick()
            Z.sup.boot[b,ii] <- max((Psi.x.J1.eval%*%tmp.J1%*%(U.J1*rnorm(length(Y)))-
                                     Psi.x.J2.eval%*%tmp.J2%*%(U.J2*rnorm(length(Y))))/asy.se)
        }

    }

    ## Compute maximum over J set for each bootstrap draw (should
    ## produce a boot.num x 1 vector)

    Z.boot <- apply(Z.sup.boot, 1, max)

    theta.star <- quantile(Z.boot, 1 - alpha, names = FALSE)

    ## Compute Lepski choice

    num.J <- length(J.x.segments.set)

    test.mat <- matrix(NA, nrow = num.J, ncol = num.J)

    for(ii in 1:nrow(J1.J2.x)){
      row.index = which(J.x.segments.set == J1.J2.x[ii, 1])
      col.index = which(J.x.segments.set == J1.J2.x[ii, 2])
      test.mat[row.index, col.index] <- (Z.sup[ii] <= 1.1 * theta.star)
    }

    test.vec <- array(NA, dim = num.J)

    for(ii in 1:(num.J - 1)){
      test.vec[ii] <- prod(test.mat[ii, (ii + 1):num.J])
    }

    ## Add 1 to segments to get dimension

    if(all(test.vec == 0 || is.na(test.vec))){
      J.seg <- max(J.x.segments.set)
    } else {
      if(any(test.vec == 1)){
        J.seg <- J.x.segments.set[min(which(test.vec == 1))]
      }else{
        J.seg <- min(J.x.segments.set)
      }
    }
    J.hat <- (J.seg + J.x.degree)^NCOL(X)

    ## Compute truncated value (second-largest element of
    ## J.x.segments.set)

    J.seg.n <- max(J.x.segments.set[-which.max(J.x.segments.set)])
    J.hat.n <- (J.seg.n + J.x.degree)^NCOL(X)

    ## Take the minimum

    J.x.seg <- min(J.seg, J.seg.n)
    J.tilde <- min(J.hat, J.hat.n)
    K.w.seg <- K.w.segments.set[which(J.x.segments.set == J.x.seg)]

    ## Return a list with various objects that might be of interest to
    ## the user

    return(list(J.tilde=J.tilde,
                J.hat=J.hat,
                J.hat.n=J.hat.n,
                J.x.seg=J.x.seg,
                K.w.seg=K.w.seg,
                theta.star=theta.star))

}

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
  K.w.segments.set <- (2^(1:L.max+K.w.smooth))-1

  ## In what follows we loop over _rows_  (makes for easy
  ##  parallelization if needed)

  test.val <- array(NA, dim = L.max)

  ## Temporary indication of where we are in the process

  pb <- progress_bar$new(format = " grid determination [:bar] :percent eta: :eta",
                         clear = TRUE,
                         width= 60,
                         total = L.max)

#  print("Determining grid")

  for(ii in 1:L.max) {

    pb$tick()

    ## Temporary indication of where we are in the process

#    print(paste("Row ",ii," of ",L.max))

    if((ii <= 2) || ((ii > 2) & (test.val[ii-2] <= 10*sqrt(NROW(X))))){

      J.x.segments <- J.x.segments.set[ii]
      K.w.segments <- K.w.segments.set[ii]

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
                               K=cbind(rep(K.w.degree,NCOL(W)),rep(K.w.segments,NCOL(W))),
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
  K.w.segments.set <- K.w.segments.set[1:L.hat.max]

  J.hat.max <- (max(J.x.segments.set) + J.x.degree)^NCOL(X)

  alpha.hat <- min(0.5, 1/J.hat.max)

  ## Return a list with various objects that might be of interest to
  ## the user

  return(list(J.x.segments.set=J.x.segments.set,
              K.w.segments.set=K.w.segments.set,
              J.hat.max=J.hat.max,
              alpha.hat=alpha.hat))

}

## Function for determining optimal spline complexity for the
## nonparametric IV estimation of Chen, Christensen and Kankanala
## (2021). Notation-wise, I try to follow their notation closely, and
## append .x and .w where needed for clarity (described further
## below).

npiv_choose_J <- function(Y,
                          X,
                          W,
                          X.eval=NULL,
                          J.x.degree=3,
                          K.w.degree=3,
                          K.w.smooth=1,
                          knots=c("uniform","quantiles"),
                          basis=c("tensor","additive","glp"),
                          eval.num=50,
                          boot.num=99,
                          check.is.fullrank=FALSE,
                          chol.pivot=FALSE,
                          lambda=sqrt(.Machine$double.eps)) {

  ## Compute \hat{J}_max and data-determined grid of J values for X and W

  tmp1 <- npiv_Jhat_max(X,
                        W,
                        J.x.degree,
                        K.w.degree,
                        K.w.smooth,
                        knots,
                        basis,
                        check.is.fullrank,
                        chol.pivot)

  ## Compute data-driven choice via bootstrap

  tmp2 <- npivJ(Y,
                X,
                W,
                X.eval,
                J.x.degree,
                K.w.degree,
                J.x.segments.set=tmp1$J.x.segments.set,
                K.w.segments.set=tmp1$K.w.segments.set,
                knots,
                basis,
                eval.num,
                boot.num,
                alpha=tmp1$alpha.hat,
                check.is.fullrank,
                chol.pivot,
                lambda)

  ## Return a list with various objects that might be of interest to
  ## the user

  return(list(J.hat.max=tmp1$J.hat.max,
              J.hat.n=tmp2$J.hat.n,
              J.hat=tmp2$J.hat,
              J.tilde=tmp2$J.tilde,
              J.x.seg=tmp2$J.x.seg,
              K.w.seg=tmp2$K.w.seg,
              theta.star=tmp2$theta.star))

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

    return(list(fitted=h,
                K.w.degree=K.w.degree,
                K.w.segments=K.w.segments,
                J.x.degree=J.x.degree,
                J.x.segments=J.x.segments,
                AIC.c=aic.c))

}
