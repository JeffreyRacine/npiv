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
                  J.x.degree=3,
                  J.x.segments.set=c(2,4,8,16,32,64)-1,
                  knots=c("uniform","quantiles"),
                  basis=c("tensor","additive","glp"),
                  eval.num=50,
                  boot.num=99,
                  alpha=0.5,
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

    J1.J2 <- expand.grid(J.x.segments.set,J.x.segments.set)
    J1.J2 <- subset(J1.J2,Var2 > Var1)

    ## In what follows we loop over _rows_ of J1.J2 (makes for easy
    ##  parallelization if needed)

    Z.sup <- numeric()
    Z.sup.boot <- matrix(NA,boot.num,NROW(J1.J2))

    for(ii in 1:NROW(J1.J2)) {

        ## Temporary indication of where we are in the process

        print(paste("Row ",ii," of ",NROW(J1.J2)))

        J.x.J1.segments <- J1.J2[ii,1]
        J.x.J2.segments <- J1.J2[ii,2]

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
                               K=cbind(rep(K.w.degree,NCOL(W)),rep(J.x.J1.segments,NCOL(W))),
                               knots=knots,
                               basis=basis)

            B.w.J2 <- prod.spline(x=W,
                               K=cbind(rep(K.w.degree,NCOL(W)),rep(J.x.J2.segments,NCOL(W))),
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

        Psi.xJ2TB.wB.wTB.w.invB.w <- t(Psi.x.J2)%*%B.w.J2%*%B.w.J2.TB.w.J2.inv%*%t(B.w.J2)
        tmp.J2 <- chol2inv(chol(Psi.xJ2TB.wB.wTB.w.invB.w%*%Psi.x.J2+diag(lambda,NCOL(Psi.x.J2)),pivot=chol.pivot))%*%Psi.xJ2TB.wB.wTB.w.invB.w
        beta.J2 <- tmp.J2%*%Y

        U.J2 <- Y-Psi.x.J2%*%beta.J2
        err.J2 <- Psi.x.J2.eval%*%tmp.J2%*%U.J2

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

        Z.sup[ii] <- max((err.J1-err.J2)/asy.se)

        ## Bootstrap the sup t-stat, store in matrix Z.sup.boot, 1
        ## column per J1/J2 combination

        for(b in 1:boot.num) Z.sup.boot[b,ii] <- max((Psi.x.J1.eval%*%tmp.J1%*%(U.J1*rnorm(length(Y)))-
                                                      Psi.x.J2.eval%*%tmp.J2%*%(U.J2*rnorm(length(Y))))/asy.se)

    }

    ## Compute maximum over J set for each bootstrap draw (should
    ## produce a boot.num x 1 vector)

    Z.boot <- apply(Z.sup.boot, 1, max)
    
    theta.star <- quantile(Z.boot, 1 - alpha, names = FALSE)
    
    ## NOTE: the J.x.segments.set and alpha will be data-determined.
    ## For now just take J.x.segments.set as given.
    
    ## Compute Lepski choice
    
    num.J <- length(J.x.segments.set)
    
    test.mat <- matrix(NA, nrow = num.J, ncol = num.J)
    
    for(ii in 1:nrow(J1.J2)){
      row.index = which(J.x.segments.set == J1.J2[ii, 1])
      col.index = which(J.x.segments.set == J1.J2[ii, 2])
      test.mat[row.index, col.index] <- (Z.sup[ii] <= 1.1 * theta.star)
    }
    
    test.vec <- array(NA, dim = num.J)
    
    for(ii in 1:(num.J - 1)){
      test.vec[ii] <- prod(test.mat[ii, (ii + 1):num.J])
    }
    
    if(any(test.vec == 1)){
      J.hat <- J.x.segments.set[min(which(test.vec == 1))]
    }else{
      J.hat <- min(J.x.segments.set)
    }
    
    ## Compute truncated value (second-largest element of
    ## J.x.segments.set)
    
    J.hat.n <- max(J.x.segments.set[-which.max(J.x.segments.set)])
    
    ## Take the minimum
    
    J.tilde <- min(J.hat, J.hat.n)

    ## Return a list with various objects that might be of interest to
    ## the user

    return(list(J.tilde=J.tilde,
                J.hat=J.hat,
                J.hat.n=J.hat.n,
                theta.star=theta.star))

}
