## UCB construction per Chen, Christensen and Kankanala (2021)
## If sieve dimension is provided, will return UCBs of Chen and
## Christensen (Quantitative Economics, 2018)

npiv.ucb <- function(Y,
                     X,
                     W,
                     X.eval=NULL,
                     deriv.index=1,
                     deriv.order=1,
                     K.w.degree=4,
                     K.w.segments=NULL,
                     K.w.smooth=1,
                     J.x.degree=3,
                     J.x.segments=NULL,
                     boot.num=99,
                     eval.num=100,
                     knots=c("uniform","quantiles"),
                     basis=c("tensor","additive","glp"),
                     X.min=NULL,
                     X.max=NULL,
                     W.min=NULL,
                     W.max=NULL,
                     alpha=0.05,
                     random.seed=42,
                     check.is.fullrank=FALSE,
                     progress=TRUE) {
  
  ## Check for regression
  
  if(all(X == W)){
    K.w.degree <- J.x.degree
    K.w.smooth <- 0
    K.w.segments <- J.x.segments
  }
  
  ## Check if dimension is provided
  
  if(is.null(K.w.segments) || is.null(J.x.segments)) {
    
    ## Chen, Christensen, Kankanala (2021) construction
    
    ## Compute \hat{J}_max and data-determined grid of J values for X and W
    
    tmp <- npiv_choose_J(Y,
                         X,
                         W,
                         X.eval,
                         J.x.degree,
                         K.w.degree,
                         K.w.smooth,
                         knots,
                         basis,
                         X.min,
                         X.max,
                         W.min,
                         W.max,
                         eval.num,
                         boot.num,
                         random.seed,
                         check.is.fullrank,
                         progress)
    
    ## Save seed prior to setting for bootstrap
    
    if(exists(".Random.seed", .GlobalEnv)) {
      
      save.seed <- get(".Random.seed", .GlobalEnv)
      exists.seed = TRUE
      
    } else {
      
      exists.seed = FALSE
      
    }
    
    ## In what follows we loop over J.x.segments.set
    
    Z.sup.boot <- matrix(NA,boot.num,length(tmp$J.x.segments.set))
    Z.sup.boot.deriv <- matrix(NA,boot.num,length(tmp$J.x.segments.set))
    
    for(ii in 1:length(tmp$J.x.segments.set)) {
      
      J.x.segments <- tmp$J.x.segments.set[ii]
      K.w.segments <- tmp$K.w.segments.set[ii]
      
      ## Generate basis functions for W for J
      
      if(K.w.degree==0) {
        B.w.J <- matrix(1,NROW(W),1)
      } else {
        B.w.J <- prod.spline(x=W,
                             K=cbind(rep(K.w.degree,NCOL(W)),rep(K.w.segments,NCOL(W))),
                             knots=knots,
                             basis=basis,
                             x.min=W.min,
                             x.max=W.max)
        if(basis!="tensor") {
          B.w.J <- cbind(1,B.w.J)
        }
      }
      
      ## Generate basis functions for X for J
      
      if(J.x.degree==0) {
        Psi.x.eval <- Psi.x <- matrix(1,NROW(X),1)
        Psi.x.deriv.eval <- Psi.x.deriv <- matrix(0,NROW(X),1)
      } else {
        Psi.x.J.eval <- Psi.x.J <- prod.spline(x=X,
                                               K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                               knots=knots,
                                               basis=basis,
                                               x.min=X.min,
                                               x.max=X.max)
        Psi.x.J.deriv.eval <- prod.spline(x=X,
                                          K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                          knots=knots,
                                          basis=basis,
                                          deriv.index=deriv.index,
                                          deriv=deriv.order,
                                          x.min=X.min,
                                          x.max=X.max)
        if(!is.null(X.eval)) {
          Psi.x.J.eval <- prod.spline(x=X,
                                      xeval=X.eval,
                                      K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                      knots=knots,
                                      basis=basis,
                                      x.min=X.min,
                                      x.max=X.max)
          Psi.x.J.deriv.eval <- prod.spline(x=X,
                                            xeval=X.eval,
                                            K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                            knots=knots,
                                            basis=basis,
                                            deriv.index=deriv.index,
                                            deriv=deriv.order,
                                            x.min=X.min,
                                            x.max=X.max)
        }
        
        if(basis!="tensor") {
          Psi.x.J <- cbind(1,Psi.x.J)
          Psi.x.J.eval <- cbind(1,Psi.x.J.eval)
          Psi.x.J.deriv <- cbind(0,Psi.x.J.deriv)
          Psi.x.J.deriv.eval <- cbind(0,Psi.x.J.deriv.eval)
        }
        
      }
      
      ## Code re-use/storage where possible... generate all objects
      ## needed to compute the t-stat vector
      
      B.w.J.TB.w.J.inv <- ginv(t(B.w.J)%*%B.w.J)
      
      Psi.xJTB.wB.wTB.w.invB.w <- t(Psi.x.J)%*%B.w.J%*%B.w.J.TB.w.J.inv%*%t(B.w.J)
      tmp.J <- ginv(Psi.xJTB.wB.wTB.w.invB.w%*%Psi.x.J)%*%Psi.xJTB.wB.wTB.w.invB.w
      beta.J <- tmp.J%*%Y
      
      U.J <- Y-Psi.x.J%*%beta.J
      hhat.J <- Psi.x.J.eval%*%beta.J
      
      ## Compute asymptotic variances and covariances
      
      D.J.inv.rho.D.J.inv <- t(t(tmp.J) * as.numeric(U.J))%*%(t(tmp.J) * as.numeric(U.J))
      asy.se.J <- sqrt(rowSums((Psi.x.J.eval%*%D.J.inv.rho.D.J.inv)*Psi.x.J.eval))
      asy.se.J.deriv <- sqrt(rowSums((Psi.x.J.deriv.eval%*%D.J.inv.rho.D.J.inv)*Psi.x.J.deriv.eval))
      
      ## Save at appropriate dimension to avoid having to compute twice
      
      if(tmp$J.x.segments.set[ii] == tmp$J.x.seg){
        
        hhat <- Psi.x.J.eval%*%beta.J
        hhat.deriv <- Psi.x.J.deriv.eval%*%beta.J
        
        asy.se <- asy.se.J
        asy.se.deriv <- asy.se.J.deriv
        
      }
      
      ## Bootstrap the sup t-stat, store in matrix Z.sup.boot and Z.sup.boot.deriv
      
      pbb <- progress_bar$new(format = "  bootstrapping [:bar] :percent eta: :eta",
                              clear = TRUE,
                              width = 60,
                              total = boot.num)
      
      ## Set seed to ensure same bootstrap draws across J
      
      set.seed(random.seed)
      
      for(b in 1:boot.num) {
        if(progress) pbb$tick()
        boot.draws <- rnorm(length(Y))
        
        if(any(asy.se.J == 0)){
          Z.sup.boot[b,ii] <- max(abs(((Psi.x.J.eval%*%tmp.J%*%(U.J*boot.draws))  / asy.se.J)[-which(asy.se.J == 0)]))
        } else {
          Z.sup.boot[b,ii] <- max(abs((Psi.x.J.eval%*%tmp.J%*%(U.J*boot.draws))  / asy.se.J))
        }
        
        if(any(asy.se.J.deriv == 0)){
          Z.sup.boot.deriv[b,ii] <- max(abs(((Psi.x.J.eval%*%tmp.J%*%(U.J*boot.draws))  / asy.se.J.deriv)[-which(asy.se.J.deriv == 0)]))
        } else {
          Z.sup.boot.deriv[b,ii] <- max(abs((Psi.x.J.deriv.eval%*%tmp.J%*%(U.J*boot.draws))  / asy.se.J.deriv))
        }
        
      }
      
    }
    
    ## Restore seed
    
    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
    
    ## Compute maximum over J set for each bootstrap draw (should
    ## produce a boot.num x 1 vector)
    
    Z.boot <- apply(Z.sup.boot, 1, max)
    Z.boot.deriv <- apply(Z.sup.boot.deriv, 1, max)
    
    z.star <- quantile(Z.boot, 1 - alpha, type = 5, names = FALSE)
    z.star.deriv <- quantile(Z.boot.deriv, 1 - alpha, type = 5, names = FALSE)
    
    ## Compute UCB critical values
    
    cv <- z.star + 0.25*log(log(length(Y)))*tmp$theta.star
    cv.deriv <- z.star.deriv + 0.25*log(log(length(Y)))*tmp$theta.star
    
    ## Store for returning
    
    J.x.segments <- tmp$J.x.seg
    K.w.segments <- tmp$K.w.seg
    
  } else {
    
    ## Chen and Christensen (2018) construction
    
    Z.sup.boot <- numeric()
    Z.sup.boot.deriv <- numeric()
    
    ## Generate basis functions for W for J
    
    if(K.w.degree==0) {
      B.w.J <- matrix(1,NROW(W),1)
    } else {
      B.w.J <- prod.spline(x=W,
                           K=cbind(rep(K.w.degree,NCOL(W)),rep(K.w.segments,NCOL(W))),
                           knots=knots,
                           basis=basis,
                           x.min=W.min,
                           x.max=W.max)
      if(basis!="tensor") {
        B.w.J <- cbind(1,B.w.J)
      }
    }
    
    ## Generate basis functions for X for J
    
    if(J.x.degree==0) {
      Psi.x.eval <- Psi.x <- matrix(1,NROW(X),1)
      Psi.x.deriv.eval <- Psi.x.deriv <- matrix(0,NROW(X),1)
    } else {
      Psi.x.J.eval <- Psi.x.J <- prod.spline(x=X,
                                             K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                             knots=knots,
                                             basis=basis,
                                             x.min=X.min,
                                             x.max=X.max)
      Psi.x.J.deriv.eval <- prod.spline(x=X,
                                        K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                        knots=knots,
                                        basis=basis,
                                        deriv.index=deriv.index,
                                        deriv=deriv.order,
                                        x.min=X.min,
                                        x.max=X.max)
      if(!is.null(X.eval)) {
        Psi.x.J.eval <- prod.spline(x=X,
                                    xeval=X.eval,
                                    K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                    knots=knots,
                                    basis=basis,
                                    x.min=X.min,
                                    x.max=X.max)
        Psi.x.J.deriv.eval <- prod.spline(x=X,
                                          xeval=X.eval,
                                          K=cbind(rep(J.x.degree,NCOL(X)),rep(J.x.segments,NCOL(X))),
                                          knots=knots,
                                          basis=basis,
                                          deriv.index=deriv.index,
                                          deriv=deriv.order,
                                          x.min=X.min,
                                          x.max=X.max)
      }
      
      if(basis!="tensor") {
        Psi.x.J <- cbind(1,Psi.x.J)
        Psi.x.J.eval <- cbind(1,Psi.x.J.eval)
        Psi.x.J.deriv <- cbind(0,Psi.x.J.deriv)
        Psi.x.J.deriv.eval <- cbind(0,Psi.x.J.deriv.eval)
      }
      
    }
    
    ## Code re-use/storage where possible... generate all objects
    ## needed to compute the t-stat vector
    
    B.w.J.TB.w.J.inv <- ginv(t(B.w.J)%*%B.w.J)
    
    Psi.xJTB.wB.wTB.w.invB.w <- t(Psi.x.J)%*%B.w.J%*%B.w.J.TB.w.J.inv%*%t(B.w.J)
    tmp.J <- ginv(Psi.xJTB.wB.wTB.w.invB.w%*%Psi.x.J)%*%Psi.xJTB.wB.wTB.w.invB.w
    beta.J <- tmp.J%*%Y
    
    U.J <- Y-Psi.x.J%*%beta.J
    hhat.J <- Psi.x.J.eval%*%beta.J
    
    ## Compute asymptotic variances and covariances
    
    D.J.inv.rho.D.J.inv <- t(t(tmp.J) * as.numeric(U.J))%*%(t(tmp.J) * as.numeric(U.J))
    asy.se.J <- sqrt(rowSums((Psi.x.J.eval%*%D.J.inv.rho.D.J.inv)*Psi.x.J.eval))
    asy.se.J.deriv <- sqrt(rowSums((Psi.x.J.deriv.eval%*%D.J.inv.rho.D.J.inv)*Psi.x.J.deriv.eval))
    
    hhat <- Psi.x.J.eval%*%beta.J
    hhat.deriv <- Psi.x.J.deriv.eval%*%beta.J
    
    asy.se <- asy.se.J
    asy.se.deriv <- asy.se.J.deriv
    
    ## Bootstrap the sup t-stat, store in matrix Z.sup.boot and Z.sup.boot.deriv
    
    pbb <- progress_bar$new(format = "  bootstrapping [:bar] :percent eta: :eta",
                            clear = TRUE,
                            width = 60,
                            total = boot.num)
    
    ## Set seed to ensure same bootstrap draws across J
    
    set.seed(random.seed)
    
    for(b in 1:boot.num) {
      if(progress) pbb$tick()
      boot.draws <- rnorm(length(Y))
      
      if(any(asy.se == 0)){
        Z.sup.boot[b] <- max(abs(((Psi.x.J.eval%*%tmp.J%*%(U.J*boot.draws))  / asy.se)[-which(asy.se == 0)]))
      } else {
        Z.sup.boot[b] <- max(abs((Psi.x.J.eval%*%tmp.J%*%(U.J*boot.draws))  / asy.se))
      }
      
      if(any(asy.se.deriv == 0)){
        Z.sup.boot.deriv[b] <- max(abs(((Psi.x.J.eval%*%tmp.J%*%(U.J*boot.draws))  / asy.se.deriv)[-which(asy.se.deriv == 0)]))
      } else {
        Z.sup.boot.deriv[b] <- max(abs((Psi.x.J.deriv.eval%*%tmp.J%*%(U.J*boot.draws))  / asy.se.deriv))
      }
      
    }
    
    ## Restore seed
    
    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
    
    cv <- quantile(Z.boot, 1 - alpha, type = 5, names = FALSE)
    cv.deriv <- quantile(Z.boot.deriv, 1 - alpha, type = 5, names = FALSE)

  }
  
  ## Compute UCBs
  
  h.lower <- hhat - cv * asy.se
  h.upper <- hhat + cv * asy.se
  
  h.lower.deriv <- hhat.deriv - cv.deriv * asy.se.deriv
  h.upper.deriv <- hhat.deriv + cv.deriv * asy.se.deriv
  
  return(list(fitted=hhat,
              deriv=hhat.deriv,
              asy.se=asy.se,
              deriv.asy.se=asy.se.deriv,
              deriv.index=deriv.index,
              deriv.order=deriv.order,
              K.w.degree=K.w.degree,
              K.w.segments=K.w.segments,
              J.x.degree=J.x.degree,
              J.x.segments=J.x.segments))
  
}
