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
                J.w.segments.set=tmp1$J.w.segments.set,
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
              J.w.seg=tmp2$J.w.seg,
              theta.star=tmp2$theta.star))
  
  
}