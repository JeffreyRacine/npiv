## We support multivariate numeric predictors. Multivariate splines
## are additive and tensor. Can be used for training or evaluation
## data. For numeric datatypes uses the B-spline. It also computes
## derivatives for the continuous variables of arbitrary order (issues
## warning when order exceeds degree of spline) with interaction if
## specified.

## Complexity can be modified via the number of knots (segments) or
## the spline degree (degree)

prodspline <- function(x,
                        z=NULL,
                        K=NULL,
                        I=NULL,
                        xeval=NULL,
                        zeval=NULL,
                        knots=c("quantiles","uniform"),
                        basis=c("additive","tensor","glp"),
                        x.min=NULL,
                        x.max=NULL,
                        deriv.index=1,
                        deriv=0) {

  basis <- match.arg(basis)
  knots <- match.arg(knots)

  if(missing(x) || missing (K)) stop(" must provide x and K")
  if(!is.matrix(K)) stop(" K must be a two-column matrix")

  ## Additive and glp models have intercept=FALSE in gsl.bs but
  ## intercept=TRUE in lm()

  gsl.intercept <- ifelse(basis=="additive" || basis=="glp", FALSE, TRUE)

  ## Care in passing (extra cast) and ensure K is a matrix of integers
  ## (K contains the spline degree [integer] for each dimension in
  ## column 1 and segments-1 for each dimension in column 2).

  x <- as.matrix(x)
  K <- round(K)

  n <- NROW(x)
  num.x <- NCOL(x)
  num.K <- nrow(K)

  if(deriv < 0) stop(" deriv is invalid")
  if(deriv > K[deriv.index,1]) warning(" deriv order too large, result will be zero")
  if(deriv.index < 1 || deriv.index > num.x) stop(" deriv.index is invalid")

  if(!is.null(z)) {
    z <- data.frame(z)
    num.z <- NCOL(z)
    num.I <- NROW(I)
    if(!is.null(zeval)) {
      zeval <- data.frame(zeval)
    }
  }

  if(is.null(xeval)) {
    xeval <- as.matrix(x)
  } else {
    xeval <- as.matrix(xeval)
    if(NCOL(x)!=NCOL(xeval)) stop(" xeval must be of the same dimension as x")
  }

  if(num.K != num.x) stop(paste(" dimension of x and K incompatible (",num.x,",",num.K,")",sep=""))
  if(!is.null(z) && (num.I != num.z)) stop(paste(" dimension of z and I incompatible (",num.z,",",num.I,")",sep=""))

  if(any(K[,1] > 0)||any(I != 0)) {

    tp <- list()

    j <- 1
    for(i in 1:num.x) {
      if(K[i,1] > 0) {
         ## nbreak is K[i,2]+1
        if(knots=="uniform") {
          knots.vec <- NULL
        } else {
          ## quantile knots
          knots.vec <- as.numeric(quantile(x[,i,drop=FALSE],probs=seq(0,1,length=(K[i,2]+1))))
          ## Correct issue of repeated knots points caused by point
          ## mass data (e.g. knots will be c(0,0,0,1,5), repeated
          ## knots will throw off gsl.bs). This adds a trivial amount
          ## to each knot and is only needed by gsl.bs(). Otherwise we
          ## retain only the unique points but then the dimension of
          ## the spline changes which can throw off predict etc. Note
          ## - there is something odd about what is produced by
          ## quantile as unique does not work as expected. 1e-20 is
          ## too small, 1e-10 works.
          knots.vec <- knots.vec + seq(0,1e-10*(max(x[,i,drop=FALSE])-min(x[,i,drop=FALSE])),length=length(knots.vec))
        }
        if((i==deriv.index)&&(deriv!=0)) {
          tp[[j]] <- predict(gsl.bs(x[,i,drop=FALSE],
                                    degree=K[i,1],
                                    nbreak=(K[i,2]+1),
                                    knots=knots.vec,
                                    deriv=deriv,
                                    x.min=if(is.null(x.min)){NULL}else{x.min[i]},
                                    x.max=if(is.null(x.max)){NULL}else{x.max[i]},
                                    intercept=gsl.intercept),newx=xeval[,i,drop=FALSE])
        } else {
          tp[[j]] <- predict(gsl.bs(x[,i,drop=FALSE],
                                    degree=K[i,1],
                                    nbreak=(K[i,2]+1),
                                    knots=knots.vec,
                                    x.min=if(is.null(x.min)){NULL}else{x.min[i]},
                                    x.max=if(is.null(x.max)){NULL}else{x.max[i]},
                                    intercept=gsl.intercept),newx=xeval[,i,drop=FALSE])
        }
        j <- j+1
      }
    }

    if(!is.null(z)) for(i in 1:num.z) {
      if(I[i] == 1) {
        if(is.null(zeval)) {
          tp[[j]] <- model.matrix(~z[,i])[,-1,drop=FALSE]
        } else {
          tp[[j]] <- model.matrix(~zeval[,i])[,-1,drop=FALSE]
        }
        j <- j+1
      }
    }

    ## When more than one element of K[,1] > 0 or I > 0 take all bases
    ## plus tensor product (all interactions), otherwise just the
    ## original bases for the one variable.

    if(NROW(tp) > 1) {
      ## First create all basis matrices for all continuous predictors
      ## (in essence, additive by default)
      P <- tp[[1]]
      for(i in 2:NROW(tp)) P <- cbind(P,tp[[i]])
      dim.P.no.tensor <- NCOL(P)
      ## Solely tensor if basis==tensor
      if(basis=="tensor") P <- tensor.prod.model.matrix(tp)
      if(basis=="glp") {
        P <- glp.model.matrix(tp)
        if(deriv!=0) {
          P.deriv <- list()
          for(i in 1:length(tp)) P.deriv[[i]] <- matrix(0,1,ncol(tp[[i]]))
          deriv.index <- deriv.index - length(which((K[1:deriv.index,1]==0)))
          while(deriv.index<=0) deriv.index <- deriv.index + 1
          P.deriv[[deriv.index]] <- matrix(NA,1,ncol(tp[[deriv.index]]))
          P[,!is.na(as.numeric(glp.model.matrix(P.deriv)))] <- 0
        }
      }
    } else {
      P <- tp[[1]]
      dim.P.no.tensor <- NCOL(P)
    }

  } else {

    ## No relevant continuous or discrete predictors.
    dim.P.no.tensor <- 0
    P <- matrix(rep(1,num.x),num.x,1)

  }

  attr(P,"dim.P.no.tensor") <- dim.P.no.tensor

  return(P)

}
