
#' Bayesian Spatial Regression with Flexible Error Distributions
#'
#' Fits a spatial regression model by Markov chain Monte Carlo (MCMC). The
#' model can include spatial dependence in the response and in the error
#' terms, heteroscedastic errors, and several symmetric or skewed error
#' distributions.
#'
#' The formula may contain one response, a first right-hand side for the
#' location model, and an optional second right-hand side for the log-scale
#' model, using the multi-part formula syntax supported by [Formula::Formula()].
#' If the second right-hand side is omitted, a constant scale is used.
#'
#' @param formula A model formula describing the response and location
#'   covariates, and optionally the scale covariates as a second right-hand
#'   side.
#' @param data An optional data frame containing the variables in the formula.
#'   If omitted, variables are looked for in the formula's environment.
#' @param W1 Optional spatial weights matrix for spatial dependence in the
#'   response. Must be a square matrix of dimension equal to the number of
#'   observations, with nonnegative entries, and must be symmetric or
#'   row-standardized.
#' @param W2 Optional spatial weights matrix for spatial dependence in the
#'   error terms. It has the same dimension and structural requirements as
#'   `W1`.
#' @param n.sim Number of MCMC draws to retain after burn-in and thinning.
#' @param n.burnin Number of initial MCMC iterations to discard.
#' @param n.thin Thinning interval for the MCMC chain.
#' @param prior A list for prior specifications. The current implementation
#'   initializes its prior settings internally; this argument is retained by
#'   the function interface.
#' @param dist Error distribution. One of `"Gaussian"`, `"Student-t"`,
#'   `"Laplace"`, `"Skew-normal"`, or `"Skew-Student-t"`.
#'
#' @return An object of class `"spreg"`, a list containing posterior draws,
#'   MCMC settings, model information, the model frame data, and auxiliary
#'   quantities used by other methods. Depending on the fitted model, the
#'   posterior draws can include `beta`, `gamma`, `rho`, `lambda`, `delta`,
#'   and `nu`; acceptance rates are also returned for selected updates.
#'
#' @details The location coefficients are represented by `beta`. The `gamma`
#'   coefficients govern the log error variance when scale covariates are
#'   supplied. The parameter `rho` describes spatial dependence in the
#'   response, and `lambda` describes spatial dependence in the error terms.
#'   The skewness parameter `delta` is returned for skewed distributions, and
#'   the degrees-of-freedom parameter `nu` is returned for Student-t-based
#'   distributions.
#'
#' @references LeSage, J. and Pace, R. K. (2009). *Introduction to Spatial
#'   Econometrics*. Chapman & Hall/CRC, Boca Raton.
#' @examples
#'
#' data(oldcol, package="spdep")
#' W2 <- W1 <- spdep::nb2mat(COL.nb)
#' fit1 <- spatialReg(CRIME ~ INC + HOVAL, W1=W1, W2=W2, data=COL.OLD,
#'                    n.sim=500, n.burnin=500, n.thin=2, dist="Gaussian")
#' fit2 <- update(fit1, dist="Skew-normal")
#'
#' ################## Model comparison using DIC and WAIC
#' DIC(fit1,fit2)
#' WAIC(fit1,fit2)
#'
#' ################## Posterior summaries of the chosen model
#' summary(fit1)
#'
#' ################## Gaussian Q-Q plot of quantile-type residuals
#' resi <- residuals(fit1)
#' dev.new()
#' qqnorm(resi, pch=20, col="blue")
#' abline(0,1,lty=3)
#'
#' ################## Computes the Geweke z-scores for convergence diagnostics
#' geweke_diagSPREG(fit1)
#'
#' ################## Computes the effective sample size
#' effectiveSize_SPREG(fit1)
#'
#'
#' @importFrom Formula Formula model.part
#' @importFrom coda mcmc as.mcmc effectiveSize geweke.diag
#' @importFrom statmod rinvgauss
#' @importFrom grDevices dev.interactive dev.new
#' @importFrom stats coef dnorm dt lm.fit median model.matrix optim pgamma pnorm pt qnorm quantile rexp rgamma rnorm runif sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom spdep nb2mat
#'
#' @export
spatialReg <- function(formula, data, W1=NULL, W2=NULL, n.sim=1000,
                       n.burnin=200, n.thin=1, prior=list(),
                       dist=c("Gaussian","Student-t","Laplace","Skew-normal","Skew-Student-t")){
  dist <- match.arg(dist)
  if(missing(data)) data <- environment(formula)
  mmf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mmf), 0)
  mmf <- mmf[c(1,m)]
  mmf$drop.unused.levels <- TRUE
  mmf[[1]] <- as.name("model.frame")
  mmf$formula <- Formula(formula)
  mmf <- eval(mmf, parent.frame())
  y <- as.matrix(model.part(Formula(formula), data = mmf, lhs=1))
  mx <- model.part(Formula(formula), data = mmf, rhs = 1, terms = TRUE)
  X <- model.matrix(mx, data = mmf)
  if(suppressWarnings(formula(Formula(formula),lhs=0,rhs=2)=="~0")){
    Z <- matrix(1,nrow(X),1)
    colnames(Z) <- "sigma^2 "
    hete <- FALSE
  }
  else{
    mz <- model.part(Formula(formula), data = mmf, rhs = 2, terms = TRUE)
    Z <- model.matrix(mz, data = mmf)
    hete <- TRUE
  }
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  argamma <- 0
  if(!is.null(W1)){
    if(!is.matrix(W1) || nrow(W1)!=ncol(W1)) stop("W1 must be a square matrix.", call.=FALSE)
    if(nrow(W1)!=n) stop("The dimensions of W1 do not match response variable.", call.=FALSE)
    if(any(W1<0)) stop("The elements of W1 must be non-negative.", call.=FALSE)
    w1ev <- eigen(W1,only.values=TRUE)$values
    if(any(abs(Im(w1ev)) > 1e-10)) stop("All eigenvalues of W1 must be real.", call. = FALSE)
    w1ev <- sort(w1ev)
    ar <- TRUE
    W1y <- W1%*%y
    arrho <- 0
  }else{
    W1 <- diag(n)
    ar <- FALSE
    W1y <- y
    w1ev <- rep(0,n)
  }
  if(!is.null(W2)){
    if(!is.matrix(W2) || nrow(W2)!=ncol(W2)) stop("W2 must be a square matrix", call.=FALSE)
    if(nrow(W2)!=n) stop("The dimensions of W2 do not match response variable", call.=FALSE)
    if(any(W2<0)) stop("The elements of W2 must be non-negative.", call.=FALSE)
    w2ev <- eigen(W2,only.values=TRUE)$values
    if(any(abs(Im(w2ev)) > 1e-10)) stop("All eigenvalues of W2 must be real.", call. = FALSE)
    w2ev <- sort(w2ev)
    em <- TRUE
    W2y <- W2%*%y
    W2X <- W2%*%X
    if(ar) W2W1y <- W2%*%W1%*%y else W2W1y <- W2y
    arlambda <- 0
  }else{
    W2 <- diag(n)
    em <- FALSE
    W2X <- X
    W2y <- y
    if(ar) W2W1y <- W1y else W2W1y <- y
    w2ev <- rep(0,n)
  }
  if(w1ev[1]<0) limr <- 0.99/w1ev[1] else limr <- -Inf
  if(w1ev[n]>0) limr <- c(limr,0.99/w1ev[n]) else limr <- c(limr,Inf)
  if(w2ev[1]<0) liml <- 0.99/w2ev[1] else liml <- -Inf
  if(w2ev[n]>0) liml <- c(liml,0.99/w2ev[n]) else liml <- c(liml,Inf)
  arnu <- 0
  prior <- priors(prior,dist,hete,em,ar)
  mu_beta <- matrix(prior$mubeta,p,1)
  V_beta <-	1/prior$sigma2beta
  if(hete){
    mu_gamma <- matrix(prior$mugamma,q,1)
    V_gamma <- 1/prior$sigma2gamma
  }else{
    alpha0 <- prior$shape/2
    rate0 <- prior$rate/2
  }
  if(ar){
    mu_rho <- prior$murho
    V_rho <- 1/prior$sigma2rho
  }
  if(em){
    mu_lambda <- prior$mulambda
    V_lambda <- 1/prior$sigma2lambda
  }
  if(dist %in% c("Skew-normal","Skew-Student-t")){
    mu_delta <- prior$mudelta
    V_delta <- 1/prior$sigma2delta
  }
  if(dist %in% c("Student-t","Skew-Student-t")){
    anu <- prior$anu
    bnu <- prior$bnu
  }
  rep <- n.sim*n.thin + n.burnin
  beta <- matrix(0,p,rep+1)
  gamma <- matrix(0,q,rep+1)
  temp <- lm.fit(y=y,x=X)
  beta[,1] <- coef(temp)
  gamma[,1] <- c(log(mean(temp$residuals^2)),rep(0,q-1))
  delta <- lambda <- rho <- matrix(0,1,rep+1)
  nu <- matrix(100,1,rep+1)
  u <- rep(1,n)
  if(dist %in% c("Student-t","Skew-Student-t")) u <- rgamma(n,shape=25,rate=25)
  if(dist=="Laplace") u <- 1/rexp(n,rate=0.5)
  V <- rep(sqrt(2/pi),n)
  keep <- seq(from=n.burnin+1, to=rep, by=n.thin)
  bar <- txtProgressBar(min=0, max=rep, initial=0, width=min(50,n.sim+n.burnin), char="+", style=3)

  for(s in 1:rep){
    sigma_inv_diag <- exp(-Z%*%gamma[,s])
    Xt <- matrix(sqrt(u*sigma_inv_diag),n,p)*(X-lambda[,s]*W2X)
    yt <- y-rho[,s]*W1y
    yt <- matrix(sqrt(u*sigma_inv_diag),n,1)*(yt-lambda[,s]*W2%*%yt-delta[,s]*V)
    Vb <- chol(crossprod(Xt) + V_beta*diag(p))
    mb <- chol2inv(Vb)%*%(crossprod(Xt,yt) + matrix(V_beta,p,1)*mu_beta)
    beta[,s+1] <- mb + backsolve(Vb,matrix(rnorm(p),p,1))
    etam <- X%*%beta[,s+1]

    et <- y-rho[,s]*W1y-etam
    et1 <- et-lambda[,s]*W2%*%et-delta[,s]*V
    if(hete){
      lpg <- function(x){
        sigma_inv_diag <- exp(-Z%*%x)
        -0.5*sum(et1^2*u*sigma_inv_diag+Z%*%x)-0.5*crossprod(sqrt(V_gamma)*(x-mu_gamma))
      }
      lpgg <- function(x) 0.5*crossprod(Z,et1^2*u*exp(-Z%*%x)-1)-V_gamma*(x-mu_gamma)
      outg <- optim(gamma[,s],lpg,lpgg,method="BFGS",control=list(fnscale=-1))
      gamma0 <- matrix(outg$par,q,1)
      Vg <- 0.5*crossprod(matrix(sqrt(et1^2*u*exp(-Z%*%gamma0)),n,q)*Z) + V_gamma*diag(q)
      gammag <- gamma0 + backsolve(chol(Vg),matrix(rnorm(q),q,1))
      R <- lpg(gammag) - lpg(gamma[,s]) -
        0.5*crossprod(gammag-gamma0,crossprod(Vg,gammag-gamma0)) +
        0.5*crossprod(gamma[,s]-gamma0,crossprod(Vg,gamma[,s]-gamma0))
      if(log(runif(1)) < R){
        gamma[,s+1] <- gammag
        if(s > n.burnin) argamma <- argamma + 1
      }else gamma[,s+1] <- gamma[,s]
    }else gamma[,s+1] <- -log(rgamma(1,shape=(alpha0+n/2),rate=(sum(et^2)/2+rate0)))
    etag <- exp(-Z%*%gamma[,s+1])

    if(ar){
      et <- y-etam-lambda[1,s]*W2y+lambda[1,s]*W2X%*%beta[,s+1]-delta[,s]*V
      yt <- W1y - lambda[1,s]*W2W1y
      lpr <- function(x){
        sum(log(1-x*w1ev)) - 0.5*sum(u*etag*(et-x*yt)^2) - 0.5*(x-mu_rho)^2*V_rho
      }
      lprg <- function(x){
        -sum(w1ev/(1-x*w1ev)) + sum((et-x*yt)*u*etag*yt) - (x-mu_rho)*V_rho
      }
      outg <- optim(rho[,s],lpr,lprg,method="L-BFGS-B",lower=limr[1],upper=limr[2],control=list(fnscale=-1))
      rho0 <- outg$par
      Vr <- 1/(sum(w1ev^2/(1-rho0*w1ev)^2) + sum(u*etag*yt^2) + V_rho)
      limr2 <- pnorm(limr,mean=rho0,sd=sqrt(Vr))
      rhog <- qnorm(limr2[1] + runif(1)*(limr2[2]-limr2[1]),mean=rho0,sd=sqrt(Vr))
      R <- lpr(rhog) - lpr(rho[,s]) -
        0.5*(rhog-rho0)^2/Vr +
        0.5*(rho[,s]-rho0)^2/Vr
      if(log(runif(1)) < R){
        rho[,s+1] <- rhog
        if(s > n.burnin) arrho <- arrho + 1
      }else rho[,s+1] <- rho[,s]
    }

    if(em){
      et <- y - rho[,s+1]*W1y - etam
      lpl <- function(x){
        sum(log(1-x*w2ev)) - 0.5*sum(u*etag*(et-x*W2%*%et-delta[,s]*V)^2) - 0.5*(x-mu_lambda)^2*V_lambda
      }
      lplg <- function(x){
        -sum(w2ev/(1-x*w2ev)) + sum((et-x*W2%*%et-delta[,s]*V)*u*etag*(W2%*%et)) - (x-mu_lambda)*V_lambda
      }
      outg <- optim(lambda[,s],lpl,lplg,method="L-BFGS-B",lower=liml[1],upper=liml[2],control=list(fnscale=-1))
      lambda0 <- outg$par
      Vl <- 1/(sum(w2ev^2/(1-lambda0*w2ev)^2) + sum(u*etag*(W2y - rho[,s+1]*W2W1y - W2X%*%beta[,s+1])^2) + V_lambda)
      liml2 <- pnorm(liml,mean=lambda0,sd=sqrt(Vl))
      lambdag <- qnorm(liml2[1] + runif(1)*(liml2[2]-liml2[1]),mean=lambda0,sd=sqrt(Vl))
      R <- lpl(lambdag) - lpl(lambda[,s]) -	0.5*(lambdag-lambda0)^2/Vl + 0.5*(lambda[,s]-lambda0)^2/Vl
      if(log(runif(1)) < R){
        lambda[,s+1] <- lambdag
        if(s > n.burnin) arlambda <- arlambda + 1
      }else lambda[,s+1] <- lambda[,s]
    }

    if(dist %in% c("Skew-normal","Skew-Student-t")){
      et1 <- y - rho[,s+1]*W1y - lambda[,s+1]*W2y + lambda[,s+1]*rho[,s+1]*W2W1y - etam + lambda[,s+1]*W2X%*%beta[,s+1]
      muV <- (delta[,s]*etag*et1)/(1 + delta[,s]^2*etag)
      VV <- 1/(u*(1 + delta[,s]^2*etag))
      lims <- pnorm(0,mean=muV,sd=sqrt(VV))
      V <- qnorm(pmin(pmax(lims + runif(n)*(1-lims),.Machine$double.eps),1-.Machine$double.eps),mean=muV,sd=sqrt(VV))
      Vde <- 1/(sum(u*etag*V^2) + V_delta)
      mde <- (sum(u*etag*V*et1) + V_delta*mu_delta)*Vde
      delta[,s+1] <- rnorm(1,mean=mde,sd=sqrt(Vde))
    }

    if(dist=="Student-t"){
      et1 <- y - rho[,s+1]*W1y - lambda[,s+1]*W2y + lambda[,s+1]*rho[,s+1]*W2W1y - etam + lambda[,s+1]*W2X%*%beta[,s+1] - delta[,s+1]*V
      u <- rgamma(n,shape=(nu[,s]+rep(1,n))/2,rate=(nu[,s]+et1^2*etag)/2)
    }
    if(dist=="Skew-Student-t"){
      et1 <- y - rho[,s+1]*W1y - lambda[,s+1]*W2y + lambda[,s+1]*rho[,s+1]*W2W1y - etam + lambda[,s+1]*W2X%*%beta[,s+1] - delta[,s+1]*V
      u <- rgamma(n,shape=(nu[,s]+rep(2,n))/2,rate=(nu[,s]+V^2+et1^2*etag)/2)
    }
    if(dist=="Laplace") u <- rinvgauss(n,mean=1/sqrt(et1^2*etag),shape=1)

    if(dist %in% c("Student-t","Skew-Student-t")){
      lpnu <- function(x) (n*x/2)*log(x/2) - n*lgamma(x/2) + (x/2)*sum(log(u)-u) + (anu-1)*log(x) - bnu*x
      lpnug <- function(x) (n/2)*(log(x/2) + 1 - digamma(x/2)) + 0.5*sum(log(u)-u) + (anu-1)/x - bnu
      outg <- optim(nu[,s],lpnu,lpnug,method="L-BFGS-B",lower=.Machine$double.eps,upper=Inf,control=list(fnscale=-1))
      nu0 <- outg$par
      Vnu <- 1/((n/4)*trigamma(nu0/2) - n/(2*nu0) + (anu-1)/nu0^2)
      limnu <- pnorm(0,mean=nu0,sd=sqrt(Vnu))
      nug <- qnorm(limnu[1] + runif(1)*(1-limnu[1]),mean=nu0,sd=sqrt(Vnu))
      R <- lpnu(nug) - lpnu(nu[,s]) - 0.5*(nug-nu0)^2/Vnu + 0.5*(nu[,s]-nu0)^2/Vnu
      if(log(runif(1)) < R){
        nu[,s+1] <- nug
        if(s > n.burnin) arnu <- arnu + 1
      }else nu[,s+1] <- nu[,s]
    }
    setTxtProgressBar(bar,s)
  }
  close(bar)
  out_ <- list(beta=matrix(beta[,keep+1],p,n.sim),gamma=matrix(gamma[,keep+1],q,n.sim))
  rm(beta,gamma)
  if(hete){
    out_$argamma <- 100*argamma/(n.sim*n.thin)
    rm(argamma)
  }
  if(em){
    out_$lambda <- matrix(lambda[,keep+1],1,n.sim)
    out_$arlambda <- 100*arlambda/(n.sim*n.thin)
    rm(lambda,arlambda)
  }
  if(ar){
    out_$rho <- matrix(rho[,keep+1],1,n.sim)
    out_$arrho <- 100*arrho/(n.sim*n.thin)
    rm(rho,arrho)
  }
  if(dist %in% c("Student-t","Skew-Student-t")){
    out_$nu <- matrix(nu[,keep+1],1,n.sim)
    out_$arnu <- 100*arnu/(n.sim*n.thin)
    rm(nu,arnu)
  }
  if(dist %in% c("Skew-normal","Skew-Student-t")){
    out_$delta <- matrix(delta[,keep+1],1,n.sim)
    rm(delta)
  }
  chains <- out_
  out_ <- list(chains=chains,n.sim=n.sim,n.burnin=n.burnin,n.thin=n.thin,n=nrow(X),
               hete=hete,em=em,ar=ar,dist=dist,colx=colnames(X),colz=colnames(Z),call=match.call(),
               data=list(y=y,X=X,Z=Z),misc=list(W1y=W1y,W2y=W2y,W2W1y=W2W1y,W2X=W2X,w1ev=w1ev,w2ev=w2ev))
  rm(W1,W2,W1y,W2X,X,Z,y,u,V,temp,mu_beta,V_beta,w1ev,w2ev)
  rm(s,keep)
  if(hete) rm(hete,mu_gamma,V_gamma) else rm(hete,alpha0,rate0)
  if(ar) rm(ar,mu_rho,V_rho)
  if(em) rm(em,mu_lambda,V_lambda)
  if(dist %in% c("Student-t","Skew-Student-t")) rm(anu,bnu)
  if(!is.null(W2y)) rm(W2y)
  if(!is.null(W2W1y)) rm(W2W1y)
  class(out_) <- "spreg"
  return(out_)
}


#' Effective Sample Sizes for Spatial Regression MCMC Draws
#'
#' Computes the effective sample size for each available parameter chain in
#' an object of class `"spreg"`.
#'
#' @param x An object of class `"spreg"`, typically returned by [spatialReg()].
#'
#' @return An object of class `"essspreg"`, a list containing effective sample
#'   sizes for the location coefficients (`beta`), scale coefficients
#'   (`gamma`), and any spatial-dependence, skewness, or degrees-of-freedom
#'   parameters present in the fitted model.
#'
#' @seealso [as.mcmc.spreg()], [coda::effectiveSize()]
#' @export
effectiveSize_SPREG <- function(x){
  # Check that the input object is of class 'spreg'
  if(!inherits(x,"spreg")) stop("Only objects of class 'spreg' are supported!",call.=FALSE)
  # Convert the mtar object into an mcmc object
  x2 <- as.mcmc(x)
  x3 <- list()
  x3$beta <- matrix(effectiveSize(x2$beta),ncol=1,byrow=TRUE)
  rownames(x3$beta) <- colnames(x2$beta)
  colnames(x3$beta) <- ""
  x3$gamma <- matrix(effectiveSize(x2$gamma),ncol=1,byrow=TRUE)
  rownames(x3$gamma) <- colnames(x2$gamma)
  colnames(x3$gamma) <- ""
  if(!is.null(x2$rho)){
     x3$rho <- matrix(effectiveSize(x2$rho))
     rownames(x3$rho) <- colnames(x2$rho)
     colnames(x3$rho) <- ""
  }
  if(!is.null(x2$lambda)){
  	 x3$lambda <- matrix(effectiveSize(x2$lambda))
     rownames(x3$lambda) <- colnames(x2$lambda)
     colnames(x3$lambda) <- ""
  }
  if(!is.null(x2$delta)){
    x3$delta <- matrix(effectiveSize(x2$delta))
    rownames(x3$delta) <- colnames(x2$delta)
    colnames(x3$delta) <- ""
  }
  if(!is.null(x2$nu)){
  	x3$nu <- matrix(effectiveSize(x2$nu))
	rownames(x3$nu) <- colnames(x2$nu)
	colnames(x3$nu) <- ""
  }
  # Assign the output class
  class(x3) <- "essspreg"
  return(x3)
}



#' Geweke Convergence Diagnostics for Spatial Regression MCMC Draws
#'
#' Computes the Geweke diagnostic z-scores for the parameter chains in an
#' object of class `"spreg"`.
#'
#' @param x An object of class `"spreg"`, typically returned by [spatialReg()].
#' @param frac1 Fraction of the chain used for the first segment in the Geweke
#'   diagnostic. Defaults to `0.1`.
#' @param frac2 Fraction of the chain used for the second segment in the Geweke
#'   diagnostic. Defaults to `0.5`.
#'
#' @return An object of class `"gdspreg"`, a list containing diagnostic
#'   z-scores for the available parameters and the two fractions used.
#'
#' @seealso [coda::geweke.diag()]
#' @export
geweke_diagSPREG <- function(x,frac1=0.1,frac2=0.5){
  # Check that the input object is of class 'spreg'
  if(!inherits(x,"spreg")) stop("Only objects of class 'spreg' are supported!",call.=FALSE)
  # Convert the mtar object into an mcmc object
  x2 <- as.mcmc(x)
  x3 <- list()
  x3$beta <- matrix(geweke.diag(x2$beta,frac1=frac1,frac2=frac2)$z,ncol=1,byrow=TRUE)
  rownames(x3$beta) <- colnames(x2$beta)
  colnames(x3$beta) <- ""
  x3$gamma <- matrix(geweke.diag(x2$gamma,frac1=frac1,frac2=frac2)$z,ncol=1,byrow=TRUE)
  rownames(x3$gamma) <- colnames(x2$gamma)
  colnames(x3$gamma) <- ""
  if(!is.null(x2$rho)){
     x3$rho <- matrix(geweke.diag(x2$rho,frac1=frac1,frac2=frac2)$z)
     rownames(x3$rho) <- colnames(x2$rho)
     colnames(x3$rho) <- ""
  }
  if(!is.null(x2$lambda)){
  	 x3$lambda <- matrix(geweke.diag(x2$lambda,frac1=frac1,frac2=frac2)$z)
     rownames(x3$lambda) <- colnames(x2$lambda)
     colnames(x3$lambda) <- ""
  }
  if(!is.null(x2$delta)){
    x3$delta <- matrix(geweke.diag(x2$delta,frac1=frac1,frac2=frac2)$z)
    rownames(x3$delta) <- colnames(x2$delta)
    colnames(x3$delta) <- ""
  }
  if(!is.null(x2$nu)){
  	x3$nu <- matrix(geweke.diag(x2$nu,frac1=frac1,frac2=frac2)$z)
	rownames(x3$nu) <- colnames(x2$nu)
	colnames(x3$nu) <- ""
  }
  # Store the fractions used in the Geweke diagnostic
  x3$frac <- c(frac1,frac2)
  x3$regim <- x$ars$nregim
  # Assign the output class
  class(x3) <- "gdspreg"
  return(x3)
}


#' Summarize a Fitted Spatial Regression Model
#'
#' Calculates posterior summaries for a fitted spatial regression model,
#' including means, standard deviations, twice one minus the posterior
#' probability of direction, and highest posterior density (HPD) credible
#' intervals.
#'
#' @param object An object of class `"spreg"`, typically returned by
#'   [spatialReg()].
#' @param credible Probability content used for the HPD credible interval.
#'   Defaults to `0.95`.
#' @param ... Additional arguments passed to or from methods (currently unused).
#'
#' @return An object of class `"summary_spreg"` containing summaries of the
#'   location and scale coefficients and, when present, spatial dependence,
#'   skewness, and degrees-of-freedom parameters.
#'
#' @seealso [print.summary_spreg()], [spatialReg()]
#' @method summary spreg
#' @export
summary.spreg <- function(object, credible=0.95,...){
  # Number of MCMC simulations
  n.sim <- object$n.sim
  # Output list
  out_ <- list()
  # Internal helper function to compute summary statistics for MCMC draws
  resumen <- function(x){
        x <- matrix(x,ifelse(is.null(nrow(x)),1,nrow(x)),ifelse(is.null(ncol(x)),length(x),ncol(x)))
        # Initialize output: mean, sd, 2(1-PD), HPD.lower, HPD.upper
        y <- matrix(0,nrow(x),5)
        # Posterior means
        y[,1] <- rowMeans(x)
        y[,2] <- apply(x,1,sd)
        # Compute 2(1-PD): twice one minus posterior probability of direction
        y[,3] <- apply(x,1,function(x) min(mean(sign(median(x))*x > 0),1-1/200000))
        y[,3] <- 2*(1 - y[,3])
        # Build high-density interval (HPD) by scanning quantile pairs
        ks <- seq(credible,1,length=n.sim*(1-credible))
        lis <- t(apply(x,1,quantile,probs=ks-credible))
        lss <- t(apply(x,1,quantile,probs=ks))
        # Pick smallest-width interval
        dif <- apply(abs(lss-lis),1,which.min)
        y[,4] <- lis[cbind(1:nrow(x),dif)]
        y[,5] <- lss[cbind(1:nrow(x),dif)]
        # Naming columns
        colnames(y) <- c("   Mean","Std.Dev"," 2(1-PD) ","HPD.Lower","HPD.Upper")
        return(y)
  }
  out_$dist <- object$dist
  out_$sample.size <- object$n
  out_$location <- resumen(object$chains$beta)
  rownames(out_$location) <- object$colx
  maxc <- max(nchar(object$colx) + nchar(round(out_$location[,1],3)) - 5)

  if(object$hete) out_$scale <- resumen(object$chains$gamma)
  else out_$scale <- resumen(exp(object$chains$gamma))
  rownames(out_$scale) <- object$colz
  maxc <- c(maxc,max(nchar(object$colz) + nchar(round(out_$scale[,1],3)) - 5))
  if(object$ar){
	out_$rho <- resumen(object$chains$rho)
	rownames(out_$rho) <- "rho"
    maxc <- c(maxc,-2 + nchar(round(out_$rho[,1],3)))
  }
  if(object$em){
	out_$lambda <- resumen(object$chains$lambda)
	rownames(out_$lambda) <- "lambda"
    maxc <- c(maxc,1 + nchar(round(out_$lambda[,1],3)))
  }
  if(object$dist %in% c("Skew-normal","Skew-Student-t")){
 	out_$delta <- resumen(object$chains$delta)
	rownames(out_$delta) <- "delta"
    maxc <- c(maxc,nchar(round(out_$delta[,1],3)))
  }
  if(object$dist %in% c("Slash","Student-t","Skew-Student-t")){
 	out_$nu <- resumen(object$chains$nu)
	rownames(out_$nu) <- "nu"
    maxc <- c(maxc,-3 + nchar(round(out_$nu[,1],3)))
  }
  out_$maxc <- maxc
  # Assign summary class
  class(out_) <- "summary_spreg"
  return(out_)
}

#' Print a Spatial Regression Model Summary
#'
#' Prints the posterior summary produced by [summary.spreg()], displaying
#' location and scale coefficients and any spatial-dependence, skewness, or
#' degrees-of-freedom parameters available in the model.
#'
#' @param x An object of class `"summary_spreg"` returned by [summary.spreg()].
#' @param digits Number of decimal places used for numeric output. The default
#'   is `max(3, getOption("digits") - 2)`.
#' @param ... Additional arguments (currently unused).
#'
#' @return The function is called for its side effect of printing the summary;
#'   it does not return a value explicitly.
#'
#' @method print summary_spreg
#' @export
print.summary_spreg <- function(x, digits=max(3, getOption("digits") - 2),...){
  # Print sample size
  message("\n\nSample size       : ",x$sample.size)
  # Print error distribution
  message("\nError Distribution: ",x$dist,"\n")
  # Print location parameters summary
  message("\nLocation coefficients\n")
  maxc <- x$maxc
  rownames(x$location)[1] <- paste0(rownames(x$location)[1],paste0(rep(" ",max(maxc)-maxc[1]),collapse=""))
  print(round(x$location,digits=digits),na.print="   |   ")
  # Print scale parameters summary
  message("\nScale coefficients\n")
  rownames(x$scale)[1] <- paste0(rownames(x$scale)[1],paste0(rep(" ",max(maxc)-maxc[2]),collapse=""))
  print(round(x$scale,digits=digits),na.print="   |   ")
  cont <- 2
  # Print spatial dependence in the response variable
  if(!is.null(x$rho)){
	message("\nSpatial dependence in the response variable")
	cont <- cont + 1
    rownames(x$rho)[1] <- paste0(rownames(x$rho)[1],paste0(rep(" ",max(maxc)-maxc[cont]),collapse=""))
	print(round(x$rho,digits=digits),na.print="   |   ")
  }
  # Print spatial dependence in the error terms
  if(!is.null(x$lambda)){
    message("\nSpatial dependence in the error terms\n")
	cont <- cont + 1
    rownames(x$lambda)[1] <- paste0(rownames(x$lambda)[1],paste0(rep(" ",max(maxc)-maxc[cont]),collapse=""))
	print(round(x$lambda,digits=digits),na.print="   |   ")
  }
  # Print skewness parameter if present
  if(x$dist %in% c("Skew-normal","Skew-Student-t")){
	message("\n\nSkewness parameter","\n")
	cont <- cont + 1
    rownames(x$delta)[1] <- paste0(rownames(x$delta)[1],paste0(rep(" ",max(maxc)-maxc[cont]),collapse=""))
	print(round(x$delta,digits=digits),na.print="   .   ",)
  }
  # Print extra parameters if present
  if(x$dist %in% c("Slash","Student-t","Skew-Student-t")){
	message("\n\nExtra parameter","\n")
	cont <- cont + 1
    rownames(x$nu)[1] <- paste0(rownames(x$nu)[1],paste0(rep(" ",max(maxc)-maxc[cont]),collapse=""))
	print(round(x$nu,digits=digits),na.print="   .   ")
  }
  message("\n\n")
}


#' Convert Spatial Regression Draws to MCMC Objects
#'
#' Converts posterior draws stored in an object of class `"spreg"` into
#' `coda` MCMC objects, preserving parameter names and the burn-in/thinning
#' metadata.
#'
#' @param x An object of class `"spreg"`, typically returned by [spatialReg()].
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `"spregmcmc"`, a list containing `coda::mcmc`
#'   objects for the available parameter groups and MCMC metadata.
#'
#' @seealso [plot.spregmcmc()], [coda::mcmc()]
#' @method as.mcmc spreg
#' @export
as.mcmc.spreg <- function(x,...){
  # Number of MCMC simulations
  n.sim <- x$n.sim
  # Initialize output list
  out <- list()
  # Store location parameters as an mcmc object
  temp <- matrix(t(x$chains$beta),n.sim,length(x$colx))
  colnames(temp) <- x$colx
  out$beta <- mcmc(temp,thin=x$n.thin,start=x$n.burnin+1)
  # Store scale parameters as an mcmc object
  if(x$hete) temp <- matrix(t(x$chains$gamma),n.sim,length(x$colz))
  else temp <- matrix(t(exp(x$chains$gamma)),n.sim,length(x$colz))
  colnames(temp) <- x$colz
  out$gamma <- mcmc(temp,thin=x$n.thin,start=x$n.burnin+1)
  # Store AR parameters as an mcmc object if present
  if(x$ar){
    temp <- matrix(x$chains$rho,n.sim,1)
	colnames(temp) <- "rho"
    out$rho <- mcmc(temp,thin=x$n.thin,start=x$n.burnin+1)
  }
  # Store EM parameters as an mcmc object if present
  if(x$em){
    temp <- matrix(x$chains$lambda,n.sim,1)
	colnames(temp) <- "lambda"
    out$lambda <- mcmc(temp,thin=x$n.thin,start=x$n.burnin+1)
  }
  # Store extra parameters as an mcmc object if present
  if(!(x$dist %in% c("Skew-normal","Gaussian","Laplace"))){
    temp <- matrix(x$chains$nu,n.sim,1)
	colnames(temp) <- "nu"
    out$nu <- mcmc(temp,thin=x$n.thin,start=x$n.burnin+1)
  }
  # Store skewness parameters as an mcmc object if present
  if(x$dist %in% c("Skew-normal","Skew-Student-t")){
    temp <- matrix(x$chains$delta,n.sim,1)
	colnames(temp) <- "delta"
    out$delta <- mcmc(temp,thin=x$n.thin,start=x$n.burnin+1)
  }
  # MCMC metadata
  out$thin <- x$n.thin
  out$start <- x$n.burnin+1
  out$end <- max(seq(x$n.burnin + 1,length.out=x$n.sim,by=x$n.thin))
  out$sz <- x$n.sim
  # Assign class and return
  class(out) <- "spregmcmc"
  return(out)
}
#' Plot MCMC Chains for a Spatial Regression Model
#'
#' Produces trace and density plots for the parameter chains stored in an
#' object of class `"spregmcmc"`. Separate graphics devices are opened for
#' the available groups of parameters.
#'
#' @param x An object of class `"spregmcmc"`, typically returned by
#'   [as.mcmc.spreg()].
#' @param trace Logical; if `TRUE`, display trace plots.
#' @param density Logical; if `TRUE`, display density plots.
#' @param smooth Logical; if `TRUE`, request smoothed density estimates where
#'   supported by the plotting method.
#' @param bwf Optional bandwidth adjustment passed to the underlying plotting
#'   method.
#' @param auto.layout Logical; whether the underlying plotting method should
#'   choose a layout automatically.
#' @param ask Logical; whether to prompt before drawing the next plot.
#' @param ... Additional graphical arguments passed to the underlying plot
#'   methods.
#'
#' @return Called for its plotting side effects.
#'
#' @method plot spregmcmc
#' @export
plot.spregmcmc <- function(x, trace=TRUE, density=TRUE, smooth=FALSE, bwf=NULL, auto.layout=TRUE, ask=dev.interactive(), ...){
    # Open a new graphics device for location parameters
    dev.new()
    message("\nLocation parameters\n")
    plot(x$beta, trace=trace, density=density, smooth=smooth,
         auto.layout=auto.layout, ask=ask, ...)
    # Open a new graphics device for scale parameters
    dev.new()
    message("\nScale parameters\n")
    plot(x$gamma, trace=trace, density=density, smooth=smooth,
         auto.layout=auto.layout, ask=ask, ...)
	# Plot AR parameters if present
	if(!is.null(x$rho)){
	  dev.new()
	  message("\nSpatial dependence in the response variable\n")
	  plot(x$rho, trace=trace, density=density, smooth=smooth,
	         auto.layout=auto.layout, ask=ask, ...)
	}
	# Plot EM parameters if present
	if(!is.null(x$lambda)){
	  dev.new()
	  message("\nSpatial dependence in the error terms\n")
	  plot(x$lambda, trace=trace, density=density, smooth=smooth,
	         auto.layout=auto.layout, ask=ask, ...)
	}
	# Plot skewness parameters if present
	if(!is.null(x$delta)){
	  dev.new()
	  message("\nSkewness parameter\n")
	  plot(x$delta, trace=trace, density=density, smooth=smooth,
	         auto.layout=auto.layout, ask=ask, ...)
	}
	# Plot extra parameters if present
	if(!is.null(x$nu)){
	  dev.new()
	  message("\nExtra parameter\n")
	  plot(x$nu, trace=trace, density=density, smooth=smooth,
	         auto.layout=auto.layout, ask=ask, ...)
	}
}

#' Deviance Information Criterion for Spatial Regression Models
#'
#' Computes the Deviance Information Criterion (DIC) for one or more fitted
#' spatial regression models of class `"spreg"`.
#'
#' @param ... One or more fitted objects of class `"spreg"` to compare.
#'
#' @return For `DIC.spreg()`, a one-column numeric matrix containing
#'   the Deviance Information Criterion (DIC) for each fitted model.
#'   Row names identify the model expressions, and the column is named
#'   `"DIC"`. Smaller values indicate a better trade-off between model
#'   fit and complexity when comparing models fitted to the same data.
#'
#' @details The criterion is calculated from the posterior mean deviance and
#'   the average deviance over the retained MCMC draws. Smaller values indicate
#'   a better trade-off between fit and model complexity, when comparing models
#'   fitted to the same response data.
#'
#' @seealso [WAIC()]
#' @method DIC spreg
#' @export
DIC.spreg <- function(...){
  # Internal function to compute minus twice the log-likelihood
  mtlogLik <- function(dist,beta,gamma,rho,lambda,delta,nu){
      # Compute residuals
      resu <- y - rho*W1y - lambda*W2y + rho*lambda*W2W1y - (X%*%beta - lambda*W2X%*%beta)
	  sigma2 <- exp(Z%*%gamma)
	  if(dist %in% c("Skew-normal","Skew-Student-t")){
		  muv <- delta*resu/(delta^2+sigma2)
		  sigma2v <- sigma2/(delta^2+sigma2)
	  	  out <- resu^2/(delta^2+sigma2)
		  if(dist=="Skew-Student-t"){
			 muv <- muv*sqrt((nu+1)/(nu+out))
			 sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pt(-x/sqrt(sigma2v),df=nu+1,lower.tail=FALSE)))
			 #out <- 2*(1 + out/nu)^(-(nu+1)/2)*gamma((nu+1)/2)*sum0/(sqrt(nu*pi*(delta^2+sigma2))*gamma(nu/2))
			 out <- 2*sum0*dt(sqrt(out),df=nu)/sqrt(delta^2+sigma2)
		  }
		  if(dist=="Skew-normal"){
			 sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pnorm(-x/sqrt(sigma2v),lower.tail=FALSE)))
			 #out <- 2*exp(-out/2)*sum0/sqrt(2*pi*(delta^2+sigma2))
			 out <- 2*sum0*dnorm(sqrt(out))/sqrt(delta^2+sigma2)
		  }
      }else{# Case of symmetric distributions
           out <- resu^2/sigma2
           # Distribution-specific likelihood kernels
           out <- switch(dist,
                         "Gaussian"={exp(-out/2)},
                         "Laplace"={sqrt(pi/2)*exp(-sqrt(out))},
                         "Student-t"={(1 + out/nu)^(-(nu+1)/2)*gamma((nu+1)/2)/(sqrt(nu/2)*gamma(nu/2))},
                         "Slash"={ifelse(out==0,nu/(nu+1),gamma((1+nu)/2)*(nu/2)*pgamma(1,shape=(1+nu)/2,rate=out/2)/((out/2)^((1+nu)/2)))})
           # Normalizing constant
           out <- out/sqrt(2*pi*sigma2)
      }
      # Convert to deviance (-2 log-likelihood)
      out <- -2*sum(log(out))-2*sum(log(1-rho*w1ev))-2*sum(log(1-lambda*w2ev))
      return(out)
  }
  # Collect fitted mtar objects
  another <- list(...)
  call. <- match.call()
  # Initialize output container
  out <- matrix(0,length(another),1)
  outnames <- vector()
  # Loop over each fitted model
  for(l in 1:(length(another))){
      n.sim <- another[[l]]$n.sim
      dist <- another[[l]]$dist
	  y <- another[[l]]$data$y
	  X <- another[[l]]$data$X
	  Z <- another[[l]]$data$Z
	  W1y <- another[[l]]$misc$W1y
	  W2y <- another[[l]]$misc$W2y
	  W2W1y <- another[[l]]$misc$W2W1y
	  W2X <- another[[l]]$misc$W2X
	  w1ev <- another[[l]]$misc$w1ev
	  w2ev <- another[[l]]$misc$w2ev
	  Dbarv <- 0
      # Loop over MCMC iterations
      for(j in 1:n.sim){
	      betaj <- another[[l]]$chains$beta[,j]
		  gammaj <- another[[l]]$chains$gamma[,j]
		  if(is.null(another[[l]]$chains$rho)) rhoj <- 0 else rhoj <- another[[l]]$chains$rho[,j]
		  if(is.null(another[[l]]$chains$lambda)) lambdaj <- 0 else lambdaj <- another[[l]]$chains$lambda[,j]
          if(is.null(another[[l]]$chains$delta)) deltaj <- 0 else deltaj <- another[[l]]$chains$delta[,j]
          if(is.null(another[[l]]$chains$nu)) nuj <- 10000 else nuj <- another[[l]]$chains$nu[,j]
          Dbarv <- Dbarv + mtlogLik(dist,betaj,gammaj,rhoj,lambdaj,deltaj,nuj)
      }
      # Compute deviance at posterior mean parameters
	  betaj <- rowMeans(another[[l]]$chains$beta)
	  gammaj <- rowMeans(another[[l]]$chains$gamma)
	  if(is.null(another[[l]]$chains$rho)) rhoj <- 0 else rhoj <- mean(another[[l]]$chains$rho)
	  if(is.null(another[[l]]$chains$lambda)) lambdaj <- 0 else lambdaj <- mean(another[[l]]$chains$lambda)
      if(is.null(another[[l]]$chains$delta)) deltaj <- 0 else deltaj <- mean(another[[l]]$chains$delta)
      if(is.null(another[[l]]$chains$nu)) nuj <- 10000 else nuj <- mean(another[[l]]$chains$nu)
      Dhat <- mtlogLik(dist,betaj,gammaj,rhoj,lambdaj,deltaj,nuj)
      # Deviance Information Criterion
      out[l] <- Dhat + 2*(Dbarv/n.sim - Dhat)
      outnames[l] <- as.character(call.[l+1])
  }
  # Assign names and return DIC values
  rownames(out) <- outnames
  colnames(out) <- "DIC"
  return(out)
}

#' Watanabe-Akaike Information Criterion for Spatial Regression Models
#'
#' Computes the Watanabe-Akaike Information Criterion (WAIC) for one or more
#' fitted spatial regression models of class `"spreg"`.
#'
#' @param ... One or more fitted objects of class `"spreg"` to compare.
#'
#' @return For `WAIC.spreg()`, a one-column numeric matrix containing
#'   the Watanabe-Akaike Information Criterion (WAIC) for each fitted
#'   model. Row names identify the model expressions, and the column
#'   is named `"WAIC"`. Smaller values indicate better expected
#'   predictive performance when comparing models fitted to the same
#'   response data.
#'
#' @details The calculation uses the likelihood evaluated across retained
#'   posterior draws. Smaller values indicate a better expected predictive
#'   performance when models are fitted to the same response data.
#'
#' @seealso [DIC()]
#' @method WAIC spreg
#' @export
WAIC.spreg <- function(...){
  # Internal function to compute the likelihood
  Lik <- function(dist,beta,gamma,rho,lambda,delta,nu){
      # Compute residuals
      resu <- y - rho*W1y - lambda*W2y + rho*lambda*W2W1y - (X%*%beta - lambda*W2X%*%beta)
	  sigma2 <- exp(Z%*%gamma)
	  if(dist %in% c("Skew-normal","Skew-Student-t")){
		  muv <- delta*resu/(delta^2+sigma2)
		  sigma2v <- sigma2/(delta^2+sigma2)
	  	  out <- resu^2/(delta^2+sigma2)
		  if(dist=="Skew-Student-t"){
			 muv <- muv*sqrt((nu+1)/(nu+out))
			 sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pt(-x/sqrt(sigma2v),df=nu+1,lower.tail=FALSE)))
			 #out <- 2*(1 + out/nu)^(-(nu+1)/2)*gamma((nu+1)/2)*sum0/(sqrt(nu*pi*(delta^2+sigma2))*gamma(nu/2))
			 out <- 2*sum0*dt(sqrt(out),df=nu)/sqrt(delta^2+sigma2)
		  }
		  if(dist=="Skew-normal"){
			 sum0 <- apply(muv,1,function(x) max(.Machine$double.xmin,pnorm(-x/sqrt(sigma2v),lower.tail=FALSE)))
			 #out <- 2*exp(-out/2)*sum0/sqrt(2*pi*(delta^2+sigma2))
			 out <- 2*sum0*dnorm(sqrt(out))/sqrt(delta^2+sigma2)
		  }
      }else{# Case of symmetric distributions
           out <- resu^2/sigma2
           # Distribution-specific likelihood kernels
           out <- switch(dist,
                         "Gaussian"={exp(-out/2)},
                         "Laplace"={sqrt(pi/2)*exp(-sqrt(out))},
                         "Student-t"={(1 + out/nu)^(-(nu+1)/2)*gamma((nu+1)/2)/(sqrt(nu/2)*gamma(nu/2))},
                         "Slash"={ifelse(out==0,nu/(nu+1),gamma((1+nu)/2)*(nu/2)*pgamma(1,shape=(1+nu)/2,rate=out/2)/((out/2)^((1+nu)/2)))})
           # Normalizing constant
           out <- out/sqrt(2*pi*sigma2)
      }
      out <- out*(prod(1-rho*w1ev)*prod(1-lambda*w2ev))^(1/nrow(y))
      return(out)
  }
  # Collect fitted mtar objects
  another <- list(...)
  call. <- match.call()
  # Initialize output container
  out <- matrix(0,length(another),1)
  outnames <- vector()
  # Loop over each fitted model
  for(l in 1:(length(another))){
      n.sim <- another[[l]]$n.sim
      dist <- another[[l]]$dist
	  y <- another[[l]]$data$y
	  X <- another[[l]]$data$X
	  Z <- another[[l]]$data$Z
	  W1y <- another[[l]]$misc$W1y
	  W2y <- another[[l]]$misc$W2y
	  W2W1y <- another[[l]]$misc$W2W1y
	  W2X <- another[[l]]$misc$W2X
	  w1ev <- another[[l]]$misc$w1ev
	  w2ev <- another[[l]]$misc$w2ev
	  a <- matrix(0,nrow(y),1)
	  b <- 0
      # Loop over MCMC iterations
      for(j in 1:n.sim){
	      betaj <- another[[l]]$chains$beta[,j]
		  gammaj <- another[[l]]$chains$gamma[,j]
		  if(is.null(another[[l]]$chains$rho)) rhoj <- 0 else rhoj <- another[[l]]$chains$rho[,j]
		  if(is.null(another[[l]]$chains$lambda)) lambdaj <- 0 else lambdaj <- another[[l]]$chains$lambda[,j]
          if(is.null(another[[l]]$chains$delta)) deltaj <- 0 else deltaj <- another[[l]]$chains$delta[,j]
          if(is.null(another[[l]]$chains$nu)) nuj <- 10000 else nuj <- another[[l]]$chains$nu[,j]
		  tempj <- matrix(Lik(dist,betaj,gammaj,rhoj,lambdaj,deltaj,nuj),nrow(y),1)
          a <- a + tempj
		  b <- b + sum(log(tempj))
      }
	  a <- sum(log(a/n.sim))
	  b <- b/n.sim
      # Watanabe-Akaike or Widely Applicable Information Criterion
      out[l] <- -2*a + 4*(a-b)
      outnames[l] <- as.character(call.[l+1])
  }
  # Assign names and return WAIC values
  rownames(out) <- outnames
  colnames(out) <- "WAIC"
  return(out)
}

#' Quantile-Type Residuals for a Spatial Regression Model
#'
#' Calculates quantile-type residuals from a fitted spatial regression model.
#' Posterior means are used for the model parameters. For non-Gaussian error
#' distributions, the residuals are transformed using the corresponding
#' distribution function; simulation is used in the skewed-distribution cases.
#'
#' @param object An object of class `"spreg"`, typically returned by
#'   [spatialReg()].
#' @param ... Additional arguments (currently unused).
#'
#' @return A one-column matrix of quantile-type residuals with class
#'   `"residuals_spreg"`.
#'
#' @seealso [spatialReg()]
#' @method residuals spreg
#' @export
residuals.spreg <- function(object,...){
	dist <- object$dist
	y <- object$data$y
	X <- object$data$X
	Z <- object$data$Z
	W1y <- object$misc$W1y
	W2y <- object$misc$W2y
	W2W1y <- object$misc$W2W1y
	W2X <- object$misc$W2X
	beta <- rowMeans(object$chains$beta)
	gamma <- rowMeans(object$chains$gamma)
	if(is.null(object$chains$rho)) rho <- 0 else rho <- mean(object$chains$rho)
	if(is.null(object$chains$lambda)) lambda <- 0 else lambda <- mean(object$chains$lambda)
	if(!is.null(object$chains$delta)) delta <- mean(object$chains$delta)
	if(!is.null(object$chains$nu)) nu <- mean(object$chains$nu)
	# Compute residuals
	resu <- y - rho*W1y - lambda*W2y + rho*lambda*W2W1y - (X%*%beta - lambda*W2X%*%beta)
	stdev <- exp(0.5*Z%*%gamma)
	if(dist=="Gaussian") resu <- resu/stdev
	if(dist=="Student-t") resu <- qnorm(pmin(pmax(pt(resu/stdev,df=nu),.Machine$double.eps),1-.Machine$double.eps))
	if(dist=="Laplace") resu <- qnorm(pmin(pmax(ifelse(resu<0,0.5*exp(resu/stdev),1-0.5*exp(-resu/stdev)),.Machine$double.eps),1-.Machine$double.eps))
	T <- 10000
	if(dist=="Skew-normal"){
      E <- matrix(rnorm(T*nrow(y)),nrow=T,ncol=nrow(y))
	  EPS <- delta*qnorm(0.5 + 0.5*runif(T)) + sweep(E,2,stdev,"*")
      resu <- colMeans(sweep(EPS,2,resu,"<="))
	  resu <- matrix(qnorm(pmin(pmax(resu,.Machine$double.eps),1-.Machine$double.eps)),nrow(y),1)
	}
	if(dist=="Skew-Student-t"){
      E <- matrix(rnorm(T*nrow(y)),nrow=T,ncol=nrow(y))
	  EPS <- (delta*qnorm(0.5 + 0.5*runif(T)) + sweep(E,2,stdev,"*"))/sqrt(matrix(rgamma(T*nrow(y),shape=nu/2,rate=nu/2),T,nrow(y)))
      resu <- colMeans(sweep(EPS,2,resu,"<="))
	  resu <- matrix(qnorm(pmin(pmax(resu,.Machine$double.eps),1-.Machine$double.eps)),nrow(y),1)
	}
	row.names(resu) <- NULL
	colnames(resu) <- ""
	class(resu) <- "residuals_spreg"
	return(resu)
}
#' Deviance Information Criterion
#'
#' Computes the Deviance Information Criterion.
#'
#' @param ... Arguments passed to methods.
#' @return The result depends on the method used. For objects of
#'   class `spreg`, `DIC.spreg()` returns a one-column numeric
#'   matrix containing the DIC values for the fitted models.
#' @export
DIC <- function(...) {
  UseMethod("DIC")
}

#' Watanabe-Akaike Information Criterion
#'
#' Computes the Watanabe-Akaike Information Criterion.
#'
#' @param ... Arguments passed to methods.
#' @return The result depends on the method used. For objects of
#'   class `spreg`, `WAIC.spreg()` returns a one-column numeric
#'   matrix containing the DIC values for the fitted models.
#' @export
WAIC <- function(...) {
  UseMethod("WAIC")
}


#' Set and Validate Prior Hyperparameters
#'
#' Sets default values and validates the hyperparameters of the prior
#' distributions used in Bayesian spatial regression models. The
#' hyperparameters are configured according to the model specification,
#' including the error distribution and the presence of spatial dependence
#' in the response variable and in the error terms.
#'
#' The function assigns default values to unspecified hyperparameters and
#' checks that the supplied values satisfy the required numerical
#' constraints. For regression location coefficients, the prior is
#' parameterized by its mean and variance. For scale coefficients, the
#' prior depends on whether heteroscedasticity is included in the model.
#' Additional hyperparameters are handled when spatial autoregressive
#' terms, skewness, or Student-t errors are present.
#'
#' @param prior A named list containing prior hyperparameters. Missing
#'   hyperparameters are assigned default values when applicable.
#' @param dist A character string specifying the error distribution.
#'   Supported values include \code{"Gaussian"}, \code{"Student-t"},
#'   \code{"Laplace"}, \code{"Skew-normal"}, and
#'   \code{"Skew-Student-t"}.
#' @param hete A logical value indicating whether the model includes
#'   covariate-dependent error variance.
#' @param em A logical value indicating whether spatial dependence is
#'   included in the error terms through the parameter \eqn{\lambda}.
#' @param ar A logical value indicating whether spatial dependence is
#'   included in the response variable through the parameter \eqn{\rho}.
#'
#' @details
#' The following prior hyperparameters are supported:
#' \describe{
#'   \item{\code{mubeta}, \code{sigma2beta}}{
#'     Mean and variance of the Gaussian prior for the location
#'     regression coefficients. Defaults are 0 and \eqn{10^6},
#'     respectively.
#'   }
#'   \item{\code{mugamma}, \code{sigma2gamma}}{
#'     Mean and variance of the Gaussian prior for the coefficients
#'     governing the log error variance, used when \code{hete = TRUE}.
#'     Defaults are 0 and \eqn{10^6}, respectively.
#'   }
#'   \item{\code{shape}, \code{rate}}{
#'     Shape and rate parameters of the Gamma prior for the error
#'     variance when \code{hete = FALSE}. Both default to
#'     \eqn{10^{-6}}.
#'   }
#'   \item{\code{murho}, \code{sigma2rho}}{
#'     Mean and variance of the Gaussian prior for the spatial
#'     autoregressive parameter \eqn{\rho}, used when \code{ar = TRUE}.
#'     Defaults are 0 and \eqn{10^6}, respectively.
#'   }
#'   \item{\code{mulambda}, \code{sigma2lambda}}{
#'     Mean and variance of the Gaussian prior for the spatial error
#'     parameter \eqn{\lambda}, used when \code{em = TRUE}.
#'     Defaults are 0 and \eqn{10^6}, respectively.
#'   }
#'   \item{\code{mudelta}, \code{sigma2delta}}{
#'     Mean and variance of the Gaussian prior for the skewness
#'     parameter \eqn{\delta}, used for the skew-normal and
#'     skew-Student-t distributions. Defaults are 0 and \eqn{10^6},
#'     respectively.
#'   }
#'   \item{\code{anu}, \code{bnu}}{
#'     Shape and rate parameters of the Gamma prior for the degrees
#'     of freedom parameter \eqn{\nu}, used for Student-t and
#'     skew-Student-t errors. The default for \code{anu} is 1.
#'     The value of \code{bnu} must be supplied and be strictly positive.
#'     The value of \code{anu} must be greater than or equal to 1.
#'   }
#' }
#'
#' All Gaussian prior variances must be finite and strictly positive.
#' Prior means must be finite numeric scalars. The Gamma prior parameters
#' must be finite positive numeric scalars, with the additional
#' constraint that \code{anu} must be at least 1.
#'
#' @return A named list containing the validated prior hyperparameters,
#'   including default values assigned to unspecified components.
#'
#' @seealso \code{\link{spatialReg}}
#'
#'
priors <- function(prior,dist,hete,em,ar){
  # Set default prior mean for location parameter
  if(is.null(prior$mubeta)) prior$mubeta <- 0
  if(!is.numeric(prior$mubeta) || length(prior$mubeta)!=1L ||
      is.na(prior$mubeta) || !is.finite(prior$mubeta)) stop("The hyperparameter 'mubeta' must be a finite numeric value.",call.=FALSE)
  # Set default prior variance for location parameter
  if(is.null(prior$sigma2beta)) prior$sigma2beta <- 1e6
  if(!is.numeric(prior$sigma2beta) || length(prior$sigma2beta)!=1L || is.na(prior$sigma2beta)
     || !is.finite(prior$sigma2beta) || prior$sigma2beta <= 0) stop("The hyperparameter 'sigma2beta' must be a finite numeric value greater than 0",call.=FALSE)
  # Set default prior mean and variance for scale parameter
  if(hete){
    if(is.null(prior$mugamma)) prior$mugamma <- 0
    if(!is.numeric(prior$mugamma) || length(prior$mugamma)!=1L ||
       is.na(prior$mugamma) || !is.finite(prior$mugamma)) stop("The hyperparameter 'mugamma' must be a finite numeric value.",call.=FALSE)
    if(is.null(prior$sigma2gamma)) prior$sigma2gamma <- 1e6
    if(!is.numeric(prior$sigma2gamma) || length(prior$sigma2gamma)!=1L || is.na(prior$sigma2gamma)
       || !is.finite(prior$sigma2gamma) || prior$sigma2gamma <= 0) stop("The hyperparameter 'sigma2gamma' must be a finite numeric value greater than 0",call.=FALSE)
  }else{
    if(is.null(prior$shape)) prior$shape <- 1/1e6
    if(!is.numeric(prior$shape) || length(prior$shape)!=1L || is.na(prior$shape)
       || !is.finite(prior$shape) || prior$shape <= 0)
      stop("The hyperparameter 'shape' must be a finite numeric value greater than 0.",call.=FALSE)
    if(is.null(prior$rate)) prior$rate <- 1/1e6
    if(!is.numeric(prior$rate) || length(prior$rate)!=1L || is.na(prior$rate)
       || !is.finite(prior$rate) || prior$rate <= 0)
      stop("The hyperparameter 'rate' must be a finite numeric value greater than 0.",call.=FALSE)
  }
  if(ar){
    # Set default prior mean for rho
    if(is.null(prior$murho)) prior$murho <- 0
    if(!is.numeric(prior$murho) || length(prior$murho)!=1L ||
      is.na(prior$murho) || !is.finite(prior$murho)) stop("The hyperparameter 'murho' must be a numeric value",call.=FALSE)
    # Set default prior variance for rho
    if(is.null(prior$sigma2rho)) prior$sigma2rho <- 1e6
    if(!is.numeric(prior$sigma2rho) || length(prior$sigma2rho)!=1L || is.na(prior$sigma2rho)
      || !is.finite(prior$sigma2rho) || prior$sigma2rho <= 0) stop("The hyperparameter 'sigma2rho' must be a finite numeric value greater than 0.",call.=FALSE)
  }
  if(em){
    # Set default prior mean for lambda
    if(is.null(prior$mulambda)) prior$mulambda <- 0
    if(!is.numeric(prior$mulambda) || length(prior$mulambda)!=1L ||
      is.na(prior$mulambda) || !is.finite(prior$mulambda)) stop("The hyperparameter 'mulambda' must be a numeric value",call.=FALSE)
    # Set default prior variance for lambda
    if(is.null(prior$sigma2lambda)) prior$sigma2lambda <- 1e6
    if(!is.numeric(prior$sigma2lambda) || length(prior$sigma2lambda)!=1L || is.na(prior$sigma2lambda)
      || !is.finite(prior$sigma2lambda) || prior$sigma2lambda <= 0) stop("The hyperparameter 'sigma2lambda' must be a finite numeric value greater than 0.",call.=FALSE)
  }
  if(dist %in% c("Skew-normal","Skew-Student-t")){
    # Set default prior mean for delta
    if(is.null(prior$mudelta)) prior$mudelta <- 0
    if(!is.numeric(prior$mudelta) || length(prior$mudelta)!=1L ||
       is.na(prior$mudelta) || !is.finite(prior$mudelta)) stop("The hyperparameter 'mudelta' must be a numeric value",call.=FALSE)
    # Set default prior variance for delta
    if(is.null(prior$sigma2delta)) prior$sigma2delta <- 1e6
    if(!is.numeric(prior$sigma2delta) || length(prior$sigma2delta)!=1L || is.na(prior$sigma2delta)
       || !is.finite(prior$sigma2delta) || prior$sigma2delta <= 0) stop("The hyperparameter 'sigma2delta' must be a finite numeric value greater than 0.",call.=FALSE)
  }
  if(dist %in% c("Student-t","Skew-Student-t")){
    # Set default prior for nu
    if(is.null(prior$anu)) prior$anu <- 1
    if(!is.numeric(prior$anu) || length(prior$anu)!=1L || is.na(prior$anu)
       || !is.finite(prior$anu) || prior$anu < 1)
       stop("The hyperparameter 'anu' must be a finite numeric value greater than or equal to 1.",call.=FALSE)
    if(is.null(prior$bnu)) prior$bnu <- 1e-6
    if(!is.numeric(prior$bnu) || length(prior$bnu)!=1L || is.na(prior$bnu)
       || !is.finite(prior$bnu) || prior$bnu <= 0)
       stop("The hyperparameter 'bnu' must be a finite numeric value greater than 0.",call.=FALSE)
  }
  # Return the validated and completed prior list
  return(prior)
}
#' Boston Census Tracts
#'
#' An \code{sf} object containing polygon geometries and associated
#' attributes for census tracts in Boston, Massachusetts.
#'
#' @format An \code{sf} object containing the census tract polygons
#'   and their associated attributes.
#'
#' @details
#' The data are read from the GeoPackage file
#' \code{boston_tracts.gpkg} distributed with the \pkg{spData} package.
#' The original geometries are used to construct the contiguity-based
#' spatial weight matrix \code{W1}. A projected version of the
#' geometries, using EPSG:26986, is used to calculate centroids for
#' constructing \code{W2}.
#'
#' @source The \pkg{spData} package.
#'
#' @seealso \code{\link{W1}}, \code{\link{W2}}
#' @keywords datasets
"boston"


#' Boston Contiguity Spatial Weight Matrix
#'
#' A row-standardized spatial weight matrix based on polygon contiguity
#' among Boston census tracts.
#'
#' @format A numeric matrix with one row and one column for each
#'   census tract in \code{\link{boston}}.
#'
#' @details
#' The matrix is constructed using \code{spdep::poly2nb()} to identify
#' polygon neighbors and \code{spdep::nb2mat()} with
#' \code{style = "W"} to obtain row-standardized weights.
#'
#' @seealso \code{\link{boston}}, \code{\link{W2}}
#' @keywords datasets
"W1"


#' Boston Four-Nearest-Neighbor Spatial Weight Matrix
#'
#' A row-standardized spatial weight matrix based on symmetrized
#' four-nearest-neighbor relationships among Boston census tract centroids.
#'
#' @format A numeric matrix with one row and one column for each
#'   census tract in \code{\link{boston}}.
#'
#' @details
#' The polygons are projected to EPSG:26986 before their centroids
#' are calculated. The four nearest neighbors are identified using
#' \code{spdep::knearneigh()}, converted to a neighbor list using
#' \code{spdep::knn2nb()}, and symmetrized using
#' \code{spdep::make.sym.nb()}. The matrix is then constructed using
#' \code{spdep::nb2mat()} with \code{style = "W"}.
#'
#' @seealso \code{\link{boston}}, \code{\link{W1}}
#' @keywords datasets
"W2"
