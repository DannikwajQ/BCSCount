#' Estimate the coefficients of the ZIP model

#' Takes in response variable,prognostic and predictive covariates of both Poisson and Zero component
#' @import credsubs
#' @import rjags

#' @param Prog.z prognostic covariate of the logistic regression (zero component)
#' @param Pred.z predictive covariate of logistic regression (zero component)
#' @param Prog.p prognostic covariate of the Poisson regression (Poisson component)
#' @param Pred.p predictive covariate of Poisson regression (Poisson component)
#' @param Y count response variable
#' @param Trt treatment covariate
#' @param offset.variable offset variable
#' @return coefficients of the Logistic and Poisson regression model
#' @export
#' @example \dontrun{
#' sun=read.csv('sunwei.csv')


#'#sun= subset(sun, sun$end== 1)

#'#y=sun$count;n=nrow(sun)
#'off.set=matrix(log(sun$month),nrow(sun),1)
#'#prog.p=cbind(sun[,c(3,5,4)])
#'#pred.p=sun[,c(5,4)]
#'#prog.z=sun[,c(3,5,4)]
#'#pred.z=sun[c(5)]

#pred.z=matrix(0,nrow(sun),1)
#'#trt=sun$treat

##########################################################
#'#nIter=15000;zeta.diff = 0;cred.level=0.8;off.set=sun$month
###########################################################

#'#JagResults=zipreg(y,trt,prog.z,prog.p,pred.p,nIter,off.set)
#'
#' }

zipreg=function (y,trt,prog.z,prog.p,pred.p,nIter,offset.variable){

  TRT=diag(trt)
  prog.z=as.matrix(prog.z)
  prog.p=as.matrix(prog.p)
  #pred.z=as.matrix(pred.z)
  pred.p=as.matrix(pred.p)

  interac.x=matrix(nrow(pred.p),ncol(pred.p))
  #interac.z=matrix(nrow(pred.z),ncol(pred.z))

  #interaction terms= product of trt and  predictive covariate
  interac.x = (TRT%*%pred.p)

  #interac.z = (TRT%*%pred.z)


  offset.variable=offset.variable

  X=cbind(1,prog.p,interac.x)
  Z=cbind(1,prog.z)

  params = c('alpha','beta','indx','indz')


  datalist1 <- list(X=X, Z=Z,y=y,n=nrow(pred.p),offset.variable=offset.variable,nKx=ncol(X),nKz=ncol(Z),Kp=2,Kz=2)

  ###JAGS for zipmodel1
  jagsmodel1 <- jags.model(textConnection(zipModel),data=datalist1,quiet=TRUE)
  ###Posterior Sample for zipmodel1
  posteriorSample<- coda.samples(jagsmodel1, params , n.thin=20, n.chains=3,n.burnin=500, n.iter=nIter, progress.bar='none')
  posteriorSample=posteriorSample[[1]]


  return(list(pos.logit = as.matrix(posteriorSample[,1:ncol(Z)]),
              pos.log =as.matrix(posteriorSample[,(ncol(Z)+1):(ncol(Z)+ncol(X))]),
              posteriorSummary=round(summary(posteriorSample)[[1]][,1],4)))



}





#' Calculates expected value of the outcome for both trt and ctrl seperately
#'
#' Takes in prognostic and predictive covariates and their corresponding coefficients of both Poisson and Zero components

#' @param design.trt design treatment =c(0,1)
#' @param design.prog.z prognostic design matrix of the logistic regression (zero component)
#' @param design.pred.z predictive design matrix of logistic regression (zero component)
#' @param design.prog.p prognostic design matrix of the Poisson regression (Poisson component)
#' @param design.pred.p predictive design matrix of Poisson regression (Poisson component)
#' @param coef.prog.z coefficients of prognostic design matrix of the logistic regression (zero component)
#' @param coef.pred.z coefficients of predictive covariate of logistic regression (zero component)
#' @param coef.prog.p coefficients of prognostic covariate of the Poisson regression (Poisson component)
#' @param coef.pred.p coefficients of predictive covariate of Poisson regression (Poisson component)
#' @return expected value of trt and ctrl
#' @export

CalRegression = function(design.trt=c(0,1),design.prog.z,design.prog.p,design.pred.p,coef.prog.z, coef.prog.p,coef.pred.p){


  design.trt <- design.trt[1]
  stopifnot(design.trt %in% c(0,1))

  if(design.trt==1){

    #trt=as.vector(rep(design.trt,nrow(design.prog.p)))

    TRT=diag(nrow(design.prog.p))

    #d.trt=  as.vector(rep(design.trt,nrow(design.prog.p)))

    prog.z=cbind(1,rep(design.trt,nrow(design.prog.p)),design.prog.z)
    prog.p=cbind(1,rep(design.trt,nrow(design.prog.p)),design.prog.p)

    #pred.z=as.matrix(design.pred.z)
    pred.p=as.matrix(design.pred.p)

    coef.prog.z=as.matrix(coef.prog.z)
    coef.prog.p=as.matrix(coef.prog.p)
    #coef.pred.z=as.matrix(coef.pred.z)
    coef.pred.p=as.matrix(coef.pred.p)


    interac.x=matrix(nrow(pred.p),ncol(pred.p))
    #interac.z=matrix(nrow(pred.z),ncol(pred.z))

    #interaction terms= product of trt and  predictive covariate
    interac.x = (TRT%*%pred.p)

    #interac.z = (TRT%*%pred.z)

    X.design=cbind(prog.p,interac.x)
    Z.design=cbind(prog.z)

    coef.pois=  cbind(coef.prog.p,coef.pred.p)
    coef.zero=  cbind(coef.prog.z)


    #design=designCtrl
    nIter = nrow(coef.pred.p)# number of iterations
    nsubj = nrow(design.prog.z)# different characteristics of patient

    ##empty matrices
    e.zero = matrix(0, nIter, nsubj)#(1-theta) = probability rv is poisson
    e.pois = matrix(0, nIter, nsubj)#(mu) mean of the poisson.

    for(isubj in 1:nsubj){
      e.zero[,isubj] = as.vector(1/(1+exp(Z.design[isubj,]%*%t(coef.zero))))
      e.pois[,isubj] = as.vector(exp(X.design[isubj,]%*%t(coef.pois)))
    }




    return(list(e.zero=e.zero, e.pois=e.pois))
  }
  else{
    trt=as.vector(rep(design.trt,nrow(design.prog.p)))
    TRT=diag(trt)
    d.trt=  as.vector(rep(design.trt,nrow(design.prog.p)))

    prog.z=as.matrix(1,d.trt,design.prog.z)
    prog.p=as.matrix(1,d.trt,design.prog.p)

    #pred.z=as.matrix(design.pred.z)
    pred.p=as.matrix(design.pred.p)

    coef.prog.z=as.matrix(coef.prog.z)
    coef.prog.p=as.matrix(coef.prog.p)
    #coef.pred.z=as.matrix(coef.pred.z)
    coef.pred.p=as.matrix(coef.pred.p)


    interac.x=matrix(nrow(pred.p),ncol(pred.p))
    #interac.z=matrix(nrow(pred.z),ncol(pred.z))

    #interaction terms= product of trt and  predictive covariate
    interac.x = (TRT%*%pred.p)

    #interac.z = (TRT%*%pred.z)

    X.design=cbind(1,d.trt,design.prog.p,interac.x)
    Z.design=cbind(1,d.trt,design.prog.z)

    coef.pois=  matrix(cbind(coef.prog.p,coef.pred.p),nIter)
    coef.zero=  matrix(cbind(coef.prog.z),nIter)


    #design=designCtrl
    nIter = nrow(coef.pred.p)# number of iterations
    nsubj = nrow(design.prog.z)# different characteristics of patient

    ##empty matrices
    e.zero = matrix(0, nIter, nsubj)#(1-theta) = probability rv is poisson
    e.pois = matrix(0, nIter, nsubj)#(mu) mean of the poisson.

    for(isubj in 1:nsubj){
      e.zero[,isubj] = as.vector(1/(1+exp(Z.design[isubj,]%*%t(coef.zero))))
      e.pois[,isubj] = as.vector(exp(X.design[isubj,]%*%t(coef.pois)))
    }

    return(list(e.zero=e.zero, e.pois=e.pois))
  }


}

#'Estimate the personalized treatment effect
#'
#' Takes the expected value of treatment and control
#'@param e.trt expected value for the outcome of the treatment
#'@param e.ctrl expected value for the outcome of the control
#'@return personalized treatment effect(PTE)
#'@export

zippte =function(e.trt,e.ctrl){

  nsubj = ncol(e.ctrl$e.pois )
  nIter = nrow(e.ctrl$e.pois)

  safety.diff  = matrix(0, nIter, nsubj)

  for(isubj in 1:nsubj){
    safety.diff[,isubj] = ((e.trt$e.zero[,isubj])*(e.trt$e.pois[,isubj]))-
      ((e.ctrl$e.zero[,isubj])*(e.ctrl$e.pois[,isubj]))
  }
  return(safety.diff)

}


#' specifies the prior and likelihood
#'
zipModel <- "model{
## Likelihood
    for(i in 1:n){
        y[i] ~ dpois(mu[i])
        mu[i] <- lambda[i]*(1-zero[i]) + 1e-10*zero[i]

        lambda[i] <- exp(mu.count[i]+log(offset.variable[i]))
        mu.count[i] <- inprod(beta[],X[i,])

        ## Zero-Inflation
        zero[i] ~ dbern(pi[i])
        pi[i] <- ilogit(mu.binary[i])
        mu.binary[i] <- inprod(alpha[], Z[i,])
    }

    ## Priors

    for(k.p in 1:Kp){
        beta[k.p] ~ dnorm(0, 0.01)
    }

    for(k.z in 1:Kz){
        alpha[k.z] ~ dnorm(0, 0.01)
    }

    for(k in (Kp+1):nKx){

        indx[k] ~ dbern(0.5)
        zz[k] <- equals(indx[k], 0)
        F[k] <- zz[k]*0.00001 + 1 - zz[k]
        eta[k] ~ dgamma(0.0001,0.0001)
        gamma[k] <- F[k]*(1/eta[k])
        beta[k] ~ dnorm(0, 1/gamma[k])
    }

    for(j in (Kz+1):nKz){

        indz[j] ~ dbern(0.5)
        zza[j] <- equals(indz[j], 0)
        Fa[j] <- zza[j]*0.00001 + 1 - zza[j]
        etaa[j] ~ dgamma(0.0001,0.0001)
        gammaa[j] <- Fa[j]*(1/etaa[j])
        alpha[j] ~ dnorm(0, 1/gammaa[j])
    }

}"

##################################################
#'credible.bands calculates the credible bands
#'

credible.bands<-function (params, design = NULL, FUN = function(x, params) {
  params %*% t(x)
}, cred.level = 0.95, method = c("asymptotic", "quantile"),
sides = c("both", "upper", "lower"), est.FUN = mean,
var.FUN = sd, point.estimate = NULL, track = numeric(0),
verbose = FALSE)
{
  params <- data.matrix(params)
  M <- nrow(params)
  if (is.null(design)) {
    N <- ncol(params)
    nonpar <- TRUE
  }
  else {
    design <- data.matrix(design)
    N <- nrow(design)
    nonpar <- FALSE
  }
  est <- var <- numeric(N)
  method <- method[1]
  stopifnot(method %in% c("asymptotic", "quantile"))
  if (method == "asymptotic") {
    m <- numeric(N)
    s <- numeric(N)
  }
  sides <- sides[1]
  stopifnot(sides %in% c("both", "upper", "lower"))
  sim.cred.band <- list(upper = rep(Inf, N), lower = rep(-Inf,
                                                         N), cred.level = cred.level, method = method, sides = sides,
                        W.crit = NA, W = NA, est = rep(NA, N), est.FUN = est.FUN,
                        var = rep(NA, N), var.FUN = var.FUN, trace = matrix(NA,
                                                                            nrow = M, ncol = length(track)))
  colnames(sim.cred.band$trace) <- track
  class(sim.cred.band) <- "sim.cred.band"
  W <- rep(-Inf, M)
  for (i in 1:N) {
    if (verbose && (i%%100 == 0)) {
      cat(i, "/", N, "\n")
    }
    if (nonpar) {
      fx <- params[, i]
    }
    else {
      fx <- FUN(design[i, , drop = FALSE], params)
    }
    if (i %in% track) {
      sim.cred.band$trace[, which(track == i)] <- fx
    }
    est[i] <- est.FUN(fx)
    var[i] <- var.FUN(fx)
    if (method == "asymptotic") {
      if (is.null(point.estimate)) {
        m[i] <- mean(fx)
        s[i] <- sd(fx)
      }
      else {
        if (nonpar) {
          point.fx <- point.estimate[i]
        }
        else {
          point.fx <- FUN(design[i, , drop = FALSE],
                          point.estimate)
        }
        m[i] <- point.fx
        s[i] <- sqrt(mean((fx - point.fx)^2))
      }
      if (sides == "both") {
        z <- abs(fx - m[i])/s[i]
      }
      else if (sides == "upper") {
        z <- (fx - m[i])/s[i]
      }
      else if (sides == "lower") {
        z <- (m[i] - fx)/s[i]
      }
    }
    else {
      if (sides == "both") {
        Fx <- to.Fx(fx)
        Gx <- 1 - to.Fx(-fx)
        z <- pmax(1 - Fx, Gx)
      }
      else if (sides == "upper") {
        z <- 1 - to.Fx(-fx)
      }
      else if (sides == "lower") {
        z <- 1 - to.Fx(fx)
      }
    }
    W <- pmax(W, z)
  }
  sim.cred.band$W <- W
  W.crit <- quantile(W, cred.level, type = 1)
  sim.cred.band$W.crit <- W.crit
  if (method == "asymptotic") {
    if (sides == "both" | sides == "upper") {
      sim.cred.band$upper <- m + W.crit * s
    }
    if (sides == "both" | sides == "lower") {
      sim.cred.band$lower <- m - W.crit * s
    }
  }
  else {
    if (nonpar) {
      if (sides %in% c("both", "upper")) {
        sim.cred.band$upper <- apply(params, 2, function(fx) {
          -quantile(-fx, prob = 1 - W.crit, type = 1)
        })
      }
      if (sides %in% c("both", "lower")) {
        sim.cred.band$lower <- apply(params, 2, function(fx) {
          quantile(fx, prob = 1 - W.crit, type = 1)
        })
      }
    }
    else {
      if (sides == "both") {
        bounds <- apply(design, 1, function(x, params) {
          fx <- FUN(t(x), params)
          upper <- -quantile(-fx, prob = 1 - W.crit,
                             type = 1)
          lower <- quantile(fx, prob = 1 - W.crit, type = 1)
          c(lower, upper)
        }, params = params)
        sim.cred.band$lower <- bounds[1, ]
        sim.cred.band$upper <- bounds[2, ]
      }
      else if (sides == "upper") {
        sim.cred.band$upper <- apply(design, 1, function(x,
                                                         params) {
          -quantile(-fx, prob = 1 - W.crit, type = 1)
        }, params = params)
      }
      else if (sides == "lower") {
        sim.cred.band$lower <- apply(design, 1, function(x,
                                                         params) {
          quantile(fx, prob = 1 - W.crit, type = 1)
        }, params = params)
      }
    }
  }
  sim.cred.band$est <- est
  sim.cred.band$var <- var
  sim.cred.band
}




####################################################
#'Constructs a credible subset pair
#'
#'cred.subgroup returns a credible subset pair over a finite set of covariate points given either a sample from the posterior of the regression surface or a function FUN(x, params) and a sample from the posterior of the parameters.

#'@param param A numeric matrix whose rows are draws from the posterior distribution of either the regression surface or the parameter vector.
#'@param design	(Optional) A numeric matrix whose rows are covariate points over which the band is to be constructed.
#'@param FUN(Optional) a function of the form function(x, params) that takes a row of design and the entire params matrix and returns a vector of the same length of x representing the regression surface.
#'@param cred.level	 the credible level.
#'@param threshold	 the value of FUN above which a covariate is included in the target subset.
#'@param method	 Either "asymptotic" (default) or "quantile"; see details.
#'@param step.down Logical (default TRUE); should the step-down procedure be used?
#'@param sides One of "both" (default), "exclusive", or "inclusive". Which bounds should be constructed?
#'@param est.FUN The function used to produce estimates of the regression surface. Default is mean.
#'@param var.FUN The function used to quantify the variability of the regression surface posterior. Default is sd.
#'@param point.estimate	 If not null, replaces the mean and sets the reference around which the standard error is computed. Useful for bootstrapping methods. Treated as a row of the params matrix.
#'@param track A numeric vector of indices indicating which rows (default none) of the design matrix should have the sample of the corresponding FUN(x, params) returned.
#'@param verbose Logical (default FALSE); print progress?
#'@return exclusive A logical vector indicating membership in the exclusive credible subset.
#'@return inclusive A logical vector indicating membership in the inclusive credible subset.
#'@return cred.level As provided.
#'@return threshold As provided.
#'@return method As provided.
#'@return step.down As provided.
#'@return sides As provided.
#'@return est Posterior estimate of the regression surface.
#'@return est.FUN As provided.
#'@return var Summary of posterior variability of the regression surface.
#'@return var.FUN As provided.
#'@return W An estimate of the extremal errors.
#'@return W.crit The critical quantile of W.
#'@return trace The posterior samples of the regression surface indicated by the track argument.
#'@details If design is NULL (default), it is taken to be the identity matrix of dimension ncol(params), so that the rows of params are treated as draws from the posterior FUN(x, params).The 'asymptotic' method assumes that the marginal posteriors of the FUN(x, params) are asymptotically normal and is usually significantly faster and less memory-intensive than the 'quantile' method, which makes no such assumption.
#'@export

cred.subgroup<-function (params, design = NULL, FUN = function(x, params) {
  params %*% t(x)
}, cred.level = 0.95, threshold = 0, method = c("asymptotic",
                                                "quantile"), step.down = TRUE, sides = c("both",
                                                                                         "exclusive", "inclusive"), est.FUN = mean, var.FUN = sd,
point.estimate = NULL, track = numeric(0), verbose = FALSE)
{
  params <- data.matrix(params)
  M <- nrow(params)
  if (is.null(design)) {
    N <- ncol(params)
    nonpar <- TRUE
  }
  else {
    design <- data.matrix(design)
    N <- nrow(design)
    nonpar <- FALSE
  }
  if (verbose) {
    cat("Computing credible subgroups over", N, "points using",
        M, "posterior draws.\n")
  }
  method <- method[1]
  stopifnot(method %in% c("asymptotic", "quantile"))
  if (method == "asymptotic") {
    m <- numeric(N)
    s <- numeric(N)
  }
  sides <- sides[1]
  stopifnot(sides %in% c("both", "exclusive", "inclusive"))
  scb.sides <- ifelse(sides == "both", "both",
                      ifelse(sides == "exclusive", "lower", "upper"))
  credsubs <- list(exclusive = rep(NA, N), inclusive = rep(NA,
                                                           N), cred.level = cred.level, threshold = threshold, method = method,
                   step.down = step.down, sides = sides, est = rep(NA, N),
                   est.FUN = est.FUN, var = rep(NA, N), var.FUN = var.FUN,
                   W = rep(NA, M), W.crit = NA, trace = matrix(NA, nrow = M,
                                                               ncol = length(track)))
  class(credsubs) <- "credsubs"
  credsubs$exclusive <- rep(FALSE, N)
  credsubs$inclusive <- rep(TRUE, N)
  test.set <- 1:N
  reject.set <- numeric(0)
  repeat {
    test.set <- setdiff(test.set, reject.set)
    if (nonpar) {
      test.par <- test.set
    }
    else {
      test.par <- 1:ncol(params)
    }
    if (length(test.set) == N) {
      scb.track <- track
    }
    else {
      scb.track <- numeric(0)
    }
    sim.cred.band <- sim.cred.band(params = params[, test.par,
                                                   drop = FALSE], design = design[test.set, , drop = FALSE],
                                   FUN = FUN, cred.level = cred.level, method = method,
                                   sides = scb.sides, est.FUN = est.FUN, var.FUN = var.FUN,
                                   point.estimate = point.estimate, track = scb.track,
                                   verbose = verbose)
    over.set <- test.set[sim.cred.band$lower > threshold]
    under.set <- test.set[sim.cred.band$upper < threshold]
    credsubs$exclusive[over.set] <- TRUE
    credsubs$inclusive[under.set] <- FALSE
    credsubs$est[test.set] <- sim.cred.band$est
    credsubs$var[test.set] <- sim.cred.band$var
    credsubs$W.crit <- sim.cred.band$W.crit
    credsubs$W <- sim.cred.band$W
    reject.set <- union(over.set, under.set)
    if (length(test.set) == N) {
      credsubs$trace <- sim.cred.band$trace
    }
    if (verbose) {
      cat(length(test.set), "hypotheses tested,",
          length(reject.set), "rejected.\n")
    }
    if (!step.down || length(reject.set) == 0 || length(test.set) ==
        length(reject.set)) {
      break
    }
  }
  credsubs
}






