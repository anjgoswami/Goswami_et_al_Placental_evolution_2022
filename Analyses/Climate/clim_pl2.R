################################################################################
##                                                                            ##
##  RPANDA : Leave-One-Out Cross-Validation for High-dimentional              ##
##           Penalized Likelihood Comparative Methods                         ##
##                                                                            ##
##   Julien Clavel - 01-08-2017                                               ##
##   require: mvMORPH, glassoFast                                             ##
##                                                                            ##
################################################################################


# Y        = traits matrix (columns: variables, rows: species)
# tree     = phylogenetic tree (an object of class 'phylo')
# model    = a function f(x, beta) beta*T(x)
# method   = "RidgeArch": Archetype (linear) Ridge penalty, "RidgeAlt": Quadratic Ridge penalty, "LASSO": Least Absolute Selection and Shrinkage Operator. "RidgeAltapprox" and "LASSOapprox" are fast approximations of the LOOCV for the Ridge quadratic and LASSO penalties
# targM    = "null", "Variance" for a diagonal unequal variance target, "unitVariance" for an equal diagonal target. Only works with "RidgeArch","RidgeAlt", and "RidgeAltapprox" methods.
# REML     = TRUE (default) or FALSE. The likelihood method used to estimate the parameters. REML must be preferred with small sample size in order to obtain unbiased estimates.
# up       = upper bound for model parameter search
# low      = lower bound for model parameter search
# tol      = lower bound tolerance for the regularizer (tuning parameter)
# starting = starting values for the parameter search. Must be a vector likes: c(model, regularization)
# sbound   = if the number of combinations for the initial grid search to find starting values is above sbound, the starting values are re-sampled to sbound

require(mvMORPH)    # >= 1.0.9
require(glassoFast) # https://github.com/JClavel/glassoFast
require(parallel)
require(BB)
#require(optimr)

fit_mvenv <- function(Y, tree, model=NULL, nparam=1, method=c("RidgeAlt","RidgeArch","RidgeAltapprox","LASSO","LASSOapprox"), targM=c("null","Variance","unitVariance"), REML=TRUE, up=NULL, low=NULL, tol=NULL, starting=NULL, SE=NULL, scale.cov=FALSE, opt="L-BFGS-B", sbound=2000, control=list(maxit=5000)){
  
  # Preliminary checks
  if(missing(tree)) stop("Please provide a phylogenetic tree of class \"phylo\" ")
  if(nrow(Y) != Ntip(tree)) stop("Length of phenotypic and phylogenetic data do not match")
  if (all(rownames(Y) %in% tree$tip.label)){
    Y <- Y[tree$tip.label,]
  }else{
    warning("rownames in Y are missing. It is assumed that they are in the same order as in the tree.")
  }
  
  # Select the model
  if(is.null(model)) stop("You must specify a function for the climate model")
  
  # Select the method
  method <- match.arg(method)[1]
  
  # Select the method
  targM <- match.arg(targM)[1]
  
  #if(!inherits(tree,"simmap") & model=="multiEnv") stop("The tree must be in simmap format")
  if(inherits(tree,"simmap")){
    betaSimmap = ncol(tree$mapped.edge)
    tree = reorderSimmap(tree, order="postorder")
    is.simmap <- TRUE
  }else{
    betaSimmap = 1
    tree = reorder(tree, order = "postorder")
    is.simmap <- FALSE
  }
  
  # Bounds for models
  if(is.null(up)){
    up <- rep(50, nparam*betaSimmap)
  }
  
  if(is.null(low)){
    low <- rep(-50, nparam*betaSimmap)
  }
  
  # Reorder the tree
  tree<-reorder.phylo(tree,"postorder")
  
  # Parameters
  if(ncol(Y)==1) stop("Only works with multivariate datasets")
  n <- nO <- nrow(Y)
  nC <- n-1
  if(REML==TRUE) n <- n-1
  
  p <- ncol(Y)
  Yest <- apply(Y,2,function(i) pic(i,tree))
  # Empirical covariance (contrasts) for starting the algorithm
  S <- t(Yest)%*%Yest / n
  
  # Identity matrix
  I <- diag(p)
  
  # Default penalty is Ridge "null"
  target <- matrix(0,p,p)
  
  # Identifying tips values
  tipsIndices <- which(tree$edge[, 2] <= Ntip(tree))
  
  # Default tolerance for the parameter search
  if(is.null(tol)){
    if(method=="RidgeArch"){
      tol = 1e-8
    }else{
      tol = 0
    } 
  }
  
  # Define the climatic function
  phyTrans <- function(phy, mod_par) {
    climTrans(phy, mod_par, model, nO, is.simmap, nparam)
  }
  
  # transformation
  transform <- function(x) x
  
  if(method!="RidgeArch"){
    upperBound <- c(up,log(1e6))
    lowerBound <- c(low,log(tol))
  }else{
    upperBound <- c(up,1)
    lowerBound <- c(low,tol)
  }
  
  # Defines index for parameters
  idx1 <- 1:(nparam*betaSimmap)
  idx2 <- max(idx1)+1
  idx3 <- idx2+1
  
  ## ------ Leave-One-Out Cross-validation / Penalty
  
  switch(method,
         "LASSOapprox"={
           # transform for the tuning parameter
           transAlpha <- function(x) exp(x)
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)             # Compute it in C
             
             # estimate noise
             if(!is.null(SE)){
               error = par[idx3]*par[idx3] 
               tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
             }
             
             Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
             Sk <- crossprod(Yk)/n                     # Compute it in C
             if(any(!is.finite(Sk))) return(1e6)
             
             #  Determinant for the phylogenetic tree
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             if(REML==TRUE){
               Ccov <- sum(log(var_contr))
             }else{
               Ccov <- sum(log(c(var_root,var_contr)))
             }
             
             # LASSO penalty
             LASSO <- glassoFast(Sk, alpha, maxIt=5000)
             G <- LASSO$w
             Gi <- LASSO$wi
             
             # LOO cross-validated log-likelihood
             LOObias <- function(Semp, Iridge, Sridge, Y, lambda){
               # Indices matrice
               Ind <- ifelse(Iridge==0,0,1);
               
               Tk <- sapply(1:nC, function(i){
                 Sk <- Y[i,]%*%t(Y[i,]) ;
                 A <-.vec((Sridge - Sk)*Ind)
                 BC <- Iridge%*%((Semp - Sk)*Ind)%*%Iridge
                 sum(A*BC)
               })
               
               bias <- (1/(2*n*(n-1))) * sum(Tk)
               
               return(bias)
             }
             
             klbias <- LOObias(Sk, Gi, G, Yk, alpha)
             
             ll <- -(1/n)*(-0.5 * (n*p*log(2*pi) + p*Ccov + n*determinant(G)$modulus + n*sum(diag(Gi%*%Sk)))) + klbias
             if (!is.finite(ll)) return(1e6)
             return(ll)
           }
         },
         "LASSO"={
           # transform for the tuning parameter
           transAlpha <- function(x) exp(x)
           require(glassoFast)
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)             # Compute it in C
             
             # estimate noise
             if(!is.null(SE)){
               error = par[idx3]*par[idx3] 
               tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
             }
             
             Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
             Sk <- crossprod(Yk)/n                     # Compute it in C
             if(any(!is.finite(Sk))) return(1e6)
             
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             # log-lik
             llik <- sapply(1:nC, function(x){
               Sk <- crossprod(Yk[-x,])/(n-1)
               LASSO <- glassoFast(Sk, alpha, maxIt=500)
               G <- LASSO$w
               Gi <- LASSO$wi
               
               Swk <- Yk[x,]%*%t(Yk[x,])
               rk <- sum(diag(Swk%*%Gi))
               determinant(G)$modulus + rk
             })
             
             # det of the phylo matrix
             if(REML==TRUE){
               Ccov <- sum(log(var_contr))
             }else{
               Ccov <- sum(log(c(var_root,var_contr)))
             }
             
             ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*mean(llik))
             if (!is.finite(ll)) return(1e6)
             return(ll)
           }
         },
         "RidgeArch"={
           # transform for the tuning parameter
           transAlpha <- function(x) (x)
           
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)
             
             # estimate noise
             if(!is.null(SE)){
               error = par[idx3]*par[idx3] 
               tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
             }
             
             # TRY SCALING??
             if(scale.cov){
               scaling = exp((1/Ntip(tr))*determinant(vcv(tr))$modulus)
               if(scaling!=0) tr$edge.length <- tr$edge.length * (1/scaling)
             }
             
             if (any(is.na(tr$edge.length))) return(1e6)   # if the tree is bad behaved
             Yk <-  apply(Y,2,function(i) pic(i,tr))
             Sk <- crossprod(Yk)/n
             if(any(!is.finite(Sk))) return(1e6)
             
             # rescale the covariance to have unit trace mean?
             #if(scale.cov) Sk <- Sk*(1/scaling)
             # Determinant for the phylogenetic tree
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             if(REML==TRUE){
               Ccov <- sum(log(var_contr))
             }else{
               Ccov <- sum(log(c(var_root,var_contr)))
             }
             
             if(targM=="Variance"){
               target <- diag(diag(Sk))
             }else if(targM=="unitVariance"){
               target <- I*mean(diag(Sk))
             }
             
             # Hoffbeck & Landgrebe 1996 parameterization
             beta <- (1 - alpha)/(n - 1)
             G <- n*beta*Sk + alpha * target
             Gi <- try(chol(G), silent = TRUE)
             if(inherits(Gi, 'try-error')) return(1e6)
             
             # log-lik
             llik <- sapply(1:nC, function(x){
               # log-lik form of Hoffbeck & Landgrebe 1996
               rk <- sum(backsolve(Gi, Yk[x,], transpose = TRUE)^2)
               (n/nC)*log(1 - beta*rk) + (rk/(1 - beta*rk))
             })
             
             ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*sum(2*log(diag(Gi))) + sum(llik))
             
             if (!is.finite(ll)) return(1e6)
             return(ll)
           }
           
         },
         "RidgeAltapprox"={
           # transform for the tuning parameter
           transAlpha <- function(x) exp(x)
           
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)             # Compute it in C
             
             # estimate noise
             if(!is.null(SE)){
               error = par[idx3]*par[idx3] #exp(par[idx3])
               tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
             }
             
             Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
             Sk <- crossprod(Yk)/n                     # Compute it in C
             if(any(!is.finite(Sk))) return(1e6)
             
             # Determinant for the phylogenetic tree
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             if(REML==TRUE){
               Ccov <- sum(log(var_contr))
             }else{
               Ccov <- sum(log(c(var_root,var_contr)))
             }
             
             # Choose an updated target matrix if it's not the usual Ridge
             if(targM=="Variance"){
               target <- diag(1/diag(Sk))
             }else if(targM=="unitVariance"){
               target <- I*(1/mean(diag(Sk)))
             }
             
             # Ridge penalty
             G <- .makePenalty(Sk,alpha,target,targM)$S
             if (any(!is.finite(G))) return(1e6)
             eig <- eigen(G, symmetric=TRUE) # we can return the vectors and inverse from .makePenalty directly
             V <- eig$vectors
             d <- eig$values
             Gi <- V%*%diag(1/d)%*%t(V)
             H <- (1/(kronecker(d,d)+alpha))
             
             
             # LOO cross-validated log-likelihood
             LOObias <- function(Semp, Iridge, Sridge, Y, lambda){
               Tk <- sapply(1:nC, function(i){
                 Sk <- Y[i,]%*%t(Y[i,]) ;
                 VSV <- .vec(t(V)%*%((Sridge - (Sk - lambda*target) - lambda*Iridge))%*%(V))
                 sum(VSV * H*.vec(t(V)%*%(-I%*%(Sk - Semp))%*%(V)))
               })
               
               bias <- (1/(2*n*(n-1))) * sum(Tk)
               
               return(bias)
             }
             
             klbias <- LOObias(Sk, Gi, G, Yk, alpha)
             
             ll <- -(1/n)*(-0.5 * (n*p*log(2*pi) + p*Ccov + n*sum(log(d)) + n*sum(diag(Gi%*%Sk)))) + klbias
             if(!is.finite(ll)) return(1e6)
             return(ll)
           }
         },
         "RidgeAlt"={
           # transform for the tuning parameter
           transAlpha <- function(x) exp(x)
           ## ---- LOOCV function
           loocv <- function(par){
             
             # parameters
             mod_par = transform(par[idx1])
             alpha = transAlpha(par[idx2])
             
             # Transform the tree
             tr <- phyTrans(tree, mod_par)             # Compute it in C
             
             if(!is.null(SE)){
               error = par[idx3]*par[idx3] #exp(par[idx3])
               tr$edge.length[tipsIndices] <- tr$edge.length[tipsIndices] + error
             }
             
             Yk <-  apply(Y,2,function(i) pic(i,tr))   # Compute it in C
             Sk <- crossprod(Yk)/n                     # Compute it in C
             if(any(!is.finite(Sk))) return(1e6)
             
             # Determinant for the phylogenetic tree
             var_pic <- pruning(tr)
             var_root <- var_pic$varRoot
             var_contr <- var_pic$varNode
             
             if(REML==TRUE){
               Ccov <- sum(log(var_contr)) 
             }else{
               Ccov <- sum(log(c(var_root,var_contr))) 
             }
             
             # Choose an updated target matrix if it's not the usual Ridge
             if(targM=="Variance"){
               target <- diag(1/diag(Sk))
             }else if(targM=="unitVariance"){
               target <- I*(1/mean(diag(Sk)))
             }
             
             # log-lik
             llik <- try(mclapply(1:nC, function(x){
               Sk <- crossprod(Yk[-x,])/(n-1)
               pen <- .makePenalty(Sk,alpha,target,targM)
               Gi <- pen$P
               G <- pen$S
               detG <- sum(log(pen$ev))
               Swk <- tcrossprod(Yk[x,]) # Yk[x,]%*%t(Yk[x,])
               rk <- sum(diag(Swk%*%Gi))
               detG + rk
             },mc.cores = getOption("mc.cores", mcores)), silent = TRUE)
             
             if(inherits(llik, 'try-error')) return(1e6)
             # det of the phylo matrix
             ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*mean(unlist(llik)))
             if(!is.finite(ll)) return(1e6)
             return(ll)
           }
           
         }     
  )
  
  # Starting values over range of parameters for quadratic ridge and LASSO (because the tuning value is between 0->Inf)
  # computationally intensive but maybe better to ensure good starting values
  if(is.null(starting)){
    
    # Here we can use the forking to spread the calculus over the grid on several cores
    # we can use a randomized search to speed up the computations of very complex models
    
    if(!is.null(SE)){
      # various errors
      guess <- c(0.001,0.01,0.1,1,10)
      error_guess = sqrt(guess)
      lowerBound = c(lowerBound,0)
      upperBound = c(upperBound,Inf)
    } 
    
    # Initiate search
    message("Initialization via grid search. Please wait...","\n")
    
    if(method=="RidgeArch"){
      range_val <- c(1e-6, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9)
    }else{
      range_val <- log(c(1e-12, 1e-9, 1e-6, 0.01, 0.1, 1, 10, 100, 1000, 10000))
    }
    
    # initial value
    mod_val <- lapply(1:(nparam*betaSimmap), function(x) c(-3,-0.5,0,0.5,3))
    
    if(!is.null(SE)){
      brute_force <- expand.grid(c(mod_val, list(range_val), list(error_guess)))
    }else{
      brute_force <- expand.grid(c(mod_val, list(range_val)))
    }
    
    # prevent the grid search to explode?
    if(nrow(brute_force)>sbound){
        brute_force <- brute_force[round(seq.int(1,nrow(brute_force),length.out = sbound)),]
    }
    
    # assess what are the best starting values
    fit_st <- apply(brute_force,1,loocv)
    start <- brute_force[which.min(fit_st),]
    
    print(start)
    tuning <- start[idx2]
    if(method=="RidgeArch"){
      cat("Best starting for the tuning: ",as.numeric(tuning))
    }else{
      cat("Best starting for the tuning: ",as.numeric(exp(tuning)))
    }
  }else{     
    start <- starting
  }
  
  # Initial guesses found we start the optimization
  message("\n", "Start optimization. Please wait...")
  
  # Optimization of the cross-validated likelihood
  if(opt=="spg"){
    estimModel <- spg(unlist(start), fn = loocv, upper=upperBound, lower=lowerBound, control=control)
  }else{
    estimModel <- optim(start, fn = loocv, method=opt, upper=upperBound, lower=lowerBound, control=control)
  }

  # Compute the scaled tree
  phy_estim <- phyTrans(tree, transform(estimModel$par[idx1]))
  if(any(is.na(phy_estim$edge.length))){
    phy_estim=tree;
    warning("Error in estimation")
  }
  
  # is SE estimated?
  if(!is.null(SE)){
    SE = (estimModel$par[idx3])*(estimModel$par[idx3])
    phy_estim$edge.length[tipsIndices] <- phy_estim$edge.length[tipsIndices] + SE
  }
  
  # TRY SCALING??
  if(scale.cov){
    scaling = exp((1/Ntip(phy_estim))*determinant(vcv(phy_estim))$modulus)
    if(scaling!=0) phy_estim$edge.length <- phy_estim$edge.length * (1/scaling)
  }
  
  # Estimated value for the model parameter
  model.par <- transform(estimModel$par[idx1])
  
  # Estimated value for the regularization parameter
  gamma <- transAlpha(estimModel$par[idx2])
  
  # Compute R
  Ytransform <-  apply(Y,2,function(i) pic(i,phy_estim))
  Snew <- crossprod(Ytransform)/n
  matMeth <- method
  if(method == "LASSOapprox"){
    matMeth <- "LASSO"
  }else if(method == "RidgeAltapprox"){
    matMeth <- "RidgeAlt"
  }
  regularizedEstimates <- .covPenalized(S=Snew, method=matMeth, targM=targM, tuning=gamma)
  
  # End
  message("Done in ", estimModel$count[1]," iterations.")
  
  # return the results
  results <- list(loocv=estimModel$value, model.par=model.par, gamma=gamma, scaled_tree=phy_estim, model="Multivariate PL-Clim", method=method, p=p, n=nO, targM=targM, R=regularizedEstimates, REML=REML, Y=Y, SE=SE, fun=loocv,  opt=estimModel, is.simmap=is.simmap, betaSimmap=betaSimmap, nparam=nparam, tree=tree, clim_fun=model)
  class(results) <- "fit_pl.rpanda"
  return(results)
}

## -------------- Miscellaneous functions
#if(is.function(model)) f<-function(x){ model((mtot+maxdiff)-x, funEnv, param) }
## function to estimate the rates changes
climTrans <- function(phy, beta, fun, nO, is.simmap=FALSE, nparam){
  
  # not ideal, but for debugging now
  # Compute the branching times
  if(is.ultrametric(phy)){
    times<-branching.times(phy)
    names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
  }else{
    # Use "phytools" called by mvMORPH
    times<-max(nodeHeights(phy))-nodeHeights(phy)[match(1:phy$Nnode+Ntip(phy),phy$edge[,1]),1]
    names(times)<-1:phy$Nnode+Ntip(phy)
  }
  
  # root age
  mtot <- max(times)
  # Set the root to zero
  times <- mtot - times
  res <- phy 
  
  if(is.simmap){
    
    # # Transforms the branch-lengths of the tree
    for (i in 1:length(phy$edge.length)) {
      
      age <- times[phy$edge[i, 1] - nO]
      # for simmap
      currentmap<-phy$maps[[i]]           # retrieve the corresponding maps on the edge "i"
      indlength<-length(currentmap)       # How many mapping there are?
      tempedge<-numeric(indlength)        # temporary vector for the mappings
      
      # loop pour traverser les "maps"
      for(betaval in 1:indlength){
        
        regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[betaval])  # retrieve the regimes within maps
        #bet <- beta[regimenumber]           
        bet <- beta[seq_len(nparam) + nparam*(regimenumber-1)]   # select the corresponding parameter for beta - betaSimmap
        bl <- currentmap[[betaval]]           # branch length under the current map
        
        int <- try(integrate(fun, lower=age, upper=age+bl, beta=bet, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
        
        # Try catch if the integrand is divergent
        if(inherits(int ,'try-error')){
          warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
          integ <- NA_real_
        } else {
          integ <- int$value
        }
        
        # tempedge <- tempedge+integ
        tempedge[betaval] <- integ 
        # on met à jour age parcequ'on va passer au maps suivant pour la lignée i
        # update "age" because we're moving to the next map for lineage i.
        age<-age+bl
      }# end of simmap branch loop
      res$edge.length[i] <- sum(tempedge)
    }
    
  }else{
      # Transforms the branch-lengths of the tree
      for (i in 1:length(phy$edge.length)) {
        bl <- phy$edge.length[i]
        age <- times[phy$edge[i, 1] - nO]
        int <- try(integrate(fun, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05, beta=beta), silent = TRUE)
        # Try catch if the integrand is divergent
        if(inherits(int ,'try-error')){
          warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
          integ <- NA_real_
        } else {
          integ <- int$value
        }
        res$edge.length[i] <- integ
      }
  }
  phy<-res
  return(phy)        
}


# Alternative penalty of van Wieringen & Peeters 2016 - Computational Statistics and Data Analysis
# see also Witten & Tibshirani 2009
.makePenalty <- function(S,lambda,target,targM){
  
  switch(targM,
         "Variance"={
           D <- (S - lambda * target)
           Alt <- D/2 + .sqM((D %*% D)/4 + lambda * diag(nrow(S)))
           AltInv <- (1/lambda)*(Alt - D)
           evalues <- eigen(Alt, symmetric=TRUE, only.values = TRUE)$values
         },
         "unitVariance"={
           eig  <- eigen(S, symmetric = TRUE)
           Q <- eig$vectors
           d <- eig$values - lambda*target[1]
           evalues <- sqrt(lambda + d^2/4) + d/2
           D1 <- diag(evalues)
           D2 <- diag(1/evalues) # Inverse
           Alt <- Q %*% D1 %*% t(Q)
           AltInv <- Q %*% D2 %*% t(Q)
         },
         "null"={
           eig  <- eigen(S, symmetric = TRUE)
           Q <- eig$vectors
           d <- eig$values
           evalues <- sqrt(lambda + d^2/4) + d/2
           D1 <- diag(evalues)
           D2 <- diag(1/evalues)
           Alt <- Q %*% D1 %*% t(Q)
           AltInv <- Q %*% D2 %*% t(Q)
         }
  )
  pen <- list(S=Alt, P=AltInv, ev=evalues)
  return(pen)
}

# Matrix square root
.sqM <- function(x){
  if(!all(is.finite(x))) return(Inf)
  eig <- eigen(x, symmetric = TRUE)
  sqrtM <- eig$vectors %*% diag(sqrt(eig$values)) %*% solve(eig$vectors)
  return(sqrtM)
}

# Compute the Regularized covariance and it's inverse
.covPenalized <- function(S, method, targM="null", tuning=0){
  
  # dim of S
  p = ncol(S)
  
  # init the target matrix
  if(method=="RidgeAlt"){
    switch(targM,
           "null"={Target <- matrix(0,p,p)},
           "Variance"={Target <- diag(1/diag(S))},
           "unitVariance"={Target <- diag(1/mean(diag(S)),p)})
  }else if(method=="RidgeArch"){
    switch(targM,
           "null"={Target <- matrix(0,p,p)},
           "Variance"={Target <- diag(diag(S))},
           "unitVariance"={Target <- diag(mean(diag(S)),p)})
  }
  
  # Construct the penalty term
  switch(method,
         "RidgeAlt"={
           pen <- .makePenalty(S,tuning,Target,targM)
           P <- pen$S
           Pi <- pen$P
         },
         "RidgeArch"={
           P <- (1-tuning)*S + tuning*Target
           eig <- eigen(P)
           V <- eig$vectors
           d <- eig$values
           Pi <- V%*%diag(1/d)%*%t(V)
         },
         "LASSO"={
           LASSO <- glassoFast(S,tuning)
           P <- LASSO$w
           Pi <- LASSO$wi
         },
         "ML"={
           P <- S
           eig <- eigen(P)
           V <- eig$vectors
           d <- eig$values
           Pi <- V%*%diag(1/d)%*%t(V)
         }
  )
  estimate <- list(R=P, Rinv=Pi)
  return(estimate)
}


# Vec operator
.vec <- function(x) as.numeric(x)

# Trace operator
.tr <- function(x) sum(diag(x))

# Transform starting values
.starting_val <- function(starting, model, method, SE=NULL){
  switch(model,
         "OU"={
           if(method=="RidgeArch"){
             start <- c(log(starting[1]), starting[2])
           }else{
             start <- c(log(starting[1]), log(starting[2]))
           }
         },
         "EB"={
           if(method=="RidgeArch"){
             start <- c(starting[1], starting[2])
           }else{
             start <- c(starting[1], log(starting[2]))
           }
         },
         "BM"={
           if(method=="RidgeArch"){
             start <- starting[1]
           }else{
             start <- log(starting[1])
           }
         },
         "lambda"={
           if(method=="RidgeArch"){
             start <- c(starting[1], starting[2])
           }else{
             start <- c(starting[1], log(starting[2]))
           }
         })
  
  if(!is.null(SE) & model!="BM") start <- c(start, log(starting[3]))
  if(!is.null(SE) & model=="BM") start <- c(start, log(starting[2]))
  return(start)
}
