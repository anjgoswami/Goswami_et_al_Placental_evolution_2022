################################################################################
##                                                                            ##
##  RPANDA : Generalized Information Criterion (GIC) for high-dimensional     ##
##           Penalized Likelihood Comparative Methods                         ##
##                                                                            ##
##   Julien Clavel - 01-08-2017                                               ##
##   require: mvMORPH, glassoFast                                             ##
##                                                                            ##
################################################################################

require(parallel)

gic_criterion3 <- function(Y, tr, model="BM", method=c("RidgeAlt","RidgeArch","LASSO","ML","RidgeAltapprox","LASSOapprox"), targM=c("null","Variance","unitVariance"), param=NULL, tuning=0, REML=TRUE, nbcores=1L, ...){
  
  # ellipsis for additional arguments
  par <- list(...)
  if(is.null(par[["fit"]])){ fit <- FALSE}else{ fit <- par$fit}
  if(is.null(par[["SE"]])){ SE <- NULL}else{ SE <- par$SE}
  
  # Select the method
  method <- match.arg(method)[1]
  if(method=="RidgeAltapprox"){
    method <- "RidgeAlt"
  }else if(method=="LASSOapprox"){
    method <- "LASSO"
  }
  
  # Select the target
  targM <- match.arg(targM)[1]
  
  # Parameters
  n <- nrow(Y)
  p <- ncol(Y)
  
  # Y must be a matrix (force coercion)
  Y <- as.matrix(Y)
  
  # check for parameters
  if(method=="ML" & p>=n) warning("The covariance matrix is singular, the log-likelihood (and the GIC) is unreliable!!")
  nC <- n
  if(REML==TRUE) n <- n-1
  
  # nb parameter
  mod.par = length(param)
  precalc <- pruning(tr, trans=FALSE)
  eig <- eigen(vcv.phylo(tr), symmetric = TRUE)
  D <- eig$vectors%*%(t(eig$vectors)*(1/sqrt(eig$values))) #tcrossprod(eig$vectors %*% diag(1/sqrt(eig$values)), eig$vectors)
  Yi <- crossprod(D,Y)
  X <- crossprod(D,matrix(1,nrow(Y)))
  beta <- pseudoinverse(X)%*%Yi
  residuals <- Yi-X%*%beta
  S <- crossprod(residuals)/n
  
  # Determinant for the phylogenetic tree
  var_root <- precalc$varRoot
  var_contr <- precalc$varNode
  
  if(REML==TRUE){
    Ccov <- sum(log(var_contr))
  }else{
    Ccov <- sum(log(c(var_root,var_contr)))
  }
  
  # Switch between targets
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
           P <- .makePenalty(S,tuning,Target,targM)$S
           eig <- eigen(P)
           V <- eig$vectors
           d <- eig$values
           Pi <- V%*%(t(V)*(1/d))
           DP <- sum(log(d))
         },
         "RidgeArch"={
           Pi <- (1-tuning)*S + tuning*Target
           eig <- eigen(Pi)
           V <- eig$vectors
           d <- eig$values
           P <- V%*%(t(V)*(1/d))
           DP <- sum(log(d))
         },
         "LASSO"={
           LASSO <- glassoFast(S,tuning)
           P <- LASSO$w
           Pi <- LASSO$wi
           DP <- as.numeric(determinant(P)$modulus)
         },
         "ML"={
           P <- S
           eig <- eigen(P)
           V <- eig$vectors
           d <- eig$values
           Pi <- V%*%(t(V)*(1/d))
           DP <- sum(log(d))
         }
  )
  
  # GIC score
  if(method=="RidgeArch"){
    # Nominal loocv
    XtX <- solve(crossprod(X))
    # hat matrix
    h <- diag(X%*%pseudoinverse(X))
    
    # check for hat score of 1 (e.g. MANOVA design)
    nloo <- 1:nC
    nloo <- nloo[!h+1e-8>=1]
    ncor = length(nloo)
    
    # First and second derivative of the functional (we can use patterned matrix to target some matrix elements)
    # We use the Kronecker-vec identity to speed up the computations
    T1 <- unlist( mclapply(nloo, function(i){
      Sk <- tcrossprod(residuals[i,]) ;
      VSV <- 0.5*(Pi - (1-tuning)*Sk - tuning*Target);
      VSV2 <- 0.5*(Pi - Sk);
      sum(VSV * 2*(P%*%VSV2%*%P))
    }, mc.cores = nbcores))
    
    df = sum(T1)/ncor
    sigma_df <- df
    
  }else if(method=="LASSO" | method=="ML"){
    # LASSO or ML
    Tf2 <- function(S, P) {
      I <- ifelse(P==0,0,1) ;
      t(.vec(S*I))%*%.vec(P%*%(S*I)%*%P)
    }
    
    sigma_df <- (1/(2*n))*sum(sapply(1:nC, function(i){
      Sk <- Yk[i,]%*%t(Yk[i,]) ;
      Tf2(Sk, Pi)})) - (1/2)*Tf2(S,Pi)
    
  }else if(method=="RidgeAlt"){
    # Alternative Ridge
    H <- (1/(0.5*(kronecker(d,d)+tuning)))
    
    # 2) First derivative of the functional
    T1 <- sapply(1:nC, function(i){
      Sk <- tcrossprod(Yk[i,]) ;
      VSV <- .vec(crossprod(V, (0.5*(P - (Sk - tuning*Target) - tuning*Pi))%*%V));
      VSV2 <- .vec(crossprod(V, (0.5*(P - Sk))%*%V));
      sum(VSV * (H*VSV2))
    })
    
    df = sum(T1)/n
    sigma_df <- df
  }
  
  # Number of parameters for the root state:
  # The Information matrix from the Hessian and gradients scores
  XtX <- solve(t(X)%*%(X))
  T2 <- sapply(nloo, function(i){
    gradient <- (X[i,])%*%t(P%*%t(Yi[i,]-X[i,]%*%beta))
    sum(gradient * (XtX%*%gradient%*%Pi))
  })
  beta_df <- sum(T2)
  
  if(!is.null(SE)) mod.par <- mod.par+1
  # LogLikelihood (minus)
  llik <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*DP + n*.tr(S%*%P))
  GIC <- 2*llik + 2*(sigma_df+beta_df+mod.par)
  
  
  # return the results
  results <- list(LogLikelihood=-llik, GIC=GIC, p=p, n=ncor, bias=sigma_df+beta_df+mod.par, bet=beta_df, sig=sigma_df, T2=T2, T1=T1, X=X, Y=Yi, res=residuals)
  class(results) <- "gic.rpanda"
  return(results)
}

# Extractor for fit_pl.rpanda 'class'
GIC3 <- function(object){
  if(!inherits(object,"fit_pl.rpanda")) stop("only works with \"fit_pl.rpanda\" class objects")
  gic_criterion3(Y=object$Y, tr=object$scaled_tree, model=object$model, method=object$method, targM=object$targM, param=object$model.par, tuning=object$gamma, REML=object$REML, fit=TRUE, SE=object$SE)
}
