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
# model    = either "BM", "OU", "EB", or "lambda"; the model of traits evolution
# method   = "RidgeArch": Archetype (linear) Ridge penalty, "RidgeAlt": Quadratic Ridge penalty, "LASSO": Least Absolute Selection and Shrinkage Operator. "RidgeAltapprox" and "LASSOapprox" are fast approximations of the LOOCV for the Ridge quadratic and LASSO penalties
# targM    = "null", "Variance" for a diagonal unequal variance target, "unitVariance" for an equal diagonal target. Only works with "RidgeArch","RidgeAlt", and "RidgeAltapprox" methods.
# REML     = TRUE (default) or FALSE. The likelihood method used to estimate the parameters. REML must be preferred with small sample size in order to obtain unbiased estimates.
# up       = upper bound for model parameter search
# low      = lower bound for model parameter search
# tol      = lower bound tolerance for the regularizer (tuning parameter)
# starting = starting values for the parameter search. Must be a vector likes: c(model, regularization)


require(mvMORPH)    # >= 1.0.9
require(glassoFast) # https://github.com/JClavel/glassoFast
require(parallel)
#require(optimr)

loocvPhylo3 <- function(Y, tree, model=c("BM","OU","EB","lambda","EnvExp","EnvLin","multiEnv", "multiEnv2","multiEnv3","EnvExp2","EnvExp3","EnvExp4","EnvExp5","EnvExp6","EnvExp7"), method=c("RidgeAlt","RidgeArch","RidgeAltapprox","LASSO","LASSOapprox"), targM=c("null","Variance","unitVariance"), REML=TRUE, up=NULL, low=NULL, tol=NULL, starting=NULL, SE=NULL, scale.height=TRUE, ncores=3L, funEnv=NULL, funEnv2=NULL, opt="L-BFGS-B", control=list(), trend=function(x) x){
  
  # Preliminary checks
  if(missing(tree)) stop("Please provide a phylogenetic tree of class \"phylo\" ")
  if(!is.ultrametric(tree) & model=="OU") stop("The method is not working with non-ultrametric trees")
  if(nrow(Y) != Ntip(tree)) stop("Length of phenotypic and phylogenetic data do not match")
  if (all(rownames(Y) %in% tree$tip.label)){
    Y <- Y[tree$tip.label,]
  }else{
    warning("rownames in Y are missing. It is assumed that they are in the same order as in the tree.")
  }
  
  # Select the model
  model <- match.arg(model)[1]
  
  # Select the method
  method <- match.arg(method)[1]
  
  # Select the method
  targM <- match.arg(targM)[1]
  
  if(!inherits(tree,"simmap") & model=="multiEnv") stop("The tree must be in simmap format")
  if(inherits(tree,"simmap")){
    betaSimmap = ncol(tree$mapped.edge)
    tree = reorderSimmap(tree, order="postorder")
  }else{
    tree = reorder(tree, order = "postorder")
  }
  
  # Bounds for models
  if(is.null(up)){
    switch(model,
           "EB"={up <- 0},
           "OU"={up <- 10},
           "lambda"={up <- 1.1},
           "EnvExp"={up <- 50},
           "EnvExp2"={up <- 50},
           "EnvExp3"={up <- 50},
           "EnvExp4"={up <- 50},
           "EnvExp5"={up <- 50},
           "EnvExp6"={up <- 50},
           "EnvExp7"={up <- 50},
           "multiEnv"={up <- rep(50, betaSimmap)},
           "multiEnv2"={ up <- rep(50, 2*betaSimmap)},
           "multiEnv3"={ up <- rep(50, 3*betaSimmap)})
  }
  
  if(is.null(low)){
    switch(model,
           "EB"={low <- -10},
           "OU"={low <- 0},
           "lambda"={low <- 1e-5},
           "EnvExp"={low <- -50},
           "EnvExp2"={low <- -50},
           "EnvExp3"={low <- -50},
           "EnvExp4"={low <- -50},
           "EnvExp5"={low <- -50},
           "EnvExp6"={low <- -50},
           "EnvExp7"={low <- -50},
           "multiEnv"={ low <- rep(-50, betaSimmap)},
           "multiEnv2"={ low <- rep(-50, 2*betaSimmap)},
           "multiEnv3"={ low <- rep(-50, 3*betaSimmap)})
  }
  names_param = NULL
  # Reorder the tree
  tree<-reorder.phylo(tree,"postorder")
  
  # Parameters
  if(ncol(Y)==1) stop("Only works with multivariate datasets")
  n <- nO <- nrow(Y)
  nC <- n-1
  if(REML==TRUE) n <- n-1
  
  # Scale the tree to unit length
  if(scale.height==TRUE) tree$edge.length <- tree$edge.length/max(branching.times(tree))
  
  p <- ncol(Y)
  Yest <- apply(Y,2,function(i) pic(i,tree))
  # Empirical covariance (contrasts) for starting the algorithm
  S <- t(Yest)%*%Yest / n
  
  # Identity matrix
  I <- diag(p)
  
  # Default penalty is Ridge "null"
  target <- matrix(0,p,p)
  
  # Identifying tips values
  tipsIndices <- which(  tree$edge[, 2] <= Ntip(tree))
  
  
  
  # Default tolerance for the parameter search
  if(is.null(tol)){
    if(method=="RidgeArch"){
      tol = 1e-8
    }else{
      tol = 0
    } 
  }
  
  # Number of cores to parallelize the computations
  mcores <- ncores
  
  ## ---- Transform phylo
  switch(model,
         "lambda"={
           phyTrans <- function(phy, lambda) {
             if(lambda==1) return(phy)
             # Compute the branching times
             if(is.ultrametric(phy)){
                 times<-branching.times(phy)
                  rootOrig <- max(times)
             }else{
                 # Use "phytools" called by mvMORPH
                 times<-max(nodeHeights(phy))-nodeHeights(phy)[match(1:phy$Nnode+Ntip(phy),phy$edge[,1]),1]
                  rootOrig <- max(times)
             }
            
             tips <- match(c(1:Ntip(phy)), phy$edge[,2])
             phy$edge.length <- phy$edge.length * lambda
             phy$edge.length[tips] <- phy$edge.length[tips] + (rootOrig * (1-lambda))
             return(phy)
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,log(1e6))
             lowerBound <- c(low,log(tol))
           }else{
             upperBound <- c(up,1)
             lowerBound <- c(low,tol)
           }
           
           idx1 <- 1
           idx2 <- 2
           idx3 <- 3
         },
         "OU"={
           phyTrans<-function(phy,alpha){
             if(alpha<=.Machine$double.eps) return(phy) # reduce to BM
             
             # Compute the branching times
             if(is.ultrametric(phy)){
                 times<-branching.times(phy)
                 names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
             }else{
                 # Use "phytools" called by mvMORPH
                 times<-max(nodeHeights(phy))-nodeHeights(phy)[match(1:phy$Nnode+Ntip(phy),phy$edge[,1]),1]
                 names(times)<-1:phy$Nnode+Ntip(phy)
             }
             
             Tmax<-times[1]
             phy2<-phy
             
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[which(names(times) == phy$edge[i, 1])]
               t1 <- max(times) - age
               t2 <- t1+bl
               phy2$edge.length[i] <- (1/(2*alpha))*exp(-2*alpha * (Tmax-t2)) * (1 - exp(-2 * alpha * t2)) -
                 (1/(2*alpha))*exp(-2*alpha * (Tmax-t1)) * (1 - exp(-2 * alpha * t1))
             }
             phy <- phy2
             return(phy)
           }
           
           
           transform <- function(x) exp(x)
           
           if(method!="RidgeArch"){
             upperBound <- c(log(up),log(1e6))
             lowerBound <- c(log(low),log(tol))
           }else{
             upperBound <- c(log(up),1)
             lowerBound <- c(log(low),tol)
           }
           
           idx1 <- 1
           idx2 <- 2
           idx3 <- 3
         },
         "EB"={
           phyTrans <- function (phy, beta)
           {
             #if(abs(beta)<=.Machine$double.eps) return(phy)
             if(beta==0) return(phy)
             
             # Compute the branching times
             if(is.ultrametric(phy)){
                 times<-branching.times(phy)
                 names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
             }else{
                 # Use "phytools" called by mvMORPH
                 times<-max(nodeHeights(phy))-nodeHeights(phy)[match(1:phy$Nnode+Ntip(phy),phy$edge[,1]),1]
                 names(times)<-1:phy$Nnode+Ntip(phy)
             }
             
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age = times[which(names(times) == phy$edge[i, 1])]
               t1 = max(times) - age
               t2 = t1+bl
               phy$edge.length[i] = (exp(beta*t2)-exp(beta*t1))/(beta)
             }
             return(phy)
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,log(1e6))
             lowerBound <- c(low,log(tol))
           }else{
             upperBound <- c(up,1)
             lowerBound <- c(low,tol)
           }
           
           idx1 <- 1
           idx2 <- 2
           idx3 <- 3
         },
         "BM"={
           phyTrans <- function(phy, beta) return(phy)
           transform <- function(x) x
           
           
           if(method!="RidgeArch"){
             upperBound <- log(1e6)
             lowerBound <- log(tol)
           }else{
             upperBound <- 1
             lowerBound <- tol
           }
           
           idx1 <- idx2 <- 1
           idx3 <- 2
         },
         "EnvExp"={
           
           # define the function for the EnvExp model
           # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
           f<-function(x){exp(beta*funEnv((mtot)-x))} # assumes ultrametric tree or tree with present day taxa
           
           phyTrans <- function(phy, beta){
             
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
             mtot<-max(times)
             # Set the root to zero
             times<-max(times)-times
             res <- phy
             
             # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
             f<-function(x){exp(beta*funEnv((mtot)-x))} # assumes ultrametric tree or tree with present day taxa
             # Transforms the branch-lengths of the tree
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[phy$edge[i, 1] - nO]
               int <- try(integrate(f, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
               # Try catch if the integrand is divergent
               if(inherits(int ,'try-error')){
                 warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                 integ <- NA_real_
               } else {
                 integ <- int$value
               }
               res$edge.length[i] <- integ
               # if(i==1)  print(integ)
             }
             phy<-res
             return(phy)        
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,log(1e6))
             lowerBound <- c(low,log(tol))
           }else{
             upperBound <- c(up,1)
             lowerBound <- c(low,tol)
           }
           
           idx1 <- 1
           idx2 <- 2
           idx3 <- 3
         },
         "EnvExp2"={
           
           phyTrans <- function(phy, beta){
             
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
             mtot<-max(times)
             # Set the root to zero
             times<-max(times)-times
             res <- phy
             
             # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
             #f<-function(x){exp(beta[1]*funEnv((mtot)-x) + beta[2]*x)} # assumes ultrametric tree or tree with present day taxa | or beta[2]*(mtot-x)?
             f<-function(x){exp(beta[1]*funEnv((mtot)-x) + beta[2]*trend(x))}
             # Transforms the branch-lengths of the tree
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[phy$edge[i, 1] - nO]
               int <- try(integrate(f, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
               # Try catch if the integrand is divergent
               if(inherits(int ,'try-error')){
                 warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                 integ <- NA_real_
               } else {
                 integ <- int$value
               }
               res$edge.length[i] <- integ
               
             }
             phy<-res
             return(phy)        
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,up,log(1e6))
             lowerBound <- c(low,low,log(tol))
           }else{
             upperBound <- c(up,up,1)
             lowerBound <- c(low,low,tol)
           }
           
           idx1 <- 1:2
           idx2 <- 3
           idx3 <- 4
         }, 
         "EnvExp3"={
           
           phyTrans <- function(phy, beta){
             
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
             mtot<-max(times)
             # Set the root to zero
             times<-max(times)-times
             res <- phy
             
             # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
             f<-function(x){exp(beta[1]*funEnv((mtot)-x) + beta[2]*trend(x) + beta[3]*funEnv((mtot)-x)*trend(x))} # assumes ultrametric tree or tree with present day taxa | or beta[2]*(mtot-x)?
             
             # Transforms the branch-lengths of the tree
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[phy$edge[i, 1] - nO]
               int <- try(integrate(f, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
               # Try catch if the integrand is divergent
               if(inherits(int ,'try-error')){
                 warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                 integ <- NA_real_
               } else {
                 integ <- int$value
               }
               res$edge.length[i] <- integ
               
             }
             phy<-res
             return(phy)        
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,up,up,log(1e6))
             lowerBound <- c(low,low,low,log(tol))
           }else{
             upperBound <- c(up,up,up,1)
             lowerBound <- c(low,low,low,tol)
           }
           
           idx1 <- 1:3
           idx2 <- 4
           idx3 <- 5
         },
         "EnvExp4"={
           
           phyTrans <- function(phy, beta){
             
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
             mtot<-max(times)
             # Set the root to zero
             times<-max(times)-times
             res <- phy
             
             # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
             f<-function(x){exp(beta[1])*exp(beta[2]*funEnv((mtot)-x) + beta[3]*x)} # assumes ultrametric tree or tree with present day taxa | or beta[2]*(mtot-x)?
             
             # Transforms the branch-lengths of the tree
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[phy$edge[i, 1] - nO]
               int <- try(integrate(f, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
               # Try catch if the integrand is divergent
               if(inherits(int ,'try-error')){
                 warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                 integ <- NA_real_
               } else {
                 integ <- int$value
               }
               res$edge.length[i] <- integ
               
             }
             phy<-res
             return(phy)        
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,up,up,log(1e6))
             lowerBound <- c(low,low,low,log(tol))
           }else{
             upperBound <- c(up,up,up,1)
             lowerBound <- c(low,low,low,tol)
           }
           
           idx1 <- 1:3
           idx2 <- 4
           idx3 <- 5
         },
         "EnvExp5"={
           
           phyTrans <- function(phy, beta){
             
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
             mtot<-max(times)
             # Set the root to zero
             times<-max(times)-times
             res <- phy
             
             # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
             f<-function(x){exp(beta[1]*funEnv((mtot)-x) + beta[2]*trend(x) + beta[3]*funEnv2((mtot)-x))} # assumes ultrametric tree or tree with present day taxa | or beta[2]*(mtot-x)?
             
             # Transforms the branch-lengths of the tree
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[phy$edge[i, 1] - nO]
               int <- try(integrate(f, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
               # Try catch if the integrand is divergent
               if(inherits(int ,'try-error')){
                 warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                 integ <- NA_real_
               } else {
                 integ <- int$value
               }
               res$edge.length[i] <- integ
               
             }
             phy<-res
             return(phy)        
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,up,up,log(1e6))
             lowerBound <- c(low,low,low,log(tol))
           }else{
             upperBound <- c(up,up,up,1)
             lowerBound <- c(low,low,low,tol)
           }
           
           idx1 <- 1:3
           idx2 <- 4
           idx3 <- 5
         },
         "EnvExp6"={
             
             phyTrans <- function(phy, beta){
                 
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
                 mtot<-max(times)
                 # Set the root to zero
                 times<-max(times)-times
                 res <- phy
                 
                 # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
                 f<-function(x){exp(beta[1]*funEnv((mtot)-x) + beta[2]*trend(x) + beta[3]*funEnv2((mtot)-x) + beta[4]*funEnv((mtot)-x)*funEnv2((mtot)-x))} # assumes ultrametric tree or tree with present day taxa | or beta[2]*(mtot-x)?
                 
                 # Transforms the branch-lengths of the tree
                 for (i in 1:length(phy$edge.length)) {
                     bl <- phy$edge.length[i]
                     age <- times[phy$edge[i, 1] - nO]
                     int <- try(integrate(f, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
                     # Try catch if the integrand is divergent
                     if(inherits(int ,'try-error')){
                         warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                         integ <- NA_real_
                     } else {
                         integ <- int$value
                     }
                     res$edge.length[i] <- integ
                     
                 }
                 phy<-res
                 return(phy)
             }
             
             transform <- function(x) x
             
             if(method!="RidgeArch"){
                 upperBound <- c(up,up,up,up,log(1e6))
                 lowerBound <- c(low,low,low,low,log(tol))
             }else{
                 upperBound <- c(up,up,up,up,1)
                 lowerBound <- c(low,low,low,low,tol)
             }
             
             idx1 <- 1:4
             idx2 <- 5
             idx3 <- 6
         },
         "EnvExp7"={ # without the trend => clim1 + clim2 + clim1*clim2
             
             phyTrans <- function(phy, beta){
                 
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
                 mtot<-max(times)
                 # Set the root to zero
                 times<-max(times)-times
                 res <- phy
                 
                 # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
                 f<-function(x){exp(beta[1]*funEnv((mtot)-x) + beta[2]*funEnv2((mtot)-x) + beta[3]*funEnv((mtot)-x)*funEnv2((mtot)-x))} # assumes ultrametric tree or tree with present day taxa | or beta[2]*(mtot-x)?
                 
                 # Transforms the branch-lengths of the tree
                 for (i in 1:length(phy$edge.length)) {
                     bl <- phy$edge.length[i]
                     age <- times[phy$edge[i, 1] - nO]
                     int <- try(integrate(f, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
                     # Try catch if the integrand is divergent
                     if(inherits(int ,'try-error')){
                         warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                         integ <- NA_real_
                     } else {
                         integ <- int$value
                     }
                     res$edge.length[i] <- integ
                     
                 }
                 phy<-res
                 return(phy)
             }
             
             transform <- function(x) x
             
             if(method!="RidgeArch"){
                 upperBound <- c(up,up,up,log(1e6))
                 lowerBound <- c(low,low,low,log(tol))
             }else{
                 upperBound <- c(up,up,up,1)
                 lowerBound <- c(low,low,low,tol)
             }
             
             idx1 <- 1:3
             idx2 <- 4
             idx3 <- 5
         },
         "EnvLin"={
           phyTrans <- function(phy, beta){
             
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
             mtot<-max(times)
             # Set the root to zero
             times<-max(times)-times
             res <- phy
             
             # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
             f<-function(x){1 + beta[1]*funEnv((mtot)-x) + beta[2]*trend(x)} # assumes ultrametric tree or tree with present day taxa
             
             # Transforms the branch-lengths of the tree
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[phy$edge[i, 1] - nO]
               int <- try(integrate(f, lower=age, upper=age+bl, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
               # Try catch if the integrand is divergent
               if(inherits(int ,'try-error')){
                 warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                 integ <- NA_real_
               } else {
                 integ <- int$value
               }
               res$edge.length[i] <- integ
             }
             phy<-res
             return(phy)        
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,up,log(1e6))
             lowerBound <- c(low,low,log(tol))
           }else{
             upperBound <- c(up,up,1)
             lowerBound <- c(low,low,tol)
           }
           
           idx1 <- 1:2
           idx2 <- 3
           idx3 <- 4
         },
         "multiEnv"={
           
           if(!inherits(tree,"simmap")) stop("The tree must be in simmap format")
           
           
           phyTrans <- function(phy, beta){
             
             # not ideal, but for debugging now (ultrametric tree only)
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
             mtot<-max(times)
             # Set the root to zero
             times<-max(times)-times
             res <- phy
             
             # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
             f<-function(x, beta){exp(beta*funEnv((mtot)-x))} # assumes ultrametric tree or tree with present day taxa
             
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
                 # bet <- beta[1]#
                 bet <- beta[regimenumber]           # select the corresponding parameter for beta
                 bl <- currentmap[[betaval]]           # branch length under the current map
                 
                 int <- try(integrate(f, lower=age, upper=age+bl, beta=bet, subdivisions=2000L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
                 
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
             
             return(res)        
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,log(1e6))
             lowerBound <- c(low,log(tol))
           }else{
             upperBound <- c(up,1)
             lowerBound <- c(low,tol)
           }
           
           idx1 <- 1:betaSimmap
           idx2 <- length(idx1)+1
           idx3 <- length(idx1)+2
           
         },
         "multiEnv3"={
             
             if(!inherits(tree,"simmap")) stop("The tree must be in simmap format")
             
             
             phyTrans <- function(phy, beta){
                 
                 # not ideal, but for debugging now (ultrametric tree only)
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
                 mtot<-max(times)
                 # Set the root to zero
                 times<-max(times)-times
                 res <- phy
                 
                 # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
                 f<-function(x, beta){exp(beta[1]*funEnv((mtot)-x) + beta[2]*funEnv2((mtot)-x) + beta[3]*funEnv((mtot)-x)*funEnv2((mtot)-x))} # assumes ultrametric tree or tree with present day taxa
                 
                 
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
                         # bet <- beta[1]#
                         bet <- beta[c(regimenumber, regimenumber+betaSimmap, regimenumber+2*betaSimmap)]
                         bl <- currentmap[[betaval]]           # branch length under the current map
                         
                         int <- try(integrate(f, lower=age, upper=age+bl, beta=bet, subdivisions=2000L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
                         
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
                 
                 return(res)
             }
             
             transform <- function(x) x
             
             if(method!="RidgeArch"){
                 upperBound <- c(up,log(1e6))
                 lowerBound <- c(low,log(tol))
             }else{
                 upperBound <- c(up,1)
                 lowerBound <- c(low,tol)
             }
             
             idx1 <- 1:(betaSimmap*3)
             idx2 <- length(idx1)+1
             idx3 <- length(idx1)+2
             
         },

         "multiEnv2"={
           
           if(!inherits(tree,"simmap")) stop("The tree must be in simmap format")
           
           
           phyTrans <- function(phy, beta){
             
             # not ideal, but for debugging now (ultrametric tree only)
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
             mtot<-max(times)
             # Set the root to zero
             times<-max(times)-times
             res <- phy
             
             # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
             #f1<-function(x, beta){exp(beta[1]*funEnv((mtot)-x) + beta[2]*x)} # assumes ultrametric tree or tree with present day taxa
             f1<-function(x, beta){exp(beta[1]*funEnv((mtot)-x) + beta[2]*funEnv2((mtot)-x))} # assumes ultrametric tree or tree with present day taxa
             
             # # Transforms the branch-lengths of the tree
             for (i in 1:length(phy$edge.length)) {
               
               age <- times[phy$edge[i, 1] - nO]
               # for simmap
               currentmap<-phy$maps[[i]]           # retrieve the corresponding maps on the edge "i"
               indlength<-length(currentmap)       # How many mapping there are?
               tempedge<-numeric(indlength)        # temporary vector for the mappings
               
               # boucle pour traverser les "maps"
               for(betaval in 1:indlength){
                 
                 regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[betaval])  # retrieve the regimes within maps
                 
                 bet <- beta[c(regimenumber, regimenumber+betaSimmap)]           # select the corresponding parameter for beta
                 bl <- currentmap[[betaval]]           # branch length under the current map
                 
                 int <- try(integrate(f1, lower=age, upper=age+bl, beta=bet, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent = TRUE)
                 
                 # Try catch if the integrand is divergent
                 if(inherits(int ,'try-error')){
                   warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                   integ <- NA_real_
                 } else {
                   integ <- int$value
                 }
                 
                 tempedge[betaval] <- integ 
                 # on met à jour age parcequ'on va passer au maps suivant pour la lignée i
                 # update "age" because we're moving to the next map for lineage i.
                 age<-age+bl
                 
               }# end of simmap branch loop
               
               res$edge.length[i] <- sum(tempedge)
             }
             
             return(res)        
           }
           
           transform <- function(x) x
           
           if(method!="RidgeArch"){
             upperBound <- c(up,log(1e6))
             lowerBound <- c(low,log(tol))
           }else{
             upperBound <- c(up,1)
             lowerBound <- c(low,tol)
           }
           
           idx1 <- 1:(betaSimmap*2)
           idx2 <- length(idx1)+1
           idx3 <- length(idx1)+2
           
         })
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
             
             if (any(is.na(tr$edge.length))) return(1e6)   # if the tree is bad behaved
             Yk <-  apply(Y,2,function(i) pic(i,tr))
             Sk <- crossprod(Yk)/n
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
             
             if(targM=="Variance"){
               target <- diag(diag(Sk))
             }else if(targM=="unitVariance"){
               target <- I*mean(diag(Sk))
             }
             
             # Hoffbeck & Landgrebe 1996 parameterization
             beta <- (1 - alpha)/(n - 1)
             G <- n*beta*Sk + alpha * target
             # Gi <- try(solve(G), silent = TRUE)
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
    
    message("Initialization via grid search. Please wait...","\n")
    if(method=="RidgeArch"){
      range_val <- c(1e-6, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9)
    }else{
      range_val <- log(c(1e-12, 1e-9, 1e-6, 0.01, 0.1, 1, 10, 100, 1000, 10000))
    }
    switch(model,
           "lambda"={
             mod_val <- c(0.2,0.5,0.8)
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val,range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val,range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             tuning <- start[2]
           },
           "OU"={
             mod_val <- log(log(2)/(max(branching.times(tree))/c(0.1,0.5,1.5,3,8)))
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val,range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val,range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             tuning <- start[2]
           },
           "EB"={
             mod_val <- -log(2)/(max(branching.times(tree))/c(0.1,0.5,1.5,3,8))
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val,range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val,range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             tuning <- start[2]
           },
           "BM"={
             mod_val = NULL
             if(!is.null(SE)){
               brute_force <- expand.grid(range_val,error_guess)
             }else{
               brute_force <- expand.grid(range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             tuning <- start[1]
           },
           "EnvExp"={
             mod_val <- c(-8,-5,-3,-1,-0.5,0,0.5,1,3,5)
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val,range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val,range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             print(start)
             tuning <- start[2]
           },
           "EnvExp2"={
             mod_val <- c(-5,-3,-1.5,-0.5,0,0.5,1,3)/2
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val, mod_val,range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val, mod_val,range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             print(start)
             tuning <- start[3]
           },
           "EnvExp3"={
             mod_val <- c(-3,-0.5,0,0.5,3)/4
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val, mod_val, mod_val, range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val, mod_val, mod_val, range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             print(start)
             tuning <- start[4]
           },
           "EnvExp4"={
             mod_val <- c(-3,-0.5,0,0.5,3)/4
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val, mod_val, mod_val, range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val, mod_val, mod_val, range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             print(start)
             tuning <- start[4]
           },
           "EnvExp5"={
             mod_val <- c(-3,-0.5,0,0.5,3)/4
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val, mod_val, mod_val, range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val, mod_val, mod_val, range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             print(start)
             tuning <- start[4]
           },
           "EnvExp6"={
               mod_val <- c(-3,-0.5,0,0.5,3)/4
               if(!is.null(SE)){
                   brute_force <- expand.grid(mod_val, mod_val, mod_val, mod_val, range_val, error_guess)
               }else{
                   brute_force <- expand.grid(mod_val, mod_val, mod_val, mod_val, range_val)
               }
               fit_st <- apply(brute_force,1,loocv)
               start <- brute_force[which.min(fit_st),]
               print(start)
               tuning <- start[5]
           },
           "EnvExp7"={
               mod_val <- c(-3,-0.5,0,0.5,3)/4
               if(!is.null(SE)){
                   brute_force <- expand.grid(mod_val, mod_val, mod_val, range_val,error_guess)
               }else{
                   brute_force <- expand.grid(mod_val, mod_val, mod_val, range_val)
               }
               fit_st <- apply(brute_force,1,loocv)
               start <- brute_force[which.min(fit_st),]
               print(start)
               tuning <- start[4]
           },
           "EnvLin"={
             mod_val <- c(-0.5,-0.1,0,0.5,1)
             if(!is.null(SE)){
               brute_force <- expand.grid(mod_val, mod_val/4,range_val,error_guess)
             }else{
               brute_force <- expand.grid(mod_val,mod_val/4,range_val)
             }
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             tuning <- start[2]
           },
           "multiEnv"={
             #param = lapply(1:betaSimmap, function(x) c(-5,-3,-1,0,1,3,5) )
             param = lapply(1:betaSimmap, function(x) c(-0.5,0,0.5) ) # particularly if the curve is unscaled may help as all the parameters are more or less on the same scale?
             
             if(!is.null(SE)){
               param[[betaSimmap+1]] = range_val
               param[[betaSimmap+2]] = error_guess
               
               brute_force <- expand.grid(param)
               
             }else{
               param[[betaSimmap+1]] = range_val
               
               brute_force <- expand.grid(param)
               
             }
             # print(brute_force)
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             
             # print(fit_st)
             print(start)
             tuning <- start[idx2]
           },
           "multiEnv2"={
             
             
             if(!is.null(SE)){
               param = c(as.list(rep(0, 2*betaSimmap)), list(range_val), list(error_guess))
               brute_force <- expand.grid(param)
               
             }else{
               param = c(as.list(rep(0, 2*betaSimmap)), list(range_val))
               brute_force <- expand.grid(param)
               
             }
             
             # param = lapply(1:betaSimmap, function(x) c(-0.5,0,0.5) ) # if the curve is unscaled may help as all the parameters are more or less on the same scale?
             # 
             # if(!is.null(SE)){
             #   param[[2*betaSimmap+1]] = range_val
             #   param[[2*betaSimmap+2]] = error_guess
             #   
             #   brute_force <- expand.grid(param)
             #   
             # }else{
             #   param[[2*betaSimmap+1]] = range_val
             #   
             #   brute_force <- expand.grid(param)
             #   
             # }
            # print(brute_force)
             fit_st <- apply(brute_force,1,loocv)
             start <- brute_force[which.min(fit_st),]
             
             # print(fit_st)
             print(start)
             tuning <- start[idx2]
           },
           "multiEnv3"={
               
               if(!is.null(SE)){
                   param = c(as.list(rep(0, 3*betaSimmap)), list(range_val), list(error_guess))
                   brute_force <- expand.grid(param)
                   
               }else{
                   param = c(as.list(rep(0, 3*betaSimmap)), list(range_val))
                   brute_force <- expand.grid(param)
                   
               }
               #print(brute_force)
               fit_st <- apply(brute_force,1,loocv)
               start <- brute_force[which.min(fit_st),]
               
               # print(fit_st)
               #print(start)
               tuning <- start[idx2]
           })
    if(method=="RidgeArch"){
      cat("Best starting for the tuning: ",as.numeric(tuning))
    }else{
      cat("Best starting for the tuning: ",as.numeric(exp(tuning)))
    }
  }else{     
    start <- starting#.starting_val(starting, model, method, SE)
  }
  # Initial guesses found we start the optimization
  message("\n", "Start optimization. Please wait...")
  
  # Optimization of the cross-validated likelihood
  estimModel <- optim(start, fn = loocv, method=opt, upper=upperBound, lower=lowerBound, control=control)
  
  # Compute the scaled tree
  phy_estim <- phyTrans(tree, transform(estimModel$par[idx1]))
  if(!is.null(SE)){
    SE = (estimModel$par[idx3])*(estimModel$par[idx3])
    phy_estim$edge.length[tipsIndices] <- phy_estim$edge.length[tipsIndices] + SE
  }
  
  # Estimated value for the model parameter
  if(model=="BM"){
    model.par <- 0
  }else{
    model.par <- transform(estimModel$par[idx1])
  }
  
  if(model=="multiEnv") names(model.par) <- names_param <- colnames(tree$mapped.edge)
  if(model=="multiEnv2") names(model.par) <- names_param <- as.vector(t(sapply(colnames(tree$mapped.edge), function(x) paste(c("beta1","beta2"), x))))
  if(model=="multiEnv3") names(model.par) <- names_param <- as.vector(t(sapply(colnames(tree$mapped.edge), function(x) paste(c("beta1","beta2","beta3"), x))))

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
  results <- list(loocv=estimModel$value, model.par=model.par, gamma=gamma, scaled_tree=phy_estim, model=model, method=method, p=p, n=nO, targM=targM, R=regularizedEstimates, REML=REML, Y=Y, SE=SE, fun=loocv, ncores=mcores, opt=estimModel, names_param=names_param)
  class(results) <- "fit_pl.rpanda"
  return(results)
}

## -------------- Miscellaneous functions

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
