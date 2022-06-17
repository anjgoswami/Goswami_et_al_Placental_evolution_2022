## Function to retrieve the rate through time curve for the climatic model

retrieve_curves <- function(object, fun, age, plot=TRUE, ...){
  
  if(missing(age)) age = max(nodeHeights(object$tree))
  if(missing(fun)) fun = object$clim_fun
  time <- seq(0, age, length.out=100)
  args <- list(...)
  
  if(object$is.simmap){
    rates <- sapply(1:object$betaSimmap, function(i) {
      fun(time, object$model.par[seq_len(object$nparam) + object$nparam*(i-1)]) * mean(diag(object$R$R))
    })
    if(plot){
        if(!is.null(args[["col"]])){
            col = args$col
            plot(-time, rev(rates[,1]), type='l', xlab="Times", 
            ylab=bquote(paste("Evolutionary rates ", sigma)), main="reconstructed rates",
            ylim=range(rates), col=col[1])
            for(i in 2:object$betaSimmap) lines(-time, rev(rates[,i]), lty=i, col=col[i])
        }else{
            plot(-time, rev(rates[,1]), type='l', xlab="Times",
            ylab=bquote(paste("Evolutionary rates ", sigma)), main="reconstructed rates",
            ylim=range(rates), ...)
            for(i in 2:object$betaSimmap) lines(-time, rev(rates[,i]), lty=i, ...)
        }
    }
    
  }else{
    rates <- fun(time, object$model.par) * mean(diag(object$R$R)) # Retrieve the rates estimates for each times steps
    if(plot){
        if(!is.null(args[["add"]])) add_plot = args$add else add_plot = FALSE
        if(add_plot) lines(-time, rev(rates), ...)
        else{
             if(is.null(args[["main"]])) plot(-time, rev(rates), type='l', xlab="Times", ylab=bquote(paste("Evolutionary rates ", sigma)), main="reconstructed rates", ...)
             else plot(-time, rev(rates), type='l', xlab="Times", ylab=bquote(paste("Evolutionary rates ", sigma)), ...)
        }
    }
  }
    
  results <- list(rates=rates, times=time, parameters=object$model.par, scaling=mean(diag(object$R$R)))
  invisible(results)
}

