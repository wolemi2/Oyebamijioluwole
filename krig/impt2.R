#####create multivariate kriging model class
#krig=setClass("krig",representation(formula="formula",design="matrix",response="matrix",d ="integer",name = "character", paramset.n = "integer", var.names = "character",sd2 = "numeric",known.covparam = "character",nugget.flag = "logical", nugget.estim = "logical",nugget = "numeric",param.n = "integer",range.n = "integer",range.names = "character",range.val = "numeric",shape.n = "integer",shape.names = "character",shape.val = "numeric",X="matrix")) 
krig=setClass("krig",representation(d="integer",n = "integer",X ="matrix",y = "matrix",p = "integer",F = "matrix",formula = "formula",trend.coef = "numeric",covar= "cpr",noise.flag = "logical",noise.var = "numeric",known.param = "character",case = "character",param.estim = "logical",method = "character",optim.method = "character",lower = "numeric",upper = "numeric",control = "list",gr = "logical",parinit = "numeric",logLik = "numeric",T = "matrix",z = "numeric",M ="matrix",vn="numeric",C="matrix"))	          
 mod=new("krig")

#########################continuation
krige=function(formula = ~1, design, response, covtype = "gauss",
    coef.trend = NULL, coef.cov = NULL, coef.var = NULL, nugget = NULL,
    nugget.estim = FALSE, noise.var = NULL, estim.method = "MLE",
    penalty = NULL, optim.method = "BFGS", lower = NULL, upper = NULL,
    parinit = NULL, multistart = 1, control = NULL, gr = TRUE,
    iso = FALSE, scaling = FALSE, knots = NULL, kernel = NULL)
{
data <- data.frame(design)
  mod@formula <- formula 
#<- drop.response(formula, data = data)
  F <- model.matrix(formula, data=data)
  
  # X <- as.matrix(design)
  X <- as.matrix(data)
  y <- as.matrix(response)
	mod@X <- X
	mod@y <- y
	mod@d <- ncol(X)
	mod@n <- nrow(X)
	mod@F <- F
	mod@p <- ncol(F)
	mod@vn=mod@noise.var <- as.numeric(noise.var)##cal from repetition and supply
#coef.var <- coef.cov <- NULL
covtype="gauss" ###other covariance like matern and exponential are possible
source("auxilliary.R",echo=TRUE)
mod@covar <- covs(covtype,d=length(nam),var.names=nam,known.covparam="None",coef.var=coef.var,coef.cov=alpha)
#######optimization of parameter set
mod@optim.method <- as.character(optim.method)
bounds <- bbound(mod@covar, design)
lower <- bounds$lower
upper <- bounds$upper
control$multistart <- multistart  #supply by user
  mod@lower <- as.numeric(lower)
  mod@upper <- as.numeric(upper)
  mod@parinit <- as.numeric(parinit)
optim.method == "BFGS"
control$pop.size <- 200
control$pop.size <- max(control$pop.size, multistart)
control$trace <- 3 ##or 0
 control$upper.alpha <- 1 - 1e-08
 mod@control <- control
 mod@gr <- as.logical(gr)
envir.LL <- new.env()
validObject(mod, complete=TRUE)
#f <- kmEstimate
###############from kmEstimate
f=function(mod,envir){
nugget <- mod@noise.var
fn <- LL   ###log likelihood defined elsewhere
fnscale <- -1
gr <- grad###gradient function also defined elsewhere
mod@control$multistart <- multistart <- 1
initList <- krig.init(mod)###defined elsewhere
lower <- mod@lower <- as.numeric(initList$lower)
upper <- mod@upper <- as.numeric(initList$upper)
parinit <- initList$par
lp <- nrow(parinit)
control <- mod@control
cat("  - parameters lower bounds : ", lower[1:(lp-1)], "\n")
 		  cat("  - parameters upper bounds : ", upper[1:(lp-1)], "\n")
 		  cat("  - variance bounds : ", c(lower[lp], upper[lp]), "\n")
		  cat("  - best initial criterion value(s) : ", initList$value, "\n")
#optimization############
BFGSargs <- c("trace", "parscale", "ndeps", "maxit", "abstol", "reltol", "REPORT", "lnm", "factr", "pgtol")
commonNames <- intersect(BFGSargs, names(control))
controlChecked <- control[commonNames]
controlChecked$REPORT <- 1
forced <- list(fnscale = fnscale)
controlChecked[names(forced)] <- forced
multistart <- control$multistart
######
mod@parinit <- parinit <- as.numeric(parinit[, 1])
mod@covar <- initList$cov[[1]]
o <- optim(par = parinit, fn = fn, gr = gr,method = "L-BFGS-B", lower =lower, upper = upper,control = controlChecked, hessian = FALSE, mod, envir=envir.LL)
mod@control$convergence <- o$convergence
mod@logLik <- as.numeric(o$value)
###########################################################################fitting

T <- envir$T
 z <- envir$z
 x <- backsolve(t(T), y, upper.tri=FALSE)		# x:=(T')^(-1)*y
M <- backsolve(t(T), F, upper.tri=FALSE)		# M:=(T')^(-1)*F
l <- lm(x ~ M-1)
beta.hat <- as.numeric(l$coef)
mod@trend.coef <-beta.hat
mod@T <- T
mod@z <- as.numeric(z)
mod@M <- M
param <- as.numeric(o$par)
lp <- length(param)
mod@covar@sd2 <- param[lp]
mod@covar <- vec1(mod@covar, param[1:(lp-1)])
mod@lower <- mod@lower[1:(lp)]    
mod@upper <- mod@upper[1:(lp)]
mod@parinit <- mod@parinit[1:(lp)]
return(mod)
}
mod <- f(mod, envir = envir.LL)
return(mod)
}
#gg2=krige(formula=f1,design,response,coef.cov=alpha,noise.var=rep(0,length(response)),covtype="gauss")
#gg=krige(formula=f1,design,response,noise.var=0)
##########################
#olu=optim(par =gg@parinit, fn =LL, gr =grad,method = "L-BFGS-B", lower =gg@lower, upper = gg@upper,hessian = FALSE,gg,envir=envir.LL)

####################################################################
################end
