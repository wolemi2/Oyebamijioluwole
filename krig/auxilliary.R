#####
#setClass("cpr",representation(d ="integer",name = "character", paramset.n = "integer", var.names = "character",sd2 = "numeric",known.covparam = "character",nugget.flag = "logical", nugget.estim = "logical",nugget = "numeric",param.n = "integer",range.n = "integer",range.names = "character",range.val = "numeric",shape.n = "integer",shape.names = "character",shape.val = "numeric"))           
 cpr=setClass("cpr",representation(d ="integer",name = "character", paramset.n = "integer", var.names = "character",sd2 = "numeric",known.covparam = "character",param.n = "integer",range.n = "integer",range.val="numeric",range.names = "character"))           
############
###########################GENERIC setting
setGeneric(name= "covMat1Mat2",def= function(obj,X1, X2) standardGeneric("covMat1Mat2"))
setGeneric(name="covMatrix",def= function(obj, X, noise.var=NULL) standardGeneric("covMatrix"))
 setGeneric(name= "bbound",def= function(obj, X) standardGeneric("bbound"))
setGeneric(name = "cov.deri",def = function(obj, X, C0, k, ...) standardGeneric("cov.deri")) 
setGeneric(name = "cov.deri2",def = function(obj, X, C0, k, ...) standardGeneric("cov.deri2"))              
setGeneric(name="vec2",def= function(obj) standardGeneric("vec2"))
setGeneric(name= "vec1",def= function(obj, param) standardGeneric("vec1"))
setGeneric(name= "paramSample",def= function(obj, n, ...) standardGeneric("paramSample"))
               
#######################
vec2.f <- function(obj){
param <- obj@range.val
return(as.numeric(param))
}
###########vect2covparam
vec1.f <- function(obj, param){
obj@range.val <- param
return(obj)
}     
setMethod("vec2",signature = "cpr",definition = vec2.f)
setMethod("vec1",signature = "cpr",definition = vec1.f)
##########covStruct
covs=function(covtype, d, known.covparam, var.names, coef.cov = NULL,coef.var = NULL){
covv <- new("cpr", d = as.integer(d), name = as.character(covtype),sd2 = as.numeric(coef.var), var.names = as.character(var.names),known.covparam = known.covparam)
covv@range.names = "alpha"
 covv@paramset.n <- as.integer(1)
covv@param.n <- as.integer(d)
covv@range.n <- as.integer(d)
#covv@range.val = coef.cov
 if (length(coef.cov)>0) covv <- vec1(covv,coef.cov)
validObject(covv)
return(covv)
}
#covs(covtype,d=length(nam),var.names=nam,known.covparam="None",coef.var=1,coef.cov=alpha)

######covParametersBounds
bbound.f=function(obj,X){
k <- obj@range.n
lower <- rep(1e-10, k)
upper <- 2 * diff(apply(X, 2, range))
upper <- as.vector(upper)
return(list(lower=lower, upper=upper))
}
setMethod("bbound", signature = "cpr",definition = function(obj, X){
if (obj@paramset.n==1) {
k <- obj@range.n
lower <- rep(1e-10, k)
upper <- 2 * diff(apply(X, 2, range))
upper <- as.vector(upper)
            #} else if (identical(obj@name, "powexp")) {              
              # coef. order : theta_1, ..., theta_d, p_1, ..., p_d
              #lower <- rep(1e-10, obj@param.n)
              #upper <- 2 * diff(apply(X, 2, range))
              #k <- obj@shape.n
              #upper <- as.vector(c(upper, rep(2, k)))
} else stop("No default values for covariance parameters bounds, the inputs 'lower' and 'upper' are required")
return(list(lower=lower, upper=upper))})
###########################


##############
#covStruct.create("gauss",d=4,var.names=c("a","b","c","d"),known.covparam="None",coef.var=1,coef.cov=alpha)
#covs("gauss",d=4,var.names=c("a","b","c","d"),known.covparam="None",coef.var=1,coef.cov=alpha)
dyn.load("CovFuns.so")  
######################################
covMatrix.cpr<- function(obj, X, noise.var=NULL) {
  d <- ncol(X)
  n <- nrow(X)
noise.var=mod@noise.var
param <- vec2(obj)
out <- .C("C_covMatrix", 
            as.double(X), as.integer(n), as.integer(d), 
            as.double(param), as.double(obj@sd2), as.character(obj@name), 
            ans = double(n * n))
  C <- matrix(out$ans, n, n)   
if (length(noise.var)>0) {
vn <- noise.var
C <- C + diag(noise.var, nrow = n)
  } else {
vn <- rep(0, n)
  }
return(list(C=C, vn=vn))
}
setMethod("covMatrix", 
          signature = "cpr", 
          definition = function(obj, X, noise.var=NULL) {
            covMatrix.cpr(obj=obj, X=X, noise.var=noise.var)
          }
)
################
covMat1Mat2.f<- function(obj, X1, X2) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  d <- ncol(X1)
 param <- vec2(obj)
#param <- params
  out <- .C("C_covMat1Mat2", 
            as.double(X1), as.integer(n1),
            as.double(X2), as.integer(n2), 
            as.integer(d),
            as.double(param), as.double(obj@sd2),as.character(obj@name),
            ans = double(n1 * n2))
  #as.double(object@sd2), as.character(object@name)
  M <- matrix(out$ans, n1, n2)
  
  #if ((!nugget.flag) | (!mod@nugget.flag)) {
    return(M)
  #} else {
   # out <- .C("C_covMat1Mat2", 
            #  as.double(X1), as.integer(n1),
             # as.double(X2), as.integer(n2), 
              #as.integer(d),
              #as.double(param), as.double(mod@nugget), "whitenoise",
              #ans = double(n1 * n2))
    #N <- matrix(out$ans, n1, n2)
    #return(M+N)
  #}
 }
setMethod("covMat1Mat2",signature = "cpr", definition = function(obj,X1, X2) {covMat1Mat2.f(obj,X1=X1,X2=X2)})
###################
cov.deri.f <- function(obj, X, C0, k) {
  ## X : n x d
n <- nrow(X)
d <- ncol(X)
##index k  starts at 0 in C language 
param <- vec2(obj)
k <- as.integer(k)
if ((k >=1) & (k <= obj@param.n)) {  # derivative with respect to alpha_k
    out <- .C("C_covMatrixDerivative", 
              as.double(X), as.integer(n), as.integer(d), 
              as.double(param), as.character(obj@name),
              as.integer(k), as.double(C0),
              ans = double(n * n))
             return(matrix(out$ans, n, n))
  } else if (k==(obj@param.n+1)) {     # derivative with respect to sigma^2
    return(C0 / obj@sd2)
  } else {
    stop("Wrong value of k")
  }
}
setMethod("cov.deri",signature = "cpr", definition = function(obj, X, C0, k) {
cov.deri.f(obj=obj, X=X, C0=C0, k=k)})
#########
cov.deri2.f <- function(obj, X, C0, k) {
  ## X : n x d
n <- nrow(X)
d <- ncol(X)
##index k  starts at 0 in C language 
param <- vec2(obj)
k <- as.integer(k)
if ((k >=1) & (k <= obj@param.n)) {  # derivative with respect to alpha_k
    out <- .C("C_covMatrixDerivative", 
              as.double(X), as.integer(n), as.integer(d), 
              as.double(param), as.character(obj@name),
              as.integer(k), as.double(C0),
              ans = double(n * n))
             return(matrix(out$ans, n, n))
  } else if (k==(obj@param.n+1)) {     # derivative with respect to sigma^2
    return(C0 / obj@sd2)
  } else {
    stop("Wrong value of k")
  }
}
setMethod("cov.deri2",signature = "cpr", definition = function(obj, X, C0, k) {
cov.deri2.f(obj=obj, X=X, C0=C0, k=k)})






setMethod("paramSample",signature = "cpr",definition = function(obj, n, lower, upper, y=NULL, type="all-sd2-nugget"){
param.n <- obj@param.n
matrixinit <- matrix(runif(n * param.n), nrow=param.n, ncol=n)
matrixinit <- lower + matrixinit*(upper - lower)})
############
compute.z <- function(x, M, beta=NULL){
  if (!is.null(beta)) {
    z <- x - M %*% beta
  } else {   
    Q <- qr.Q(qr(M))
    H <- Q %*% t(Q)
    z <- x - H %*% x
  }
  return(as.numeric(z))
}

#############other auxilliary functions
#source("atta.R",echo=FALSE)
source("derivative.R",echo=FALSE)
source("kmNuggets.init.R",echo=FALSE)
source("LogLikFun.R",echo=FALSE)

###################################################END
eq2.36_Sigma <- function(H, Sigma, d){
betahat <- betahat_mult_Sigma(H,Sigma,d)
sqrt((1/det(Sigma)) / det(quad.form.inv(Sigma,H))) * 
  exp( -0.5*quad.form.inv(Sigma,   d-H %*% betahat))
}

