#############LogLikDerivative
grad=function (param,mod, envir)
{
#else if (identical(mod@case, "LLconcentration_beta")) {
        nparam <- length(param)
        mod@covar <- vec1(mod@covar, param[1:(nparam -
            1)])
        mod@covar@sd2 <- sigma2 <- param[nparam]
        logLik.derivative <- matrix(0, nparam, 1)
        C <- envir$C
        T <- envir$T
        vn <- envir$vn
        z <- envir$z
        x <- backsolve(T, z)
        Cinv <- chol2inv(T)
        for (k in 1:(nparam)) {
            gradC.k <- cov.deri(mod@covar,
                X = mod@X, C0 = C - diag(vn, nrow = nrow(C)),
                k = k)
            term1 <- -t(x) %*% gradC.k %*% x
            term2 <- sum(Cinv * gradC.k)
            logLik.derivative[k] <- -0.5 * (term1 + term2)
        }
return(logLik.derivative)
}
#############################

grad2=function (param,modd,envir=NULL)
{     param=modd@covar@sd2
nparam <- length(param)
modd@covar@param.n=nparam
sigma2 <- param
betahat=function(F,Sigmainv,Y)
{out <- as.vector(solve(t(F) %*%Sigmainv%*%F, crossprod(crossprod(Sigmainv,F), Y)))
names(out) <- colnames(F)
return(out)
}
logLik.derivative <- matrix(0, nparam, 1)
C <- as.matrix(nearPD(modd@C)$mat)#modd@C
Y=modd@y
F=modd@F
#T <- envir$T
vn <- modd@vn
#z <- envir$z
#x <- backsolve(T, z)
Cinv <- chol2inv(C)
beta=betahat(F,Cinv,Y)
        for (k in 1:(nparam)) {
            gradC.k <- cov.deri2(modd@covar,
                X = modd@X, C0 = C - diag(vn, nrow = nrow(C)),
                k = k)
            term1 <- -t(Y-F%*%beta)%*%Cinv %*% gradC.k %*%Cinv%*%(Y-F%*%beta)
            term2 <- sum(Cinv * gradC.k)
            logLik.derivative[k] <- -0.5 * (term1 + term2)
        }
return(logLik.derivative)
}
