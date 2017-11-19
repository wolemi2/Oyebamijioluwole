###################Multivariate kriging
#############likelihood equation for multivariate matrix
LL2<- function(mod, ...) {
logLik <-  -0.5*(mod@n*log(2*pi) + 2*sum(log(diag(mod@T))) + t(mod@z)%*%mod@z)     
return(logLik)
}
##########VAR.MATRIX
nobs=rep(NA,lb)
for(i in 1:lb){
nobs[i]=c(nrow(mm[[i]]))
nob=sum(nobs)
}
var.mat=function(model){
Sigma <- matrix(NA,nob,nob)
for (i in 1:lb) {
for (j in 1:lb) {
ii <- which(i ==model[[i]]@T)
jj <- which(j == model[[j]]@T)
Sigma[ii, jj]=sum(model[[i]]@T,model[[j]]@T)/2
}}
return(Sigma)
}
###########
olu=list()
for (i in 1:lb) {
#for (j in 1:lb) {
olu[[i]]=c(model[[i]]@T)
olu2=c(olu)
}
###################

output_type1=c("temp","rain","humidity")
output_type2=c("temp","rain","humidity")##could be different

var.matrix2=function (x1, x2 = x1,Bo,Mo){
type1=types(x1)
type2=types(x2)
#Sigma <- matrix(NA, nrow(x1), nrow(x2))
rownames(Sigma) <- rownames(x1)
colnames(Sigma) <- rownames(x2)
for (i in output_type) {
for (j in output_type) {
ni <- which(i == output_type1)
nj <- which(j == output_type2)
ii <- which(i == type1)
jj <- which(j == type2)
p=1:ninput
x1=as.matrix(x1);x2=as.matrix(x2)
Sigma[ii, jj]=CovMat1Mat2(mod@covar,XX[[i]],XX[[j]],nugget.flag=FALSE)
#Sigma[ii, jj] <- ss(Bo[, , ni], Bo[, , nj]) * Mo[ni,nj] * corr.matrix(x2[jj,p],x1[ii,p], pos.def.matrix = solve(solve(Bo[,, i])/2 + solve(Bo[, , j])/2))
}}
return(Sigma)
}
###########
make_M <- function(vec){  #returns an 'M' matrix from a vector.
jj_M <- diag(alpha)
jj_M[upper.tri(jj_M,diag=TRUE )] <- vec
jj_M[lower.tri(jj_M,diag=FALSE)] <- 0
jj_M <- jj_M + t(jj_M)
diag(jj_M) <- diag(jj_M)/2
return(jj_M)
}
#######################################
sigmahat=rep(NA,lb)
for(i in 1:lb){
sigmahat[i]=model[[i]]@covar@sd2
}

opt_M=function(mm,dd,Mo,Bo){
jj <- sigmahat
jjM <- diag(Mo)
corre <- Mo/ sqrt(kronecker(jjM,t(jjM)))###find new correlatn matrix
M <- sqrt(kronecker(jj,t(jj)))*corre    ##obtain covariance matrix  
start_vec <- M[upper.tri(M,diag=TRUE)]
opt_vec <- optim(par=start_vec,fn=f, ...)
return(make_M(opt_vec$par))
}



######################################################
model.matrix(f1, data=as.data.frame(design))
regressor.multi(design,func=form[[1]])

     optimal_params(toy_expt,toy_LoF,toy_mhp,option='b',control=list(maxit=1))
     
gg2=km(formula=f1,design,response,coef.cov=alpha,noise.var=rep(0,length(response)),covtype="gauss")


##alphas
1/sqrt(2*diag(Bo[[1]]))
0.2416841 0.2913583 0.6565322 0.2746175 
1/sqrt(2*diag(Bo[[2]]))
0.3230853 0.3247849 0.2349205 0.2255320
1/sqrt(2*diag(Bo[[3]]))
0.3230853 0.3247849 0.2286947 0.2302656
Mo
         temp rain humidity
temp      1.0 -0.7      0.5
rain     -0.7  2.0      1.4
humidity  0.5  1.4      3.0
$M
               temp      rain  humidity
temp      2.0863437 -1.169590 0.7765303
rain     -1.1695903  2.676182 2.0880970
humidity  0.7765303  2.088097 3.4682640

1/(2*(model[[3]]@covar@range.val)^2)





$B
, , temp

         a    b    c    d
a 10.61014 0.00 0.00 0.00
b  0.00000 5.89 0.00 0.00
c  0.00000 0.00 1.16 0.00
d  0.00000 0.00 0.00 6.63

, , rain

     a    b        c    d
a 4.79 0.00  0.00000 0.00
b 0.00 4.74  0.00000 0.00
c 0.00 0.00 11.38632 0.00
d 0.00 0.00  0.00000 9.83

, , humidity

     a    b    c        d
a 4.79 0.00 0.00  0.00000
b 0.00 4.74 0.00  0.00000
c 0.00 0.00 9.56  0.00000
d 0.00 0.00 0.00 11.81837

###################################################
vv=function (x1, x2 = x1, hp, ...)
{
    stopifnot(is.mdm(x1))
    stopifnot(is.mdm(x2))
    stopifnot(compatible(x1, x2))
    stopifnot(compatible(x1, hp))
    t1 <- types(x1)
    t2 <- types(x2)
    Sigma <- matrix(NA, nrow(xold(x1)), nrow(xold(x2)))
    rownames(Sigma) <- rownames(x1)
    colnames(Sigma) <- rownames(x2)
    B <- B(hp)
    for (i in levels(t1)) {
        for (j in levels(t2)) {
            ni <- which(i == levels(t1))
            nj <- which(j == levels(t2))
            ii <- which(i == t1)
            jj <- which(j == t2)
            if (FALSE) {
                print(i)
                print(j)
                print("-------")
            }
            Sigma[ii, jj] <- corr.matrix(xold(x2)[jj, , drop = FALSE],
                xold(x1)[ii, , drop = FALSE], pos.def.matrix = solve(solve(B[,
                  , i])/2 + solve(B[, , j])/2), ...)
        }
    }
    return(Sigma)


alpha=1/(sqrt(2*diag(Bo[[i]])))
sig=1/(sqrt(2*Sigma1))
ola=model[[1]]
ola@covar@sd2=1
ola@covar@range.val=1/sqrt(2*diag(Bo[[i]]))
corr.matrix(xold=ola@X,scales=diag(bb))
covMatrix(ola@covar,ola@X)
bb=Bo[[1]]

