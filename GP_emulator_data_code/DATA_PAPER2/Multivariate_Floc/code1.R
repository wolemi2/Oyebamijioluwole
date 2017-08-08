###########Multivariate emulator for floc data

library(MASS)
library(Matrix)
library(emulator)
library(multivator)
library(abind)
#set.seed(1)
setwd("./Multivariate_Floc")
dirname=getwd()
##load the scaled inputs between[0,1]
X=read.csv("X.csv")## matrix of inputs
Y=read.csv("Y.csv")##load the matrix of logarithms of output data
X=as.matrix(X)
Y=as.matrix(Y)
id=c(9,28:31,33:36)
X=X[,id]
nam=c("MumHET","s","o2","no2","no3","y1","y2","y3","y4")
###########################################EMULATOR FITTING
## model matrix
#form= ~1 + y4 + y3 + y2 + y1 + no3 + no2 + o2 + s + MumHET 
form=~1+MumHET+s+o2+no2+no3+y1+y2+y3+y4
H<-model.matrix(form,data=as.data.frame(X))
## Define number of columns of output and model matrices
k<-dim(Y)[2]
n<-dim(Y)[1]
m<-dim(H)[2]
ninp=dim(X)[2]
id2=ninp+1
############## Function giving the log marginal likelihood
tau=1
likfun2<-function(theta){
R=corr.matrix(X,scales=theta[-id2])#corr matrix for output#2
A=R+(theta[id2]*diag(tau,dim(R)))##add random nugget/noise
if(any(Re(eigen(A,TRUE,TRUE)$values) < 0)){
A=as.matrix(nearPD(A)$mat)
} else {
A=A
}
Omega=quad.form.inv(A,H)#t(H)%*%Ainv%*%H
dA=(determinant(A)$modulus[1])###return log ofdeterminant
dOmega=(determinant(Omega)$modulus[1])##return log of determinant
sig=(quad.form.inv(A,Y) - quad.form(quad.form.inv(quad.form.inv(A,H),t(solve(A,H))), Y))
S=determinant(sig)$modulus[1]
-0.5*k*dA-0.5*k*dOmega-.5*(n-m)*S
}
#########initialize the parameter values
theta=c(0.7287402,0.0000100,0.0000100,1.9737827,0.0000100,7.1322561,4.5044932,0.6049189,1.2098281,0.1294618)
initial_scale=theta#rep(.2,(ninp+1))#theta
lower=rep(1e-5,length(initial_scale))
upper=2*(apply(X,2,max)-apply(X,2,min))
opt<-optim(fn=likfun2,par=initial_scale,method="L-BFGS-B",control=list(fnscale= -1,trace=3),lower=lower,upper=upper)
save(opt,file=paste(getwd(),"opt",sep="/"))

R=corr.matrix(X,scales=opt$par[-id2])#corr matrix for output#2
A=R+(opt$par[id2]*diag(tau,dim(R)))
if(class(try(chol(A), silent = T)) == "try-error") # positive definite
{A=as.matrix(nearPD(A)$mat)
Ainv=chol2inv(chol(A))###
} else {
Ainv=chol2inv(chol(A))###
}
model=list(X=X,Y=Y,theta=opt$par,form=form,Ainv=Ainv,nam=nam)
dirname=getwd()
save(model,file=paste(dirname,"model",sep="/"))###spa
##########################PREDICTION FUNCTIONS
pred=function(model,newdata){	
X<- model$X
Y<- model$Y
theta1=model$theta[id2]
theta<- model$theta[-id2]
form<- model$form
Ainv<- model$Ainv
A=solve(Ainv)
n2<- dim(newdata)[1]
n<- dim(X)[1]
X0<- newdata
colnames(X0)=model$nam
######model matrix 
H<-model.matrix(form,data=as.data.frame(X))
H0<-model.matrix(form,as.data.frame(X0))
A01=corr.matrix(X0,X,scales=theta)###cross correlation
A00=corr.matrix(X0,X0,scales=theta)##test point correlation
A00=A00
Omega=quad.form.inv(A,H)#t(H)%*%Ainv%*%H
betahat=solve(Omega,crossprod(H, solve(A, Y)))##betahat=ginv(t(H)%*%Ainv%*%H)%*%(t(H)%*%Ainv%*%Y)
mu_star=(H0%*%betahat)+crossprod(A01,Ainv)%*%(Y-H%*%betahat)#H0%*%betahat+t(A01)%*%Ainv%*%(Y-H%*%betahat)
c_star=A00-(t(A01)%*%Ainv%*%A01)+(((H0-(t(A01)%*%Ainv%*%H))%*%solve(Omega)%*%t(H0-(t(A01)%*%Ainv%*%H))))
Sigma=(t(Y-H%*%betahat)%*%Ainv%*%(Y-H%*%betahat))/(n-m)
#######Predictive variance per output
if(class(try(chol(c_star), silent = T)) == "try-error") # positive definite
{c_star=as.matrix(nearPD(c_star)$mat)
} else {
c_star=c_star###
}
svar=list()
for(i in 1:n2){
svar[[i]]=(diag(c_star)[i]*diag(Sigma))
}
K0=abind(svar,along=2)
out=list(mu=mu_star,K=t(K0))
return(out)
}
##########Validation results
nout=4
nres=20##number of validation results
load("model")###R object
##load the matrix of logarithms of the 20 validation datasets
val=list()
out=list()
dens=list()
for(k in 1:nres){
val[[k]]=read.csv(paste(sub("k",k,"Y_k.csv"),sep="/"))
}
ini=as.matrix(read.csv("XX_ini.csv"))##matrix of initial parameters
#j=1###time
for(k in 1:nres){
j=1
inp=as.matrix(read.csv(sub("k",k,"XX_k.csv")))###matrix of input parameters
input=matrix(c(inp[j,],ini[k,]),nrow=1)
input=input[,id,drop=FALSE]
colnames(input)=nam=c("MumHET","s","o2","no2","no3","y1","y2","y3","y4")
##iterative use of emulator (dynamically)
nsim=200##number of Monte carlo simulation
npoints=1#number of points to predict per time point
mstar=list();result=list();mstar2=list()
sims=matrix(NA,nrow=npoints,ncol=nsim)
ntime=dim(val[[k]])[1]
dens0=array(NA,c(ntime,nsim,nout))

m1=pred(model,newdata=input)
mu=m1$mu###predictions
varr=m1$K#variance
mstar=mvrnorm(nsim,c(mu),diag(c(varr),ncol=nout,nrow=nout))
V1=varr
V2=0
mstar2[[j]]=mstar
result[[j]]=cbind(mu,V1+V2)
ntime=dim(val[[k]])[1]
dens0[j,,]=mstar

for(j in 2:ntime){
g2a=mstar2[[j-1]]
inp2=matrix(rep(inp[j,],nsim),nrow=nsim,byrow=TRUE)
input=cbind(inp2,g2a)
input=input[,id,drop=FALSE]
colnames(input)=nam
m1=pred(model,newdata=input)
p1=array(m1$mu,c(npoints,nsim,nout))
p2=array(m1$K,c(npoints,nsim,nout))
mu=t(matrix(apply(p1,3,rowMeans)))
dens0[j,,]=p1[1,,]
mu2=array(rep(mu,nsim),c(npoints,nout,nsim));mu2=aperm(mu2,c(1,3,2))
varr=apply(p2,c(1,3),mean)#expectation of variance
#varr2=apply((p1-mu2)^2,3,rowMeans)#variance of expectation
varr2=diag(var(p1[1,,]))/nsim#variance of expectation
V1=varr;V2=varr2
mstar=mvrnorm(nsim,c(mu),diag(c(V1+V2),ncol=nout,nrow=nout))
mstar2[[j]]=mstar
result[[j]]=cbind(mu,V1+V2)
}

res=abind(result,along=1)
emu=res[,1:4]
vv=res[,5:8]
out[[k]]=list(emu,vv)
dens[[k]]=dens0
}
#Compute the multivariate Proportion of variance (P) and Root Mean squared Error (RMSE) Table 1
vv2=list()
for(k in 1:20){
vv2[[k]]=out[[k]][[1]]###emulator
}
lammp=abind(val,along=1)
emu=abind(vv2,along=1)
P=rep(NA,nout)
RMSE=P
for(i in 1:nout){
P[i]=1-(sum((lammp[,i]-emu[,i])^2)/sum((lammp[,i]- mean(lammp[,1]))^2))##propotion
RMSE[i]=sqrt(mean((lammp[,i]-emu[,i])^2))##RMSE
}
P
RMSE
save(val,file="val")#lammps
save(out,file="out")##emu +
save(dens,file="dens")##
###############################################################PLOTS
##FIGURE 8 IN THE MANUSCRIPTS
#k=10
for(k in 1:20){
emu2=out[[k]][[1]]###emulator
vv=out[[k]][[2]]###emulator variance
lammp=val[[k]]####LAMMPS simulated values
ntime2=nrow(lammp)
dev=sqrt(vv)
###################
tit=c(expression("ln(Floc eqv. diameter)  metre"),expression("ln(Fractal dimension) "),"ln(Number of particles)","ln(Total mass)  kg")
pow=c(1,1,1,1)
dirn3b=sub("k",k,paste("Density_revision","k",sep="_"))
dir.create(dirn3b)
dirn=paste(getwd(),dirn3b,sep="/")

pos=c("topleft","topright","topleft","topleft")
for(i in 1:(nout)){
dirname=file.path(dirn,paste(paste(dirn3b,"_"),i,".jpeg",sep=""))
jpeg(file=dirname,quality=100)

bw=c(.5,0.5,.5,.5)
band=bw[i]
dat=dens[[k]]
dd=array(NA,c(512,nsim,2))
for(r in 1:nsim){
#gap=density(dat[j,,i],bw=band)
gap=density(dat[,r,i],bw=band)
dd[,r,1]=gap$x
dd[,r,2]=gap$y
#band=gap$bw
}
emu=apply(dd,c(1,3),median)
ssd=apply(dd,c(1,3),sd)
#
lammps=(lammp[,i]*pow[i])
d1=density(lammps,bw=band)
h0=emu[,1]###x
h00=emu[,2]
h1=emu[,2]+(2*ssd[,2])
h2=emu[,2]-(2*ssd[,2])
x1=min(h0,d1$x)
x2 = max(h0,d1$x)
y1 = min(h1,h2,d1$y)
y2 = max(h1,h2,d1$y)
plot(d1,type="l",main=tit[i],xlab="",cex.lab=1.5,cex.axis=1.6,cex.main=1.7,ylim=c(y1,y2)+c(0,0),xlim=c(x1,x2))
lines(d1$x,h00,col="green")
lines(d1$x,h1, col="red",lwd=3,lty=2)
lines(d1$x,h2, col="red",lwd=3,lty=2)
legend(pos[i],c("Predicted","Observed","95% C.I"),
fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
}

##FIGURE 7 IN THE MANUSCRIPTS
#k=10
emu2=out[[k]][[1]]###emulator
vv=out[[k]][[2]]###emulator variance
lammp=val[[k]]####LAMMPS simulated values
ntime2=nrow(lammp)
dev=sqrt(vv)
dirn2=sub("k",k,paste("Linear_revision","k",sep="_"))
#tit=c(expression("ln(floc eqv. diameter)   metre"),"ln(fractal dimension)","ln(number of particle)",expression("ln(floc total mass)  kg"))
dir.create(dirn2)
dirn=paste(getwd(),dirn2,sep="/")
for(i in 1:(nout)){
dirname=file.path(dirn,paste(paste(dirn2,"_"),i,".jpeg",sep=""))
jpeg(file=dirname,quality=100)
mytitle=tit[i]
lammps=lammp[,i]*pow[i]
emu=emu2[,i]*pow[i]
sd=dev[,i]*pow[i]

U=emu+(2*sd)
L=emu-(2*sd)
#time=(1:ntime2)*10^4/(24*60*60)###days
time=(1:ntime2)*30*275/(24*60*60)###days
x1 = min(U,lammps,emu,L)
x2 = max(U,lammps,emu,L)
par(mar=c(5.1,5.4,4.1,2.1))
plot(time,lammps,main=mytitle,ylab=expression("Output"),lty=1,col="black",xlab="time (day)",cex.lab=1.7,cex.axis=1.7,cex.main=1.7,xlim=range(time),ylim=c(x1,x2)+c(0,.03))
lines(time,L, col="red",lwd=3,lty=2)
lines(time,U, col="red",lwd=3,lty=2)
lines(time,emu,col="green",lwd=3)
legend(1,(x2+.03),c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
}
}
