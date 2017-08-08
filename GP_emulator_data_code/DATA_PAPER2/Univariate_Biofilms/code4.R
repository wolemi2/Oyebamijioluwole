###########Univariate emulator for biofilms data
library(doParallel)
library(abind)
setwd("./Univariate_Biofilms")
dirname=getwd()
nout=4
#########################################################################
##Histograms in Figure 4
bb3=read.csv("All_Biofilms.csv")##read in all Biofilms data
bb3=log(bb3)
tit=c(expression("ln(Average height)  metre"),expression("ln(Surface roughness)  metre"),"ln(Segregation indices)","ln(Species diversity indices)")
dirn3="Histogram"
for(i in 1:nout){
dirname=file.path(paste(dirn3,i,".jpeg",sep=""))
jpeg(file=dirname,quality=100)
hist(bb3[,i],main="",xlab=tit[i],cex.main=1.5,cex.lab=1.5,cex.axis=1.5,col="pink")
dev.off()
}
nam=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YHET","YAOB","YNOB","YEPS","Y1","Do2","Dnh4","Dno2","Dno3","Ds","diffT","s","o2","no2","no3","nh4","y1","y2","y3","y4")
###########################################EMULATOR FITTING
f=list(
~1 + y1 + no3 + no2 + o2 + s + diffT + Ds + Dno3 + Dno2 + 
    Dnh4 + Do2 + Y1 + YEPS + YNOB + YAOB + YHET + bEPS + bNOB + 
    bAOB + bHET + etaHET + MumNOB + MumAOB + MumHET + Ko2NOB + 
    Kno2NOB + Ko2AOB + Knh4AOB + Kno3HET + Kno2HET + Ko2HET + 
    KsHET
,
~1 + y2 + no3 + no2 + o2 + s + diffT + Ds + Dno3 + Dno2 + 
    Dnh4 + Do2 + Y1 + YEPS + YNOB + YAOB + YHET + bEPS + bNOB + 
    bAOB + bHET + etaHET + MumNOB + MumAOB + MumHET + Ko2NOB + 
    Kno2NOB + Ko2AOB + Knh4AOB + Kno3HET + Kno2HET + Ko2HET + 
    KsHET
,
~1 + y3 + no3 + no2 + o2 + s + diffT + Ds + Dno3 + Dno2 + 
    Dnh4 + Do2 + Y1 + YEPS + YNOB + YAOB + YHET + bEPS + bNOB + 
    bAOB + bHET + etaHET + MumNOB + MumAOB + MumHET + Ko2NOB + 
    Kno2NOB + Ko2AOB + Knh4AOB + Kno3HET + Kno2HET + Ko2HET + 
    KsHET
,
~1 + y4 + no3 + no2 + o2 + s + diffT + Ds + Dno3 + Dno2 + 
    Dnh4 + Do2 + Y1 + YEPS + YNOB + YAOB + YHET + bEPS + bNOB + 
    bAOB + bHET + etaHET + MumNOB + MumAOB + MumHET + Ko2NOB + 
    Kno2NOB + Ko2AOB + Knh4AOB + Kno3HET + Kno2HET + Ko2HET + 
    KsHET)

#####
Y=read.csv("Y.csv")##read matrix of training data
inp=read.csv("X.csv")##read matrix of input data
X=list()
for(i in 1:nout){
vv=inp[,c(1:32,(32+i))]
X[[i]]=vv
}
nug=c(0.0001,0.0001,.0001,0.0001)
modd=foreach(i=1:nout,.verbose=FALSE,.errorhandling="stop") %dopar% {
library(DiceKriging)
km(formula=f[[i]],X[[i]][,-32],response=Y[,i],covtype="exp",nugget=nug[i])
}
dirname=getwd()
save(modd,file=paste(dirname,"modd",sep="/"))###spa
n <- 1000
nvar=27+5;nout=4
#nam=c(para,nam2,nam3)
nam=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YHET","YAOB","YNOB","YEPS","Y1","Do2","Dnh4","Dno2","Dno3","Ds","diffT","s","o2","no2","no3","nh4","y1","y2","y3","y4")
X1 <- data.frame(matrix(runif((nvar) * n), nrow = n))
X2 <- data.frame(matrix(runif((nvar) * n), nrow = n))
library(sensitivity)
sens=list()
for(i in 1:nout){
sens[[i]]=sobolGP(modd[[i]],X1=X1,X2=X2,type="UK")
}
save(sens,file=paste(dirname,"sens",sep="/"))###spa
#############################################CROSS-VALIDATION
##################
dirname=getwd()
nout=4
nres=20##number of validation results
load("modd")
##load the matrix of logarithms of the 20 validation datasets
val=list()
out=list()
setwd("../")
for(k in 1:nres){
#val[[k]]=read.csv(paste(sub("k",k,"Univariate_Biofilms/Y_k.csv"),sep="/"))
val[[k]]=read.csv(paste(sub("k",k,paste("Univariate_Biofilms","Y_k.csv",sep="/")),sep="/"))
}
ini=as.matrix(read.csv("Univariate_Biofilms/XX_ini.csv"))##matrix of initial parameters
setwd(dirname)

X=list()
for(k in 1:nres){
inp=as.matrix(read.csv(sub("k",k,"XX_k.csv")))###matrix of input parameters
inp=inp[1:19,]

##iterative use of emulator (dynamically)
j=1
input=matrix(c(inp[j,],ini[k,]),nrow=1)
colnames(input)=nam

nsim=1000##number of Monte carlo simulation
npoints=1#number of points to predict per time point
mstar=matrix(NA,nrow=nsim,ncol=nout);result=list();mstar2=list()
sims=matrix(NA,nrow=npoints,ncol=nsim)
mu=varr=rep(NA,nout)
for(i in 1:nout){
vv=matrix(input[,c(1:32,(32+i))],nrow=1)
colnames(vv)=nam[c(1:32,(32+i))]
X[[i]]=vv
m1=predict(modd[[i]],newdata=X[[i]][,-32,drop=FALSE],type="UK",se.compute=TRUE)
mu[i]=m1$mean###predictions
varr[i]=(m1$sd)^2#variance
mstar[,i]=rnorm(nsim,mu[i],sqrt(varr[i]))
}
V1=varr
V2=0
mstar2[[j]]=mstar
result[[j]]=c(mu,V1+V2)
p1=p2=matrix(NA,nrow=nsim,ncol=nout)
ntime=dim(val[[k]])[1]
for(j in 2:ntime){
g1=matrix(rep(t(inp[j,]),nsim),nrow=npoints*nsim,ncol=32,byrow=TRUE)
g2a=mstar2[[j-1]]
input=cbind(g1,g2a)
colnames(input)=nam
for(i in 1:nout){
vv=matrix(input[,c(1:32,(32+i))],nrow=nsim)
colnames(vv)=nam[c(1:32,(32+i))]
X[[i]]=vv
m1=predict(modd[[i]],newdata=X[[i]][,-32,drop=FALSE],type="UK",se.compute=TRUE)
p1[,i]=m1$mean###predictions
p2[,i]=(m1$sd)^2#variance
}
mu=apply(p1,2,mean)
mu2=t(array(rep(mu,nsim),c(nout,nsim)))
varr=apply(p2,2,mean)#expectation of variance
varr2=apply((p1-mu2)^2,2,mean)#variance of expectation
V1=varr;V2=varr2
V=sqrt(V1+V2)
for(i in 1:nout){
mstar[,i]=rnorm(nsim,mu[i],V[i])
}
mstar2[[j]]=mstar
result[[j]]=c(mu,V1+V2)
}
res=abind(result,along=0)
emu=res[,1:4]
vv=res[,5:8]
out[[k]]=list(emu,vv)
}
#Compute the univariate Proportion of variance (P) and Root Mean squared Error (RMSE) Table 1
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

