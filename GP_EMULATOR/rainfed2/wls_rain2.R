##############################wls
#source("rainfed2/MM.R")
if(out=="current"){
load("rainfed2/masked/crop_err4c")
load("rainfed2/masked/N1")
T=N1
} else {
load("rainfed2/crop_N")
load("rainfed2/unmasked/crop_err4c")
T=crop_N
}
##########
v4=sub("man",man,(sub("j",j,"p_j_man")))
load(paste("rainfed2/pred",v4,sep="/"))
crop_pred4b=eval(parse(text = as.name(v4)))

library(mlegp)
m=i
w=j##############correction
###NEW CROP ANALYSIS REL TO 1ST DEC" 2nd STAGE EMUL
##########
m=1
s1=array(crop_pred4b[[m]],c(59199,4*5))
l=186
e1=crop_err4c[[m]]
######
#E=(e1[,w,,man])
#E=array(E,c(186,4*5))
#Y2=array(s1,c(59199,4*5))
#m1=prcomp(Y2[T[[m]],],scale=TRUE,center=TRUE)
#l_cereal=m1$rotation ###for 84 by 84
#y_star=prediction[T[[m]],m]
#Gamma=Y2[T[[m]],]%*%l_cereal ###for
#x_i=t(Gamma[,1:4])%*%Y2[T[[m]],]
#x_star=t(Gamma[,1:4])%*%y_star
#X=t(x_i)
#X_star=t(x_star)

E=(e1[,w,,man])
E=array(E,c(186,4*5))
Y2=array(s1,c(59199,4*5))
m1=prcomp(t(Y2[T[[m]],]),scale=TRUE,center=TRUE)
l_cereal=m1$rotation
##ada=(Y2[T[[m]],]-m1$center)/matrix(rep(m1$scale,20),ncol=20)
X=m1$x[,1:4]#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=predict(m1,t(as.matrix(y_star)))[,1:4]

###############fit GP emulator
residual=t(E)
ind=which(residual[m,]!=0) 
u=which(duplicated(X))
library(snowfall)
sfInit(parallel=TRUE,cpus=4)
sfLibrary(mlegp)
if(length(u)==0){
fitPC1 =mlegp(X=X,Z=residual[,ind],constantMean=1,parallel=TRUE)
}else{
fitPC1 =mlegp(X=X[-u,],Z=residual[-u,ind],constantMean=1,parallel=TRUE)}
sfStop()
id=length(ind)
pred=rep(0,id)
for(r in 1:id){
pred[r]=predict(fitPC1[[r]],X_star)
}
predY=rep(0,186)
predY[ind]=pred
E_crop1=predY
#######################################
m=2
s1=array(crop_pred4b[[m]],c(59199,4*5))
l=186
e1=crop_err4c[[m]]
######
E=(e1[,w,,man])
E=array(E,c(186,4*5))
Y2=array(s1,c(59199,4*5))
m1=prcomp(t(Y2[T[[m]],]),scale=TRUE,center=TRUE)
l_cereal=m1$rotation
##ada=(Y2[T[[m]],]-m1$center)/matrix(rep(m1$scale,20),ncol=20)
X=m1$x[,1:4]#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=predict(m1,t(as.matrix(y_star)))[,1:4]
###############fit GP emulator
residual=t(E)
ind=which(residual[m,]!=0) 
u=which(duplicated(X))
library(snowfall)
sfInit(parallel=TRUE,cpus=4)
sfLibrary(mlegp)
if(length(u)==0){
fitPC2 =mlegp(X=X,Z=residual[,ind],constantMean=1,parallel=TRUE)
}else{
fitPC2 =mlegp(X=X[-u,],Z=residual[-u,ind],constantMean=1,parallel=TRUE)}
sfStop()
id=length(ind)
pred=rep(0,id)
for(r in 1:id){
pred[r]=predict(fitPC2[[r]],X_star)
}
predY=rep(0,186)
predY[ind]=pred
E_crop2=predY
#######################################
m=3
s1=array(crop_pred4b[[m]],c(59199,4*5))
e1=crop_err4c[[m]]
######
E=(e1[,w,,man])
E=array(E,c(186,4*5))
Y2=array(s1,c(59199,4*5))
m1=prcomp(t(Y2[T[[m]],]),scale=TRUE,center=TRUE)
l_cereal=m1$rotation
##ada=(Y2[T[[m]],]-m1$center)/matrix(rep(m1$scale,20),ncol=20)
X=m1$x[,1:4]#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=predict(m1,t(as.matrix(y_star)))[,1:4]
##t(l_cereal[,1:4])%*%(as.matrix(y_star-m1$center)/m1$scale)
###############fit GP emulator
residual=t(E)
ind=which(residual[m,]!=0) 
u=which(duplicated(X))
library(snowfall)
sfInit(parallel=TRUE,cpus=4)
sfLibrary(mlegp)
if(length(u)==0){
fitPC3 =mlegp(X=X,Z=residual[,ind],constantMean=1,parallel=TRUE)
}else{
fitPC3 =mlegp(X=X[-u,],Z=residual[-u,ind],constantMean=1,parallel=TRUE)}
sfStop()
id=length(ind)
pred=rep(0,id)
for(r in 1:id){
pred[r]=predict(fitPC3[[r]],X_star)
}
predY=rep(0,186)
predY[ind]=pred
E_crop3=predY
#######################################
m=4
s1=array(crop_pred4b[[m]],c(59199,4*5))
e1=crop_err4c[[m]]
######
E=(e1[,w,,man])
E=array(E,c(186,4*5))
Y2=array(s1,c(59199,4*5))
m1=prcomp(t(Y2[T[[m]],]),scale=TRUE,center=TRUE)
l_cereal=m1$rotation
##ada=(Y2[T[[m]],]-m1$center)/matrix(rep(m1$scale,20),ncol=20)
X=m1$x[,1:4]#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=predict(m1,t(as.matrix(y_star)))[,1:4]
##t(l_cereal[,1:4])%*%(as.matrix(y_star-m1$center)/m1$scale)
###############fit GP emulator
residual=t(E)
ind=which(residual[m,]!=0) 
u=which(duplicated(X))
library(snowfall)
sfInit(parallel=TRUE,cpus=4)
sfLibrary(mlegp)
if(length(u)==0){
fitPC4 = mlegp(X=X[,1:4],Z=residual[,ind],constantMean=1,parallel=TRUE)
}else{
fitPC4 =mlegp(X=X[-u,],Z=residual[-u,ind],constantMean=1,parallel=TRUE)}
sfStop()
id=length(ind)
pred=rep(0,id)
for(r in 1:id){
pred[r]=predict(fitPC4[[r]],X_star)
}
predY=rep(0,186)
predY[ind]=pred
E_crop4=predY
#######################################
#E_crop2=E_crop2*(-1)
#E_crop3=E_crop3*(-1)

#E_crop4=E_crop4*(-1)
E_crop=list(E_crop1,E_crop2,E_crop3,E_crop4)
#detach("package:MASS", unload=TRUE)
source("rainfed2/aggregate3.R",echo=FALSE)
#rm(E_crop1,E_crop2,E_crop3,E_crop4,crop_err4c,X,x_star,Gamma,x_i,Y2,y_star,predY,fitPC,residual,m1,E,e1,l_cereal)
