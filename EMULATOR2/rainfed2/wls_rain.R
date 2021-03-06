############################wls
load("rainfed2/crop_N")
load("rainfed2/crop_err4b")
load("rainfed2/crop_pred4b")
library(MASS)
m=i
w=j##############correction
###NEW CROP ANALYSIS REL TO 1ST DEC" 2nd STAGE EMUL
#######################################
m=1
l=length(crop_N[[m]])
e1=array(crop_err4b[[m]],c(length(crop_N[[m]]),8,16,7))
s1=array(crop_pred4b[[m]],c(length(crop_N[[m]]),8,16,7))
##########PCA
E=(e1[,w,,man])
Y2=s1[,w,,man]
#E=matrix(E,nrow=length(crop_N[[m]]))
#Y2=matrix(Y2,nrow=length(crop_N[[m]]))
m1=prcomp(Y2,scale=F,center=F)
l_cereal=m1$rotation ###for 16 by 16
y_star=prediction[crop_N[[m]],m]
Gamma=Y2%*%l_cereal ###for 16648 by 16648
x_i=t(Gamma[,1:4])%*%Y2
x_star=t(Gamma[,1:4])%*%y_star
lambda=(m1$sdev)^2/sum((m1$sdev)^2)
#lambda=(m1$sdev)^2
d=round(sqrt(colSums(lambda*(matrix(x_star,nrow=4,ncol=16)-x_i)^2)))
#d=sort(d)
w_i=(1/d^2)/sum((1/d^2))
#w_i=(1/d^2)+.5
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop1=E[,ind] 
}else {E_crop1=rowMeans(E[,ind])}} else {
X=t(x_i)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%x_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop1=pmax(EE,(r1))}

#########
m=2
l=length(crop_N[[m]])
e1=array(crop_err4b[[m]],c(length(crop_N[[m]]),8,16,7))
s1=array(crop_pred4b[[m]],c(length(crop_N[[m]]),8,16,7))
##########PCA
E=(e1[,w,,man])
Y2=s1[,w,,man]
#E=matrix(E,nrow=length(crop_N[[m]]))
#Y2=matrix(Y2,nrow=length(crop_N[[m]]))
m1=prcomp(Y2,scale=F,center=F)
l_cereal=m1$rotation ###for 16 by 16
y_star=prediction[crop_N[[m]],m]
Gamma=Y2%*%l_cereal ###for 16648 by 16648
x_i=t(Gamma[,1:4])%*%Y2
x_star=t(Gamma[,1:4])%*%y_star
lambda=(m1$sdev)^2/sum((m1$sdev)^2)
d=round(sqrt(colSums(lambda*(matrix(x_star,nrow=4,ncol=16)-x_i)^2)))
w_i=(1/d^2)/sum((1/d^2))
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop2=E[,ind] 
}else {E_crop2=rowMeans(E[,ind])}
} else {
X=t(x_i)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%x_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop2=pmax(EE,(r1))}
###################
m=3
l=length(crop_N[[m]])
e1=array(crop_err4b[[m]],c(length(crop_N[[m]]),8,16,7))
s1=array(crop_pred4b[[m]],c(length(crop_N[[m]]),8,16,7))
##########PCA
E=(e1[,w,,man])
Y2=s1[,w,,man]
#E=matrix(E,nrow=length(crop_N[[m]]))
#Y2=matrix(Y2,nrow=length(crop_N[[m]]))
m1=prcomp(Y2,scale=F,center=F)
l_cereal=m1$rotation ###for 16 by 16
y_star=prediction[crop_N[[m]],m]
Gamma=Y2%*%l_cereal ###for 16648 by 16648
x_i=t(Gamma[,1:4])%*%Y2
x_star=t(Gamma[,1:4])%*%y_star
lambda=(m1$sdev)^2/sum((m1$sdev)^2)
d=round(sqrt(colSums(lambda*(matrix(x_star,nrow=4,ncol=16)-x_i)^2)))
w_i=(1/d^2)/sum((1/d^2))
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop3=E[,ind] 
}else {E_crop3=rowMeans(E[,ind])}
} else {
X=t(x_i)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%x_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop3=pmax(EE,(r1))}
############
m=4
l=length(crop_N[[m]])
e1=array(crop_err4b[[m]],c(length(crop_N[[m]]),8,16,7))
s1=array(crop_pred4b[[m]],c(length(crop_N[[m]]),8,16,7))
##########PCA
E=(e1[,w,,man])
Y2=s1[,w,,man]
E=matrix(E,nrow=length(crop_N[[m]]))
Y2=matrix(Y2,nrow=length(crop_N[[m]]))
m1=prcomp(Y2,scale=F,center=F)
l_cereal=m1$rotation ###for 16 by 16
y_star=prediction[crop_N[[m]],m]
Gamma=Y2%*%l_cereal ###for 16648 by 16648
x_i=t(Gamma[,1:4])%*%Y2
x_star=t(Gamma[,1:4])%*%y_star
lambda=(m1$sdev)^2/sum((m1$sdev)^2)
d=round(sqrt(colSums(lambda*(matrix(x_star,nrow=4,ncol=16)-x_i)^2)))
w_i=(1/d^2)/sum((1/d^2))
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop4=E[,ind] 
}else {E_crop4=rowMeans(E[,ind])}
} else {
X=t(x_i)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%x_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop4=pmax(EE,(r1))}
E_crop=list(E_crop1,E_crop2,E_crop3,E_crop4)
#rm(E_crop1,E_crop2,E_crop3,E_crop4,crop_err4b,crop_pred4b,beta,E_star,X,x_star,W,Gamma,x_i,Y2,y_star)

