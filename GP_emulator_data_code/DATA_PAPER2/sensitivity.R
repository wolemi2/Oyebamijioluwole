setwd("C:/Users/user/Desktop/review/")

load("C:\\Users\\user\\Desktop\\review\\Univariate_Floc\\sens")
nam=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YHET","YAOB","YNOB","YEPS","Y1","Do2","Dnh4","Dno2","Dno3","Ds","diffT","s","o2","no2","no3","nh4")
library(sensitivity)
dat=rbind(sens[[1]]$S$mean,sens[[2]]$S$mean,sens[[3]]$S$mean,sens[[4]]$S$mean)
gap1=dat
dat1=dat[,1:32]
rm(sens)
######biofilm
load("C:\\Users\\user\\Documents\\Univariate_Biofilms/sens")
dat=rbind(sens[[1]]$S$mean,sens[[2]]$S$mean,sens[[3]]$S$mean,sens[[4]]$S$mean)
gap2=dat
dat2=dat[,1:32]
dat2[,32]=0######nh4 for biofilms
dat=rbind(dat1,dat2)
dat[which(dat<0)]=0
par(mar=c(5.1,4.5,1.9,.5)+0.0)
nam4=c("Floc equivalent diameter","Floc fractal dimension","Floc total number of particles","Floc total mass","Biofilm average height","Biofilm surface roughness","Biofilm segregation indices","Biofilm species diversity")
xvals=barplot((dat),ylab="Sensitivity indices",main="Sensitivity analysis",offset=0,cex.lab=1.6,col=c("green","blue","pink","red","grey","yellow","blueviolet","brown"),density=c(40,60,90,150,180,235,270,300),cex.main=1.8,space=0.3,angle=c(45,90,135,180,225),cex.axis=1.6,cex.names=1.6,args.legend=list(x=c(3,8),y=c(2.5,3.5),bty="n",cex=2.2),legend=paste(nam4))
text(xvals,par("usr")[3]-.01,srt=45,adj=1,labels=nam,xpd=TRUE,cex=1.6)
###

dat1[which(dat1<0)]=0
dat2[which(dat2<0)]=0
r1=rowSums(dat1)
r2=rowSums(dat2)
ss=array(NA,dim(dat1))
for(i in 1:nrow(dat1)){
ss[i,]=dat1[i,]/r1[i]
}
ss2=array(NA,dim(dat2))
for(i in 1:nrow(dat2)){
ss2[i,]=dat2[i,]/r2[i]
}
sss=rbind(ss,ss2)
sss[which(is.na(sss))]=0
par(mar=c(5.1,4.5,1.9,.5)+0.0)
nam4=c("Floc equivalent diameter","Floc fractal dimension","Floc total number of particle","Floc total mass","Biofilm average height","Biofilm surface roughness","Biofilm segregation indices","Biofilm species diversity")
xvals=barplot((sss),ylab="Sensitivity indices",main="Sensitivity analysis",offset=0,cex.lab=1.6,col=c("green","blue","pink","red","grey","yellow","blueviolet","brown"),density=c(40,60,90,150,180,235,270,300),cex.main=1.8,space=0.3,angle=c(45,90,135,180,225),cex.axis=1.6,cex.names=1.6,args.legend=list(x=c(3,8),y=c(1.5,3),bty="n",cex=2.2),legend=paste(nam4))
text(xvals,par("usr")[3]-.01,srt=45,adj=1,labels=nam,xpd=TRUE,cex=1.6)
