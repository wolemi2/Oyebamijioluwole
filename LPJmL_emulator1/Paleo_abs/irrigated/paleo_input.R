###############Compute paleo inputs
#path2= "C:\\Users\\olu\\Dropbox\\LPJmL_emulator1\\Paleo_abs\\anthropocene2"
path2=".\\anthropocene2"
##############################################NEW
#rcp=”RCP3PD”
data2= list.files(path=path2,pattern="GENIE")

output="output"
cld2=paste(path2,data2[[1]],sep="/")
v3=paste(output,data2[[1]],sep="_")
v3=sub("output_cld_pat_CMIP3",output,paste(v3))
pre2=paste(path2,data2[[2]],sep="/")
tmp2=paste(path2,data2[[3]],sep="/")
#wet2=paste(path2,data2[[2]],sep="/")

w=sort(c(seq(1,12,12),seq(2,12,12),seq(12,12,12)))
s=sort(c(seq(6,12,12),seq(7,12,12),seq(8,12,12)))
sp=sort(c(seq(3,12,12),seq(4,12,12),seq(5,12,12)))
au=sort(c(seq(9,12,12),seq(10,12,12),seq(11,12,12)))


################
load("anthropocene/aye2b")
library("rhdf5")
my_list=list(cld2,pre2,tmp2)
h=c("CLD_MONTHLY","PRECIP_MONTHLY","TMP_MONTHLY")
#olu=c(h5read(my_list[[m]],name=h[m],start=c(skip,1,1),count=c(12*length(1),-1,-1)))
m=1
d=j-1

if(d==0){
dat_paleo=matrix(0,nrow=59199,ncol=12)
wdays=inp4
} else {

library(doParallel)
cl <- makeCluster(12)
registerDoParallel(cl)
datan<-{foreach(i=((12*d)-11):(12*d),.combine =rbind) %do%{
h5read(cld2,name=h[m],index=list(i,1:280,1:720))
}}
datan=c(datan)
gg=which(datan>9.96921e+35)
datan[gg]=NaN
e=array(datan,c(12,280*720))
X1=cbind(c(colMeans(e[s,])),c(colMeans(e[w,])),c(colMeans(e[sp,])),c(colMeans(e[au,])))
X2c=X1[ind,]
vinp1=X2c
###
m=2
datan<-{foreach(i=((12*d)-11):(12*d),.combine =rbind) %do%{
h5read(pre2,name=h[m],index=list(i,1:280,1:720))
}}
datan=c(datan)
gg=which(datan>9.96921e+35)
datan[gg]=NaN
e=array(datan,c(12,280*720))
wwet=e
X1=cbind(c(colMeans(e[s,])),c(colMeans(e[w,])),c(colMeans(e[sp,])),c(colMeans(e[au,])))
X2c=X1[ind,]
vinp2=X2c
###
m=3
datan<-{foreach(i=((12*d)-11):(12*d),.combine =rbind) %do%{
h5read(tmp2,name=h[m],index=list(i,1:280,1:720))
}}
datan=c(datan)
gg=which(datan>9.96921e+35)
datan[gg]=NaN
e=array(datan,c(12,280*720))
X1=cbind(c(colMeans(e[s,])),c(colMeans(e[w,])),c(colMeans(e[sp,])),c(colMeans(e[au,])))
X2c=X1[ind,]
vinp3=X2c
#####
m=4
aye2d=cbind(c(rowMeans(aye2b[,s])),c(rowMeans(aye2b[,w])),c(rowMeans(aye2b[,sp])),c(rowMeans(aye2b[,au])))

vinp4b=(inp2+vinp2)
vinp4b[which(vinp4b<0)]=0###set negative abs precip to zero

eee=round((aye2d*vinp4b*.01)^.45)##convert absol prep to wetdays
eee[which(eee>31)]=31
wdays=eee
dat_paleo=cbind(vinp1,vinp2,vinp3)#####paleo_change##############

}
##########################END
