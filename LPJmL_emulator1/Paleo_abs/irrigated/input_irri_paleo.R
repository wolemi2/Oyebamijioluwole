##########climgen relative to baseline (2001-2010) average
climgen_path=path
##############################################NEW
rcp="RCP3PD"
if(rcp=="other") {
 data1<- list.files(path=climgen_path,pattern=".nc")
} else {
data1<- list.files(path=climgen_path,pattern=sub("rcp",rcp,"rcp"))}

output="output"
cld=paste(climgen_path,data1[[1]],sep="/")
v3=paste(output,data1[[1]],sep="_")
v3=sub("output_cld_pat_CMIP3",output,paste(v3))
pre=paste(climgen_path,data1[[2]],sep="/")
tmp=paste(climgen_path,data1[[3]],sep="/")
wet=paste(climgen_path,data1[[4]],sep="/")
#############################
new_co2 <- list.files(path=climgen_path,pattern = ".txt")
####################### 
n=j
npixel=59199
load("irrigated/k3")
library(ncdf)
#load("irrigated/N")
load("irrigated/grid_out")
###################################################

#########################################
library(ncdf)
grid_len=npixel
cld=open.ncdf(cld)
kenny=cld[8]$dim$MONTH$vals
pre=open.ncdf(pre)
tmp=open.ncdf(tmp)
wet=open.ncdf(wet)
npixel=59199
w=sort(c(seq(1,120,12),seq(2,120,12),seq(12,120,12)))
s=sort(c(seq(6,120,12),seq(7,120,12),seq(8,120,12)))
sp=sort(c(seq(3,120,12),seq(4,120,12),seq(5,120,12)))
au=sort(c(seq(9,120,12),seq(10,120,12),seq(11,120,12)))

library(doParallel)
cl <- makeCluster(length(n))
registerDoParallel(cl)
tt<-{foreach(z=1:length(n),.combine =cbind) %do%{
T=matrix(0,nrow=59199,ncol=length(n))
T=k3+((z-1)*201600)}}

ind=c(tt)
my_list=list(cld,pre,tmp,wet)
load("irrigated/k3")

#####################################
m=1
h=c("CLD_MONTHLY","PRECIP_MONTHLY","TMP_MONTHLY","WET_MONTHLY")
skip=12*(2001-kenny[1])+1
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(1),-1,-1)))
e=array(olu,dim=c(120,length(1),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(1),201600,4))
X2a=aperm(X2,c(2,1,3))
X2b=array(X2a,c(201600*length(1),4))
X2c=X2b[ind[1:59199],]
if(j==0) {
  X2d=rbind(X2c,X2c,X2c,X2c,X2c,X2c,X2c,X2c)
} else {
  X2d=X2c
}
inp1=X2d
###########
m=2
h=c("CLD_MONTHLY","PRECIP_MONTHLY","TMP_MONTHLY","WET_MONTHLY")
skip=12*(2001-kenny[1])+1
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(1),-1,-1)))
e=array(olu,dim=c(120,length(1),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(1),201600,4))
X2a=aperm(X2,c(2,1,3))
X2b=array(X2a,c(201600*length(1),4))
X2c=X2b[ind[1:59199],]
if(j==0) {
  X2d=rbind(X2c,X2c,X2c,X2c,X2c,X2c,X2c,X2c)
} else {
  X2d=X2c
}

inp2=X2d
###########
m=3
h=c("CLD_MONTHLY","PRECIP_MONTHLY","TMP_MONTHLY","WET_MONTHLY")
skip=12*(2001-kenny[1])+1
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(1),-1,-1)))
e=array(olu,dim=c(120,length(1),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(1),201600,4))
X2a=aperm(X2,c(2,1,3))
X2b=array(X2a,c(201600*length(1),4))
X2c=X2b[ind[1:59199],]
if(j==0) {
  X2d=rbind(X2c,X2c,X2c,X2c,X2c,X2c,X2c,X2c)
} else {
  X2d=X2c
}

inp3=X2d
#########
m=4
h=c("CLD_MONTHLY","PRECIP_MONTHLY","TMP_MONTHLY","WET_MONTHLY")
skip=12*(2001-kenny[1])+1
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(1),-1,-1)))
e=array(olu,dim=c(120,length(1),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(1),201600,4))
X2a=aperm(X2,c(2,1,3))
X2b=array(X2a,c(201600*length(1),4))
X2c=X2b[ind[1:59199],]
if(j==0) {
  X2d=rbind(X2c,X2c,X2c,X2c,X2c,X2c,X2c,X2c)
} else {
  X2d=X2c
}

inp4=X2d
datt=cbind(inp1,inp2,inp3)######initial
rm(h,m,my_list)
##############################################################
source("irrigated/paleo_input.R")
dat=(dat_paleo+datt)
DAT=cbind(dat,wdays)

#dat=dat[,c(1:4,9:12,17:20,25:28,5:8,13:16,21:24,29:32)]
load("irrigated/other_inp")
soil=other_inp[[1]][1:59199,3]
lat.lpj=other_inp[[1]][1:59199,4]

###########################################################END