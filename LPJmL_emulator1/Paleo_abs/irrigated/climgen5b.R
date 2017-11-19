
################################
if(j==0) {

j=1
source("irrigated/input_irri_paleo.R",echo=TRUE)
if(fertilization==0){
fert=3:4
} else {
if(fertilization==1) {
fert=1:2
} else {
fert=1:4
}}

npixel=59199
load("irrigated/k3")
library(ncdf)
#load("irrigated/ini")
load("irrigated/crop_N")
load("irrigated/grid_out")
load("irrigated/crop_form4")
load("irrigated/coef4")
m=i
###################################################
nam=c("scld" ,   "wcld"  ,  "spcld" ,  "acld" ,"spre",   
 "wpre" ,   "sppre" ,  "apre" ,   "stmp" ,   "wtmp" ,   "sptmp"  , "atmp",   
"swet" ,   "wwet" ,   "spwet" ,  "awet" ,"soil", "lat.lpj","LA") 
################################
w=1
LA=rep(man,npixel)
input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)

#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)


prediction1=prediction

################
#####include PCA analysis
#source("irrigated/wls_irri.R",echo=F)

################################
j=2
w=2
source("irrigated/input_irri_paleo.R",echo=TRUE)
input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)

#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction2=prediction

################
#####include PCA analysis
#source("irrigated/wls_irri.R",echo=F)
################################
j=3
w=3
source("irrigated/input_irri_paleo.R",echo=TRUE)

input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)

#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction3=prediction


################
#####include PCA analysis
################################
j=4
w=4
source("irrigated/input_irri_paleo.R",echo=TRUE)

input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)

#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction4=prediction


################
#####include PCA analysis
################################
j=5
w=5

source("irrigated/input_irri_paleo.R",echo=TRUE)
input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)

#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction5=prediction

################
#####include PCA analysis
################################
j=6
w=6
source("irrigated/input_irri_paleo.R",echo=TRUE)
input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)

#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction6=prediction

################
#####include PCA analysis
################################
j=7
w=7
source("irrigated/input_irri_paleo.R",echo=TRUE)

source("irrigated/input_irri_paleo.R",echo=TRUE)
input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)

#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction7=prediction

################
#####include PCA analysis
################################
j=8
w=8
source("irrigated/input_irri_paleo.R",echo=TRUE)

source("irrigated/input_irri_paleo.R",echo=TRUE)
input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)
##########################
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction8=prediction
#########################################
j=9
w=9
source("irrigated/input_irri_paleo.R",echo=TRUE)
source("irrigated/input_irri_paleo.R",echo=TRUE)
input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction9=prediction

################
#####include PCA analysis
#source("irrigated/wls_irri.R",echo=F)


j=10
w=10
source("irrigated/input_irri_paleo.R",echo=TRUE)

input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction10=prediction


j=11
w=11
source("irrigated/input_irri_paleo.R",echo=TRUE)

input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction11=prediction



#################
if(resolution=="LOW"){
crop=cbind(prediction1,prediction2,prediction3,prediction4,prediction5,prediction6,prediction7,prediction8,prediction9,prediction10,prediction11)
crop=as.matrix(crop)
crop=array(crop,dim=c(npixel,2*ncol(crop1),10))
j=0
source("irrigated/aggregate.R",echo=F)
} else {
crop=cbind(prediction1,prediction2,prediction3,prediction4,prediction5,prediction6,prediction7,prediction8,prediction9,prediction10,prediction11)
wal=cbind(prediction1,prediction2,prediction3,prediction4,prediction5,prediction6,prediction7,prediction8,prediction9,prediction10,prediction11)

crop=as.matrix(crop)
crop=array(crop,dim=c(npixel,3,11))

j=0
source("irrigated/net.R",echo=F)}
###############################END
} else {
if(fertilization==0){
fert=3:4
} else {
if(fertilization==1) {
fert=1:2
} else {
fert=1:4
}}
j=j
source("irrigated/input_irri_paleo.R",echo=TRUE)
npixel=59199
#load("irrigated/k3")
library(ncdf)
#load("irrigated/ini")
load("irrigated/crop_N")
load("irrigated/grid_out")
load("irrigated/crop_form4")
load("irrigated/coef4")
m=i
###################################################
nam=c("scld" ,   "wcld"  ,  "spcld" ,  "acld" ,"spre",   
 "wpre" ,   "sppre" ,  "apre" ,   "stmp" ,   "wtmp" ,   "sptmp"  , "atmp",   
"swet" ,   "wwet" ,   "spwet" ,  "awet" ,"soil", "lat.lpj","LA") 
w=1
LA=rep(man,npixel)
input=cbind(DAT,soil,lat.lpj,LA)
colnames(input)=nam
input=as.data.frame(input)

################cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
#ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
ras=crop_form4[[m]]
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.matrix(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction1=prediction

#################################################
if(resolution=="LOW"){
crop=prediction1
crop=as.matrix(crop)
crop=array(crop,dim=c(59199,ncol(crop),1))
source("irrigated/aggregate.R",echo=F)
} else {
crop=prediction1
crop=as.matrix(crop)
crop=array(crop,dim=c(59199,ncol(crop),1))
source("irrigated/net.R",echo=F)
}}