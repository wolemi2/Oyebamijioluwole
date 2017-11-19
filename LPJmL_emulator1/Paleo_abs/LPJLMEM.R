#########################################################
## THE MAIN EMULATOR CODE
### DO NOT ALTER ANY THING BELOW THIS LINE
#########################################################

i=1:3
rcp="other"
fertilization=0;man=1   
#########################
new_dir=getwd()
source("configuration_option.R",echo=F)

#j=12-j#######reverse the decadal data
#j=11-j#######reverse the decadal data



if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:11
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed/climgen5a.R",echo=TRUE)
} else {
source("irrigated/climgen5b.R",echo=TRUE)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
