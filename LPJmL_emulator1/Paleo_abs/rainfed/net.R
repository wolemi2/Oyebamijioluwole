############netcdf output
crp=aperm(crop,c(1,3,2))
load("rainfed/k3")
if(j==0) {
   w=1:11
} else {
   w=j
}
load("rainfed/lonlat")
xvals <- Lon
yvals <- Lat
zvals=w
nz=length(zvals)
ny <- length(yvals)
nx <- length(xvals)
new_dir=getwd()
dirname <- paste(getwd(),"output_result",sep="/")
setwd(dirname)

xdim <- dim.def.ncdf("LONGITUDE","degrees",xvals )
ydim <- dim.def.ncdf("LATITUDE", "degrees",yvals )
tdim <- dim.def.ncdf("TIME","decade",zvals)
#---------
# Make var
#---------
mv <- 1.e30 # missing value
var_cereal <- var.def.ncdf("cereal","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal cereal")
var_rice <- var.def.ncdf("rice","gC/m^2",dim=list(xdim,ydim,tdim), mv,longname="mean decadal rice")
var_maize <- var.def.ncdf("maize","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal maize")
#var_oil_max <- var.def.ncdf("oil_max","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadaloil_max")
#var_grd <- var.def.ncdf("grd","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal groundnut")


#var_ini_cereal <- var.def.ncdf("ini_cereal","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean baseline ini_cereal")
#var_ini_rice <- var.def.ncdf("ini_rice","gC/m^2",dim=list(xdim,ydim,tdim), mv,longname="mean baseline ini_rice")
#var_ini_maize <- var.def.ncdf("ini_maize","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean baseline ini_maize")
#var_ini_oil_max <- var.def.ncdf("ini_oil_max","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean baseline ini_oil_max")
#var_ini_grd <- var.def.ncdf("ini_grd","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean baseline groundnut")



#---------------------
# Make new output file
#---------------------
crop=create.ncdf(paste(as.name(paste("crop",v3,sep="_")),"nc",sep="."),list(var_cereal,var_rice,var_maize))

dat <- array(NA,dim=c(nx*ny,nz))
dat[k3,]=crp[,,1]  
dat=array(dat,dim=c(280,720,length(w)))
dat=aperm(dat,c(2,1,3))

dat2 <- array(NA,dim=c(nx*ny,nz))
dat2[k3,]=crp[,,2]  
dat2=array(dat2,dim=c(280,720,length(w)))
dat2=aperm(dat2,c(2,1,3))

dat3 <- array(NA,dim=c(nx*ny,nz))
dat3[k3,]=crp[,,3]  
dat3=array(dat3,dim=c(280,720,length(w)))
dat3=aperm(dat3,c(2,1,3))


put.var.ncdf(crop,var_cereal,dat, start=c(1,1,1), count=c(nx,ny,nz))
put.var.ncdf(crop,var_rice,dat2, start=c(1,1,1), count=c(nx,ny,nz))
put.var.ncdf(crop,var_maize,dat3, start=c(1,1,1), count=c(nx,ny,nz))

close.ncdf(crop)
setwd(new_dir)
