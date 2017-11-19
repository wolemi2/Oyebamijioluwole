##Configuration options for the LPJML PALEO emulator
# Select "j" which represent time slice from options (j=c(0,1,2,3,4,5,6,7,8,9,10,11)
# Select the "output resolution" type  by uncomment the line
# select "crop types" from options (irrigated or rainfed)
# NOTE: j==0 implies output all the 11 decades
#Requred: ClimGEN netcdf climate data to run the emulator in the format below 
##############indicate path to ClimGEN climate data
path=".\\anthropocene"

j=2  ###choose time point(1:11); (j=0 for ALL output ONCE)


###########################choose output resolution
#resolution="LOW"
resolution="HIGH"

##########################choose crop types
#type="irrigated"
type="rainfed"
#########################################################END