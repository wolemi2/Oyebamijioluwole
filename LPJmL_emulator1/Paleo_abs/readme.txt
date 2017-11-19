LPJmL CROP YIELD EMULATOR "LPJLMEM"
DESCRIPTION
This version of LPJmL crop emulator produces crop yield change (gC/M^2) and its baseline values for 3 different crops: Temperate cereal, rice and  maize, rainfed and irrigated. All the outputs are decadal values. The emulator can provide its outputs either as a “LOW” or “HIGH” resolution. The HIGH resolution is a gridded scale of 0.5 by 0.5 degrees while the LOW resolution provides the results at Country level. The emulator is run from ClimGEN paleo data. The output format will be an Excel "csv" data for the LOW resolution and Netcdf file for the HIGH resolution output. 

All required data and R scripts for the emulator are in the directories rainfed and irrigated subfolder. All the outputs from running the emulator go into the subfolder “output_results”. The main R source code is "LPJMLEM.R". There are some options to select and set in the configuration file before the emulator can be run to produce the desire results.

CONFIGURATION FILE
This emulator can run for various configuration options in order to change the default options. There are 3 basic options namely (time-step, output resolutions and crop type. Any desire output result can be set in the configuration file "configuration_option" before running the main code. 

Major options to set in the configuration file:
(i) Select "j" from options (j=c(1,2,3,4,5,6,7,8,9,10,11); decadal levels

(ii) Select the "resolution" by uncomment the line; this represents output format (LOW/HIGH) resolution.

(iii)Select "type" by uncomment the line; this option represents crop types, it can either be a rainfed or irrigated crop.

RUN THE CODE
The main source code is "LPJMLEM.R". The emulator can be run from either the Window or Linux platform. 
 (A) Put the climate data in the subfolder “CLIMGEN”. These are precipitation, temperature, wetdays frequency and cloud cover.
 (B) Set the configuration file.
 (C) Run the emulator: Require ncdf package
(i) To run from Window machine:
 At R console………
> setwd(“path to/EMU”)
> source(“LPJMLEM.R”,echo=TRUE)

(ii) To run from Linux terminal:
$ cd path to/EMU
$ R CMD BATCH LPJMLEM.R result.txt 

 (D) Obtain the results from subfolder “output_results”.

OUTPUT RESULTS
For each run of the emulator, the Netcdf or CSV output data produces decadal yield. For instance, any typical CSV output data for 1 time point will have 11 columns of data. COLUMN 1== Country; COLUMNS 2:6== change in yield for cereal, rice and maize.
COLUMN 7:11 are the baseline yield average(2005-2014) for ini_cereal, ini_rice, ini_maize, Ini_oil and ini_grd respectively.
Also, for option j=0, there will be 11X8 columns accordingly for each decade.
The Netcdf are self-explanatory for each corresponding results.

NOTE: This emulator will produce output correspond to the following time-steps for combinations of any option set above for different "j"
j=1 change in crop yield between (2005-2014) and (2015-2024)
j=2  ............................(2005-2014) and (2025-2034)
j=3  ............................(2005-2014) and (2035-2044)
j=4  ............................(2005-2014) and (2045-2054)
j=5  ............................(2005-2014) and (2055-2064)
j=6  ............................(2005-2014) and (2065-2074)
j=7  ............................(2005-2014) and (2075-2084)
j=8  ............................(2005-2014) and (2085-2094)
j==0, ALL 8 decade results at once.


CONTACT: oluwole.oyebamiji@open.ac.uk or wolemi2@yahoo.com for further details.
