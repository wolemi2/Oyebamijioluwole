Title: Gaussian process emulation of an individual-based model simulation of microbial communities
O.K. Oyebamiji, D.J. Wilkinson, P.G. Jayathilake, T.P. Curtis, S.P. Rushton, B. Li, P. Gupta

This folder contains data and code to reproduce results in the manuscripts. There are 4 subfolders namely:
(i) Multivariate_Floc - This subfolder contains the data and code for fitting the multivariate emulator for the floc and generating Figures 7 & 8 and (/rho and RMSE multivariate portion  Table 1) in the manuscript.
(a) code1.R R functions for fitting emulator and performing crossvalidation for multivariate floc emulation
(b) X.csv - A 2000 X 12 matrix of sampled input parameters (floc training data)
(c) Y.csv - A 2000 X 4 matrix of sampled characterized outputs  (floc training data)
(d) XX.csv - A 20 X 8 matrix of input parameters for 20 validation points (floc test data)
(e) XX_ini.csv - A 20 X 4 matrix of initial input parameters for 20 validation points (floc test data)
(f) (Y_1.csv -Y_20.csv) -  Each is a matrix of 4 simulated outputs of varying dimension for 20 different validation points (floc test data). 

(ii) Multivariate_Biofilms - This folder contain the data and code for fitting the
 multivariate emulator for the biofilms and generating Figures 9 & 10  and (/rho and RMSE multivariate portion Table 1) in the manuscript.
(a) code3.R - R functions for fitting emulator and performing crossvalidation for mutivariate biofilms emulation
(b) X.csv - A 1500 X 12 matrix of sampled input parameters (biofilms training data)
(c) Y.csv - A 1500 X 4 matrix of sampled characterized outputs  (biofilms training data)
(d) XX.csv - A 20 X 8 matrix of input parameters for 20 validation points (biofilms test data)
(e) XX_ini.csv - A 20 X 4 matrix of initial input parameters for 20 validation points (biofilms test data)
(f) (Y_1.csv -Y_20.csv) -  Each is a matrix of 4 simulated outputs of varying dimension for 20 different validation points (biofilms  test data). 

(iii) Univariate_Floc - This folder contain the data and code for fitting the univariate emulator for the floc and (/rho and RMSE of univariate portion Table 1) in the manuscript.
(a) code2.R - R functions for fitting univariate emulator, performing sensitivity analysis and crossvalidation of univariate floc emulation and plotting Histograms in Figure 4
(b) X.csv - A 2500 X 36 matrix of sampled input parameters ie (32 parameters +4 initial values)  (floc training data)
(c) Y.csv - A 2500 X 4 matrix of sampled characterized outputs  (floc training data)
(d) XX_i.csv - A n X 32 matrix of input parameters for 20 validation points (i=1:20, floc test data)
(e) Y_i.csv - A n X 4 matrix of output for 20 validation points (i =1:20, floc test data)
(f) All_Floc.csv - A 15487 X 4 matrix containing all the 4 simulated floc outputs data before sampling
(g) XX_ini.csv - A 20 X 4 matrix of initial input parameters for 20 validation points (floc test data)

(iv) Univariate_Biofilms - This folder contains the data and code for fitting the univariate emulator for the biofilms
(a) code4.R - R functions for fitting univariate emulator, performing sensitivity analysis and crossvalidation of univariate biofilms emulation and plotting Histograms in Figure 5
(b) X.csv - A 2500 X 36 matrix of sampled input parameters ie (32 parameters +4 initial values)  (biofilms training data)
(c) Y.csv - A 2500 X 4 matrix of sampled characterized outputs  (biofilms training data)
(d) XX_i.csv - A n X 32 matrix of input parameters for 20 validation points (i =1:20, biofilms test data)
(e) Y_i.csv - A n X 4 matrix of output for 20 validation points (i =1:20, biofilms test data)
(f) All_Biofilms.csv - A 5320 X 4 matrix containing all the 4 simulated biofilms outputs data before sampling.
(g) XX_ini.csv - A 20 X 4 matrix of initial input parameters for 20 validation points (biofilms test data)

Others:
Sensitivity.R - R codes for generating the plot of sensitivity indices in Figure 11.
NOTE: Figures 1, 2, 3 and 6 in the manuscripts are made outside R.
###################################################################
Dr. O.K. Oyebamiji
School of Mathematics & Statistics
Newcastle University
Newcastle upon Tyne
NE1 7RU
UK
E-mail: wolemi2@yahoo.com; Oluwole.Oyebamiji@ncl.ac.uk