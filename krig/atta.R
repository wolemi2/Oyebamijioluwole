##################
#cpr=setClass("cpr", 		
         #representation(
                    #    d = "integer",            	## (spatial) dimension
                     #   name = "character",             ## "powexp"
                      #  paramset.n = "integer",         ## number of parameters sets 
                        ##   gauss, exp : 1;  powexp : 2
                      #  var.names = "character",  	## e.g.  c("Lat", "Long") length d
			## s.d. of dor the non-nugget part of error
                      #  sd2 = "numeric",       		## variance (stationarity)
			## nugget part
                      #  known.covparam = "character",   ## known covariance parameters (except nugget): "All" or "Known"
                        #nugget.flag = "logical",  	## logical : is there a nugget effect ?
                        #nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
                        #nugget = "numeric",    		## nugget (variance)
  			## total number of parameters (except sigma and nugget)
                      #  param.n = "integer",            ## range.n + shape.n
  			## range part 
                       # range.n = "integer",            ## number of distinct range parms
                       # range.names = "character",	## their name (usually "theta")
                       # range.val = "numeric",          ## their values
  			## shape part, if any 
                        #shape.n = "integer",            ## number of distinct shape parms
                        #shape.names = "character",      ## their name ("p", "nu", "alpha", etc.)
                        #shape.val = "numeric"           ## their values
                        #),
         validity = function(obj) {
           
           covset <- c("gauss", "exp", "matern3_2", "matern5_2", "powexp")
           if (!is.element(obj@name, covset)) {
             cat("The list of available covariance functions is:\n", covset, "\n")
             return("invalid character string for 'covtype' argument")
           }
           
           if (!identical(obj@sd2, numeric(0))) {
             if (obj@sd2 < 0) {
               return("The model variance should be non negative")
             }
           }
           
           if (length(obj@nugget) > 1L) {
             return("'nugget' must be a single non-negative number. For heteroskedastic noisy observations, use 'noise.var' instead.")
           }
           
           if (!identical(obj@nugget, numeric(0))) {
             if (obj@nugget < 0) {
               return("The nugget effect should be non negative")
             }
           }
           
           if (!identical(obj@range.val, numeric(0))) {
             if (min(obj@range.val) < 0) {
               return("The range parameters must have positive values")
             }
             if (length(obj@range.val) != obj@d) {
               return("Incorrect number of range parameters")
             }
           }
           
           if (!identical(obj@shape.val, numeric(0))) {
             if (min(obj@shape.val) < 0) {
               return("The shape parameters must have positive values")
             }
             if (length(obj@shape.val) != obj@d) {
               return("Incorrect number of shape parameters")
             }
             if (identical(obj@name, "powexp") && (max(obj@shape.val) > 2)) {
               return("The exponents must be <= 2 for a Power-Exponential covariance")
             }
           }
           TRUE
         }
         #)

#cpr2=new("cpr")
setMethod("covMatrix", 
          signature = "cpr", 
          definition = function(obj, X, noise.var=NULL) {
            covMatrix.f(obj=obj, X=X, noise.var=noise.var)
          }
)


setMethod("covMat1Mat2", 
          signature = "cpr", 
          definition = function(obj, X1, X2, nugget.flag=FALSE) {
            covMat1Mat2.f(obj=obj, X1=X1, X2=X2, nugget.flag=nugget.flag)
          }
)


setMethod("vec2",signature = "cpr",definition = vec2.f)

setMethod("vec1",signature = "cpr",definition = vec1.f)


setMethod("bbound", 
          signature = "cpr", 
          definition = function(obj, X){
            if (obj@paramset.n==1) {
              k <- obj@range.n
              lower <- rep(1e-10, k)
              upper <- 2 * diff(apply(X, 2, range))
              upper <- as.vector(upper)
            } else if (identical(obj@name, "powexp")) {              
              # coef. order : theta_1, ..., theta_d, p_1, ..., p_d
              lower <- rep(1e-10, obj@param.n)
              upper <- 2 * diff(apply(X, 2, range))
              k <- obj@shape.n
              upper <- as.vector(c(upper, rep(2, k)))
            } else stop("No default values for covariance parameters bounds, the inputs 'lower' and 'upper' are required")
            return(list(lower=lower, upper=upper))
          }
)


setMethod("cov.deri",
          signature = "cpr", 
          definition = function(obj, X, C0, k) {
            cov.deri.f(obj=obj, X=X, C0=C0, k=k)
          }
)

######################################################
