require("ape")
require("rexpokit")
require("cladoRcpp")


#######################################################
# LaplacesDemon functions for Bayesian analysis with BioGeoBEARS
#######################################################


#######################################################
# plotlp
#######################################################
#' Plot LaplacesDemon plots after excluding some burnin
#'
#' Takes the objects output by a \code{\link[LaplacesDemon]{LaplacesDemon}} MCMC search, and plots LnP, parameter values, etc.
#'
#' @param Fit The \code{\link[LaplacesDemon]{LaplacesDemon}} output object.
#' @param MyData The \code{\link[LaplacesDemon]{LaplacesDemon}} input data.
#' @param PDF Plot to a PDF? (Recommended!)
#' @param Parms The parameters
#' @param burnfract The burnin fraction
#' @param ... Additional arguments to standard functions
#' @return The plot
#' @export
#' @seealso \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
plotlp <- function(Fit, MyData, PDF=TRUE, Parms=NULL, burnfract=0.0, ...)
	{
	defaults='
	PDF=TRUE
	Parms=NULL
	burnfract=0.1
	'
	
	
	if (burnfract > 0)
		{
		# Edit the burnin for the plot
		newFit = change_burnin(Fit, MyData, burnfract=0.1, newsamps=NULL)
		} else {
		newFit = Fit
		}

	plot_demonoid2(x=newFit, Data=MyData, PDF=PDF, Parms=Parms)
	
	}




# Re-parameterization
# Transform continuous-sampled parameter on -Inf, +Inf to
# a boundedly-sampled parameter
# 

#######################################################
# transform_with_logistic
#######################################################
#' Transform continuously-sampled parameter on -Inf, +Inf to a boundedly-sampled parameter
#'
#' \code{\link[LaplacesDemon]{LaplacesDemon}} likes to run its MCMC sampling on a simple number line.  Thus, the likelihood
#' function etc. should transform the numbers into the range desired, e.g. 0-1.
#'
#' @param paramval The input value in the range [-Inf, +Inf]
#' @param minval The minimum rescaled value (default 0)
#' @param maxval The maximum rescaled value (default 1)
#' @return \code{paramval}
#' @export
#' @seealso \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite LaplacesDemon_Tutorial
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
#' transform_with_logistic(paramval=-10, minval=0, maxval=5)
#' # NOTE that when minval=0, maxval=1, you have logit!
#' 
#' # Transform regular-space parameters to logistic space
#' invlogit(0.5)
#' transform_with_logistic(paramval=0.5, minval=0, maxval=1)
#' transform_with_logistic(paramval=0.5, minval=0, maxval=5)
#' 
#' invlogit(0.3)
#' transform_with_logistic(paramval=0.3, minval=0, maxval=1)
#' transform_with_logistic(paramval=0.3, minval=0, maxval=5)
#' 
#' # Transform logistic-space parameters back to regular space
#' logit(0.5744425)
#' invtransform_with_logistic(transformed_paramval=0.5744425, minval=0, maxval=1)
#' invtransform_with_logistic(transformed_paramval=3.112297, minval=0, maxval=5)
#' 
#' logit(0.5744425)
#' invtransform_with_logistic(transformed_paramval=0.5744425, minval=0, maxval=1)
#' invtransform_with_logistic(transformed_paramval=2.872213, minval=0, maxval=5)
#' 
#' # These should transform, then undo the transform
#' invtransform_with_logistic(transformed_paramval=transform_with_logistic(
#' paramval=0.3, minval=0, maxval=1), minval=0, maxval=1)
#' invtransform_with_logistic(transformed_paramval=transform_with_logistic(
#' paramval=0.3, minval=0, maxval=5), minval=0, maxval=5)
#' 
transform_with_logistic <- function(paramval, minval=0, maxval=1)
	{
	transformed_paramval = minval + maxval * (exp(paramval) / (exp(paramval) + 1))
	
	return(transformed_paramval)
	}





# 
# # Transform regular-space parameters to logistic space
# invlogit(0.5)
# transform_with_logistic(paramval=0.5, minval=0, maxval=1)
# transform_with_logistic(paramval=0.5, minval=0, maxval=5)
# 
# invlogit(0.3)
# transform_with_logistic(paramval=0.3, minval=0, maxval=1)
# transform_with_logistic(paramval=0.3, minval=0, maxval=5)
# 
# # Transform logistic-space parameters back to regular space
# logit(0.5744425)
# invtransform_with_logistic(transformed_paramval=0.5744425, minval=0, maxval=1)
# invtransform_with_logistic(transformed_paramval=3.112297, minval=0, maxval=5)
# 
# logit(0.5744425)
# invtransform_with_logistic(transformed_paramval=0.5744425, minval=0, maxval=1)
# invtransform_with_logistic(transformed_paramval=2.872213, minval=0, maxval=5)
# 
# # These should transform, then undo the transform
# invtransform_with_logistic(transformed_paramval=transform_with_logistic(
#' paramval=0.3, minval=0, maxval=1), minval=0, maxval=1)
# invtransform_with_logistic(transformed_paramval=transform_with_logistic(
#' paramval=0.3, minval=0, maxval=5), minval=0, maxval=5)
#
#######################################################
# invtransform_with_logistic
#######################################################
#' Transform rescaled values back to continuously-sampled parameter on -Inf, +Inf to a boundedly-sampled parameter
#'
#' \code{\link[LaplacesDemon]{LaplacesDemon}} likes to run its MCMC sampling on a simple number line.  Thus, the likelihood
#' function etc. should transform the numbers into the range desired, e.g. 0-1. invtransform_with_logistic does the reverse transform.
#'
#' @param transformed_paramval The input value in the range [minval, maxval]
#' @param minval The minimum rescaled value (default 0)
#' @param maxval The maximum rescaled value (default 1)
#' @return \code{paramval}
#' @export
#' @seealso \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite LaplacesDemon_Tutorial
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
#' # Transform regular-space parameters to logistic space
#' invlogit(0.5)
#' transform_with_logistic(paramval=0.5, minval=0, maxval=1)
#' transform_with_logistic(paramval=0.5, minval=0, maxval=5)
#' 
#' invlogit(0.3)
#' transform_with_logistic(paramval=0.3, minval=0, maxval=1)
#' transform_with_logistic(paramval=0.3, minval=0, maxval=5)
#' 
#' # Transform logistic-space parameters back to regular space
#' logit(0.5744425)
#' invtransform_with_logistic(transformed_paramval=0.5744425, minval=0, maxval=1)
#' invtransform_with_logistic(transformed_paramval=3.112297, minval=0, maxval=5)
#' 
#' logit(0.5744425)
#' invtransform_with_logistic(transformed_paramval=0.5744425, minval=0, maxval=1)
#' invtransform_with_logistic(transformed_paramval=2.872213, minval=0, maxval=5)
#' 
#' # These should transform, then undo the transform
#' invtransform_with_logistic(transformed_paramval=transform_with_logistic(
#' paramval=0.3, minval=0, maxval=1), minval=0, maxval=1)
#' invtransform_with_logistic(transformed_paramval=transform_with_logistic(
#' paramval=0.3, minval=0, maxval=5), minval=0, maxval=5)
#' 
invtransform_with_logistic <- function(transformed_paramval, minval=0, maxval=1)
	{
	# transformed_paramval = minval + maxval * (exp(paramval) / (exp(paramval) + 1))
	# (paramval-minval)/maxval = (exp(transformed_paramval) / (exp(transformed_paramval) + 1))
	# log((paramval-minval)/maxval) = log(exp(transformed_paramval)) - log((exp(transformed_paramval) + 1))
	# log((paramval-minval)/maxval) = transformed_paramval - (transformed_paramval - 1)
	# log((paramval-minval)/maxval) = transformed_paramval

	# Error check (as in logit())
    if ({
        any(transformed_paramval < minval)
    } || {
        any(transformed_paramval > maxval)
    }) 
        stop("transformed_paramval must be in [minval, maxval].")

	
	# Value must be between 0 and 1
	val01 = (transformed_paramval-minval)/maxval
	
	paramval = log(val01 / (1-val01))
	
	return(paramval)
	}





#######################################################
# reparam_LapDem_output
#######################################################
#' Re-parameterize LaplacesDemon MCMC output
#'
#' \code{\link[LaplacesDemon]{LaplacesDemon}} likes to run its MCMC sampling on a simple number line.  Thus, the likelihood
#' function etc. should transform the numbers into the range desired, e.g. 0-1.
#'
#' This function transforms the \code{\link[LaplacesDemon]{LaplacesDemon}} output
#'
#' @param Fit The \code{\link[LaplacesDemon]{LaplacesDemon}} output object.
#' @param MyData The \code{\link[LaplacesDemon]{LaplacesDemon}} input data.
#' @param transfun The function to use for the transformation.
#' @return \code{Fit} The transformed MCMC output.
#' @export
#' @seealso \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite LaplacesDemon_Tutorial
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
reparam_LapDem_output <- function(Fit, MyData, transfun=transform_with_logistic)
	{
	require(LaplacesDemon)
	
	# Safety
	origFit = Fit
	
	# Get the data
	BioGeoBEARS_model_object = MyData$BioGeoBEARS_model_object
	params_table = BioGeoBEARS_model_object@params_table
	freeTF = params_table$type == "free"

	
	

	for (i in 1:ncol(Fit$Posterior1))
		{
		minval = params_table$min[freeTF][i]
		maxval = params_table$max[freeTF][i]
		

		#######################################################
		# Posterior1
		#######################################################
		# Mark the 0.0s
		zerosTF = Fit$Posterior1[, i] == 0.0
		
		# Change the rest to regularspace
		Fit$Posterior1[, i] = transfun(paramval=Fit$Posterior1[, i], minval, maxval)
		
		# Convert 0.0s back to 0.0
		Fit$Posterior1[zerosTF, i] = 0.0


		#######################################################
		# Posterior2
		#######################################################
		if ( !is.null(nrow(Fit$Posterior2)) && (nrow(Fit$Posterior2) >= 1)) 
			{
			# Mark the 0.0s
			zerosTF = Fit$Posterior2[, i] == 0.0
			
			# Change the rest to regularspace
			Fit$Posterior2[, i] = transfun(paramval=Fit$Posterior2[, i], minval, maxval)
			
			# Convert 0.0s back to 0.0
			Fit$Posterior2[zerosTF, i] = 0.0
			} else {
			
			# If it's 1 row...
			for (j in 1:length(Fit$Posterior2))
				{
				zerosTF = Fit$Posterior2[i] == 0.0
				
				# Change the rest to regularspace
				Fit$Posterior2[i] = transfun(paramval=Fit$Posterior2[i], minval, maxval)
				
				# Convert 0.0s back to 0.0
				if (zerosTF == TRUE)
					{
					Fit$Posterior2[i] = 0.0
					}
				
				}
			
			}
	}


	#######################################################
	# Re-do the summary, from the parameter-space posterior sample
	#######################################################
	
	
	#######################################################
	# re-do the summary, copying from LaplacesDemon()
	#######################################################
	
	# Posterior1 is the thinned samples, including everything
	LIV = length(Fit$Initial.Values)
	thinned = Fit$Posterior1
	thinned.rows <- nrow(thinned)
	Thinning = Fit$Thinning
	


    Dev <- matrix(Fit$Deviance, ncol=1)
    Mon <- Fit$Monitor
    Num.Mon <- ncol(Mon)
    
	thinned2 = Fit$Posterior2
	
	if ( !is.null(nrow(thinned2)) && (nrow(thinned2) >= 1)) 
		{
		HD <- BMK.Diagnostic(thinned2, batches = 10)
		Ind <- 1 * (HD > 0.5)
		BurnIn <- thinned.rows
		batch.list <- seq(from = 1, to = nrow(thinned2), by = floor(nrow(thinned2)/10))
		for (i in 1:9) 
			{
			if (sum(Ind[, i:9]) == 0)
				{
				BurnIn <- batch.list[i] - 1
				break
				}
			}
	
		Stat.at <- BurnIn + 1
		#rm(batch.list, HD, Ind, thinned2)
    	}
    
    cat("Assessing Thinning and ESS\n")
    acf.temp <- matrix(1, trunc(10 * log10(thinned.rows)), LIV)
    ESS1 <- Rec.Thin <- rep(1, LIV)
    for (j in 1:LIV) {
        temp0 <- acf(thinned[, j], lag.max = nrow(acf.temp), 
            plot = FALSE)
        acf.temp[, j] <- abs(temp0$acf[2:{
            nrow(acf.temp) + 1
        }, , 1])
        ESS1[j] <- ESS(thinned[, j])
        Rec.Thin[j] <- which(acf.temp[, j] <= 0.1)[1] * Thinning
    }
    Rec.Thin <- ifelse(is.na(Rec.Thin), nrow(acf.temp), Rec.Thin)
    ESS2 <- ESS(Dev)
    ESS3 <- ESS(Mon)
 
	if ( !is.null(nrow(thinned2)) && (nrow(thinned2) >= 1)) 
		{
		if (Stat.at < thinned.rows)
			{
			ESS4 <- ESS(thinned[Stat.at:thinned.rows, ])
			ESS5 <- ESS(Dev[Stat.at:thinned.rows, ])
			ESS6 <- ESS(Mon[Stat.at:thinned.rows, ])
		    }
		}
    cat("Creating Summaries\n")
    Num.Mon <- ncol(Mon)
    Summ1 <- matrix(NA, LIV, 7, dimnames = list(MyData$parm.names, 
        c("Mean", "SD", "MCSE", "ESS", "LB", "Median", "UB")))
    Summ1[, 1] <- colMeans(thinned)
    Summ1[, 2] <- apply(thinned, 2, sd)
    Summ1[, 3] <- 0
    Summ1[, 4] <- ESS1
    Summ1[, 5] <- apply(thinned, 2, quantile, c(0.025), na.rm = TRUE)
    Summ1[, 6] <- apply(thinned, 2, quantile, c(0.5), na.rm = TRUE)
    Summ1[, 7] <- apply(thinned, 2, quantile, c(0.975), na.rm = TRUE)
    for (i in 1:ncol(thinned)) {
        temp <- try(MCSE(thinned[, i]), silent = TRUE)
        if (!inherits(temp, "try-error")) 
            Summ1[i, 3] <- temp
        else Summ1[i, 3] <- MCSE(thinned[, i], method = "sample.variance")
    }
    Deviance <- rep(NA, 7)
    Deviance[1] <- mean(Dev)
    Deviance[2] <- sd(as.vector(Dev))
    temp <- try(MCSE(as.vector(Dev)), silent = TRUE)
    if (inherits(temp, "try-error")) 
        temp <- MCSE(as.vector(Dev), method = "sample.variance")
    Deviance[3] <- temp
    Deviance[4] <- ESS2
    Deviance[5] <- as.numeric(quantile(Dev, probs = 0.025, na.rm = TRUE))
    Deviance[6] <- as.numeric(quantile(Dev, probs = 0.5, na.rm = TRUE))
    Deviance[7] <- as.numeric(quantile(Dev, probs = 0.975, na.rm = TRUE))
    Summ1 <- rbind(Summ1, Deviance)
    for (j in 1:Num.Mon) {
        Monitor <- rep(NA, 7)
        Monitor[1] <- mean(Mon[, j])
        Monitor[2] <- sd(as.vector(Mon[, j]))
        temp <- try(MCSE(as.vector(Mon[, j])), silent = TRUE)
        if (inherits(temp, "try-error")) 
            temp <- MCSE(Mon[, j], method = "sample.variance")
        Monitor[3] <- temp
        Monitor[4] <- ESS3[j]
        Monitor[5] <- as.numeric(quantile(Mon[, j], probs = 0.025, 
            na.rm = TRUE))
        Monitor[6] <- as.numeric(quantile(Mon[, j], probs = 0.5, 
            na.rm = TRUE))
        Monitor[7] <- as.numeric(quantile(Mon[, j], probs = 0.975, 
            na.rm = TRUE))
        Summ1 <- rbind(Summ1, Monitor)
        rownames(Summ1)[nrow(Summ1)] <- MyData$mon.names[j]
    }
    
	if ( !is.null(nrow(thinned2)) && (nrow(thinned2) >= 1)) 
		{

		Summ2 <- matrix(NA, LIV, 7, dimnames = list(MyData$parm.names, 
			c("Mean", "SD", "MCSE", "ESS", "LB", "Median", "UB")))
		if (Stat.at < thinned.rows)
			{
			thinned2 <- matrix(thinned[Stat.at:thinned.rows, ], thinned.rows - 
				Stat.at + 1, ncol(thinned))
			Dev2 <- matrix(Dev[Stat.at:thinned.rows, ], thinned.rows - 
				Stat.at + 1, ncol(Dev))
			Mon2 <- matrix(Mon[Stat.at:thinned.rows, ], thinned.rows - 
				Stat.at + 1, ncol(Mon))
			Summ2[, 1] <- colMeans(thinned2)
			Summ2[, 2] <- apply(thinned2, 2, sd)
			Summ2[, 3] <- 0
			Summ2[, 4] <- ESS4
			Summ2[, 5] <- apply(thinned2, 2, quantile, c(0.025), 
				na.rm = TRUE)
			Summ2[, 6] <- apply(thinned2, 2, quantile, c(0.5), na.rm = TRUE)
			Summ2[, 7] <- apply(thinned2, 2, quantile, c(0.975), 
				na.rm = TRUE)
			for (i in 1:ncol(thinned2))
				{
				temp <- try(MCSE(thinned2[, i]), silent = TRUE)
				if (!inherits(temp, "try-error")) 
					Summ2[i, 3] <- temp
				else Summ2[i, 3] <- MCSE(thinned2[, i], method = "sample.variance")
				}
			Deviance <- rep(NA, 7)
			Deviance[1] <- mean(Dev2)
			Deviance[2] <- sd(as.vector(Dev2))
			temp <- try(MCSE(as.vector(Dev2)), silent = TRUE)
			if (inherits(temp, "try-error")) 
				temp <- MCSE(as.vector(Dev2), method = "sample.variance")
			Deviance[3] <- temp
			Deviance[4] <- ESS5
			Deviance[5] <- as.numeric(quantile(Dev2, probs = 0.025, 
				na.rm = TRUE))
			Deviance[6] <- as.numeric(quantile(Dev2, probs = 0.5, 
				na.rm = TRUE))
			Deviance[7] <- as.numeric(quantile(Dev2, probs = 0.975, 
				na.rm = TRUE))
			Summ2 <- rbind(Summ2, Deviance)
			for (j in 1:Num.Mon)
				{
				Monitor <- rep(NA, 7)
				Monitor[1] <- mean(Mon2[, j])
				Monitor[2] <- sd(as.vector(Mon2[, j]))
				temp <- try(MCSE(as.vector(Mon[, j])), silent = TRUE)
				if (inherits(temp, "try-error")) 
					temp <- MCSE(as.vector(Mon[, j]), method = "sample.variance")
				Monitor[3] <- temp
				Monitor[4] <- ESS6[j]
				Monitor[5] <- as.numeric(quantile(Mon2[, j], probs = 0.025, 
					na.rm = TRUE))
				Monitor[6] <- as.numeric(quantile(Mon2[, j], probs = 0.5, 
					na.rm = TRUE))
				Monitor[7] <- as.numeric(quantile(Mon2[, j], probs = 0.975, 
					na.rm = TRUE))
				Summ2 <- rbind(Summ2, Monitor)
				rownames(Summ2)[nrow(Summ2)] <- MyData$mon.names[j]
				}
			}
		Fit$Summary2 = Summ2
    	}
    
    Fit$Summary1 = Summ1
    
    # Update Deviance, DIC1, DIC2 (actually shouldn't change)
    Fit$Deviance = as.vector(Dev)
    
    Fit$DIC1 = c(mean(as.vector(Dev)), var(as.vector(Dev))/2, 
            mean(as.vector(Dev)) + var(as.vector(Dev))/2)
    
	if ( !is.null(nrow(thinned2)) && (nrow(thinned2) >= 1)) 
		{
		if (Stat.at < thinned.rows)
			{
			Fit$DIC2 = c(mean(as.vector(Dev2)), var(as.vector(Dev2))/2, mean(as.vector(Dev2)) + var(as.vector(Dev2))/2)
			} else {
			Fit$DIC2 = rep(NA, 3)
			}
		}
    
    #######################################################
	# Find out how much of each parameter equals 0.0
	#######################################################
	zeros_TF = Fit$Posterior1 == 0.0
	Posterior1_iszero = colSums(zeros_TF) / nrow(Fit$Posterior1)
	Fit$Posterior1_iszero = Posterior1_iszero
        
	if ( !is.null(nrow(Fit$Posterior2)) && (nrow(Fit$Posterior2) >= 1)) 
		{
		zeros_TF = Fit$Posterior2 == 0.0
		Posterior2_iszero = colSums(zeros_TF) / nrow(Fit$Posterior2)
	    Fit$Posterior2_iszero = Posterior2_iszero
		}
    
    
    return(Fit)
	}








#######################################################
# reparam_LapDem_output
#######################################################
#' Change the burnin fraction when re-parameterizing LaplacesDemon MCMC output
#'
#' \code{\link[LaplacesDemon]{LaplacesDemon}} likes to run its MCMC sampling on a simple number line.  Thus, the likelihood
#' function etc. should transform the numbers into the range desired, e.g. 0-1.
#'
#' This function changes the burnin fraction
#'
#' @param Fit The \code{\link[LaplacesDemon]{LaplacesDemon}} output object.
#' @param MyData The \code{\link[LaplacesDemon]{LaplacesDemon}} input data.
#' @param burnfract The new burnin fraction
#' @param newsamps If you want to subset the samples.
#' @return \code{Fit} The modified MCMC output.
#' @export
#' @seealso \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite LaplacesDemon_Tutorial
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
change_burnin <- function(Fit, MyData, burnfract=0.1, newsamps=NULL)
	{
	require(LaplacesDemon)
	
	# Safety
	origFit = Fit
	
	# Get the data
	BioGeoBEARS_model_object = MyData$BioGeoBEARS_model_object
	params_table = BioGeoBEARS_model_object@params_table
	freeTF = params_table$type == "free"

	
	# Cut the number of rows by burnfract
	rowstart = round(burnfract * nrow(Fit$Posterior1)) + 1
	rowend = nrow(Fit$Posterior1)

	# The NULL is NOT transferred to seq successfully
	if (!is.null(newsamps))
		{
		rownums = seq(from=rowstart, to=rowend, length.out=newsamps)
		} else {
		rownums = seq(from=rowstart, to=rowend)
		}
	
	# Cut the tables down
	Fit$Posterior1 = Fit$Posterior1[rownums, ]
	Fit$Monitor = Fit$Monitor[rownums, ]
	Fit$Deviance = Fit$Deviance[rownums]
	
	# Revise summary also
	Fit$Thinned.Samples = nrow(Fit$Posterior1)
	Fit$Iterations = nrow(Fit$Posterior1)
	
	
	# Maybe also do Posterior2 (should be unnecessary)
	dothis = FALSE
	if (dothis && !is.null(nrow(Fit$Posterior2)) && (nrow(Fit$Posterior2) >= 1))
		{
		# Cut the number of rows by burnfract
		rowstart = round(burnfract * nrow(Fit$Posterior2)) + 1
		rowend = nrow(Fit$Posterior2)
		# The NULL is transferred to seq successfully
		rownums = seq(start=rowstart, end=rowend, length.out=newsamps)	

		Fit$Posterior2 = Fit$Posterior2[rownums, ]
		}


	#######################################################
	# Re-do the summary, from the parameter-space posterior sample
	#######################################################
	
	
	#######################################################
	# re-do the summary, copying from LaplacesDemon()
	#######################################################
	
	# Posterior1 is the thinned samples, including everything
	LIV = length(Fit$Initial.Values)
	thinned = Fit$Posterior1
	thinned.rows <- nrow(thinned)
	Thinning = Fit$Thinning
	


    Dev <- matrix(Fit$Deviance, ncol=1)
    Mon <- Fit$Monitor
    Num.Mon <- ncol(Mon)
    
	thinned2 = Fit$Posterior2
	
	if ( !is.null(nrow(thinned2)) && (nrow(thinned2) >= 1)) 
		{
		HD <- BMK.Diagnostic(thinned2, batches = 10)
		Ind <- 1 * (HD > 0.5)
		BurnIn <- thinned.rows
		batch.list <- seq(from = 1, to = nrow(thinned2), by = floor(nrow(thinned2)/10))
		for (i in 1:9) 
			{
			if (sum(Ind[, i:9]) == 0)
				{
				BurnIn <- batch.list[i] - 1
				break
				}
			}
	
		Stat.at <- BurnIn + 1
		#rm(batch.list, HD, Ind, thinned2)
    	}
    
    cat("Assessing Thinning and ESS\n")
    acf.temp <- matrix(1, trunc(10 * log10(thinned.rows)), LIV)
    ESS1 <- Rec.Thin <- rep(1, LIV)
    for (j in 1:LIV) {
        temp0 <- acf(thinned[, j], lag.max = nrow(acf.temp), 
            plot = FALSE)
        acf.temp[, j] <- abs(temp0$acf[2:{
            nrow(acf.temp) + 1
        }, , 1])
        ESS1[j] <- ESS(thinned[, j])
        Rec.Thin[j] <- which(acf.temp[, j] <= 0.1)[1] * Thinning
    }
    Rec.Thin <- ifelse(is.na(Rec.Thin), nrow(acf.temp), Rec.Thin)
    ESS2 <- ESS(Dev)
    ESS3 <- ESS(Mon)
 
	if ( !is.null(nrow(thinned2)) && (nrow(thinned2) >= 1)) 
		{
		if (Stat.at < thinned.rows)
			{
			ESS4 <- ESS(thinned[Stat.at:thinned.rows, ])
			ESS5 <- ESS(Dev[Stat.at:thinned.rows, ])
			ESS6 <- ESS(Mon[Stat.at:thinned.rows, ])
		    }
		}
    cat("Creating Summaries\n")
    Num.Mon <- ncol(Mon)
    Summ1 <- matrix(NA, LIV, 7, dimnames = list(MyData$parm.names, 
        c("Mean", "SD", "MCSE", "ESS", "LB", "Median", "UB")))
    Summ1[, 1] <- colMeans(thinned)
    Summ1[, 2] <- apply(thinned, 2, sd)
    Summ1[, 3] <- 0
    Summ1[, 4] <- ESS1
    Summ1[, 5] <- apply(thinned, 2, quantile, c(0.025), na.rm = TRUE)
    Summ1[, 6] <- apply(thinned, 2, quantile, c(0.5), na.rm = TRUE)
    Summ1[, 7] <- apply(thinned, 2, quantile, c(0.975), na.rm = TRUE)
    for (i in 1:ncol(thinned)) {
        temp <- try(MCSE(thinned[, i]), silent = TRUE)
        if (!inherits(temp, "try-error")) 
            Summ1[i, 3] <- temp
        else Summ1[i, 3] <- MCSE(thinned[, i], method = "sample.variance")
    }
    Deviance <- rep(NA, 7)
    Deviance[1] <- mean(Dev)
    Deviance[2] <- sd(as.vector(Dev))
    temp <- try(MCSE(as.vector(Dev)), silent = TRUE)
    if (inherits(temp, "try-error")) 
        temp <- MCSE(as.vector(Dev), method = "sample.variance")
    Deviance[3] <- temp
    Deviance[4] <- ESS2
    Deviance[5] <- as.numeric(quantile(Dev, probs = 0.025, na.rm = TRUE))
    Deviance[6] <- as.numeric(quantile(Dev, probs = 0.5, na.rm = TRUE))
    Deviance[7] <- as.numeric(quantile(Dev, probs = 0.975, na.rm = TRUE))
    Summ1 <- rbind(Summ1, Deviance)
    for (j in 1:Num.Mon) {
        Monitor <- rep(NA, 7)
        Monitor[1] <- mean(Mon[, j])
        Monitor[2] <- sd(as.vector(Mon[, j]))
        temp <- try(MCSE(as.vector(Mon[, j])), silent = TRUE)
        if (inherits(temp, "try-error")) 
            temp <- MCSE(Mon[, j], method = "sample.variance")
        Monitor[3] <- temp
        Monitor[4] <- ESS3[j]
        Monitor[5] <- as.numeric(quantile(Mon[, j], probs = 0.025, 
            na.rm = TRUE))
        Monitor[6] <- as.numeric(quantile(Mon[, j], probs = 0.5, 
            na.rm = TRUE))
        Monitor[7] <- as.numeric(quantile(Mon[, j], probs = 0.975, 
            na.rm = TRUE))
        Summ1 <- rbind(Summ1, Monitor)
        rownames(Summ1)[nrow(Summ1)] <- MyData$mon.names[j]
    }
    
	if ( !is.null(nrow(thinned2)) && (nrow(thinned2) >= 1)) 
		{

		Summ2 <- matrix(NA, LIV, 7, dimnames = list(MyData$parm.names, 
			c("Mean", "SD", "MCSE", "ESS", "LB", "Median", "UB")))
		if (Stat.at < thinned.rows)
			{
			thinned2 <- matrix(thinned[Stat.at:thinned.rows, ], thinned.rows - 
				Stat.at + 1, ncol(thinned))
			Dev2 <- matrix(Dev[Stat.at:thinned.rows, ], thinned.rows - 
				Stat.at + 1, ncol(Dev))
			Mon2 <- matrix(Mon[Stat.at:thinned.rows, ], thinned.rows - 
				Stat.at + 1, ncol(Mon))
			Summ2[, 1] <- colMeans(thinned2)
			Summ2[, 2] <- apply(thinned2, 2, sd)
			Summ2[, 3] <- 0
			Summ2[, 4] <- ESS4
			Summ2[, 5] <- apply(thinned2, 2, quantile, c(0.025), 
				na.rm = TRUE)
			Summ2[, 6] <- apply(thinned2, 2, quantile, c(0.5), na.rm = TRUE)
			Summ2[, 7] <- apply(thinned2, 2, quantile, c(0.975), 
				na.rm = TRUE)
			for (i in 1:ncol(thinned2))
				{
				temp <- try(MCSE(thinned2[, i]), silent = TRUE)
				if (!inherits(temp, "try-error")) 
					Summ2[i, 3] <- temp
				else Summ2[i, 3] <- MCSE(thinned2[, i], method = "sample.variance")
				}
			Deviance <- rep(NA, 7)
			Deviance[1] <- mean(Dev2)
			Deviance[2] <- sd(as.vector(Dev2))
			temp <- try(MCSE(as.vector(Dev2)), silent = TRUE)
			if (inherits(temp, "try-error")) 
				temp <- MCSE(as.vector(Dev2), method = "sample.variance")
			Deviance[3] <- temp
			Deviance[4] <- ESS5
			Deviance[5] <- as.numeric(quantile(Dev2, probs = 0.025, 
				na.rm = TRUE))
			Deviance[6] <- as.numeric(quantile(Dev2, probs = 0.5, 
				na.rm = TRUE))
			Deviance[7] <- as.numeric(quantile(Dev2, probs = 0.975, 
				na.rm = TRUE))
			Summ2 <- rbind(Summ2, Deviance)
			for (j in 1:Num.Mon)
				{
				Monitor <- rep(NA, 7)
				Monitor[1] <- mean(Mon2[, j])
				Monitor[2] <- sd(as.vector(Mon2[, j]))
				temp <- try(MCSE(as.vector(Mon[, j])), silent = TRUE)
				if (inherits(temp, "try-error")) 
					temp <- MCSE(as.vector(Mon[, j]), method = "sample.variance")
				Monitor[3] <- temp
				Monitor[4] <- ESS6[j]
				Monitor[5] <- as.numeric(quantile(Mon2[, j], probs = 0.025, 
					na.rm = TRUE))
				Monitor[6] <- as.numeric(quantile(Mon2[, j], probs = 0.5, 
					na.rm = TRUE))
				Monitor[7] <- as.numeric(quantile(Mon2[, j], probs = 0.975, 
					na.rm = TRUE))
				Summ2 <- rbind(Summ2, Monitor)
				rownames(Summ2)[nrow(Summ2)] <- MyData$mon.names[j]
				}
			}
		Fit$Summary2 = Summ2
    	}
    
    Fit$Summary1 = Summ1
    
    # Update Deviance, DIC1, DIC2 (actually shouldn't change)
    Fit$Deviance = as.vector(Dev)
    
    Fit$DIC1 = c(mean(as.vector(Dev)), var(as.vector(Dev))/2, 
            mean(as.vector(Dev)) + var(as.vector(Dev))/2)
    
	if ( !is.null(nrow(thinned2)) && (nrow(thinned2) >= 1)) 
		{
		if (Stat.at < thinned.rows)
			{
			Fit$DIC2 = c(mean(as.vector(Dev2)), var(as.vector(Dev2))/2, mean(as.vector(Dev2)) + var(as.vector(Dev2))/2)
			} else {
			Fit$DIC2 = rep(NA, 3)
			}
		}
    
    #######################################################
	# Find out how much of each parameter equals 0.0
	#######################################################
	zeros_TF = Fit$Posterior1 == 0.0
	Posterior1_iszero = colSums(zeros_TF) / nrow(Fit$Posterior1)
	Fit$Posterior1_iszero = Posterior1_iszero
        
	if ( !is.null(nrow(Fit$Posterior2)) && (nrow(Fit$Posterior2) >= 1)) 
		{
		zeros_TF = Fit$Posterior2 == 0.0
		Posterior2_iszero = colSums(zeros_TF) / nrow(Fit$Posterior2)
	    Fit$Posterior2_iszero = Posterior2_iszero
		}
    
    
    return(Fit)
	}









#######################################################
# Bayesian models, for LaplacesDemon
#######################################################
# params = c(0.1, 0.1, 0.3)
# Calculate the log posterior probability


#######################################################
# LapDem_Model_calc_logPP
#######################################################
#' The function calculating posterior probability of parameter values
#'
#' What it says
#'
#' @param params the parameter values
#' @param MyData The \code{\link[LaplacesDemon]{LaplacesDemon}} input data.
#' @return \code{Modelout} The MCMC output.
#' @export
#' @seealso \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite LaplacesDemon_Tutorial
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
LapDem_Model_calc_logPP <- function(params, MyData)
	{
	require(LaplacesDemon) 	# for logit, invlogit
	
	orig_params = params
	
	# Get the data
	tr = MyData$tr
	tipranges = MyData$tipranges
	BioGeoBEARS_model_object = MyData$BioGeoBEARS_model_object
	
	params_table = BioGeoBEARS_model_object@params_table
	freeTF = params_table$type == "free"
	
	# NO: Convert the params (range -Inf -> +Inf) to the min/max values
	# convert anything outside the range, to min/max, using the interval function
	
	params
	
	for (p in 1:length(params))
		{
		# Check for parameters being set to 0 by RJMCMC
		if (params[p] == 0.0)
			{
			params[p] = 0.0
			} else {
		
			# Inverse logit of param
			# params[p] = invlogit(params[p])
			
			# Reset things outside the boundary to (min-max)
			minval = params_table$min[freeTF][p]
			maxval = params_table$max[freeTF][p]
			
			params[p] = transform_with_logistic(paramval=params[p], minval, maxval)
				
			}
		}

	params

	# NO, INTERVAL() CRASHES with mapply or anything else
	# 		d = 0.1
	# 		e = 0.1
	# 		j = 0.3
	# 		params = c(d,e,j)
	# 		#params = c(d,e)
	# 		
	# 		params_table = BioGeoBEARS_model_object@params_table
	# 		freeTF = params_table$type == "free"
	# 	
	# 		
	# 		# Re-parameterize/constrain with LaplacesDemon:::interval()	
	# 		params_not_zero_TF = params != 0.0
	# 		minvals = params_table$min[freeTF][params_not_zero_TF]
	# 		maxvals = params_table$min[freeTF][params_not_zero_TF]
	# 		mapply(
	
	# Input the current parameters into the BioGeoBEARS model object
	BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object, params)
	
	# Update the dependent parameters accordingly
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	# Input the basics for calc_loglike_for_optim()
	tip_condlikes_of_data_on_each_state = MyData$tip_condlikes_of_data_on_each_state
	print_optim = MyData$print_optim
	areas_list = MyData$areas_list
	states_list = MyData$states_list
	force_sparse = MyData$force_sparse
	cluster_already_open = MyData$cluster_already_open

	
	
	### Parameters
	params
	
	### Log of Prior Densities (these are all uniform)
	#beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
	#sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
	params_table = BioGeoBEARS_model_object@params_table
	freeTF = params_table$type == "free"
	rownums = (1:nrow(params_table))[freeTF]
	LnPrior = 0
	for (rn in rownums)
		{
		# Assume a Unif(min, max) prior for all variables
		LnPrior = LnPrior + dunif(x=params_table$est[rn], min=params_table$min[rn], max=params_table$max[rn], log=TRUE)
		
		# Assume a very flat, but normal, prior, for all variables
		# (This will avoid -- NAH, doesn't work cause variables are in regular flat space, not logit space
		# NO DOESNT WORK
		#LnPrior = LnPrior + rnorm(x=params_table$est[rn], min=params_table$min[rn], max=params_table$max[rn], log=TRUE)
		
		
		# Try a beta prior with shape1 & shape2 close to 1.00001
		plotthis = '
		x = seq(0,1,0.01)
		y=dbeta(x, shape1=1.01, shape2=1.01)
		plot(x,y)
		dbeta(x, shape1=1.01, shape2=1.01, log=TRUE)
		'
		
# 		if ( (rownames(params_table)[rn] == "d") || (rownames(params_table)[rn] == "e"))
# 			{
# 			# beta prior for d and e
# 			LnPrior = LnPrior + dbeta(x=params_table$est[rn], shape1=1.000001, shape2=1.000001, log=TRUE)
# 			} else {
# 			# unif prior for j etc.
# 			LnPrior = LnPrior + dunif(x=params_table$est[rn], min=params_table$min[rn], max=params_table$max[rn], log=TRUE)
# 			}
# 		
		
		}
	
	
	### Log-Likelihood
	# mus are the coefficients of the betas
	# LL is the sum of the logs of the deviance/densities under the model
	#mu <- tcrossprod(Data$X, t(beta))
	#LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
	LnL = calc_loglike_for_optim(params, BioGeoBEARS_model_object, phy=tr, tip_condlikes_of_data_on_each_state, print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
	
	
	### Log-Posterior
	LnP <- LnL + LnPrior
	
	### yhat = data simulated under these model parameters, for posterior predictive sampling/testing
	### (I think) 
	# yhat = rnorm(length(mu), mu, sigma)
	yhat = NA
	


	# Convert the params back to (range -Inf -> +Inf) to the min/max values
# 	for (p in 1:length(params))
# 		{
# 		# Scale to logit
# 		minval = params_table$min[freeTF][p]
# 		maxval = params_table$max[freeTF][p]
# 		
# 		
# 		if (params[p] == 0.0)
# 			{
# 			# Convert to -70 in logit space, to keep parameter "off the bottom" of precision
# 			params[p] = -70.0
# 			} else {
# 			# Don't subtract minval, as that results in 0, which results in -Inf under logit
# 			params[p] = (params[p] - minval) / (maxval-minval)
# 			
# 			# logit of param
# 			params[p] = logit(params[p])
# 			}
# 		}

		
	### Model output to return
	Modelout <- list(LP=LnP, Dev=-2*LnL, Monitor=c(LnP,LnL,LnPrior), yhat=yhat, parm=orig_params)
	return(Modelout)
	}




# Calculate the log posterior probability
# This version attempts rescaling rather than re-writing parameters...

#######################################################
# LapDem_Model_calc_logPP_orig
#######################################################
#' The function calculating posterior probability of parameter values
#'
#' What it says
#'
#' @param params the parameter values
#' @param MyData The \code{\link[LaplacesDemon]{LaplacesDemon}} input data.
#' @return \code{Modelout} The MCMC output.
#' @export
#' @seealso \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite LaplacesDemon_Tutorial
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
LapDem_Model_calc_logPP_orig <- function(params, MyData)
	{
	require(LaplacesDemon) 	# for logit, invlogit
	
	# Get the data
	tr = MyData$tr
	tipranges = MyData$tipranges
	BioGeoBEARS_model_object = MyData$BioGeoBEARS_model_object
	
	params_table = BioGeoBEARS_model_object@params_table
	freeTF = params_table$type == "free"
	
	# Convert the params (range -Inf -> +Inf) to the min/max values
	for (p in 1:length(params))
		{
		# Check for parameters being set to 0 by RJMCMC
		if (params[p] == 0.0)
			{
			params[p] = 0.0
			} else {
		
			# Inverse logit of param
			params[p] = invlogit(params[p])
			
			# Scale to (min-max)
			minval = params_table$min[freeTF][p]
			maxval = params_table$max[freeTF][p]
			
			params[p] = minval + (maxval-minval) * params[p]
			}
		}


	
	# Input the current parameters into the BioGeoBEARS model object
	BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object, params)
	
	# Update the dependent parameters accordingly
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	# Input the basics for calc_loglike_for_optim()
	tip_condlikes_of_data_on_each_state = MyData$tip_condlikes_of_data_on_each_state
	print_optim = MyData$print_optim
	areas_list = MyData$areas_list
	states_list = MyData$states_list
	force_sparse = MyData$force_sparse
	cluster_already_open = MyData$cluster_already_open

	
	
	### Parameters
	params
	
	### Log of Prior Densities (these are all uniform)
	#beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
	#sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
	params_table = BioGeoBEARS_model_object@params_table
	freeTF = params_table$type == "free"
	rownums = (1:nrow(params_table))[freeTF]
	LnPrior = 0
	for (rn in rownums)
		{
		# Assume a Unif(min, max) prior for all variables
		#LnPrior = LnPrior + dunif(x=params_table$est[rn], min=params_table$min[rn], max=params_table$max[rn], log=TRUE)
		
		# Assume a very flat, but normal, prior, for all variables
		# (This will avoid -- NAH, doesn't work cause variables are in regular flat space, not logit space
		# NO DOESNT WORK
		#LnPrior = LnPrior + rnorm(x=params_table$est[rn], min=params_table$min[rn], max=params_table$max[rn], log=TRUE)
		
		
		# Try a beta prior with shape1 & shape2 close to 1.00001
		plotthis = '
		x = seq(0,1,0.01)
		y=dbeta(x, shape1=1.01, shape2=1.01)
		plot(x,y)
		dbeta(x, shape1=1.01, shape2=1.01, log=TRUE)
		'
		
		if ( (rownames(params_table)[rn] == "d") || (rownames(params_table)[rn] == "e"))
			{
			# beta prior for d and e
			LnPrior = LnPrior + dbeta(x=params_table$est[rn], shape1=1.000001, shape2=1.000001, log=TRUE)
			} else {
			# unif prior for j etc.
			LnPrior = LnPrior + dunif(x=params_table$est[rn], min=params_table$min[rn], max=params_table$max[rn], log=TRUE)
			}
		
		
		}
	
	
	### Log-Likelihood
	# mus are the coefficients of the betas
	# LL is the sum of the logs of the deviance/densities under the model
	#mu <- tcrossprod(Data$X, t(beta))
	#LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
	LnL = calc_loglike_for_optim(params, BioGeoBEARS_model_object, phy=tr, tip_condlikes_of_data_on_each_state, print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
	
	
	### Log-Posterior
	LnP <- LnL + LnPrior
	
	### yhat = data simulated under these model parameters, for posterior predictive sampling/testing
	### (I think) 
	# yhat = rnorm(length(mu), mu, sigma)
	yhat = NA
	


	# Convert the params back to (range -Inf -> +Inf) to the min/max values
	for (p in 1:length(params))
		{
		# Scale to logit
		minval = params_table$min[freeTF][p]
		maxval = params_table$max[freeTF][p]
		
		
		if (params[p] == 0.0)
			{
			# Convert to -70 in logit space, to keep parameter "off the bottom" of precision
			params[p] = -70.0
			} else {
			# Don't subtract minval, as that results in 0, which results in -Inf under logit
			params[p] = (params[p] - minval) / (maxval-minval)
			
			# logit of param
			params[p] = logit(params[p])
			}
		}

		
	### Model output to return
	Modelout <- list(LP=LnP, Dev=-2*LnL, Monitor=c(LnP,LnL,LnPrior), yhat=yhat, parm=params)
	return(Modelout)
	}












#######################################################
# PLotting functions
#######################################################
#LaplacesDemon:::plot.demonoid <- function (x, BurnIn = 0, Data = NULL, PDF = FALSE, Parms = NULL, 
#   ...) 


#######################################################
# plot_demonoid2
#######################################################
#' Plot LaplacesDemon MCMC output
#'
#' Modified version of \code{\link[LaplacesDemon]{plot.demonoid}} to plot to PDF.
#'
#' @param x The MCMC output
#' @param BurnIn The burnin amount
#' @param Data The \code{\link[LaplacesDemon]{LaplacesDemon}} input data.
#' @param PDF If TRUE, plot to PDF.
#' @param Parms The parameters.
#' @param ... Additional arguments to standard functions
#' @return Nothing.
#' @export
#' @seealso \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite LaplacesDemon_Tutorial
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
plot_demonoid2 <- function (x, BurnIn = 0, Data = NULL, PDF = FALSE, Parms = NULL, 
    ...) 
{
	greyval1 = 0.75
	greyval1 = rgb(greyval1, greyval1, greyval1)
	greyval2 = 0.35
	greyval2 = rgb(greyval2, greyval2, greyval2)
	
    if (missing(x)) 
        stop("The x argument is required.")
    if (class(x) != "demonoid") 
        stop("x must be of class demonoid.")
    if (is.null(Data)) 
        stop("The Data argument is NULL.")
    if (BurnIn >= nrow(x$Posterior1)) 
        BurnIn <- 0
    Stat.at <- BurnIn + 1
    if (is.null(Parms)) {
        Posterior <- x$Posterior1
    }
    else {
        Parms <- sub("\\[", "\\\\[", Parms)
        Parms <- sub("\\]", "\\\\]", Parms)
        Parms <- sub("\\.", "\\\\.", Parms)
        if (length(grep(Parms[1], colnames(x$Posterior1))) == 
            0) 
            stop("Parameter in Parms does not exist.")
        keepcols <- grep(Parms[1], colnames(x$Posterior1))
        if (length(Parms) > 1) {
            for (i in 2:length(Parms)) {
                if (length(grep(Parms[i], colnames(x$Posterior1))) == 
                  0) 
                  stop("Parameter in Parms does not exist.")
                keepcols <- c(keepcols, grep(Parms[i], colnames(x$Posterior1)))
            }
        }
        Posterior <- as.matrix(x$Posterior1[, keepcols])
        colnames(Posterior) <- colnames(x$Posterior1)[keepcols]
    }
    if (PDF == TRUE) {
        pdf("LaplacesDemon.Plots.pdf")
        par(mfrow = c(3, 3))
    }
    else {
        par(mfrow = c(3, 3), ask = TRUE)
    }
    for (j in 1:ncol(Posterior)) {
        plot(Stat.at:x$Thinned.Samples, Posterior[Stat.at:x$Thinned.Samples, 
            j], type = "l", xlab = "Generations", ylab = "Value", 
            main = colnames(Posterior)[j])
        panel.smooth(Stat.at:x$Thinned.Samples, Posterior[Stat.at:x$Thinned.Samples, 
            j], pch = "")

		###########################################
		# Plot the parameter values
		###########################################
		numbreaks = 30
        vals = pretty(Posterior[Stat.at:x$Thinned.Samples, j])
        breakpoints = seq(from=min(vals), to=max(vals), length.out=numbreaks)
        
        # Get values for density plot
        histres2 = hist(x=Posterior[Stat.at:x$Thinned.Samples, j], breaks=100, plot=FALSE)

        histres = hist(x=Posterior[Stat.at:x$Thinned.Samples, j], breaks=breakpoints, freq=FALSE, add=FALSE, xlab="Value", ylab="Density", 
        	col=greyval1, border=greyval2, main=colnames(Posterior)[j])
        lines(density(unlist(mapply(rep,histres2$mids,histres2$counts)), adjust=1), col="black", lwd=1.5)
        #polygon(density(Posterior[Stat.at:x$Thinned.Samples, j]), col = "white", border = "black")

        abline(v = 0, col = "red", lty = 2)
        if (!is.constant(Posterior[Stat.at:x$Thinned.Samples, 
            j])) {
            z <- acf(Posterior[Stat.at:x$Thinned.Samples, j], 
                plot = FALSE)
            se <- 1/sqrt(length(Posterior[Stat.at:x$Thinned.Samples, 
                j]))
            plot(z$lag, z$acf, ylim = c(min(z$acf, -2 * se), 
                1), type = "h", main = colnames(Posterior)[j], 
                xlab = "Lag", ylab = "Correlation")
            abline(h = (2 * se), col = "red", lty = 2)
            abline(h = (-2 * se), col = "red", lty = 2)
        }
        else {
            plot(0, 0, main = paste(colnames(Posterior)[j], "is a constant."))
        }
    }
    plot(Stat.at:length(x$Deviance), x$Deviance[Stat.at:length(x$Deviance)], 
        type = "l", xlab = "Iterations", ylab = "Value", main = "Deviance")
    panel.smooth(Stat.at:length(x$Deviance), x$Deviance[Stat.at:length(x$Deviance)], 
        pch = "")
    
    
    
    #plot(density(x$Deviance[Stat.at:length(x$Deviance)]), xlab = "Value", 
    #    main = "Deviance")
    #polygon(density(x$Deviance[Stat.at:length(x$Deviance)]), 
    #    col = "black", border = "black")
	# Get values for density plot

	###########################################
	# Plot the Deviance values
	###########################################
	numbreaks = 30
	vals = pretty(x$Deviance[Stat.at:length(x$Deviance)])
	breakpoints = seq(from=min(vals), to=max(vals), length.out=numbreaks)
	histres2 = hist(x=x$Deviance[Stat.at:length(x$Deviance)], breaks=100, plot=FALSE)
	histres = hist(x=x$Deviance[Stat.at:length(x$Deviance)], breaks=breakpoints, freq=FALSE, add=FALSE, xlab="Value", ylab="Density", 
		col=greyval1, border=greyval2, main="Deviance")
	lines(density(unlist(mapply(rep,histres2$mids,histres2$counts)), adjust=1), col="black", lwd=1.5)

    
    
    
    abline(v = 0, col = "red", lty = 2)
    if (!is.constant(x$Deviance[Stat.at:length(x$Deviance)])) {
        z <- acf(x$Deviance[Stat.at:length(x$Deviance)], plot = FALSE)
        se <- 1/sqrt(length(x$Deviance[Stat.at:length(x$Deviance)]))
        plot(z$lag, z$acf, ylim = c(min(z$acf, -2 * se), 1), 
            type = "h", main = "Deviance", xlab = "Lag", ylab = "Correlation")
        abline(h = (2 * se), col = "red", lty = 2)
        abline(h = (-2 * se), col = "red", lty = 2)
    }
    else {
        plot(0, 0, main = "Deviance is a constant.")
    }
    if (is.vector(x$Monitor)) {
        J <- 1
        nn <- length(x$Monitor)
    }
    else if (is.matrix(x$Monitor)) {
        J <- ncol(x$Monitor)
        nn <- nrow(x$Monitor)
    }
    for (j in 1:J) {
        plot(Stat.at:nn, x$Monitor[Stat.at:nn, j], type = "l", 
            xlab = "Iterations", ylab = "Value", main = Data$mon.names[j])
        panel.smooth(Stat.at:nn, x$Monitor[Stat.at:nn, j], pch = "")
        #plot(density(x$Monitor[Stat.at:nn, j]), xlab = "Value", 
        #    main = Data$mon.names[j])
        #polygon(density(x$Monitor[Stat.at:nn, j]), col = "black", 
        #    border = "black")

		###########################################
		# Plot the LP, LnL, LnPrior values
		###########################################
		numbreaks = 30
		vals = pretty(x$Monitor[Stat.at:nn, j])
		breakpoints = seq(from=min(vals), to=max(vals), length.out=numbreaks)
		histres2 = hist(x=x$Monitor[Stat.at:nn, j], breaks=100, plot=FALSE)
		histres = hist(x=x$Monitor[Stat.at:nn, j], breaks=breakpoints, freq=FALSE, add=FALSE, xlab="Value", ylab="Density", 
			col=greyval1, border=greyval2, main=Data$mon.names[j])
		lines(density(unlist(mapply(rep,histres2$mids,histres2$counts)), adjust=1), col="black", lwd=1.5)
		
		
        
        abline(v = 0, col = "red", lty = 2)
        if (!is.constant(x$Monitor[Stat.at:nn, j])) {
            z <- acf(x$Monitor[Stat.at:nn, j], plot = FALSE)
            se <- 1/sqrt(length(x$Monitor[Stat.at:nn, j]))
            plot(z$lag, z$acf, ylim = c(min(z$acf, -2 * se), 
                1), type = "h", main = Data$mon.names[j], xlab = "Lag", 
                ylab = "Correlation")
            abline(h = (2 * se), col = "red", lty = 2)
            abline(h = (-2 * se), col = "red", lty = 2)
        }
        else {
            plot(0, 0, main = paste(Data$mon.names[j], "is a constant."))
        }
    }
    if (nrow(x$CovarDHis) > 1) {
        if (x$Algorithm %in% c("Adaptive Hamiltonian Monte Carlo", 
            "Hamiltonian Monte Carlo with Dual-Averaging", "No-U-Turn Sampler")) {
            plot(x$CovarDHis[, 1], type = "l", xlab = "Adaptations", 
                ylab = expression(epsilon))
        }
        else {
            Diff <- abs(diff(x$CovarDHis))
            adaptchange <- matrix(NA, nrow(Diff), 3)
            for (i in 1:nrow(Diff)) {
                adaptchange[i, 1:3] <- as.vector(quantile(Diff[i, 
                  ], probs = c(0.025, 0.5, 0.975)))
            }
            plot(adaptchange[, 2], ylim = c(min(adaptchange), 
                max(adaptchange)), type = "l", col = "red", xlab = "Adaptations", 
                ylab = "Absolute Difference", main = "Proposal Variance", 
                sub = "Median=Red, 95% Bounds=Gray")
            lines(adaptchange[, 1], col = "gray")
            lines(adaptchange[, 3], col = "gray")
        }
    }
    if (PDF == TRUE) 
        dev.off()
}





