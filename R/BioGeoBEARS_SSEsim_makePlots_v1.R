

plot_simDEC_DECJ_inferences <- function(simdir=NULL, sim_params_Rdata_table="/simdata/BGB/dej_params.Rdata", modelnull="DEC", modelalt="DEC+J", makePDFs=TRUE, openPDFs=TRUE)
	{
	defaults='
	simdir = "/simdata/BGB/ps0011_sim003/"
	sim_params_Rdata_table="/simdata/BGB/dej_params.Rdata"
	modelnull="DEC"
	modelalt="DEC+J"
	makePDFs=TRUE
	openPDFs=TRUE
	'
	
	# If you just want the stats, not the plots
	if (makePDFs == FALSE)
		{
		juststats = TRUE
		openPDFs = FALSE
		} else {
		juststats = FALSE
		}
	
	
	#######################################################
	# Get whatever directory your simulations are in
	#######################################################
	if (is.null(simdir))
		{
		cat("\n\nWARNING: no simulation directory specified, using default on Nick Matzke's computer.\n\n")
		simdir = "/simdata/BGB/ps0011_sim003/"
		}
	

	# Set the working directory to the simulation directory
	setwd(simdir)
	
	# Process the directory to get the row of the simulation parameters
	words = strsplit(simdir, "/")[[1]]
	dir_txt = words[length(words)]

	# Get the row of parameter values
	words = strsplit(dir_txt, "_")[[1]]
	word = words[1]
	param_iter = as.numeric(gsub(pattern="ps", replacement="", x=word))
	param_iter


	# Get the parameters for this run
	load(file="/simdata/BGB/dej_params.Rdata")
	nums = 1:nrow(dej_params)
	dej_params = cbind(nums, dej_params)
	dim(dej_params)
	param_vals = dej_params[param_iter,]
	param_vals




	#######################################################
	# Read saved results
	#######################################################
	# Read the results
	resDEC_fn = "DEC_inf.Rdata"
	load(resDEC_fn)
	resDEC = res

	resDECJ_fn = "DECJ_inf.Rdata"
	load(resDECJ_fn)
	resDECj = res


	SSEsim_results_fn = "SSEsim_results_processed.Rdata"
	load(SSEsim_results_fn)
	SSEsim_results_processed




	# Inputs for graphics and stats
	tr = read.tree(resDEC$inputs$trfn)
	tipranges = getranges_from_LagrangePHYLIP(resDEC$inputs$geogfn)
	BioGeoBEARS_run_object = resDEC$inputs



	# Results table
	restable = NULL
	teststable = NULL

	
	if (makePDFs==TRUE)
		{
		pdffn = paste(dir_txt, "_DEC_vs_DECj_SSE.pdf", sep="")
		pdf(pdffn, width=8.5, height=11)
		} # end makePDFs
	
	#######################################################
	# Plot ancestral states - DEC
	#######################################################
	analysis_titletxt = paste("BioGeoBEARS ", modelnull, " on ", dir_txt, " d=", param_vals$d, " e=", param_vals$e, " j=", param_vals$j, " brate=", param_vals$brate, " drate=", param_vals$drate, " bexp=", param_vals$b_exp, " dexp=", param_vals$d_exp, sep="")

	# Setup
	results_object = resDEC
	scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

	# States
	res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, juststats=juststats)
	row.names(res1) = modelnull
	
	# Pie chart
	plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, juststats=juststats)

	#######################################################
	# Plot ancestral states - DECJ
	#######################################################
	analysis_titletxt = paste("BioGeoBEARS ", modelalt, " on ", dir_txt, " d=", param_vals$d, " e=", param_vals$e, " j=", param_vals$j, " brate=", param_vals$brate, " drate=", param_vals$drate, " bexp=", param_vals$b_exp, " dexp=", param_vals$d_exp, sep="")

	# Setup
	results_object = resDECj
	scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

	# States
	res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, juststats=juststats)
	row.names(res2) = modelalt
	
	# Pie chart
	plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, juststats=juststats)

	if (makePDFs==TRUE)
		{
		dev.off()
		if (openPDFs==TRUE)
			{
			cmdstr = paste("open ", pdffn, sep="")
			system(cmdstr)
			}
		}


	#######################################################
	# Stats
	#######################################################
	# We have to extract the log-likelihood differently, depending on the 
	# version of optim/optimx
	if (BioGeoBEARS_run_object$use_optimx == TRUE)
		{
		# Using optimx() results
		if (packageVersion("optimx") < 2013)
			{
			# optimx 2012
			LnL_2 = as.numeric(resDEC$optim_result$fvalues)
			LnL_1 = as.numeric(resDECj$optim_result$fvalues)
			} else {
			# optimx 2013
			LnL_2 = as.numeric(resDEC$optim_result$value)
			LnL_1 = as.numeric(resDECj$optim_result$value)
			} # end optimx 2012 vs. 2013
		} else {
		# Using optim() results
		LnL_2 = as.numeric(resDEC$optim_result$value)
		LnL_1 = as.numeric(resDECj$optim_result$value)
		} # end optim vs. optimx

	numparams1 = 3
	numparams2 = 2
	stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
	row.names(stats) = paste(modelnull, "_v_", modelalt, sep="")
	stats

	res1
	res2

	rbind(res1, res2)
	tmp_tests = conditional_format_table(stats)

	restable = rbind(restable, res1, res2)
	teststable = rbind(teststable, tmp_tests)


	restable
	teststable

	cat("\nParameter values for this simulation:\n")
	print(param_vals)
	
	
	# Assemble output list
	simDEC_DECJ_inf_results = NULL	
	simDEC_DECJ_inf_results$SSEsim_results_processed = SSEsim_results_processed
	simDEC_DECJ_inf_results$res1 = res1
	simDEC_DECJ_inf_results$res2 = res2
	simDEC_DECJ_inf_results$tests = tmp_tests
	
	return(simDEC_DECJ_inf_results)
	}




hist_sim_v_inf_vals <- function(inferred_vals, trueval)
	{
	h1 = hist(inferred_vals, plot=FALSE)
	xticks = c(0, pretty(c(0, h1$mids), n=1))
	yticks = c(0, pretty(c(0,h1$counts), n=1))
	xlims = c(0, max(xticks))
	ylims = c(0, max(yticks))
	plot(h1$mids, h1$counts, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n", xlim=xlims, ylim=ylims)
	plot(h1, add=TRUE)
	
	abline(v=trueval, lwd=2, col="blue", lty="dashed")
	box()
	# x-axis tick marks
	axis(side=1, at=xticks, labels=TRUE, tick=TRUE)
	# y-axis tick marks (not really needed for histograms)
	#axis(side=1, at=xticks, labels=TRUE, tick=TRUE)
	
	return(h1)
	}
















#######################################################
# Make plots of DEC vs. DEC+J parameter inference and model choice
#######################################################



plot_inf_accuracy_vs_SSEsims_v2 <- function(plot_inputs, xwidth_manual=NULL)
	{
	

defaults='
titletxt1 = bquote(paste("100 SSE simulations (", lambda, "=0.3, ", mu, "=0.3, ", alpha, "=1, ", omega, "=-1)", sep=""))
titletxt2 = paste("DEC (white) vs. DEC+J (grey) inference", sep="")
pdffn = "DEC_DECJ_inf_boxplots_SSE_b03d03bx1dx-1.pdf"

simtype = "SSE"
segwidth_truth = 0.05
plot_inputs$segwidth_truth = segwidth_truth

simtype = "SSE"
brate_plot = 0.3
drate_plot = 0.3
b_exp_plot = 1
d_exp_plot = -1

plot_inputs = NULL
plot_inputs$titletxt1 = NULL
plot_inputs$pdffn = pdffn
plot_inputs$doPDF = TRUE
plot_inputs$siminf_stats_good = siminf_stats_good
plot_inputs$simtype = simtype
plot_inputs$brate_plot = brate_plot
plot_inputs$drate_plot = drate_plot
plot_inputs$b_exp_plot = b_exp_plot
plot_inputs$d_exp_plot = d_exp_plot
plot_inputs$uniq_params_txt = uniq_params_txt
plot_inputs$segwidth_truth = segwidth_truth

'# end defaults

titletxt1 = plot_inputs$titletxt1
pdffn = plot_inputs$pdffn
doPDF = plot_inputs$doPDF
siminf_stats_good = plot_inputs$siminf_stats_good
brate_plot = plot_inputs$brate_plot
drate_plot = plot_inputs$drate_plot
b_exp_plot = plot_inputs$b_exp_plot
d_exp_plot = plot_inputs$d_exp_plot
uniq_params_txt = plot_inputs$uniq_params_txt
segwidth_truth = plot_inputs$segwidth_truth
simtype = plot_inputs$simtype

if (is.null(titletxt1))
	{
	titletxt1 = bquote(paste(.(simtype), " simulations (", lambda, "=", .(brate_plot), ", ",  mu, "=", .(drate_plot), ", ", alpha, "=", .(b_exp_plot), ", ", omega, "=", .(d_exp_plot), ")", sep=""))
	}



#######################################################
# Spacer between models
#######################################################
xspace_value = 0.5

if (is.null(xwidth_manual))
	{
	xwidth = uniq_params_num + (uniq_params_num-1) * xspace_value
	} else {
	xwidth = xwidth_manual
	}

#######################################################
# Make plots of DEC vs. DEC+J parameter inference and model choice
#######################################################
if(plot_inputs$doPDF == TRUE)
	{
	pdf(pdffn, height=11, width=8.5)
	}

par(oma = c(5,4,4,1))
par(mar=c(0,2,1,1))
nums = 1:5
mat = matrix(data=nums, nrow=5, ncol=5, byrow=FALSE)
mat

layout(mat)

i=1
ymax = 0.25
yaxis_ticks = c(0, 0.05, 0.1, 0.15, 0.2, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("d")), side=2, las=1, line=4)

# Plot the page title
mtext(text=titletxt1, side=3, line=0, cex=1.2, outer=TRUE)
#mtext(text=titletxt2, side=3, line=-0.5, cex=0.8, outer=TRUE)

# Plot a legend
xpos = 1
points(x=xpos, y=0.925*ymax, pch=18, cex=2, col="black")
segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=0.925*ymax, y1=0.925*ymax, lwd=2, col="black")
text(x=xpos+1.5*segwidth_truth, y=0.925*ymax, labels="true parameter value", pos=4, offset=1)
points(x=xpos, y=0.825*ymax, pch=22, cex=3, col="black", bg="white")
text(x=xpos+1.5*segwidth_truth, y=0.825*ymax, labels="DEC inference", pos=4, offset=1)
points(x=xpos, y=0.725*ymax, pch=22, cex=3, col="black", bg="lightgray")
text(x=xpos+1.5*segwidth_truth, y=0.725*ymax, labels="DEC+J inference", pos=4, offset=1)



# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	dtrue = unique(siminf_stats_sim$d)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$d_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$d_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=dtrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=dtrue, y1=dtrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through d truth/inference plots



i=1
ymax = 0.25
yaxis_ticks = c(0, 0.05, 0.1, 0.15, 0.2, ymax)

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("e")), side=2, las=1, line=4)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	etrue = unique(siminf_stats_sim$e)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$e_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$e_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=etrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=etrue, y1=etrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through e truth/inference plots



i=1
ymax = 0.5
yaxis_ticks = pretty(c(0, ymax))
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("j")), side=2, las=1, line=4)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$j_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$j_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=jtrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=jtrue, y1=jtrue, lwd=2, col="black")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through j truth/inference plots





#######################################################
# Accuracy
#######################################################
i=1
ymax = 1
yaxis_ticks =c(0, 1/15, 0.25, 0.5, 0.75, ymax)
yaxis_ticks_txt =c("0", "1/15", "0.25", "0.5", "0.75", "1")

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks_txt, las=1)
mtext(text="Proportion correct", side=2, las=3, line=3.5, cex=1)

# Line representing random guess among states
abline(h=1/15, lty="dotted", col="black", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	
	# Calculate accuracy
	MLstate_TF_ancstate_accuracies_DEC_subset = ancstate_accuracies_DEC_subset[TF,99:147]
	MLstate_accuracy_DEC = apply(X=MLstate_TF_ancstate_accuracies_DEC_subset, MARGIN=1, FUN=mean, na.rm=TRUE)

	MLstate_TF_ancstate_accuracies_DECJ_subset = ancstate_accuracies_DECJ_subset[TF,99:147]
	MLstate_accuracy_DECJ = apply(X=MLstate_TF_ancstate_accuracies_DECJ_subset, MARGIN=1, FUN=mean, na.rm=TRUE)
	
	
	jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=MLstate_accuracy_DEC, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=MLstate_accuracy_DECJ, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	#points(x=i, y=jtrue, pch="*", cex=4, col="darkblue")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through accuracy boxplots








#######################################################
# Model choice -- with AICc
#######################################################
i=1
ymax = 130
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 50, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote(paste(Delta, "AICc", sep=""))
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3.5, cex=1)

# Dashed red line at the significance level
abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	LnL_DEC = siminf_stats_sim$LnL_DEC
	LnL_DECJ = siminf_stats_sim$LnL_DECJ
	
	# Boxplot of the log-likelihood advantage of DEC+J
	LnL_advantage = LnL_DECJ-LnL_DEC
	# Min on this plot is 0.1
	# Also a few optimization issues at high d, e values
	# (could be datasets with no signal left; functions will misfire if LnL DEC>LnL DEC+J)
	TF = LnL_advantage < 0
	LnL_DECJ[TF] = LnL_DEC[TF]
	
	if (sum(TF) > 0)
		{
		cat("Note: In ", sum(TF), " cases, LnL DEC > LnL DEC+J.\n", sep="")
		# This may be very minor, or may indicate a failure of optimx ML routine
		# to find the true optimum. For now, in these cases we are setting LnL 
		# DEC+J to equal LnL DEC.
		}
	
	# Calculate DeltaAICc 
	AICc_DEC = getAICc(LnL=LnL_DEC, numparams=2, samplesize=100)
	AICc_DECJ = getAICc(LnL=LnL_DECJ, numparams=3, samplesize=100)
	delta_AICc_vals = AICc_DEC - AICc_DECJ

	# Plot AICc weights
	# DEC or DEC+J
	xpos = extra_xspace + i + 0.0
	
	# Plot one box (grey or white)
	# AICc weights look stupid
	jtrue = unique(siminf_stats_sim$j)
	if (jtrue > 0)
		{
		# width makes fatter boxes
		b1 = boxplot(x=delta_AICc_vals, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n", col="lightgray", border="black", width=3)
		} else {
		# width makes fatter boxes
		b1 = boxplot(x=delta_AICc_vals, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n", col="lightgray", border="black", width=3)
		}
	
	# Add "true model:" text
	if (i==1)
		{
		txtplot = "True model:"
		mtext(side=1, at=-1, text=txtplot, cex=0.55, las=1, line=0)
		}
	
	if (jtrue > 0)
		{
		txtplot = "+J"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		} else {
		txtplot = "DEC"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		}
	
		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot


if(plot_inputs$doPDF == TRUE)
	{
	dev.off()
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)
	}



	
	} # end plot_inf_accuracy_vs_SSEsims_v2()
























#######################################################
# Do a subset plot for publication
#######################################################

plot_inf_accuracy_vs_SSEsims_v3<- function(plot_inputs, xwidth_manual=NULL)
	{
	

defaults='
titletxt1 = bquote(paste("100 SSE simulations (", lambda, "=0.3, ", mu, "=0.3, ", alpha, "=1, ", omega, "=-1)", sep=""))
titletxt2 = paste("DEC (white) vs. DEC+J (grey) inference", sep="")
pdffn = "DEC_DECJ_inf_boxplots_SSE_b03d03bx1dx-1.pdf"

simtype = "SSE"
segwidth_truth = 0.05
plot_inputs$segwidth_truth = segwidth_truth

simtype = "SSE"
brate_plot = 0.3
drate_plot = 0.3
b_exp_plot = 1
d_exp_plot = -1

plot_inputs = NULL
plot_inputs$titletxt1 = NULL
plot_inputs$pdffn = pdffn
plot_inputs$doPDF = TRUE
plot_inputs$siminf_stats_good = siminf_stats_good
plot_inputs$simtype = simtype
plot_inputs$brate_plot = brate_plot
plot_inputs$drate_plot = drate_plot
plot_inputs$b_exp_plot = b_exp_plot
plot_inputs$d_exp_plot = d_exp_plot
plot_inputs$uniq_params_txt = uniq_params_txt
plot_inputs$segwidth_truth = segwidth_truth

'# end defaults

titletxt1 = plot_inputs$titletxt1
pdffn = plot_inputs$pdffn
doPDF = plot_inputs$doPDF
siminf_stats_good = plot_inputs$siminf_stats_good
brate_plot = plot_inputs$brate_plot
drate_plot = plot_inputs$drate_plot
b_exp_plot = plot_inputs$b_exp_plot
d_exp_plot = plot_inputs$d_exp_plot
uniq_params_txt = plot_inputs$uniq_params_txt
segwidth_truth = plot_inputs$segwidth_truth
simtype = plot_inputs$simtype

if (is.null(titletxt1))
	{
	titletxt1 = bquote(paste(.(simtype), " simulations (", lambda, "=", .(brate_plot), ", ",  mu, "=", .(drate_plot), ", ", alpha, "=", .(b_exp_plot), ", ", omega, "=", .(d_exp_plot), ")", sep=""))
	}



#######################################################
# Spacer between models
#######################################################
xspace_value = 0.5

if (is.null(xwidth_manual))
	{
	xwidth = uniq_params_num + (uniq_params_num-1) * xspace_value
	} else {
	xwidth = xwidth_manual
	}

#######################################################
# Make plots of DEC vs. DEC+J parameter inference and model choice
#######################################################
if(plot_inputs$doPDF == TRUE)
	{
	pdf(pdffn, height=6, width=6)
	}

par(oma = c(5,4,4,1))
par(mar=c(0,2,1,1))
nums = 1:4
mat = matrix(data=nums, nrow=4, ncol=4, byrow=FALSE)
mat

layout(mat)

i=1
ymax = 0.20
yaxis_ticks = c(0, 0.05, 0.1, 0.15, 0.2, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("d")), side=2, las=1, line=4)

# Plot the page title
mtext(text=titletxt1, side=3, line=0, cex=1.2, outer=TRUE)
#mtext(text=titletxt2, side=3, line=-0.5, cex=0.8, outer=TRUE)

# Plot a legend
xpos = 1
points(x=xpos, y=0.925*ymax, pch=18, cex=2, col="black")
segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=0.925*ymax, y1=0.925*ymax, lwd=2, col="black")
text(x=xpos+1.5*segwidth_truth, y=0.925*ymax, labels="true parameter value", pos=4, offset=1)
points(x=xpos, y=0.725*ymax, pch=22, cex=3, col="black", bg="white")
text(x=xpos+1.5*segwidth_truth, y=0.725*ymax, labels="DEC inference", pos=4, offset=1)
points(x=xpos, y=0.525*ymax, pch=22, cex=3, col="black", bg="lightgray")
text(x=xpos+1.5*segwidth_truth, y=0.525*ymax, labels="DEC+J inference", pos=4, offset=1)



# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	dtrue = unique(siminf_stats_sim$d)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$d_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$d_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=dtrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=dtrue, y1=dtrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through d truth/inference plots




# Cut "e"




i=1
ymax = 0.5
yaxis_ticks = pretty(c(0, ymax))
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("j")), side=2, las=1, line=4)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$j_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$j_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=jtrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=jtrue, y1=jtrue, lwd=2, col="black")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through j truth/inference plots





#######################################################
# Accuracy
#######################################################
i=1
ymax = 1
yaxis_ticks =c(1/15, 0.25, 0.5, 0.75, ymax)
yaxis_ticks_txt =c("1/15", "0.25", "0.5", "0.75", "1")

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks_txt, las=1)
mtext(text="Proportion correct", side=2, las=3, line=3.5, cex=1)

# Line representing random guess among states
abline(h=1/15, lty="dotted", col="black", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	
	# Calculate accuracy
	MLstate_TF_ancstate_accuracies_DEC_subset = ancstate_accuracies_DEC_subset[TF,99:147]
	MLstate_accuracy_DEC = apply(X=MLstate_TF_ancstate_accuracies_DEC_subset, MARGIN=1, FUN=mean, na.rm=TRUE)

	MLstate_TF_ancstate_accuracies_DECJ_subset = ancstate_accuracies_DECJ_subset[TF,99:147]
	MLstate_accuracy_DECJ = apply(X=MLstate_TF_ancstate_accuracies_DECJ_subset, MARGIN=1, FUN=mean, na.rm=TRUE)
	
	
	jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=MLstate_accuracy_DEC, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=MLstate_accuracy_DECJ, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	#points(x=i, y=jtrue, pch="*", cex=4, col="darkblue")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through accuracy boxplots








#######################################################
# Model choice -- with AICc
#######################################################
i=1
ymax = 130
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 50, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote(paste(Delta, "AICc", sep=""))
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3.5, cex=1)

# Dashed red line at the significance level
abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	LnL_DEC = siminf_stats_sim$LnL_DEC
	LnL_DECJ = siminf_stats_sim$LnL_DECJ
	
	# Boxplot of the log-likelihood advantage of DEC+J
	LnL_advantage = LnL_DECJ-LnL_DEC
	# Min on this plot is 0.1
	# Also a few optimization issues at high d, e values
	# (could be datasets with no signal left; functions will misfire if LnL DEC>LnL DEC+J)
	TF = LnL_advantage < 0
	LnL_DECJ[TF] = LnL_DEC[TF]
	
	if (sum(TF) > 0)
		{
		cat("Note: In ", sum(TF), " cases, LnL DEC > LnL DEC+J.\n", sep="")
		# This may be very minor, or may indicate a failure of optimx ML routine
		# to find the true optimum. For now, in these cases we are setting LnL 
		# DEC+J to equal LnL DEC.
		}
	
	# Calculate DeltaAICc 
	AICc_DEC = getAICc(LnL=LnL_DEC, numparams=2, samplesize=100)
	AICc_DECJ = getAICc(LnL=LnL_DECJ, numparams=3, samplesize=100)
	delta_AICc_vals = AICc_DEC - AICc_DECJ

	# Plot AICc weights
	# DEC or DEC+J
	xpos = extra_xspace + i + 0.0
	
	# Plot one box (grey or white)
	# AICc weights look stupid
	jtrue = unique(siminf_stats_sim$j)
	if (jtrue > 0)
		{
		# width makes fatter boxes
		b1 = boxplot(x=delta_AICc_vals, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n", col="lightgray", border="black", width=3)
		} else {
		# width makes fatter boxes
		b1 = boxplot(x=delta_AICc_vals, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n", col="lightgray", border="black", width=3)
		}
	
	# Add "true model:" text
	if (i==1)
		{
		txtplot = "True model:"
		mtext(side=1, at=-0.25, text=txtplot, cex=0.55, las=1, line=0)
		}
	
	if (jtrue > 0)
		{
		txtplot = "+J"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		} else {
		txtplot = "DEC"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		}
	
		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot


if(plot_inputs$doPDF == TRUE)
	{
	dev.off()
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)
	}



	
	} # end plot_inf_accuracy_vs_SSEsims_v2()






















#######################################################
# Make plots of Tree Statistics for these simulations,
# to show variability in the simulated trees
#######################################################



#######################################################
# Make plots of keys statistics for simulated tree/event histories
#######################################################

plot_SSEsims_treestats_pt1 <- function(plot_inputs, xwidth_manual=NULL, sizes=c(0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4))
	{
	

defaults='
titletxt1 = bquote(paste("100 SSE simulations (", lambda, "=0.3, ", mu, "=0.3, ", alpha, "=1, ", omega, "=-1)", sep=""))
titletxt2 = paste("fossil (white) vs. observed (grey) tree stats", sep="")
pdffn = "simtree_stats_SSE_b03d03bx1dx-1.pdf"

simtype = "SSE"
segwidth_truth = 0.05
plot_inputs$segwidth_truth = segwidth_truth

simtype = "SSE"
brate_plot = 0.3
drate_plot = 0.3
b_exp_plot = 1
d_exp_plot = -1

plot_inputs = NULL
plot_inputs$titletxt1 = NULL
plot_inputs$pdffn = pdffn
plot_inputs$doPDF = TRUE
plot_inputs$siminf_stats_good = siminf_stats_good
plot_inputs$sim_tipstates_good = sim_tipstates_good
plot_inputs$ancstate_accuracies_DEC_good = ancstate_accuracies_DEC_good
plot_inputs$ancstate_accuracies_DECJ_good = ancstate_accuracies_DECJ_good
plot_inputs$simevent_counts_good = simevent_counts_good
plot_inputs$treeheights_incl_fossil_good = treeheights_incl_fossil_good

plot_inputs$simtype = simtype
plot_inputs$brate_plot = brate_plot
plot_inputs$drate_plot = drate_plot
plot_inputs$b_exp_plot = b_exp_plot
plot_inputs$d_exp_plot = d_exp_plot
plot_inputs$uniq_params_txt = uniq_params_txt
plot_inputs$segwidth_truth = segwidth_truth

sizes = c(0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4)
'# end defaults


# Load the SSE parameters, subset the data

titletxt1 = plot_inputs$titletxt1
pdffn = plot_inputs$pdffn
doPDF = plot_inputs$doPDF
siminf_stats_good = plot_inputs$siminf_stats_good
sim_tipstates_good = plot_inputs$sim_tipstates_good
ancstate_accuracies_DEC_good = plot_inputs$ancstate_accuracies_DEC_good
ancstate_accuracies_DECJ_good = plot_inputs$ancstate_accuracies_DECJ_good
simevent_counts_good = plot_inputs$simevent_counts_good
treeheights_incl_fossil_good = plot_inputs$treeheights_incl_fossil_good

brate_plot = plot_inputs$brate_plot
drate_plot = plot_inputs$drate_plot
b_exp_plot = plot_inputs$b_exp_plot
d_exp_plot = plot_inputs$d_exp_plot
uniq_params_txt = plot_inputs$uniq_params_txt
segwidth_truth = plot_inputs$segwidth_truth
simtype = plot_inputs$simtype

TF1 = siminf_stats_good$brate == brate_plot
TF2 = siminf_stats_good$drate == drate_plot
TF3 = siminf_stats_good$b_exp == b_exp_plot
TF4 = siminf_stats_good$d_exp == d_exp_plot
TF = (TF1 + TF2 + TF3 + TF4) == 4


# Do the subsetting
siminf_stats_subset = dfnums_to_numeric(siminf_stats_good[TF, ])
sim_tipstates_subset = sim_tipstates_good[TF, ]
ancstate_accuracies_DEC_subset = ancstate_accuracies_DEC_good[TF,]
ancstate_accuracies_DECJ_subset = ancstate_accuracies_DECJ_good[TF,]
simevent_counts_subset = simevent_counts_good[TF,]
treeheights_incl_fossil_subset = treeheights_incl_fossil_good[TF]


if (is.null(titletxt1))
	{
	titletxt1 = bquote(paste("Tree stats for ", .(simtype), " simulations (", lambda, "=", .(brate_plot), ", ",  mu, "=", .(drate_plot), ", ", alpha, "=", .(b_exp_plot), ", ", omega, "=", .(d_exp_plot), ")", sep=""))
	}
titletxt1


#######################################################
# Spacer between models
#######################################################
xspace_value = 0.5

if (is.null(xwidth_manual))
	{
	xwidth = uniq_params_num + (uniq_params_num-1) * xspace_value
	} else {
	xwidth = xwidth_manual
	}

#######################################################
# Make plots of Tree Statistics
#######################################################
if(plot_inputs$doPDF == TRUE)
	{
	pdf(pdffn, height=11, width=8.5)
	}

par(oma = c(5,4,4,1))
par(mar=c(0,2,1,1))
nums = 1:5
mat = matrix(data=nums, nrow=5, ncol=5, byrow=FALSE)
mat

layout(mat)


#######################################################
# Tree heights, all and observed
#######################################################

i=1
ymax = 150
yaxis_ticks = c(0, 50, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote("Tree height"), side=2, las=3, line=3)

# Plot the page title
mtext(text=titletxt1, side=3, line=0, cex=1.2, outer=TRUE)
#mtext(text=titletxt2, side=3, line=-0.5, cex=0.8, outer=TRUE)

# Plot a legend
xpos = 1
#points(x=xpos, y=0.925*ymax, pch=18, cex=2, col="black")
#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=0.925*ymax, y1=0.925*ymax, lwd=2, col="black")
#text(x=xpos+1.5*segwidth_truth, y=0.925*ymax, labels="true parameter value", pos=4, offset=1)
points(x=xpos, y=0.925*ymax, pch=22, cex=3, col="black", bg="white")
text(x=xpos+1.5*segwidth_truth, y=0.925*ymax, labels="with fossils", pos=4, offset=1)
points(x=xpos, y=0.825*ymax, pch=22, cex=3, col="black", bg="gray")
text(x=xpos+1.5*segwidth_truth, y=0.825*ymax, labels="observed tree", pos=4, offset=1)



# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	observed_tree_ages = ancstate_accuracies_DEC_subset[TF,1]
	observed_tree_ages[c(1:5, 96:100)]
	actual_tree_ages = treeheights_incl_fossil_subset[TF]
	actual_tree_ages[c(1:5, 96:100)]
	
	
	#dtrue = unique(siminf_stats_sim$d)
	
	# Actual tree age (with fossils)
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=actual_tree_ages, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Observed tree ages (excluding fossils, extinct nodes, etc.)
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=observed_tree_ages, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	#points(x=xpos, y=dtrue, pch=18, cex=2, col="black")
	#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=dtrue, y1=dtrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through d truth/inference plots



#######################################################
# Number of tips, fossil and observed
#######################################################

i=1
ymax = 400
yaxis_ticks = c(0, 200, ymax)

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text="# species", side=2, las=3, line=3)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	simevent_counts_sim
	

	# Line representing number of species observed
	#abline(h=50, lty="dotted", col="black", lwd=0.5)
	axis(side=2, at=50, labels="50 species\nobserved", las=1, cex=0.5)
	#mtext(side=2, at=50, text="observed: 50", las=1, line=1, cex=0.55, lty=1)


	#etrue = unique(siminf_stats_sim$e)
	
	# Number of tips in true tree, including fossils
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=50+simevent_counts_sim$num_fossil_tips, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Number of tips in observed tree (always 50)
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=50+simevent_counts_sim$num_fossil_tips-simevent_counts_sim$num_fossil_tips, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	#points(x=xpos, y=etrue, pch=18, cex=2, col="black")
	#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=etrue, y1=etrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through e truth/inference plots



#######################################################
# Tip range sizes
#######################################################
# Range sizes
#sizes = c(0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4)
#sizes


foo <- function(val, sizes)
	{
	sizes[val]
	}



i=1
ymax = 4
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 1, 2, 3, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote("range size"), side=2, las=3, line=3)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	sim_tipstates_sim = sim_tipstates_subset[TF,]
	sim_tipstates_sim2 = sim_tipstates_sim - 1
	
	sim_tipsizes = sapply(X=sim_tipstates_sim2, FUN=foo, sizes=sizes)
	sim_tipsizes2 = matrix(data=sim_tipsizes, nrow=nrow(sim_tipstates_sim2), ncol=ncol(sim_tipstates_sim2))

	meanSizes = rowMeans(sim_tipsizes2, na.rm=FALSE)


	
	#jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.0
	b1 = boxplot(x=meanSizes, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	#points(x=xpos, y=jtrue, pch=18, cex=2, col="black")
	#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=jtrue, y1=jtrue, lwd=2, col="black")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through j truth/inference plots





#######################################################
# Number of range-change events
#######################################################

i=1
ymax = 400
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 200, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote("# events (total)"), side=2, las=3, line=3)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$y_actual + simevent_counts_sim$s_actual + simevent_counts_sim$v_actual + simevent_counts_sim$j_actual + simevent_counts_sim$d_actual + simevent_counts_sim$e_actual + simevent_counts_sim$a_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$y_observed + simevent_counts_sim$s_observed + simevent_counts_sim$v_observed + simevent_counts_sim$j_observed + simevent_counts_sim$d_observed + simevent_counts_sim$e_observed + simevent_counts_sim$a_observed
	
	
	#jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	#points(x=xpos, y=jtrue, pch=18, cex=2, col="black")
	#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=jtrue, y1=jtrue, lwd=2, col="black")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through j truth/inference plots





#######################################################
# Num events - cladogenetic, nonsympatric
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text="# events (clad.)", side=2, las=3, line=3, cex=1)

# Line representing random guess among states
#abline(h=1/15, lty="dotted", col="black", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$s_actual + simevent_counts_sim$v_actual + simevent_counts_sim$j_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$s_observed + simevent_counts_sim$v_observed + simevent_counts_sim$j_observed
		
	#jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	#points(x=i, y=jtrue, pch="*", cex=4, col="darkblue")




	# Plot one box (grey or white)
	# AICc weights look stupid
	siminf_stats_sim = siminf_stats_subset[TF,]
	jtrue = unique(siminf_stats_sim$j)
	xpos = extra_xspace + i+0.0
	
	# Add "true model:" text
	if (i==1)
		{
		txtplot = "True model:"
		mtext(side=1, at=-1, text=txtplot, cex=0.55, las=1, line=0)
		}
	
	if (jtrue > 0)
		{
		txtplot = "+J"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		} else {
		txtplot = "DEC"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		}
	



	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through accuracy boxplots







par(oma = c(5,4,4,1))
par(mar=c(0,2,1,1))
nums = 1:6
mat = matrix(data=nums, nrow=6, ncol=6, byrow=FALSE)
mat

layout(mat)



#######################################################
# Number of "d" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (d)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]


	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$d_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$d_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")





		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot





#######################################################
# 2nd page (counts of e, y, s, v, j)
#######################################################


titletxt1_cont = bquote(paste("Tree stats for ", .(simtype), " simulations, cont. (", lambda, "=", .(brate_plot), ", ",  mu, "=", .(drate_plot), ", ", alpha, "=", .(b_exp_plot), ", ", omega, "=", .(d_exp_plot), ")", sep=""))

titletxt1_cont



#######################################################
# Number of "e" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (e)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)


# Plot the page title
mtext(text=titletxt1_cont, side=3, line=0, cex=1.2, outer=TRUE)


# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$e_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$e_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot






#######################################################
# Number of "y" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (y)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$y_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$y_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot









#######################################################
# Number of "s" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (s)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$s_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$s_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot










#######################################################
# Number of "v" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (v)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$v_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$v_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot






#######################################################
# Number of "j" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (j)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$j_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$j_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")



	# Plot one box (grey or white)
	# AICc weights look stupid
	siminf_stats_sim = siminf_stats_subset[TF,]
	jtrue = unique(siminf_stats_sim$j)
	xpos = extra_xspace + i+0.0
	
	# Add "true model:" text
	if (i==1)
		{
		txtplot = "True model:"
		mtext(side=1, at=-1, text=txtplot, cex=0.55, las=1, line=0)
		}
	
	if (jtrue > 0)
		{
		txtplot = "+J"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		} else {
		txtplot = "DEC"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		}
	

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot





if(plot_inputs$doPDF == TRUE)
	{
	dev.off()
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)
	}



	
	} # end plot_SSEsims_treestats_pt1()













