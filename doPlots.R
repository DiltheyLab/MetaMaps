library("RColorBrewer")

args = commandArgs(trailingOnly=TRUE)

capitalize <- function(x)
{
	if(nchar(x) > 0)
	{
		x2 <- x
		substr(x2, 1, 1) <- toupper(substr(x2, 1, 1))
		x2
	}
	else
	{
		x
	}
}

simulationsDirectory <- "databases/miniSeq_100/simulations_logNormal"

if(length(args) > 0)
{
	simulationsDirectory <- args[[1]]
}

cat(paste("doPlots.R\n\tDirectory: ", simulationsDirectory, "\n\n", paste = ""))

prefix <- paste(simulationsDirectory, "/_forPlot_", sep = "")

barplotPal <- brewer.pal(n = 5, name = "PuBu")

pdf(paste(simulationsDirectory, "/plots_frequencies.pdf", sep = ""), height = 10, width = 10)
library(ggplot2)
freqFile <- paste(prefix, "frequencies_xy", sep = "")
pBefore <- par(mfrow=c(4,4), oma=c(0,3,6,0), mar = c(2,2,2,1)) 
freqD <- read.delim(freqFile, header = T, stringsAsFactors = F)
freq_methods <- c("Kraken-Dist", "Bracken-Dist", "MetaMap-EM-Dist", "MetaMap-U-Dist")
freq_levels <- c("definedGenomes", "species", "genus", "family")
freqD[["freqTarget"]] <- as.numeric(freqD[["freqTarget"]])
freqD[["freqIs"]] <- as.numeric(freqD[["freqIs"]])
idx_minmax <- which(freqD[["level"]] %in% freq_levels)
pC <- list()
pC[["Unclassified"]] <- "#00FF00"
# pC[["NotLabelledAtLevel"]] <- "#009900"
# pC[["Mappable"]] <- barplotPal[[1]]
# pC[["species"]] <- barplotPal[[2]]
# pC[["genus"]] <- barplotPal[[3]]
# pC[["family"]] <- barplotPal[[4]]
# pC[["Other"]] <- barplotPal[[5]]
pC[["Mappable"]] <- "blue"
pC[["species"]] <- "blue"
pC[["genus"]] <- "blue"
pC[["family"]] <- "blue"
pC[["Other"]] <- "blue"
pC[["novelAdmixture"]] <- "orange"

pCex <- list()
pCex[["Unclassified"]] <- 2
# pCex[["NotLabelledAtLevel"]] <- 2
pCex[["Mappable"]] <- 2
pCex[["species"]] <- 2
pCex[["genus"]] <- 2
pCex[["family"]] <- 2
pCex[["Other"]] <- 2


pS <- list()
pS[["Unclassified"]] <- 19
# pS[["NotLabelledAtLevel"]] <- 19
pS[["Mappable"]] <- 3
pS[["species"]] <- 4
pS[["genus"]] <- 5
pS[["family"]] <- 6
pS[["Other"]] <- 7


pL <- list()
pL[["Unclassified"]] <- "Unclassified"
# pL[["NotLabelledAtLevel"]] <- "Unlabelled"
pL[["Mappable"]] <- "Genome"
pL[["species"]] <- "Species"
pL[["genus"]] <- "Genus"
pL[["family"]] <- "Family"
pL[["Other"]] <- ">Family"

frequency_plotTypes <-  c()
if(any(freqD[["variety"]] == "fullDB"))
{
	frequency_plotTypes <- c(frequency_plotTypes, "complete")
}
if(any(freqD[["variety"]] != "fullDB"))
{
	frequency_plotTypes <- c(frequency_plotTypes, "incompleteSummary")
}
varieties_nonFullDB <- sort(unique(freqD[["variety"]][freqD[["variety"]] != "fullDB"]))
frequency_plotTypes <- c(frequency_plotTypes, varieties_nonFullDB)
for(plotType in frequency_plotTypes)
{
	frequency_varities <- c("Summary")

	for(l in freq_levels)
	{
		xAxis_values <- c()
		yAxis_values <- c()
		for(iteration in c("defineAxes", "realDeal"))
		{
			for(m in freq_methods)
			{
				l_lookup <- l
				if((l == "definedGenomes") && (m == "MetaMap-U-Dist"))
				{
					l_lookup <- "definedAndHypotheticalGenomes"
				}

				indices <- c()
				if(plotType == "complete")
				{
					pC[["Mappable"]] <- "blue"			
					indices <- which((freqD[["method"]] == m) & (freqD[["level"]] == l_lookup) & (freqD[["variety"]] == "fullDB"))
				}
				else if(plotType == "incompleteSummary")
				{
					pC[["Mappable"]] <- "gray"			
					indices <- which((freqD[["method"]] == m) & (freqD[["level"]] == l_lookup) & (freqD[["variety"]] != "fullDB"))
				}
				else
				{
					pC[["Mappable"]] <- "gray"			
					indices <- which((freqD[["method"]] == m) & (freqD[["level"]] == l_lookup) & (freqD[["variety"]] == plotType))				
				}
				
				if(iteration == "defineAxes")
				{
					xAxis_values <- c(xAxis_values, freqD[["freqTarget"]][indices])
					yAxis_values <- c(yAxis_values, freqD[["freqIs"]][indices])
				}
				else
				{
					xAxis_min <- min(c(xAxis_values, yAxis_values))
					xAxis_max <- max(c(xAxis_values, yAxis_values))
					xAxis_half <- mean(c(xAxis_min, xAxis_max))
					xAxis_twoThirds <- xAxis_min + (2/5) * (xAxis_max - xAxis_min)

					yAxis_min <- min(c(xAxis_values, yAxis_values))
					yAxis_max <- max(c(xAxis_values, yAxis_values))
					yAxis_half <- mean(c(yAxis_min, yAxis_max))
					yAxis_twoThirds <- yAxis_min + (4/5) * (yAxis_max - yAxis_min)
					yAxis_twoThirds_2 <- yAxis_min + (3/5) * (yAxis_max - yAxis_min)
							
					# cat(paste(m, "@", l, ", ", "fullDB", ", r = ", cor(frTarget, frIs), sep = ""), "\n")
					
					pointColours <- c()
					pointCexs <- c()
					pointSymbols <- c()
					for(i in indices)
					{
						taxonID <- freqD[["taxonID"]][[i]]
						pointLevel <- freqD[["taxonIDCategory"]][[i]]
						mappable <- freqD[["isMappable"]][[i]]
						proportionNovelTotal <- freqD[["proportionNovelTotal"]][[i]]
						pointColour <- ""
						pointCex <- 0
						pointSymbol <- 0
						
						if(taxonID %in% names(pC))
						{
							pointColour <- pC[[taxonID]]
							pointCex <- pCex[[taxonID]]
							pointSymbol <- pS[[taxonID]]
						}
						else if(mappable == 1)
						{
							pointColour <- pC[["Mappable"]]
							pointCex <- pCex[["Mappable"]]
							pointSymbol <- pS[["Mappable"]]
						}
						else if(pointLevel %in% names(pC))
						{
							pointColour <- pC[[pointLevel]]
							pointCex <- pCex[[pointLevel]]
							pointSymbol <- pS[[pointLevel]]
						}
						else
						{
							pointColour <- pC[["Other"]]
							pointCex <- pCex[["Other"]]
							pointSymbol <- pS[["Other"]]
							cat("Unknown point level: ", pointLevel, "\n")
						}
						
						if(proportionNovelTotal > 0)
						{
							pointColour <- pC[["novelAdmixture"]]				
						}				

						pointColours <- c(pointColours, pointColour)
						pointCexs <- c(pointCexs, pointCex)
						pointSymbols <- c(pointSymbols, pointSymbol)
					}
						
					if(!((l == "definedGenomes") && ((m == "Kraken-Dist") || (m == "Bracken-Dist"))))
					{
						plot(0, 0, col = "white", main = "", cex.main = 0.7, xaxt = "n", yaxt = "n", xlim = c(xAxis_min, xAxis_max), ylim = c(yAxis_min, yAxis_max))		
						frTarget <- freqD[["freqTarget"]][indices]
						frIs <- freqD[["freqIs"]][indices]
						taxonIDs <- freqD[["taxonID"]][indices]
						taxonLabels <- freqD[["taxonLabel"]][indices]		
						stopifnot(length(indices) > 0)				
						points(frTarget, frIs,  pch = pointSymbols, cex = pointCexs, col = pointColours)
						stopifnot(length(frTarget) == length(pointColours))
						r <- cor(frTarget, frIs)
						L1 <- sum(abs(frTarget - frIs))
						text(xAxis_twoThirds, yAxis_twoThirds, labels = paste("r**2 = ", sprintf("%.3f", r**2), sep = ""), cex = 1.5)
						text(xAxis_twoThirds, yAxis_twoThirds_2, labels = paste("L1 = ", sprintf("%.3f", L1), sep = ""), cex = 1.5)
						axis(2, xaxs = "i", labels = F)			
						axis(1, xaxs = "i", labels = F)
						lines(c(xAxis_min, xAxis_max), c(yAxis_min, yAxis_max), col = "gray")
					}	
					else
					{
						plot(0, 0, col = "white", main = "", cex.main = 0.7, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlim = c(xAxis_min, xAxis_max), ylim = c(yAxis_min, yAxis_max), axes = F)
						if(m == "Kraken-Dist")
						{
							legend("topright", legend = c(pL, recursive = T, use.names = F),  col = c(pC, recursive = T, use.names = F), pch = as.vector(pS), cex = 1.5)
						}
					}
					
					if(m == freq_methods[[1]])
					{
						axis(2, xaxs = "i")		
						axis(2, at = xAxis_half, tick = F, line = 1, labels = c("Inferred"), cex.axis = 1)				
						axis(2, at = yAxis_half, tick = F, line = 2.5, labels = c(capitalize(l)), cex.axis = 1.8)	
					}
					else
					{
					}
					
					if(!((l == "definedGenomes") && ((m == "Kraken-Dist") || (m == "Bracken-Dist"))))
					{
						axis(1, xaxs = "i")			
						axis(1, at = xAxis_half, tick = F, line = 1, labels = c("Truth"), cex.axis = 1)	
					}
					else
					{
					}
					
					if(l == freq_levels[[1]])
					{
						title(main = m, cex.main = 1.8, line = 1)		
					}
				}
				# frTarget_idx_005 <- which(frIs <= 0.05)
				# plot(frTarget[frTarget_idx_005], frIs[frTarget_idx_005], main = paste("[Zoom] Frequencies ", m, ", ", l, ", ", v, sep = ""), xlab = "True frequency", ylab = "Inferred frequency", cex.main = 0.7)
			}
		}
	}
	if(plotType == "complete")
	{
		mtext("Inference with complete DB", side = 3, line = 2.2, outer = TRUE, cex = 1.5)
	}
	else if(plotType == "incompleteSummary")
	{
		mtext("Inference with incomplete DB (summary)", side = 3, line = 2.2, outer = TRUE, cex = 1.5)	
	}
	else
	{
		mtext(paste("Inference for ", plotType, sep = ""), side = 3, line = 2.2, outer = TRUE, cex = 1.5)		
	}				
}

dev.off()

par(pBefore)
# q()

# par(mfrow=c(1,1)) 


pdf(paste(simulationsDirectory, "/plots_reads.pdf", sep = ""))


barPlotFile <- paste(prefix, "barplots_fullDB", sep = "")
barplotD <- read.delim(barPlotFile, header = T)
#pdf(paste(prefix, "plots.pdf", sep = ""))
for(rL in unique(barplotD[["readLevel"]]))
{
	for(v in unique(barplotD[["variety"]]))
	{
		if((rL == "ALL") && (v == "fullDB"))	
		{
			barPlotVector <- c()
			callRatesVector <- c()
			methodNames <- unique(barplotD[["method"]])
			for(m in methodNames)
			{
				barPlotVector_thisMethod <- c()

				for(l in c("absolute", "species", "genus", "family"))
				{
					indices <- which((barplotD[["readLevel"]] == rL) & (barplotD[["variety"]] == v) & (barplotD[["level"]] == l) & (barplotD[["method"]] == m))
					stopifnot(length(indices) == 1)

					callRate <- barplotD[["callRate"]][indices[[1]]]
					accuracy <- barplotD[["accuracy"]][indices[[1]]]
					barPlotVector_thisMethod <- c(barPlotVector_thisMethod, accuracy)
					if(l == "absolute")
					{
						callRatesVector <- c(callRatesVector, callRate)
					}						
				}
				
				runningSum <- barPlotVector_thisMethod[[1]]
				for(i in 2:length(barPlotVector_thisMethod))
				{
					stopifnot(runningSum <= barPlotVector_thisMethod[[i]])
					runningSum_new <- barPlotVector_thisMethod[[i]]										
					barPlotVector_thisMethod[[i]] <- barPlotVector_thisMethod[[i]] - runningSum
					runningSum <- runningSum_new
				}
				barPlotVector <- c(barPlotVector, barPlotVector_thisMethod)
			}
			#print(barPlotVector)
			matrix_for_barplot <- matrix(barPlotVector, nrow = 4, ncol = length(methodNames))
			colnames(matrix_for_barplot) <- methodNames
			bpPos <- barplot(matrix_for_barplot, col = barplotPal, main = paste("Reads - complete DB (", rL, " reads against ", v, ")", sep = ""), ylab = "Accuracy", ylim = c(0, 1.3), axes = F)
			axis(1, at = bpPos, tick = F, labels = methodNames)			
			axis(2, at = c(0, 0.2, 0., 0.6, 0.8, 1), tick = T)

			points(bpPos, rep(1.2, length(bpPos)), cex = 5 * callRatesVector, pch = 19, col = "#777777")
			text(bpPos, rep(1.1, length(bpPos)), labels = paste(sprintf("%.f", 100*callRatesVector), "%", sep = ""), adj = c(0.5, 0.5), cex = 1)
			
			legend("bottomright", legend = rev(c("Genome", "Species", "Genus", "Family")), fill = rev(barplotPal[1:4]))			

			axis(2, at = c(1.2), tick = F, labels = c("Call rate"), cex.axis = 1, line = -1)
			
		}
	}
}
# dev.off()


byReadLengthFile <- paste(prefix, "byReadLength_fullDB", sep = "")
byReadLengthD <- read.delim(byReadLengthFile, header = T, stringsAsFactors = F)
stopifnot(all(byReadLengthD[["variety"]] == "fullDB"))
methodNames_reads <- unique(byReadLengthD[["method"]])
if(length(methodNames_reads) > 0)
{
	rL_l_min <- min(byReadLengthD[["readLength"]])
	rL_l_max <- max(byReadLengthD[["readLength"]])
	rL_accuracy_min <- min(byReadLengthD[["accuracyAvg"]])
	rL_accuracy_max <- max(byReadLengthD[["accuracyAvg"]])
	for(mI in 1:length(methodNames_reads))
	{
		m <- methodNames_reads[[mI]]
		plot(0, 0, col = "white", main = paste("Reads assignment v/s read length: ", m, sep = ""), xlab = "Read length", xlim = c(rL_l_min, rL_l_max), ylim = c(0, 1), ylab = "")
		plotLevels <- c("absolute", "species", "genus", "family")
		for(lI in 1:length(plotLevels))
		{
			l <- plotLevels[[lI]]
			indices <- which((byReadLengthD[["method"]] == m) & (byReadLengthD[["readCategory"]] == "ALL") & (byReadLengthD[["evaluationLevel"]] == l))	
			indices_validAccuracy <- which((byReadLengthD[["method"]] == m) & (byReadLengthD[["readCategory"]] == "ALL") & (byReadLengthD[["evaluationLevel"]] == l) & (byReadLengthD[["accuracyAvg"]] >= 0))	
			stopifnot(length(indices) > 0)
			stopifnot(length(indices_validAccuracy) > 0)
			if(lI == 1)
			{
				points(byReadLengthD[["readLength"]][indices], byReadLengthD[["callRateAvg"]][indices], col = "gray", pch = 2)
				lines(byReadLengthD[["readLength"]][indices], byReadLengthD[["callRateAvg"]][indices], col = "gray")
			}
			
			points(byReadLengthD[["readLength"]][indices_validAccuracy], byReadLengthD[["accuracyAvg"]][indices_validAccuracy], col = barplotPal[[lI]], pch = 19)		
		}
		legend("bottomright", legend = rev(c("Accuracy @ Genome", "Accuracy @ Species", "Accuracy @ Genus", "Accuracy @ Family", "Call rate")),  col = rev(c(barplotPal[1:4], "gray")), lty = c(0, 0, 0, 0, 1), pch = rev(c(19, 19, 19, 19, 2)))
		
	}
}

barPlotFile_byLevel <- paste(prefix, "barplots_readCategory", sep = "")
barplotD_byLevel <- read.delim(barPlotFile_byLevel, header = T, stringsAsFactors = F)
# pdf(paste(prefix, "plots.pdf", sep = ""))

attachedToFile_byLevel <- paste(prefix, "barplots_attachedTo", sep = "")
attachedToD <- read.delim(attachedToFile_byLevel, header = T, stringsAsFactors = F)

rCs <- c("truthLeafInDB", "novel_to_species", "novel_to_genus", "novel_to_superkingdom")
rCs <- rCs[rCs %in% barplotD_byLevel[["readCategory"]]]
rCs_labels <- list()
rCs_labels[["truthLeafInDB"]] <- "Read source in DB"
rCs_labels[["novel_to_species"]] <- "Reads to novel species node"
rCs_labels[["novel_to_genus"]] <- "Reads to novel genus node"
rCs_labels[["novel_to_superkingdom"]] <- "Reads to novel superkingdom node"

rCs_to_evaluationLevels <- list()
rCs_to_evaluationLevels[["truthLeafInDB"]] <- ""
rCs_to_evaluationLevels[["novel_to_species"]] <- "species"
rCs_to_evaluationLevels[["novel_to_genus"]] <- "genus"
rCs_to_evaluationLevels[["novel_to_superkingdom"]] <- "superkingdom"

evaluationLevels <- unique(barplotD_byLevel[["evaluationLevel"]])
evaluationLevels_ordered <- c("species", "genus", "family", "superkingdom")
evaluationLevels_ordered_assignment_names <- list()
evaluationLevels_ordered_assignment_names[["species"]] <- "Species/Strain"
evaluationLevels_ordered_assignment_names[["genus"]] <- "Genus"
evaluationLevels_ordered_assignment_names[["family"]] <- "Family"
evaluationLevels_ordered_assignment_names[["superkingdom"]] <- "Superkingdom"
evaluationLevels_to_i <- list()
for(i in 1:length(evaluationLevels_ordered))
{
	evaluationLevels_to_i[[evaluationLevels_ordered[[i]]]] <- i
}
barplotPal_byLevel <- barplotPal[2:5]

for(rC in rCs)
{
	stopifnot(rC %in% labels(rCs_labels))
	stopifnot(rC %in% labels(rCs_to_evaluationLevels))
	
	iMs <- sort(unique(barplotD_byLevel[["method"]][which(barplotD_byLevel[["readCategory"]] == rC)]))
	callRates_by_Method <- c()
	callRates_by_Method_list <- list()
	accuracies_by_Method_byLevel <- list()
	for(eL in evaluationLevels)
	{
		accuracies_by_Method_byLevel[[eL]] <- list()
	}
	for(iM in iMs)
	{
		callRates_iM <- c()
		for(eL in evaluationLevels)
		{
			relevantIndices <- which((barplotD_byLevel[["readCategory"]] == rC) & (barplotD_byLevel[["method"]] == iM) & (barplotD_byLevel[["evaluationLevel"]] == eL))
			stopifnot(length(relevantIndices) == 1)
			callRates_iM <- c(callRates_iM, barplotD_byLevel[["callRateAvg"]][relevantIndices[[1]]])
			accuracies_by_Method_byLevel[[iM]][[eL]] <- barplotD_byLevel[["accuracyAvg"]][relevantIndices[[1]]]
		}
		callRate <- callRates_iM[[1]]
		if(!all(abs(callRates_iM - callRate) <= 1e-5))
		{
			print("Call rate mismatch")
			print(callRate)
			print(callRates_iM)
		}
		stopifnot(all(abs(callRates_iM - callRate) <= 1e-5))
		callRates_by_Method <- c(callRates_by_Method, callRate)
		callRates_by_Method_list[[iM]] <- callRate
	}
	
	vector_for_barplot <- c()
	for(eL in c("absolute", evaluationLevels_ordered))
	{
		stopifnot(length(iMs) > 0)
		for(iM in iMs)
		{
			stopifnot(eL %in% names(accuracies_by_Method_byLevel[[iM]]))
			vector_for_barplot <- c(vector_for_barplot, accuracies_by_Method_byLevel[[iM]][[eL]])
			#cat(eL, " ", iM, " ", length(accuracies_by_Method_byLevel[[iM]][[eL]]), "\n")
		}
	}
	#print(vector_for_barplot)
	stopifnot(length(vector_for_barplot) > 1)
	matrix_for_barplot <- matrix(vector_for_barplot, ncol = 1 + length(evaluationLevels_ordered), nrow = length(iMs))
	# colnames(matrix_for_barplot) <- c("CR", evaluationLevels_ordered)
	colorVector <- c(rep(barplotPal[[1]], length(iMs)))
	for(i in 1:length(evaluationLevels_ordered))
	{
		thisEvaluationLevel <- evaluationLevels_ordered[[i]]
		thisNovelLevel <- rCs_to_evaluationLevels[[rC]]
		if(nchar(thisNovelLevel) > 0)
		{
			stopifnot(thisEvaluationLevel %in% names(evaluationLevels_to_i))
			stopifnot(thisNovelLevel %in% names(evaluationLevels_to_i))
			if(evaluationLevels_to_i[[thisEvaluationLevel]] < evaluationLevels_to_i[[thisNovelLevel]])
			{
				colorVector <- c(colorVector, rep("#EEEEEE", length(iMs)))
			}
			else
			{
				colorVector <- c(colorVector, rep(barplotPal_byLevel[[i]], length(iMs)))
			}
		}
		else
		{
			colorVector <- c(colorVector, rep(barplotPal_byLevel[[i]], length(iMs)))
		}
	}
	
	pMar <- par()$mar
	par(mar = pMar + c(4,0,0,0)) 
	bpPos <- barplot(matrix_for_barplot, beside = T, main = paste("Reads - all experiments (compl./incompl.): ", rCs_labels[[rC]]), col = colorVector, ylim = c(0, 1.45), axes = F)

	axis(1, at = colMeans(bpPos), tick = F, labels = c("Absolute", sapply(evaluationLevels_ordered, function(x){capitalize(x)}, USE.NAMES = F)), pos = 0)
	axis(1, at = bpPos, tick = F, labels = rep(iMs, 1 + length(evaluationLevels_ordered)), las = 2, pos = -0.1, cex.axis = 0.7)
	axis(2, at = c(0, 0.5, 1), tick = T)
	
	par(mar = pMar)
	
	attachedToProportions_methods <- c("Metamap-U-Reads", "Metamap-EM-Reads", "Kraken-Reads")

	plotCircles_x_names <- c("Call rate", sapply(evaluationLevels_ordered, function(x){paste("Assn.: ", evaluationLevels_ordered_assignment_names[[x]])}, USE.NAMES = F) )
	axis(3, at = colMeans(bpPos), tick = F, labels = plotCircles_x_names, cex.axis = 0.7, pos = 1.35)
	
	
	plotCircles_realized_y_values <- c()
	for(iMi in 1:length(attachedToProportions_methods))
	{
		iM <- attachedToProportions_methods[[iMi]]
		
		plotCircles_x_values <- c(callRates_by_Method_list[[iM]])
		for(eLi in 1:length(evaluationLevels_ordered))
		{
			thisEvaluationLevel <- evaluationLevels_ordered[[eLi]]
			indices <- which((attachedToD[["method"]] == iM) & (attachedToD[["readCategory"]] == rC))
			if(!(length(indices) == 1))
			{
				# print(indices)
				cat("method = ", iM, ", readCategory = ", rC, "\n")
			}
			stopifnot(length(indices) == 1)
			columnName_attachedTo <- paste("attachedTo_", thisEvaluationLevel, sep = "")
			if(!(columnName_attachedTo %in% names(attachedToD)))
			{
				cat("Missing column: ", columnName_attachedTo, "\n")
			}
			stopifnot(columnName_attachedTo %in% names(attachedToD))
			
			plotCircles_x_values <- c(plotCircles_x_values, attachedToD[[columnName_attachedTo]][[indices[[1]]]])
		}
		stopifnot(length(plotCircles_x_names) == length(plotCircles_x_values))
		
		points_x <- colMeans(bpPos)
		yPos <- 1.38 - (iMi - 1) * 0.1
		plotCircles_realized_y_values <- c(plotCircles_realized_y_values, yPos)
		points_y <- rep(yPos, length(points_x))
		
		points(points_x, points_y, pch = 19, cex = plotCircles_x_values * 1.9, col = c("black", barplotPal_byLevel))		
		
		text(points_x, points_y - 0.05, labels = paste(sprintf("%.f", 100*plotCircles_x_values), "%", sep = ""), adj = c(0.5, 0.5), cex = 0.5)
	}
	
	axis(2, at = plotCircles_realized_y_values, tick = F, labels = attachedToProportions_methods, las = 2, cex.axis = 0.7, line = -2.4)
	axis(2, at = c(0.5), tick = F, labels = c("Accuracy"), line = 2)

}

par(mfrow=c(1,1)) 

dev.off()
