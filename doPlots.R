library("RColorBrewer")

simulationsDirectory <- "databases/miniSeq_100/simulations_logNormal"

prefix <- paste(simulationsDirectory, "/_forPlot_", sep = "")

pdf(paste(simulationsDirectory, "/plots.pdf", sep = ""))


barPlotFile <- paste(prefix, "barplots_fullDB", sep = "")
barplotD <- read.delim(barPlotFile, header = T)
barplotPal <- brewer.pal(n = 5, name = "PuBu")
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

				for(l in c("mappingTarget", "species", "genus", "family"))
				{
					indices <- which((barplotD[["readLevel"]] == rL) & (barplotD[["variety"]] == v) & (barplotD[["level"]] == l) & (barplotD[["method"]] == m))
					stopifnot(length(indices) == 1)

					callRate <- barplotD[["callRate"]][indices[[1]]]
					accuracy <- barplotD[["accuracy"]][indices[[1]]]
					barPlotVector_thisMethod <- c(barPlotVector_thisMethod, accuracy)
					if(l == "mappingTarget")
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
			barplot(matrix_for_barplot, col = barplotPal, main = paste("Reads - complete DB (", rL, " reads against ", v, ")", sep = ""))
			legend("topright", legend = rev(c("Genome", "Species", "Genus", "Family")), fill = rev(barplotPal[1:4]))			
		}
	}
}
# dev.off()


barPlotFile_byLevel <- paste(prefix, "barplots_readCategory", sep = "")
barplotD_byLevel <- read.delim(barPlotFile_byLevel, header = T, stringsAsFactors = F)
# pdf(paste(prefix, "plots.pdf", sep = ""))

rCs <- c("truthLeafInDB", "novel_to_superkingdom")
rCs_labels <- list()
rCs_labels[["truthLeafInDB"]] <- "Read source in DB"
rCs_labels[["novel_to_superkingdom"]] <- "Reads to novel superkingdom node"

rCs_to_evaluationLevels <- list()
rCs_to_evaluationLevels[["truthLeafInDB"]] <- ""
rCs_to_evaluationLevels[["novel_to_superkingdom"]] <- "superkingdom"

evaluationLevels <- unique(barplotD_byLevel[["evaluationLevel"]])
evaluationLevels_ordered <- c("species", "genus", "family", "superkingdom")
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
		stopifnot(all(callRates_iM == callRate))
		callRates_by_Method <- c(callRates_by_Method, callRate)
	}
	
	vector_for_barplot <- callRates_by_Method
	for(eL in evaluationLevels_ordered)
	{
		for(iM in iMs)
		{
			stopifnot(eL %in% names(accuracies_by_Method_byLevel[[iM]]))
			vector_for_barplot <- c(vector_for_barplot, accuracies_by_Method_byLevel[[iM]][[eL]])
		}
	}
	matrix_for_barplot <- matrix(vector_for_barplot, ncol = 1 + length(evaluationLevels_ordered), nrow = length(iMs))
	# colnames(matrix_for_barplot) <- c("CR", evaluationLevels_ordered)
	colorVector <- c(rep("gray", length(iMs)))
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
	bpPos <- barplot(matrix_for_barplot, beside = T, main = paste("Reads - incomplete DB: ", rCs_labels[[rC]]), col = colorVector)
	axis(1, at = colMeans(bpPos), tick = F, labels = c("Read call rate", evaluationLevels_ordered), pos = 0)
	axis(1, at = bpPos, tick = F, labels = rep(iMs, 1 + length(evaluationLevels_ordered)), las = 2, pos = -0.1, cex.axis = 0.7)
	par(mar = pMar)
}

freqFile <- paste(prefix, "frequencies_xy", sep = "")
freqD <- read.delim(freqFile, header = T)
for(v in unique(freqD[["variety"]]))
{
	for(l in unique(freqD[["level"]]))
	{
		for(m in sort(unique(freqD[["method"]])))
		{
			indices <- which((freqD[["variety"]] == v) & (freqD[["level"]] == l) & (freqD[["method"]] == m))
			if(length(indices) > 0)
			{
				frTarget <- freqD[["freqTarget"]][indices]
				frIs <- freqD[["freqIs"]][indices]
				taxonIDs <- freqD[["taxonID"]][indices]
				taxonLabels <- freqD[["taxonLabel"]][indices]
				cat(paste(m, "@", l, ", ", v, ", r = ", cor(frTarget, frIs), sep = ""), "\n")
				if((v == "fullDB") && (((m == "MetaMap-EM-Dist") && (l == "definedGenomes")) || ((m == "MetaMap-U-Dist") && (l == "definedAndHypotheticalGenomes"))))
				{
					plot(frTarget, frIs, main = paste("Frequencies ", m, ", ", l, ", ", v, sep = ""), xlab = "True frequency", ylab = "Inferred frequency", cex.main = 0.7)
					frTarget_idx_005 <- which(frIs <= 0.05)
					plot(frTarget[frTarget_idx_005], frIs[frTarget_idx_005], main = paste("[Zoom] Frequencies ", m, ", ", l, ", ", v, sep = ""), xlab = "True frequency", ylab = "Inferred frequency", cex.main = 0.7)
				}
			}
		}
	}
}
par(mfrow=c(1,1)) 

dev.off()

