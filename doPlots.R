prefix <- "_forPlot_"
freqFile <- paste(prefix, "frequencies_xy", sep = "")
freqD <- read.delim(freqFile, header = T)
pdf(paste(prefix, "plots.pdf", sep = ""))
for(v in unique(freqD[["variety"]]))
{
	for(l in unique(freqD[["level"]]))
	{
		for(m in unique(freqD[["method"]]))
		{
			indices <- which((freqD[["variety"]] == v) & (freqD[["level"]] == l) & (freqD[["method"]] == m))
			if(length(indices) > 0)
			{
				frTarget <- freqD[["freqTarget"]][indices]
				frIs <- freqD[["freqIs"]][indices]
				taxonIDs <- freqD[["taxonID"]][indices]
				taxonLabels <- freqD[["taxonLabel"]][indices]
				cat(paste(m, "@", l, ", ", v, ", r = ", cor(frTarget, frIs), sep = ""), "\n")
				if((v == "fullDB") && (l == "mappingTarget") && ((m == "MetaMap-EM-Dist") || (m == "MetaMap-U-Dist")))
				{
					plot(frTarget, frIs, main = paste(m, "@", l, ", ", v, sep = ""))
					frTarget_idx_005 <- which(frIs <= 0.05)
					plot(frTarget[frTarget_idx_005], frIs[frTarget_idx_005], main = paste(m, "@", l, ", ", v, sep = ""))
				}
			}
		}
	}
}

library("RColorBrewer")

barPlotFile <- paste(prefix, "barplots_fullDB", sep = "")
barplotD <- read.delim(barPlotFile, header = T)
barplotPal <- brewer.pal(n = 4, name = "PuBu")
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
			barplot(matrix_for_barplot, col = barplotPal, main = paste("Correctly classified: ", rL, " reads / ", v, sep = ""))
			legend("topright", legend = rev(c("Strain", "Species", "Genus", "Family")), fill = rev(barplotPal))			
		}
	}
}
dev.off()


barPlotFile_byLevel <- paste(prefix, "barplots_readCategory", sep = "")
barplotD_byLevel <- read.delim(barPlotFile_byLevel, header = T, stringsAsFactors = F)
pdf(paste(prefix, "plots.pdf", sep = ""))
rCs <- unique(barplotD_byLevel[["readCategory"]])
rCs_labels <- list()
rCs_labels[["truthLeafInDB"]] <- "Read source in DB"
evaluationLevels <- unique(barplotD_byLevel[["evaluationLevel"]])
evaluationLevels_ordered <- c("mappingTarget", "species", "genus")
barplotPal_byLevel <- brewer.pal(n = length(evaluationLevels_ordered), name = "PuBu")
for(rC in c("truthLeafInDB"))
{
	iMs <- unique(barplotD_byLevel[["method"]][which(barplotD_byLevel[["readCategory"]] == rC)])
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
		colorVector <- c(colorVector, rep(barplotPal_byLevel[[i]], length(iMs)))
	}
	bpPos <- barplot(matrix_for_barplot, beside = T, main = rCs_labels[[rC]], col = colorVector)
	axis(1, at = colMeans(bpPos), tick = F, labels = c("CR", evaluationLevels_ordered), pos = -1)
	
	axis(1, at = bpPos, tick = F, labels = rep(iMs, 1 + length(evaluationLevels_ordered)), las = 2)
	
	
}
