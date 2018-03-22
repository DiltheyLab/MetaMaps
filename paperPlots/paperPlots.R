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

barplotPal <- brewer.pal(n = 5, name = "PuBu")

captions <- list()
captions[["Kraken-Reads"]] <- "Kraken"
captions[["Metamap-EM-Reads"]] <- "MetaMap-Complete"
captions[["Metamap-U-Reads"]] <- "MetaMap-Unknown"

colourByMethod <- list()
colourByMethod[["Kraken-Reads"]] <- "darkred"
colourByMethod[["Metamap-EM-Reads"]] <- "darkblue"
colourByMethod[["Metamap-U-Reads"]] <- "lightblue"


lineStyleByLevel <- list()
lineStyleByLevel[["species"]] <- "solid"
lineStyleByLevel[["absolute"]] <- "dashed"

dotStyleByLevel <- list()
dotStyleByLevel[["species"]] <- 1
dotStyleByLevel[["absolute"]] <- 2


translateCaptions <- function(n)
{
	stopifnot(n %in% names(captions))
	captions[[n]]
}

readLengthPlot <- function(d1, t1)
{
	pdf("readLengthPlot.pdf", width = 12, height = 6)
	par(mfrow=c(1,2)) 
	

	readLengthFile <- paste(d1, "/_forPlot_byReadLength_fullDB", sep = "")
	
	byReadLengthD <- read.delim(readLengthFile, header = T, stringsAsFactors = F)
	stopifnot(all(byReadLengthD[["variety"]] == "fullDB"))
	methodNames_reads <- c("Metamap-EM-Reads", "Metamap-U-Reads", "Kraken-Reads")

	rL_l_min <- min(byReadLengthD[["readLength"]])
	rL_l_max <- max(byReadLengthD[["readLength"]])
	rL_accuracy_min <- min(byReadLengthD[["accuracyAvg"]])
	rL_accuracy_max <- max(byReadLengthD[["accuracyAvg"]])
		
	plot(0, 0, col = "white", main = paste("Read assignment in ", t1), xlab = "Read length bin", ylab = "Proportion reads correctly assigned", xlim = c(rL_l_min, rL_l_max), ylim = c(0, 1))
	
	legend_titles <- c()
	legend_colours <- c()
	legend_symbols <- c()
	legend_lineStyles <- c()
	
	if(length(methodNames_reads) > 0)
	{
		for(mI in 1:length(methodNames_reads))
		{
			m <- methodNames_reads[[mI]]
			plotLevels <- c("absolute", "species")
			if(m == "Kraken-Reads")
			{
				plotLevels <- c("species")
			}
			for(lI in 1:length(plotLevels))
			{
				l <- plotLevels[[lI]]
				indices <- which((byReadLengthD[["method"]] == m) & (byReadLengthD[["readCategory"]] == "ALL") & (byReadLengthD[["evaluationLevel"]] == l))	
				indices_validAccuracy <- which((byReadLengthD[["method"]] == m) & (byReadLengthD[["readCategory"]] == "ALL") & (byReadLengthD[["evaluationLevel"]] == l) & (byReadLengthD[["accuracyAvg"]] >= 0))	
				stopifnot(length(indices) > 0)
				stopifnot(length(indices_validAccuracy) > 0)
				if(lI == 1)
				{
					#points(byReadLengthD[["readLength"]][indices], byReadLengthD[["callRateAvg"]][indices], col = "gray", pch = 2)
					#lines(byReadLengthD[["readLength"]][indices], byReadLengthD[["callRateAvg"]][indices], col = "gray")
				}
				
				x_per_point <- byReadLengthD[["readLength"]][indices_validAccuracy]
				n_per_point <- byReadLengthD[["Ntotal"]][indices_validAccuracy]
				p_per_point <- byReadLengthD[["accuracyAvg"]][indices_validAccuracy]
				
				lines(x_per_point, p_per_point, col = colourByMethod[[m]], lty = lineStyleByLevel[[l]], lwd = 1)		
				points(byReadLengthD[["readLength"]][indices_validAccuracy], byReadLengthD[["accuracyAvg"]][indices_validAccuracy], col = colourByMethod[[m]], pch = dotStyleByLevel[[l]])		
				
						
				confidence_interval_size <- c()
				for(j in 1:length(p_per_point))
				{
					confidence_interval_size <- c(confidence_interval_size, qnorm(0.975)*sqrt((1/n_per_point[[j]])*p_per_point[[j]]*(1-p_per_point[[j]])))
					lines(rep(x_per_point[[j]], 2), p_per_point[[j]] + c(1, -1) * confidence_interval_size[[j]], col = colourByMethod[[m]])
				}
				
				legend_titles <- c(legend_titles, paste(translateCaptions(m), " ", l, sep = ""))
				legend_colours <- c(legend_colours, colourByMethod[[m]])
				legend_symbols <- c(legend_symbols, dotStyleByLevel[[l]])
				legend_lineStyles <- c(legend_lineStyles, lineStyleByLevel[[l]])				
			}
			
			# legend("bottomright", legend = rev(c("Accuracy @ Genome", "Accuracy @ Species", "Accuracy @ Genus", "Accuracy @ Family", "Call rate")),  col = rev(c(barplotPal[1:4], "gray")), lty = c(0, 0, 0, 0, 1), pch = rev(c(19, 19, 19, 19, 2)))
			
		}
		
		legend("bottomright", legend = legend_titles, col = legend_colours, lty = legend_lineStyles, pch = legend_symbols, bg = "white")			

								
	}
	
	dev.off()
	
}

twoReadPlots <- function(d1, t1, d2, t2)
{
	pdf("readAccuracyPlot.pdf", width = 12, height = 6)
	par(mfrow=c(1,2)) 
	dirs <- c(d1, d2)
	titles <- c(t1, t2)
	for(dI in 1:length(dirs)) # todo
	{
		barPlotFile <- paste(dirs[[dI]], "/_forPlot_barplots_fullDB", sep = "")
		cat("Reading ", barPlotFile, "\n")
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
					methodNames <- as.character(unique(barplotD[["method"]]))
					
					for(m in methodNames)
					{
						barPlotVector_thisMethod <- c()

						for(l in c("absolute", "species", "genus", "family"))
						{
							indices <- which((barplotD[["readLevel"]] == rL) & (barplotD[["variety"]] == v) & (barplotD[["level"]] == l) & (barplotD[["method"]] == m))
							# print(indices)
							stopifnot(length(indices) == 1)

							callRate <- barplotD[["callRate"]][indices[[1]]]
							accuracy <- barplotD[["accuracy"]][indices[[1]]]
							if((m == "Kraken-Reads") && (l == "absolute"))
							{
								accuracy <- 0
							}
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
					# colnames(matrix_for_barplot) <- methodNames
					plotTitle <- paste("Reads - complete DB (", rL, " reads against ", v, ")", sep = "")
					plotTitle <- titles[[dI]]
					bpPos <- barplot(matrix_for_barplot, col = barplotPal, main = plotTitle, ylab = "Accuracy", ylim = c(0, 1.3), axes = F)
					axis(1, at = bpPos, tick = F, labels = lapply(methodNames, translateCaptions), cex.axis = 0.8)			
					axis(2, at = c(0, 0.2, 0., 0.6, 0.8, 1), tick = T)

					points(bpPos, rep(1.2, length(bpPos)), cex = 5 * callRatesVector, pch = 19, col = "#777777")
					text(bpPos, rep(1.1, length(bpPos)), labels = paste(sprintf("%.f", 100*callRatesVector), "%", sep = ""), adj = c(0.5, 0.5), cex = 1)
					
					if(dI == 2)
					{
						legend("bottomright", legend = rev(c("Genome", "Species", "Genus", "Family")), fill = rev(barplotPal[1:4]))			
					}
					
					axis(2, at = c(1.2), tick = F, labels = c("Call rate"), cex.axis = 1, line = -1)
					
				}
			}
		}
	}
	
	dev.off()
	
	
}

readLengthPlot("../databases/miniSeq+H/simulations_i100_specifiedFrequencies", "i100")
# twoReadPlots("../databases/miniSeq+H/simulations_i100_specifiedFrequencies", "i100", "../databases/miniSeq+H/simulations_p25_logNormal", "p25")