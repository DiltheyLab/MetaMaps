pdf("HMPplot.pdf", width = 10, height = 5)
sources_nonTruth_order <- c("Kraken", "Bracken", "MetaMap-EM", "MetaMap-U")
#for(v in c("PacBio", "Nanopore"))
for(v in c("PacBio"))
{
	D <- read.delim(paste("_HMP_distributions_", v, ".txt", sep = ""), stringsAsFactors = F)

	sources <- unique(D[["Source"]])
	stopifnot("truth" %in% sources)
	sources_nonTruth <- sources[sources != "truth"]
	stopifnot(all(sources_nonTruth_order %in% sources_nonTruth))
	stopifnot(all(sources_nonTruth %in% sources_nonTruth_order))
	sources_nonTruth <- sources_nonTruth_order
	for(analysisLevel in c("species", "genus"))
	{
		D_level <- D[D[["Level"]] == analysisLevel,]
		taxonIDs <- unique(D_level[["taxonID"]])
		taxonIDs_labels <- c()
		for(tI in taxonIDs)
		{
			i <- which((D_level[["Source"]] == "truth") & (D_level[["taxonID"]] == tI))
			stopifnot(length(i) == 1)	
			tI_label <- D_level[["taxonLabel"]][[i]]
			taxonIDs_labels <- c(taxonIDs_labels, tI_label)
		}
		
		taxonIDs <- taxonIDs[order(taxonIDs_labels)]
		
		taxonIDs_true <- c()
		taxonID_true_labels <- c()
		for(tI in taxonIDs)
		{
			i <- which((D_level[["Source"]] == "truth") & (D_level[["taxonID"]] == tI))
			stopifnot(length(i) == 1) 
			if(D_level[["F"]][[i]] != 0)
			{
				taxonIDs_true <- c(taxonIDs_true, tI)
				tLabel <- D_level[["taxonLabel"]][[i]]
				taxonID_true_labels <- c(taxonID_true_labels, tLabel)
			}
		}	
		
		
		legend_for_plot <- c()
		colors_for_legend <- c("blue", "lightblue", "cyan", "red", "pink")
		for(S in c("GoldStandard", sources_nonTruth))
		{
			legend_for_plot <- c(legend_for_plot, S)
		}
		stopifnot(length(legend_for_plot) == length(colors_for_legend))
		
		realizedSums <- list()
		for(S in c("truth", sources_nonTruth))
		{
			if(!(abs(1 - sum(D_level[["F"]][(D_level[["Source"]] == S)])) <= 1e-5))
			{
				cat("Error - truth doesn't add up\n")
				cat("Source ", S, "\n")
				cat(analysisLevel, "\n")
				cat("sum(D_level[[F]][(D_level[[Source]] == S)]) = ", sum(D_level[["F"]][(D_level[["Source"]] == S)]), "\n")
			}
			stopifnot(abs(1 - sum(D_level[["F"]][(D_level[["Source"]] == S)])) <= 1e-5)
			realizedSums[[S]] <- 0
		}
		
		vector_for_plot <- c()	
		for(tI in taxonIDs_true)
		{
			for(S in c("truth", sources_nonTruth))
			{
				i <- which((D_level[["Source"]] == S) & (D_level[["taxonID"]] == tI))
				stopifnot(length(i) == 1)
				f <- D_level[["F"]][[i]]
				vector_for_plot <- c(vector_for_plot, f)
				realizedSums[[S]] <- realizedSums[[S]] + f
			}
		}	
		
		for(S in c("truth", sources_nonTruth))
		{
			missingS <- 1 - realizedSums[[S]]
			vector_for_plot <- c(vector_for_plot, missingS)
		}
		
		
		matrix_for_barplot <- matrix(vector_for_plot, nrow = 1 + length(sources_nonTruth))
		colnames(matrix_for_barplot) <- 1:(length(taxonIDs_true)+1)

		barplot(matrix_for_barplot, beside = T, col = colors_for_legend, las=2, main = paste(v, ": Abundance estimates [", analysisLevel, "]", sep = ""), cex.axis = 1.4, cex.names = 1.7) 
		
		legend("topright", legend = legend_for_plot, fill = colors_for_legend, cex = 1.0)	
		
		cat(paste("Labels ", v, " @ ", analysisLevel, "\n", sep = ""))
		for(lI in 1:length(taxonID_true_labels))
		{
			cat(paste("\t", lI, " ", taxonID_true_labels[[lI]], "\n", sep = ""))
		}
		cat(paste("\t", length(taxonID_true_labels)+1, " ", "Other", "\n", sep = ""))
		
		xyPlot_max <- max(D_level[["F"]])
		
		
		plot(0, 0, xlim = c(0, xyPlot_max), ylim = c(0, xyPlot_max), main = paste(v, ": Abundance estimates [", analysisLevel, "]", sep = ""), col = "white", xlab = "Truth", ylab = "Estimation")

		labels_for_xyLegend <- c()
		for(Si in 1:length(sources_nonTruth))
		{
			S <- sources_nonTruth[[Si]]
			x_values <- c()
			y_values <- c()
			for(taxonID in taxonIDs)
			{
				x_value_i <- which((D_level[["Source"]] == "truth") & (D_level[["taxonID"]] == taxonID))
				y_value_i <- which((D_level[["Source"]] == S)  & (D_level[["taxonID"]] == taxonID))
				stopifnot(length(x_value_i) == 1)
				stopifnot(length(y_value_i) == 1)
				x_values <- c(x_values, D_level[["F"]][[x_value_i]])
				y_values <- c(y_values, D_level[["F"]][[y_value_i]])
			
			}					
			points(x_values, y_values, col = colors_for_legend[[Si+1]])
			
			r <- sprintf("%.3f", cor(x_values, y_values))
			# labels_for_xyLegend <- c(labels_for_xyLegend, paste(sources_nonTruth[[Si]], " - r = ", r, sep = ""))
			labels_for_xyLegend <- c(labels_for_xyLegend, paste(sources_nonTruth[[Si]], sep = ""))
		}
		lines(c(0, xyPlot_max), c(0, xyPlot_max), col = "gray")
		
		legend("topright", legend = labels_for_xyLegend, fill = colors_for_legend[2:length(colors_for_legend)])
		# cat(paste(c(paste("\n", v, " labels: ", analysisLevel, sep = ""), taxonID_true_labels, "", recursive = T), collapse = "\n"))
	}

}

dev.off()
