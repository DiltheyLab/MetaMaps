D <- read.delim("_HMP_distributions.txt", stringsAsFactors = F)

pdf("HMPplot.pdf", width = 10, height = 5)
sources <- unique(D[["Source"]])
stopifnot("truth" %in% sources)
sources_nonTruth <- sources[sources != "truth"]
for(analysisLevel in c("species", "genus"))
{
	D_level <- D[D[["Level"]] == analysisLevel,]
	taxonIDs <- unique(D_level[["taxonID"]])
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

	barplot(matrix_for_barplot, beside = T, col = colors_for_legend, las=2, main = paste("Abundance estimates [", analysisLevel, "]", sep = ""), cex.axis = 1.4, cex.names = 1.7) 
	
	legend("topright", legend = legend_for_plot, fill = colors_for_legend, cex = 1.0)	
	
	cat(paste(c(paste("\nLabels: ", analysisLevel, sep = ""), taxonID_true_labels, "", recursive = T), collapse = "\n"))
}

dev.off()