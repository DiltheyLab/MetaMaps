pdf("HMP_and_Zymo_plot.pdf", width = 10, height = 10)
par(mfrow=c(2,1))
#layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE), widths=c(18,11))
  
sources_nonTruth_order <- c("Kraken", "Bracken", "MetaMap-EM", "MetaMap-U")
sources_nonTruth_order <- c("Kraken", "Bracken", "MetaMap-EM")
sources_nonTruth_order <- c("MetaMap-EM", "Bracken", "Centrifuge")
sources_nonTruth_order <- c("MetaMap-EM", "Bracken", "Centrifuge", "Megan")

target_width <- 19
for(what in c("HMP", "Zymo"))
{
	#for(v in c("PacBio", "Nanopore"))
	for(v in c("PacBio"))
	{
		stopifnot(v == "PacBio")
		#D <- read.delim(paste("_HMP_distributions_", v, ".txt", sep = ""), stringsAsFactors = F)
		D <- read.delim(paste("_HMP_distributions_", what, ".txt", sep = ""), stringsAsFactors = F)

		sources <- unique(D[["Source"]])
		stopifnot("truth" %in% sources)
		sources_nonTruth <- sources[sources != "truth"]
		if(!all(sources_nonTruth_order %in% sources_nonTruth))
		{
			print(unique(sources_nonTruth))
		}
		stopifnot(all(sources_nonTruth_order %in% sources_nonTruth))
		# stopifnot(all(sources_nonTruth %in% sources_nonTruth_order))
		sources_nonTruth <- sources_nonTruth_order
		#for(analysisLevel in c("species", "genus"))
		for(analysisLevel in c("genus"))
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
			colors_for_legend <- c("blue", "lightblue", "cyan", "red")
			colors_for_legend <- c("blue", "red", "lightblue")
			colors_for_legend <- c("gray", "blue", "firebrick2", "orange", "lightpink3")
			for(S in c("Truth", sources_nonTruth))
			{
				if(S == "MetaMap-EM")
				{ 
					S <- "MetaMaps"
				}	
				if(S == "Megan")
				{ 
					S <- "MEGAN-LR"
				}				
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
			L1_by_method <- list()	
			for(m in sources_nonTruth)
			{
				L1_by_method[[m]] <- 0
			}

			for(tI in taxonIDs_true)
			{
				for(S in c("truth", sources_nonTruth))
				{
					i <- which((D_level[["Source"]] == S) & (D_level[["taxonID"]] == tI))
					stopifnot(length(i) == 1)
					f <- D_level[["F"]][[i]]
					vector_for_plot <- c(vector_for_plot, f)
					realizedSums[[S]] <- realizedSums[[S]] + f
					if(S != "truth")
					{
						i_truth <- which((D_level[["Source"]] == "truth") & (D_level[["taxonID"]] == tI))
						f_truth <- D_level[["F"]][[i_truth]]
						L1_by_method[[S]] <- L1_by_method[[S]] + abs(f_truth - f)
					}
				}
			}	
			

			for(S in c("truth", sources_nonTruth))
			{
				missingS <- 1 - realizedSums[[S]]
				vector_for_plot <- c(vector_for_plot, missingS)
				if(S != "truth")
				{
					L1_by_method[[S]] <- L1_by_method[[S]] + missingS
				}
			}

			vector_for_plot <- c(vector_for_plot, 0)
			for(S in sources_nonTruth)
			{
				vector_for_plot <- c(vector_for_plot, L1_by_method[[S]])
			}

			missing_entries_horizontal <- target_width - (length(taxonIDs_true) + 2)
			stopifnot(missing_entries_horizontal >= 0)
			print(missing_entries_horizontal)
			
			colnames_for_barplot <- c(1:(length(taxonIDs_true)+1), "L1")
			
			if(missing_entries_horizontal >= 1)
			{
				for(i in 1:missing_entries_horizontal)
				{
					for(S in c("truth", sources_nonTruth))
					{
						vector_for_plot <- c(vector_for_plot, 0)				
					}
					colnames_for_barplot <- c(colnames_for_barplot, "")
				}
			}

			

			matrix_for_barplot <- matrix(vector_for_plot, nrow = 1 + length(sources_nonTruth))
			matrix_for_barplot_old <- matrix_for_barplot
			colnames(matrix_for_barplot) <- colnames_for_barplot

			mainT <- paste(v, ": Abundance estimates [", analysisLevel, "]", sep = "")
			mainT <- paste(what, " composition", sep = "")
			values_in_plot <- sort(as.vector(matrix_for_barplot), decreasing = T)
			thirdBiggest <- values_in_plot[[3]]
			index_values_bigger <- which(matrix_for_barplot > thirdBiggest)
			#print(c("thirsBiggest", thirdBiggest, length(index_values_bigger)))
			#print(index_values_bigger)
			for(biggerIdx in index_values_bigger)
			{
				oldValue <- matrix_for_barplot[[biggerIdx]]
				newValue <- (oldValue - thirdBiggest) / 10 + thirdBiggest
				matrix_for_barplot[[biggerIdx]] <- newValue
			}
			barplot_return <- barplot(matrix_for_barplot, beside = T, col = colors_for_legend, las=2, main = mainT, cex.axis = 1.4, cex.names = 1.7, cex.main = 2)
			
			for(biggerIdx in index_values_bigger)
			{
				# print(barplot_return[[biggerIdx]])
				text(barplot_return[[biggerIdx]], matrix_for_barplot[[biggerIdx]], labels = sprintf("%.2f", matrix_for_barplot_old[[biggerIdx]]), adj = c(0.5, 0))
			}
			
			# print(barplot_return)
			
			if(what == "HMP")
			{
				legend("topleft", legend = legend_for_plot, fill = colors_for_legend, cex = 1.0)	
			}
			cat(paste(what, " labels ", v, " @ ", analysisLevel, "\n", sep = ""))
			for(lI in 1:length(taxonID_true_labels))
			{
				cat(paste("\t", lI, " ", taxonID_true_labels[[lI]], "\n", sep = ""))
			}
			cat(paste("\t", length(taxonID_true_labels)+1, " ", "Other", "\n", sep = ""))
			
			xyPlot_max <- max(D_level[["F"]])
			
			if(1 == 0)
			{
				
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
				
				legend("topright", legend = labels_for_xyLegend, fill = colors_for_legend[2:length(colors_for_legend)], cex = 1.5)
				# cat(paste(c(paste("\n", v, " labels: ", analysisLevel, sep = ""), taxonID_true_labels, "", recursive = T), collapse = "\n"))
			}
		}

	}
}

dev.off()
