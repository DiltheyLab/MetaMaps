args = commandArgs(trailingOnly=TRUE)
file_prefix <- ""
if(length(args) >= 1)
{
	stopifnot(length(args) >= 1)
	file_prefix <- args[[1]]
} else
{
	file_prefix <- "/scratch/tmp/MetaMap/hmp_set7"
	file_prefix <- "tmp/hmp7_2_miniSeq+H"
}


# Rscript plotUnknownResults.R /scratch/tmp/hmp_set7_run2/withAll.wimpAnalysis.unknown

file_out <- paste(file_prefix, ".U.plots.pdf", sep = "")
idty <- read.delim(paste(file_prefix, ".U.lengthAndIdentitiesPerTaxonID", sep = ""), stringsAsFactors = F, na.strings = "", header = T)
shifted <- read.delim(paste(file_prefix, ".U.shiftedHistogramsPerTaxonID", sep = ""), stringsAsFactors = F, na.strings = "", header = T)
potDistr <- read.delim(paste(file_prefix, ".U.WIMP", sep = ""), stringsAsFactors = F, na.strings = "", header = T)
potDistr <- potDistr[potDistr[[1]] == "definedGenomes",]

pdf(file_out)
par(mfrow=c(3,3))
taxonIDs <- unique(potDistr[["taxonID"]])

potDistr[["taxonID_directIndirect"]] <- paste(potDistr[["taxonID"]], potDistr[["directIndirect"]], sep = "_")
for(taxonID in taxonIDs)
{
	for(dI in c("direct", "indirect"))
	{
		idty_local <- idty[(idty[["taxonID"]] == taxonID) & (idty[["directIndirect"]] == dI),]
		
		identities <- idty_local[["Identity"]]		
		lengths <- idty_local[["Length"]]
		
		taxonName <- unique(potDistr[["Name"]][potDistr[["taxonID"]] == taxonID])[[1]]
		n_reads <- length(identities)	
		
		f <- n_reads / dim(idty)[[1]]
		
		plotIdentifier <- paste(taxonID, dI, sep = "_")
		if(f >= 0.001)
		{
			
	
			f_both_pot <- potDistr[["PotFrequency"]][potDistr[["taxonID"]] == taxonID][[1]]
			f_direct <- potDistr[["EMFrequency_direct"]][potDistr[["taxonID"]] == taxonID][[1]]
			f_indirect <- potDistr[["EMFrequency_indirect"]][potDistr[["taxonID"]] == taxonID][[1]]
		
			hist(lengths, main = paste("Read lengths", plotIdentifier), cex.main = 0.7)

			if(any((shifted[["directIndirect"]] == dI) & (shifted[["taxonID"]] == taxonID)))
			{
				barplot_shifted_y <- shifted[["P"]][(shifted[["directIndirect"]] == dI) & (shifted[["taxonID"]] == taxonID)]
				names(barplot_shifted_y) <- shifted[["identity"]][(shifted[["directIndirect"]] == dI) & (shifted[["taxonID"]] == taxonID)]
				barplot_shifted_y <- barplot_shifted_y[order(as.numeric(names(barplot_shifted_y)))]
				barplot(barplot_shifted_y, main = taxonName, cex.main = 0.7)
			}
			else
			{
				plot(0,0)
			}
			
			hist(identities, main = paste(plotIdentifier, " ", n_reads, " r.; ", f_both_pot, sep = ""), cex.main = 0.7)
			
			if(f_both_pot > f)
			{
				cat(paste(plotIdentifier, ": ", f, " (from idty) vs ", f_both_pot, " (from distr)\n", sep = ""))
			}
		}
	
	}
}

dev.off()


cat(paste("\nProduced file ", file_out, "\n\n", sep = ""))

