args = commandArgs(trailingOnly=TRUE)

prefix <- "/data/projects/phillippy/projects/MetaMap/tmp/hmp7_2_miniSeq+H"
minimumPlotFreq <- 0.001
if(length(args) > 0)
{
	prefix <- args[[1]]
}

if(length(args) > 1)
{
	minimumPlotFreq <- args[[2]]
}

data_LandI <- read.delim(paste(prefix, ".EM.lengthAndIdentitiesPerMappingUnit", sep = ""))
data_LandI <- data_LandI[data_LandI[["AnalysisLevel"]] == "EqualCoverageUnit",]

data_coverage <- read.delim(paste(prefix, ".EM.contigCoverage", sep = ""))

countsPerUnit <- table(data_LandI[["ID"]])
countsPerUnit <- sort(countsPerUnit, d = T)
freqPerUnit <- countsPerUnit/sum(countsPerUnit)

fn_output <- paste(prefix, '.identitiesAndCoverage.pdf', sep = "")
fn_output_coveragePlots <- paste(prefix, '.coveragePerContig_MetaMapComplete.pdf', sep = "")

plotLabels <- c()
identities_per_label <- list()
lengths_per_label <- list()
plotTaxonIDs <- list()
taxonID_2_mappingUnits <- list()
for(i in 1:length(countsPerUnit))
{
	iLabel <- names(countsPerUnit)[[i]]
	iCount <- countsPerUnit[[i]]
	iFreq <- freqPerUnit[[i]]
	if(iFreq >= minimumPlotFreq)
	{
		if(!(length(regmatches(iLabel, regexec("kraken:taxid\\|(x?\\d+)\\|", iLabel))[[1]]) == 2))
		{
			cat("Can't match ", iLabel, "\n")
		}
		stopifnot(length(regmatches(iLabel, regexec("kraken:taxid\\|(x?\\d+)\\|", iLabel))[[1]]) == 2)
		taxonID <- regmatches(iLabel, regexec("kraken:taxid\\|(x?\\d+)\\|", iLabel))[[1]][[2]]
		plotTaxonIDs[[taxonID]] <- 1
		
		if(!(taxonID %in% names(taxonID_2_mappingUnits)))
		{
			taxonID_2_mappingUnits[[taxonID]] <- list()
		}
		taxonID_2_mappingUnits[[taxonID]][[iLabel]] <- 1
		
		#iIdentities <- data_LandI[["Identity"]][data_LandI[["ID"]] == iLabel]
		#iLengths <- data_LandI[["Length"]][data_LandI[["ID"]] == iLabel]
		#identities_per_label[[iLabel]] <- iIdentities
		#lengths_per_label[[iLabel]] <- iLengths
		#plotLabels <- c(plotLabels, iLabel)
	}
}

# pdf(fn_output_coveragePlots, width = 15, height = 5)
# for(taxonID in names(plotTaxonIDs))
# {

# }
# dev.off()



pdf(fn_output, width = 12, height = 8)
#par(mfrow=c(1,3),oma=c(0,0,2,0))
par(mar = c(5, 3, 5, 3), oma = c(0,0,2,0))
m <- rbind(c(1,2,3), c(4, 4, 4))
layout(m)

allIdentityDensities <- c()
allIdentity_min_not0 <- c()
for(doWhat in c("limits", "plot"))
{
	for(taxonID in names(plotTaxonIDs))
	{
		indices_coverage_taxonID <- which( data_coverage[["taxonID"]] == taxonID )
		stopifnot(length(indices_coverage_taxonID) > 0)
		
		taxonLabel <- as.character(data_coverage[["equalCoverageUnitLabel"]][[indices_coverage_taxonID[[1]]]])
				reads_count <- 0
		allReads_lengths <- c()
		allReads_identities <- c()
		allWindows_coverages <- c()
		for(mappingUnit in names(taxonID_2_mappingUnits[[taxonID]]))
		{	
			reads_count <- reads_count + countsPerUnit[[mappingUnit]]
			coverage_indices_mappingUnit <- which(data_coverage[["contigID"]] == mappingUnit)
			stopifnot(length(coverage_indices_mappingUnit) > 0)
			allWindows_coverages <- c(allWindows_coverages, data_coverage[["readCoverage"]][coverage_indices_mappingUnit])
			
			LandI_indices_mappingUnit <- which(data_LandI[["ID"]] == mappingUnit)
			stopifnot(length(LandI_indices_mappingUnit) > 0)		
			
			allReads_lengths <- c(allReads_lengths, data_LandI[["Length"]][LandI_indices_mappingUnit])
			allReads_identities <- c(allReads_identities, data_LandI[["Identity"]][LandI_indices_mappingUnit])
			
		}
		stopifnot(length(allReads_lengths) == reads_count)
		
		if(doWhat == "plot")
		{
			histogram_length <- hist(allReads_lengths, plot = F)
			plot(histogram_length, main = "Read length histogram", xlim = c(0, max(data_LandI[["Length"]])), xlab = "Read length")
		}
		
		# min_identity <- 
		vector_identities <- rep(0, 101)
		names(vector_identities) <- 0:100	
		allReads_identities_100 <- round(allReads_identities*100)+1
		identitiesTable <- table(allReads_identities_100)
		identitiesTable <- identitiesTable/sum(identitiesTable)
		for(idt in names(identitiesTable))
		{
			vector_identities[[as.integer(idt)]] <- identitiesTable[[idt]]		
		}	
		
		if(doWhat == "plot")
		{
			barplot(vector_identities[min(allIdentity_min_not0):length(vector_identities)], xlab = "Identity", ylab = "Density", main = paste("Read identities"), ylim = c(0, max(allIdentityDensities)))		
			histogram_coverage <- hist(allWindows_coverages, plot = F)	
			plot(histogram_coverage, main = "Genome window coverage histogram", xlab = "Coverage")
			title(paste("MetaMaps mapping summary for ", taxonLabel, " (taxon ID ", taxonID, ") - ", reads_count, " mapped reads assigned", sep = ""), outer=TRUE, cex.main = 1.5)			
		}
		else
		{
			allIdentityDensities <- c(allIdentityDensities, vector_identities)
			allIdentity_min_not0 <- c(allIdentity_min_not0, min(which(vector_identities != 0)))
		}
		
		if(doWhat == "plot")
		{
			indices_coverage_taxonID <- which( data_coverage[["taxonID"]] == taxonID )
			stopifnot(length(indices_coverage_taxonID) > 0)
			
			taxonLabel <- as.character(data_coverage[["equalCoverageUnitLabel"]][[indices_coverage_taxonID[[1]]]])
			
			reads_count <- 0	 
			allWindows_coverages <- c()
			allWindows_coverages_colors <- c()
			mappingUnits <- names(taxonID_2_mappingUnits[[taxonID]])	
			mappingUnits <- unique(data_coverage[["contigID"]][data_coverage[["taxonID"]] == taxonID])
			for(mappingUnitI in 1:length(mappingUnits))
			{	
				reads_count <- reads_count + countsPerUnit[[mappingUnit]]	
				mappingUnit <- mappingUnits[[mappingUnitI]]
				# reads_count <- reads_count + countsPerUnit[[mappingUnit]]
				coverage_indices_mappingUnit <- which(data_coverage[["contigID"]] == mappingUnit)
				stopifnot(length(coverage_indices_mappingUnit) > 0)
				coverages_thisMappingUnit <- data_coverage[["readCoverage"]][coverage_indices_mappingUnit]
				allWindows_coverages <- c(allWindows_coverages, coverages_thisMappingUnit)
				thisMappingUnit_color <- "blue"
				if((mappingUnitI %% 2) == 0)
				{
					thisMappingUnit_color <- "red"
				}
				allWindows_coverages_colors <- c(allWindows_coverages_colors, rep(thisMappingUnit_color, length(coverages_thisMappingUnit)))
			}
			plot(1:length(allWindows_coverages), allWindows_coverages, col = allWindows_coverages_colors, main = paste("Genome-wide coverage over all contigs for ", taxonLabel, " (taxon ID ", taxonID, ") - ", reads_count, " mapped reads assigned",  sep = ""), xlab = "Coordinate concatenated genome (1000s)", ylab = "Coverage", pch = 20, cex.main = 1)		
		}
	}
}
dev.off()


cat("\nGenerated file: ", fn_output, "\n")

