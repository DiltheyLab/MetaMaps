mappings <- "/scratch/tmp/MetaMap/hmp_set7.EM"
data_LandI <- read.delim(paste(mappings, ".lengthAndIdentitiesPerMappingUnit", sep = ""))
data_LandI <- data_LandI[data_LandI[["AnalysisLevel"]] == "EqualCoverageUnit",]

countsPerUnit <- table(data_LandI[["ID"]])
countsPerUnit <- sort(countsPerUnit, d = T)
freqPerUnit <- countsPerUnit/sum(countsPerUnit)

pdf("_plottedStatistics.pdf")

plotLengthIdentity_rows <- 4
plotLengthIdentity_cols <- 3
op <- par(mfrow=c(plotLengthIdentity_rows,plotLengthIdentity_cols)) 
plotLabels <- c()
identities_per_label <- list()
lengths_per_label <- list()
for(i in 1:length(countsPerUnit))
{
	iLabel <- names(countsPerUnit)[[i]]
	iCount <- countsPerUnit[[i]]
	iFreq <- freqPerUnit[[i]]
	if(iFreq >= 0.001)
	{
		iIdentities <- data_LandI[["Identity"]][data_LandI[["ID"]] == iLabel]
		iLengths <- data_LandI[["Length"]][data_LandI[["ID"]] == iLabel]
		identities_per_label[[iLabel]] <- iIdentities
		lengths_per_label[[iLabel]] <- iLengths
		plotLabels <- c(plotLabels, iLabel)
	}
}

for(label in plotLabels)
{
	hist(identities_per_label[[label]], main = label)
}
missingPlots <- length(plotLabels) %% (plotLengthIdentity_rows * plotLengthIdentity_cols)
for(i in 1:missingPlots)
{
	plot(0,0)
}	

for(label in plotLabels)
{
	hist(lengths_per_label[[label]], main = label)
}
for(i in 1:missingPlots)
{
	plot(0,0)
}	

dev.off()




