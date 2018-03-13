prefix <- "/data/projects/phillippy/projects/MetaMap/tmp/hmp7_2_miniSeq+H"

fn_observed_EM <- paste(prefix, "", sep = "")


fn_shifted_U <- paste(prefix, ".U.shiftedHistogramsPerTaxonID", sep = "")
fn_observed_U <- paste(prefix, ".U.lengthAndIdentitiesPerTaxonID", sep = "")

D_shifted_U <- read.delim(fn_shifted_U, stringsAsFactors = F, header = T)
D_observed_U <- read.delim(fn_observed_U, stringsAsFactors = F)
D_observed_U[["identityIndex"]] <- round(D_observed_U[["Identity"]]*100)+1


fn_output <- paste(prefix, '.identities_MetaMapUnknown.pdf', sep = "")
pdf(fn_output, width = 20, height = 5)

one_taxonID_direct <- D_shifted_U[["taxonID"]][D_shifted_U[["directIndirect"]] == "direct"]
one_taxonID_direct <- one_taxonID_direct[[1]]

for(plotTaxonID in c(1654))
{
	par(mfrow=c(1,4),oma=c(0,0,2,0))
	
	idx_baseline_U <- which(D_shifted_U[["taxonID"]] == one_taxonID_direct)
	idx_shifted_U <- which(D_shifted_U[["taxonID"]] == plotTaxonID)
	idx_observed_U <- which(D_observed_U[["taxonID"]] == plotTaxonID)
	
	vector_baseline <- rep(0, 101)
	names(vector_baseline) <- 0:100
	for(idxI in idx_baseline_U)
	{
		idt <- as.integer(round(D_shifted_U[["identity"]][[idxI]]))+1
		stopifnot(idt <= 100)
		vector_baseline[[idt]] <- D_shifted_U[["P"]][[idxI]]
	}
	barplot(vector_baseline, xlab = "Identity", ylab = "Density", main = "H_D - Read identities")
	
	taxonName <- unique(D_observed_U[["taxonName"]][idx_observed_U])	
	stopifnot(length(taxonName) == 1)
	taxonName <- taxonName[[1]]
	
	vector_expected <- rep(0, 101)
	names(vector_expected) <- 0:100
	for(idxI in idx_shifted_U)
	{
		idt <- as.integer(round(D_shifted_U[["identity"]][[idxI]]))+1
		stopifnot(idt <= 100)
		vector_expected[[idt]] <- D_shifted_U[["P"]][[idxI]]
	}
	barplot(vector_expected, xlab = "Identity", ylab = "Density", main = paste("[H_D x H_(t;l)] - Expected identities for ", plotTaxonID, " (0 = expected unmapped)", sep = ""))
	
	vector_expected_no0 <- vector_expected
	vector_expected_no0[[1]] <- 0
	vector_expected_no0 <- vector_expected_no0 / sum(vector_expected_no0)
	barplot(vector_expected_no0, xlab = "Identity", ylab = "Density", main = paste("Trunacted [H_D x H_(t;l)] for ", plotTaxonID, sep = ""))
	
	vector_observed <- rep(0, 101)
	names(vector_observed) <- 0:100	
	observedTable <- table(D_observed_U[["identityIndex"]][idx_observed_U])
	observedTable <- observedTable/sum(observedTable)
	for(idt in names(observedTable))
	{
		vector_observed[[as.integer(idt)]] <- observedTable[[idt]]		
	}	
	barplot(vector_observed, xlab = "Identity", ylab = "Density", main = paste("Observed identities of reads assigned to ", plotTaxonID, sep = ""))
	
	#plot(0:100, vector_expected, type = "s")
		
	#vector_observed <- rep(0, 100)
	#names(vector_observed) <- 0:100

	#plot(D_shifted_U[["identity"]][idx_shifted_U], D_shifted_U[["P"]][idx_shifted_U], main = "Expected identities")
	
	title(paste("MetaMap-Unknown for ", taxonName, " (taxon ID ", plotTaxonID, ")", sep = ""), outer=TRUE)
}
dev.off()
