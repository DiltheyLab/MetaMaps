/*
 * fU.h
 *
 *  Created on: Jul 10, 2017
 *      Author: diltheyat
 */

#ifndef META_FU_H_
#define META_FU_H_

#include <string>
#include <assert.h>
#include <vector>
#include <map>
#include <fstream>
#include <utility>

#include "fEM.h"
#include "util.h"

#include "fU_helper.h"


namespace meta
{

void cleanF_U(std::pair<std::map<std::string, double>, std::map<std::string, double>>& f, const std::pair<std::map<std::string, size_t>, std::map<std::string, size_t>>& assignedReads, size_t distributedReads);
std::string getReadIDFromReadLines(const std::vector<std::string>& readLines);

void produceShiftedHistograms(std::string output_shifted_frequencies_file, identityManager& iM, const std::pair<std::map<std::string, double>, std::map<std::string, double>>& f)
{	
	std::ofstream strout_histogram(output_shifted_frequencies_file);
	assert(strout_histogram.is_open());

	strout_histogram 	<< "taxonID" << "\t"
						<< "directIndirect" << "\t"
						<< "identity" << "\t"
						<< "P" << "\n";				
	for(auto tF : f.first)
	{
		if(tF.second > 1e-5)
		{
			std::map<int, double> nodeHistogram = iM.getHistogramForNode(tF.first, true);
			double hE_sum = 0;
			for(auto hE : nodeHistogram)
			{
				strout_histogram 	<< tF.first << "\t"
									<< "direct" << "\t"
									<< hE.first << "\t"
									<< hE.second << "\n";
				hE_sum += hE.second;
			}
			assert(abs(1 - hE_sum) <= 1e-3);
		}
	}
	
	for(auto tF : f.second)
	{
		std::map<int, double> nodeHistogram = iM.getHistogramForNode(tF.first, false);
		double hE_sum = 0;
		for(auto hE : nodeHistogram)
		{
			strout_histogram 	<< tF.first << "\t"
								<< "indirect" << "\t"
								<< hE.first << "\t"
								<< hE.second << "\n";
			hE_sum += hE.second;
		}
		assert(abs(1 - hE_sum) <= 1e-3);
	}
		
	

}

void printGenusLevelSummary(const taxonomy& T, std::pair<std::map<std::string, double>, std::map<std::string, double>> frequencies)
{
	std::map<std::string, double> combinedF;
	
	for(auto fE : frequencies.first)
	{
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(fE.first, {"genus"});
		std::string fE_genus = upwardByLevel.at("genus");
		std::string fE_genus_name = (fE_genus != "Undefined") ? T.getNode(fE_genus).name.scientific_name : fE_genus;
		
		if(combinedF.count(fE_genus_name) == 0)
		{
			combinedF[fE_genus_name] = 0;
		}
		combinedF.at(fE_genus_name) += fE.second;
	}
	
	for(auto fE : frequencies.second)
	{
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(fE.first, {"genus"});
		std::string fE_genus = upwardByLevel.at("genus");
		std::string fE_genus_name = (fE_genus != "Undefined") ? T.getNode(fE_genus).name.scientific_name : fE_genus;
		
		if(combinedF.count(fE_genus_name) == 0)
		{
			combinedF[fE_genus_name] = 0;
		}
		combinedF.at(fE_genus_name) += fE.second;
	}
	
	std::vector<std::string> f_keys;
	for(auto fE : combinedF)
	{
			f_keys.push_back(fE.first);
	}
	std::sort(f_keys.begin(), f_keys.end());
	
	std::cout << "Genus-level summary:\n";
	for(auto fK : f_keys)
	{		
		double f = combinedF.at(fK);
		if(f >= 0.001)
		{
			std::cout << fK << "\t" << f << "\n";
		}
	}
	std::cout << "\n\n";
}

void producePotFile_U(std::string outputFN, const taxonomy& T, std::tuple<std::map<std::string, double>, std::map<std::string, double>, std::map<std::string, double>> frequencies, std::tuple<std::map<std::string, size_t>, std::map<std::string, size_t>> readCount, size_t mappableReads, const std::set<std::string>& relevantTaxonIDs_mappable)
{
	//double initial_f_sum = 0;
	std::set<std::string> combinedKeys;
	for(auto f : std::get<0>(frequencies))
	{
		combinedKeys.insert(f.first);
		//initial_f_sum += f.second;
	}
	for(auto f : std::get<1>(frequencies))
	{
		combinedKeys.insert(f.first);
		//initial_f_sum += f.second;		
	}
	for(auto f : std::get<2>(frequencies))
	{
		combinedKeys.insert(f.first);
		//initial_f_sum += f.second;
	}

	//assert(abs(1 - initial_f_sum) <= 1e-3);
	
	for(auto f : std::get<0>(readCount))
	{
		combinedKeys.insert(f.first);
	}
	for(auto f : std::get<1>(readCount))
	{
		combinedKeys.insert(f.first);
	}

	std::map<std::string, std::tuple<std::map<std::string, double>, std::map<std::string, double>, std::map<std::string, double>>> frequencies_perLevel;
	std::map<std::string, std::pair<std::map<std::string, size_t>, std::map<std::string, size_t>>> readCount_perLevel;

	std::map<std::string, std::set<std::string>> combinedKeys_perLevel;

	std::set<std::string> targetLevels = getRelevantLevelNames();
	
	for(auto taxonID : combinedKeys)
	{
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(taxonID, targetLevels);
		upwardByLevel["definedAndHypotheticalGenomes"] = taxonID;
		upwardByLevel["definedGenomes"] = taxonID;

		for(auto uN : upwardByLevel)
		{
			std::string level = uN.first;
			std::string levelValue = uN.second;
			if(level == "definedGenomes")
			{
				if(! relevantTaxonIDs_mappable.count(levelValue))
					continue;
			}
			combinedKeys_perLevel[level].insert(levelValue);

			if(std::get<0>(frequencies_perLevel[level]).count(levelValue) == 0)
			{
				std::get<0>(frequencies_perLevel[level])[levelValue] = 0;
				std::get<1>(frequencies_perLevel[level])[levelValue] = 0;
				std::get<2>(frequencies_perLevel[level])[levelValue] = 0;
				readCount_perLevel[level].first[levelValue] = 0;
				readCount_perLevel[level].second[levelValue] = 0;
			}

			if(std::get<0>(frequencies).count(taxonID))
			{
				std::get<0>(frequencies_perLevel.at(level)).at(levelValue) += std::get<0>(frequencies).at(taxonID);
			}

			if(std::get<1>(frequencies).count(taxonID))
			{
				std::get<1>(frequencies_perLevel.at(level)).at(levelValue) += std::get<1>(frequencies).at(taxonID);
			}

			if(std::get<2>(frequencies).count(taxonID))
			{
				std::get<2>(frequencies_perLevel.at(level)).at(levelValue) += std::get<2>(frequencies).at(taxonID);
			}


			if(std::get<0>(readCount).count(taxonID))
			{
				readCount_perLevel.at(level).first.at(levelValue) += std::get<0>(readCount).at(taxonID);
			}

			if(std::get<1>(readCount).count(taxonID))
			{
				readCount_perLevel.at(level).second.at(levelValue) += std::get<1>(readCount).at(taxonID);
			}
		}
	}


	std::ofstream strout_frequencies(outputFN);
	assert(strout_frequencies.is_open());

	strout_frequencies << "AnalysisLevel" << "\t" <<  "taxonID" << "\t" << "Name" << "\t" << "readsDirectlyAssigned_inDB" << "\t" << "readsDirectlyAssigned_potentiallyNovel" << "\t" << "frDirect" << "\t" << "frIndirect" << "\t" <<  "frFromUnmapped" << "\t" << "Absolute" << "\t" << "PotFrequency" << "\n";

	for(auto l : combinedKeys_perLevel)
	{
		std::string levelName = l.first;
		double levelFreqSum = 0;
		size_t levelReadSum = 0;
		for(auto taxonID : l.second)
		{
			if(taxonID == "Undefined")
				continue;

			std::string taxonIDName = T.getNode(taxonID).name.scientific_name;
			
			size_t taxonID_combinedReads = (readCount_perLevel.at(levelName).first.at(taxonID) + readCount_perLevel.at(levelName).second.at(taxonID));
			double taxonID_combinedFreq = (std::get<0>(frequencies_perLevel.at(levelName)).at(taxonID) + std::get<1>(frequencies_perLevel.at(levelName)).at(taxonID) + std::get<2>(frequencies_perLevel.at(levelName)).at(taxonID));
			
			strout_frequencies <<
					levelName << "\t" <<
					taxonID << "\t" <<
					taxonIDName << "\t" <<
					readCount_perLevel.at(levelName).first.at(taxonID) << "\t" <<
					readCount_perLevel.at(levelName).second.at(taxonID) << "\t" <<
					std::get<0>(frequencies_perLevel.at(levelName)).at(taxonID) << "\t" <<
					std::get<1>(frequencies_perLevel.at(levelName)).at(taxonID) << "\t" <<
					std::get<2>(frequencies_perLevel.at(levelName)).at(taxonID) << "\t" <<
					taxonID_combinedReads << "\t" <<
					taxonID_combinedFreq << "\n";

			levelReadSum += taxonID_combinedReads;					
			levelFreqSum += taxonID_combinedFreq;
		}
				
		long long levelReadsUnclassified = mappableReads - levelReadSum;
		assert(levelReadsUnclassified >= 0);
		
		
		assert((levelFreqSum >= 0) && (levelFreqSum <= (1+1e-3)));
		if(levelFreqSum > 1)
			levelFreqSum = 1;
		
		double levelFreqUnclassified = 1 - levelFreqSum;
		//levelFreqUnclassified = 0;
		
		strout_frequencies <<
				levelName << "\t" <<
				0 << "\t" <<
				"Unclassified" << "\t" <<
				0 << "\t" <<
				0 << "\t" <<
				0 << "\t" <<
				0 << "\t" <<
				0 << "\t" <<
				levelReadsUnclassified << "\t" <<
				levelFreqUnclassified << "\n";		
	}

	strout_frequencies.close();
}

class oneMappingLocation_U
{
public:
	std::string taxonID;
	std::string contigID;
	double identity;
	size_t readLength;
	size_t start;
	size_t stop;
	double p;
	double l;
	bool direct;
};

oneMappingLocation_U getBestMapping_U(const std::vector<oneMappingLocation_U>& locations)
{
	assert(locations.size() > 0);
	double maxP;
	size_t whichI_maxP;
	for(size_t i = 0; i < locations.size(); i++)
	{
		const auto& l = locations.at(i);
		if((i == 0) || (l.p > maxP))
		{
			whichI_maxP = i;
			maxP = l.p;
		}
	}
	return locations.at(whichI_maxP);
}


std::vector<oneMappingLocation_U> getMappingLocations_U(identityManager& iM, const std::map<std::string, std::vector<std::string>>& indirectUpwardNodes, const std::pair<std::map<std::string, double>, std::map<std::string, double>>& f, const std::vector<std::string>& readLines)
{
	assert(readLines.size() > 0);

	bool verbose = false;
	std::set<std::string> directMappingsTaxonIDs;
	for(auto line : readLines)
	{
		std::vector<std::string> line_fields = split(line, " ");
		std::string contigID = line_fields.at(5);
		std::string contig_taxonID = extractTaxonId(contigID);
		directMappingsTaxonIDs.insert(contig_taxonID);
	}
	
	if(directMappingsTaxonIDs.count("1063"))
	{
		// verbose = true;
	}
		
	std::string readID;
	long long readLength = -1;
		
	for(auto line : readLines)
	{
		std::vector<std::string> line_fields = split(line, " ");
		std::string contigID = line_fields.at(5);
		std::string contig_taxonID = extractTaxonId(contigID);

		if(readID.length() == 0)
		{
			readID = line_fields.at(0);
		}
		else
		{
			assert(readID == line_fields.at(0));
		}
		
		
		if(readLength == -1)
			readLength = std::stoi(line_fields.at(1));
	}
	
	if(verbose)
	{
		std::cout << "\ngetMappingLocations_U(..): Analyze read " << readID << ", length " << readLength << "\n";
	}
			
	std::map<std::string, oneMappingLocation_U> bestDirectMappings;
	std::map<std::string, oneMappingLocation_U> bestInDirectMappings;
	
	for(auto line : readLines)
	{
		std::vector<std::string> line_fields = split(line, " ");

		std::string contigID = line_fields.at(5);
		size_t contigID_start = std::stoull(line_fields.at(7));
		size_t contigID_stop = std::stoull(line_fields.at(8));

		std::string contig_taxonID = extractTaxonId(contigID);

		if(readLength == -1)
			readLength = std::stoi(line_fields.at(1));

		// todo we might want to use field 12 (mashmap-internal estimate) again
		double identity = std::stod(line_fields.at(9))/100.0;
		assert(identity >= 0);
		assert(identity <= 1);
		int identityInt = int((identity * 100) + 0.5);

		if(verbose)
			std::cout << "\tLine with identity " << identityInt << "\n";
				
		oneMappingLocation_U l;
		l.taxonID = contig_taxonID;
		l.contigID = contigID;
		l.start = contigID_start;
		l.stop = contigID_stop;
		l.p = 0;
		l.readLength = readLength;
		l.identity = identity;
		l.direct = true;
		l.l = f.first.at(contig_taxonID) * iM.getIdentityP(identityInt, contig_taxonID, readLength, true);

		if((bestDirectMappings.count(contig_taxonID) == 0) || (bestDirectMappings.at(contig_taxonID).identity < l.identity))
		{
			bestDirectMappings[contig_taxonID] = l;
		}
		
		if(verbose)
			std::cout << "\t\tDirect to " << contig_taxonID << ": " << f.first.at(contig_taxonID) << " x " << iM.getIdentityP(identityInt, contig_taxonID, readLength, true) << "\n";
		
		//std::cout << "\nTaxon " << contig_taxonID << " indirect-upward: " << indirectUpwardNodes.at(contig_taxonID).size() << "\n" << std::flush;
		
		for(auto indirectTaxon : indirectUpwardNodes.at(contig_taxonID))
		{
			oneMappingLocation_U lI;
			lI.taxonID = indirectTaxon;
			lI.contigID = "";
			lI.start = 0;
			lI.stop = 0;
			lI.p = 0;
			lI.readLength = readLength;
			lI.identity = identity;
			lI.direct = false;
			lI.l = f.second.at(indirectTaxon) * iM.getIdentityP(identityInt, indirectTaxon, readLength, false);
			
			if((bestInDirectMappings.count(indirectTaxon) == 0) || (bestInDirectMappings.at(indirectTaxon).identity < lI.identity))
			{
				bestInDirectMappings[indirectTaxon] = lI;
			}
						
			//if(verbose || (indirectTaxon == "1063"))
			if(1 == 0)
			{
				std::cout << "\t\tDirect to " << contig_taxonID << ": " << f.first.at(contig_taxonID) << " x " << iM.getIdentityP(identityInt, contig_taxonID, readLength, true) << "\n";				
				std::cout << "\t\tIndirect to " << indirectTaxon << ": " << f.second.at(indirectTaxon) << " x " << iM.getIdentityP(identityInt, indirectTaxon, readLength, false) << "\n";
			}		
		}
	}
	
	std::vector<oneMappingLocation_U> mappingLocations;			
	for(auto l : bestDirectMappings)
	{
			mappingLocations.push_back(l.second);
	}
	for(auto l : bestInDirectMappings)
	{
			mappingLocations.push_back(l.second);
	}

	if(verbose)
	{
		std::cout << "\tPost-normalization I:\n";
		for(auto l : mappingLocations)
		{
			std::cout << "\t\t" << l.taxonID << (l.direct ? " direct " : " indirect ") << " " << l.p << "\n";
		}
		std::cout << "\n" << std::flush;
	}
	
	double l_sum = 0;
	for(const auto& l : mappingLocations)
	{
		l_sum += l.l;
	}
	
	if(verbose)
	{
		std::cout << "\tPost-normalization II:\n";
		for(auto l : mappingLocations)
		{
			std::cout << "\t\t" << l.taxonID << (l.direct ? " direct " : " indirect ") << " " << l.p << "\n";
		}
		std::cout << "\n" << std::flush;
	}
	
	
	assert(l_sum > 0);
	for(auto& l : mappingLocations)
	{
		l.p = l.l/l_sum;
	}
	
	if(verbose)
	{
		std::cout << "\tPost-normalization III:\n";
		for(auto l : mappingLocations)
		{
			std::cout << "\t\t" << l.taxonID << (l.direct ? " direct " : " indirect ") << " " << l.p << "\n";
		}
		std::cout << "\n" << std::flush;
	}
		
	if(verbose)
	{
		std::cout << "\tPost-normalization:\n";
		for(auto l : mappingLocations)
		{
			std::cout << "\t\t" << l.taxonID << (l.direct ? " direct " : " indirect ") << " " << l.p << "\n";
		}
		std::cout << "\n" << std::flush;
		
		assert(3 == 4);
	}
	
	std::string printTaxonID = "1063";
	printTaxonID = "";
	if(bestDirectMappings.count(printTaxonID) || bestInDirectMappings.count(printTaxonID))
	{
		std::cout << "Have mappings to " << printTaxonID << "\n";
		for(auto mL : mappingLocations)
		{
			int identityInt = int((mL.identity * 100) + 0.5);
				
			std::cout << "\tOne mapping " << (mL.direct ? "DIRECT" : "INDIRECT") << "\n";
			std::cout << "\t\tTaxon ID: " << mL.taxonID << "\n";
			std::cout << "\t\tIdentity: " << mL.identity << "\n";
			std::cout << "\t\tl       : " << iM.getIdentityP(identityInt, mL.taxonID, mL.readLength, mL.direct) << "\n";
			std::cout << "\t\tf       : " << (mL.direct ? f.first.at(mL.taxonID) : f.second.at(mL.taxonID)) << "\n";
			std::cout << "\t\tP       : " << mL.p << "\n";
		}
	}

	return mappingLocations;
}

std::pair<int, int> getMinMaxIdentities(std::string mappedFile)
{
	int minIdentity = -1;
	int maxIdentity = -1;
	
	// size_t processedRead = 0;
	std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines) -> void
	{
		for(const auto& line : readLines)
		{
			std::vector<std::string> line_fields = split(line, " ");
			// todo consider 12
			double identity = std::stod(line_fields.at(9))/100.0;
			assert(identity >= 0);
			assert(identity <= 1);
			int identityInt = int((identity * 100) + 0.5);
			assert(identityInt >= 0);
			assert(identityInt <= 100);
			
			if((minIdentity == -1) || (identityInt < minIdentity))
			{
				minIdentity = identityInt;
			}
			if((maxIdentity == -1) || (identityInt > maxIdentity))
			{
				maxIdentity = identityInt;
			}	
		}		
	};
	callBackForAllReads(mappedFile, processOneRead);	
	
	assert(maxIdentity > 1);
	
	return std::make_pair(minIdentity, maxIdentity);
}

void doU(std::string DBdir, std::string mappedFile, size_t minimumReadsPerBestContig)
{
	unsigned int round_first_unknown = 5;

	taxonomy T(DBdir+"/taxonomy");

	std::set<std::string> taxonIDsInMappings = getTaxonIDsFromMappingsFile(mappedFile);

	std::set<std::string> mappableTaxonIDs = getDirectlyMappableTaxonIDs(DBdir);
	
	std::string fn_fittedLengthAndIdentities = mappedFile + ".EM.lengthAndIdentitiesPerMappingUnit";
	if(! fileExists(fn_fittedLengthAndIdentities))
	{
		std::cerr << "\n\nERROR: File " << fn_fittedLengthAndIdentities << " not existing.\n\nThis file is generated automatically by the EM step. Run the EM step first.\n\n";
	}
	
	std::pair<int, int> identity_minmax_inAlignments = getMinMaxIdentities(mappedFile);

	identityAndReadLengthHistogram iAndL;
	iAndL.readFromEMOutput(fn_fittedLengthAndIdentities, identity_minmax_inAlignments, minimumReadsPerBestContig);

	std::string fn_tree_selfSimilarities = DBdir + "/selfSimilarities.txt";
	treeAdjustedIdentities tAI;
	tAI.readFromFile(fn_tree_selfSimilarities, taxonIDsInMappings, T);

	identityManager iM(iAndL, tAI);

	// mapping stats
	std::map<std::string, size_t> mappingStats = getMappingStats(mappedFile);
	size_t nTotalReads = mappingStats.at("TotalReads");	
	size_t nTooShort = mappingStats.at("ReadsTooShort");
	size_t nUnmapped = mappingStats.at("ReadsNotMapped");
	size_t nMapped = mappingStats.at("ReadsMapped");
	assert(nTotalReads == (nTooShort + nUnmapped + nMapped));
	
	size_t nReadsMappable = nTotalReads - nTooShort;
	assert(nReadsMappable <= nTotalReads);
	
	std::vector<size_t> unmappedReadsLengths = getUnmappedReadsStats(mappedFile);
	assert(unmappedReadsLengths.size() == nUnmapped); // can remove later
	
	// Todo consider removing the following
	std::cout << "\nnTotalReads: " << nTotalReads << "\n" << std::flush;	
	std::cout << "Reads too short: " << nTooShort << "\n" << std::flush;	
	std::cout << "Reads mappable: " << nReadsMappable << "\n" << std::flush;	
	std::cout << "Unmapped (but long enough) reads: " << nUnmapped << "\n\n" << std::flush;
	
	// set up nodes to map to

	std::set<std::string> relevantTaxonIDs_direct = taxonIDsInMappings;
	std::set<std::string> relevantTaxonIDs_indirect;

	std::map<std::string, std::vector<std::string>> indirectUpwardNodes;
	for(auto tI : taxonIDsInMappings)
	{
		std::vector<std::string> upwardTaxonIDs = T.getUpwardNodes(tI);
		indirectUpwardNodes[tI] = std::vector<std::string>();
		for(auto uTI : upwardTaxonIDs)
		{
			if(tAI.nodeForIndirectAttachment(uTI))
			{
				relevantTaxonIDs_indirect.insert(uTI);
				indirectUpwardNodes[tI].push_back(uTI);
				
				/*
				if((uTI == "67753") || (uTI == "1914297"))
				{
					std::cout << "uTI:" << uTI << "\n";
					std::cout << "tI:" << tI << "\n";
					std::cout << std::flush;
				}
				*/
			}
		}
	}

	// set up initial frequency distribution

	std::pair<std::map<std::string, double>, std::map<std::string, double>> f;
	size_t combined_n_taxonIDs = relevantTaxonIDs_direct.size() + relevantTaxonIDs_indirect.size();
	for(auto nI : relevantTaxonIDs_direct)
	{
		f.first[nI] = 1.0/(double)combined_n_taxonIDs;
	}
	for(auto nI : relevantTaxonIDs_indirect)
	{
		f.second[nI] = 1.0/(double)combined_n_taxonIDs;
	}

	// set up U-EM

	auto normalize_f = [](std::pair<std::map<std::string, double>, std::map<std::string, double>>& fToN)
	{
		double f_sum = 0;
		for(auto tF : fToN.first)
		{
			f_sum += tF.second;
		}
		for(auto tF : fToN.second)
		{
			f_sum += tF.second;
		}

		// todo remove
		std::cout << "f_sum: " << f_sum << "\n" << std::flush;
		
		for(auto& tF : fToN.first)
		{
			tF.second /= f_sum;
		}
		for(auto& tF : fToN.second)
		{
			tF.second /= f_sum;
		}
		
		double f_sum_after = 0;
		for(auto tF : fToN.first)
		{
			f_sum_after += tF.second;
		}
		for(auto tF : fToN.second)
		{
			f_sum_after += tF.second;
		}	
		assert(abs(1 - f_sum_after) <= 1e-3);
	};
	
	auto normalize_f_triplet = [](std::tuple<std::map<std::string, double>, std::map<std::string, double>, std::map<std::string, double>>& fToN)
	{
		double f_sum = 0;
		for(auto tF : std::get<0>(fToN))
		{
			f_sum += tF.second;
		}
		for(auto tF : std::get<1>(fToN))
		{
			f_sum += tF.second;
		}
		for(auto tF : std::get<2>(fToN))
		{
			f_sum += tF.second;
		}

		// todo remove
		std::cout << "f_sum: " << f_sum << "\n" << std::flush;

		for(auto& tF : std::get<0>(fToN))
		{
			tF.second /= f_sum;
		}
		for(auto& tF : std::get<1>(fToN))
		{
			tF.second /= f_sum;
		}
		for(auto& tF : std::get<2>(fToN))
		{
			tF.second /= f_sum;
		}

		double f_sum_after = 0;
		for(auto tF : std::get<0>(fToN))
		{
			f_sum_after += tF.second;
		}
		for(auto tF : std::get<1>(fToN))
		{
			f_sum_after += tF.second;
		}
		for(auto tF : std::get<2>(fToN))
		{
			f_sum_after += tF.second;
		}
		assert(abs(1 - f_sum_after) <= 1e-3);
	};

	std::cout << "Starting EM-U..." << std::endl;
	double ll_lastIteration;
	size_t EMiteration = 0;
	bool continueEM = true;
	while(continueEM)
	{
		double ll_thisIteration_mapped = 0;
		std::pair<std::map<std::string, double>, std::map<std::string, double>> f_nextIteration = f;
		for(auto& fNextEntry : f_nextIteration.first)
		{
			fNextEntry.second = 0;
		}
		for(auto& fNextEntry : f_nextIteration.second)
		{
			fNextEntry.second = 0;
		}
		
		std::map<std::string, size_t> bestMappings_perTaxonID;

		size_t allReads_has_1654 = 0;
		size_t allReads_has_209 = 0;
		size_t allReads_has_1654_and_1301 = 0;

		size_t processedRead = 0;
		int printedReads = 0;
		std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines) -> void
		{
			processedRead++;
			if((processedRead % 10000) == 0)
			{
				std::cout << "\r EM-U round " << EMiteration << ", read << " << processedRead << " / " << mappingStats.at("ReadsMapped") << "   " << std::flush;
			}			
			// todo reinstate
			// std::cout << "\r EM round " << EMiteration << ", mapped read << " << processedRead << " / " << mappingStats.at("ReadsMapped") << "   " << std::flush;
					
			assert(readLines.size() > 0);
			std::vector<oneMappingLocation_U> mappingLocations = getMappingLocations_U(iM, indirectUpwardNodes, f, readLines);

			oneMappingLocation_U bestMapping = getBestMapping_U(mappingLocations);
			if(bestMappings_perTaxonID.count(bestMapping.taxonID) == 0)
			{
				bestMappings_perTaxonID[bestMapping.taxonID] = 0;
			}
			bestMappings_perTaxonID.at(bestMapping.taxonID)++;
			
			std::set<std::string> observedTaxonIDs_direct;
			std::set<std::string> observedTaxonIDs_indirect;
			
			bool has_1654 = false;
			bool has_209 = false;
			double l_read = 0;
			double p_read = 0;
			for(const auto& mL : mappingLocations)
			{
				l_read += mL.l;
				p_read += mL.p;
				if(mL.direct)
				{
					assert(f_nextIteration.first.count(mL.taxonID));					
					f_nextIteration.first.at(mL.taxonID) += mL.p;
					assert(observedTaxonIDs_direct.count(mL.taxonID) == 0);
					observedTaxonIDs_direct.insert(mL.taxonID);
				}
				else
				{
					assert(f_nextIteration.second.count(mL.taxonID));										
					f_nextIteration.second.at(mL.taxonID) += mL.p;
					assert(observedTaxonIDs_indirect.count(mL.taxonID) == 0);	
					observedTaxonIDs_indirect.insert(mL.taxonID);	

					if(mL.taxonID == "1654")
						has_1654 = true;
					
					if(mL.taxonID == "209")
						has_209 = true;
				}
			}
			
			if(has_1654 && 0)
			{
				// assert(EMiteration == 0);
				
				allReads_has_1654++;
				if(has_209)
				{
					//allReads_has_1654_and_1301++;
				}
			}
			
			
			if(has_209 && 0)
			{
				allReads_has_209++;
				
				/// todo
				if(printedReads < 15)
				{
					std::cout << "Round " << EMiteration << " read ID " << getReadIDFromReadLines(readLines) << "\n";
					for(const auto& mL : mappingLocations)
					{
						std::cout << "\t" << mL.taxonID << " " << mL.p << " " << mL.identity << "\n";
					}
					std::cout << "\n" << std::flush;
					printedReads++;
				}				
			}
			
			assert(abs(1 - p_read) <= 1e-3);
			ll_thisIteration_mapped += log(l_read);
		};

		callBackForAllReads(mappedFile, processOneRead);
		std::cout << "\n";
		
		double f_sum_preNormalization = 0;
		for(auto tF : f_nextIteration.first)
		{
			f_sum_preNormalization += tF.second;
		}
		for(auto tF : f_nextIteration.second)
		{
			f_sum_preNormalization += tF.second;
		}
		if(!abs(nMapped - f_sum_preNormalization) <= 1e-2)
		{
			std::cerr << "nMapped: " << nMapped << "\n";
			std::cerr << "f_sum_preNormalization: " << f_sum_preNormalization << "\n";
			std::cerr << std::flush;
		}
		assert(abs(nMapped - f_sum_preNormalization) <= 1e-2);

		// todo
		/*
		std::cout << "allReads_has_1654" << "\t" << allReads_has_1654 << "\n";
		std::cout << "allReads_has_209" << "\t" << allReads_has_209 << "\n";
		std::cout << "allReads_has_1654_and_1301" << "\t" << allReads_has_1654_and_1301 << "\n" << std::flush;
		// assert(2 == 3);
			
		std::string taxonID_print = "1654";

		std::cout << "Intermediate stats round " << EMiteration << "\n";
		std::cout << "1654 pre unmapped: " << f_nextIteration.second.at("1654") << "\n";
		std::cout << "209 pre unmapped: " << f_nextIteration.second.at("209") << "\n";
		std::cout << "52773 pre unmapped: " << f_nextIteration.first.at("52773") << "\n";
		*/
		
		double ll_thisIteration_unmapped = 0;
		if(0 && EMiteration >= round_first_unknown)
		{
			if(EMiteration > round_first_unknown)
			{
				// assert(2 == 3);
			}
			
			std::cout << "\tEM round " << EMiteration << ", " << unmappedReadsLengths.size() << " unmapped reads.\n" << std::flush;
			
			std::map<std::string, double> taxonIDs_assignedUnmapped;
			bool firstRead = true;
			for(auto oneUnmappedReadLength : unmappedReadsLengths)
			{
				std::map<std::string, double> unmapped_p; 
				double unmapped_p_sum = 0;

				for(auto taxonIDf : f.second)
				{
					double p_unmapped = taxonIDf.second * iM.getIdentityP(0, taxonIDf.first, oneUnmappedReadLength, false);
					unmapped_p[taxonIDf.first] = p_unmapped;
					unmapped_p_sum += p_unmapped;
				}

				for(auto& tID : unmapped_p)
				{
					tID.second /= unmapped_p_sum;
				}
				
				// todo
				/*
				if(unmapped_p.count(taxonID_print) && firstRead)
				{
					std::cout << "One read of length " << oneUnmappedReadLength << "\n";
					
					std::string maxP;
					double maxP_value;
					for(auto tID : unmapped_p)
					{		
						if(tID.second > 0.001)
						{
							std::cout << "\t" << tID.first << " " << tID.second << " [" << f.second.at(tID.first) << " x " << iM.getIdentityP(0, tID.first, oneUnmappedReadLength, false) << "]\n";
						}
						if((maxP.length() == 0) || (maxP_value < tID.second))
						{
							maxP = tID.first;
							maxP_value = tID.second;
						}
					}
					std::cout << "\t MAX " << maxP << " " << maxP_value << "\n";
					std::cout << std::flush;
				}
				*/

				for(auto tID : unmapped_p)
				{
					f_nextIteration.second.at(tID.first) += tID.second;
					if(taxonIDs_assignedUnmapped.count(tID.first) == 0)
					{
						taxonIDs_assignedUnmapped[tID.first] = 0;
					}
					taxonIDs_assignedUnmapped.at(tID.first) += tID.second;
				}
				
				ll_thisIteration_unmapped += log(unmapped_p_sum);
				
				firstRead = false;
			}
			
			printSorted(taxonIDs_assignedUnmapped, "Unmapped read assignment:", 10);
		}

		// todo remove
		/*
		std::cout << "1654 post unmapped: " << f_nextIteration.second.at("1654") << "\n";
		std::cout << "209 post unmapped: " << f_nextIteration.second.at("209") << "\n";
		std::cout << "52773 post unmapped: " << f_nextIteration.first.at("52773") << "\n";
		*/
		normalize_f(f_nextIteration);
		/*
		std::cout << "1654: " << f_nextIteration.second.at("1654") << "\n" << std::flush;
		std::cout << "209: " << f_nextIteration.second.at("209") << "\n" << std::flush;
		std::cout << "52773: " << f_nextIteration.first.at("52773") << "\n" << std::flush;
		std::cout << "\n";
		*/
		
		double ll_thisIteration = ll_thisIteration_mapped + ll_thisIteration_unmapped;

		std::cout << "\tLog likelihood: " << ll_thisIteration << std::endl;
		std::cout << "\t\tcontribution mapped reads  : " << ll_thisIteration_mapped << std::endl;
		std::cout << "\t\tcontribution unmapped reads: " << ll_thisIteration_unmapped << std::endl;
		
		//if((EMiteration > 0) && (EMiteration != round_first_unknown))
		if(EMiteration > 0)
		{
			double ll_diff = ll_thisIteration - ll_lastIteration;
			assert(ll_diff >= 0);

			double ll_relative = ll_thisIteration/ll_lastIteration;

			std::cout << "\tImprovement: " << ll_diff << std::endl;
			std::cout << "\tRelative   : " << ll_relative << std::endl;
		
			double ll_relative_imp = 1 - ll_relative;
			
			// todo 
			if((ll_diff < 20) || (EMiteration > 30))
			{
				// continueEM = false;
			}		
			
			if(ll_diff <= 1)
			{
				continueEM = false; 
			}
			
			if(ll_relative_imp < 0.0001)
			{
				continueEM = false;
			}
						
		/*
			if(abs(1-ll_relative) < 0.01)
			{
				continueEM = false;
			}
		*/
		}
		
		if(EMiteration > 10)
		{
			// continueEM = false;
		}
		
		if(0 && (EMiteration == (round_first_unknown - 1)))
		{
			// todo
			/*
			std::cout << "Before renormalization:\n";
			for(auto& taxonIDf : f_nextIteration.second)
			{
				if(taxonIDf.second > 0.001)
					std::cout << "\t" << taxonIDf.first << ": " << taxonIDf.second << "\n";
			}
			std::cout << "\n";
			*/
			
			for(auto& taxonIDf : f_nextIteration.second)
			{
				double p_unmapped = iM.getIdentityP(0, taxonIDf.first, -1, false, true);
				double p_mapped = 1 - p_unmapped;
				assert(p_mapped >= 0);
				assert(p_mapped <= 1);
				double increaseBy  = (p_mapped > 0) ? (1 / p_mapped) : 1;
				
				if(bestMappings_perTaxonID.count(taxonIDf.first) && (bestMappings_perTaxonID.at(taxonIDf.first) >= 10))
				{
					double newF = taxonIDf.second * increaseBy;
					// std::cout << "Re-normalize " << taxonIDf.first << " (p unmapped = " << p_unmapped << ") by factor " << increaseBy << " from " << taxonIDf.second << " to " << newF << "\n";
					taxonIDf.second = newF;
				}
			}
				
			normalize_f(f_nextIteration);
			
			/*
			std::cout << "After renormalization:\n";
			for(auto& taxonIDf : f_nextIteration.second)
			{
				if(taxonIDf.second > 0.001)				
					std::cout << "\t" << taxonIDf.first << ": " << taxonIDf.second << "\n";
			}
			std::cout << "\n" << std::flush;
			*/
		}
		
		double f_sum_postNormalization  = 0;
		for(auto tF : f_nextIteration.first)
		{
			f_sum_postNormalization += tF.second;
		}
		for(auto tF : f_nextIteration.second)
		{
			f_sum_postNormalization += tF.second;
		}	
		//std::cerr << "f_sum_postNormalization: " << f_sum_postNormalization << "\n" << std::flush;
		assert(abs(1 - f_sum_postNormalization) <= 1e-3);
		
		// std::cout << "For round " << EMiteration << "\n";
		// printGenusLevelSummary(T, f_nextIteration);
		
		f = f_nextIteration;
		EMiteration++;
		ll_lastIteration = ll_thisIteration;
		
		//std::string taxonID_print = "1654";
		/*
		std::cout << "\nNext iteration " << taxonID_print << ":\n";
		std::cout << "\tDirect:   " << (f_nextIteration.first.count(taxonID_print) ? f_nextIteration.first.at(taxonID_print) : 0) << "\n";
		std::cout << "\tInDirect: " << (f_nextIteration.second.count(taxonID_print) ? f_nextIteration.second.at(taxonID_print) : 0) << "\n";
		std::cout << "\n" << std::flush;
		*/
		//assert(1 == 4);
	}

	// output file names

	std::string output_assigned_reads_and_identities = mappedFile + ".U.lengthAndIdentitiesPerTaxonID";
	std::string output_pot_frequencies = mappedFile + ".U.WIMP";
	std::string output_reads_taxonID = mappedFile + ".U.reads2Taxon";
	std::string output_shifted_frequencies_file = mappedFile + ".U.shiftedHistogramsPerTaxonID";

	std::ofstream strout_reads_identities(output_assigned_reads_and_identities);
	assert(strout_reads_identities.is_open());
	strout_reads_identities << "taxonID" << "\t" << "directIndirect" << "\t" << "taxonName" << "\t" << "Identity" << "\t" << "Length" << "\n";

	std::ofstream strout_reads_taxonIDs(output_reads_taxonID);
	assert(strout_reads_taxonIDs.is_open());
	
	std::pair<std::map<std::string, size_t>, std::map<std::string, size_t>> assignedReads;

	std::function<void(const std::vector<std::string>&)> processOneRead_final = [&](const std::vector<std::string>& readLines) -> void
	{
		assert(readLines.size() > 0);
		std::vector<oneMappingLocation_U> mappingLocations = getMappingLocations_U(iM, indirectUpwardNodes, f, readLines);

		std::string readID = getReadIDFromReadLines(readLines);
		
		oneMappingLocation_U bestMapping = getBestMapping_U(mappingLocations);

		std::map<std::string, size_t>& assignedReads_dirIndir = (bestMapping.direct) ? assignedReads.first : assignedReads.second;
		if(assignedReads_dirIndir.count(bestMapping.taxonID) == 0)
		{
			assignedReads_dirIndir[bestMapping.taxonID] = 0;
		}
		assignedReads_dirIndir.at(bestMapping.taxonID)++;

		strout_reads_identities <<
				bestMapping.taxonID << "\t" <<
				((bestMapping.direct) ? "direct" : "indirect") << "\t" <<
				T.getNode(bestMapping.taxonID).name.scientific_name << "\t" <<
				bestMapping.identity << "\t" <<
				bestMapping.readLength << "\n";
				
		strout_reads_taxonIDs << readID << "\t" << bestMapping.taxonID << "\n";
				
	};

	callBackForAllReads(mappedFile, processOneRead_final);

	strout_reads_identities.close();
	strout_reads_taxonIDs.close();
	
	cleanF_U(f, assignedReads, mappingStats.at("ReadsMapped"));

	std::tuple<std::map<std::string, double>, std::map<std::string, double>, std::map<std::string, double>> frequencies_triplet;
	if(unmappedReadsLengths.size())
	{
		std::cerr << "Analysis of unmapped reads:\n";
		
		double proportion_mapped = (double)nMapped / (double)nReadsMappable;
		double proportion_toDistribute = (double)nUnmapped / (double)nReadsMappable;
		assert(abs((proportion_mapped + proportion_toDistribute) - 1) <= 1e-3);

		std::cerr << "\t" << "proportion_mapped" << ": " << proportion_mapped << "\n";
		std::cerr << "\t" << "proportion_toDistribute" << ": " << proportion_toDistribute << "\n";

		// find average proportion of expected unmapped reads, at given unmapped read lengths
		std::map<std::string, double> proportion_unmapped_averagedEmpiricalReadLengths;
		std::cerr << "\t" << "Unmapped for - expected\n";
		
		for(auto tF : f.second)
		{
			double proportion_unmapped_sum = 0;
			for(auto oneUnmappedReadLength : unmappedReadsLengths)
			{
				proportion_unmapped_sum += iM.getIdentityP(0, tF.first, oneUnmappedReadLength, false);
			}
			double proportion_unmapped_average = proportion_unmapped_sum / (double)unmappedReadsLengths.size();
			assert((proportion_unmapped_average >= 0) && (proportion_unmapped_average <= 1));
			proportion_unmapped_averagedEmpiricalReadLengths[tF.first] = proportion_unmapped_average;
			std::cerr << "\t\t" << "Taxon ID " << tF.first << ": indirectly " << tF.second << "; unmapped for " << tF.first << ": " << proportion_unmapped_average << "\n";
		}

		
		std::map<std::string, double> totalFrequencies_before;
		for(auto tF : f.first)
		{
			if(totalFrequencies_before.count(tF.first) == 0)
			{
				totalFrequencies_before[tF.first] = 0;
			}
			totalFrequencies_before.at(tF.first) += tF.second;
		}	
		for(auto tF : f.second)
		{
			if(totalFrequencies_before.count(tF.first) == 0)
			{
				totalFrequencies_before[tF.first] = 0;
			}
			totalFrequencies_before.at(tF.first) += tF.second;
		}	
		
		std::cerr << "\t" << "Would like to add:\n";
		
		// find how many unmapped reads we might want to add
		double allTaxa_wantAdditionalReads = 0;
		std::map<std::string, double> perTaxon_wantAdditionalReads;
		for(auto tF : f.second)
		{
			double currentIndirectFrequency = tF.second;
			assert((currentIndirectFrequency >= 0) && (currentIndirectFrequency <= 1));
			double approximate_read_number = nMapped * currentIndirectFrequency;
			double expected_proportion_mapped = 1 - proportion_unmapped_averagedEmpiricalReadLengths.at(tF.first);
			double wouldLikeTotal = (1.0/expected_proportion_mapped) * approximate_read_number;
			double wouldLikeToAdd = wouldLikeTotal - approximate_read_number;
			assert(wouldLikeToAdd >= 0);
			allTaxa_wantAdditionalReads += wouldLikeToAdd;
			perTaxon_wantAdditionalReads[tF.first] = wouldLikeToAdd;
			
			std::cerr << "\t\t" << "Taxon ID " << tF.first  << ": " << approximate_read_number << " reads indirectly assigned; expected mapped: " << expected_proportion_mapped << "; desireed total: " << wouldLikeTotal << "; now would like to add: " << wouldLikeToAdd << "\n";
		}

		
		// scale the number of reads we can add
		double addReads_scalingFactor = 1;
		if(allTaxa_wantAdditionalReads > nUnmapped)
		{
			addReads_scalingFactor = nUnmapped /  allTaxa_wantAdditionalReads;
		}
		
		std::cerr << "\t" << "Total reads we'd like to add: " << allTaxa_wantAdditionalReads << "; can spread (max.): " << nUnmapped << "; scaling factor: " << addReads_scalingFactor << "\n";

		double leaveUnassigned = nUnmapped - (allTaxa_wantAdditionalReads * addReads_scalingFactor);
		if(!(leaveUnassigned >= -1e-5))
		{
			std::cerr << "Warning: leaveUnassigned = " << leaveUnassigned << "\n" << std::flush;
		}
		assert(leaveUnassigned >= -1e-5);
		if(leaveUnassigned < 0)
		{
			leaveUnassigned = 0; 
		}
		double leaveUnassignedProp = leaveUnassigned / nReadsMappable;
		assert((leaveUnassignedProp >= 0) && (leaveUnassignedProp <= 1));

		std::cerr << "\t" << "Leave unassigned: " << leaveUnassigned << "; prop: " << leaveUnassignedProp << "\n" << std::flush;
		
		// now re-assign reads across taxa
		for(auto& tF : f.first)
		{
			std::get<0>(frequencies_triplet)[tF.first] = tF.second * nMapped;
		}

		for(auto& tF : f.second)
		{
			std::get<1>(frequencies_triplet)[tF.first] = tF.second * nMapped;
			std::get<2>(frequencies_triplet)[tF.first] = addReads_scalingFactor * perTaxon_wantAdditionalReads.at(tF.first);
		}

		// check that the total number of reads is OK
		double total_reads_post_unmapped = 0;
		for(auto tF : std::get<0>(frequencies_triplet))
		{
			total_reads_post_unmapped += tF.second;
		}
		for(auto tF : std::get<1>(frequencies_triplet))
		{
			total_reads_post_unmapped += tF.second;
		}
		for(auto tF : std::get<2>(frequencies_triplet))
		{
			total_reads_post_unmapped += tF.second;
		}

		assert(abs((total_reads_post_unmapped + leaveUnassigned) - nReadsMappable ) <= 1e-3);

		normalize_f_triplet(frequencies_triplet);

				
		std::map<std::string, double> totalFrequencies_after;
		double checkSum_postUnassignedRemoval = 0;
		for(auto& tF : std::get<0>(frequencies_triplet))
		{
			tF.second *= (1 - leaveUnassignedProp);
			checkSum_postUnassignedRemoval += tF.second;
			if(totalFrequencies_after.count(tF.first) == 0)
			{
				totalFrequencies_after[tF.first] = 0;
			}
			totalFrequencies_after.at(tF.first) += tF.second;
		}
		for(auto& tF : std::get<1>(frequencies_triplet))
		{
			tF.second *= (1 - leaveUnassignedProp);
			checkSum_postUnassignedRemoval += tF.second;
			if(totalFrequencies_after.count(tF.first) == 0)
			{
				totalFrequencies_after[tF.first] = 0;
			}
			totalFrequencies_after.at(tF.first) += tF.second;			
		}
		for(auto& tF : std::get<2>(frequencies_triplet))
		{
			tF.second *= (1 - leaveUnassignedProp);
			checkSum_postUnassignedRemoval += tF.second;
			if(totalFrequencies_after.count(tF.first) == 0)
			{
				totalFrequencies_after[tF.first] = 0;
			}
			totalFrequencies_after.at(tF.first) += tF.second;			
		}
		assert(abs(1 - (checkSum_postUnassignedRemoval + leaveUnassignedProp)) <= 1e-3);
		
		assert(totalFrequencies_before.size() == totalFrequencies_after.size());
		std::cerr << "\tTotal freq change summary\n";
		for(auto tF : totalFrequencies_before)
		{
			std::cerr << "\t\tTaxon ID " << tF.first << ": " << totalFrequencies_before.at(tF.first) << " -> " << totalFrequencies_after.at(tF.first) << "\n";
		}
		std::cerr << std::flush;
		std::cout << std::flush;
	}
	else
	{
		for(auto tF : f.first)
		{
			std::get<0>(frequencies_triplet)[tF.first] = tF.second;
		}
		for(auto tF : f.second)
		{
			std::get<1>(frequencies_triplet)[tF.first] = tF.second;
		}

		normalize_f_triplet(frequencies_triplet);
	}


	/*
	for(auto oneUnmappedReadLength : unmappedReadsLengths)
	{
		std::map<std::string, double> unmapped_p;
		// double unmapped_p_sum = 0;

		double best_unmapped_P;
		std::string best_unmapped_P_which;
		for(auto taxonIDf : f.second)
		{
			double p_unmapped = taxonIDf.second * iM.getIdentityP(0, taxonIDf.first, oneUnmappedReadLength, false);
			if((best_unmapped_P_which.length() == 0) || (best_unmapped_P < p_unmapped))
			{
				best_unmapped_P = p_unmapped;
				best_unmapped_P_which = taxonIDf.first;
			}
		}

		assignedReads.second[best_unmapped_P_which]++;
	}
	*/



	producePotFile_U(output_pot_frequencies, T, frequencies_triplet, assignedReads, nReadsMappable, mappableTaxonIDs);
	
	produceShiftedHistograms(output_shifted_frequencies_file, iM, f);
}

void cleanF_U(std::pair<std::map<std::string, double>, std::map<std::string, double>>& f, const std::pair<std::map<std::string, size_t>, std::map<std::string, size_t>>& assignedReads, size_t distributedReads)
{
	double minFreq = 0.9 * (1.0/(double)distributedReads);
	std::set<std::string> taxonID_to_delete;
	
	std::map<std::string, double> f_combined;
	
	for(auto fI : f.first)
	{
		if(f_combined.count(fI.first) == 0)
		{
			f_combined[fI.first] = 0;
		}
		f_combined.at(fI.first) += fI.second;
	}
	for(auto fI : f.second)
	{
		if(f_combined.count(fI.first) == 0)
		{
			f_combined[fI.first] = 0;
		}
		f_combined.at(fI.first) += fI.second;
	}	
	
	for(auto fI : f_combined)
	{
		if((fI.second < minFreq) && (assignedReads.first.count(fI.first) == 0) && (assignedReads.second.count(fI.first) == 0))
		{
			taxonID_to_delete.insert(fI.first);
		}
	}
	
	for(auto taxonID : taxonID_to_delete)
	{
		f.first.erase(taxonID);
		f.second.erase(taxonID);
	}
	
	double f_sum = 0;
	for(auto fI : f.first)
	{
		f_sum += fI.second;
	}
	for(auto fI : f.second)
	{
		f_sum += fI.second;
	}	
	assert(f_sum > 0);
	for(auto& fI : f.first)
	{
		fI.second /= f_sum;
	}
	for(auto& fI : f.second)
	{
		fI.second /= f_sum;
	}

	double f_sum_2 = 0;
	for(auto fI : f.first)
	{
		f_sum_2 += fI.second;
	}
	for(auto fI : f.second)
	{
		f_sum_2 += fI.second;
	}	
	assert(abs(1 - f_sum_2) <= 1e-3);
}


std::string getReadIDFromReadLines(const std::vector<std::string>& readLines)
{
	assert(readLines.size());
	std::string readID;
	for(unsigned int lineI = 0; lineI < readLines.size(); lineI++)
	{
		std::vector<std::string> line_fields = split(readLines.at(lineI), " ");	
		if(lineI == 0)
		{
			readID = line_fields.at(0);
		}
		else
		{
			assert(readID == line_fields.at(0));
		}
	}
	
	return readID;
}		

}
#endif /* META_FU_H_ */
