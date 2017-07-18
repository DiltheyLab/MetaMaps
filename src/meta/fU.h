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


void producePotFile_U(std::string outputFN, const taxonomy& T, std::pair<std::map<std::string, double>, std::map<std::string, double>> frequencies, std::pair<std::map<std::string, size_t>, std::map<std::string, size_t>> readCount, size_t mappableReads)
{
	std::set<std::string> combinedKeys;
	for(auto f : frequencies.first)
	{
		combinedKeys.insert(f.first);
	}
	for(auto f : frequencies.second)
	{
		combinedKeys.insert(f.first);
	}
	for(auto f : readCount.first)
	{
		combinedKeys.insert(f.first);
	}
	for(auto f : readCount.second)
	{
		combinedKeys.insert(f.first);
	}

	std::map<std::string, std::pair<std::map<std::string, double>, std::map<std::string, double>>> frequencies_perLevel;
	std::map<std::string, std::pair<std::map<std::string, size_t>, std::map<std::string, size_t>>> readCount_perLevel;

	std::map<std::string, std::set<std::string>> combinedKeys_perLevel;

	std::set<std::string> targetLevels = getRelevantLevelNames();
	
	for(auto taxonID : combinedKeys)
	{
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(taxonID, targetLevels);
		for(auto uN : upwardByLevel)
		{
			std::string level = uN.first;
			std::string levelValue = uN.second;
			combinedKeys_perLevel[level].insert(levelValue);

			if(frequencies_perLevel[level].first.count(levelValue) == 0)
			{
				frequencies_perLevel[level].first[levelValue] = 0;
				frequencies_perLevel[level].second[levelValue] = 0;
				readCount_perLevel[level].first[levelValue] = 0;
				readCount_perLevel[level].second[levelValue] = 0;
			}

			if(frequencies.first.count(taxonID))
			{
				frequencies_perLevel.at(level).first.at(levelValue) += frequencies.first.at(taxonID);
			}

			if(frequencies.second.count(taxonID))
			{
				frequencies_perLevel.at(level).second.at(levelValue) += frequencies.second.at(taxonID);
			}

			if(readCount.first.count(taxonID))
			{
				readCount_perLevel.at(level).first.at(levelValue) += readCount.first.at(taxonID);
			}

			if(readCount.second.count(taxonID))
			{
				readCount_perLevel.at(level).second.at(levelValue) += readCount.second.at(taxonID);
			}
		}
	}


	std::ofstream strout_frequencies(outputFN);
	assert(strout_frequencies.is_open());

	strout_frequencies << "AnalysisLevel" << "\t" <<  "taxonID" << "\t" << "Name" << "\t" << "Absolute_direct" << "\t" << "Absolute_indirect" << "\t" << "EMFrequency_direct" << "\t" << "EMfrequency_indirect" << "\t" <<  "Absolute" << "\t" << "EMFrequency" << "\t" << "PotFrequency" << "\n";

	for(auto l : combinedKeys_perLevel)
	{
		std::string levelName = l.first;
		double levelFreqSum = 0;
		size_t levelReadSum = 0;
		for(auto taxonID : l.second)
		{
			std::string taxonIDName = (taxonID != "Undefined") ? T.getNode(taxonID).name.scientific_name : taxonID;
			
			size_t taxonID_combinedReads = (readCount_perLevel.at(levelName).first.at(taxonID) + readCount_perLevel.at(levelName).second.at(taxonID));
			double taxonID_combinedFreq = (frequencies_perLevel.at(levelName).first.at(taxonID) + frequencies_perLevel.at(levelName).second.at(taxonID));
			
			strout_frequencies <<
					levelName << "\t" <<
					((taxonID != "Undefined") ? taxonID : "0") << "\t" <<
					taxonIDName << "\t" <<
					readCount_perLevel.at(levelName).first.at(taxonID) << "\t" <<
					readCount_perLevel.at(levelName).second.at(taxonID) << "\t" <<
					frequencies_perLevel.at(levelName).first.at(taxonID) << "\t" <<
					frequencies_perLevel.at(levelName).second.at(taxonID) << "\t" <<
					taxonID_combinedReads << "\t" <<
					taxonID_combinedFreq << "\t" << 
					taxonID_combinedFreq << "\n";

			levelReadSum += taxonID_combinedReads;					
			levelFreqSum += taxonID_combinedFreq;
		}
				
		long long levelReadsUnclassified = mappableReads - levelReadSum;
		assert(levelReadsUnclassified >= 0);
		
		assert((levelFreqSum >= 0) && (levelFreqSum <= 1));
		double levelFreqUnclassified = 1 - levelFreqSum;
		
		strout_frequencies <<
				levelName << "\t" <<
				0 << "\t" <<
				"Unclassified" << "\t" <<
				0 << "\t" <<
				levelReadsUnclassified << "\t" <<
				0 << "\t" <<
				levelFreqUnclassified << "\t" <<
				levelReadsUnclassified << "\t" <<
				levelFreqUnclassified << "\t" << 
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

	std::vector<oneMappingLocation_U> mappingLocations;

	std::string readID;
	long long readLength = -1;
	for(auto line : readLines)
	{
		std::vector<std::string> line_fields = split(line, " ");

		std::string contigID = line_fields.at(5);
		size_t contigID_start = std::stoull(line_fields.at(7));
		size_t contigID_stop = std::stoull(line_fields.at(8));

		std::string contig_taxonID = extractTaxonId(contigID);


		if(readLength == -1)
			readLength = std::stoi(line_fields.at(1));

		if(readID.length() == 0)
			readID = line_fields.at(0);
		else
			assert(readID == line_fields.at(0));

		double identity = std::stod(line_fields.at(12))/100.0;
		assert(identity >= 0);
		assert(identity <= 1);
		int identityInt = int((identity * 100) + 0.5);

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
		mappingLocations.push_back(l);

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
			mappingLocations.push_back(l);
		}
	}

	double l_sum = 0;
	for(const auto& l : mappingLocations)
	{
		l_sum += l.l;
	}
	assert(l_sum > 0);
	for(auto& l : mappingLocations)
	{
		l.p = l.l/l_sum;
	}

	return mappingLocations;
}

void doU(std::string DBdir, std::string mappedFile, size_t minimumReadsPerBestContig)
{
	unsigned int round_first_unknown = 0;

	taxonomy T(DBdir+"/taxonomy");

	std::set<std::string> taxonIDsInMappings = getTaxonIDsFromMappingsFile(mappedFile);

	std::string fn_fittedLengthAndIdentities= mappedFile + ".EM.lengthAndIdentitiesPerMappingUnit";
	if(! fileExists(fn_fittedLengthAndIdentities))
	{
		std::cerr << "\n\nERROR: File " << fn_fittedLengthAndIdentities << " not existing.\n\nThis file is generated automatically by the EM step. Run the EM step first.\n\n";
	}

	identityAndReadLengthHistogram iAndL;
	iAndL.readFromEMOutput(fn_fittedLengthAndIdentities, minimumReadsPerBestContig);

	std::string fn_tree_selfSimilarities = DBdir + "/selfSimilarities.txt";
	treeAdjustedIdentities tAI;
	tAI.readFromFile(fn_tree_selfSimilarities, taxonIDsInMappings, T);

	identityManager iM(iAndL, tAI);

	// mapping stats
	std::map<std::string, size_t> mappingStats = getMappingStats(mappedFile);
	size_t nTotalReads = mappingStats.at("TotalReads");	
	size_t nTooShort = mappingStats.at("ReadsTooShort");
	size_t nUnmapped = mappingStats.at("ReadsNotMapped");
	
	size_t nReadsMappable = nTotalReads - nTooShort;
	assert(nReadsMappable <= nTotalReads);
	
	std::vector<size_t> unmappedReadsLengths = getUnmappedReadsStats(mappedFile);
	assert(unmappedReadsLengths.size() == nUnmapped); // can remove later
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
			}
		}
	}

	// set up initial frequency distribution

	std::pair<std::map<std::string, double>, std::map<std::string, double>> f;
	for(auto nI : relevantTaxonIDs_direct)
	{
		f.first[nI] = 1.0/(double)relevantTaxonIDs_direct.size();
	}
	for(auto nI : relevantTaxonIDs_indirect)
	{
		f.second[nI] = 1.0/(double)relevantTaxonIDs_indirect.size();
	}

	// set up U-EM

	std::cout << "Starting EM..." << std::endl;	
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

		size_t processedRead = 0;		
		std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines) -> void
		{
			processedRead++;
			std::cout << "\r EM round " << EMiteration << ", mapped read << " << processedRead << " / " << mappingStats.at("ReadsMapped") << "   " << std::flush;
					
			assert(readLines.size() > 0);
			std::vector<oneMappingLocation_U> mappingLocations = getMappingLocations_U(iM, indirectUpwardNodes, f, readLines);

			double l_read = 0;
			for(const auto& mL : mappingLocations)
			{
				l_read += mL.l;
				if(mL.direct)
				{
					assert(f_nextIteration.first.count(mL.taxonID));					
					f_nextIteration.first.at(mL.taxonID) += mL.p;
				}
				else
				{
					assert(f_nextIteration.second.count(mL.taxonID));										
					f_nextIteration.second.at(mL.taxonID) += mL.p;
				}
			}
			ll_thisIteration_mapped += log(l_read);
		};

		callBackForAllReads(mappedFile, processOneRead);
		std::cout << "\n";
		
		double ll_thisIteration_unmapped = 0;
		if(EMiteration >= round_first_unknown)
		{
			std::cout << "\tEM round " << EMiteration << ", " << unmappedReadsLengths.size() << " unmapped reads.\n" << std::flush;

			
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

				for(auto tID : unmapped_p)
				{
					f_nextIteration.second.at(tID.first) += tID.second;
				}

				ll_thisIteration_unmapped += log(unmapped_p_sum);
			}
		}

		double f_sum = 0;
		for(auto tF : f_nextIteration.first)
		{
			f_sum += tF.second;
		}
		for(auto tF : f_nextIteration.second)
		{
			f_sum += tF.second;
		}

		for(auto& tF : f_nextIteration.first)
		{
			tF.second /= f_sum;
		}
		for(auto& tF : f_nextIteration.second)
		{
			tF.second /= f_sum;
		}

		double ll_thisIteration = ll_thisIteration_mapped + ll_thisIteration_unmapped;

		std::cout << "\tLog likelihood: " << ll_thisIteration << std::endl;
		std::cout << "\t\tcontribution mapped reads  : " << ll_thisIteration_mapped << std::endl;
		std::cout << "\t\tcontribution unmapped reads: " << ll_thisIteration_unmapped << std::endl;
		
		if((EMiteration > 0) && (EMiteration != round_first_unknown))
		{
			double ll_diff = ll_thisIteration - ll_lastIteration;
			assert(ll_diff >= 0);

			double ll_relative = ll_thisIteration/ll_lastIteration;

			std::cout << "\tImprovement: " << ll_diff << std::endl;
			std::cout << "\tRelative   : " << ll_relative << std::endl;
		
			if(abs(1-ll_relative) < 0.01)
			{
				continueEM = false;
			}
		}
		f = f_nextIteration;
		EMiteration++;
		ll_lastIteration = ll_thisIteration;
	}

	// output file names

	std::string output_assigned_reads_and_identities = mappedFile + ".U.lengthAndIdentitiesPerTaxonID";
	std::string output_pot_frequencies = mappedFile + ".U.WIMP";


	std::ofstream strout_reads_identities(output_assigned_reads_and_identities);
	assert(strout_reads_identities.is_open());
	strout_reads_identities << "taxonID" << "\t" << "directIndirect" << "\t" << "taxonName" << "\t" << "Identity" << "\t" << "Length" << "\n";

	std::pair<std::map<std::string, size_t>, std::map<std::string, size_t>> assignedReads;

	std::function<void(const std::vector<std::string>&)> processOneRead_final = [&](const std::vector<std::string>& readLines) -> void
	{
		assert(readLines.size() > 0);
		std::vector<oneMappingLocation_U> mappingLocations = getMappingLocations_U(iM, indirectUpwardNodes, f, readLines);

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
	};

	callBackForAllReads(mappedFile, processOneRead_final);

	strout_reads_identities.close();
	
	for(auto oneUnmappedReadLength : unmappedReadsLengths)
	{
		std::map<std::string, double> unmapped_p;
		double unmapped_p_sum = 0;

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

	cleanF_U(f, assignedReads, mappingStats.at("ReadsMapped"));

	producePotFile_U(output_pot_frequencies, T, f, assignedReads, nReadsMappable);

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
}

}
#endif /* META_FU_H_ */
