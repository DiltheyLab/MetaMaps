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

#include "fEM.h"
#include "util.h"

#include "fU_helper.h"


namespace meta
{


void producePotFile_U(std::string outputFN, const taxonomy& T, std::pair<std::map<std::string, double>, std::map<std::string, double>> frequencies, std::pair<std::map<std::string, size_t>, std::map<std::string, size_t>> readCount)
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

	std::map<std::string, std::set<std::string>> combinedKeys_perLevel;

	std::map<std::string, std::map<std::string, double>> f_direct_per_level;
	for(auto freq_per_node : frequencies.first)
	{
		std::string nodeID = freq_per_node.first;
		assert(T.knowNode(nodeID));
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(nodeID);
		upwardByLevel["EqualCoverageUnit"] = nodeID;
		for(auto uN : upwardByLevel)
		{
			if(f_per_level[uN.first].count(uN.second) == 0)
				f_per_level[uN.first][uN.second] = 0;

			f_per_level[uN.first][uN.second] += freq_per_node.second;
			combinedKeys_perLevel[uN.first].insert(uN.second);

			assert(f_per_level[uN.first][uN.second] >= 0);
			assert(f_per_level[uN.first][uN.second] <= 1);
		}
	}
	std::map<std::string, std::map<std::string, size_t>> rC_per_level;
	for(auto count_per_node : readCount)
	{
		std::string nodeID = count_per_node.first;
		assert(T.knowNode(nodeID));
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(nodeID);
		upwardByLevel["EqualCoverageUnit"] = nodeID;
		for(auto uN : upwardByLevel)
		{
			if(f_per_level[uN.first].count(uN.second) == 0)
				rC_per_level[uN.first][uN.second] = 0;

			rC_per_level[uN.first][uN.second] += count_per_node.second;
			combinedKeys_perLevel[uN.first].insert(uN.second);
		}
	}


	std::ofstream strout_frequencies(outputFN);
	assert(strout_frequencies.is_open());

	strout_frequencies << "AnalysisLevel" << "\t" <<  "ID" << "\t" << "Name" << "\t" <<  "Absolute" << "\t" <<  "PotFrequency" << "\n";

	for(auto l : combinedKeys_perLevel)
	{
		std::string levelName = l.first;

		for(auto taxonID : l.second)
		{
			std::string taxonIDName = T.getNode(taxonID).name.scientific_name;
			double f = (f_per_level[levelName].count(taxonID)) ? f_per_level[levelName][taxonID] : 0;
			size_t rC = (rC_per_level[levelName].count(taxonID)) ? rC_per_level[levelName][taxonID] : 0;

			strout_frequencies << levelName << "\t" << taxonID << "\t" << taxonIDName << "\t" << rC << "\t" << f << "\n";
		}
		strout_frequencies << levelName << "\t" << 0 << "\t" <<  "Unmapped" << "\t" << nUnmapped << "\t" << 0 << "\n";
		strout_frequencies << levelName << "\t" << 0 << "\t" <<  "TooShort" << "\t" << nTooShort << "\t" << 0 << "\n";

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
		l.l = f.first.at(contig_taxonID) * iM.getIdentityP(identityInt, contig_taxonID, true);
		mappingLocations.push_back(l);

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
			lI.l = f.second.at(indirectTaxon) * iM.getIdentityP(identityInt, indirectTaxon, false);
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
	taxonomy T(DBdir+"/taxonomy");

	std::set<std::string> taxonIDsInMappings = getTaxonIDsFromMappingsFile(mappedFile);

	std::string fn_fittedLengthAndIdentities= mappedFile + ".lengthAndIdentitiesPerMappingUnit";
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
	size_t nUnmapped = mappingStats.at("ReadsNotMapped");
	size_t nTooShort = mappingStats.at("ReadsTooShort");

	// set up nodes to map to

	std::set<std::string> relevantTaxonIDs_direct = taxonIDsInMappings;
	std::set<std::string> relevantTaxonIDs_indirect;

	std::map<std::string, std::vector<std::string>> indirectUpwardNodes;
	for(auto tI : taxonIDsInMappings)
	{
		std::vector<std::string> upwardTaxonIDs = T.getUpwardNodes(tI);
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

	// set up EM

	double ll_lastIteration;
	size_t EMiteration = 0;
	bool continueEM = true;
	while(continueEM)
	{
		double ll_thisIteration = 0;
		std::pair<std::map<std::string, double>, std::map<std::string, double>> f_nextIteration = f;
		for(auto& fNextEntry : f_nextIteration.first)
		{
			fNextEntry.second = 0;
		}
		for(auto& fNextEntry : f_nextIteration.second)
		{
			fNextEntry.second = 0;
		}

		std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines) -> void
		{
			assert(readLines.size() > 0);
			std::vector<oneMappingLocation_U> mappingLocations = getMappingLocations_U(iM, indirectUpwardNodes, f, readLines);

			double l_read = 0;
			for(const auto& mL : mappingLocations)
			{
				l_read += mL.l;
				if(mL.direct)
				{
					f_nextIteration.first.at(mL.taxonID) += mL.p;
				}
				else
				{
					f_nextIteration.second.at(mL.taxonID) += mL.p;
				}
			}
			ll_thisIteration += log(l_read);
		};

		callBackForAllReads(mappedFile, processOneRead);

		double ll_diff = ll_thisIteration - ll_lastIteration;
		assert(ll_diff >= 0);

		double ll_improvement = ll_thisIteration / ll_lastIteration;

		if((EMiteration > 0) && (ll_improvement < 1.1))
		{
			continueEM = false;
		}

		f = f_nextIteration;
		EMiteration++;
		ll_lastIteration = ll_thisIteration;
	}





}

}
#endif /* META_FU_H_ */
