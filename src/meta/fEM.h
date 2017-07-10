/*
 * fEM.h
 *
 *  Created on: Jul 7, 2017
 *      Author: diltheyat
 */

#ifndef META_FEM_H_
#define META_FEM_H_

#include <map>
#include <set>
#include <utility>
#include <string>
#include <fstream>
#include <assert.h>
#include <boost/regex.hpp>

#include "util.h"
#include "taxonomy.h"

namespace meta
{
std::string extractTaxonId(std::string line);
std::map<std::string, std::map<std::string, size_t>> loadRelevantTaxonInfo(std::string DBdir, const std::set<std::string>& taxonIDs);
std::set<std::string> getTaxonIDsFromMappingsFile(std::string mappedFile);
void callBackForAllReads(std::string mappedFile, std::function<void(const std::vector<std::string>&)>& callbackFunction);

class oneMappingLocation
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
};

void producePotFile(std::string outputFN, const taxonomy& T, std::map<std::string, double> frequencies, std::map<std::string, size_t> readCount, size_t nUnmapped, size_t nTooShort)
{
	std::set<std::string> combinedKeys;
	for(auto f : frequencies)
	{
		combinedKeys.insert(f.first);
	}
	for(auto f : readCount)
	{
		combinedKeys.insert(f.first);
	}

	std::map<std::string, std::set<std::string>> combinedKeys_perLevel;

	std::map<std::string, std::map<std::string, double>> f_per_level;
	for(auto freq_per_node : frequencies)
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

oneMappingLocation getBestMapping(const std::vector<oneMappingLocation>& locations)
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
std::vector<oneMappingLocation> getMappingLocations(const std::map<std::string, std::map<std::string, size_t>>& taxonInfo, const std::map<std::string, double>& f, const std::vector<std::string>& readLines)
{
	assert(readLines.size() > 0);

	std::set<std::string> saw_contigIDs;
	std::set<std::string> saw_taxonIDs;

	std::vector<oneMappingLocation> mappingLocations;

	std::string readID;
	long long readLength = -1;
	for(auto line : readLines)
	{
		std::vector<std::string> line_fields = split(line, " ");

		std::string contigID = line_fields.at(5);
		size_t contigID_start = std::stoull(line_fields.at(7));
		size_t contigID_stop = std::stoull(line_fields.at(8));

		std::string contig_taxonID = extractTaxonId(contigID);

		assert(taxonInfo.at(contig_taxonID).count(contigID));

		double mappingQuality = std::stod(line_fields.at(13));
		assert(mappingQuality >= 0);
		assert(mappingQuality <= 1);

		if(readLength == -1)
			readLength = std::stoi(line_fields.at(1));

		if(readID.length() == 0)
			readID = line_fields.at(0);
		else
			assert(readID == line_fields.at(0));

		double identity = (std::stoi(line_fields.at(12)))/100.0;

		saw_taxonIDs.insert(contig_taxonID);
		saw_contigIDs.insert(contigID);

		oneMappingLocation l;
		l.taxonID = contig_taxonID;
		l.contigID = contigID;
		l.start = contigID_start;
		l.stop = contigID_stop;
		l.p = mappingQuality;
		l.readLength = readLength;
		l.identity = identity;
		mappingLocations.push_back(l);
	}

	std::map<std::string, size_t> mappingLocations_per_taxonID;
	for(auto targetTaxonID : saw_taxonIDs)
	{
		size_t possibleMappingLocations = 0;
		for(auto possibleContig : taxonInfo.at(targetTaxonID))
		{
			std::string contigID = possibleContig.first;
			size_t contigLength = possibleContig.second;

			if(contigLength >= readLength)
			{
				possibleMappingLocations += (contigLength - readLength + 1);
			}
			else
			{
				if(saw_contigIDs.count(contigID))
				{
					possibleMappingLocations++;
				}
			}
		}

		mappingLocations_per_taxonID[targetTaxonID] = possibleMappingLocations;
		assert(possibleMappingLocations > 0);
	}

	double p_allLocations = 0;
	for(auto& mL : mappingLocations)
	{
		mL.l = f.at(mL.taxonID) * (1/(double)mappingLocations_per_taxonID.at(mL.taxonID)) * mL.p;
		p_allLocations += mL.l ;
	}

	assert(p_allLocations > 0);
	for(auto& mL : mappingLocations)
	{
		mL.p = mL.l / p_allLocations;
	}

	return mappingLocations;
}

std::map<std::string, size_t> getMappingStats(std::string mappedFile)
{
	std::map<std::string, size_t> forReturn;
	std::ifstream statsStream (mappedFile + ".meta");
	assert(statsStream.is_open());
	std::string line;
	while(statsStream.good())
	{
		std::getline(statsStream, line);
		eraseNL(line);
		if(line.length() == 0)
		{
			std::vector<std::string> line_fields = split(line, " ");
			assert(line_fields.size() == 2);
			forReturn[line_fields.at(0)] = std::stoull(line_fields.at(1));
		}
	}
	return forReturn;
}

void doEM(std::string DBdir, std::string mappedFile)
{
	std::set<std::string> relevantTaxonIDs = getTaxonIDsFromMappingsFile(mappedFile);

	// Get numbers of unmapped, short reads
	// todo
	std::map<std::string, size_t> mappingStats = getMappingStats(mappedFile);
	size_t nUnmapped = mappingStats.at("ReadsNotMapped");
	size_t nTooShort = mappingStats.at("ReadsTooShort");

	// get genome info
	std::map<std::string, std::map<std::string, size_t>> taxonInfo = loadRelevantTaxonInfo(DBdir, relevantTaxonIDs);

	// load taxonomy
	taxonomy T(DBdir+"/taxonomy");

	// initialize frequency distribution
	std::map<std::string, double> f;
	for(auto t : relevantTaxonIDs)
	{
		f[t] = 1/(double)relevantTaxonIDs.size();
	}

	double ll_lastIteration;
	size_t EMiteration = 0;
	bool continueEM = true;
	while(continueEM)
	{
		double ll_thisIteration = 0;
		std::map<std::string, double> f_nextIteration = f;
		for(auto& fNextEntry : f_nextIteration)
		{
			fNextEntry.second = 0;
		}

		std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines) -> void
		{
			assert(readLines.size() > 0);
			std::vector<oneMappingLocation> mappingLocations = getMappingLocations(taxonInfo, f, readLines);

			double l_read = 0;
			for(const auto& mL : mappingLocations)
			{
				l_read += mL.l;
				f_nextIteration.at(mL.taxonID) += mL.p;
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

	std::string output_pot_frequencies = mappedFile + ".WIMP";
	std::string output_recalibrated_mappings = mappedFile + ".postEM";
	std::string output_assigned_reads_and_identities = mappedFile + ".lengthAndIdentitiesPerMappingUnit";
	std::string output_contig_coverage = mappedFile + ".contigCoverage";

	std::ofstream strout_reads_identities(output_assigned_reads_and_identities);
	assert(strout_reads_identities.is_open());
	strout_reads_identities << "AnalysisLevel" << "\t" << "ID" << "\t" << "readI" << "\t" << "Identity" << "\t" << "Length" << "\n";

	std::ofstream strout_recalibrated_mappings(output_recalibrated_mappings);
	assert(strout_recalibrated_mappings.is_open());

	size_t coverage_windowSize = 1000;
	std::map<std::string, std::map<std::string, std::vector<size_t>>> coverage_per_contigID;
	std::map<std::string, std::map<std::string, size_t>> size_last_window;

	std::map<std::string, size_t> reads_per_taxonID;
	size_t runningReadI = 0;
	std::function<void(const std::vector<std::string>&)> processOneRead_final = [&](const std::vector<std::string>& readLines) -> void
	{
		assert(readLines.size() > 0);
		std::vector<oneMappingLocation> mappingLocations = getMappingLocations(taxonInfo, f, readLines);

		assert(readLines.size() == mappingLocations.size());
		for(unsigned int lineI = 0; lineI < readLines.size(); lineI++)
		{
			double finalMappingQuality = mappingLocations.at(lineI).p;
			std::vector<std::string> line_fields = split(readLines.at(lineI), " ");
			line_fields.at(13) = std::to_string(finalMappingQuality);
			strout_recalibrated_mappings << join(line_fields, " ") << "\n";
		}

		oneMappingLocation bestMapping = getBestMapping(mappingLocations);
		strout_reads_identities << "EqualCoverageUnit" << "\t" << bestMapping.contigID << "\t" << runningReadI << "\t" << bestMapping.identity << "\t" << bestMapping.readLength << "\n";
		if(reads_per_taxonID.count(bestMapping.taxonID) == 0)
		{
			reads_per_taxonID[bestMapping.taxonID] = 0;
		}
		reads_per_taxonID[bestMapping.taxonID]++;

		if(coverage_per_contigID[bestMapping.taxonID].count(bestMapping.contigID) == 0)
		{
			size_t contigLength = taxonInfo.at(bestMapping.taxonID).at(bestMapping.contigID);
			size_t n_windows = contigLength / coverage_windowSize;
			if(n_windows == 0)
			{
				n_windows++;
				size_last_window[bestMapping.taxonID][bestMapping.contigID] = contigLength;
			}
			else
			{
				if((n_windows * coverage_windowSize) != contigLength)
				{
					n_windows++;
					size_last_window[bestMapping.taxonID][bestMapping.contigID] = contigLength - (n_windows * coverage_windowSize);
				}
			}
			coverage_per_contigID[bestMapping.taxonID][bestMapping.contigID].resize(n_windows, 0);
		}

		size_t bestMappingStopPos = (bestMapping.stop >= taxonInfo.at(bestMapping.taxonID).at(bestMapping.contigID)) ? (taxonInfo.at(bestMapping.taxonID).at(bestMapping.contigID) - 1) : bestMapping.stop;
		for(size_t contigPos = bestMapping.start; contigPos <= bestMappingStopPos; contigPos += coverage_windowSize)
		{
			size_t windowI = contigPos / coverage_windowSize;
			assert(windowI >= 0);
			assert(windowI < coverage_per_contigID.at(bestMapping.taxonID).at(bestMapping.contigID).size());

			size_t windowStart = windowI * coverage_windowSize;
			size_t windowStop = (windowI+1) * coverage_windowSize - 1;
			if(windowStop > taxonInfo.at(bestMapping.taxonID).at(bestMapping.contigID))
			{
				windowStop = taxonInfo.at(bestMapping.taxonID).at(bestMapping.contigID) - 1;
			}
			assert(windowStop >= windowStart);

			size_t read_window_overlap = overlap(windowStart, windowStop, bestMapping.start, bestMappingStopPos);
			assert(read_window_overlap > 0);
			assert(read_window_overlap <= coverage_windowSize);
			coverage_per_contigID.at(bestMapping.taxonID).at(bestMapping.contigID).at(windowI) += read_window_overlap;
		}

		runningReadI++;
	};

	callBackForAllReads(mappedFile, processOneRead_final);
	strout_reads_identities.close();

	producePotFile(output_pot_frequencies, T, f, reads_per_taxonID, nUnmapped, nTooShort);

	std::ofstream strout_coverage(output_contig_coverage);
	strout_coverage << "taxonID" << "\t" << "equalCoverageUnitLabel" << "\t" << "contigID" << "\t" << "start" << "\t" << "stop" << "\t" << "nBases" << "\t" << "readCoverage" << "\n";

	for(auto taxonData : coverage_per_contigID)
	{
		for(auto contigData : taxonData.second)
		{
			for(size_t windowI = 0; windowI < contigData.second.size(); windowI++)
			{
				size_t windowLength = coverage_windowSize;
				if((windowI == (contigData.second.size()-1)) && (size_last_window.at(taxonData.first).at(contigData.first) != 0))
				{
					windowLength = size_last_window.at(taxonData.first).at(contigData.first);
				}

				size_t nBases = coverage_per_contigID.at(taxonData.first).at(contigData.first).at(windowI);
				strout_coverage <<
						taxonData.first << "\t" <<
						T.getNode(taxonData.first).name.scientific_name << "\t" <<
						contigData.first << "\t" <<
						nBases << "\t" <<
						(double)nBases/(double)windowLength << "\n";
			}
		}
	}


}

void callBackForAllReads(std::string mappedFile, std::function<void(const std::vector<std::string>&)>& callbackFunction)
{
	std::ifstream mappingsStream (mappedFile);
	assert(mappingsStream.is_open());

	std::string runningReadID;
	std::vector<std::string> runningReadLines;

	std::string line;
	while(mappingsStream.is_open())
	{
		std::getline(mappingsStream, line);
		eraseNL(line);
		if(line.length() == 0)
			continue;

		size_t firstSpacePos = line.find(" ");
		assert(firstSpacePos != std::string::npos);

		std::string readID = line.substr(0, firstSpacePos);

		if(runningReadID != readID)
		{
			callbackFunction(runningReadLines);

			runningReadID = readID;
			runningReadLines.clear();
		}

		runningReadLines.push_back(line);
	}

	if(runningReadLines.size())
	{
		callbackFunction(runningReadLines);
	}
}

std::map<std::string, std::map<std::string, size_t>> loadRelevantTaxonInfo(std::string DBdir, const std::set<std::string>& taxonIDs)
{
	std::map<std::string, std::map<std::string, size_t>> forReturn;

	std::string fn_taxons = DBdir + "/taxonInfo.txt";
	std::ifstream tS (fn_taxons);
	if(! tS.is_open())
	{
		std::cerr << "Could not open file " << fn_taxons << " -- perhaps you have specified an incomplete DB?" << std::endl;
		exit(1);
	}
	assert(tS.is_open());

	std::string line;
	while(tS.good())
	{
		std::getline(tS, line);
		eraseNL(line);
		if(line.length() == 0)
			continue;
		std::vector<std::string> line_fields = split(line, " ");
		assert(line_fields.size() == 2);

		std::string taxonID = line_fields.at(0);

		std::vector<std::string> contigs = split(line_fields.at(1), ";");

		for(auto c : contigs)
		{
			std::vector<std::string> contigFields = split(c, "=");
			assert(contigFields.size() == 2);
			assert(forReturn[taxonID].count(contigFields.at(0)) == 0);
			forReturn[taxonID][contigFields.at(0)] = std::stoull(contigFields.at(1));
		}
	}

	return forReturn;
}

std::set<std::string> getTaxonIDsFromMappingsFile(std::string mappedFile)
{
	std::set<std::string> forReturn;
	std::ifstream mappingsStream(mappedFile);
	assert(mappingsStream.is_open());
	std::string line;
	while(mappingsStream.good())
	{
		std::getline(mappingsStream, line);
		eraseNL(line);
		std::vector<std::string> line_fields = split(line, " ");
		std::string mappingTarget = line_fields.at(5);
		std::string taxonID = extractTaxonId(mappingTarget);
		forReturn.insert(taxonID);
	}
	mappingsStream.close();
	return forReturn;
}

boost::regex matchTaxonID{"kraken:taxid\\|(x?\\d+)"};
boost::smatch matchBrackets;
std::string extractTaxonId(std::string line)
{
	assert(boost::regex_search(line, matchBrackets, matchTaxonID));
	std::string taxonID = matchBrackets[0];
	return taxonID;
}

}

#endif /* META_FEM_H_ */
