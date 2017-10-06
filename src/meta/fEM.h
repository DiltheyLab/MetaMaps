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

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>

namespace meta
{
std::string extractTaxonId(std::string line);
std::map<std::string, std::map<std::string, size_t>> loadRelevantTaxonInfo(std::string DBdir, const std::set<std::string>& taxonIDs);
std::set<std::string> getTaxonIDsFromMappingsFile(std::string mappedFile);
void callBackForAllReads(std::string mappedFile, std::function<void(const std::vector<std::string>&)>& callbackFunction);
std::set<std::string> getRelevantLevelNames();
void cleanF(std::map<std::string, double>& f, const std::map<std::string, size_t>& reads_per_taxonID, size_t distributedReads);
std::map<std::string, std::vector<size_t>> get_NS_per_window(std::string DBdir, size_t windowSize, const std::map<std::string, std::map<std::string, std::vector<size_t>>>& coverage_per_contigID);

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

void producePotFile(std::string outputFN, const taxonomy& T, std::map<std::string, double> frequencies, std::map<std::string, size_t> readCount, size_t nTotalReads, size_t nUnmapped, size_t nTooShort)
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

	std::set<std::string> targetLevels = getRelevantLevelNames();
	
	std::map<std::string, std::set<std::string>> combinedKeys_perLevel;

	std::map<std::string, std::map<std::string, double>> f_per_level;
	for(auto freq_per_node : frequencies)
	{
		std::string nodeID = freq_per_node.first;
		assert(T.knowNode(nodeID));
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(nodeID, targetLevels);
		upwardByLevel["mappingTarget"] = nodeID;
		
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
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(nodeID, targetLevels);
		upwardByLevel["mappingTarget"] = nodeID;
		for(auto uN : upwardByLevel)
		{
			if(f_per_level[uN.first].count(uN.second) == 0)
				rC_per_level[uN.first][uN.second] = 0;

			rC_per_level[uN.first][uN.second] += count_per_node.second;
			combinedKeys_perLevel[uN.first].insert(uN.second);
		}
	}

	long long nMappable = nTotalReads - nTooShort;
	assert(nMappable > 0);
	long long nMapped = nMappable - nUnmapped;
	assert(nMapped >= 0);
	
	std::ofstream strout_frequencies(outputFN);
	assert(strout_frequencies.is_open());

	strout_frequencies << "AnalysisLevel" << "\t" <<  "taxonID" << "\t" << "Name" << "\t" <<  "Absolute" <<  "\t" <<  "EMFrequency" << "\t" <<  "PotFrequency" << "\n";

	double freq_unmapped = (double)nUnmapped/(double)nMappable;
	assert(freq_unmapped >= 0);
	assert(freq_unmapped <= 1);

	for(auto l : combinedKeys_perLevel)
	{
		std::string levelName = l.first;

		size_t sum_reads_assigned_thisLevel = 0;
		double sum_assigned_thisLevel = 0;
		for(auto taxonID : l.second)
		{
			double f = (f_per_level[levelName].count(taxonID)) ? f_per_level[levelName][taxonID] : 0;
			size_t rC = (rC_per_level[levelName].count(taxonID)) ? rC_per_level[levelName][taxonID] : 0;
			
			sum_assigned_thisLevel += f;
			sum_reads_assigned_thisLevel += rC;
		}		
		
		assert(sum_assigned_thisLevel >= 0);
		assert(abs(1-sum_assigned_thisLevel) <= 1e-3);
		if(sum_assigned_thisLevel > 1)
			sum_assigned_thisLevel = 1;
		
		long long nReads_unassigned_level = nMapped - sum_reads_assigned_thisLevel;
		assert(nReads_unassigned_level >= 0);
		
		double additional_freq_unassigned = (1 - freq_unmapped) * (1 - sum_assigned_thisLevel);
		
		assert(additional_freq_unassigned >= 0);
		assert(additional_freq_unassigned <= 1);
		
		for(auto taxonID : l.second)
		{
			if(f_per_level[levelName].count(taxonID))
			{
				f_per_level[levelName].at(taxonID) /= sum_assigned_thisLevel;	
			}
		}
		
		double total_thisLevel_freq_unassigned = freq_unmapped + additional_freq_unassigned;
		assert(total_thisLevel_freq_unassigned >= 0);
		assert(total_thisLevel_freq_unassigned <= 1);
		
		double l_freq = 0;
		double l_freq_2 = 0;
		for(auto taxonID : l.second)
		{
			std::string taxonIDName = (taxonID != "Undefined") ? T.getNode(taxonID).name.scientific_name : taxonID;
			double f = (f_per_level[levelName].count(taxonID)) ? f_per_level[levelName][taxonID] : 0;
			size_t rC = (rC_per_level[levelName].count(taxonID)) ? rC_per_level[levelName][taxonID] : 0;

			double taxonID_rescaled_f = (f * (1-total_thisLevel_freq_unassigned));
			strout_frequencies <<
				levelName << "\t" <<
				((taxonID != "Undefined") ? taxonID : "0") << "\t" <<
				taxonIDName << "\t" <<
				rC << "\t" <<
				f << "\t" <<
				taxonID_rescaled_f << "\n";
				
			l_freq += taxonID_rescaled_f;
			l_freq_2 += f;
		}
		
		strout_frequencies << levelName << "\t" << 0 << "\t" <<  "Unclassified" << "\t" << nReads_unassigned_level << "\t" << 0 << "\t" << additional_freq_unassigned << "\n";
		strout_frequencies << levelName << "\t" << 0 << "\t" <<  "Unmapped" << "\t" << nUnmapped << "\t" << 0 << "\t" << freq_unmapped << "\n";
		strout_frequencies << levelName << "\t" << 0 << "\t" <<  "TooShort" << "\t" << nTooShort << "\t" << 0 << "\t" << 0 << "\n";
		
		l_freq += (freq_unmapped+additional_freq_unassigned);
		
		assert(abs(1 - l_freq) <= 1e-3);
		assert(abs(1 - l_freq_2) <= 1e-3);
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
	double mQ_sum = 0;
	for(auto line : readLines)
	{
		std::vector<std::string> line_fields = split(line, " ");

		std::string contigID = line_fields.at(5);
		size_t contigID_start = std::stoull(line_fields.at(7));
		size_t contigID_stop = std::stoull(line_fields.at(8));

		std::string contig_taxonID = extractTaxonId(contigID);

		if(! taxonInfo.count(contig_taxonID))
		{
			std::cerr << "Unknown taxonID " << contig_taxonID << "; please check that your mappings file was mapped against the database now specified." << std::endl;
			exit(1);
		}
		assert(taxonInfo.at(contig_taxonID).count(contigID));

		double mappingQuality = std::stod(line_fields.at(13));
		assert(mappingQuality >= 0);
		assert(mappingQuality <= 1);
		mQ_sum += mappingQuality;
		
		if(readLength == -1)
			readLength = std::stoi(line_fields.at(1));

		if(readID.length() == 0)
			readID = line_fields.at(0);
		else
			assert(readID == line_fields.at(0));

		// todo perhaps back to 12
		double identity = std::stod(line_fields.at(9))/100.0;
		assert(identity >= 0);
		assert(identity <= 1);

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
	assert(abs(1 - mQ_sum) <= 1e-3);

	std::map<std::string, size_t> mappingLocations_per_taxonID;
	for(auto targetTaxonID : saw_taxonIDs)
	{
		size_t possibleMappingLocations = 0;
		for(auto possibleContig : taxonInfo.at(targetTaxonID))
		{
			std::string contigID = possibleContig.first;
			size_t contigLength = possibleContig.second;

			if((long long)contigLength >= readLength)
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
		p_allLocations += mL.l;
	}

	assert(p_allLocations > 0);
	for(auto& mL : mappingLocations)
	{
		assert(mL.l >= 0);
		assert(mL.l <= 1);		
		mL.p = mL.l / p_allLocations;
	}
	
	double postNormalization_sum = 0;
	for(const auto& mL : mappingLocations)
	{
		postNormalization_sum += mL.p;
	}
	assert(abs(1 - postNormalization_sum) <= 1e-3);
	
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
		if(line.length() != 0)
		{
			std::vector<std::string> line_fields = split(line, " ");
			assert(line_fields.size() == 2);
			forReturn[line_fields.at(0)] = std::stoull(line_fields.at(1));
		}
	}
	return forReturn;
}

std::vector<size_t> getUnmappedReadsStats(std::string mappedFile)
{
	std::vector<size_t> forReturn;
	std::ifstream statsStream (mappedFile + ".meta.unmappedReadsLengths");
	assert(statsStream.is_open());
	std::string line;
	while(statsStream.good())
	{
		std::getline(statsStream, line);
		eraseNL(line);
		if(line.length() != 0)
		{
			size_t l = std::stoull(line);
			forReturn.push_back(l);
		}
	}
	return forReturn;
}

void doEM(std::string DBdir, std::string mappedFile, size_t minimumReadsPerBestContig)
{
	std::set<std::string> relevantTaxonIDs = getTaxonIDsFromMappingsFile(mappedFile);

	// Get numbers of unmapped, short reads
	std::map<std::string, size_t> mappingStats = getMappingStats(mappedFile);
	size_t nUnmapped = mappingStats.at("ReadsNotMapped");
	size_t nTooShort = mappingStats.at("ReadsTooShort");
	size_t nTotalReads = mappingStats.at("TotalReads");

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

	std::cout << "Starting EM..." << std::endl;	
	double ll_lastIteration;
	size_t EMiteration = 0;
	bool continueEM = true;
	while(continueEM)
	{
		double f_sum = 0;
		for(auto fE : f)
		{
			f_sum += fE.second;
		}
		assert(abs(1 - f_sum) <= 1e-3);
		
		std::cout << "EM round " << EMiteration << std::endl;
		
		double ll_thisIteration = 0;
		std::map<std::string, double> f_nextIteration = f;
		for(auto& fNextEntry : f_nextIteration)
		{
			fNextEntry.second = 0;
		}

		size_t processedRead = 0;
		std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines) -> void
		{
			processedRead++;
			std::cout << "\r EM round " << EMiteration << ", read << " << processedRead << " / " << mappingStats.at("ReadsMapped") << "   " << std::flush;
			
			assert(readLines.size() > 0);
			std::vector<oneMappingLocation> mappingLocations = getMappingLocations(taxonInfo, f, readLines);

			double l_read = 0;
			double mappingLocations_sum = 0;
			for(const auto& mL : mappingLocations)
			{
				l_read += mL.l;
				assert(f_nextIteration.count(mL.taxonID));
				f_nextIteration.at(mL.taxonID) += mL.p;
				mappingLocations_sum += mL.p;
			}
			assert(abs(1 - mappingLocations_sum) < 1e-3);
			ll_thisIteration += log(l_read);
		};

		callBackForAllReads(mappedFile, processOneRead);
		std::cout << "\n";
		std::cout << "\tLog likelihood: " << ll_thisIteration << std::endl;
		
		double sum_f_nextIteration = 0;
		for(auto fNextIterationE : f_nextIteration)
		{
			sum_f_nextIteration += fNextIterationE.second;
		}
		
		for(auto& fNextIterationE : f_nextIteration)
		{
			fNextIterationE.second /= sum_f_nextIteration;
		}
		
		double sum_postNormalization = 0;
		for(const auto& fNextIterationE : f_nextIteration)
		{
			sum_postNormalization += fNextIterationE.second;
		}
		assert(abs(1 - sum_postNormalization) <= 1e-3);
		
		if(EMiteration > 0) 
		{
			double ll_diff = ll_thisIteration - ll_lastIteration;
			assert(ll_diff >= 0);

			double ll_relative = ll_thisIteration/ll_lastIteration;

			std::cout << "\tImprovement: " << ll_diff << std::endl;
			std::cout << "\tRelative   : " << ll_relative << std::endl;
		
			if((ll_diff < 20) || (EMiteration > 30))
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
		
		f = f_nextIteration;
		EMiteration++;
		ll_lastIteration = ll_thisIteration;
	}

	std::string output_pot_frequencies = mappedFile + ".EM.WIMP";
	std::string output_recalibrated_mappings = mappedFile + ".EM";
	std::string output_assigned_reads_and_identities = mappedFile + ".EM.lengthAndIdentitiesPerMappingUnit";
	std::string output_reads_taxonID = mappedFile + ".EM.reads2Taxon";
	std::string output_contig_coverage = mappedFile + ".EM.contigCoverage";
	std::string output_evidence_unknownSpecies = mappedFile + ".EM.evidenceUnknownSpecies";

	std::ofstream strout_reads_identities(output_assigned_reads_and_identities);
	assert(strout_reads_identities.is_open());
	strout_reads_identities << "AnalysisLevel" << "\t" << "ID" << "\t" << "readI" << "\t" << "Identity" << "\t" << "Length" << "\n";

	std::ofstream strout_recalibrated_mappings(output_recalibrated_mappings);
	assert(strout_recalibrated_mappings.is_open());

	std::ofstream strout_reads_taxonIDs(output_reads_taxonID);
	assert(strout_reads_taxonIDs.is_open());
	
	size_t coverage_windowSize = 1000;
	std::map<std::string, std::map<std::string, std::vector<size_t>>> coverage_per_contigID;
	std::map<std::string, std::map<std::string, std::vector<size_t>>> coverage_per_contigID_countReads;
	std::map<std::string, std::map<std::string, size_t>> size_last_window;

	std::map<std::string, size_t> reads_per_taxonID;
	size_t runningReadI = 0;
	std::map<std::string, std::vector<double>> identities_per_taxonID;
	std::function<void(const std::vector<std::string>&)> processOneRead_final = [&](const std::vector<std::string>& readLines) -> void
	{
		assert(readLines.size() > 0);
		std::vector<oneMappingLocation> mappingLocations = getMappingLocations(taxonInfo, f, readLines);

		std::string readID;
		assert(readLines.size() == mappingLocations.size());
		for(unsigned int lineI = 0; lineI < readLines.size(); lineI++)
		{
			double finalMappingQuality = mappingLocations.at(lineI).p;
			std::vector<std::string> line_fields = split(readLines.at(lineI), " ");
			readID = line_fields.at(0);			
			line_fields.at(13) = std::to_string(finalMappingQuality);
			strout_recalibrated_mappings << join(line_fields, " ") << "\n";
		}

		oneMappingLocation bestMapping = getBestMapping(mappingLocations);

		strout_reads_identities << "EqualCoverageUnit" << "\t" << bestMapping.contigID << "\t" << runningReadI << "\t" << bestMapping.identity << "\t" << bestMapping.readLength << "\n";
		strout_reads_taxonIDs << readID << "\t" << bestMapping.taxonID << "\n";
		identities_per_taxonID[bestMapping.taxonID].push_back(bestMapping.identity);
		
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
				else
				{
					size_last_window[bestMapping.taxonID][bestMapping.contigID] = coverage_windowSize;
				}
			}
			coverage_per_contigID[bestMapping.taxonID][bestMapping.contigID].resize(n_windows, 0);
			coverage_per_contigID_countReads[bestMapping.taxonID][bestMapping.contigID].resize(n_windows, 0);
		}

		size_t bestMappingStopPos = (bestMapping.stop >= taxonInfo.at(bestMapping.taxonID).at(bestMapping.contigID)) ? (taxonInfo.at(bestMapping.taxonID).at(bestMapping.contigID) - 1) : bestMapping.stop;
		for(size_t contigPos = bestMapping.start; contigPos <= bestMappingStopPos; contigPos += coverage_windowSize)
		{
			size_t windowI = contigPos / coverage_windowSize;
			assert(windowI >= 0);
			assert(windowI < coverage_per_contigID.at(bestMapping.taxonID).at(bestMapping.contigID).size());
			assert(windowI < coverage_per_contigID_countReads.at(bestMapping.taxonID).at(bestMapping.contigID).size());

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
			coverage_per_contigID_countReads.at(bestMapping.taxonID).at(bestMapping.contigID).at(windowI)++;
		}

		runningReadI++;
	};
	
	std::cout << "Outputting mappings with adjusted alignment qualities." << std::endl;
	callBackForAllReads(mappedFile, processOneRead_final);
	strout_reads_identities.close();
	strout_reads_taxonIDs.close();

	cleanF(f, reads_per_taxonID, mappingStats.at("ReadsMapped"));
	producePotFile(output_pot_frequencies, T, f, reads_per_taxonID, nTotalReads, nUnmapped, nTooShort);

	std::ofstream strout_coverage(output_contig_coverage);
	strout_coverage << "taxonID" << "\t" << "equalCoverageUnitLabel" << "\t" << "contigID" << "\t" << "start" << "\t" << "stop" << "\t" << "nBases" << "\t" << "readCoverage" << "\n";

	std::map<std::string, std::string> contigID_2_taxon;
	for(auto taxonData : coverage_per_contigID)
	{
		for(auto contigData : taxonData.second)
		{
			for(size_t windowI = 0; windowI < contigData.second.size(); windowI++)
			{
				size_t windowLength = coverage_windowSize;
				if((windowI == (contigData.second.size()-1)))
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

				if(contigID_2_taxon.count(contigData.first))
				{
					assert(contigID_2_taxon.at(contigData.first) == taxonData.first);
				}
				else
				{
					contigID_2_taxon[contigData.first] = taxonData.first;
				}

			}
		}
	}

	// find taxon ID with highest median identity
	std::string taxonID_highestMedianIdentity;
	double highestMedianIdentity;
	double highestMedian_oneThird;
	double highestMedian_twoThirds;
	double highestMedian_oneThird_cumulativeP;
	// double highestMedian_twoThirds_cumulativeP;

	for(auto identitiesEntry : identities_per_taxonID)
	{
		std::string taxonID = identitiesEntry.first;
		std::vector<double> identities = identitiesEntry.second;
		if(identities.size() >= minimumReadsPerBestContig)
		{
			std::sort(identities.begin(), identities.end());
			double medianIdentity = identities.at(identities.size()/2);
			if((taxonID_highestMedianIdentity.length() == 0) || (medianIdentity > highestMedianIdentity))
			{
				highestMedianIdentity = medianIdentity;
				taxonID_highestMedianIdentity = taxonID;

				highestMedian_oneThird = identities.at(identities.size() * (1.0/3.0));
				highestMedian_twoThirds = identities.at(identities.size() * (2.0/3.0));
				assert(highestMedian_oneThird <= highestMedian_twoThirds);

				size_t highestMedian_oneThird_cumulativeN = 0;
				size_t highestMedian_twoThirds_cumulativeN = 0;

				for(auto idty : identities)
				{
					if(idty <= highestMedian_oneThird)
					{
						highestMedian_oneThird_cumulativeN++;
					}
					if(idty <= highestMedian_twoThirds)
					{
						highestMedian_twoThirds_cumulativeN++;
					}
				}

				assert(highestMedian_twoThirds_cumulativeN <= highestMedian_oneThird_cumulativeN);
				highestMedian_oneThird_cumulativeP = (double)highestMedian_oneThird_cumulativeN / (double)identities.size();
				// highestMedian_twoThirds_cumulativeP = (double)highestMedian_twoThirds_cumulativeN / (double)identities.size();
			}
		}
	}

	// todo
	size_t minimum_notManyNs_eitherSide = -1;
	std::map<std::string, std::vector<size_t>> Ns_per_genomeWindow = get_NS_per_window(DBdir, coverage_windowSize, coverage_per_contigID);
	std::map<std::string, std::vector<bool>> use_genomeWindow;
	std::map<std::string, size_t> genomeWindows;
	std::map<std::string, size_t> genomeWindows_usable;
	std::map<std::string, size_t> genomeWindows_usable_nReads;
	std::map<std::string, double> genomeWindows_usable_coverageZero;

	std::map<std::string, size_t> genomeWindows_coverageIsZero;
	std::map<std::string, double> genomeWindows_coverageIsZero_expected;
	std::map<std::string, double> genomeWindows_coverageIsZero_p;

	for(auto contigData : Ns_per_genomeWindow)
	{
		std::string taxonID = contigID_2_taxon.at(contigData.first);

		std::map<std::string, std::vector<size_t>> runningBases_notManyNs_forward;
		std::map<std::string, std::vector<size_t>> runningBases_notManyNs_backward;

		runningBases_notManyNs_forward[contigData.first].resize(contigData.second.size(), 0);
		runningBases_notManyNs_backward[contigData.first].resize(contigData.second.size(), 0);

		use_genomeWindow[contigData.first].resize(contigData.second.size(), false);

		size_t running_forward = 0;
		for(size_t windowI = 0; windowI < contigData.second.size(); windowI++)
		{
			runningBases_notManyNs_forward.at(contigData.first).at(windowI) = running_forward;
			size_t windowI_Ns = contigData.second.at(windowI);
			size_t windowI_length = coverage_windowSize;
			if((windowI == (contigData.second.size()-1)))
			{
				windowI_length = size_last_window.at(taxonID).at(contigData.first);
			}

			double propN = (double)windowI_Ns/(double)windowI_length;

			if(propN <= 0.02)
			{
				running_forward += windowI_length;
			}
			else
			{
				running_forward = 0;
			}
		}

		size_t running_backward = 0;
		for(long long windowI = (contigData.second.size() - 1); windowI >= 0; windowI--)
		{
			runningBases_notManyNs_backward.at(contigData.first).at(windowI) = running_backward;
			size_t windowI_Ns = contigData.second.at(windowI);
			size_t windowI_length = coverage_windowSize;
			if((windowI == ((long long)contigData.second.size()-1)))
			{
				windowI_length = size_last_window.at(taxonID).at(contigData.first);
			}

			double propN = (double)windowI_Ns/(double)windowI_length;

			if(propN <= 0.02)
			{
				running_backward += windowI_length;
			}
			else
			{
				running_backward = 0;
			}
		}

		size_t use_windows = 0;
		size_t use_windows_nReads = 0;
		size_t use_windows_zeroCoverage = 0;
		for(size_t windowI = 0; windowI < contigData.second.size(); windowI++)
		{
			if((runningBases_notManyNs_forward.at(contigData.first).at(windowI) >= minimum_notManyNs_eitherSide)
				&& (runningBases_notManyNs_backward.at(contigData.first).at(windowI) >= minimum_notManyNs_eitherSide)
			)
			{
				use_genomeWindow[contigData.first].at(windowI) = true;
				use_windows++;
				use_windows_nReads += coverage_per_contigID_countReads.at(taxonID).at(contigData.first).at(windowI);
				if(coverage_per_contigID_countReads.at(taxonID).at(contigData.first).at(windowI) == 0)
				{
					use_windows_zeroCoverage++;
				}
			}
		}

		if(genomeWindows_usable.count(taxonID) == 0)
		{
			genomeWindows[taxonID] = 0;
			genomeWindows_usable[taxonID] = 0;
			genomeWindows_usable_nReads[taxonID] = 0;
			genomeWindows_usable_coverageZero[taxonID] = 0;
		}

		genomeWindows.at(taxonID) += contigData.second.size();
		genomeWindows_usable.at(taxonID) += use_windows;
		genomeWindows_usable_nReads.at(taxonID) += use_windows_nReads;
		genomeWindows_usable_coverageZero.at(taxonID) += use_windows_zeroCoverage;

	}

	{
		double oneThird_p = highestMedian_oneThird_cumulativeP;
		std::ofstream evidenceNewSpeciesStream(output_evidence_unknownSpecies);
		evidenceNewSpeciesStream << "taxonID"
				<< "\t" << "nReads"
				<< "\t" << "propBottomThirdReadIdentities"
				<< "\t" << "expectedPropBottomThirdReadIdentities"
				<< "\t" << "pValue_BottomThirdReadIdentities"
				<< "\t" << "coverageWindows_totalGenome"
				<< "\t" << "coverageWindows_usable"
				<< "\t" << "coverageWindows_usable_averageCoverage"
				<< "\t" << "coverageWindows_usable_coverageIsZero"
				<< "\t" << "coverageWindows_usable_coverageIsZero_expected"
				<< "\t" << "coverageWindows_usable_coverageIsZero_P" << "\n";

		boost::math::chi_squared chiSq_oneDf(1);

		for(auto identitiesEntry : identities_per_taxonID)
		{
			std::string taxonID = identitiesEntry.first;

			// identities p-value
			std::string propBottomThirdReadIdentities_str = "NA";
			std::string pValue_identities_str = "NA";
			std::string oneThird_p_str = "NA";
			std::vector<double> identities = identitiesEntry.second;
			if(taxonID_highestMedianIdentity.length())
			{
				size_t observed_n_OneThird = 0;
				for(auto idty : identities)
				{
					if(idty <= highestMedian_oneThird)
					{
						observed_n_OneThird++;
					}
				}
				assert(identities.size());

				size_t observed_n_nonOneThird = identities.size() - observed_n_OneThird;

				double expected_n_oneThird =  oneThird_p * identities.size();
				double expected_n_nonOneThird = identities.size() - expected_n_oneThird;

				oneThird_p_str = std::to_string(oneThird_p);
				double testStatistic = pow(observed_n_OneThird - expected_n_oneThird, 2)/expected_n_oneThird
						+ pow(observed_n_nonOneThird - expected_n_nonOneThird, 2)/expected_n_nonOneThird;

				propBottomThirdReadIdentities_str = std::to_string(((double)observed_n_nonOneThird / (double)identities.size()));
				pValue_identities_str = std::to_string(boost::math::cdf(chiSq_oneDf, testStatistic));
			}

			// coverage p-value
			std::string usableWindows_averageCoverage_str = "NA";
			std::string usableWindows_coverageZero_expected_str = "NA";
			std::string usableWindows_coverageZero_p_str = "NA";

			if(genomeWindows_usable.at(taxonID) > 0)
			{
				double usableWindows_averageCoverage = (double)genomeWindows_usable_nReads.at(taxonID)/(double)genomeWindows_usable.at(taxonID);
				usableWindows_averageCoverage_str = std::to_string(usableWindows_averageCoverage);

				double windows_0_probability = boost::math::pdf(boost::math::poisson_distribution<>(usableWindows_averageCoverage), 0);

				usableWindows_coverageZero_expected_str = std::to_string(genomeWindows_usable.at(taxonID) * windows_0_probability);
				double windows_0_pValue = 1;
				if(genomeWindows_usable_coverageZero.at(taxonID) > 0)
				{
					double p_belowObservedValue = boost::math::cdf(boost::math::binomial_distribution<>(genomeWindows_usable.at(taxonID), windows_0_probability), genomeWindows_usable_coverageZero.at(taxonID) - 1);
					assert((p_belowObservedValue >= 0) && (p_belowObservedValue <= 1));
					windows_0_pValue = 1 - p_belowObservedValue;
				}

				usableWindows_coverageZero_p_str = std::to_string(windows_0_pValue);
			}

			evidenceNewSpeciesStream << taxonID
					<< "\t" << identities.size()
					<< "\t" << propBottomThirdReadIdentities_str
					<< "\t" << oneThird_p_str
					<< "\t" << pValue_identities_str
					<< "\t" << genomeWindows.at(taxonID)
					<< "\t" << genomeWindows_usable.at(taxonID)
					<< "\t" << usableWindows_averageCoverage_str
					<< "\t" << genomeWindows_usable_coverageZero.at(taxonID)
					<< "\t" << usableWindows_coverageZero_expected_str
					<< "\t" << usableWindows_coverageZero_p_str << "\n";
		}
	}
}

void cleanF(std::map<std::string, double>& f, const std::map<std::string, size_t>& reads_per_taxonID, size_t distributedReads)
{
	double minFreq = 0.9 * (1.0/(double)distributedReads);
	std::set<std::string> taxonID_to_delete;
	
	for(auto fI : f)
	{
		if((fI.second < minFreq) && (reads_per_taxonID.count(fI.first) == 0))
		{
			taxonID_to_delete.insert(fI.first);
		}
	}
	
	for(auto taxonID : taxonID_to_delete)
	{
		f.erase(taxonID);
	}
	
	double f_sum = 0;
	for(auto fI : f)
	{
		f_sum += fI.second;
	}
	assert(f_sum > 0);
	for(auto& fI : f)
	{
		fI.second /= f_sum;
	}
}
void callBackForAllReads(std::string mappedFile, std::function<void(const std::vector<std::string>&)>& callbackFunction)
{
	std::ifstream mappingsStream (mappedFile);
	assert(mappingsStream.is_open());

	std::string runningReadID;
	std::vector<std::string> runningReadLines;

	std::string line;
	while(mappingsStream.good())
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
			if(runningReadLines.size())
			{
				callbackFunction(runningReadLines);
			}
			
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
	size_t lineNumber = 0;
	while(mappingsStream.good())
	{
		std::getline(mappingsStream, line);
		eraseNL(line);
		lineNumber++;
		if(line.length())
		{
			std::vector<std::string> line_fields = split(line, " ");
			if(line_fields.size() < 6)
			{
				std::cerr << "File " << mappedFile << " has weird format - is this a mappings file generated by MetaMap? Offending line: " << lineNumber << std::endl;
				exit(1);
			}
			std::string mappingTarget = line_fields.at(5);
			std::string taxonID = extractTaxonId(mappingTarget);
			forReturn.insert(taxonID);
		}
	}
	mappingsStream.close();
	return forReturn;
}

boost::regex matchTaxonID{"kraken:taxid\\|(x?\\d+)"};
boost::smatch matchBrackets;
std::string extractTaxonId(std::string line)
{
	assert(boost::regex_search(line, matchBrackets, matchTaxonID));
	std::string taxonID = matchBrackets[1];
	return taxonID;
}

std::set<std::string> getRelevantLevelNames()
{
	return std::set<std::string>({"species", "genus", "family", "order", "phylum", "superkingdom"});
}

std::map<std::string, std::vector<size_t>> get_NS_per_window(std::string DBdir, size_t windowSize, const std::map<std::string, std::map<std::string, std::vector<size_t>>>& coverage_per_contigID)
{
	std::map<std::string, std::vector<size_t>> forReturn;
	std::string fn_windows = DBdir + "/contigNstats_windowSize_" + std::to_string(windowSize) + ".fa";
	std::ifstream windowStream;
	windowStream.open(fn_windows.c_str());
	assert(windowStream.is_open());
	std::string line;
	while(windowStream.good())
	{
		std::getline(windowStream, line);
		eraseNL(line);
		if(line.length() != 0)
		{
			std::vector<std::string> line_fields = split(line, " ");
			assert(line_fields.size() == 3);
			std::string taxonID = line_fields.at(0);
			std::string contigID = line_fields.at(1);
			if(coverage_per_contigID.count(taxonID) && coverage_per_contigID.at(taxonID).count(contigID))
			{
				std::vector<std::string> n_fields = split(line_fields.at(2), ";");
				assert(n_fields.size() == coverage_per_contigID.at(taxonID).at(contigID).size());
				std::vector<size_t> n_fields_int;
				for(auto e : n_fields)
				{
					n_fields_int.push_back(std::stoull(e));
				}
				assert(forReturn.count(contigID) == 0);
				forReturn[contigID] = n_fields_int;
			}
		}
	}

	for(auto taxonEntries : coverage_per_contigID)
	{
		for(auto contigEntries : taxonEntries.second)
		{
			assert(forReturn.count(contigEntries.first));
		}
	}

	return forReturn;
}

}

#endif /* META_FEM_H_ */
