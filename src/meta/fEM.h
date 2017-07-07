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
#include <string>
#include <fstream>
#include <assert.h>
#include <boost/regex.hpp>

#include "util.h"

std::string extractTaxonId(std::string line);
std::map<std::string, std::map<std::string, size_t>> loadRelevantTaxonInfo(std::string DBdir, const std::set<std::string>& taxonIDs);
std::set<std::string> getTaxonIDsFromMappingsFile(std::string mappedFile);
void callBackForAllReads(std::string mappedFile, std::function<void(const std::vector<std::string>&)>& callbackFunction);

void doEM(std::string DBdir, std::string mappedFile)
{
	std::set<std::string> relevantTaxonIDs = getTaxonIDsFromMappingsFile(mappedFile);

	// get genome info
	std::map<std::string, std::map<std::string, size_t>> taxonInfo = loadRelevantTaxonInfo(DBdir, relevantTaxonIDs);

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

		std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines ) -> void
		{
			assert(readLines.size() > 0);

			std::set<std::string> saw_contigIDs;
			std::map<std::string, std::vector<double>> mappingQualities_perTaxon;
			std::map<std::string, std::vector<double>> mappingProbabilities_perTaxon;

			std::string readID;
			size_t readLength = -1;
			double p_thisRead = 0;

			for(auto line : readLines)
			{
				std::vector<std::string> line_fields = split(line, " ");

				std::string contigID = line_fields.at(5);
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

				saw_contigIDs.insert(contigID);
				mappingQualities_perTaxon[contig_taxonID].push_back(mappingQuality);
			}

			double p_allLocations = 0;
			for(auto targetTaxon : mappingQualities_perTaxon)
			{
				std::string targetTaxonID = targetTaxon.first;
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

				assert(possibleMappingLocations > 0);
				for(auto mappingQuality : targetTaxon.second)
				{
					double p_thisMapping = f.at(targetTaxonID) * (1/(double)possibleMappingLocations) * mappingQuality;
					mappingProbabilities_perTaxon[targetTaxonID].push_back(p_thisMapping);
					p_allLocations += p_thisMapping;
					p_thisRead += p_thisMapping;
				}
			}

			assert(p_allLocations > 0);
			for(auto targetTaxon : mappingProbabilities_perTaxon)
			{
				for(auto& mappingP : targetTaxon.second)
				{
					mappingP /= p_allLocations;
					f_nextIteration.at(targetTaxon.first) += mappingP;
				}
			}

			assert(p_thisRead > 0);
			assert(p_thisRead < 1);
			ll_thisIteration += log(p_thisRead);
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


#endif /* META_FEM_H_ */
