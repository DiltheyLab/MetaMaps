/*
 * fU_helper.h
 *
 *  Created on: Jul 10, 2017
 *      Author: diltheyat
 */

#ifndef META_FU_HELPER_H_
#define META_FU_HELPER_H_

#include <string>
#include <assert.h>
#include <vector>
#include <map>
#include <fstream>

#include "fEM.h"
#include "util.h"

namespace meta
{


class identityAndReadLengthHistogram
{
protected:
	int minimumIdentity;
	std::map<int, double> identityHistogram;
	std::map<size_t, double> readLengthHistogram;

public:

	std::set<int> getIdentityKeys() const
	{
		std::set<int> forReturn;
		for(auto hI : identityHistogram)
		{
			forReturn.insert(hI.first);
		}
		return forReturn;
	}

	int getIdentityMinimum() const
	{
		return minimumIdentity;
	}

	double getIdentityP(int idty) const
	{
		assert(idty <= 100);
		assert(idty >= 0);
		if(identityHistogram.count(idty))
		{
			return identityHistogram.at(idty);
		}
		else
		{
			return 0;
		}
	}

	void readFromEMOutput(std::string fn, size_t minimumReadsPerContig)
	{
		minimumIdentity = -1;
		std::ifstream f (fn);
		assert(f.is_open());

		std::map<std::string, std::vector<double>> identitiesPerUnit;
		std::map<std::string, std::vector<size_t>> lengthPerUnit;

		std::string line;
		while(f.good())
		{
			std::getline(f, line);
			eraseNL(line);
			assert(line.length());
			std::vector<std::string> line_fields = split(line, "\t");

			std::string contigID = line_fields.at(1);
			double identity = std::stod(line_fields.at(3));
			size_t length = std::stoull(line_fields.at(4));

			identitiesPerUnit[contigID].push_back(identity);
			lengthPerUnit[contigID].push_back(length);
		}

		double highestMedian;
		std::string highestMedianContigID;
		for(auto contigData : identitiesPerUnit)
		{
			if(contigData.second.size() > minimumReadsPerContig)
			{
				std::vector<double> identities = contigData.second;
				std::sort(identities.begin(), identities.end());
				double medianIdentity = identities.at(identities.size()/2);
				if((highestMedianContigID.length() == 0) || (medianIdentity > highestMedian))
				{
					highestMedian = medianIdentity;
					highestMedianContigID = contigData.first;
				}
			}
		}

		if(highestMedianContigID.length() == 0)
		{
			std::cerr << "\n\nERROR: Cannot fit read length and identity distribution from file " << fn << "\n";
			std::cerr << "\nThe most likely explanation is that no contig in the fitted data has more than " << minimumReadsPerContig << " (this is a parameter that can be changed) reads assigned to it - abort.\n\n";
			std::cerr << std::flush;
			exit(1);
		}

		std::vector<double> identities_bestContig = identitiesPerUnit.at(highestMedianContigID);
		std::vector<size_t> lengths_bestContig = lengthPerUnit.at(highestMedianContigID);

		std::map<int, size_t> identityHistogram_int;
		for(auto i : identities_bestContig)
		{
			int iI = int((i*100) + 0.5);
			if(identityHistogram_int.count(iI) == 0)
			{
				identityHistogram_int[iI] = 0;
			}
			assert(iI >= 0);
			assert(iI <= 100);
			identityHistogram_int.at(iI)++;
			if((minimumIdentity == -1) || (iI < minimumIdentity))
			{
				minimumIdentity = iI;
			}
		}

		identityHistogram.clear();
		for(auto i_and_p : identityHistogram_int)
		{
			identityHistogram[i_and_p.first] = (double)i_and_p.second/(double)identities_bestContig.size();
		}

		std::map<size_t, size_t> readLengthHistogram_int;
		for(auto l : lengths_bestContig)
		{
			size_t l_1000 = 1000 * int((l / 1000) + 0.5);
			if(readLengthHistogram_int.count(l_1000) == 0)
			{
				readLengthHistogram_int[l_1000] = 0;
			}
			readLengthHistogram_int.at(l_1000)++;
		}

		readLengthHistogram.clear();
		for(auto lI : readLengthHistogram_int)
		{
			readLengthHistogram[lI.first] = (double)lI.second / lengths_bestContig.size();
		}


		int notDefined = 0;
		for(int i = minimumIdentity; i < 100; i++)
		{
			if(identityHistogram.count(i) == 0)
			{
				notDefined++;
			}
		}

		double newP_notDefined = 1.0/(double)notDefined;
		double newSum = 0;
		for(int i = minimumIdentity; i < 100; i++)
		{
			if(identityHistogram.count(i) == 0)
			{
				identityHistogram[i] = newP_notDefined;
			}
			newSum += identityHistogram.at(i);
		}
		for(auto& p : identityHistogram)
		{
			p.second /= newSum;
		}
	}

	double getReadLengthP(size_t readLength)
	{
		std::vector<size_t> readLengths_vec;
		readLengths_vec.reserve(readLengthHistogram.size());
		for(auto rL : readLengthHistogram)
		{
			readLengths_vec.push_back(rL.first);
		}
		std::sort(readLengths_vec.begin(), readLengths_vec.end());
		if(readLengths_vec.size() > 1)
		{
			assert(readLengths_vec.at(0) < readLengths_vec.at(1));
		}
		if(readLength < readLengths_vec.at(0))
		{
			return readLengthHistogram.at(readLengths_vec.at(0));
		}
		else if(readLength >= readLengths_vec.at(readLengths_vec.size()-1))
		{
			return readLengthHistogram.at(readLengths_vec.at(readLengths_vec.size()-1));
		}
		else
		{
			bool found = false;
			for(size_t i = 0; i < (readLengths_vec.size() - 1); i++)
			{
				size_t rL_i = readLengths_vec.at(i);
				size_t rL_ip1 = readLengths_vec.at(i+1);
				if((readLength >= rL_i) && (readLength < rL_ip1))
				{
					found = true;
					long long readLength_diff = rL_ip1 - rL_i;
					assert(readLength_diff > 0);
					double weight_right = (readLength - rL_i)/(double)readLength_diff;
					assert((weight_right >= 0) && (weight_right <= 1));
					double weight_left = 1 - weight_right;
					double combinedP = readLengthHistogram.at(rL_i) * weight_left + readLengthHistogram.at(rL_ip1) * weight_right;
					return combinedP;
				}
			}
			assert(found);
			return 0;
		}
	}
};

class treeAdjustedIdentities
{
public:
	std::map<std::string, std::map<size_t, std::map<int, double>>> D;

	std::vector<size_t> getTwoClosestReadLenghts(std::string taxonID, size_t targetReadLength)
	{
		std::vector<size_t> readLengths_vec;
		readLengths_vec.reserve(D.at(taxonID).size());
		for(auto rL : D.at(taxonID))
		{
			readLengths_vec.push_back(rL.first);
		}
		std::sort(readLengths_vec.begin(), readLengths_vec.end());
		assert(readLengths_vec.size());

		if(targetReadLength < readLengths_vec.at(0))
		{
			return std::vector<size_t>({readLengths_vec.at(0)});
		}
		else if(targetReadLength >= readLengths_vec.at(readLengths_vec.size()-1))
		{
			return std::vector<size_t>({readLengths_vec.at(readLengths_vec.size()-1)});
		}
		else
		{
			bool found = false;
			for(size_t i = 0; i < (readLengths_vec.size() - 1); i++)
			{
				size_t rL_thisI = readLengths_vec.at(i);
				size_t rl_nextI = readLengths_vec.at(i+1);
				if((targetReadLength >= rL_thisI) && (targetReadLength < rl_nextI))
				{
					return std::vector<size_t>({rL_thisI, rl_nextI});
					found = true;
				}
			}
			assert(found);
			return std::vector<size_t>();
		}
	}

	bool nodeForIndirectAttachment(std::string taxonID)
	{
		return (D.count(taxonID) > 0);
	}

	void readFromFile(std::string fn, const std::set<std::string>& mappings_taxonIDs, const taxonomy& T)
	{
		std::set<std::string> relevantTaxonIDs;
		for(auto tI : mappings_taxonIDs)
		{
			relevantTaxonIDs.insert(tI);
			std::vector<std::string> upwardNodes = T.getUpwardNodes(tI);
			for(auto tI2 : upwardNodes)
			{
				relevantTaxonIDs.insert(tI2);
			}
		}

		std::ifstream f (fn);
		assert(f.is_open());

		std::string line;
		while(f.good())
		{
			std::getline(f, line);
			eraseNL(line);
			assert(line.length());
			std::vector<std::string> line_fields = split(line, "\t");

			std::string nodeID = line_fields.at(0);
			size_t readLength = std::stoull(line_fields.at(1));
			int idty = int(100 * std::stod(line_fields.at(2)));
			double p = std::stod(line_fields.at(3));

			if(relevantTaxonIDs.count(nodeID))
			{
				D[nodeID][readLength][idty] = p;
			}
		}
	}
};

std::map<size_t, double> readReadLengthHistogram(std::string fn)
{
	std::map<size_t, double> forReturn;
	std::ifstream f (fn);
	assert(f.is_open());

	double p_sum = 0;
	std::string line;
	while(f.good())
	{
		std::getline(f, line);
		eraseNL(line);
		assert(line.length());
		std::vector<std::string> line_fields = split(line, "\t");
		assert(line_fields.size() == 0);
		size_t readLength = std::stoull(line_fields.at(0));
		double p = std::stod(line_fields.at(1));
		assert(forReturn.count(readLength) == 0);
		forReturn[readLength] = p;
		p_sum += p;
	}
	assert(abs(1 - p_sum) <= 1e-4);
	return forReturn;
}

class identityManager
{
protected:
	identityAndReadLengthHistogram& iH;
	treeAdjustedIdentities& tAI;

	std::map<std::string, std::map<int, double>> indirectMapping_cache;
	std::map<std::string, std::map<size_t, std::map<int, double>>> indirectMapping_cache_granularReads;

	bool granularReadIdentities;

public:
	identityManager(identityAndReadLengthHistogram& iH, treeAdjustedIdentities& tAI) : iH(iH), tAI(tAI)
	{
		granularReadIdentities = false;
	}

	double getIdentityP(int identity, std::string taxonID, size_t readLength, bool directlyAttached)
	{
		if(directlyAttached)
		{
			double p = iH.getIdentityP(identity);
			if(p == 0)
			{
				return 1e-4;
			}
			else
			{
				return p;
			}
		}
		else
		{
			assert(tAI.D.count(taxonID));
			if(granularReadIdentities)
			{
				// todo perhaps reconsider
				readLength = size_t(readLength/(double)100 + 0.5) * 100;

				if(indirectMapping_cache_granularReads[taxonID][readLength].count(identity))
				{
					return indirectMapping_cache_granularReads.at(taxonID).at(readLength).at(identity);
				}
				else
				{
					std::map<int, double> histogramForNode = getShiftedIdentityHistogramForNode_oneReadLength(taxonID, readLength);
					double forReturn =(histogramForNode.count(identity)) ? histogramForNode.at(identity) : 0;
					indirectMapping_cache_granularReads.at(taxonID).at(readLength)[identity] = forReturn;
					return forReturn;
				}
			}
			else
			{
				if(indirectMapping_cache[taxonID].count(identity))
				{
					return indirectMapping_cache.at(taxonID).at(identity);
				}
				else
				{
					std::map<int, double> histogramForNode = getShiftedIdentityHistogramForNode(taxonID);
					double forReturn = (histogramForNode.count(identity)) ? histogramForNode.at(identity) : 0;
					indirectMapping_cache.at(taxonID)[identity] = forReturn;
					return forReturn;
				}
			}
		}
	}

	std::map<int, double> getShiftedIdentityHistogramForNode_oneReadLength(std::string taxonID, size_t readLength) const
	{
		assert(tAI.D.count(taxonID));

		std::vector<size_t> closestReadLengths = tAI.getTwoClosestReadLenghts(taxonID, readLength);
		assert((closestReadLengths.size() == 1) || (closestReadLengths.size() == 2));
		if(closestReadLengths.size() == 1)
		{
			return getConvolutedHistogram(iH, tAI.D.at(taxonID).at(closestReadLengths.at(0)));
		}
		else
		{
			size_t existing_rL1 = closestReadLengths.at(0);
			size_t existing_rL2  = closestReadLengths.at(1);
			assert(existing_rL1 < existing_rL2);

			long long readLength_diff = existing_rL2 - existing_rL1;
			assert(readLength_diff > 0);
			double weight_right = (readLength - existing_rL1)/(double)readLength_diff;
			assert((weight_right >= 0) && (weight_right <= 1));
			double weight_left = 1 - weight_right;

			std::map<int, double> h1 = getConvolutedHistogram(iH, tAI.D.at(taxonID).at(existing_rL1));
			std::map<int, double> h2 = getConvolutedHistogram(iH, tAI.D.at(taxonID).at(existing_rL2));

			std::set<int> combinedIdentities;
			for(auto i : h1)
			{
				combinedIdentities.insert(i.first);
			}
			for(auto i : h2)
			{
				combinedIdentities.insert(i.first);
			}
			std::map<int, double> forReturn;

			for(auto iV : combinedIdentities)
			{
				double v1 = (h1.count(iV)) ? h1.at(iV) : 0;
				double v2 = (h2.count(iV)) ? h2.at(iV) : 0;
				double vAvg = weight_left * v1 + weight_right * v2;
				forReturn[iV] = vAvg;
			}

			double fR_sum = 0;
			for(auto fR : forReturn)
			{
				fR_sum += fR.second;
			}
			assert(abs(1 - fR_sum) <= 1e-3);

			return forReturn;
		}
	}


	std::map<int, double> getShiftedIdentityHistogramForNode(std::string taxonID) const
	{
		assert(tAI.D.count(taxonID));
		std::map<int, double> forReturn;

		double forReturn_sum = 0;
		for(auto rLI : tAI.D.at(taxonID))
		{
			size_t readLength = rLI.first;
			double readLength_P = iH.getReadLengthP(readLength);

			std::set<int> keys_iH = iH.getIdentityKeys();

			for(auto k1 : keys_iH)
			{
				double p1 = iH.getIdentityP(k1);
				for(auto k2P: rLI.second)
				{
					int k2 = k2P.first;
					double p2 = k2P.first;

					double newK = ((double)k1/100.0) * ((double)k2/100.0);
					assert(newK >= 0);
					assert(newK <= 1);
					int newK_int = int((newK*100) + 0.5);

					double newP = readLength_P * p1 * p2;
					assert(newP > 0);
					assert(newP <= 1);

					if(newK_int < iH.getIdentityMinimum())
					{
						newK_int = 0;
					}

					if(forReturn.count(newK_int) == 0)
					{
						forReturn[newK_int] = 0;
					}
					forReturn.at(newK_int) += newP;
					forReturn_sum += newP;
				}
			}

		}

		assert(forReturn_sum > 0);
		for(auto& fR : forReturn)
		{
			fR.second /= forReturn_sum;
		}

		return forReturn;
	}

	static std::map<int, double> getConvolutedHistogram(const identityAndReadLengthHistogram& iH, std::map<int, double> additionalHistogram)
	{
		std::map<int, double> forReturn;
		std::set<int> keys_iH = iH.getIdentityKeys();
		for(auto k1 : keys_iH)
		{
			double p1 = iH.getIdentityP(k1);
			for(auto k2P: additionalHistogram)
			{
				int k2 = k2P.first;
				double p2 = k2P.first;

				double newK = ((double)k1/100.0) * ((double)k2/100.0);
				assert(newK >= 0);
				assert(newK <= 1);
				int newK_int = int((newK*100) + 0.5);
				double newP = p1 * p2;

				if(newK_int < iH.getIdentityMinimum())
				{
					newK_int = 0;
				}

				if(forReturn.count(newK_int) == 0)
				{
					forReturn[newK_int] = 0;
				}
				forReturn.at(newK_int) += newP;
			}
		}
		return forReturn;
	}

};


}

#endif /* META_FU_HELPER_H_ */
