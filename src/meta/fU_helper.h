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
	int minimum;
	std::map<int, double> h;

public:

	std::set<int> getKeys() const
	{
		std::set<int> forReturn;
		for(auto hI : h)
		{
			forReturn.insert(hI.first);
		}
		return forReturn;
	}

	int getMinimum() const
	{
		return minimum;
	}

	double getP(int idty) const
	{
		assert(idty <= 100);
		assert(idty >= 0);
		if(h.count(idty))
		{
			return h.at(idty);
		}
		else
		{
			return 0;
		}
	}

	void readFromEMOutput(std::string fn)
	{
		minimum = -1;
		std::ifstream f (fn);
		assert(f.is_open());

		std::string line;
		assert(f.good());
		std::getline(f, line);
		eraseNL(line);
		int minIdty = int(100 * std::stod(line) + 0.5);
		while(f.good())
		{
			std::getline(f, line);
			eraseNL(line);
			assert(line.length());
			std::vector<std::string> line_fields = split(line, "\t");
			assert(line_fields.size() == 2);
			int idty = int(100 * std::stod(line_fields.at(0)) + 0.5);
			double p = std::stod(line_fields.at(1));
			assert(idty >= 0);
			assert(idty <= 100);
			assert(p >= 0);
			assert(p <= 1);
			h[idty] = p;
			if((minimum == -1) || (idty < minimum))
			{
				minimum = idty;
			}
		}

		assert(minimum == minIdty);

		int notDefined = 0;
		for(int i = minimum; i < 100; i++)
		{
			if(h.count(i) == 0)
			{
				notDefined++;
			}
		}
		double newP_notDefined = 1.0/(double)notDefined;
		double newSum = 0;
		for(int i = minimum; i < 100; i++)
		{
			if(h.count(i) == 0)
			{
				h[i] = newP_notDefined;
			}
			newSum += h.at(i);
		}

		for(auto& p : h)
		{
			p.second /= newSum;
		}
	}
};

class treeAdjustedIdentities
{
public:
	std::map<std::string, std::map<size_t, std::map<int, double>>> D;

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
	identityHistogram& iH;
	treeAdjustedIdentities& tAI;
	std::map<size_t, double> readLengthHistogram;

	std::map<std::string, std::map<int, double>> indirectMapping_cache;

public:
	identityManager(identityHistogram& iH, treeAdjustedIdentities& tAI, std::map<size_t, double> readLengthHistogram) : iH(iH), tAI(tAI)
	{

	}

	double getP(int identity, std::string taxonID, bool directlyAttached)
	{
		if(directlyAttached)
		{
			double p = iH.getP(identity);
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
			if(indirectMapping_cache[taxonID].count(identity))
			{
				return indirectMapping_cache.at(taxonID).at(identity);
			}
			else
			{
				std::map<int, double> histogramForNode = getShiftedIdentityHistogramForNode(taxonID);
				return (histogramForNode.count(identity)) ? histogramForNode.at(identity) : 0;
			}
		}
	}

	std::map<int, double> getShiftedIdentityHistogramForNode(std::string taxonID)
	{
		assert(tAI.D.count(taxonID));
		std::map<int, double> forReturn;

		double forReturn_sum = 0;
		for(auto rLI : tAI.D.at(taxonID))
		{
			size_t readLength = rLI.first;
			double readLength_P = getReadLengthP(readLength);

			std::set<int> keys_iH = iH.getKeys();

			for(auto k1 : keys_iH)
			{
				double p1 = iH.getP(k1);
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

					if(newK_int < iH.getMinimum())
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

	static std::map<int, double> getConvolutedHistogram(const identityHistogram& iH, std::map<int, double> additionalHistogram)
	{
		std::map<int, double> forReturn;
		std::set<int> keys_iH = iH.getKeys();
		for(auto k1 : keys_iH)
		{
			double p1 = iH.getP(k1);
			for(auto k2P: additionalHistogram)
			{
				int k2 = k2P.first;
				double p2 = k2P.first;

				double newK = ((double)k1/100.0) * ((double)k2/100.0);
				assert(newK >= 0);
				assert(newK <= 1);
				int newK_int = int((newK*100) + 0.5);
				double newP = p1 * p2;

				if(newK_int < iH.getMinimum())
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
					double readLength_diff = rL_ip1 - rL_i;
					assert(readLength_diff > 0);
					double weight_right = (readLength - rL_i)/readLength_diff;
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


}

#endif /* META_FU_HELPER_H_ */
