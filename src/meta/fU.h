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


class oneMappingLocation_U
{
public:
	std::string readID;
	std::string taxonID;
	double originalIdentity;
	size_t minimizerUnion;
	size_t minimizerIntersection;
	long long readLength;
	double p;
	double mapQ;
	bool direct;
};

std::vector<oneMappingLocation_U> getMappingLocations_U(const identityManager& iM, const std::map<std::string, std::vector<std::string>>& indirectUpwardNodes, const std::map<std::string, int> indirectUpwardNodes_nSourceGenomes, const std::vector<std::string>& readLines)
{
	assert(readLines.size() > 0);

	bool verbose = false;

	std::string readID;
	long long readLength = -1;

	for(auto line : readLines)
	{
		std::vector<std::string> line_fields = split(line, " ");

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
	assert(readLength != -1);

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
		// size_t contigID_start = std::stoull(line_fields.at(7));
		// size_t contigID_stop = std::stoull(line_fields.at(8));

		std::string contig_taxonID = extractTaxonId(contigID);

		double identity = (std::stod(line_fields.at(9)))/100.0;
		int intersection_size = std::stoi(line_fields.at(10));
		int sketch_size = std::stoi(line_fields.at(11));
		assert(intersection_size <= sketch_size);

		assert(identity >= 0);
		assert(identity <= 1);

		{
			oneMappingLocation_U l;
			l.readID = readID;
			l.taxonID = contig_taxonID;
			l.p = 0;
			l.mapQ = 0;
			l.originalIdentity = identity;
			l.minimizerUnion = sketch_size;
			l.minimizerIntersection = intersection_size;
			l.readLength = readLength;
			l.direct = true;

			if((bestDirectMappings.count(contig_taxonID) == 0) || (bestDirectMappings.at(contig_taxonID).originalIdentity < l.originalIdentity))
			{
				bestDirectMappings[contig_taxonID] = l;
			}
		}

		for(auto indirectTaxon : indirectUpwardNodes.at(contig_taxonID))
		{
			oneMappingLocation_U lI;
			lI.readID = readID;
			lI.taxonID = indirectTaxon;
			lI.p = 0;
			lI.mapQ = 0;
			lI.readLength = readLength;
			lI.originalIdentity = identity;
			lI.minimizerUnion = sketch_size;
			lI.minimizerIntersection = intersection_size;
			lI.direct = false;

			if((bestInDirectMappings.count(indirectTaxon) == 0) || (bestInDirectMappings.at(indirectTaxon).originalIdentity < lI.originalIdentity))
			{
				bestInDirectMappings[indirectTaxon] = lI;
			}

			//if(verbose || (indirectTaxon == "1063"))
			//if(1 == 0)
			//{
			//	std::cout << "\t\tDirect to " << contig_taxonID << ": " << f.first.at(contig_taxonID) << " x " << iM.getIdentityP(identityInt, contig_taxonID, readLength, true) << "\n";
			//	std::cout << "\t\tIndirect to " << indirectTaxon << ": " << f.second.at(indirectTaxon) << " x " << iM.getIdentityP(identityInt, indirectTaxon, readLength, false) << "\n";
			//}
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

	return mappingLocations;
}



void compute_U_mappingQualities(std::vector<oneMappingLocation_U>& mappingLocations, const identityManager& iM, int kmerSize)
{
	bool verbose = false;
	assert(mappingLocations.size());

	for(auto& mL : mappingLocations)
	{
		mL.mapQ = 0;
	}
	
	int max_int_identity = -1;
	for(auto& mL : mappingLocations)
	{
		if(mL.direct)
		{
			int this_mL_int_identity = ceil(mL.originalIdentity * 100);
			if((max_int_identity == -1) || (this_mL_int_identity > max_int_identity))
			{
				max_int_identity = this_mL_int_identity;
			}
		}
		mL.p = 0;
	}
	
	std::string readID = mappingLocations.at(0).readID;
	assert((max_int_identity > 0) && (max_int_identity <= 100)); 
	
	if(readID == "m131217_210252_42161_c100527102550000001823083509281480_s1_p0/4505/6895_13844")
	{
		verbose = false;
	}

	if(verbose)
	{
		std::cerr << "Read " << mappingLocations.at(0).readID << ", maximum direct mapping identity (int) " << max_int_identity << "\n" << std::flush;
	}

	int iM_maxReadIdentity = iM.getMaximumReadIdentity();
	if(!(max_int_identity <= iM_maxReadIdentity)) 
	{
		std::cerr << "max_int_identity: " << max_int_identity << "\n";
		std::cerr << "iM_maxReadIdentity: " << iM_maxReadIdentity << "\n";
		std::cerr << std::flush;
	}
	assert(max_int_identity <= iM_maxReadIdentity);
	for(int readIdentity = max_int_identity; readIdentity <= iM_maxReadIdentity; readIdentity++)
	{
		double readIdentity_P = iM.getReadIdentityP(readIdentity);
		assert((readIdentity_P > 0) && (readIdentity_P <= 1));

		for(auto& mL : mappingLocations)
		{
			int n_kmers = mL.readLength - kmerSize + 1;

			double thisReadIdentity_thisML_l = 0;
			if(mL.direct)
			{
				/*
				std::cerr << "kmerSize" << ": " << kmerSize << "\n";
				std::cerr << "n_kmers" << ": " << n_kmers << "\n";
				std::cerr << "readIdentity" << ": " << readIdentity << "\n";
				std::cerr << "mL.minimizerUnion" << ": " << mL.minimizerUnion << "\n";
				std::cerr << "mL.minimizerIntersection" << ": " << mL.minimizerIntersection << "\n" << std::flush;
				*/
				
				double directL = mapWrap::likelihood_observed_set_sizes(kmerSize, n_kmers, readIdentity/100.0, mL.minimizerUnion, mL.minimizerIntersection);
				thisReadIdentity_thisML_l += directL;
				
				if(verbose)
				{
					std::cerr << "Direct " << readIdentity << "\n";
					std::cerr << "\t" << "directL (observed " << mL.originalIdentity << " v/s assumed " << readIdentity << ") " << ": " << directL << "\n";
					std::cerr << "\t" << "readIdentity_P" << ": " << readIdentity_P << "\n";
					std::cerr << "\t" << "(readIdentity_P * directL)" << ": " << readIdentity_P * directL << "\n" << std::flush;
				}
			}
			else  
			{
				if(verbose)
				{
					std::cerr << "Indirect " << readIdentity << "\n";
				}
				double indirect_acrossDistribution = 0;
				std::map<int, double> node_shiftDistribution = iM.getOriginalUHistogramForNode_oneReadLength(mL.taxonID, mL.readLength);
				double p_sum_generateMapping = 0;
				for(int count_p_generateMapping = 1; count_p_generateMapping >= 0; count_p_generateMapping--)
				{
					for(auto shiftE : node_shiftDistribution) 
					{
							
						if(verbose)
						{
							// std::cerr << "\t\t" << "shift" << " " << shiftE.first << " " << shiftE.second << "\n";
						}	
						
						if(shiftE.first == 0)  
						{
							continue;
						}
						int thisShift_identityLoss = 100 - shiftE.first;
						assert(thisShift_identityLoss >= 0);
						
						// int thisShift_identity = readIdentity - thisShift_identityLoss;
						double thisShift_identity = (readIdentity/100.0) * (shiftE.first/100.0);
						if(!(thisShift_identity > 0)) 
						{
							//std::cerr << "thisShift_identity" << ": " << thisShift_identity << "\n";
							//std::cerr << "readIdentity" << ": " << readIdentity << "\n";
							//std::cerr << "shiftE.first" << ": " << shiftE.first << "\n" << std::flush;		
							for(auto sD : node_shiftDistribution)
							{
								// std::cerr << "\t" << sD.first << " " << sD.second << "\n" << std::flush;
							}
						}
						assert(thisShift_identity > 0);
						
						if(verbose)
						{
							// std::cerr << "\t\t\teffective identity: " << thisShift_identity << "\n";
						} 
						
						if(thisShift_identity > (iM.getMinimumReadIdentity()/100.0))
						{	
							/*
							std::cerr << "kmerSize" << ": " << kmerSize << "\n";
							std::cerr << "n_kmers" << ": " << n_kmers << "\n";
							std::cerr << "thisShift_identity" << ": " << thisShift_identity << "\n";
							std::cerr << "mL.minimizerUnion" << ": " << mL.minimizerUnion << "\n";
							std::cerr << "mL.minimizerIntersection" << ": " << mL.minimizerIntersection << "\n" << std::flush;
							*/
							
							double indirectL = mapWrap::likelihood_observed_set_sizes(kmerSize, n_kmers, thisShift_identity, mL.minimizerUnion, mL.minimizerIntersection);
							
							if(count_p_generateMapping)
							{
								p_sum_generateMapping += (shiftE.second * indirectL);
							}
							else
							{
								assert(p_sum_generateMapping > 0);
								indirect_acrossDistribution += ((shiftE.second * indirectL)/p_sum_generateMapping);
							}
							if(verbose)
							{						
								//std::cerr << "\t\t\t\t" << "indirectL (observed " << mL.originalIdentity << " v/s assumed " << thisShift_identity << ") " << indirectL << "\n";
								//std::cerr << "\t\t\t\t" << "readIdentity_P: " << readIdentity_P << "\n";
								//std::cerr << "\t\t\t\t" << "(readIdentity_P * shiftE.second * indirectL): " << readIdentity_P * shiftE.second * indirectL << "\n" << std::flush;
							}
						}
					}
				}
				thisReadIdentity_thisML_l += indirect_acrossDistribution;
				
				if(verbose)
				{
					std::cerr << "\t" << "indirect_acrossDistribution" << ": " << indirect_acrossDistribution << "\n" << std::flush;
				}				
			}
			// assert(thisReadIdentity_thisML_l > 0); // perhaps too strong
			mL.mapQ += (readIdentity_P * thisReadIdentity_thisML_l); 
		}
	}
	
	double mL_p_sum = 0;
	for(auto mL : mappingLocations)
	{
		mL_p_sum += mL.mapQ;
	}	
	assert(mL_p_sum > 0);
	
	double mL_after = 0;
	double mL_after_direct = 0;
	double mL_after_indirect = 0;
	for(auto& mL : mappingLocations)
	{
		mL.mapQ /= mL_p_sum;
		mL_after += mL.mapQ;
		if(mL.direct)
		{
			mL_after_direct += mL.mapQ;
		}
		else
		{
			mL_after_indirect += mL.mapQ;
		}
	}	
	
	assert(abs(1 - mL_after) <= 1e-3);
	
	if(verbose)
	{
		std::cerr << "p_direct: " << mL_after_direct << "\n";
		std::cerr << "p_indirect: " << mL_after_indirect << "\n";
		std::cerr << std::flush;
	}
	if(verbose)
	{
		assert(1 == 0);
	}
	/*
	std::cerr << "newRead II\n";
	for(auto& mL : mappingLocations)
	{
		std::cerr << mL.p << "\n" << std::flush;
	}
	*/	
	
}

void generateUnknownMapQFile(std::string DBdir, std::string mappedFile, const identityManager& iM, const taxonomy& T, const skch::Parameters& parameters)
{
	std::map<std::string, size_t> mappingStats = getMappingStats(mappedFile);
	std::map<std::string, std::string> mappingParameters = getMappingParameters(mappedFile);

	std::set<std::string> taxonIDsInMappings = getTaxonIDsFromMappingsFile(mappedFile);

	std::set<std::string> mappableTaxonIDs = getDirectlyMappableTaxonIDs(DBdir);

	// std::pair<int, int> identity_minmax_inAlignments = getMinMaxIdentities(mappedFile);

	int kMerSize = std::stoi(mappingParameters.at("kmerSize"));

	std::set<std::string> relevantTaxonIDs_direct = taxonIDsInMappings;
	std::set<std::string> relevantTaxonIDs_indirect;

	std::map<std::string, std::vector<std::string>> indirectUpwardNodes;
	std::map<std::string, int> indirectUpwardNodes_nSourceGenomes;
	for(auto tI : taxonIDsInMappings)
	{
		std::vector<std::string> upwardTaxonIDs = T.getUpwardNodes(tI);
		indirectUpwardNodes[tI] = std::vector<std::string>();
		for(auto uTI : upwardTaxonIDs)
		{
			if(iM.tAI.nodeForIndirectAttachment(uTI))
			{
				relevantTaxonIDs_indirect.insert(uTI);
				indirectUpwardNodes[tI].push_back(uTI);
				if(indirectUpwardNodes_nSourceGenomes.count(uTI) == 0)
				{
					indirectUpwardNodes_nSourceGenomes[uTI] = iM.tAI.indirectAttachmentNode_sources(uTI);
				}
			}
		}
	}

	std::string fn_mappings_unknownMapQ = mappedFile + ".mapQ_U";
	std::ofstream outputStr_mappings_unknownMapQ;
	outputStr_mappings_unknownMapQ.open(fn_mappings_unknownMapQ.c_str());
	assert(outputStr_mappings_unknownMapQ.is_open());

	size_t processedRead = 0;

	std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines) -> void
	{
		processedRead++;
		if((processedRead % 10000) == 0)
		{
			std::cout << "\r U mapping quality calculation read << " << processedRead << " / " << mappingStats.at("ReadsMapped") << "   " << std::flush;
		}

		std::vector<oneMappingLocation_U> mappingLocations = getMappingLocations_U(iM, indirectUpwardNodes, indirectUpwardNodes_nSourceGenomes, readLines);
		
		compute_U_mappingQualities(mappingLocations, iM, kMerSize);
		
		double sum_mapQ = 0;
		for(auto location : mappingLocations)
		{
			outputStr_mappings_unknownMapQ << location.readID << " " << location.taxonID << " " << location.direct << " " << location.mapQ << " " << location.originalIdentity << "\n";
			sum_mapQ += location.mapQ;
		}
		assert(abs(1 - sum_mapQ) <= 1e-3);
	};

	callBackForAllReads(mappedFile, processOneRead, parameters);

	outputStr_mappings_unknownMapQ.close();

	std::cerr << "Generated " << fn_mappings_unknownMapQ << "\n" << std::flush;

	/*

	if(lines.size() > 0)
	{
		std::set<std::string> readIDs;
		std::set<int> read_lengths;

		std::vector<unsigned int> fields_per_line;
		std::vector<std::pair<unsigned int, unsigned int>> observed_set_sizes;
		std::vector<double> identities;

		double max_identity = -1;
		for(const std::string& line : lines)
		{
			std::vector<std::string> fields = split(line, " ");

			assert((fields.size() == 14) || (fields.size() == 15) || (fields.size() == 12));
			fields_per_line.push_back(fields.size());

			readIDs.insert(fields.at(0));
			read_lengths.insert(std::stoi(fields.at(1)));

			double identity = (std::stod(fields.at(9)))/100.0;
			int intersection_size = std::stoi(fields.at(10));
			int sketch_size = std::stoi(fields.at(11));
			assert(intersection_size <= sketch_size);
			if(identity > max_identity)
			{
				max_identity = identity;
			}

			identities.push_back(identity);
			observed_set_sizes.push_back(std::make_pair(sketch_size, intersection_size));
		}

		double identities_sum = 0;
		for(auto idty : identities)
		{
			identities_sum += idty;
		}
		//double mean_identity = identities_sum / (double)identities.size();

		// convert to true max identity
		assert(readIDs.size() == 1);
		assert(read_lengths.size() == 1);
		assert(max_identity != -1);
		max_identity = exp(-(1-max_identity));

		int read_length = *(read_lengths.begin());
		assert(read_length > param.kmerSize);

		int n_kmers = read_length - param.kmerSize + 1;

		std::vector<double> likelihoods;
		double likelihood_sum = 0;
		assert(observed_set_sizes.size());
		for(auto p : observed_set_sizes)
		{
			double likelihood = likelihood_observed_set_sizes(param.kmerSize, n_kmers, max_identity, p.first, p.second);
			assert(likelihood >= 0);
			assert(likelihood <= 1);
			likelihoods.push_back(likelihood);
			likelihood_sum += likelihood;
		}
		if((!likelihood_sum))
		{
			std::cerr << "WARNING!\n";
			std::cerr << "\t" << "likelihood_sum" << ": " << likelihood_sum << std::endl;
			std::cerr << "\t" << "max_identity" << ": " << max_identity << std::endl;
			std::cerr << "\t" << "readID" << ": " << *(readIDs.begin()) << std::endl;

			for(auto p : observed_set_sizes)
			{
				likelihood_observed_set_sizes(param.kmerSize, n_kmers, max_identity, p.first, p.second, true);
			}

			for(auto l : lines)
			{
				std::cerr << l << "\n" << std::flush;
			}

			std::cerr << "========= END ==========\n" << std::endl;
		}
		assert(likelihood_sum > 0);
		//std::cout << "max_identity: " << max_identity << "\n" << std::flush;
		//std::cout << "n_kmers: " << n_kmers < "\n";
		for(unsigned int lI = 0; lI < likelihoods.size(); lI++)
		{
			likelihoods.at(lI) = likelihoods.at(lI) / likelihood_sum;
			assert(likelihoods.at(lI) >= 0);
			assert(likelihoods.at(lI) <= 1);
			//std::cout << "\t" << observed_set_sizes.at(lI).first << " " << observed_set_sizes.at(lI).second << " " << likelihoods.at(lI) << "\n" << std::flush;
		}

		for(unsigned int lineI = 0; lineI < lines.size(); lineI++)
		{
			double mappingQuality = likelihoods.at(lineI);
			double reportedIdentity = identities.at(lineI);
			float correctedIdentity = exp(-(1-reportedIdentity));

			//std::string extension = "binomialIdentity=" + std::to_string(int(correctedIdentity*100+0.5)) + ";mappingQuality=" + std::to_string(mappingQuality);
			//lines.at(lineI) += ((fields_per_line.at(lineI) <= 12) ? " " : ";") + extension;

			std::stringstream toAdd;
			toAdd << " " << correctedIdentity*100 << " " << mappingQuality;
			lines.at(lineI) += toAdd.str();
		}
	}

	*/
}


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

void produceEM2U(std::string mappedFile, const taxonomy& T)
{
	std::ifstream resultsPerRead_EM;
	std::ifstream resultsPerRead_U;
	resultsPerRead_EM.open(mappedFile + ".EM.reads2Taxon");
	resultsPerRead_U.open(mappedFile + ".U.reads2Taxon");
	assert(resultsPerRead_EM.is_open());
	assert(resultsPerRead_U.is_open());
	
	std::map<std::string, std::map<std::string, size_t>> EM2U_details;
	std::map<std::string, std::map<std::string, size_t>> EM2U_taxonLevel;
	
	std::string line_EM;
	std::string line_U;
	while(resultsPerRead_EM.good())
	{
		assert(resultsPerRead_U.good());
		
		std::getline(resultsPerRead_EM, line_EM);
		std::getline(resultsPerRead_U, line_U);
		
		eraseNL(line_EM);
		eraseNL(line_U);
		
		if(line_EM.length() == 0)
		{
			assert(line_U.length() == 0);
			continue;
		}	
		
		std::vector<std::string> fields_EM = split(line_EM, "\t");
		std::vector<std::string> fields_U = split(line_U, "\t");
		
		assert(fields_EM.size() == 2);
		assert(fields_U.size() == 2);

		assert(fields_EM.at(0) == fields_U.at(0));
		
		if(fields_EM.at(1) == "0")
		{
			continue;
		}
		
		if((EM2U_details.count(fields_EM.at(1)) == 0) || (EM2U_details.at(fields_EM.at(1)).count(fields_U.at(1)) == 0))
		{
			EM2U_details[fields_EM.at(1)][fields_U.at(1)] = 0;
		}
		EM2U_details.at(fields_EM.at(1)).at(fields_U.at(1))++;
		
		std::string U_level = T.getNode(fields_U.at(1)).rank;
		if(fields_EM.at(1) == fields_U.at(1))
		{
			U_level = "identical";
		}
		if((EM2U_taxonLevel.count(fields_EM.at(1)) == 0) || (EM2U_taxonLevel.at(fields_EM.at(1)).count(U_level) == 0))
		{
			EM2U_taxonLevel[fields_EM.at(1)][U_level] = 0;
		}
		EM2U_taxonLevel.at(fields_EM.at(1)).at(U_level)++;		
	}
	assert(!resultsPerRead_U.good());
	
	std::string fn_EM2U_details = mappedFile + ".EM2U.details";
	std::string fn_EM2U_summary = mappedFile + ".EM2U.summary";
	std::ofstream stream_EM2U_details;
	std::ofstream stream_EM2U_summary;
	stream_EM2U_details.open(fn_EM2U_details);
	stream_EM2U_summary.open(fn_EM2U_summary);
	
	for(auto outerKey : EM2U_details)
	{
		for(auto innerKey : outerKey.second)
		{
			stream_EM2U_details << outerKey.first << "\t" << innerKey.first << "\t" << innerKey.second << "\n";
		}
	}	
	
	for(auto outerKey : EM2U_taxonLevel)
	{
		for(auto innerKey : outerKey.second)
		{
			stream_EM2U_summary << outerKey.first << "\t" << innerKey.first << "\t" << innerKey.second << "\n";
		}
	}		
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
	std::map<std::string, double> classifiedAt_freq_absolute;
	std::map<std::string, size_t> classifiedAt_reads_absolute;
	
	for(auto taxonID : combinedKeys)
	{
		std::map<std::string, std::string> upwardByLevel = T.getUpwardNodesByRanks(taxonID, targetLevels);
		upwardByLevel["definedAndHypotheticalGenomes"] = taxonID;
		upwardByLevel["definedGenomes"] = taxonID;

		{
			std::string taxonID_level;

			if(relevantTaxonIDs_mappable.count(taxonID))
			{
				taxonID_level = "definedGenomes";
			}		
			else
			{
				taxonID_level = T.getNode(taxonID).rank;
			}
			
			double combinedF = 0;
			if(std::get<0>(frequencies).count(taxonID))
				combinedF += std::get<0>(frequencies).at(taxonID);

			if(std::get<1>(frequencies).count(taxonID))
				combinedF += std::get<1>(frequencies).at(taxonID);
			
			if(std::get<2>(frequencies).count(taxonID))
				combinedF += std::get<2>(frequencies).at(taxonID);
			
			size_t combinedReads = 0;
			if(std::get<0>(readCount).count(taxonID))
				combinedReads += std::get<0>(readCount).at(taxonID);
			
			if(std::get<1>(readCount).count(taxonID))
				combinedReads += std::get<1>(readCount).at(taxonID);
			
			if(classifiedAt_freq_absolute.count(taxonID_level) == 0)
			{
				classifiedAt_freq_absolute[taxonID_level] = 0;
				classifiedAt_reads_absolute[taxonID_level] = 0;
			}
			classifiedAt_freq_absolute.at(taxonID_level) += combinedF;
			classifiedAt_reads_absolute.at(taxonID_level) += combinedReads;
				
		}
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

			double combinedF = 0;
			if(std::get<0>(frequencies).count(taxonID))
			{
				std::get<0>(frequencies_perLevel.at(level)).at(levelValue) += std::get<0>(frequencies).at(taxonID);
				combinedF += std::get<0>(frequencies).at(taxonID);
			}

			if(std::get<1>(frequencies).count(taxonID))
			{
				std::get<1>(frequencies_perLevel.at(level)).at(levelValue) += std::get<1>(frequencies).at(taxonID);
				combinedF += std::get<1>(frequencies).at(taxonID);
			}

			if(std::get<2>(frequencies).count(taxonID))
			{
				std::get<2>(frequencies_perLevel.at(level)).at(levelValue) += std::get<2>(frequencies).at(taxonID);
				combinedF += std::get<2>(frequencies).at(taxonID);
			}


			size_t combinedReads = 0;
			if(std::get<0>(readCount).count(taxonID))
			{
				readCount_perLevel.at(level).first.at(levelValue) += std::get<0>(readCount).at(taxonID);
				combinedReads += std::get<0>(readCount).at(taxonID);
			}

			if(std::get<1>(readCount).count(taxonID))
			{
				readCount_perLevel.at(level).second.at(levelValue) += std::get<1>(readCount).at(taxonID);
				combinedReads += std::get<1>(readCount).at(taxonID);
			}
			
		}
	}
	
	std::string classifiedAt = outputFN + ".absoluteClassifiedAt";
	std::ofstream strout_absolute(classifiedAt);
	strout_absolute << "Level" << "\t" << "f" << "\t" << "nReads" << "\n";
	for(auto lV : classifiedAt_freq_absolute)
	{
		strout_absolute << lV.first << "\t" << lV.second << "\t" << classifiedAt_reads_absolute.at(lV.first) << "\n";
	}
	
	strout_absolute.close();

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

oneMappingLocation_U getHighestProbabilityMapping_U(const std::vector<oneMappingLocation_U>& locations)
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



std::pair<int, int> getMinMaxIdentities(std::string mappedFile, const skch::Parameters& parameters)
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
	callBackForAllReads(mappedFile, processOneRead, parameters);
	
	assert(maxIdentity > 1);
	
	return std::make_pair(minIdentity, maxIdentity);
}

void normalize_f(std::pair<std::map<std::string, double>, std::map<std::string, double>>& fToN)
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

	// to do remove
	// std::cout << "f_sum: " << f_sum << "\n" << std::flush;

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
}

void normalize_f_triplet(std::tuple<std::map<std::string, double>, std::map<std::string, double>, std::map<std::string, double>>& fToN)
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

	// to do remove
	// std::cout << "f_sum: " << f_sum << "\n" << std::flush;

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
}



void doU(std::string DBdir, std::string mappedFile, size_t minimumReadsPerBestContig, const skch::Parameters& parameters)
{
	// unsigned int round_first_unknown = 5;

	taxonomy T(DBdir+"/taxonomy");

	std::set<std::string> taxonIDsInMappings = getTaxonIDsFromMappingsFile(mappedFile);
	std::set<std::string> mappableTaxonIDs = getDirectlyMappableTaxonIDs(DBdir);
	
	std::string fn_fittedLengthAndIdentities = mappedFile + ".EM.lengthAndIdentitiesPerMappingUnit";
	if(! fileExists(fn_fittedLengthAndIdentities))
	{
		std::cerr << "\n\nERROR: File " << fn_fittedLengthAndIdentities << " not existing.\n\nThis file is generated automatically by the EM step. Run the EM step first.\n\n";
	}
	
	std::pair<int, int> identity_minmax_inAlignments = getMinMaxIdentities(mappedFile, parameters);

	identityAndReadLengthHistogram iAndL;
	iAndL.readFromEMOutput(fn_fittedLengthAndIdentities, identity_minmax_inAlignments, minimumReadsPerBestContig);

	std::string fn_tree_selfSimilarities = DBdir + "/selfSimilarities.txt";
	treeAdjustedIdentities tAI;
	tAI.readFromFile(fn_tree_selfSimilarities, taxonIDsInMappings, T);

	identityManager iM(iAndL, tAI);
	generateUnknownMapQFile(DBdir, mappedFile, iM, T, parameters);

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
	std::map<std::string, int> indirectUpwardNodes_nSourceGenomes;
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
				if(indirectUpwardNodes_nSourceGenomes.count(uTI) == 0)
				{
					indirectUpwardNodes_nSourceGenomes[uTI] = tAI.indirectAttachmentNode_sources(uTI);
				}
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


	auto get_mappings_with_P = [&](const std::pair<std::map<std::string, double>, std::map<std::string, double>>& f, const std::vector<std::string>& readLines, std::vector<oneMappingLocation_U>& mappings_with_P_ret) -> double
	{
		assert(readLines.size() > 0);
		mappings_with_P_ret.clear();

		std::set<std::string> thisRead_direct;
		std::set<std::string> thisRead_indirect;

		double l_read = 0;
		for(auto rL : readLines)
		{
			assert(rL.length());
			std::vector<std::string> rL_fields = split(rL, " ");
			if(!(rL_fields.size() == 5))
			{
				std::cerr << "Weird line for parsing:\n\t'" << rL << "'\n" << std::flush;
			}
			assert(rL_fields.size() == 5);

			oneMappingLocation_U thisMapping;
			thisMapping.taxonID = rL_fields.at(1);
			thisMapping.direct = std::stoi(rL_fields.at(2));
			thisMapping.mapQ = std::stod(rL_fields.at(3));
			thisMapping.originalIdentity = std::stod(rL_fields.at(4));
			assert((thisMapping.mapQ >= 0) && (thisMapping.mapQ <= 1));

			double l;
			if(thisMapping.direct)
			{
				assert(thisRead_direct.count(thisMapping.taxonID) == 0);
				thisRead_direct.insert(thisMapping.taxonID);
				l = f.first.at(thisMapping.taxonID) * thisMapping.mapQ;
			}
			else
			{
				assert(thisRead_indirect.count(thisMapping.taxonID) == 0);
				thisRead_indirect.insert(thisMapping.taxonID);
				l = f.second.at(thisMapping.taxonID) * thisMapping.mapQ;
			}
			l_read += l;
			thisMapping.p = l;
			mappings_with_P_ret.push_back(thisMapping);
		}
		if(l_read <= 0)
		{
			std::cerr << "l_read: " << l_read << "\n";
			for(auto rL : readLines)
			{
				std::cerr << "\t" << rL << "\n";
			}
			std::cerr << "Source file: " << mappedFile +".mapQ_U" << "\n";
			std::cerr << std::flush;
		}
		assert(l_read > 0);

		for(auto& m : mappings_with_P_ret)
		{
			m.p /= l_read;
		}

		return l_read;
	};

	auto sumF = [](const std::pair<std::map<std::string, double>, std::map<std::string, double>>& f) -> double {
		double S = 0;
		for(auto tF : f.first)
		{
			S += tF.second;
		}
		for(auto tF : f.second)
		{
			S += tF.second;
		}
		return S;
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

		//std::map<std::string, size_t> bestMappings_perTaxonID;

		size_t processedRead = 0;
		// int printedReads = 0;
		std::function<void(const std::vector<std::string>&)> processOneRead = [&](const std::vector<std::string>& readLines) -> void
		{
			processedRead++;
			if((processedRead % 10000) == 0)
			{
				std::cout << "\r EM-U round " << EMiteration << ", read << " << processedRead << " / " << mappingStats.at("ReadsMapped") << "   " << std::flush;
			}
			std::vector<oneMappingLocation_U> mappings_with_P;
			double l_read = get_mappings_with_P(f, readLines, mappings_with_P);

			/*
			std::string bestMapping_taxonID;
			double bestMapping_l;
			bool bestMapping_isDirect;

			if((bestMapping_taxonID.length() == 0) || (bestMapping_l < l))
			{
				bestMapping_taxonID = taxonID;
				bestMapping_isDirect = direct;
				bestMapping_l = l;
			}

			assert(readLines.size() > 0);
			double l_read = 0;
			std::string bestMapping_taxonID;

			double bestMapping_l;
			bool bestMapping_isDirect;

			for(auto rL : readLines)
			{
				std::vector<std::string> rL_fields = split(rL, " ");
				assert(rL_fields.size() == 4);
				std::string taxonID = rL_fields.at(1);
				bool direct = std::stoi(rL_fields.at(2));
				double mapQ = std::stod(rL_fields.at(3));
				assert((mapQ >= 0) && (mapQ <= 1));

				double l;
				if(direct)
				{
					assert(thisRead_direct.count(taxonID) == 0);
					l = f.first.at(taxonID) * mapQ;
					thisRead_direct[taxonID] = l;
					l_read += l;
				}
				else
				{
					assert(thisRead_indirect.count(taxonID) == 0);
					l = f.second.at(taxonID) * mapQ;
					thisRead_indirect[taxonID] = l;
					l_read += l;
				}

				if((bestMapping_taxonID.length() == 0) || (bestMapping_l < l))
				{
					bestMapping_taxonID = taxonID;
					bestMapping_isDirect = direct;
					bestMapping_l = l;
				}
			}
			*/

			assert(l_read > 0);
			ll_thisIteration_mapped += log(l_read);

			double p_read = 0;
			for(auto& l : mappings_with_P)
			{
				if(l.direct)
				{
					f_nextIteration.first.at(l.taxonID) += l.p;
				}
				else
				{
					f_nextIteration.second.at(l.taxonID) += l.p;

				}
				p_read += l.p;
			}

			assert(abs(1 - p_read) <= 1e-3);

			/*
			if(bestMappings_perTaxonID.count(bestMapping_taxonID) == 0)
			{
				bestMappings_perTaxonID[bestMapping_taxonID] = 0;
			}
			bestMappings_perTaxonID.at(bestMapping_taxonID)++;
			*/
		};

		callBackForAllReads(mappedFile+".mapQ_U", processOneRead, parameters);
		std::cout << "\n";

		double f_sum_preNormalization = sumF(f_nextIteration);
		assert(abs(nMapped - f_sum_preNormalization) <= 1e-2);

		normalize_f(f_nextIteration);

		double ll_thisIteration = ll_thisIteration_mapped;

		std::cout << "\tLog likelihood: " << ll_thisIteration << std::endl;
		std::cout << "\t\tcontribution mapped reads  : " << ll_thisIteration_mapped << std::endl;
		std::cout << "\t\tcontribution unmapped reads: " << "NA" << std::endl;

		if(EMiteration > 0)
		{
			double ll_diff = ll_thisIteration - ll_lastIteration;
			assert(ll_diff >= 0);

			double ll_relative = ll_thisIteration/ll_lastIteration;

			std::cout << "\tImprovement: " << ll_diff << std::endl;
			std::cout << "\tRelative   : " << ll_relative << std::endl;

			double ll_relative_imp = 1 - ll_relative;

			if((ll_diff <= 1) && (ll_relative_imp < 0.0001))
			{
				continueEM = false;
			}
		}

		if(EMiteration > 3)
		{
			// continueEM = false; // todo
		}

		double f_sum_postNormalization  = sumF(f_nextIteration);
		assert(abs(1 - f_sum_postNormalization) <= 1e-3);

		f = f_nextIteration;
		EMiteration++;
		ll_lastIteration = ll_thisIteration;
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
		std::string readID = getReadIDFromReadLines(readLines);

		std::vector<oneMappingLocation_U> mappings_with_P;
		get_mappings_with_P(f, readLines, mappings_with_P);

		oneMappingLocation_U bestMapping = getHighestProbabilityMapping_U(mappings_with_P);

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
				bestMapping.originalIdentity << "\t" <<
				bestMapping.readLength << "\n";

		strout_reads_taxonIDs << readID << "\t" << bestMapping.taxonID << "\n";
	};

	callBackForAllReads(mappedFile+".mapQ_U", processOneRead_final, parameters);

	// long-enough-but-unmapped reads are set to unassigned
	std::vector<std::string> readIDs_notMapped_despiteLongEnough = getUnmappedReadsIDs(mappedFile);
	for(auto readID : readIDs_notMapped_despiteLongEnough)
	{
		strout_reads_taxonIDs << readID << "\t" << 0 << "\n";
	}

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

	produceEM2U(mappedFile, T);
	
	std::cout << "\n\nDone. Frequencies output in " << "\n\t" << output_pot_frequencies << "\n" << std::flush;

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
