/*
 * mapWrap.h
 *
 *  Created on: Jul 6, 2017
 *      Author: diltheyat
 */

#ifndef MAP_MAPWRAP_H_
#define MAP_MAPWRAP_H_

#include <string>
#include <assert.h>
#include <boost/filesystem.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <fstream>
#include <zlib.h>
#include <utility>
#include <boost/math/distributions/binomial.hpp>

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"
#include "meta/util.h"

class mapWrap {
protected:
	void unifyFiles(std::string unifiedOutputFN, const skch::Parameters& param, std::vector<std::string> mappingFiles, std::vector<std::string> querySequences)
	{
		bool paranoid = true;
		std::set<std::string> processedReadIDs;

		std::ofstream combinedOutput(unifiedOutputFN);
		assert(combinedOutput.is_open());

		std::vector<std::ifstream*> openFiles;
		for(auto resultsFile : mappingFiles)
		{
			std::ifstream* i = new std::ifstream(resultsFile);
			assert(i->is_open());
			openFiles.push_back(i);
		}

		std::string lineBuffer;
		lineBuffer.reserve(100000);

		auto queryOpenFileForReadData = [&](int fI, const std::string& readID) -> std::vector<std::string>
		{
			std::vector<std::string> forReturn;
			if(openFiles.at(fI)->good())
			{
				std::streampos posBeforeStart = openFiles.at(fI)->tellg();
				std::string line;
				while(openFiles.at(fI)->good())
				{
					std::getline(*(openFiles.at(fI)), line);
					size_t pos_first_space = line.find(' ');
					if(pos_first_space == std::string::npos)
					{
						break;
					}
					else
					{
						std::string thisLine_readID = line.substr(0, pos_first_space);
						if(paranoid && (processedReadIDs.count(thisLine_readID) != 0))
						{
							std::cerr << "Seems that read ID " << thisLine_readID << " has already been processed - this target ID " << readID << "\n" << std::endl;
							exit(1);
						}
						if(thisLine_readID == readID)
						{
							forReturn.push_back(line);
							posBeforeStart = openFiles.at(fI)->tellg();
						}
						else
						{
							break;
						}
					}
				}
				openFiles.at(fI)->seekg(posBeforeStart);
				return forReturn;
			}
			else
			{
				return forReturn;
			}
		};

		std::ofstream metaLengths(unifiedOutputFN+".meta.unmappedReadsLengths");
		assert(metaLengths.is_open());

		size_t totalInputReads = 0;
		size_t inputReads_mapped = 0;
		size_t inputReads_tooShort = 0;
		size_t inputReads_notMapped = 0;
		std::string lastProcessedReadID;
		for(const auto& qSF : querySequences)
		{
			//Open the file using kseq
			FILE *file = fopen(qSF.c_str(), "r");
			gzFile fp = gzdopen(fileno(file), "r");
			kseq_t *seq = kseq_init(fp);

			//size of sequence
			skch::offset_t len;

			while ((len = kseq_read(seq)) >= 0)
			{
				totalInputReads++;

				//Is the read too short?
				if(len < param.windowSize || len < param.kmerSize || len < param.minReadLength)
				{
					inputReads_tooShort++;
					continue;
				}
				else
				{
					std::string readID = seq->name.s;
					std::vector<std::string> combinedLines;
					for(unsigned int fI = 0; fI < openFiles.size(); fI++)
					{
						std::vector<std::string> file_lines = queryOpenFileForReadData(fI, readID);
						combinedLines.insert(combinedLines.end(), file_lines.begin(), file_lines.end());
					}

					if(combinedLines.size() == 0)
					{
						inputReads_notMapped++;
						int len_int = len;
						metaLengths << len_int << "\t" << readID << "\n";
						//metaLengths << (int)len << "\t" << readID < "\n";
					}
					else
					{
						inputReads_mapped++;
					}
					addMappingQualities(param, combinedLines);
					for(auto l : combinedLines)
					{
						combinedOutput << l << "\n";
					}
					if(paranoid)
					{
						processedReadIDs.insert(readID);
					}
					lastProcessedReadID = readID;
				}
			}

			//Close the input file
			kseq_destroy(seq);
			gzclose(fp);
		}

		size_t reads_mappable = totalInputReads - inputReads_tooShort;
		for(size_t fileI = 0; fileI < openFiles.size(); fileI++)
		{
			std::ifstream* iS = openFiles.at(fileI);
			if(iS->good() && (reads_mappable != 0))
			{
				std::cerr << "Error: output file " << fileI << " / " << openFiles.size() << " not completely processed.\n" << std::flush;
				std::cerr << "\tFile: " << mappingFiles.at(fileI) << std::endl;
				std::cerr << "\tLast processed read ID: " << lastProcessedReadID << std::endl;
				exit(1);
			}
			iS->close();
			delete(iS);
		}

		std::ofstream metaOutput(unifiedOutputFN+".meta");
		assert(metaOutput.is_open());
		metaOutput << "TotalReads" << " " << totalInputReads << "\n";
		metaOutput << "ReadsTooShort" << " " << inputReads_tooShort << "\n";
		metaOutput << "ReadsMapped" << " " << inputReads_mapped << "\n";
		metaOutput << "ReadsNotMapped" << " " << inputReads_notMapped << "\n";
		metaOutput.close();

		metaLengths.close();

		combinedOutput.close();
		
		for(auto resultsFile : mappingFiles)
		{
			boost::filesystem::wpath F(resultsFile);
			assert(boost::filesystem::remove(F));
		}		
		
		std::ofstream parametersStream(unifiedOutputFN+".parameters");
		assert(parametersStream.is_open());
		
		parametersStream << "kmerSize" << " " << param.kmerSize << "\n";
		parametersStream << "windowSize" << " " << param.windowSize << "\n";
		parametersStream << "minReadLength" << " " << param.minReadLength << "\n";
		parametersStream << "alphabetSize" << " " << param.alphabetSize << "\n";
		parametersStream << "referenceSize" << " " << param.referenceSize << "\n";
		parametersStream << "percentageIdentity" << " " << param.percentageIdentity << "\n";
		parametersStream << "p_value" << " " << param.p_value << "\n";
		parametersStream << "refSequences" << " " << param.refSequences << "\n";
		parametersStream << "querySequences" << " " << param.querySequences << "\n";
		parametersStream << "outFileName" << " " << param.outFileName << "\n";
		parametersStream << "reportAll" << " " << param.reportAll << "\n";
		parametersStream << "index" << " " << param.index << "\n";
		parametersStream << "maximumMemory" << " " << param.maximumMemory << "\n";
			
	}

	void addMappingQualities(const skch::Parameters& param, std::vector<std::string>& lines)
	{
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
	}

	double likelihood_observed_set_sizes(int k, int n_kmers, double identity, int sketch_size, int intersection_size, bool verbose = false)
	{
		assert(intersection_size <= sketch_size);
		double p_kMer_survival = std::pow(identity, k);
		double E_surviving_kMers = p_kMer_survival * n_kmers;
		double E_surviving_kMers_int = std::round(E_surviving_kMers);
		double E_union_size = n_kmers + (n_kmers - E_surviving_kMers_int);
		double E_intersection_size = E_surviving_kMers_int;
		double likelihood = boost::math::pdf(boost::math::binomial_distribution<>(sketch_size, E_intersection_size/E_union_size), intersection_size);
		
		if(verbose)
		{
			std::cerr << "OUTPUT: likelihood = " << likelihood << "\n";
			std::cerr << "k" << ": " << k << std::endl;
			std::cerr << "n_kmers" << ": " << n_kmers << std::endl;
			std::cerr << "identity" << ": " << identity << std::endl;
			std::cerr << "sketch_size" << ": " << sketch_size << std::endl;
			std::cerr << "intersection_size" << ": " << intersection_size << std::endl;
			std::cerr << "E_intersection_size" << ": " << E_intersection_size << std::endl;
			std::cerr << "E_union_size" << ": " << E_union_size << std::endl;			
			std::cerr << std::endl;
		}
		
		return likelihood;
	}

public:
	mapWrap()
	{

	}

	void createIndex(skch::Parameters param, std::string prefix, size_t maximumMemory)
	{
		assert(param.outFileName.length() == 0);
		param.outFileName = "";

		std::ofstream indexBuilt (prefix + ".index");
		assert(indexBuilt.is_open());
		indexBuilt << 0 << "\n";
		indexBuilt.close();

		std::string fn_serialize_arguments = prefix + ".arguments";
		std::ofstream arguments_serialization_ostream(fn_serialize_arguments.c_str());
		if(! arguments_serialization_ostream.is_open())
		{
			throw std::runtime_error("Cannot open file " + fn_serialize_arguments + " for serialization.");
		}
		{
			boost::archive::text_oarchive arguments_archive(arguments_serialization_ostream);
			arguments_archive & param;
			arguments_serialization_ostream.close();
		}
		std::vector<std::string> generatedIndexFiles;
		std::function<void(skch::Sketch*, size_t)> indexStoreFunction = [&](skch::Sketch* referSketch, size_t N)
		{
			std::string outputFn = prefix + "." + std::to_string(N);
			std::ofstream ostream(outputFn.c_str());
			if(! ostream.is_open())
			{
			  throw std::runtime_error("Cannot open file " + outputFn + " for serialization.");
			}
			boost::archive::binary_oarchive archive(ostream);
			archive & (*referSketch);
			ostream.close();
			std::cout << "Stored state in file " << outputFn << "\n" << std::flush;
			generatedIndexFiles.push_back(outputFn);
		};

		skch::Sketch referSketch(param, maximumMemory, &indexStoreFunction);

		std::ofstream indexBuilt2 (prefix + ".index");
		assert(indexBuilt2.is_open());
		indexBuilt2 << 1 << "\n";
		for(auto f : generatedIndexFiles)
		{
			indexBuilt2 << f << "\n";
		}
		indexBuilt2.close();
	}

	void mapDirectly(skch::Parameters param, size_t maximumMemory)
	{
		std::string prefix = param.outFileName;
		std::vector<std::string> individualMappingFiles;

		std::function<void(skch::Sketch*, size_t)> mapFunction = [&](skch::Sketch* referSketch, size_t N)
		{
			std::string outputFn = prefix + "." + std::to_string(N);
			skch::Parameters localP = param;
			localP.outFileName = outputFn;

			skch::Map mapper = skch::Map(localP, *referSketch);

			individualMappingFiles.push_back(outputFn);
		};

		param.outFileName = "";
		skch::Sketch referSketch(param, maximumMemory, &mapFunction);

		unifyFiles(prefix, param, individualMappingFiles, param.querySequences);
	}

	void mapAgainstIndex(skch::Parameters param, std::string indexPrefix)
	{
		std::string fn_indexOK = indexPrefix + ".index";

		std::vector<std::string> indexFiles;
		{
			std::vector<std::string> fn_indexOK_lines;
			std::ifstream fn_indexOK_stream (fn_indexOK.c_str());
			assert(fn_indexOK_stream.is_open());
			std::string line;
			while(fn_indexOK_stream.good())
			{
				std::getline(fn_indexOK_stream, line);
				eraseNL(line);
				if(line.length())
					fn_indexOK_lines.push_back(line);
			}
			fn_indexOK_stream.clear();
			if((fn_indexOK_lines.size() == 0) || (fn_indexOK_lines.at(0) != "1"))
			{
				std::cerr << "The file " << fn_indexOK << " does not indicate that index " << indexPrefix << " was built successfully, abort." << std::endl;
				exit(1);
			}
			if(fn_indexOK_lines.size() < 2)
			{
				std::cerr << "Index " << indexPrefix << " was built successfully, but no index files present?" << std::endl;
				exit(1);
			}
			indexFiles = std::vector<std::string>(fn_indexOK_lines.begin()+1, fn_indexOK_lines.end());
		}

		std::string fn_serialize_arguments = indexPrefix + ".arguments";
		std::ifstream arguments_serialization_istream(fn_serialize_arguments.c_str());
		if(! arguments_serialization_istream.is_open())
		{
			std::cerr << "Expected file " << fn_serialize_arguments << " not found - have you supplied a valid index?" << std::endl;
			exit(1);
		}
		
		boost::archive::text_iarchive arguments_archive(arguments_serialization_istream);
		skch::Parameters restoredParameters;
		arguments_archive >> restoredParameters;
		arguments_serialization_istream.close();

		skch::Parameters useParameters;
		useParameters.alphabetSize = restoredParameters.alphabetSize;
		useParameters.kmerSize = restoredParameters.kmerSize;
		useParameters.minReadLength = restoredParameters.minReadLength;
		useParameters.p_value = restoredParameters.p_value;
		useParameters.percentageIdentity = restoredParameters.percentageIdentity;
		useParameters.windowSize = restoredParameters.windowSize;
		useParameters.referenceSize = restoredParameters.referenceSize;
		useParameters.reportAll = param.reportAll;
  
		std::cout << "Parameters restored from index " << fn_serialize_arguments << "\n";
		std::cout << "\t" << "- alphabetSize: " << useParameters.alphabetSize << "\n";
		std::cout << "\t" << "- kmerSize: " << useParameters.kmerSize << "\n";
		std::cout << "\t" << "- minReadLength: " << useParameters.minReadLength << "\n";
		std::cout << "\t" << "- p_value: " << useParameters.p_value << "\n";
		std::cout << "\t" << "- percentageIdentity: " << useParameters.percentageIdentity << "\n";
		std::cout << "\t" << "- windowSize: " << useParameters.windowSize << "\n";
		std::cout << "\n" << std::flush;

		std::vector<std::string> outputFiles;
		std::string outputPrefix = param.outFileName;
		size_t indexFileI = 0;
		for(auto indexFile : indexFiles)
		{
			std::ifstream sketch_serialization_istream(indexFile);
			if(! sketch_serialization_istream.is_open())
			{
					std::cerr << "Cannot open file " << indexFile << " for reading -- invalid index " << indexPrefix << std::endl;
			}
			boost::archive::binary_iarchive sketch_archive(sketch_serialization_istream);

			skch::Sketch referSketchI;

			sketch_archive >> referSketchI;

			sketch_serialization_istream.close();

			useParameters.outFileName = outputPrefix + "." + std::to_string(indexFileI);
			useParameters.querySequences = param.querySequences;
			skch::Map mapper = skch::Map(useParameters, referSketchI);

			outputFiles.push_back(useParameters.outFileName);
			indexFileI++;
		}

		unifyFiles(outputPrefix, useParameters, outputFiles, param.querySequences);
	}
};

#endif /* MAP_MAPWRAP_H_ */
