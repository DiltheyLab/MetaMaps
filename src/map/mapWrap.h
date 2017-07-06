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
#include <boost/archive/text_oarchive.hpp>
#include <fstream>
#include <zlib.h>

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"
#include "meta/util.h"

class mapWrap {
protected:
	void unifyFiles(std::string unifiedOutputFN, skch::Parameters param, std::vector<std::string> mappingFiles, std::vector<std::string> querySequences)
	{
		std::cerr << "UNIFY!!" << std::endl;

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

		auto queryOpenFileForReadData = [&](int fI, const std::string& readID) -> std::string
		{
			lineBuffer.clear();
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
						if(thisLine_readID == readID)
						{
							lineBuffer.append(line+"\n");
							posBeforeStart = openFiles.at(fI)->tellg();
						}
						else
						{
							break;
						}
					}
				}
				openFiles.at(fI)->seekg(posBeforeStart);
				return lineBuffer;
			}
			else
			{
				return "";
			}

		};

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
				//Is the read too short?
				if(len < param.windowSize || len < param.kmerSize || len < param.minReadLength)
				{
					continue;
				}
				else
				{
					std::string readID = seq->name.s;
					for(unsigned int fI = 0; fI < openFiles.size(); fI++)
					{
						combinedOutput << queryOpenFileForReadData(fI, readID);
					}
				}
			}

			//Close the input file
			kseq_destroy(seq);
			gzclose(fp);
		}

		for(auto iS : openFiles)
		{
			iS->close();
			delete(iS);
		}

		combinedOutput.close();
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
		boost::archive::text_oarchive arguments_archive(arguments_serialization_ostream);
		arguments_archive & param;
		arguments_serialization_ostream.close();

		std::vector<std::string> generatedIndexFiles;
		std::function<void(skch::Sketch*, size_t)> indexStoreFunction = [&](skch::Sketch* referSketch, size_t N)
		{
			std::string outputFn = prefix + "." + std::to_string(N);
			std::ofstream ostream(outputFn.c_str());
			if(! ostream.is_open())
			{
			  throw std::runtime_error("Cannot open file " + outputFn + " for serialization.");
			}
			boost::archive::text_oarchive archive(ostream);
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
			boost::archive::text_iarchive sketch_archive(sketch_serialization_istream);

			skch::Sketch referSketchI;

			sketch_archive >> referSketchI;

			sketch_serialization_istream.close();

			std::cout << "Restored archive from " << indexFile << "\n";

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
