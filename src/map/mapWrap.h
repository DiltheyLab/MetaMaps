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

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"

class mapWrap {
protected:
	void unifyFiles(std::vector<std::string> mappingFiles, std::vector<std::string> querySequences)
	{
		assert( 1 == 0 );
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

		skch::Sketch referSketch(param, maximumMemory, &mapFunction);
		unifyFiles(individualMappingFiles, param.querySequences);
	}


};

#endif /* MAP_MAPWRAP_H_ */
