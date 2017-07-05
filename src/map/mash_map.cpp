/**
 * @file    mash_map.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <assert.h>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/parseCmdArgs.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"

//External includes
#include "common/argvparser.hpp"
#include <boost/filesystem.hpp>

void mapAgainstPrefix(const skch::Parameters& parameters, std::string prefix)
{
	skch::Parameters useParameters = parameters;

	// check that parameters are compatible
	{
		std::string fn_serialize_arguments = prefix + ".arguments";
		std::ifstream arguments_serialization_istream(fn_serialize_arguments.c_str());
		if(! arguments_serialization_istream.is_open())
		{
		  throw std::runtime_error("Expected file " + fn_serialize_arguments + " not found - have you supplied a valid prefix?");
		}
		boost::archive::text_iarchive arguments_archive(arguments_serialization_istream);
		skch::Parameters restoredParameters;
		arguments_archive >> restoredParameters;
		arguments_serialization_istream.close();

		useParameters.alphabetSize = restoredParameters.alphabetSize;
		useParameters.kmerSize = restoredParameters.kmerSize;
		useParameters.minReadLength = restoredParameters.minReadLength;
		useParameters.p_value = restoredParameters.p_value;
		useParameters.percentageIdentity = restoredParameters.percentageIdentity;
		useParameters.windowSize = restoredParameters.windowSize;
		useParameters.referenceSize = restoredParameters.referenceSize;

		std::cout << "Parameters restored from index " << fn_serialize_arguments << "\n";
		std::cout << "\t" << "- alphabetSize: " << useParameters.alphabetSize << "\n";
		std::cout << "\t" << "- kmerSize: " << useParameters.kmerSize << "\n";
		std::cout << "\t" << "- minReadLength: " << useParameters.minReadLength << "\n";
		std::cout << "\t" << "- p_value: " << useParameters.p_value << "\n";
		std::cout << "\t" << "- percentageIdentity: " << useParameters.percentageIdentity << "\n";
		std::cout << "\t" << "- windowSize: " << useParameters.windowSize << "\n";
		std::cout << "\n" << std::flush;
	}

	std::string outputPrefix = parameters.outFileName;
	int read_index_files = 0;
	int runningIndexNumber = 1;
	while(boost::filesystem::exists( prefix + "." + std::to_string(runningIndexNumber)))
	{
		read_index_files++;

		std::ifstream sketch_serialization_istream(prefix + "." + std::to_string(runningIndexNumber));
		if(! sketch_serialization_istream.is_open())
		{
				throw std::runtime_error("Cannot open file " + prefix + "." + std::to_string(runningIndexNumber) + " for reading");
		}
		boost::archive::text_iarchive sketch_archive(sketch_serialization_istream);

		skch::Sketch referSketchI;

		sketch_archive >> referSketchI;

		sketch_serialization_istream.close();

		std::cout << "Restored archive from " << prefix << "." << std::to_string(runningIndexNumber) << std::endl;

		useParameters.outFileName = outputPrefix + "." + std::to_string(runningIndexNumber);
		skch::Map mapper = skch::Map(useParameters, referSketchI);

		runningIndexNumber++;
	}

	if(read_index_files == 0)
		throw std::runtime_error("No index files read for prefix " + prefix);
}

int main(int argc, char** argv)
{
	{
		/*
		{
			std::ifstream sketch_serialization_istream("index.1");
			if(! sketch_serialization_istream.is_open())
			{
					throw std::runtime_error("Cannot open file index.1 for reading");
			}
			boost::archive::text_iarchive sketch_archive(sketch_serialization_istream);

			skch::Sketch referSketch;

			sketch_archive >> referSketch;

		}
		*/

		skch::Parameters parameters;        //sketching and mapping parameters

		parameters.refSequences.push_back("C:/Temp/jh1_long.fasta");
		parameters.querySequences.push_back("C:/Temp/reads.fasta");
		parameters.outFileName = "mappingResults";

	    parameters.referenceSize = skch::CommonFunc::getReferenceSize(parameters.refSequences);
	    parameters.alphabetSize = 4;
		parameters.reportAll = true;
		parameters.kmerSize = 16;

		parameters.p_value = 1e-03;

		parameters.minReadLength = 50;
		parameters.percentageIdentity = 85;

		parameters.windowSize = 50;

		/*
		parameters.windowSize = skch::Stat::recommendedWindowSize(
				parameters.p_value,
				parameters.kmerSize,
				parameters.alphabetSize,
				parameters.percentageIdentity,
				parameters.minReadLength,
				parameters.referenceSize);
		*/


		std::cout << "Referenze size / window size: " << parameters.referenceSize << " / " << parameters.windowSize << "\n" << std::flush;

		skch::Sketch referSketch(parameters, "index", 500 * 1024);

		mapAgainstPrefix(parameters, "index");

		assert( 5 == 4);


		std::string fn_serialize_sketch = "sketch_serialized";

		{
			skch::Sketch referSketch(parameters);

			skch::Map mapper = skch::Map(parameters, referSketch);

			std::cout << "INFO, skch::main, mapping results saved in : " << parameters.outFileName << std::endl;

			std::cout << "\n\nMapping done, results in " << parameters.outFileName << "\n\n" << std::flush;

			std::ofstream sketch_serialization_ostream(fn_serialize_sketch.c_str());
			if(! sketch_serialization_ostream.is_open())
			{
					throw std::runtime_error("Cannot open file " + fn_serialize_sketch + " for serialization");
			}
			boost::archive::text_oarchive sketch_archive(sketch_serialization_ostream);
			sketch_archive & referSketch;
		}

		{
			std::ifstream sketch_serialization_istream(fn_serialize_sketch.c_str());
			if(! sketch_serialization_istream.is_open())
			{
					throw std::runtime_error("Cannot open file " + fn_serialize_sketch + " for reading");
			}
			boost::archive::text_iarchive sketch_archive(sketch_serialization_istream);

			skch::Sketch referSketch;

			sketch_archive >> referSketch;

			skch::Map mapper = skch::Map(parameters, referSketch);

			std::cout << "INFO, skch::main, mapping results saved in : " << parameters.outFileName << std::endl;

			std::cout << "\n\nMapping done, results in " << parameters.outFileName << "\n\n" << std::flush;


		}

		std::cout << "\n\nDone.\n\n" << std::flush;
	}
	assert( 2 == 4);

  CommandLineProcessing::ArgvParser cmd;

  //Setup command line options
  skch::initCmdParser(cmd);

  //Parse command line arguements   
  skch::Parameters parameters;        //sketching and mapping parameters

  skch::parseandSave(argc, argv, cmd, parameters);   

  auto t0 = skch::Time::now();

  //Build the sketch for reference
  skch::Sketch referSketch(parameters);

  std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
  std::cout << "INFO, skch::main, Time spent sketching the reference : " << timeRefSketch.count() << " sec" << std::endl;

  //Map the sequences in query file
  t0 = skch::Time::now();

  skch::Map mapper = skch::Map(parameters, referSketch);

  std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;
  std::cout << "INFO, skch::main, Time spent mapping the query : " << timeMapQuery.count() << " sec" << std::endl;

  std::cout << "INFO, skch::main, mapping results saved in : " << parameters.outFileName << std::endl;

}
