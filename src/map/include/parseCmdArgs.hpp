/**
 * @file    parseCmdArgs.hpp
 * @brief   Functionality related to command line parsing for indexing and mapping
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PARSE_CMD_HPP 
#define PARSE_CMD_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <exception>
#include <stdexcept>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/commonFunc.hpp"

//More includes

//External includes
#include "common/argvparser.hpp"

namespace skch
{

  /**
   * @brief           Initialize the command line argument parser 
   * @param[out] cmd  command line parser object
   */
  void initCmdParser(CommandLineProcessing::ArgvParser &cmd, std::string mode)
  {
    cmd.setIntroductoryDescription("Approximate read mapper based on Jaccard similarity");
    cmd.setIntroductoryDescription("");

    cmd.setHelpOption("h", "help", "Print this help page");

    if(mode == "mapAgainstIndex")
    {
        cmd.defineOption("index", "output prefix for ", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("index","i");

        cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("query","q");

        cmd.defineOption("output", "output file", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("output","o");

        cmd.defineOption("threads", "count of threads for parallel execution [default : 1]", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("threads","t");

        cmd.defineOption("all", "report all the mapping locations for a read, default is to consider few best ones");
    }
    else if(mode == "mapDirectly" || mode == "index")
    {
        cmd.defineOption("reference", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("reference", "r");

        cmd.defineOption("kmer", "kmer size <= 16 [default 16 (DNA)]", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("kmer","k");

        cmd.defineOption("pval", "p-value cutoff, used to determine window/sketch sizes [default e-03]", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("pval","p");

        cmd.defineOption("maxmemory", "maximum memory, in GB [default e-03]", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("maxmemory","mm");

        cmd.defineOption("window", "window size [default : computed using pvalue cutoff]\n\
    P-value is not considered if a window value is provided. Lower window size implies denser sketch",
            ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("window","w");

        cmd.defineOption("minReadLen", "minimum read length to map [default : 2000]", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("minReadLen","m");

        cmd.defineOption("perc_identity", "threshold for identity [default : 80]", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("perc_identity","pi");

        cmd.defineOption("threads", "count of threads for parallel execution [default : 1]", ArgvParser::OptionRequiresValue);
        cmd.defineOptionAlternative("threads","t");

        if(mode == "index")
        {
            cmd.defineOption("index", "output prefix for index", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
            cmd.defineOptionAlternative("index","i");
        }
        if(mode == "mapDirectly")
        {
            cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
            cmd.defineOptionAlternative("query","q");

            cmd.defineOption("all", "report all the mapping locations for a read, default is to consider few best ones");

			cmd.defineOption("output", "output file", ArgvParser::OptionRequiresValue);
			cmd.defineOptionAlternative("output","o");
        }
    }

    if(mode == "classify")
    {
        cmd.defineOption("DB", "Path to DB", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
        cmd.defineOption("mappings", "Path to mappings file", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
        cmd.defineOption("minreads", "Minimum number of reads per contig to be considered for fitting identity and length for the 'Unknown' functionality", ArgvParser::OptionRequiresValue);		
    }
	
    if(mode == "classifyU")
    {
        cmd.defineOption("DB", "Path to DB", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
        cmd.defineOption("mappings", "Path to mappings file", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
        cmd.defineOption("minreads", "Minimum number of reads per contig to be considered for fitting identity and length for the 'Unknown' functionality", ArgvParser::OptionRequiresValue);
    }	
  }


  /**
   * @brief                   Parse the file which has list of reference or query files
   * @param[in]   fileToRead  File containing list of ref/query files 
   * @param[out]  fileList    List of files will be saved in this vector   
   */
  template <typename VEC>
    void parseFileList(std::string &fileToRead, VEC &fileList)
    {
      std::string line;

      std::ifstream in(fileToRead);

      if (in.fail())
      {
        std::cerr << "ERROR, skch::parseFileList, Could not open " << fileToRead << "\n";
        exit(1);
      }

      while (std::getline(in, line))
      {
        fileList.push_back(line);
      }
    }

  /**
   * @brief                     validate the reference and query file(s)
   * @param[in] querySequences  vector containing query file names
   * @param[in] refSequences    vector containing reference file names
   */
  template <typename VEC>
    void validateInputFiles(VEC &querySequences, VEC &refSequences)
    {
      //Open file one by one
      for(auto &e : querySequences)
      {
		std::vector<std::string> e_parts = split(e, ",");
		for(auto e_part : e_parts)
		{
			std::ifstream in(e_part);
			
			if (in.fail())
			{
			  std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e_part << "\n";
			  exit(1);
			}
		}
      }

      for(auto &e : refSequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << "\n";
          exit(1);
        }
      }
    }

  /**
   * @brief                   Print the parsed cmd line options
   * @param[in]  parameters   parameters parsed from command line
   */
  void printCmdOptions(skch::Parameters &parameters, std::string mode)
  {
	if(mode == "mapAgainstIndex")
	{
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
		std::cout << "Index = " << parameters.index << std::endl;
		std::cout << "Query = " << parameters.querySequences << std::endl;
		std::cout << "Report all = " << parameters.reportAll << std::endl;
		std::cout << "Threads = " << parameters.threads << std::endl;
		std::cout << "Mapping output file = " << parameters.outFileName << std::endl;
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
	}
	else if(mode == "index")
	{
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
		std::cout << "Reference = " << parameters.refSequences << std::endl;
		std::cout << "Kmer size = " << parameters.kmerSize << std::endl;
		std::cout << "Window size = " << parameters.windowSize << std::endl;
		std::cout << "Min. read length >= " << parameters.minReadLength << std::endl;
		std::cout << "Alphabet = " << (parameters.alphabetSize == 4 ? "DNA" : "AA") << std::endl;
		std::cout << "P-value = " << parameters.p_value << std::endl;
		std::cout << "Percentage identity threshold = " << parameters.percentageIdentity << std::endl;
		std::cout << "Target max. memory = " << parameters.maximumMemory << std::endl;
		std::cout << "Index output prefix = " << parameters.index << std::endl;
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
	}
	else if(mode == "mapDirectly")
	{
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
		std::cout << "Reference = " << parameters.refSequences << std::endl;
		std::cout << "Query = " << parameters.querySequences << std::endl;
		std::cout << "Kmer size = " << parameters.kmerSize << std::endl;
		std::cout << "Window size = " << parameters.windowSize << std::endl;
		std::cout << "Min. read length >= " << parameters.minReadLength << std::endl;
		std::cout << "Alphabet = " << (parameters.alphabetSize == 4 ? "DNA" : "AA") << std::endl;
		std::cout << "P-value = " << parameters.p_value << std::endl;
		std::cout << "Percentage identity threshold = " << parameters.percentageIdentity << std::endl;
		std::cout << "Report all = " << parameters.reportAll << std::endl;
		std::cout << "Target max. memory = " << parameters.maximumMemory << std::endl;
		std::cout << "Threads = " << parameters.threads << std::endl;
		std::cout << "Mapping output file = " << parameters.outFileName << std::endl;
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
	}
	else if(mode == "classify")
	{
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
		std::cout << "DB = " << parameters.DB << std::endl;
		std::cout << "Mappings = " << parameters.mappingsForClassification << std::endl;
		std::cout << "Min. reads for 'U' = " << parameters.minimumReadsForU << std::endl;
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
	}
	else if(mode == "classifyU")
	{
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
		std::cout << "DB = " << parameters.DB << std::endl;
		std::cout << "Mappings = " << parameters.mappingsForClassification << std::endl;
		std::cout << "Min. reads for 'U' = " << parameters.minimumReadsForU << std::endl;
		std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
	}	
	else
	{
		throw std::runtime_error("Invalid code path");
	}
  }

  /**
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  sketch parameters are saved here
   */
void parseandSave(int argc, char** argv,  CommandLineProcessing::ArgvParser &cmd, skch::Parameters &parameters, std::string mode)
{
	int result = cmd.parse(argc, argv);

	//Make sure we get the right command line args
	if (result != ArgvParser::NoParserError)
	{
		std::cout << cmd.parseErrorDescription(result) << "\n";
		exit(1);
	}

	if((mode == "mapAgainstIndex") || (mode == "index"))
	{
		if(!cmd.foundOption("index"))
		{
			std::cerr << "Please provide index" << std::endl;
			exit(1);
		}
		else
		{
			std::stringstream str;
			std::string index;
			str << cmd.optionValue("index");
			str >> index;
			parameters.index = index;
		}
	}

	if((mode == "mapDirectly") || (mode == "index"))
	{
		if (!cmd.foundOption("reference"))
		{
			std::cerr << "Provide reference file (s)\n";
			exit(1);
		}
		else
		{
			std::stringstream str;
			std::string ref;
			str << cmd.optionValue("reference");
			str >> ref;
			parameters.refSequences.push_back(ref);
		}

		parameters.referenceSize = skch::CommonFunc::getReferenceSize(parameters.refSequences);
		parameters.alphabetSize = 4;

		if(cmd.foundOption("maxmemory"))
		{
			std::stringstream str;
			str << cmd.optionValue("maxmemory");
			size_t mm;
			str >> mm;
			parameters.maximumMemory = pow(1024,3) * mm;
		}
		else
		{
			parameters.maximumMemory = 0;
		}

		if(cmd.foundOption("kmer"))
		{
			std::stringstream str;
			str << cmd.optionValue("kmer");
			str >> parameters.kmerSize;
		}
		else
		{
			if(parameters.alphabetSize == 4)
				parameters.kmerSize = 16;
			else
				parameters.kmerSize = 5;
		}

		if(cmd.foundOption("pval"))
		{
			std::stringstream str;
			str << cmd.optionValue("pval");
			str >> parameters.p_value;
		}
		else
		{
			parameters.p_value = 1e-03;
		}

		if(cmd.foundOption("minReadLen"))
		{
			std::stringstream str;
			str << cmd.optionValue("minReadLen");
			str >> parameters.minReadLength;
		}
		else
		{
			parameters.minReadLength = 2000;
		}

		if(cmd.foundOption("perc_identity"))
		{
			std::stringstream str;
			str << cmd.optionValue("perc_identity");
			str >> parameters.percentageIdentity;
		}
		else
		{
			parameters.percentageIdentity = 80;
		}


		if(cmd.foundOption("window"))
		{
			std::stringstream str;
			str << cmd.optionValue("window");
			str >> parameters.windowSize;

			//Re-estimate p value
			int s = parameters.minReadLength * 2 / parameters.windowSize;
			parameters.p_value = skch::Stat::estimate_pvalue (s, parameters.kmerSize, parameters.alphabetSize,
			parameters.percentageIdentity,
			parameters.minReadLength, parameters.referenceSize);
		}
		else
		{
			//Compute optimal window size
			parameters.windowSize = skch::Stat::recommendedWindowSize(parameters.p_value,
			parameters.kmerSize, parameters.alphabetSize,
			parameters.percentageIdentity,
			parameters.minReadLength, parameters.referenceSize);
		}
	}

	if((mode == "mapDirectly") || (mode == "mapAgainstIndex"))
	{
		if(!cmd.foundOption("query"))
		{
			std::cerr << "Provide query file (s)\n";
			exit(1);
		}
		else
		{
			std::stringstream str;
			std::string query;
			str << cmd.optionValue("query");
			str >> query;
			parameters.querySequences.push_back(query);
		}

		if(cmd.foundOption("all"))
		{
			parameters.reportAll = true;
		}
		else
		{
			parameters.reportAll = false;
		}

    if(cmd.foundOption("threads"))
    {
			std::stringstream str;
			str << cmd.optionValue("threads");
			str >> parameters.threads;
    }
    else
    {
      parameters.threads = 1;
    }

		if(!cmd.foundOption("output"))
		{
			std::cerr << "Provide output file\n";
			exit(1);
		}
		else
		{
			std::stringstream str;
			std::string output;
			str << cmd.optionValue("output");
			str >> output;
			parameters.outFileName = output;
		}
	}

	if(mode == "classify")
	{
		if(!cmd.foundOption("DB"))
		{
			std::cerr << "Provide path to DB.\n";
			exit(1);
		}
		else
		{
			std::stringstream str;
			str << cmd.optionValue("DB");
			str >> parameters.DB;
		}

		if(!cmd.foundOption("mappings"))
		{
			std::cerr << "Provide path to mappings.\n";
			exit(1);
		}
		else
		{
			std::stringstream str;
			str << cmd.optionValue("mappings");
			str >> parameters.mappingsForClassification;
		}

		if(!cmd.foundOption("minreads"))
		{
			parameters.minimumReadsForU = 10000;
		}
		else
		{
			std::stringstream str;
			str << cmd.optionValue("minreads");
			str >> parameters.minimumReadsForU;
		}
	}
	

	if(mode == "classifyU")
	{
		if(!cmd.foundOption("DB"))
		{
			std::cerr << "Provide path to DB.\n";
			exit(1);
		}
		else
		{
			std::stringstream str;
			str << cmd.optionValue("DB");
			str >> parameters.DB;
		}

		if(!cmd.foundOption("mappings"))
		{
			std::cerr << "Provide path to mappings.\n";
			exit(1);
		}
		else
		{
			std::stringstream str;
			str << cmd.optionValue("mappings");
			str >> parameters.mappingsForClassification;
		}

		if(!cmd.foundOption("minreads"))
		{
			parameters.minimumReadsForU = 10000;
		}
		else
		{
			std::stringstream str;
			str << cmd.optionValue("minreads");
			str >> parameters.minimumReadsForU;
		}
	}	

	printCmdOptions(parameters, mode);

	//Check if files are valid
	validateInputFiles(parameters.querySequences, parameters.refSequences);
}
}


#endif
