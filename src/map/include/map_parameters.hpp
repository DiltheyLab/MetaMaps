/**
 * @file    map_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_CONFIG_HPP 
#define SKETCH_CONFIG_HPP

#include <vector>

//Switch to enable timing of L1 and L2 stages for each read
//Timings are reported in a file
#define ENABLE_TIME_PROFILE_L1_L2 0

#include <boost/serialization/vector.hpp>


namespace boost {
namespace serialization {
class access;
}
}

namespace skch
{
  /**
   * @brief   configuration parameters for building sketch
   *          expected to be initialized using command line arguments
   */


  class Parameters
  {
	private:
		friend class boost::serialization::access;

		template<typename Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & kmerSize;
			ar & windowSize;
			ar & minReadLength;
			ar & alphabetSize;
			ar & referenceSize;
			ar & percentageIdentity;
			ar & p_value;
			ar & threads;
			ar & refSequences;
			ar & querySequences;
			ar & outFileName;
			ar & reportAll;
			ar & index;
			ar & maximumMemory;
		}

  public:
    int kmerSize;                                     //kmer size for sketching
    int windowSize;                                   //window size used for sketching 
    int minReadLength;                                //minimum read length which code maps
    int alphabetSize;                                 //alphabet size
    uint64_t referenceSize;                           //Approximate reference size
    float percentageIdentity;                         //user defined threshold for good similarity
    double p_value;                                   //user defined threshold for p value
    int threads;                                      //execution thread count
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
    bool reportAll;                                   //Report all alignments if this is true
    //const static int L2slideJump = 5;                 //parameter to jump read during L2 evaluation, 
                                                      //  1 is the most accurate

    std::string index;
    size_t maximumMemory;
    std::string DB;
	std::string mappingsForClassification;
	size_t minimumReadsForU;
	
	Parameters()
	{
		kmerSize = 0;
		windowSize = 0;
		minReadLength = 0;
		alphabetSize = 0;
		referenceSize = 0;
		percentageIdentity = 0;
		p_value = 0;
		threads = 0;
		reportAll = false;
		maximumMemory = 0;
		minimumReadsForU = 0;
	}
  };
}

#endif
