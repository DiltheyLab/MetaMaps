/**
 * @file    winSketch.hpp
 * @brief   routines to index the reference 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WIN_SKETCH_HPP 
#define WIN_SKETCH_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <cmath>
#include <cassert>
#include <zlib.h>  
#include <fstream>

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/kseq.h"
#include "common/murmur3.h"
#include "common/prettyprint.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/unordered_map.hpp>

namespace boost {
namespace serialization {
class access;
}
}


KSEQ_INIT(gzFile, gzread)

namespace skch
{
  /**
   * @class     skch::Sketch
   * @brief     sketches and indexes the reference (subject sequence)
   * @details  
   *            1.  Minimizers are computed in streaming fashion
   *                Computing minimizers is using double ended queue which gives
   *                O(reference size) complexity
   *                Algorithm described here:
   *                https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html
   *
   *            2.  Index hashes into appropriate format to enable fast search at L1 mapping stage
   */

	/*
	class Sketch;

	class Map {
		Map(const skch::Parameters &p, const skch::Sketch &refsketch, PostProcessResultsFn_t f = nullptr);
	};

	*/

    class Sketch
    {
	private:
		friend class boost::serialization::access;

		template<typename Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & percentageThreshold;
			ar & freqThreshold;
			ar & metadata;
			ar & sequencesByFileInfo;
			ar & minimizerPosLookupIndex;
			ar & minimizerIndex;
			ar & minimizerFreqHistogram;
		}

      //private members
    
      //algorithm parameters
      // skch::Parameters &param;

      //Ignore top % most frequent minimizers while lookups
      float percentageThreshold = 0.001;

      //Minimizers that occur this or more times will be ignored (computed based on percentageThreshold)
      int freqThreshold = std::numeric_limits<int>::max();

      public:

      typedef std::vector< MinimizerInfo > MI_Type;
      using MIIter_t = MI_Type::const_iterator;

      //Keep sequence length, name that appear in the sequence (for printing the mappings later)
      std::vector< ContigInfo > metadata;

      /*
       * Keep the information of what sequences come from what file#
       * Example [a, b, c] implies 
       *  file 0 contains 0 .. a-1 sequences
       *  file 1 contains a .. b-1 
       *  file 2 contains b .. c-1
       */
      std::vector< int > sequencesByFileInfo;

      //Index for fast seed lookup
      /*
       * [minimizer #1] -> [pos1, pos2, pos3 ...]
       * [minimizer #2] -> [pos1, pos2...]
       * ...
       */
      using MI_Map_t = std::unordered_map< MinimizerMapKeyType, MinimizerMapValueType >;
      MI_Map_t minimizerPosLookupIndex;

      private:

      /**
       * Keep list of minimizers, sequence# , their position within seq , here while parsing sequence 
       * Note : position is local within each contig
       * Hashes saved here are non-unique, ordered as they appear in the reference
       */
      MI_Type minimizerIndex;

      //Frequency histogram of minimizers
      //[... ,x -> y, ...] implies y number of minimizers occur x times
      std::map<int, int> minimizerFreqHistogram;

      public:

      /**
       * @brief   constructor
       *          also builds, indexes the minimizer table
       */

      //Make the default constructor private, non-accessible
      Sketch()
      {

      }

      /*
      Sketch(skch::Parameters &p)
      {
           this->build(p);
           this->index();
           this->computeFreqHist();
      }
	  */

      Sketch(skch::Parameters &p, size_t maximumMemory, std::function<void(Sketch*, size_t)>* oneIndexAction = nullptr)
      {
    	  assert(p.outFileName.length() == 0);
          build_and_store_index(p, maximumMemory, oneIndexAction);
      }

      private:

      size_t getMemoryOf(size_t hashes, size_t minimizers)
      {
    	  size_t estimated_number_buckets = hashes / 10;
    	  size_t memory_hash_table =
    			  estimated_number_buckets*(sizeof(size_t) + sizeof(void*)) + // pointers into hash table buckets
    			  hashes * sizeof(void*) + // from one bucket element to the next
				  hashes * sizeof(MinimizerMapValueType) + // each bucket element is a vector, and has a size
				  minimizers * sizeof(MinimizerMetaData); // total number of vector elements in all buckets
    	  memory_hash_table *= 1.2;

    	  size_t memory_vector = sizeof(MI_Type) + minimizers * sizeof(MinimizerInfo);

    	  return memory_hash_table + memory_vector;
      }

      void build_and_store_index(skch::Parameters &param, size_t maximumMemory, std::function<void(Sketch*, size_t)>* oneIndexAction)
      {
    	  std::cout << "Parameters used:\n";
    	  std::cout << "\t" << "- alphabetSize: " << param.alphabetSize << "\n";
    	  std::cout << "\t" << "- kmerSize: " << param.kmerSize << "\n";
    	  std::cout << "\t" << "- minReadLength: " << param.minReadLength << "\n";
    	  std::cout << "\t" << "- p_value: " << param.p_value << "\n";
    	  std::cout << "\t" << "- percentageIdentity: " << param.percentageIdentity << "\n";
    	  std::cout << "\t" << "- windowSize: " << param.windowSize << "\n";
    	  std::cout << "\t" << "- maximumMemory: ~" << maximumMemory/std::pow(1024,3) << " GB\n";
    	  std::cout << "\n" << std::flush;

          seqno_t seqCounter = 0;
          MI_Type thisContig_minimizerIndex;
          size_t runningHashes = 0;
          size_t runningMinimizers = 0;
    	  int runningIndexNumber = 1;

    	  size_t currentIndex_counter_sawSequences = 0;

    	  std::string outputFn = param.outFileName;

    	  std::vector<std::string> generatedFiles;

		  auto processCurrentState = [&](int N)
		  {
			  assert(metadata.size() == currentIndex_counter_sawSequences);
			  std::cerr << "\nCall storeCurrentState with " << currentIndex_counter_sawSequences << "\n\n" << std::flush;

			  this->computeFreqHist();

			  (*oneIndexAction)(this, N);

			  /*
			  if(mappingMode)
			  {
				  (*oneIndexAction)(this, N);
				  //std::cout << "Mapping for fragment " << N << std::endl;
				  //std::string thisN_outputFn = outputFn + "." + std::to_string(N);
				  //skch::Parameters param_thisN = param;
				  //param_thisN.outFileName = thisN_outputFn;
				  // skch::Map* mapper = new skch::Map(param_thisN, *(this));
				  //generatedFiles.push_back(thisN_outputFn);
			  }
			  else
			  {
				  std::string outputFn = prefix + "." + std::to_string(N);
				  std::ofstream ostream(outputFn.c_str());
				  if(! ostream.is_open())
				  {
					  throw std::runtime_error("Cannot open file " + outputFn + " for serialization.");
				  }
				  boost::archive::text_oarchive archive(ostream);
				  archive & (*this);
				  ostream.close();
				  std::cout << "Stored state in file " << outputFn << "\n" << std::flush;
				  generatedFiles.push_back(outputFn);
			  }

			  */
		  };

          for(const auto &fileName : param.refSequences)
          {
            //Open the file using kseq
            FILE *file = fopen(fileName.c_str(), "r");
            gzFile fp = gzdopen(fileno(file), "r");
            kseq_t *seq = kseq_init(fp);

            //size of sequence
            offset_t len;

            while ((len = kseq_read(seq)) >= 0)
            {

              // std::cerr << "Adding " << seq->name.s << " with length " << len << std::endl;

              //Is the sequence too short?
              if(len < param.windowSize || len < param.kmerSize)
              {
                metadata.push_back( ContigInfo{seq->name.s, (offset_t)seq->seq.l} );

                seqCounter++;
                currentIndex_counter_sawSequences++;
              }
              else
              {
                thisContig_minimizerIndex.clear();

                skch::CommonFunc::addMinimizers(thisContig_minimizerIndex, seq->seq.s, len, param.kmerSize, param.windowSize, param.alphabetSize, currentIndex_counter_sawSequences);

                size_t thisContig_wouldAdd_hashes = 0;
                size_t thisContig_wouldAdd_minimizers = thisContig_minimizerIndex.size();

                std::set<MinimizerMapKeyType> uniqueNovelHashes;
                for(const auto &e : thisContig_minimizerIndex)
                {
                	if((uniqueNovelHashes.count(e.hash) == 0) && (minimizerPosLookupIndex.count(e.hash) == 0))
                	{
                		thisContig_wouldAdd_hashes++;
                		uniqueNovelHashes.insert(e.hash);
                	}
                }

                size_t ifAdd_total_hashes = runningHashes + thisContig_wouldAdd_hashes;
                size_t ifAdd_total_minimizers = runningMinimizers + thisContig_wouldAdd_minimizers;

                size_t memory_after_add = getMemoryOf(ifAdd_total_hashes, ifAdd_total_minimizers);
				
				/*
                std::cerr << "thisContig_wouldAdd_hashes" << ": " << thisContig_wouldAdd_hashes << "\n";
                std::cerr << "thisContig_wouldAdd_minimizers" << ": " << thisContig_wouldAdd_minimizers << "\n";
                std::cerr << "ifAdd_total_hashes" << ": " << ifAdd_total_hashes << "\n";
                std::cerr << "ifAdd_total_minimizers" << ": " << ifAdd_total_minimizers << "\n";
                std::cerr << "memory_after_add" << ": " << memory_after_add << "\n";
                std::cerr << "maximumMemory" << ": " << maximumMemory << "\n\n";
				*/
				
                if((maximumMemory > 0) && (memory_after_add > maximumMemory))
                {
                	processCurrentState(runningIndexNumber);

                	minimizerIndex.clear();
                	minimizerPosLookupIndex.clear();
                	metadata.clear();

                	runningHashes = 0;
                	runningMinimizers = 0;
                	currentIndex_counter_sawSequences = 0;
                	runningIndexNumber++;

                    std::set<MinimizerMapKeyType> uniqueNovelHashes;
                    for(auto &e : thisContig_minimizerIndex)
                    {
                    	uniqueNovelHashes.insert(e.hash);
                    	e.seqId = 0;
                    }
                    thisContig_wouldAdd_hashes = uniqueNovelHashes.size();

                    ifAdd_total_hashes = runningHashes + thisContig_wouldAdd_hashes;
                    ifAdd_total_minimizers = runningMinimizers + thisContig_wouldAdd_minimizers;

                    memory_after_add = getMemoryOf(ifAdd_total_hashes, ifAdd_total_minimizers);
                    if(memory_after_add > maximumMemory)
                    {
                    	std::cerr << "Can't index file " << fileName << " within current memory limits - contig " << seq->name.s << " is too large. I estimate I'd need " << (memory_after_add/1024) << " kb, but I am instructed to only use " << (maximumMemory/1024) << " kb." << std::endl;
                    	throw std::runtime_error("Can't index file " + fileName + " within current memory limits - contig " + seq->name.s + " is too large");
                    }

                }

				for(const auto &e : thisContig_minimizerIndex)
				{
				  // [hash value -> info about minimizer]
				  minimizerPosLookupIndex[e.hash].push_back(
					  MinimizerMetaData{e.seqId, e.wpos, e.strand});
				}

				minimizerIndex.insert(minimizerIndex.end(), thisContig_minimizerIndex.begin(), thisContig_minimizerIndex.end());

				metadata.push_back( ContigInfo{seq->name.s, (offset_t)seq->seq.l} );

				runningHashes = ifAdd_total_hashes;
				runningMinimizers = ifAdd_total_minimizers;
				currentIndex_counter_sawSequences++;
				seqCounter++;
				
				std::cerr << "Added " << seq->name.s << " with length " << len << "; est. memory ~" << (memory_after_add/pow(1024,3)) << " GB" << std::endl;
              }
            }

            sequencesByFileInfo.push_back(seqCounter);

            kseq_destroy(seq);
            gzclose(fp); //close the file handler
            fclose(file);

          }

          processCurrentState(runningIndexNumber);


          std::cout << "INFO, skch::Sketch::build, minimizers picked from reference = " << minimizerIndex.size() << std::endl;

          std::cout << "\nIndex construction DONE, wrote " <<  runningIndexNumber << " files.\n\n" << std::flush;
      }

      /**
       * @brief     build the sketch table
       * @details   compute and save minimizers from the reference sequence(s)
       *            assuming a fixed window size
       */

      /*
      void build(skch::Parameters &param)
      {

        //sequence counter while parsing file
        seqno_t seqCounter = 0;

        for(const auto &fileName : param.refSequences)
        {

#ifdef DEBUG
        std::cout << "INFO, skch::Sketch::build, building minimizer index for " << fileName << std::endl;
#endif

          //Open the file using kseq
          FILE *file = fopen(fileName.c_str(), "r");
          gzFile fp = gzdopen(fileno(file), "r");
          kseq_t *seq = kseq_init(fp);


          //size of sequence
          offset_t len;

          while ((len = kseq_read(seq)) >= 0) 
          {
            //Save the sequence name
            metadata.push_back( ContigInfo{seq->name.s, (offset_t)seq->seq.l} );

            //Is the sequence too short?
            if(len < param.windowSize || len < param.kmerSize)
            {
#ifdef DEBUG
              std::cout << "WARNING, skch::Sketch::build, found an unusually short sequence relative to kmer and window size" << std::endl;
#endif
              seqCounter++;
              continue;  
            }
            else
            {
              skch::CommonFunc::addMinimizers(this->minimizerIndex, seq, param.kmerSize, param.windowSize, param.alphabetSize, seqCounter);
            }

            seqCounter++;
          }

          sequencesByFileInfo.push_back(seqCounter);

          kseq_destroy(seq);  
          gzclose(fp); //close the file handler 
          fclose(file);
        }

        std::cout << "INFO, skch::Sketch::build, minimizers picked from reference = " << minimizerIndex.size() << std::endl;

      }

	*/
      /**
       * @brief   build the index for fast lookups using minimizer table
       */
      /*
      void index()
      {
        //Parse all the minimizers and push into the map
        for(auto &e : minimizerIndex)
        {
          // [hash value -> info about minimizer]
          minimizerPosLookupIndex[e.hash].push_back( 
              MinimizerMetaData{e.seqId, e.wpos, e.strand});
        }

        std::cout << "INFO, skch::Sketch::index, unique minimizers = " << minimizerPosLookupIndex.size() << std::endl;
      }
	*/

      /**
       * @brief   report the frequency histogram of minimizers using position lookup index
       *          and compute which high frequency minimizers to ignore
       */
      void computeFreqHist()
      {
        //1. Compute histogram

    	if(minimizerPosLookupIndex.size() > 0)
    	{
        for(auto &e : this->minimizerPosLookupIndex)
          this->minimizerFreqHistogram[e.second.size()] += 1;

        std::cout << "INFO, skch::Sketch::computeFreqHist, Frequency histogram of minimizers = " <<  *this->minimizerFreqHistogram.begin() <<  " ... " << *this->minimizerFreqHistogram.rbegin() << std::endl;

        //2. Compute frequency threshold to ignore most frequent minimizers

        int64_t totalUniqueMinimizers = this->minimizerPosLookupIndex.size();
        int64_t minimizerToIgnore = totalUniqueMinimizers * percentageThreshold / 100;

        int64_t sum = 0;

        //Iterate from highest frequent minimizers
        for(auto it = this->minimizerFreqHistogram.rbegin(); it != this->minimizerFreqHistogram.rend(); it++)
        {
          sum += it->second; //add frequency
          if(sum < minimizerToIgnore)
          {
            this->freqThreshold = it->first;
            //continue
          }
          else if(sum == minimizerToIgnore)
          {
            this->freqThreshold = it->first;
            break;
          }
          else
          {
            break;
          }
        }

        if(this->freqThreshold != std::numeric_limits<int>::max())
          std::cout << "INFO, skch::Sketch::computeFreqHist, With threshold " << this->percentageThreshold << "\%, ignore minimizers occurring >= " << this->freqThreshold << " times during lookup." << std::endl;
        else
          std::cout << "INFO, skch::Sketch::computeFreqHist, With threshold " << this->percentageThreshold << "\%, consider all minimizers during lookup." << std::endl;
    	}
      }

      public:

      /**
       * @brief               search hash associated with given position inside the index
       * @details             if MIIter_t iter is returned, than *iter's wpos >= winpos
       * @param[in]   seqId
       * @param[in]   winpos
       * @return              iterator to the minimizer in the index
       */
      MIIter_t searchIndex(seqno_t seqId, offset_t winpos) const
      {
        std::pair<seqno_t, offset_t> searchPosInfo(seqId, winpos);

        /*
         * std::lower_bound --  Returns an iterator pointing to the first element in the range
         *                      that is not less than (i.e. greater or equal to) value.
         */
        MIIter_t iter = std::lower_bound(this->minimizerIndex.begin(), this->minimizerIndex.end(), searchPosInfo, cmp);

        return iter;
      }

      /**
       * @brief     Return end iterator on minimizerIndex
       */
      MIIter_t getMinimizerIndexEnd() const
      {
        return this->minimizerIndex.end();
      }

      int getFreqThreshold() const
      {
        return this->freqThreshold;
      }

      private:

      /**
       * @brief     functor for comparing minimizers by their position in minimizerIndex
       * @details   used for locating minimizers with the required positional information
       */
      struct compareMinimizersByPos
      {
        typedef std::pair<seqno_t, offset_t> P;

        bool operator() (const MinimizerInfo &m, const P &val)
        {
          return ( P(m.seqId, m.wpos) < val);
        }

        bool operator() (const P &val, const MinimizerInfo &m)
        {
          return (val < P(m.seqId, m.wpos) );
        }
      } cmp;

    }; //End of class Sketch
} //End of namespace skch

#endif
