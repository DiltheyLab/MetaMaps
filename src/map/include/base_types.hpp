/**
 * @file    base_types.hpp
 * @brief   Critical type defintions for mapping algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef BASE_TYPES_MAP_HPP 
#define BASE_TYPES_MAP_HPP

#include <tuple>
#include <boost/serialization/serialization.hpp>
#include <vector>

namespace boost {
namespace serialization {
class access;
}
}

namespace skch
{
  typedef uint32_t hash_t;    //hash type
  typedef int offset_t;       //position within sequence
  typedef int seqno_t;        //sequence counter in file
  typedef int32_t strand_t;   //sequence strand 

  //C++ timer
  typedef std::chrono::high_resolution_clock Time;

  //Information about each minimizer
  class MinimizerInfo
  {
  private:
	    friend class boost::serialization::access;


	    template<typename Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & hash;
			ar & seqId;
			ar & wpos;
			ar & strand;
		}
  public:
    hash_t hash;                              //hash value
    seqno_t seqId;                            //sequence or contig id
    offset_t wpos;                            //First (left-most) window position when the minimizer is saved
    strand_t strand;                          //strand information

    //Lexographical less than comparison
    bool operator <(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand) 
        < std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    //Lexographical equality comparison
    bool operator ==(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand) 
        == std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    bool operator !=(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand) 
        != std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    static bool equalityByHash(const MinimizerInfo& x, const MinimizerInfo& y) {
      return x.hash == y.hash;
    }

    static bool lessByHash(const MinimizerInfo& x, const MinimizerInfo& y) {
      return x.hash < y.hash;
    }

  };

  //Type for map value type used for
  //L1 stage lookup index
  class MinimizerMetaData
  {
  private:
	    friend class boost::serialization::access;


	    template<typename Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & seqId;
			ar & wpos;
			ar & strand;
		}
  public:

    seqno_t seqId;          //sequence or contig id
    offset_t wpos;          //window position (left-most window)
    strand_t strand;        //strand information

    bool operator <(const MinimizerMetaData& x) const {
      return std::tie(seqId, wpos, strand) 
        < std::tie(x.seqId, x.wpos, x.strand);
    }
  };

  typedef hash_t MinimizerMapKeyType;
  typedef std::vector<MinimizerMetaData> MinimizerMapValueType;

  //Metadata recording for contigs in the reference DB
  class ContigInfo
  {
  private:
	    friend class boost::serialization::access;

	    template<typename Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & name;
			ar & len;
		}
  public:
    std::string name;       //Name of the sequence
    offset_t len;           //Length of the sequence
  };

  //Label tags for strand information
  enum strnd : strand_t
  {
    FWD = 1,  
    REV = -1
  };  

  //Final mapping result
  struct MappingResult
  {
    offset_t queryLen;                //length of the query sequence
    offset_t refStartPos;             //start position of the mapping on reference
    offset_t refEndPos;               //end pos
    seqno_t refSeqId;                 //internal sequence id of the reference contig
    seqno_t querySeqId;               //internal sequence id of the query sequence
    float nucIdentity;                //calculated identity
    float nucIdentityUpperBound;      //upper bound on identity (90% C.I.)
    int sketchSize;                   //sketch size
    int conservedSketches;            //count of conserved sketches
    strand_t strand;                  //strand
    float mappedRegionComplexity;     //estimated entropy in the mapped region on reference
    std::string queryName;            //name of query sequence
  };

  typedef std::vector<MappingResult> MappingResultsVector_t;

  //Container to save copy of kseq object
  struct InputSeqContainer
  {
    seqno_t seqCounter;                 //sequence counter
    offset_t len;                       //sequence length             
    std::string seq;                    //sequence string
    std::string seqName;                //sequence id

    /*
     * @brief               constructor
     * @param[in] kseq_seq  complete read or reference sequence
     * @param[in] kseq_id   sequence id name
     * @param[in] len       length of sequence
     */
    InputSeqContainer(const char * kseq_seq, const char * kseq_id, offset_t len, seqno_t seqcount)
    {
      this->seq = std::string{kseq_seq, std::size_t(len)};
      this->seqName = std::string{kseq_id};
      this->len = len;
      this->seqCounter = seqcount; 
    }
  };

  //Output type of map function
  struct MapModuleOutput
  {
    MappingResultsVector_t readMappings;  //read mapping coordinates
    std::string qseqName;                 //query sequence id

    //Function to erase all output mappings
    void reset()
    {
      this->readMappings.clear();
    }
  };

  //Information about query sequence during L1/L2 mapping
  template <typename MinimizerVec>
    struct QueryMetaData
    {
      char *seq;                          //query sequence pointer 
      seqno_t seqCounter;                 //query sequence counter
      offset_t len;                       //length of this query sequence
      int sketchSize;                     //sketch size
      MinimizerVec minimizerTableQuery;   //Vector of minimizers in the query 
    };

}

#endif
