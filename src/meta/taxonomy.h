/*
 * taxonomy.h
 *
 *  Created on: Jul 5, 2017
 *      Author: diltheyat
 */

#ifndef META_TAXONOMY_H_
#define META_TAXONOMY_H_

#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <boost/regex.hpp>

#include "util.h"

namespace meta
{

class oneTaxonName {
public:
	std::string scientific_name;
	std::string genbank_common_name;

};
class oneTaxonNode {
public:
	std::string id;
	std::string parent_id;
	std::string rank;
	oneTaxonName name;
	std::set<std::string> children;
};

class taxonomy {
protected:

public:
	std::map<std::string, oneTaxonNode> T;

	bool knowNode(std::string nodeID) const
	{
		return(T.count(nodeID) > 0);
	}

	
	std::string getFirstNonXNode(std::string nodeID)
	{
		if(nodeID.find("x") == std::string::npos)
		{
				return nodeID;
		}
		else
		{
			std::string runningNodeID = nodeID;
			while(runningNodeID.find("x") != std::string::npos)
			{
				runningNodeID = T.at(runningNodeID).parent_id;
			}	
			assert(nodeID.find("x") == std::string::npos);
			return runningNodeID;
		}
	}
	
	std::map<std::string, std::string> getUpwardNodesByRanks(std::string nodeID, std::set<std::string> targetRanks = std::set<std::string>() ) const
	{
		std::map<std::string, std::string> forReturn;
		
		std::vector<std::string> uN = getUpwardNodes(nodeID);
		for(auto n : uN)
		{
			std::string rank = T.at(n).rank;
			
			if(targetRanks.size())
			{
				if(!targetRanks.count(rank))
					continue;
			}
			
			if(rank != "no rank")
			{
				if(forReturn.count(rank))
				{
					std::cerr << "Node " << nodeID << " has multiple entries for rank " << rank << std::endl;
					exit(1);
				}
				forReturn[rank] = n;
			}
		}
		
		
		for(auto tL : targetRanks)
		{
			if(forReturn.count(tL) == 0)
			{
				forReturn[tL] = "Undefined";
			}
		}
		
		return forReturn;
	}

	std::vector<std::string> getUpwardNodes(std::string nodeID, bool includeStartNode = true) const
	{
		assert(T.count(nodeID));
		std::vector<std::string> u;

		if(includeStartNode)
			u.push_back(nodeID);

		while(nodeID != "1")
		{
			nodeID = T.at(nodeID).parent_id;
			u.push_back(nodeID);
		}

		return u;
	}

	const oneTaxonNode& getNode(std::string nodeID) const
	{
		assert(T.count(nodeID));
		return T.at(nodeID);
	}

	taxonomy(std::string dir)
	{
		std::string fn_names = dir + "/names.dmp";
		std::string fn_nodes = dir + "/nodes.dmp";
		std::string fn_merged = dir + "/merged.dmp";

		std::string line;

		boost::regex re_space("\\s*\\|\\s*");

		std::map<std::string, oneTaxonName> names;
		{
			std::ifstream namesStream;
			namesStream.open(fn_names.c_str());
			if(! namesStream.is_open())
			{
				std::cerr << "Cannot open file " << fn_names << " -- is '" << dir << "' a valid NCBI taxonomy?" << std::endl;
				exit(1);
			}

			while(namesStream.good())
			{
				std::getline(namesStream, line);
				eraseNL(line);
				if(line.length() == 0)
				{
					continue;
				}
				
				line = boost::regex_replace(line, re_space, "|");				
				std::vector<std::string> line_components = split(line, "|");

				
				assert(line_components.at(0).length());
				assert(line_components.at(1).length());

				std::string id = line_components.at(0);
				std::string name = line_components.at(1);
				std::string type = line_components.at(3);
				
				if(type == "scientific name")
				{
					names[id].scientific_name = name;
				}
				else if(type == "genbank common name")
				{
					names[id].genbank_common_name = name;
				}
			}
		}

		{
			std::ifstream nodesStream;
			nodesStream.open(fn_nodes.c_str());
			if(! nodesStream.is_open())
			{
				std::cerr << "Cannot open file " << fn_nodes << " -- is '" << dir << "' a valid NCBI taxonomy?" << std::endl;
				exit(1);
			}

			while(nodesStream.good())
			{
				std::getline(nodesStream, line);
				eraseNL(line);
				if(line.length() == 0)
				{
					continue;
				}
				
				line = boost::regex_replace(line, re_space, "|");				
				std::vector<std::string> line_components = split(line, "|");

				assert(line_components.at(0).length());

				std::string id = line_components.at(0);
				std::string parent = line_components.at(1);
				std::string rank = line_components.at(2);
				
				assert(rank.length());

				if(names.count(id) == 0)
				{
					std::cerr << "No name for taxon ID " << id << " in taxonomy directory " << dir << std::endl;
					exit(1);
				}

				oneTaxonNode thisNode;
				thisNode.name = names.at(id);
				thisNode.id = id;
				thisNode.parent_id = parent;
				thisNode.rank = rank;

				assert(T.count(id) == 0);
				T[id] = thisNode;
			}

			for(const auto& n : T)
			{
				const oneTaxonNode& thisNode = n.second;
				assert(thisNode.id == n.first);
				if(thisNode.parent_id != "1")
				{
					assert(T.count(thisNode.parent_id));
					T.at(thisNode.parent_id).children.insert(thisNode.id);
				}
			}
		}
		
		std::cout << "Read taxonomy from " << dir << " -- have " << T.size() << " nodes." << std::endl;
	}
};

}
#endif /* META_TAXONOMY_H_ */
