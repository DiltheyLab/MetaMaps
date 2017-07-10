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

#include "fEM.h"
#include "util.h"


namespace meta
{

void doU(std::string DBdir, std::string mappedFile)
{
	std::set<std::string> relevantTaxonIDs = getTaxonIDsFromMappingsFile(mappedFile);

	std::string fn_fittedLengthAndIdentities= mappedFile + ".lengthAndIdentitiesPerMappingUnit";
	if(! fileExists(fn_fittedLengthAndIdentities))
	{
		std::cerr << "\n\nERROR: File " << fn_fittedLengthAndIdentities << " not existing.\n\nThis file is generated automatically by the EM step. Run the EM step first.\n\n";
	}

	identityHistogram iH;
	treeAdjustedIdentities tAI;
	std::map<size_t, double> readLengthHistogram;


}

}
#endif /* META_FU_H_ */
