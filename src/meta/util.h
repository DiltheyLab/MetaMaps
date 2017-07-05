/*
 * util.h
 *
 *  Created on: Jul 5, 2017
 *      Author: diltheyat
 */

#ifndef META_UTIL_H_
#define META_UTIL_H_

#include <string>
#include <vector>
#include <iostream>
#include <assert.h>

void errEx(const std::string& s)
{
	std::cerr << s << std::endl;
	exit(1);
}

void eraseNL(std::string& s)
{
	if (!s.empty() && s[s.length()-1] == '\n') {
	    s.erase(s.length()-1);
	}
}

std::vector<std::string> split(std::string input, std::string delimiter)
{
	std::vector<std::string> output;
	if(input.length() == 0)
	{
		return output;
	}

	if(delimiter == "")
	{
		output.reserve(input.size());
		for(unsigned int i = 0; i < input.length(); i++)
		{
			output.push_back(input.substr(i, 1));
		}
	}
	else
	{
		if(input.find(delimiter) == std::string::npos)
		{
			output.push_back(input);
		}
		else
		{
			int s = 0;
			int p = input.find(delimiter);

			do {
				output.push_back(input.substr(s, p - s));
				s = p + delimiter.size();
				p = input.find(delimiter, s);
			} while (p != (int)std::string::npos);
			output.push_back(input.substr(s));
		}
	}

	return output;
}



#endif /* META_UTIL_H_ */
