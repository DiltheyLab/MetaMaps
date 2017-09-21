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

template<typename T>
void printSorted(const std::map<std::string, T>& in, std::string title = "")
{
	std::vector<std::string> keys;
	for(const auto& e : in)
	{
		keys.push_back(e.first);
	}
	
	std::sort(keys.begin(), keys.end(), [&](const std::string& a, const std::string& b){
		return (in.at(a) < in.at(b));
	});
	
	std::reverse(keys.begin(), keys.end());
	
	if(title.length())
		std::cout << title << "\n";
	
	for(std::string k : keys)
	{
		if(title.length())
			std::cout << "\t";
		
		std::cout << k << ": " << in.at(k) << "\n";
	}
	
	std::cout << std::flush;
	
}

void errEx(const std::string& s)
{
	std::cerr << s << std::endl;
	exit(1);
}

std::string join(std::vector<std::string> parts, std::string delim)
{
	if(parts.size() == 0)
		return "";

	std::string ret = parts.at(0);

	for(unsigned int i = 1; i < parts.size(); i++)
	{
		ret.append(delim);
		ret.append(parts.at(i));
	}

	return ret;
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


size_t overlap_eins_larger(size_t eins_left, size_t eins_right, size_t zwei_left, size_t zwei_right)
{
	assert(eins_left < eins_right);
	assert(zwei_left < zwei_right);
	size_t eins_length = eins_right - eins_left + 1;
	size_t zwei_length = zwei_right - zwei_left + 1;
	assert(eins_length > 0);
	assert(zwei_length > 0);
	assert(eins_length >= zwei_length);

	if((eins_left <= zwei_left) and (eins_right >= zwei_right))
	{
		return zwei_length;
	}
	else if((zwei_left >= eins_left) and (zwei_left <= eins_right))
	{
		assert(zwei_right >= eins_right);
		return (eins_right - zwei_left + 1);
	}
	else if((zwei_right >= eins_left) and (zwei_right <= eins_right))
	{
		assert(zwei_left <= eins_left);
		return (zwei_right - eins_left + 1);
	}
	else
	{
		return 0;
	}
}

size_t overlap(size_t eins_left, size_t eins_right, size_t zwei_left, size_t zwei_right)
{
	assert(eins_left < eins_right);
	assert(zwei_left < zwei_right);
	size_t eins_length = eins_right - eins_left + 1;
	size_t zwei_length = zwei_right - zwei_left + 1;
	assert(eins_length > 0);
	assert(zwei_length > 0);

	size_t O;
	if(eins_length > zwei_length)
	{
		O = overlap_eins_larger(eins_left, eins_right, zwei_left, zwei_right);
		assert(O >= 0);

	}
	else
	{
		O = overlap_eins_larger(zwei_left, zwei_right, eins_left, eins_right);
		assert(O >= 0);
	}

	return O;
}

bool fileExists(std::string filepath)
{
	std::ifstream ifile(filepath);
	if(ifile.good())
	{
		ifile.close();
		return true;
	}
	else
	{
		return false;
	}
}


#endif /* META_UTIL_H_ */
