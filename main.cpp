#include <iostream>
#include <numeric>
#include "fuzzyMatch.h"

unsigned char randNucleotide()
{
	constexpr unsigned char DNA[4] = {'A','C','G','T'};
	return DNA[rand()%4];
}

int main()
{
	constexpr size_t seqLength = 10000000;
	constexpr size_t numTargets = 10000;
	constexpr size_t targetLength = 23;
	constexpr size_t mismatches = 4;
	const std::string filter = "XXXXXXXXXXXXXXXXXXXX*GG";
	char * sequence = new char[seqLength];
	for(size_t i = 0; i < seqLength; ++i)
	{
		sequence[i] = randNucleotide();
	}
	std::vector<std::string> targetStrings;
	
	for(size_t i = 0; i < numTargets; i++)
	{
		std::string target;
		for(uint_fast8_t j = 0; j < targetLength;++j)
			target += randNucleotide();
		//targets will always have the filter
		for(uint_fast8_t j = 0; j < filter.size();j++)
		{
			if((filter[j]=='A')|(filter[j]=='C')|(filter[j]=='G')|(filter[j]=='T'))
			{
				target[targetLength-1-j]=filter[j];
			}
		}
		targetStrings.push_back(target);
	}
	std::vector<std::string> sequences;
	sequences.push_back(std::string(sequence));
	auto matches = match(sequences,targetStrings, targetLength,mismatches,filter);

	delete[] sequence;
	size_t totalMatches = 0;
	for(size_t i = 0; i < matches.size(); ++i)
	{
		// std::cout << targetStrings[i] <<std::endl;
		// for(int j = 0; j < matches[i].size();j++)
		// {
		// 	int similarity = 0;
		// 	for(int k = 0; k < targetLength; k++)
		// 		if(matches[i][j][k]==targetStrings[i][k])
		// 			similarity++;
		// 	std::cout << matches[i][j] <<"\t"<<similarity<< std::endl;
		// }
		// std::cout << std::endl;
		// if(matches1[i].size()!=matches[i].size())
		// {
		// 	std::cout << i << "\t" << matches1[i].size() << "\t" << matches[i].size() << std::endl;
		// 	std::cout << matches[i].at(0) << std::endl;
		// 	std::cout << targetStrings[i] << std::endl;
		// 	if(matches1[i].size()==1)
		// 	{
		// 		std::cout << targetStrings[i] << std::endl;
		// 		//std::cout << matches1[0] << std::endl;
		// 	}
		// }
	 	totalMatches += matches[i].size();
	}
	std::cout << totalMatches << std::endl;
	
	return 0;
}