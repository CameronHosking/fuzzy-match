#include <iostream>
#include <numeric>
#include "fuzzyMatch.h"

uint32_t trivialGetSimilarity(const char *a, const char* b)
{
    uint32_t i = 0;
    uint32_t similarity = 0;
    while(a[i])
    {
        if(a[i]==b[i])
        {
            similarity++;
        }
        i++;
    }
    return similarity;    
}

unsigned char randNucleotide()
{
	constexpr unsigned char DNA[4] = {'A','C','G','T'};
	return DNA[rand()%4];
}

int main()
{
	constexpr size_t seqLength = 1000000;
	constexpr size_t numTargets = 1000000;
	constexpr size_t targetLength = 23;
	DNA4::setLength(targetLength);
	constexpr size_t mismatches = 4;
	const std::string filterString = "x*XXXXXXX*XXXXXXXXXXGXG";
	char * sequence = new char[seqLength+1];
	for(size_t i = 0; i < seqLength; ++i)
	{
		sequence[i] = randNucleotide();
	}
	sequence[seqLength] = '\0';
	std::vector<std::string> targetStrings;
	
	for(size_t i = 0; i < numTargets; i++)
	{
		std::string target;
		for(uint_fast8_t j = 0; j < targetLength;++j)
			target += randNucleotide();
		//targets will always have the filter
		for(int j = (int)filterString.size()-1; j >= 0;j--)
		{
			if((filterString[j]=='A')|(filterString[j]=='C')|(filterString[j]=='G')|(filterString[j]=='T'))
			{
				target[j]=filterString[j];
			}
		}
		targetStrings.push_back(target);
		
	}
	stop_watch s;
	Filter filter(filterString);
	DNA4Set d(filter);
	DNA4Set newSet(filter);

	s.start();
	doForTargetsInSequence(sequence,filter,AddToSet(d));
	s.stop();
	std::cout << "inserted " << d.size() << " targets in "<< s << std::endl;
	s.start();
	doForTargetsInSequence(sequence,filter,AddToSetIfExistingInOtherSet(newSet,d));
	s.stop();
	std::cout << "found and inserted " << newSet.size() << " targets in "<< s <<  std::endl;
	s.start();
	uint32_t size = d.size();
	doForTargetsInSequence(sequence,filter,RemoveFromSet(d));
	s.stop();
	std::cout << "removed " << size - d.size() << " targets in "<< s << " " << d.size() << " remain"<< std::endl;

	std::vector<std::string> sequences;
	sequences.push_back(std::string(sequence));
	s.start();
	OffTargetFinder offTargetFinder(stringsToDNA4(targetStrings),mismatches,filter,1500000000);
	s.stop();
	std::cout << "finished indexing in " << s <<std::endl;
	s.start();
	doForTargetsInSequence(sequence,filter,FindIfOffTarget(offTargetFinder,0));
	s.stop();
	std::cout << "comparisons took " << s <<std::endl;
	auto &offTargets = offTargetFinder.getOffTargets();

	delete[] sequence;
	size_t totalMatches = 0;
	uint32_t mismatchFrequencies[mismatches+1];
	for(uint32_t i = 0; i < mismatches+1; ++i)
	{
		mismatchFrequencies[i] = 0;
	}
	for(size_t i = 0; i < offTargets.size(); ++i)
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
		for(uint m = 0; m <= mismatches; m++)
		{
			for(size_t j = 0; j < offTargets[i][m].size(); j++)
		{
				mismatchFrequencies[m]++;
			}
			totalMatches += offTargets[i][m].size();
		}

	}
	for(uint32_t i = 0; i < mismatches+1; ++i)
	{
		std::cout << mismatchFrequencies[i] << " offtargets with distance " << i << std::endl;
	}
	std::cout << totalMatches << std::endl;
	
	return 0;
}