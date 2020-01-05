#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <x86intrin.h>    //AVX/SSE Extensions


//represents up to a length 32 strand of DNA
struct DNA4{
	uint32_t letters[4];
	uint32_t& operator[](size_t i){return letters[i];}
};


uint32_t similarity(DNA4 a, DNA4 b)
{
	return __builtin_popcount ((a[0]&b[0])+(a[1]&b[1])+(a[2]&b[2])+(a[3]&b[3]));
}

DNA4 createDNA4(const char *s, uint_fast8_t length)
{
	uint32_t As = 0;
	uint32_t Cs = 0;
	uint32_t Gs = 0;
	uint32_t Ts = 0;
	for(int i = 0; i < length; i++)
	{
		char c = s[i];
		As <<= 1;
		Cs <<= 1;
		Gs <<= 1;
		Ts <<= 1;
		if(c=='A'||c=='a')
			As |= 1;
		else if(c=='C'||c=='c')
			Cs |= 1;
		else if(c=='G'||c=='g')
			Gs |= 1;
		else if(c=='T'||c=='t')
			Ts |= 1;
	}
	return DNA4{As,Cs,Gs,Ts};
}

DNA4 createDNA4(const std::string &s, uint_fast8_t length)
{
	return createDNA4(s.data(),length);
}


std::vector<std::vector<std::string>> match(const char * sequenceString, size_t sequenceLength, const std::vector<std::string> &targetStrings, uint_fast8_t targetLengths, uint_fast8_t minimumMatches)
{
	std::vector<DNA4> targets;
	std::vector<std::vector<std::string>> matches;
	for(std::string s: targetStrings)
	{
		targets.push_back(createDNA4(s,targetLengths));
		matches.push_back(std::vector<std::string>());
	}
	
	DNA4 sequence = createDNA4(sequenceString,targetLengths-1);
	size_t currentPos = targetLengths-1;
	int *matched1 = new int[targetStrings.size()/2+1];
	int *matched2 = new int[targetStrings.size()/2+1];
	while(currentPos < sequenceLength)
	{
		int numberOfMatches1 = 0;
		int numberOfMatches2 = 0;
		//add next character
		char c = sequenceString[currentPos];
		sequence[0] <<= 1;
		sequence[1] <<= 1;
		sequence[2] <<= 1;
		sequence[3] <<= 1;
		if(c=='A'||c=='a')
			sequence[0] |= 1;
		else if(c=='C'||c=='c')
			sequence[1] |= 1;
		else if(c=='G'||c=='g')
			sequence[2] |= 1;
		else if(c=='T'||c=='t')
			sequence[3] |= 1;
		
		//compare all the targets with this sequence
		for(int i = 0; i < targets.size();i+=2)
		{
			if(similarity(sequence,targets[i])>=minimumMatches)
			{
				matched1[numberOfMatches1] = i;
				numberOfMatches1++;
			}
			if(similarity(sequence,targets[i+1])>=minimumMatches)
			{
				matched2[numberOfMatches2] = i;
				numberOfMatches2++;
			}
		}
		
		for(int i = 0; i < numberOfMatches1;++i)
		{
			matches[matched1[i]].push_back("");//std::string(sequenceString + currentPos - targetLengths + 1,targetLengths));
		}
		
		for(int i = 0; i < numberOfMatches2;++i)
		{
			matches[matched2[i]].push_back("");//std::string(sequenceString + currentPos - targetLengths + 1,targetLengths));
		}
		currentPos++;
	}
	return matches;
}

char randNucleotide()
{
	constexpr char DNA[4] = {'A','C','G','T'};
	return DNA[rand()%4];
}

int main()
{
	constexpr size_t seqLength = 1000000;
	constexpr size_t numTargets = 1000;
	constexpr size_t targetLength = 23;
	constexpr size_t minSimilarity = 15;
	
	char * sequence = new char[seqLength];
	size_t As = 0;
	for(int i = 0; i < seqLength; ++i)
	{
		sequence[i] = randNucleotide();
	}
	std::vector<std::string> targetStrings;
	
	for(int i = 0; i < numTargets; i++)
	{
		std::string target;
		for(int j = 0; j < targetLength;++j)
			target += randNucleotide();
		targetStrings.push_back(target);
	}

	auto matches = match(sequence,seqLength,targetStrings,targetLength,minSimilarity);
	size_t totalMatches = 0;
	for(int i = 0; i < matches.size(); ++i)
	{
		/*std::cout << targetStrings[i] <<std::endl;
		for(int j = 0; j < matches[i].size();j++)
		{
			int similarity = 0;
			for(int k = 0; k < targetLength; k++)
				if(matches[i][j][k]==targetStrings[i][k])
					similarity++;
			std::cout << matches[i][j] <<"\t"<<similarity<< std::endl;
		}
		std::cout << std::endl;*/
		totalMatches += matches[i].size();
	}
	std::cout << totalMatches << std::endl;
	
	return 0;
}