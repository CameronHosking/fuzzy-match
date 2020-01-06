#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>


//represents up to a length 32 strand of DNA
struct DNA4{
	uint32_t letters[4];
	const uint32_t& operator[](size_t i) const {return letters[i];}
	uint32_t& operator[](size_t i){return letters[i];}
};

constexpr char as[128] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
								
constexpr char cs[128] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
constexpr char gs[128] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
constexpr char ts[128] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0 };

void addCharacter(DNA4 &s, char c)
{
 	s[0] <<= 1;
	s[1] <<= 1;
	s[2] <<= 1;
	s[3] <<= 1;
	s[0] |= as[c];
	s[1] |= cs[c];
	s[2] |= gs[c];
	s[3] |= ts[c];
}

uint32_t similarity(const DNA4 &a,const DNA4 &b)
{
	return __builtin_popcount ((a[0]&b[0])+(a[1]&b[1])+(a[2]&b[2])+(a[3]&b[3]));
}

DNA4 createDNA4(const char *s, uint_fast8_t length)
{
	DNA4 ret{0,0,0,0};
	for(int i = 0; i < length; i++)
	{
		addCharacter(ret,s[i]);
	}
	return ret;
}

DNA4 createDNA4(const std::string &s, uint_fast8_t length)
{
	return createDNA4(s.data(),length);
}

std::vector<DNA4> stringsToDNA4(const std::vector<std::string> &targetStrings)
{
	std::vector<DNA4> targets;
	for(std::string s: targetStrings)
	{
		targets.push_back(createDNA4(s,s.size()));
	}
	return targets;
}
std::vector<std::vector<std::string>> match(const char * sequenceString, size_t sequenceLength, const std::vector<std::string> &targetStrings, uint_fast8_t targetLengths, uint_fast8_t mismatches)
{
	std::vector<DNA4> targets = stringsToDNA4(targetStrings);
	auto matches = std::vector<std::vector<std::string>>(targets.size());
	uint_fast8_t minimumMatches = targetLengths - mismatches;
	
	DNA4 sequence = createDNA4(sequenceString,targetLengths-1);
	size_t currentPos = targetLengths-1;
	int *matched1 = new int[targetStrings.size()/2+1];
	int *matched2 = new int[targetStrings.size()/2+1];
	while(currentPos < sequenceLength)
	{
		int numberOfMatches1 = 0;
		int numberOfMatches2 = 0;
		//add next character
		addCharacter(sequence,sequenceString[currentPos]);
		
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
				matched2[numberOfMatches2] = i+1;
				numberOfMatches2++;
			}
		}
		if(targets.size()%2)
		{
			if(similarity(sequence,targets[targets.size()-1])>=minimumMatches)
			{
				matched1[numberOfMatches1] = targets.size()-1;
				numberOfMatches1++;
			}
		}
		for(int i = 0; i < numberOfMatches1;++i)
		{
			matches[matched1[i]].push_back(std::string(sequenceString + currentPos - targetLengths + 1,targetLengths));
		}
		
		for(int i = 0; i < numberOfMatches2;++i)
		{
			matches[matched2[i]].push_back(std::string(sequenceString + currentPos - targetLengths + 1,targetLengths));
		}
		currentPos++;
	}
	return matches;
}

std::vector<std::vector<std::string>> match(const char * sequenceString, size_t sequenceLength, const std::vector<std::vector<std::vector<std::vector<std::pair<uint32_t,DNA4> > > > > &brackets, uint_fast8_t targetLength, uint32_t numTargets, uint_fast8_t mismatches)
{
	
	auto matches = std::vector<std::vector<std::string>>(numTargets);
	uint_fast8_t minimumMatches = targetLength - mismatches;
	
	DNA4 sequence = createDNA4(sequenceString,targetLength-1);
	size_t currentPos = targetLength-1;
	int *matched1 = new int[numTargets/2+1];
	int *matched2 = new int[numTargets/2+1];
	DNA4 As = createDNA4(std::string(targetLength,'A'),targetLength);
	DNA4 Cs = createDNA4(std::string(targetLength,'C'),targetLength);
	DNA4 Gs = createDNA4(std::string(targetLength,'G'),targetLength);
	while(currentPos < sequenceLength)
	{
		//std::cout << currentPos << std::endl;
		int numberOfMatches1 = 0;
		int numberOfMatches2 = 0;
		//add next character
		addCharacter(sequence,sequenceString[currentPos]);
		int numAs = similarity(As,sequence);
		int numCs = similarity(Cs,sequence);
		int numGs = similarity(Gs,sequence);
		//std::cout << numAs << "\t" << numCs << "\t" << numGs << std::endl;
		auto &targets = brackets[numAs][numCs][numGs];
		if(targets.size()==0)
		{
			currentPos++;
			continue;
		}
		//compare all the targets with this sequence
		for(int i = 0; i < targets.size()-1;i+=2)
		{
			//std::cout << targets.size()-1 << std::endl;
			if(similarity(sequence,targets[i].second)>=minimumMatches)
			{
				matched1[numberOfMatches1] = targets[i].first;
				numberOfMatches1++;
			}
			if(similarity(sequence,targets[i+1].second)>=minimumMatches)
			{
				matched2[numberOfMatches2] = targets[i+1].first;
				numberOfMatches2++;
			}

		}
		if(targets.size()%2)
		{
			if(similarity(sequence,targets[targets.size()-1].second)>=minimumMatches)
			{
				matched1[numberOfMatches1] = targets[targets.size()-1].first;
				numberOfMatches1++;
			}
		}
			
		for(int i = 0; i < numberOfMatches1;++i)
		{
			matches[matched1[i]].push_back(std::string(sequenceString + currentPos - targetLength + 1,targetLength));
		}
		for(int i = 0; i < numberOfMatches2;++i)
		{
			matches[matched2[i]].push_back(std::string(sequenceString + currentPos - targetLength + 1,targetLength));
		}
		currentPos++;
	}
	return matches;
}

std::vector<std::vector<std::vector<std::vector<std::pair<uint32_t,DNA4> > > > > getMatchBrackets(const std::vector<DNA4> &targets, int targetLength,int mismatches)
{
	std::vector<std::vector<std::vector<std::vector<std::pair<uint32_t,DNA4> > > > > brackets;
	for(int As = 0; As <= targetLength;As++)
	{
		brackets.push_back(std::vector<std::vector<std::vector<std::pair<uint32_t,DNA4> > > >());
		for(int Cs = 0; Cs <= targetLength - As;Cs++)
		{
			brackets[As].push_back(std::vector<std::vector<std::pair<uint32_t,DNA4> > >());
			for(int Gs = 0; Gs <= targetLength - As - Cs; Gs++)
			{
				brackets[As][Cs].push_back(std::vector<std::pair<uint32_t,DNA4> >());
				//Ts is no greater than targetLength - As - Cs - Gs
			}
		}
	}
	
	DNA4 As = createDNA4(std::string(targetLength,'A'),targetLength);
	DNA4 Cs = createDNA4(std::string(targetLength,'C'),targetLength);
	DNA4 Gs = createDNA4(std::string(targetLength,'G'),targetLength);
	
	for(int i = 0; i < targets.size();++i)
	{
		const DNA4 &target = targets[i];
		int numAs = similarity(As,target);
		int numCs = similarity(Cs,target);
		int numGs = similarity(Gs,target);
		
		for(int aOffset = -mismatches; aOffset <= mismatches; aOffset++)
		{
			if((numAs + aOffset < 0) | (numAs + aOffset > targetLength)) continue;
			for(int cOffset = -mismatches - ((aOffset < 0)?aOffset:0); cOffset <= mismatches - ((aOffset > 0)?aOffset:0); cOffset++)
			{
				if((numCs + cOffset < 0) | (numCs + cOffset + numAs + aOffset> targetLength)) continue;
				for(int gOffset = - mismatches - ((aOffset < 0)?aOffset:0) - ((cOffset < 0)?cOffset:0); gOffset <=  mismatches - ((aOffset > 0)?aOffset:0) - ((cOffset > 0)?cOffset:0); gOffset++)
				{
					if((numGs + gOffset < 0) | (numGs + gOffset + numCs + cOffset + numAs + aOffset> targetLength)) continue;
					//std::cout << numAs+aOffset << "\t" << numCs+cOffset << "\t" << numGs+gOffset << std::endl;
					brackets[numAs+aOffset][numCs+cOffset][numGs+gOffset].push_back(std::pair<uint32_t,DNA4>(i,target));
				}
			}
		}
	}
	/*
	for(int As = 0; As <= targetLength;As++)
	{
		for(int Cs = 0; Cs <= targetLength - As;Cs++)
		{
			for(int Gs = 0; Gs <= targetLength - As - Cs; Gs++)
			{
				if(brackets[As][Cs][Gs].size() > 0)
				{
					std::cout << As << "\t" << Cs << "\t" << Gs << "\t" << targetLength - As - Cs - Gs <<std::endl; 
					std::cout<<brackets[As][Cs][Gs].size() << std::endl;
				}
				//Ts is no greater than targetLength - As - Cs - Gs
			}
		}
	}*/
	std::cout << "test" << std::endl;
	return brackets;
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
	constexpr size_t mismatches = 8;
	
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

	auto targetBrackets = getMatchBrackets(stringsToDNA4(targetStrings),targetLength,mismatches);
	auto matches = match(sequence,seqLength,targetBrackets,targetLength,numTargets,mismatches);

	//auto matches = match(sequence,seqLength,targetStrings,targetLength,mismatches);

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