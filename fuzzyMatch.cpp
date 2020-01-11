#include <cstdlib>
#include <string>
#include <vector>
#include <array>
#include <iostream>


//represents up to a length 32 strand of DNA
struct DNA4{
	uint32_t letters[4];
	const uint32_t& operator[](size_t i) const {return letters[i];}
	uint32_t& operator[](size_t i){return letters[i];}
	bool operator==(const DNA4 & o) const
	{
		return o[0]==letters[0]&o[1]==letters[1]&o[2]==letters[2]&o[3]==letters[3];
	}
	bool operator!=(const DNA4 &o)const{return !operator==(o);}
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

uint32_t getSimilarity(const DNA4 &a,const DNA4 &b)
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

std::string DNA4ToString(const DNA4 &a, uint_fast8_t length)
{
	std::string s(length,'_');
	for(uint_fast8_t i=0;i<length;i++)
	{
		if((a[0]>>i)&1)
		{
			s[length-i-1] = 'A';
		}
		else if((a[1]>>i)&1)
		{
			s[length-i-1] = 'C';
		}
		else if((a[2]>>i)&1)
		{
			s[length-i-1] = 'G';
		}
		else if((a[3]>>i)&1)
		{
			s[length-i-1] = 'T';
		}
	}
	return s;
}

std::vector<DNA4> stringsToDNA4(const std::vector<std::string> &targetStrings)
{
	std::vector<DNA4> targets;
	targets.reserve(targetStrings.size());
	for(std::string s: targetStrings)
	{
		targets.push_back(createDNA4(s,s.size()));
	}
	return targets;
}

//storage of target buckets
std::vector<std::pair<DNA4,uint32_t>> buckets[32*32*32];

inline std::vector<std::pair<DNA4,uint32_t>> &getTargets(uint32_t a, uint32_t c, uint32_t g)
{
	return buckets[a*32*32+c*32+g];
}

std::vector<std::vector<std::string>> match(const char * sequenceString, size_t sequenceLength, const std::vector<std::string> &targetStrings, uint_fast8_t targetLengths, uint_fast8_t mismatches)
{
	std::vector<DNA4> targets = stringsToDNA4(targetStrings);
	auto matches = std::vector<std::vector<std::string>>(targets.size());
	if(targets.size()==0)
	{
		return matches;
	}
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
		for(int i = 0; i < targets.size()-1;i+=2)
		{
			if(getSimilarity(sequence,targets[i])>=minimumMatches)
			{
				matched1[numberOfMatches1] = i;
				numberOfMatches1++;
			}
			if(getSimilarity(sequence,targets[i+1])>=minimumMatches)
			{
				matched2[numberOfMatches2] = i+1;
				numberOfMatches2++;
			}
		}
		if(targets.size()%2)
		{
			if(getSimilarity(sequence,targets[targets.size()-1])>=minimumMatches)
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
	delete[] matched1;
	delete[] matched2;
	return matches;
}

std::vector<std::vector<std::string>> match(const char * sequenceString, size_t sequenceLength, uint_fast8_t targetLength, const std::vector<std::pair<std::array<uint32_t,33>, std::vector<std::pair<DNA4,uint32_t>>>> &targetTargetSimilarities, uint_fast8_t mismatches)
{
	std::cout << "starting to match" << std::endl;
	auto matches = std::vector<std::vector<std::string>>(targetTargetSimilarities.size());
	uint_fast8_t minimumMatches = targetLength - mismatches;
	
	DNA4 sequence = createDNA4(sequenceString,targetLength-1);
	size_t currentPos = targetLength-1;
	int *matched1 = new int[targetTargetSimilarities.size()/2+1];
	int *matched2 = new int[targetTargetSimilarities.size()/2+1];
	bool specificTargetMatched = new bool[targetTargetSimilarities.size()];
	//masks off the the bits that aren't in the target length
	uint32_t mask = (~0u)>>(32-targetLength);
	while(currentPos < sequenceLength)
	{
		//std::cout << currentPos << std::endl;
		int numberOfMatches1 = 0;
		int numberOfMatches2 = 0;
		//add next character
		addCharacter(sequence,sequenceString[currentPos]);

		int numAs = __builtin_popcount(sequence[0]&mask);
		int numCs = __builtin_popcount(sequence[1]&mask);
		int numGs = __builtin_popcount(sequence[2]&mask);

		//std::cout << numAs << "\t" << numCs << "\t" << numGs << std::endl;
		const std::pair<DNA4,uint32_t> *target = buckets[numAs*32*32+numCs*32+numGs].data();
		const std::pair<DNA4,uint32_t> *end = target + buckets[numAs*32*32+numCs*32+numGs].size();	
		uint_fast8_t maxSimilarity = 11;
		if(target==end)
		{
			currentPos++;
			continue;
		}
		uint totalChecked = 0;
		//compare all the targets with this sequence
		while(target < end-1)
		{
			//std::cout << targets.size()-1 << std::endl;
			uint_fast8_t similarity1 = getSimilarity(sequence,target->first);
			uint_fast8_t similarity2 = getSimilarity(sequence,(target+1)->first);
			totalChecked+=2;
			if(similarity1>maxSimilarity||similarity2>maxSimilarity)
			{
				if(similarity2>similarity1)
				{
					target = target+1;
					similarity1 = similarity2;
				}
				maxSimilarity = similarity1;
				auto &t = targetTargetSimilarities[target->second];

				uint_fast8_t lowestSimilarity = similarity1<mismatches?0:similarity1-mismatches;
				uint_fast8_t highestSimilarity = similarity1+mismatches>targetLength?targetLength:similarity1+mismatches;
				target = t.second.data() + t.first[23-highestSimilarity];
				end = t.second.data() + t.first[23-lowestSimilarity+1];
				if(currentPos%10000==0)
					std::cout << (int)similarity1 << "\t" << end-target << "\t"<<totalChecked<<std::endl;
				numberOfMatches1 = 0;
				numberOfMatches2 = 0;
				continue;
			}

			if(similarity1>=minimumMatches)
			{
				matched1[numberOfMatches1] = target->second;
				numberOfMatches1++;
			}
			if(similarity2>=minimumMatches)
			{
				matched2[numberOfMatches2] = (target+1)->second;
				numberOfMatches2++;
			}
			
			target+=2;
		}
		if(currentPos%10000==0)
			std::cout << "\t\t" << totalChecked << std::endl;
		//with an odd number of targets we need to check the last value 
		if(target == end - 1)
		{
			if(getSimilarity(sequence,(end-1)->first)>=minimumMatches)
			{
				matched1[numberOfMatches1] = (end-1)->second;
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
	
	delete[] matched1;
	delete[] matched2;
	return matches;
}

//returns a vector which contains an array and a vector for each target. The array stores the start of each mismatch bucket in the vector.
//so to access targets with mismatch 10-16 to target i you would access v[i].second[v[i].first[10]] up to v[i].second[v[i].first[17]]
//or std::pair<DNA4,uint32_t> * it = v[i].second.data(); for(int i = v[i].first[10]; i < v[i].first[17]; i++){//do something with it[i]} 
std::vector<std::pair<std::array<uint32_t,33>, std::vector<std::pair<DNA4,uint32_t>>>> getTargetSimilarities(const std::vector<DNA4> &targets, uint_fast8_t maxMismatch)
{
	std::vector<std::array<std::vector<std::pair<DNA4,uint32_t>>,33>> mismatches(targets.size());
	uint_fast8_t targetLength = getSimilarity(targets[0],targets[0]);

	for(int i = 0; i < targets.size();i++)
	{
		//include the similarity of the target with itself to make later operaions more simple
		mismatches[i][0].push_back(std::pair<DNA4,uint32_t>(targets[i],i));

		//store the similarity of each target with this one
		for(int j = i+1; j < targets.size();j++)
		{
			uint_fast8_t similarity = getSimilarity(targets[i],targets[j]);
			//we only go into this array if similarity is at least 12 so the minimum similarity = 12 - maxMismatch
			if(similarity>=11-maxMismatch)
			{
				uint_fast8_t mismatch = targetLength - similarity;
				mismatches[i][mismatch].push_back(std::pair<DNA4,uint32_t>(targets[j],j));
				mismatches[j][mismatch].push_back(std::pair<DNA4,uint32_t>(targets[i],i));
			}
				
		}
		}

	std::vector<std::pair<std::array<uint32_t,33>, std::vector<std::pair<DNA4,uint32_t>>>> packedMismatches(targets.size());
	for(int i = 0; i < targets.size();i++)
	{
		uint32_t total = 0;
		auto &targetPair = packedMismatches[i];
		for(int mismatch = 0; mismatch <= 32; mismatch++)
		{
			targetPair.first[mismatch] = total;
			total += mismatches[i][mismatch].size();
		}
		targetPair.second.reserve(total);
		for(int j = 0; j <= targetLength; j++)
		{
			targetPair.second.insert(targetPair.second.end(),mismatches[i][j].begin(),mismatches[i][j].end());
		}
	}

	return packedMismatches;
}

void clearTargetBuckets()
{
	for(int i = 0; i < 32*32*32; i++)
	{
		buckets[i].clear();
	}
}

bool isNewTargetSet(const std::vector<DNA4> &targets, int targetLength,int mismatches)
{
	static struct BucketInput
	{
		std::vector<DNA4> targets;
		uint32_t targetLength;
		uint_fast8_t mismatches;
	}currentBucketInput;
	
	if(currentBucketInput.targets.size()==targets.size())
	{
		bool sameInput = true;
		for(int i = 0; i < currentBucketInput.targets.size(); i++)
		{
			if(currentBucketInput.targets[i]!=targets[i])
			{
				sameInput = false;
				break;
			}
		}
		sameInput &= currentBucketInput.targetLength == targetLength;
		sameInput &= currentBucketInput.mismatches == mismatches;
		if(sameInput)
			return false;
	}

	currentBucketInput.targets=targets;
	currentBucketInput.targetLength = targetLength;
	currentBucketInput.mismatches = mismatches;
	return true;
}

void fillTargetBuckets(const std::vector<DNA4> &targets, int targetLength,int mismatches)
{
	// std::vector<std::vector<std::vector<std::vector<std::pair<DNA4,uint32_t> > > > > brackets;
	// int i = 0;
	// for(int As = 0; As <= targetLength;As++)
	// {
	// 	//std::cout << As << "\t" << i << std::endl;
	// 	//length of A[n] is ((k-n+1)^2+(k-n+1))/2 or (n^2 -(3+2k)n + (k^2 + 3k + 2))/2 
	// 	//A[n] starts at A[0] + n(n^2 - (6+3k)n + 3k^2 + 12k + 11)/6 
	// 	//total length of A arrays is (k^3 + 6k^2 + 11k + 6)/6
	// 	//brackets.push_back(std::vector<std::vector<std::vector<std::pair<DNA4,uint32_t> > > >());
	// 	for(int Cs = 0; Cs <= targetLength - As;Cs++)
	// 	{
	// 		//length of C[n] is k-A-n+1
	// 		//C[n] starts at C[0] + (k-A+1)n - (n-1)n/2
	// 		//total length of C arrays is ((k-A+1)^2+(k-A+1))/2
	// 		//brackets[As].push_back(std::vector<std::vector<std::pair<DNA4,uint32_t> > >());
	// 		//length of G[n] is 1
	// 		//G[n] starts at G[0] + n
	// 		//total length of G arrays is k-A-C+1
	// 		for(int Gs = 0; Gs <= targetLength - As - Cs; Gs++)
	// 		{

	// 			//brackets[As][Cs].push_back(std::vector<std::pair<DNA4,uint32_t> >());
	// 			//buckets[32*32*As+32*Cs+Gs].push_back()
	// 			//Ts is no greater than targetLength - As - Cs - Gs
	// 			i++;
	// 		}
	// 	}
	// }
	//std::cout << targetLength << "\t" << i << std::endl;
	
	if(isNewTargetSet(targets,targetLength,mismatches))
	{
		clearTargetBuckets();
	}

	uint32_t mask = (~0u)>>(32-targetLength);
	for(int i = 0; i < targets.size();++i)
	{
		const DNA4 &target = targets[i];
		int numAs = __builtin_popcount(target[0]&mask);
		int numCs = __builtin_popcount(target[1]&mask);
		int numGs = __builtin_popcount(target[2]&mask);
		
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
					//brackets[numAs+aOffset][numCs+cOffset][numGs+gOffset].push_back(std::pair<DNA4,uint32_t>(i,target));
					buckets[32*32*(numAs+aOffset)+32*(numCs+cOffset)+(numGs+gOffset)].push_back(std::pair<DNA4,uint32_t>(target,i));
				}
			}
		}
	}
	double totalWork;
	for(int As = 0; As <= targetLength;As++)
	{
		for(int Cs = 0; Cs <= targetLength - As;Cs++)
		{
			for(int Gs = 0; Gs <= targetLength - As - Cs; Gs++)
			{
				totalWork += buckets[32*32*As+32*Cs+Gs].size()*buckets[32*32*As+32*Cs+Gs].size();
				if(buckets[32*32*As+32*Cs+Gs].size() > 0)
				{
					//std::cout << As << "\t" << Cs << "\t" << Gs << "\t" << targetLength - As - Cs - Gs <<std::endl; 
					//std::cout<<brackets[As][Cs][Gs].size() << std::endl;

				}
				//Ts is no greater than targetLength - As - Cs - Gs
			}
		}
	}
	std::cout << totalWork << std::endl;
}

char randNucleotide()
{
	constexpr char DNA[4] = {'A','C','G','T'};
	return DNA[rand()%4];
}

int main()
{
	constexpr size_t seqLength = 1000000;
	constexpr size_t numTargets = 10000;
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
	fillTargetBuckets(stringsToDNA4(targetStrings),targetLength,mismatches);
	auto matches = match(sequence,seqLength,targetLength, getTargetSimilarities(stringsToDNA4(targetStrings),mismatches),mismatches);

	//auto matches = match(sequence,seqLength,targetStrings,targetLength,mismatches);
	delete[] sequence;
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
		// if(matches1[i].size()!=matches[i].size())
		// {
		// 	std::cout << i << "\t" << matches1[i].size() << "\t" << matches[i].size() << std::endl;
		// }
		totalMatches += matches[i].size();
	}
	std::cout << totalMatches << std::endl;
	
	return 0;
}