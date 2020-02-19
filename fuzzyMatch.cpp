#include <cstdlib>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <numeric>


//represents up to a length 32 strand of DNA
struct DNA4{
	uint64_t ACandGT[2];
	const uint64_t& operator[](size_t i) const {return ACandGT[i];}
	uint64_t& operator[](size_t i){return ACandGT[i];}
	bool operator==(const DNA4 & o) const
	{
		return o[0]==ACandGT[0]&o[1]==ACandGT[1];
	}
	bool operator!=(const DNA4 &o)const{return !operator==(o);}
	uint32_t As(uint32_t length) const {return __builtin_popcountll(ACandGT[1]>>32&&((1ULL<<length)-1));}
	uint32_t Cs(uint32_t length) const {return __builtin_popcountll(ACandGT[0]<<64-length);}
	uint32_t Gs(uint32_t length) const {return __builtin_popcountll(ACandGT[1]>>32&&((1ULL<<length)-1));}
	uint32_t Ts(uint32_t length) const {return __builtin_popcountll(ACandGT[1]<<64-length);}

	//removes characters that are in a marked position in mask
	//eg if this = ACCGG and mask = __GG_ 
	//then this becomes AC__G
	void subtractMask(const DNA4 &maskDNA4)
	{
		uint64_t mask = maskDNA4[0] | maskDNA4[1];
		mask = ~(mask << 32 | mask >> 32 | mask);
		ACandGT[0] &= mask;
		ACandGT[1] &= mask;		
	}
};

template<uint32_t pow>
uint32_t pow2()
{
	return 2ul*pow2<pow-1>();
}
template<>
uint32_t pow2<0ul>()
{
	return 1ul;
}

constexpr uint64_t acs[128] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0x100000000ULL,0,1ULL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0x100000000ULL,0,1ULL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
								
constexpr uint64_t gts[128] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0x100000000ULL,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0x100000000ULL,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0 };

void addCharacter(DNA4 &s, char c)
{
 	s[0] <<= 1;
	s[1] <<= 1;
	s[0] &= ~(1ULL|(1ULL<<32));
	s[1] &= ~(1ULL|(1ULL<<32));
	s[0] |= acs[c];
	s[1] |= gts[c];
}

uint32_t getSimilarity(const DNA4 &a,const DNA4 &b)
{
	return __builtin_popcountll ((a[0]&b[0])|(a[1]&b[1]));
}

DNA4 createDNA4(const char *s, uint_fast8_t length)
{
	uint64_t ac = 0;
	uint64_t gt = 0;
	char c;
	// DNA4 ret{0,0};
	// for(int i = 0; i < length; ++i)
	// 	addCharacter(ret,s[i]);
	// return ret;
	switch(length)
	{
		case 32: {c = s[length-32]; ac |= acs[c] << 31; gt |= gts[c] << 31;} 
		case 31: {c = s[length-31]; ac |= acs[c] << 30; gt |= gts[c] << 30;}
		case 30: {c = s[length-30]; ac |= acs[c] << 29; gt |= gts[c] << 29;}
		case 29: {c = s[length-29]; ac |= acs[c] << 28; gt |= gts[c] << 28;} 
		case 28: {c = s[length-28]; ac |= acs[c] << 27; gt |= gts[c] << 27;}
		case 27: {c = s[length-27]; ac |= acs[c] << 26; gt |= gts[c] << 26;}
		case 26: {c = s[length-26]; ac |= acs[c] << 25; gt |= gts[c] << 25;}
		case 25: {c = s[length-25]; ac |= acs[c] << 24; gt |= gts[c] << 24;}
		case 24: {c = s[length-24]; ac |= acs[c] << 23; gt |= gts[c] << 23;}
		case 23: {c = s[length-23]; ac |= acs[c] << 22; gt |= gts[c] << 22;}
		case 22: {c = s[length-22]; ac |= acs[c] << 21; gt |= gts[c] << 21;}
		case 21: {c = s[length-21]; ac |= acs[c] << 20; gt |= gts[c] << 20;}
		case 20: {c = s[length-20]; ac |= acs[c] << 19; gt |= gts[c] << 19;}
		case 19: {c = s[length-19]; ac |= acs[c] << 18; gt |= gts[c] << 18;}
		case 18: {c = s[length-18]; ac |= acs[c] << 17; gt |= gts[c] << 17;}
		case 17: {c = s[length-17]; ac |= acs[c] << 16; gt |= gts[c] << 16;}
		case 16: {c = s[length-16]; ac |= acs[c] << 15; gt |= gts[c] << 15;}
		case 15: {c = s[length-15]; ac |= acs[c] << 14; gt |= gts[c] << 14;}
		case 14: {c = s[length-14]; ac |= acs[c] << 13; gt |= gts[c] << 13;}
		case 13: {c = s[length-13]; ac |= acs[c] << 12; gt |= gts[c] << 12;}
		case 12: {c = s[length-12]; ac |= acs[c] << 11; gt |= gts[c] << 11;}
		case 11: {c = s[length-11]; ac |= acs[c] << 10; gt |= gts[c] << 10;}
		case 10: {c = s[length-10]; ac |= acs[c] << 9; gt |= gts[c] << 9;}
		case 9: {c = s[length-9]; ac |= acs[c] << 8; gt |= gts[c] << 8;}
		case 8: {c = s[length-8]; ac |= acs[c] << 7; gt |= gts[c] << 7;}
		case 7: {c = s[length-7]; ac |= acs[c] << 6; gt |= gts[c] << 6;}
		case 6: {c = s[length-6]; ac |= acs[c] << 5; gt |= gts[c] << 5;}
		case 5: {c = s[length-5]; ac |= acs[c] << 4; gt |= gts[c] << 4;}
		case 4: {c = s[length-4]; ac |= acs[c] << 3; gt |= gts[c] << 3;}
		case 3: {c = s[length-3]; ac |= acs[c] << 2; gt |= gts[c] << 2;}
		case 2: {c = s[length-2]; ac |= acs[c] << 1; gt |= gts[c] << 1;}
		case 1: {c = s[length-1]; ac |= acs[c]; gt |= gts[c];}
	}
	return DNA4{ac,gt};
}
std::vector<DNA4> stringsToDNA4(const std::vector<std::string> &targetStrings)
{
	std::vector<DNA4> targets;
	targets.reserve(targetStrings.size());
	for(std::string s: targetStrings)
{
		targets.push_back(createDNA4(s.c_str(),s.size()));
	}
	return targets;
}

std::string DNA4ToString(const DNA4 &a, uint_fast8_t length)
{
	std::string s(length,'_');
	for(uint_fast8_t i=0;i<length;i++)
	{
		if((a[0]>>(32+i))&1)
		{
			s[length-i-1] = 'A';
		}
		else if((a[0]>>i)&1)
		{
			s[length-i-1] = 'C';
		}
		else if((a[1]>>(32+i))&1)
		{
			s[length-i-1] = 'G';
		}
		else if((a[1]>>i)&1)
		{
			s[length-i-1] = 'T';
		}
	}
	return s;
}

struct Location
{
	size_t seqID;
	size_t positionInSeq;
};

//storage of target buckets
std::vector<std::pair<DNA4,uint32_t>> buckets[32*32*32];

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

	//store targets by acg content	
	// std::vector<std::pair<DNA4,uint32_t>> ACGcontent[32*32*32];
	// for(int i = 0; i < targets.size();++i)
	// {
	// 	const DNA4 &target = targets[i];
	// 	ACGcontent[target.As(targetLength)*32*32+target.Cs(targetLength)*32+target.Gs(targetLength)].push_back(std::pair<DNA4,uint32_t>(target,i));
	// }

	int numRecordedCharacters = getSimilarity(targets[0],targets[0]);
	for(int i = 0; i < targets.size();++i)
	{
		const DNA4 &target = targets[i];
		int numAs = target.As(targetLength);
		int numCs = target.Cs(targetLength);
		int numGs = target.Gs(targetLength);
		
		for(int a = std::max(numAs-mismatches,0); a <= std::min(numAs+mismatches,numRecordedCharacters); a++)
		{
			int aSurplus = std::max(a - numAs,0);
			int aDeficit = std::max(numAs - a,0);
			for(int c = std::max(numCs-mismatches+aDeficit,0); c <= std::min(numCs+mismatches-aSurplus,numRecordedCharacters); c++)
			{
				int acSurplus = aSurplus + std::max(c - numCs,0);
				int acDeficit = aDeficit + std::max(numCs - c,0);
				for(int g = std::max(numGs-mismatches+acDeficit,0); g <= std::min(numGs+mismatches-acSurplus,numRecordedCharacters); g++)
				{
					//std::cout << numAs+aOffset << "\t" << numCs+cOffset << "\t" << numGs+gOffset << std::endl;
					//brackets[numAs+aOffset][numCs+cOffset][numGs+gOffset].push_back(std::pair<DNA4,uint32_t>(i,target));
					buckets[32*32*a+32*c+g].push_back(std::pair<DNA4,uint32_t>(target,i));
				}
			}
		}
	}
	//double totalWork;
	// for(int As = 0; As <= targetLength;As++)
	// {
	// 	for(int Cs = 0; Cs <= targetLength - As;Cs++)
	// 	{
	// 		for(int Gs = 0; Gs <= targetLength - As - Cs; Gs++)
	// 		{
		
	// 			for(int a = std::max(As-mismatches,0); a <= std::min(As+mismatches,targetLength); a++)
	// 			{
	// 				int aSurplus = std::max(a - As,0);
	// 				int aDeficit = std::max(As - a,0);
	// 				for(int c = std::max(Cs-mismatches+aDeficit,0); c <= std::min(Cs+mismatches-aSurplus,targetLength); c++)
	// 				{
	// 					int acSurplus = aSurplus + std::max(c - Cs,0);
	// 					int acDeficit = aDeficit + std::max(Cs - c,0);
	// 					uint32_t total = 0;
	// 					auto &bucket = buckets[32*32*As+32*Cs+Gs];
	// 					for(int g = std::max(Gs-mismatches+acDeficit,0); g <= std::min(Gs+mismatches-acSurplus,targetLength); g++)
	// 					{
	// 						total += ACGcontent[a*32*32+c*32+g].size();
	// 					}
	// 					bucket.reserve(total);
	// 					for(int g = std::max(Gs-mismatches+acDeficit,0); g <= std::min(Gs+mismatches-acSurplus,targetLength); g++)
	// 					{
	// 						bucket.insert(bucket.end(),ACGcontent[a*32*32+c*32+g].begin(),ACGcontent[a*32*32+c*32+g].end());
	// 					}
	// 				}
	// 			}

	// 			//totalWork += buckets[32*32*As+32*Cs+Gs].size()*buckets[32*32*As+32*Cs+Gs].size();
	// 			// if(buckets[32*32*As+32*Cs+Gs].size() > 100)
	// 			// {
	// 			// 	std::cout << As << "\t" << Cs << "\t" << Gs << "\t" << targetLength - As - Cs - Gs << "\t"; 
	// 			// 	std::cout<<buckets[32*32*As+32*Cs+Gs].size() << std::endl;

	// 			// }
	// 			// if(ACGcontent[32*32*As+32*Cs+Gs].size()>0)
	// 			// {
	// 			// 	std::cout << As << "\t" << Cs << "\t" << Gs << "\t" << targetLength - As - Cs - Gs << "\t" << ACGcontent[32*32*As+32*Cs+Gs].size() << std::endl; 				
	// 			// }
	// 			//Ts is no greater than targetLength - As - Cs - Gs
	// 		}
	// 	}
	// }
	//std::cout << totalWork << std::endl;
}

inline std::vector<std::pair<DNA4,uint32_t>> &getTargets(uint32_t a, uint32_t c, uint32_t g)
{
	return buckets[a*32*32+c*32+g];
}
uint_fast8_t sensibleMinSimilarity(uint64_t numTargets, uint_fast8_t targetSize, uint_fast8_t mismatches)
{

}
//returns a vector which contains an array and a vector for each target. The array stores the start of each mismatch bucket in the vector.
//so to access targets with mismatch 10-16 to target i you would access v[i].second[v[i].first[10]] up to v[i].second[v[i].first[17]]
//or std::pair<DNA4,uint32_t> * it = v[i].second.data(); for(int i = v[i].first[10]; i < v[i].first[17]; i++){//do something with it[i]} 
std::vector<std::pair<std::array<uint32_t,33>, std::vector<std::pair<DNA4,uint32_t>>>> getTargetSimilarities(const std::vector<DNA4> &targets,uint_fast8_t minSimilarityToStore)
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
			//we only go into this array if similarity is at least the minimum similarity
			if(similarity>=minSimilarityToStore)
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

std::vector<std::vector<Location>> simpleMatch(const std::vector<std::string> &sequences, const std::vector<std::string> &targetStrings, uint_fast8_t targetLengths, uint_fast8_t mismatches, std::string requiredMatch = "")
{
	DNA4 filterSeq = createDNA4(requiredMatch.data(),requiredMatch.size());
	bool filterExists = requiredMatch.size() > 0;
	uint_fast8_t numFilterChars = getSimilarity(filterSeq,filterSeq);

	auto targets = stringsToDNA4(targetStrings);
	auto matches = std::vector<std::vector<Location>>(targets.size());
	if(targets.size()==0)
	{
		return matches;
	}
	uint_fast8_t minimumMatches = targetLengths - mismatches;
	//auto similarities = std::array<uint64_t,24>();
	//similarities.fill(0);
	int *matched1 = new int[targetStrings.size()/2+1];
	int *matched2 = new int[targetStrings.size()/2+1];
	for(size_t seqID = 0; seqID < sequences.size();++seqID)
	{
		size_t sequenceLength = sequences[seqID].size();
		const char * sequenceString = sequences[seqID].data();
		auto sequence = createDNA4(sequenceString,targetLengths-1);
		size_t currentPos = targetLengths-1;
		while(currentPos < sequenceLength)
		{
			int numberOfMatches1 = 0;
			int numberOfMatches2 = 0;
			//add next character
			addCharacter(sequence,sequenceString[currentPos]);
			
			if(filterExists)
			{
				if(getSimilarity(sequence,filterSeq)<numFilterChars)
				{
					currentPos++;
					continue;
				}
			}
			//compare all the targets with this sequence
			for(int i = 0; i < targets.size()-1;i+=2)
			{
				//similarities[getSimilarity(sequence,targets[i])]++;
				//similarities[getSimilarity(sequence,targets[i+1])]++;
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
				//similarities[getSimilarity(sequence,targets[targets.size()-1])]++;
				if(getSimilarity(sequence,targets[targets.size()-1])>=minimumMatches)
				{
					matched1[numberOfMatches1] = targets.size()-1;
					numberOfMatches1++;
				}
			}
			for(int i = 0; i < numberOfMatches1;++i)
			{
				matches[matched1[i]].push_back(Location{seqID,currentPos});
			}
			
			for(int i = 0; i < numberOfMatches2;++i)
			{
				matches[matched2[i]].push_back(Location{seqID,currentPos});
			}
			currentPos++;
		}
	}
	
	// uint64_t sum = std::accumulate(similarities.begin(),similarities.end(),0ULL);
	// uint64_t atLeast = sum;
	// for(int i = 0; i < similarities.size();i++)
	// {
	// 	std::cout << i << '\t' << similarities[i]/(double)sum << '\t' << atLeast/(double)sum << std::endl;
	// 	atLeast -= similarities[i];
	// }
	delete[] matched1;
	delete[] matched2;
	return matches;
}

std::vector<std::vector<Location> > match(const std::vector<std::string> &sequences, const std::vector<std::string> &targetStrings, uint_fast8_t targetLength, uint_fast8_t mismatches, std::string requiredMatch = "")
{
	DNA4 filterSeq = createDNA4(requiredMatch.data(),requiredMatch.size());
	uint_fast8_t numFilterChars = getSimilarity(filterSeq,filterSeq);

	//hacky solution to approximate when it is best to use the simple solution
	//instead of the target-target similarity solution
	//this is only the case for target length = 23 and large sequence strings (> 1 million characters)
	//if the same set of targets is to be used on many sets of sequences consider 
	//precomputing target buckets and targetTargetSimilarities
	size_t totalSequenceLength = 0;
	for(const std::string &s: sequences) 
		totalSequenceLength += s.size();
	int overheadCost = mismatches;
	overheadCost += numFilterChars/2;
	if(totalSequenceLength<2000000) overheadCost++;
	if(targetStrings.size() < 2000) overheadCost++;
	if(targetStrings.size() < 200) overheadCost++;
	if(targetStrings.size() < 100) overheadCost++;
	if(targetStrings.size() < 50) overheadCost+=2;
	if(overheadCost>6)
		return simpleMatch(sequences,targetStrings,targetLength,mismatches,requiredMatch);

	bool filterExists = numFilterChars > 0;
	auto DNA4TargetStrings = stringsToDNA4(targetStrings);
	if(filterExists)
	{
		for(DNA4 &t: DNA4TargetStrings)
		{
			t.subtractMask(filterSeq);
		}
	}

	uint_fast8_t importantChars = getSimilarity(DNA4TargetStrings[0],DNA4TargetStrings[0]);

	//specifically chosen for length 23
	//TODO calculate optimum for strings of length n
	uint_fast8_t minimumSimilarityToCheckBuckets = 12;
	if(mismatches>4) //otherwise memory exceeds 1.5GB
		minimumSimilarityToCheckBuckets = 13;
	if(targetStrings.size()<2000)//with fewer targets we are less likely to come accross a high match by chance so we need to accept lower similarity
		minimumSimilarityToCheckBuckets = 11;

	fillTargetBuckets(DNA4TargetStrings,targetLength,mismatches);
	const std::vector<std::pair<std::array<uint32_t,33>, std::vector<std::pair<DNA4,uint32_t>>>> targetTargetSimilarities = getTargetSimilarities(DNA4TargetStrings,minimumSimilarityToCheckBuckets - mismatches);

	std::cout << "starting to match" << std::endl;
	auto matches = std::vector<std::vector<Location> >(targetTargetSimilarities.size());
	uint_fast8_t minimumMatches = importantChars - mismatches;

	int *matched1 = new int[targetTargetSimilarities.size()/2+1];
	int *matched2 = new int[targetTargetSimilarities.size()/2+1];
	bool specificTargetMatched = new bool[targetTargetSimilarities.size()];
	for(size_t seqID = 0; seqID < sequences.size();++seqID)
	{
		size_t sequenceLength = sequences[seqID].size();
		const char * sequenceString = sequences[seqID].data();
		DNA4 rawSequence = createDNA4(sequenceString,targetLength-1);
		size_t currentPos = targetLength-1;
		while(currentPos < sequenceLength)
		{
			//std::cout << currentPos << std::endl;
			int numberOfMatches1 = 0;
			int numberOfMatches2 = 0;
			//add next character
			addCharacter(rawSequence,sequenceString[currentPos]);
			DNA4 sequence = rawSequence;
			if(filterExists)
			{
				if(getSimilarity(sequence,filterSeq)<numFilterChars)
				{
					currentPos++;
					continue;
				}
				
				sequence.subtractMask(filterSeq);
			}

			int numAs = sequence.As(targetLength);
			int numCs = sequence.Cs(targetLength);
			int numGs = sequence.Gs(targetLength);

			//std::cout << numAs << "\t" << numCs << "\t" << numGs << std::endl;
			const std::pair<DNA4,uint32_t> *target = buckets[numAs*32*32+numCs*32+numGs].data();
			const std::pair<DNA4,uint32_t> *end = target + buckets[numAs*32*32+numCs*32+numGs].size();	
			uint_fast8_t maxSimilarity = minimumSimilarityToCheckBuckets - 1;
			if(target==end)
			{
				currentPos++;
				continue;
			}
			//uint totalChecked = 0;
			//compare all the targets with this sequence
			while(target < end-1)
			{
				uint_fast8_t similarity1 = getSimilarity(sequence,target->first);
				uint_fast8_t similarity2 = getSimilarity(sequence,(target+1)->first);
				//totalChecked+=2;

				if(similarity1>maxSimilarity||similarity2>maxSimilarity)
				{
					if(similarity2>similarity1)
					{
						target = target+1;
						similarity1 = similarity2;
					}
					maxSimilarity = similarity1;
					auto &t = targetTargetSimilarities[target->second];

					uint_fast8_t lowestSimilarity = std::max(similarity1-mismatches,0);
					uint_fast8_t highestSimilarity = std::min(similarity1+mismatches,(int)importantChars);
					target = t.second.data() + t.first[importantChars-highestSimilarity];
					end = t.second.data() + t.first[importantChars-lowestSimilarity+1];

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
			//std::cout << "\t\t" << totalChecked << std::endl;
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
				matches[matched1[i]].push_back(Location{seqID,currentPos});
			}
			for(int i = 0; i < numberOfMatches2;++i)
			{
				matches[matched2[i]].push_back(Location{seqID,currentPos});
			}
			currentPos++;
		}
	}
	delete[] matched1;
	delete[] matched2;
	return matches;
}

char randNucleotide()
{
	constexpr char DNA[4] = {'A','C','G','T'};
	return DNA[rand()%4];
}

int main()
{
	constexpr size_t seqLength = 10000000;
	constexpr size_t numTargets = 10000;
	constexpr size_t targetLength = 23;
	constexpr size_t mismatches = 4;
	const std::string filter = "";
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
		//targets will always have the filter
		for(int j = 0; j < filter.size();j++)
			target[targetLength-1-j]=filter[j];
		targetStrings.push_back(target);
	}
	std::vector<std::string> sequences;
	sequences.push_back(std::string(sequence));
	auto matches = match(sequences,targetStrings, targetLength,mismatches,filter);

	delete[] sequence;
	size_t totalMatches = 0;
	for(int i = 0; i < matches.size(); ++i)
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