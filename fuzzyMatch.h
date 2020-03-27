#include <cstdlib>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <memory>
#include <iterator>
#include <algorithm>
#include <bitset>
#include "stopwatch.h"

//represents up to a length 32 strand of DNA
struct DNA4{
	uint64_t ACandGT[2];
	const uint64_t& operator[](size_t i) const {return ACandGT[i];}
	uint64_t& operator[](size_t i){return ACandGT[i];}
	bool operator==(const DNA4 & o) const
	{
		return (o[0]==ACandGT[0])&&(o[1]==ACandGT[1]);
	}
	//this is not a lexigraphical sort
	bool operator>(const DNA4 &o)const{
		if(ACandGT[0]!=o[0])
			return ACandGT[0]>o[0];
		return ACandGT[1]>o[1];
	}
	bool operator<(const DNA4 &o)const{
		if(ACandGT[0]!=o[0])
			return ACandGT[0]<o[0];
		return ACandGT[1]<o[1];
	}
	bool operator>=(const DNA4 &o)const{return !operator<(o);}
	bool operator<=(const DNA4 &o)const{return !operator>(o);}
	bool operator!=(const DNA4 &o)const{return !operator==(o);}
	uint32_t As(uint32_t length) const {return __builtin_popcountll(ACandGT[1]>>32&&((1ULL<<length)-1));}
	uint32_t Cs(uint32_t length) const {return __builtin_popcountll(ACandGT[0]<<(64-length));}
	uint32_t Gs(uint32_t length) const {return __builtin_popcountll(ACandGT[1]>>32&&((1ULL<<length)-1));}
	uint32_t Ts(uint32_t length) const {return __builtin_popcountll(ACandGT[1]<<(64-length));}

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

	void addWildcards(const char *s, uint_fast8_t length, char wildcardChar)
	{
		uint64_t wildcardPositions=0;
		for(uint_fast8_t i = 0; i < length; i++)
		{
			wildcardPositions <<= 1;
			if(s[i]==wildcardChar)
			{
				wildcardPositions |= 1ULL;
			}
		}
		//duplicate over the two 32 bit sections of the 64bit
		wildcardPositions = wildcardPositions | (wildcardPositions << 32);
		ACandGT[0] |= wildcardPositions;
		ACandGT[1] |= wildcardPositions;
	}
	
	uint_fast8_t getNumMatchableChars() const
	{
		return __builtin_popcount(uint32_t((ACandGT[0]>>32)|(ACandGT[1]>>32)|ACandGT[0]|ACandGT[1]));
	}
	bool matchableCharacterAtPosition(uint32_t i,uint_fast8_t length) const
	{
		uint64_t mask = (1ULL|(1ULL<<32)<<(length-i-1));
		return (ACandGT[0]|ACandGT[1])&mask;
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
uint32_t nChooseK(uint32_t n,uint32_t k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( uint32_t i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

template <typename T>
//returns a vector of all combinations (n choose k) of elements stored at begin
std::vector<std::vector<T>> getAllCombinations(T* begin, uint32_t n, uint32_t k)
{
	std::vector<std::vector<T> > allCombinations;
	std::string bitmask(k, 1); // K leading 1's
    bitmask.resize(n, 0); // N-K trailing 0's
 
    // print integers and permute bitmask
    do {
		std::vector<T> combination;
        for (uint32_t i = 0; i < n; ++i) // [0..N-1] integers
        {
            if (bitmask[i])
			{
				combination.push_back(begin[i]);
			}
        }
		std::sort(combination.begin(),combination.end());
        allCombinations.push_back(combination);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
	return allCombinations;
}

constexpr uint64_t acs[256] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0x100000000ULL,0,1ULL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0x100000000ULL,0,1ULL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
								
constexpr uint64_t gts[256] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
								0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0x100000000ULL,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0x100000000ULL,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

void addCharacter(DNA4 &s, unsigned char c)
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
	unsigned char c;
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

struct TargetBucket
{
	DNA4 *begin_DNA;
	uint32_t *begin_position;
	uint64_t size;
};

class TargetContainer
{
	public:
	TargetContainer(const std::vector<DNA4> &targets, uint32_t targetLength, uint32_t mismatches, const DNA4 &filter = DNA4{0,0}, uint64_t maxIndexSize = ~0ULL)
	:stringLength(targetLength),buckets_DNA4(nullptr),buckets_positions(nullptr),bucketsForLastHash(nullptr),numberOfVariableChars(stringLength - filter.getNumMatchableChars()),mismatches(mismatches),numberOfTargets(targets.size())
	{
		if(filter == DNA4{0,0}) 
			hasFilter = false;
		numberOfDivisions = optimumNumberOfDivisions(maxIndexSize);
		if(numberOfDivisions==0)
			return;
		totalBuckets = numberOfDivisions*arrangementsPerDivision*bucketsPerArrangement;
		totalHashmaps = numberOfDivisions*arrangementsPerDivision;
		buckets_positions = new std::vector<uint32_t>[totalBuckets];
		//buckets = new std::vector<std::pair<DNA4,uint32_t> >[totalBuckets];
		bucketsForLastHash = new uint32_t[totalHashmaps];
		createHashHelpers(numberOfDivisions,filter);
		putTargetsInHashmaps(targets);
		buckets_DNA4 = new std::vector<DNA4>[totalBuckets];
		for(uint i = 0; i < totalBuckets;++i)
		{
			std::vector<uint32_t> &bucket = buckets_positions[i];
			buckets_DNA4[i] = std::vector<DNA4>(bucket.size());
			for(uint j = 0; j < bucket.size(); j++)
			{
				buckets_DNA4[i][j] = targets[buckets_positions[i][j]];
			}
		}

	}

	~TargetContainer()
	{
		delete[] buckets_DNA4;
		delete[] buckets_positions;
		delete[] bucketsForLastHash;
	}
	//returns the optimum number of divisions to use, if it is not worth using divisions zero is returned
	uint32_t optimumNumberOfDivisions(uint64_t maxIndexSize)
	{
		//std::cout << "Divisions\tdivisionSize\tmismatchesPerDivision\tarrangementsPerDivision\ttargetsPerBucket\ttotalWork\ttotalSize"<<std::endl;
		//calculate the optimum division size
		uint64_t minimumTotalWork = numberOfTargets;
		uint32_t bestNumberOfDivisions = 0;
		for(uint32_t divisions = 1; divisions <= numberOfVariableChars; ++divisions)
		{
			divisionSize = numberOfVariableChars/divisions;//round down
			mismatchesPerDivision = mismatches/divisions;//round down
			arrangementsPerDivision = nChooseK(divisionSize,mismatchesPerDivision);
			bucketsPerArrangement = std::pow(4,divisionSize-mismatchesPerDivision);//where 4 is the size of the alphabet
			uint64_t totalSize = arrangementsPerDivision*divisions*(bucketsPerArrangement*sizeof(std::vector<std::pair<DNA4,uint32_t> >) + numberOfTargets*(sizeof(std::pair<DNA4,uint32_t>)));
			//overhead to generate a hash for each arrangement within each division is approximately = numberOfVariableChars
			//each arrangement contains on average targets/bucketsPerArrangement
			double totalWork = (numberOfVariableChars*10+numberOfTargets/bucketsPerArrangement)*(arrangementsPerDivision*divisions);
			if(totalWork < minimumTotalWork&&totalSize<maxIndexSize)
			{
				minimumTotalWork = totalWork;
				bestNumberOfDivisions = divisions;
			}
			//std::cout << divisions <<'\t'<< divisionSize <<'\t'<< mismatchesPerDivision <<'\t'<< arrangementsPerDivision <<'\t'<< numberOfTargets/bucketsPerArrangement <<'\t'<< totalWork << '\t'<< totalSize << std::endl;
		}
		if(bestNumberOfDivisions==0)
			return 0;
		divisionSize = numberOfVariableChars/bestNumberOfDivisions;//round down
		mismatchesPerDivision = mismatches/bestNumberOfDivisions;//round down
		arrangementsPerDivision = nChooseK(divisionSize,mismatchesPerDivision);
		bucketsPerArrangement = std::pow(4,divisionSize-mismatchesPerDivision);//where 4 is the size of the alphabet
		//std::cout << bestNumberOfDivisions << std::endl;
		
		return bestNumberOfDivisions;
	}

	void createHashHelpers(uint32_t numberOfDivisions,const DNA4 &filter)
	{
		std::unique_ptr<uint32_t[]> positions(new uint32_t[numberOfVariableChars]);
		uint32_t ithMatchableCharacter = 0;
		for(uint32_t i = 0; i < numberOfVariableChars;++i)
		{
			if(!filter.matchableCharacterAtPosition(ithMatchableCharacter,stringLength))
			{
				positions[i] = stringLength - ithMatchableCharacter - 1;
			}
			else
			{
				i--;
			}	
			ithMatchableCharacter++;
		}
		//positions now contains the position of each matchable character counting from the end
		//if there was no filter or all the filter characters were at the end of the string then
		//the array will just contain the numbers numberOfVariableChars-1 to 0 in decreasing order

		for(uint32_t i = 0; i < numberOfDivisions;++i)
		{
			auto combinationsInDivision = getAllCombinations<uint32_t>(positions.get()+i*divisionSize,divisionSize,divisionSize-mismatchesPerDivision);
			hashHelpers.insert(hashHelpers.end(),combinationsInDivision.begin(),combinationsInDivision.end());
		}	
	}

	void getBuckets(const DNA4 &string, TargetBucket* bucketsForString)
	{
		slowestHashInTheWorld(string);
		for(uint32_t j = 0; j < totalHashmaps;j++)
		{
			uint32_t bucket = bucketsForLastHash[j];
			bucketsForString[j] = TargetBucket{buckets_DNA4[bucket].data(),buckets_positions[bucket].data(),buckets_positions[bucket].size()};
		}
	}

	void slowestHashInTheWorld(const DNA4 &string) const
	{
		//bitmask of positions of each character
	#pragma GCC diagnostic ignored "-Wnarrowing"
		uint32_t chars[4] = {string[0]>>32,string[1]>>32,(uint32_t)string[0],uint32_t(string[1])};//the high 32 bits of DNA4 and the low 32 bits
	#pragma GCC diagnostic pop
		for(uint32_t i = 0; i < totalHashmaps; ++i)
		{
			uint32_t hash = 0;
			const std::vector<uint32_t> &relevantPositions = hashHelpers[i];
			for(uint32_t j = 0; j < relevantPositions.size();j++)
			{
				uint32_t pos = relevantPositions[j];
				hash <<= 2;
				//only uses 3 of the four character of the alphabet, the fourth is represented by 0's
				//this will hash strings that contain none alphabet characters as if the had the 4th character there instead
				hash |= ((chars[0]>>pos)&1UL) | ((chars[1]>>pos)&1UL)*2 | ((chars[2]>>pos)&1UL)*3;
			}
			bucketsForLastHash[i] = bucketsPerArrangement*i + hash;
		}
	}

	void putTargetsInHashmaps(const std::vector<DNA4> &targets)
	{
		for(uint32_t i = 0; i < targets.size();++i)
		{
			DNA4 target = targets[i];
			slowestHashInTheWorld(target);
			for(uint32_t j = 0; j < totalHashmaps;j++)
			{
				buckets_positions[bucketsForLastHash[j]].push_back(i);
			}
		}
	}

	uint32_t numberOfHashmaps() const {return totalHashmaps;}

	operator bool(){return numberOfDivisions!=0;}

	private:
	uint32_t stringLength;
	//std::vector<std::pair<DNA4,uint32_t> > *buckets;
	std::vector<DNA4> *buckets_DNA4;
	std::vector<uint32_t> *buckets_positions;
	mutable uint32_t * bucketsForLastHash;
	std::vector<std::vector<uint32_t> > hashHelpers;
	bool hasFilter;
	DNA4 filter;
	uint32_t numberOfVariableChars;
	uint32_t mismatches;
	uint32_t numberOfTargets;
	uint32_t numberOfDivisions;
	uint32_t divisionSize;
	uint32_t mismatchesPerDivision;
	uint32_t arrangementsPerDivision;
	uint64_t bucketsPerArrangement;
	uint32_t totalHashmaps;
	uint64_t totalBuckets;


};

struct offtarget
{
	uint32_t seqID;
	size_t positionInSeq;
	DNA4 offTargetSequence;
	uint32_t mismatches;
	bool operator!=(const offtarget &o) { return o.seqID!=seqID||o.positionInSeq!=positionInSeq;}
	bool operator==(const offtarget &o) { return !operator!=(o);}
};

std::vector<std::vector<offtarget>> simpleMatch(const std::vector<std::string> &sequences, const std::vector<std::string> &targetStrings, uint_fast8_t targetLengths, uint_fast8_t mismatches, std::string requiredMatch = "")
{
	DNA4 filterSeq = createDNA4(requiredMatch.data(),requiredMatch.size());
	bool filterExists = requiredMatch.size() > 0;
	filterSeq.addWildcards(requiredMatch.data(),requiredMatch.size(),'*');
	uint_fast8_t numFilterChars = filterSeq.getNumMatchableChars();

	auto targets = stringsToDNA4(targetStrings);
	auto matches = std::vector<std::vector<offtarget>>(targets.size());
	if(targets.size()==0)
	{
		return matches;
	}
	
	if(filterExists)
	{
		for(DNA4 &t: targets)
		{
			t.subtractMask(filterSeq);
		}
	}
	stop_watch s;
	s.start();
	std::cout << "Using simple method "<< std::endl;
	uint_fast8_t minimumMatches = targets.at(0).getNumMatchableChars() - mismatches;
	if(targets.at(0).getNumMatchableChars() < mismatches) minimumMatches = 0;//avoid underflow
	//auto similarities = std::array<uint64_t,24>();
	//similarities.fill(0);
	int *matched1 = new int[targetStrings.size()/2+1];
	int *matched2 = new int[targetStrings.size()/2+1];
	uint64_t numComparisons = 0;
	for(uint32_t seqID = 0; seqID < sequences.size();++seqID)
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
			numComparisons+= targets.size();
			//compare all the targets with this sequence
			for(size_t i = 0; i < targets.size()-1;i+=2)
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
				matches[matched1[i]].push_back(offtarget{seqID,currentPos,sequence,targets[matched1[i]].getNumMatchableChars()-getSimilarity(sequence,targets[matched1[i]])});
			}
			
			for(int i = 0; i < numberOfMatches2;++i)
			{
				matches[matched2[i]].push_back(offtarget{seqID,currentPos,sequence,targets[matched2[i]].getNumMatchableChars()-getSimilarity(sequence,targets[matched2[i]])});
			}
			currentPos++;
		}
	}

	delete[] matched1;
	delete[] matched2;
	s.stop();
	std::cout<< numComparisons << " comparisons took " << s << std::endl;
	return matches;
}

std::vector<std::vector<offtarget> > match(const std::vector<std::string> &sequences, const std::vector<std::string> &targetStrings, uint_fast8_t targetLength, uint_fast8_t mismatches, std::string requiredMatch = "", uint64_t maxIndexSize = ~0ULL)
{
	auto matches = std::vector<std::vector<offtarget>>(targetStrings.size());
	if(targetStrings.size()==0||sequences.size()==0)
	{
		return matches;
	}
	uint_fast8_t targetLengths = targetStrings[0].size();

	DNA4 filterSeq = createDNA4(requiredMatch.data(),requiredMatch.size());
	filterSeq.addWildcards(requiredMatch.data(),requiredMatch.size(),'*');
	uint_fast8_t numFilterChars = filterSeq.getNumMatchableChars();

	auto DNA4TargetStrings = stringsToDNA4(targetStrings);
	bool filterExists = numFilterChars>0;
	if(filterExists)
	{
		for(DNA4 &t: DNA4TargetStrings)
		{
			t.subtractMask(filterSeq);
		}
	}
	stop_watch s;
	s.start();
	TargetContainer targetContainer(DNA4TargetStrings, targetLengths, mismatches, filterSeq, maxIndexSize);
	if(!targetContainer)
	{
		return simpleMatch(sequences,targetStrings,targetLength,mismatches,requiredMatch);
	}
	TargetBucket * bucketsArray = new TargetBucket[targetContainer.numberOfHashmaps()];
	s.stop();
	std::cout << "indexing targets took "<< s << std::endl;
	s.start();
	uint_fast8_t minimumMatches = DNA4TargetStrings.at(0).getNumMatchableChars() - mismatches;
	if(DNA4TargetStrings.at(0).getNumMatchableChars() < mismatches) minimumMatches = 0;//avoid underflow
	uint64_t comparisons = 0;
	uint64_t naiveComparisons = 0;
	for(uint32_t seqID = 0; seqID < sequences.size();++seqID)
	{
		uint32_t sequenceLength = sequences[seqID].size();
		const char * sequenceString = sequences[seqID].data();
		auto sequence = createDNA4(sequenceString,targetLengths-1);
		size_t currentPos = targetLengths-1;
		while(currentPos < sequenceLength)
		{
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
			naiveComparisons += targetStrings.size();
			targetContainer.getBuckets(sequence,bucketsArray);
			for(size_t i = 0; i < targetContainer.numberOfHashmaps(); ++i)
			{
				TargetBucket bucket = bucketsArray[i];
				comparisons += bucket.size;
				for(uint32_t targetNum = 0; targetNum<bucket.size; targetNum++)
				{
					uint32_t similarity = getSimilarity(sequence,bucket.begin_DNA[targetNum]);
					if(similarity>=minimumMatches)
					{
						auto & targetMatches = matches[bucket.begin_position[targetNum]];
						//make sure we haven't added this match from another hashmap already
						if(targetMatches.size()==0||targetMatches.back().positionInSeq!=currentPos||targetMatches.back().seqID!=seqID)
						{
							targetMatches.push_back(offtarget{seqID,currentPos,sequence,bucket.begin_DNA[targetNum].getNumMatchableChars()-similarity});
						}
					}
				}
			}
			currentPos++;
		}
	}
	s.stop();
	std::cout<< "comparisons took " << s << std::endl;
	std::cout << "performed " << comparisons << " comparisons; naive would perform " << naiveComparisons <<" comparisons" << std::endl;
	delete[] bucketsArray;
	return matches;
}

