#pragma once

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
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//represents up to a length 32 strand of DNA
struct DNA4{

	static void setLength(uint32_t length)
	{
		DNA4::length = length;
		DNA4::lowLengthMask = (1ULL<<length)-1;
		DNA4::highLengthMask = lowLengthMask << 32;
		DNA4::fullMask = DNA4::lowLengthMask | DNA4::highLengthMask;
	}

	static uint32_t getLength(){return length;}
	uint64_t ACandGT[2];

	DNA4(const char *s)
	{
		uint64_t ac = 0;
		uint64_t gt = 0;
		unsigned char c;

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
		ACandGT[0] = ac;
		ACandGT[1] = gt;
	}

	DNA4(const std::string &s)
	:DNA4(s.data())
	{}

	DNA4()
	:ACandGT{0,0}
	{}
	
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
	uint32_t numAs() const {return __builtin_popcountll(ACandGT[0]&highLengthMask);}
	uint32_t numCs() const {return __builtin_popcountll(ACandGT[0]&lowLengthMask);}
	uint32_t numGs() const {return __builtin_popcountll(ACandGT[1]&highLengthMask);}
	uint32_t numTs() const {return __builtin_popcountll(ACandGT[1]&lowLengthMask);}

	uint32_t As() const {return ACandGT[0]>>32&lowLengthMask;}
	uint32_t Cs() const {return ACandGT[0]&lowLengthMask;}
	uint32_t Gs() const {return ACandGT[1]>>32&lowLengthMask;}
	uint32_t Ts() const {return ACandGT[1]&lowLengthMask;}

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

	bool addCharacter(unsigned char c)
	{
		ACandGT[0] <<= 1;
		ACandGT[1] <<= 1;
		uint64_t mask = (~(1ULL|(1ULL<<32)))&DNA4::fullMask;
		ACandGT[0] &= mask;
		ACandGT[1] &= mask;
		ACandGT[0] |= acs[c];
		ACandGT[1] |= gts[c];
		return c;
	}
	void addWildcards(const char *s, char wildcardChar)
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
	bool matchableCharacterAtPosition(uint32_t i) const
	{
		uint64_t mask = (1ULL|(1ULL<<32)<<(length-i-1));
		return (ACandGT[0]|ACandGT[1])&mask;
	}
	uint32_t getSimilarity(const DNA4 &b) const
	{
		return __builtin_popcountll ((ACandGT[0]&b[0])|(ACandGT[1]&b[1]));
	}
	std::string toString() const
	{
		std::string s(DNA4::length,'N');
		for(uint_fast8_t i=0;i<DNA4::length;i++)
		{
			if((ACandGT[0]>>(32+i))&1)
			{
				s[DNA4::length-i-1] = 'A';
			}
			else if((ACandGT[0]>>i)&1)
			{
				s[DNA4::length-i-1] = 'C';
			}
			else if((ACandGT[1]>>(32+i))&1)
			{
				s[DNA4::length-i-1] = 'G';
			}
			else if((ACandGT[1]>>i)&1)
			{
				s[DNA4::length-i-1] = 'T';
			}
		}
		return s;
	}
	private:
	static uint_fast8_t length;
	static uint64_t lowLengthMask;
	static uint64_t highLengthMask;
	static uint64_t fullMask;
};

uint64_t DNA4::highLengthMask = ~0ULL;
uint64_t DNA4::lowLengthMask = ~0ULL;
uint64_t DNA4::fullMask = ~0ULL;
uint_fast8_t DNA4::length = 0;

class Filter
{
	public:
	Filter(std::string s)
	:filterSeq((std::string(DNA4::getLength()-s.size(),'x')+s).data()),numFilterChars(0)
	{
		if(s.size()==0)
			return;
		if(s.size()!=DNA4::getLength())
		{
			std::cout << "WARNING: filter string \'" << s << "\' is " << s.size() << " characters. Target length is " << DNA4::getLength() << " characters." << std::endl;
			if(s.size() < DNA4::getLength())
			{
				std::cout << "Pad string to the complete length." << std::endl;
				std::cout << "eg "<< s << std::string(DNA4::getLength()-s.size(),'x') << " if the filter is at the start of the target," << std::endl;
				std::cout << "or "<< std::string(DNA4::getLength()-s.size(),'x') << s << " if the filter is at the end of the target." << std::endl;
			}
			else
			{
				std::cout << "Remember to set the length of targets by calling DNA4::setLength(val) before calling any other functions." << std::endl;
			}
		}
		filterSeq.addWildcards(s.data(),'*');
		numFilterChars = filterSeq.getNumMatchableChars();
	}

	bool passes(const DNA4 & seq) const
	{
		return filterSeq.getSimilarity(seq)==numFilterChars;
	}

	//adds filter characters to target in the positions where filter has matchable characters
	void addFilterCharsToTarget(DNA4 & seq) const
	{
		seq[0]|=filterSeq[0];
		seq[1]|=filterSeq[1];
	}

	bool exists() const {return numFilterChars;}

	uint32_t numberOfFilterChars() const {return numFilterChars;}

	DNA4 asDNA4() const {return filterSeq;}

	private:
	DNA4 filterSeq;
	uint_fast8_t numFilterChars;
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

std::vector<DNA4> stringsToDNA4(const std::vector<std::string> &targetStrings)
{
	std::vector<DNA4> targets;
	targets.reserve(targetStrings.size());
	for(std::string s: targetStrings)
	{
		targets.push_back(DNA4(s.c_str()));
	}
	return targets;
}

class DNA4Set
{
	public:
	//the internal hash function used is based upon the filter
	//inserting targets that do not pass the filter will be less efficient but still correct
	DNA4Set(Filter filter = std::string())
	:numTargets(0),buckets(nullptr),mask(0),offset(0)
	{
		//find postion and size of longest run of variable characters
		calculateGoodFilter(filter.asDNA4());
		buckets = new std::vector<DNA4>[mask+1];
	}

	DNA4Set(const DNA4Set& o)
	:numTargets(o.numTargets),buckets(new std::vector<DNA4>[o.mask+1]),mask(o.mask),offset(o.offset)
	{
		copyBuckets(o);
	}
 
	// Move constructor
	DNA4Set(DNA4Set&& o)
	:numTargets(o.numTargets),buckets(o.buckets),mask(o.mask),offset(o.offset)
	{
		o.buckets = nullptr; 
	}
 
	// Copy assignment
	DNA4Set& operator=(const DNA4Set& o)
	{
		// Self-assignment detection
		if (&o == this)
			return *this;
 
		//if we need a different number of buckets
		if(mask!=o.mask)
		{
			delete[] buckets;
			mask = o.mask;
			buckets = new std::vector<DNA4>[o.mask+1];
		}
		numTargets = o.numTargets;
		offset = o.offset;
		copyBuckets(o);

		return *this;
	}
 
	// Move assignment
	DNA4Set& operator=(DNA4Set&& o)
	{
		// Self-assignment detection
		if (&o == this)
			return *this;
 
		delete[] buckets;

		numTargets = o.numTargets;
		buckets = o.buckets;
		o.buckets = nullptr;
		mask = o.mask;
		offset = o.offset;

		return *this;
	}

	~DNA4Set()
	{
		delete[] buckets;
	}

	bool contains(DNA4 target)
	{
		std::vector<DNA4> & bucket = getBucket(target);
		for(uint32_t i = 0; i < bucket.size();++i)
		{
			if(bucket[i]==target)
				return true;
		}
		return false;
	}

	bool insert(DNA4 target)
	{
		std::vector<DNA4> & bucket = getBucket(target);
		for(uint32_t i = 0; i < bucket.size();++i)
		{
			if(bucket[i]==target)
				return false;
		}
		numTargets++;
		bucket.push_back(target);
		return true;
	}

	bool remove(DNA4 target)
	{
		std::vector<DNA4> & bucket = getBucket(target);
		for(uint32_t i = 0; i < bucket.size();++i)
		{
			if(bucket[i]==target)
			{
				numTargets--;
				bucket.erase(bucket.begin()+i);
				return true;
			}
		}
		return false;
	}

	size_t size() const {return numTargets;}

	std::vector<DNA4> getAllTargets() const
	{
		std::vector<DNA4> targets(numTargets);
		size_t total = 0;
		for(uint32_t b = 0; b < mask+1; ++b)
		{
			auto &bucket = buckets[b];
			for(uint32_t i = 0; i < bucket.size(); ++i)
			{
				targets[total] = bucket[i];
				total++;
			}
		}
		return targets;
	}

	private:

	inline void copyBuckets(const DNA4Set& o)
	{
		for(uint32_t i = 0; i < mask+1;++i)
		{
			buckets[i] = o.buckets[i];
		}
	}

	inline std::vector<DNA4> & getBucket(const DNA4 & target)
	{
		return buckets[((target[0]|target[1])>>offset)&mask];
	}

	void calculateGoodFilter(const DNA4 &filter)
	{
		int longestRun = 0;
		int endOfLongestRun = 0;
		int start = -1;
		uint32_t H1 = filter.As();
		uint32_t H2 = filter.Gs();
		uint32_t L1 = filter.Cs();
		uint32_t L2 = filter.Ts();
		//all are the same or both sets of highs and lows are different
		uint32_t equallyVariablePositions = (~((H1^L1) | (H2^L2) | (H1^H2))) |
											((H1^H2) & (L1^L2));
		for(uint32_t i = 0;i < DNA4::getLength();++i)
		{
			//this is another variable position
			if(equallyVariablePositions&(1ULL<<i))
			{
				if(int(i-start) > longestRun)
				{
					endOfLongestRun = i;
					longestRun = i-start;
				}
			}
			else
			{
				start=i;
			}
		}
		offset = endOfLongestRun-longestRun + 1;
		longestRun = std::min(longestRun,20);
		mask = (1ULL << longestRun)-1;
	}

	size_t numTargets;
	std::vector<DNA4> * buckets;
	uint64_t mask;
	uint32_t offset;
};

class OffTargetFinder
{

	struct TargetBucket
	{
		DNA4 *begin_DNA;
		uint32_t *begin_position;
		uint64_t size;
	};

	class TargetContainer
	{
		public:
		TargetContainer(uint32_t mismatches, Filter filter)
		:buckets_DNA4(nullptr),buckets_positions(nullptr),bucketsForLastHash(nullptr),filter(filter),numberOfVariableChars(DNA4::getLength() - filter.numberOfFilterChars()),mismatches(mismatches),numTargets(0), good(false)
		{}

		~TargetContainer()
		{
			delete[] buckets_DNA4;
			delete[] buckets_positions;
			delete[] bucketsForLastHash;
		}

		bool addTargets(const std::vector<DNA4> &targets,uint64_t maxIndexSize = ~0ULL)
		{
			numTargets = targets.size();
			uint32_t numberOfDivisions = optimumNumberOfDivisions(maxIndexSize);
			if(!good)
				return false;
			uint32_t totalBuckets = numberOfDivisions*arrangementsPerDivision*bucketsPerArrangement;
			totalHashmaps = numberOfDivisions*arrangementsPerDivision;
			buckets_positions = new std::vector<uint32_t>[totalBuckets];
			bucketsForLastHash = new uint32_t[totalHashmaps];
			createHashHelpers(numberOfDivisions,mismatchesPerDivision);	
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
			return true;
		}

		//returns the optimum number of divisions to use, if it is not worth using divisions zero is returned
		uint32_t optimumNumberOfDivisions(uint64_t maxIndexSize)
		{
			if(numberOfVariableChars<mismatches)
				return 0;
			//std::cout << "Divisions\tdivisionSize\tmismatchesPerDivision\tarrangementsPerDivision\ttargetsPerBucket\ttotalWork\ttotalSize"<<std::endl;
			//calculate the optimum division size
			uint64_t minimumTotalWork = numTargets;
			uint32_t bestNumberOfDivisions = 0;
			for(uint32_t divisions = 1; divisions <= numberOfVariableChars; ++divisions)
			{
				divisionSize = numberOfVariableChars/divisions;//round down
				uint32_t mismatchesPerDivision = mismatches/divisions;//round down
				arrangementsPerDivision = nChooseK(divisionSize,mismatchesPerDivision);
				bucketsPerArrangement = std::pow(4,divisionSize-mismatchesPerDivision);//where 4 is the size of the alphabet
				uint64_t totalSize = arrangementsPerDivision*divisions*(bucketsPerArrangement*sizeof(std::vector<std::pair<DNA4,uint32_t> >) + numTargets*(sizeof(std::pair<DNA4,uint32_t>)));
				//overhead to generate a hash for each arrangement within each division is approximately = numberOfVariableChars
				//each arrangement contains on average targets/bucketsPerArrangement
				double totalWork = (numberOfVariableChars*10+numTargets/bucketsPerArrangement)*(arrangementsPerDivision*divisions);
				if(totalWork < minimumTotalWork&&totalSize<maxIndexSize)
				{
					minimumTotalWork = totalWork;
					bestNumberOfDivisions = divisions;
				}
				//std::cout << divisions <<'\t'<< divisionSize <<'\t'<< mismatchesPerDivision <<'\t'<< arrangementsPerDivision <<'\t'<< numberOfTargets/bucketsPerArrangement <<'\t'<< totalWork << '\t'<< totalSize << std::endl;
			}
			if(bestNumberOfDivisions==0)
				return 0;
			good = true;
			divisionSize = numberOfVariableChars/bestNumberOfDivisions;//round down
			mismatchesPerDivision = mismatches/bestNumberOfDivisions;//round down
			arrangementsPerDivision = nChooseK(divisionSize,mismatchesPerDivision);
			bucketsPerArrangement = std::pow(4,divisionSize-mismatchesPerDivision);//where 4 is the size of the alphabet
			//std::cout << bestNumberOfDivisions << std::endl;
			
			return bestNumberOfDivisions;
		}

		void createHashHelpers(uint32_t numberOfDivisions, uint32_t mismatchesPerDivision)
		{
			std::unique_ptr<uint32_t[]> positions(new uint32_t[numberOfVariableChars]);
			uint32_t ithMatchableCharacter = 0;
			for(uint32_t i = 0; i < numberOfVariableChars;++i)
			{
				if(!filter.asDNA4().matchableCharacterAtPosition(ithMatchableCharacter))
				{
					positions[i] = DNA4::getLength() - ithMatchableCharacter - 1;
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
					//this will hash strings that contain non alphabet characters as if the had the 4th character there instead
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

		uint32_t numberOfTargets() const {return numTargets;}

		uint32_t numberOfHashmaps() const {return totalHashmaps;}

		operator bool(){return good;}

		private:
		std::vector<DNA4> *buckets_DNA4;
		std::vector<uint32_t> *buckets_positions;
		mutable uint32_t * bucketsForLastHash;
		std::vector<std::vector<uint32_t> > hashHelpers;
		Filter filter;
		uint32_t numberOfVariableChars;
		uint32_t mismatches;
		uint32_t numTargets;
		bool good;
		uint32_t divisionSize;
		uint32_t mismatchesPerDivision;
		uint32_t arrangementsPerDivision;
		uint64_t bucketsPerArrangement;
		uint32_t totalHashmaps;
	};
	public:
	struct offtarget
	{
		uint32_t seqID;
		size_t positionInSeq;
		DNA4 offTargetSequence;
		uint32_t mismatches;
		bool exactMatchTo(DNA4 t)const {return offTargetSequence.getSimilarity(t)==DNA4::getLength();}
		bool operator!=(const offtarget &o) { return o.seqID!=seqID||o.positionInSeq!=positionInSeq;}
		bool operator==(const offtarget &o) { return !operator!=(o);}
		std::string asString() const {return offTargetSequence.toString();}
	};

	OffTargetFinder(const std::vector<DNA4> &targets, uint_fast8_t mismatches, Filter filter, uint64_t maxIndexSize = ~0ULL)
	:offTargets(targets.size()),
		targetContainer(mismatches,filter),
		targets(targets),
		matched1(nullptr),matched2(nullptr),
		bucketsArray(nullptr)
	{
		if(targets.size()==0)
		{
			return;
		}
		for(uint32_t i = 0; i < targets.size();i++)
		{
			offTargets[i].first = targets[i];
		}
		if(filter.exists())
		{
			for(DNA4 &t: this->targets)
			{
				filter.addFilterCharsToTarget(t);
			}
		}
		targetContainer.addTargets(this->targets,maxIndexSize);
		minSimilarity = DNA4::getLength() -  mismatches;
		if(DNA4::getLength() < mismatches) minSimilarity = 0;//avoid underflow

		if(targetContainer)
		{
			std::cout << "using index method" << std::endl;
			bucketsArray = new TargetBucket[targetContainer.numberOfHashmaps()];
		}
		else
		{
			std::cout << "Using simple method"<< std::endl;
			matched1 = new uint32_t[targets.size()/2+1];
			matched2 = new uint32_t[targets.size()/2+1];
		}
	}

	~OffTargetFinder()
	{
		delete[] matched1;
		delete[] matched2;
		delete[] bucketsArray;
	}

	void findMatches(DNA4 & seq, size_t position, uint32_t seqID)
	{
		if(targetContainer)
		{
			findMatchesWithIndex(seq, position, seqID);
		}
		else
		{
			findMatchesSimple(seq, position, seqID);
		}
	}

	std::vector<std::pair<DNA4,std::vector<offtarget>>> & getOffTargets()
	{
		return offTargets;
	}

	private:
	void findMatchesWithIndex(DNA4 & seq, size_t position, uint32_t seqID)
	{
		//naiveComparisons += targetContainer.numberOfTargets();
		targetContainer.getBuckets(seq,bucketsArray);
		for(size_t i = 0; i < targetContainer.numberOfHashmaps(); ++i)
		{
			TargetBucket bucket = bucketsArray[i];
			//std::cout << bucket.size << std::endl;
			for(uint32_t targetNum = 0; targetNum<bucket.size; targetNum++)
			{
				uint32_t similarity = seq.getSimilarity(bucket.begin_DNA[targetNum]);
				
				if(similarity>=minSimilarity)
				{
					//std::cout << similarity <<std::endl;
					auto & targetMatches = offTargets[bucket.begin_position[targetNum]].second;
					//make sure we haven't added this match from another hashmap already
					if(targetMatches.size()==0||targetMatches.back().positionInSeq!=position||targetMatches.back().seqID!=seqID)
					{
						targetMatches.push_back(offtarget{seqID,position,seq, DNA4::getLength()-similarity});
					}
				}
			}
		}
	}

	void findMatchesSimple(DNA4 & seq, size_t position, uint32_t seqID)
	{
		int numberOfMatches1 = 0;
		int numberOfMatches2 = 0;

		//compare all the targets with this sequence
		for(size_t i = 0; i < targets.size()-1;i+=2)
		{
			if(seq.getSimilarity(targets[i])>=minSimilarity)
			{
				matched1[numberOfMatches1] = i;
				numberOfMatches1++;
			}
			if(seq.getSimilarity(targets[i+1])>=minSimilarity)
			{
				matched2[numberOfMatches2] = i+1;
				numberOfMatches2++;
			}
		}
		if(targets.size()%2)
		{
			if(seq.getSimilarity(targets[targets.size()-1])>=minSimilarity)
			{
				matched1[numberOfMatches1] = targets.size()-1;
				numberOfMatches1++;
			}
		}

		for(int i = 0; i < numberOfMatches1;++i)
		{
			offTargets[matched1[i]].second.push_back(offtarget{seqID,position,seq,DNA4::getLength()-seq.getSimilarity(targets[matched1[i]])});
		}	
		for(int i = 0; i < numberOfMatches2;++i)
		{
			offTargets[matched2[i]].second.push_back(offtarget{seqID,position,seq,DNA4::getLength()-seq.getSimilarity(targets[matched2[i]])});
		}
	}

	std::vector<std::pair<DNA4,std::vector<offtarget>>> offTargets;
	TargetContainer targetContainer;
	std::vector<DNA4> targets;
	uint32_t * matched1;
	uint32_t * matched2;
	TargetBucket * bucketsArray;
	uint_fast8_t minSimilarity;
};

struct FindIfOffTarget{
	OffTargetFinder &offTargetFinder;
	uint32_t sequenceID;
	FindIfOffTarget(OffTargetFinder &offTargetFinder,uint32_t sequenceID)
	:offTargetFinder(offTargetFinder),sequenceID(sequenceID)
	{}
	void doAction(DNA4 &t,size_t position)
	{
		offTargetFinder.findMatches(t,position,sequenceID);
	}
};

struct AddToSet{
	DNA4Set &s;
	AddToSet(DNA4Set & set)
	:s(set){}
	void doAction(DNA4 &t,size_t){s.insert(t);}
};

struct RemoveFromSet{
	DNA4Set &s;
	RemoveFromSet(DNA4Set & set)
	:s(set){}
	void doAction(DNA4 &t,size_t){s.remove(t);}
};

struct AddToSetIfExistingInOtherSet{
	DNA4Set &s;
	DNA4Set &otherSet;
	AddToSetIfExistingInOtherSet(DNA4Set &set, DNA4Set &otherSet)
	:s(set),otherSet(otherSet){}
	void doAction(DNA4 &t,size_t)
	{
		if(otherSet.contains(t))
			s.insert(t);
	}
};

template<class Action>
void doForTargetsInSequence(const char * sequence, Filter filter, Action &&action)
{
	for(size_t i = 0; i < DNA4::getLength();++i)
	{
		if(sequence[i]=='\0')
		{
		return;
		}
	}
	DNA4 sequenceDNA4 = DNA4(sequence);
	sequence+=DNA4::getLength();
	size_t currentPos = 0;
	do
	{
		//check that it matches the filter
		if(filter.passes(sequenceDNA4))
		{
			action.doAction(sequenceDNA4,currentPos);
		}
		++currentPos;
	}while(sequenceDNA4.addCharacter(*(sequence++)));
		}

template<class Action>
void doForTargetsInSequence(const std::string &sequence, Filter filter, Action &&action)
{
	doForTargetsInSequence(sequence.data(),filter,action);
}

template<class Action>
void doForTargetsInSequence(std::istream &seqStream, Filter filter, Action &&action)
{
	constexpr size_t buffSize = 4096;
	char buf[buffSize+1];

	seqStream.read(buf,buffSize);
	if(seqStream.gcount()<DNA4::getLength())
		return;
	buf[seqStream.gcount()] = '\0';
	DNA4 sequenceDNA4 = DNA4(buf);
	size_t currentPos = 0;
	char *c = &buf[DNA4::getLength()-1];
	do
	{	
		do
		{
			//check that it matches the filter
			if(filter.passes(sequenceDNA4))
			{
				action.doAction(sequenceDNA4,currentPos);
			}
			currentPos++;
		}while(*(++c)&&sequenceDNA4.addCharacter(*c));
		c = buf;
		seqStream.read(buf,buffSize);
		buf[seqStream.gcount()] = '\0';
		//adds character since this was skipped above via short circuiting
		//if gcount is 0 null will be added but we exit if this is the case anyway
		sequenceDNA4.addCharacter(*c);
	} while (seqStream.gcount()>0);
}