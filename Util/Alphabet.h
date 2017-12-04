//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Alphabet.h - Abstraction of the alphabet
// that is used in the suffix array, fm-index,
// etc data structures
//
#ifndef ALPHABET_H
#define ALPHABET_H

//#define USE_SSE 1
//#define ALPHACOUNT_VALIDATE 1
#include <utility>
#include <stdint.h>
#include <limits>
#include <math.h>

#if USE_SSE
#include <emmintrin.h>
#include <xmmintrin.h>
#endif

#include <iostream>
#include <iterator>
#include <algorithm>
#include "Util.h"

//
// Constants
//
// TODO: Add amino acid alphabet here
// affected files
// Utils/Pileup.cpp
// SGA/kmer-count.cpp
// DNA parameters
const uint8_t DNA_ALPHABET_SIZE = 5;
const char DNA_RAW_ALPHABET[DNA_ALPHABET_SIZE] = {'A', 'C', 'G', 'T', '$'};
const char DNA_RANK_ALPHABET[DNA_ALPHABET_SIZE] = {'$', 'A', 'C', 'G', 'T'};
const uint8_t DNA_ALPHABET_DEDUCED_SIZE = 4;
// Amino acids parameters
const uint8_t AA_ALPHABET_SIZE = 21;
const char AA_RAW_ALPHABET[AA_ALPHABET_SIZE] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '$'};
const char AA_RANK_ALPHABET[AA_ALPHABET_SIZE] = {'$', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
const uint8_t AA_ALPHABET_REDUCED_SIZE = 20;
typedef uint64_t BaseCount;

/*
class DNAConst  {
  // TODO: make these variable private and return with functions
 public: 
  DNAConst()  {
    ALPHABET_SIZE = 5;
    ALPHABET_DEDUCED_SIZE = 4;
    char alphabet[5] = {'A', 'C', 'G', 'T', '$'};
    for(int i = 0; i < ALPHABET_SIZE; ++ i)
      ALPHABET[i] = alphabet [i];
    char rankAlphabet[5] = {'$', 'A', 'C', 'G', 'T'};    
    for(int i = 0; i < ALPHABET_SIZE; ++ i)
      RANK_ALPHABET[i] = rankAlphabet [i];
  }
  static uint8_t ALPHABET_SIZE;
  static uint8_t ALPHABET_DEDUCED_SIZE ;
  static char ALPHABET[5];
  static char RANK_ALPHABET[5];  
};

class AAConst  {
  // TODO: make these variable private and return with functions
 public: 
  AAConst()  {
    ALPHABET_SIZE = 21;
    ALPHABET_DEDUCED_SIZE = 20;
    char alphabet[21] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '$'};
    for(int i = 0; i < ALPHABET_SIZE; ++ i)
      ALPHABET[i] = alphabet [i];
    char rankAlphabet[21] = {'$', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
    for(int i = 0; i < ALPHABET_SIZE; ++ i)
      RANK_ALPHABET[i] = rankAlphabet [i];
  }
  static uint8_t ALPHABET_SIZE;
  static uint8_t ALPHABET_DEDUCED_SIZE ;
  static char ALPHABET[21];
  static char RANK_ALPHABET[21];  
};
*/



namespace DNA_ALPHABET
{
    static const uint8_t s_dnaLexoRankLUT[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
        0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };

    inline static uint8_t getBaseRank(char b)
    {
        return s_dnaLexoRankLUT[static_cast<uint8_t>(b)];
    }
    
    inline char getBase(size_t idx)
    {
        assert(idx < DNA_ALPHABET_DEDUCED_SIZE);
        return DNA_RAW_ALPHABET[idx];
    }

    static const uint8_t size = 4;
};

namespace BWT_DNA_ALPHABET
{
    static const uint8_t s_bwtLexoRankLUT[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
        0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };

    static const uint8_t size = 5;

    inline static uint8_t getRank(char b)
    {
        return s_bwtLexoRankLUT[static_cast<uint8_t>(b)];
    }
    
    inline char getChar(size_t idx)
    {
        return DNA_RANK_ALPHABET[idx];
    }
};

// IUPAC ambiguity alphabet
namespace IUPAC
{
    // Returns true if c is [ACGT]
    bool isUnambiguous(char c);

    // Returns true if c is a valid ambiguity code
    bool isAmbiguous(char c);

    // Returns true if c is a valid symbol in this alphabet
    bool isValid(char c);

    // Returns a string defining the possible unambiguous bases for each symbol
    // in the alphabet
    std::string getPossibleSymbols(char c);
};

// IUPAC amino acid ambiguity alphabet
namespace IUPAC_AA
{
    // Returns true if c is [ARNDCQEGHILKMFPSTWYV]
    bool isUnambiguous(char c);

    // Returns true if c is a valid ambiguity code
    bool isAmbiguous(char c);

    // Returns true if c is a valid symbol in this alphabet
    bool isValid(char c);

    // Returns a string defining the possible unambiguous bases for each symbol
    // in the alphabet
    std::string getPossibleSymbols(char c);
};

//
// A simple class holding the count for each base of a DNA string (plus the terminator)
// Note that this uses RANK_ALPHABET
// TODO: has many DNA specific function, check if we need to modify them
template<typename Storage, const bool isAminoAcid>
class AlphaCount
{
    public:
        //
        inline AlphaCount()
        {
            clear();
        }

        inline void clear()
        {
            memset(m_counts, 0, (isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE) * sizeof(Storage));
        }

        //
        inline void set(char b, Storage v)
        {
            m_counts[getBaseRank(b)] = v;
        }

        //
        inline void setByIdx(size_t i, Storage v)
        {
            m_counts[i] = v;
        }


        // 
        inline void increment(char b)
        {
            int br = getBaseRank(b);
#ifdef ALPHACOUNT_VALIDATE
            assert(m_counts[br] != maxValue);
#endif
            m_counts[br]++;
        }

        //
        inline void add(char b, Storage v)
        {
            int br = getBaseRank(b);
#ifdef ALPHACOUNT_VALIDATE
            assert(m_counts[br] + v <= maxValue);
#endif
            m_counts[br] += v;
        }

        //
        inline void subtract(char b, Storage v)
        {
            int br = getBaseRank(b);
#ifdef ALPHACOUNT_VALIDATE
            assert(m_counts[br] > v);
#endif            
            m_counts[br] -= v;
        }

        // 
        inline Storage get(char b) const
        {
            return m_counts[getBaseRank(b)];
        }

        //
        inline Storage getByIdx(const int i) const
        {
            return m_counts[i];
        }

        // Return the base for index i
        static char getBase(size_t i)
        {
            return DNA_RANK_ALPHABET[i];
        }
        
        // Return the maximum possible count for a symbol
        static size_t getMaxValue()
        {
            return maxValue;
        }

        // Swap the (A,T) and (C,G) entries, which turns the AlphaCount
        // into the AlphaCount for the complemented sequence
        inline void complement()
        {
            Storage tmp;

            // A,T
            tmp = m_counts[4];
            m_counts[4] = m_counts[1];
            m_counts[1] = tmp;

            // C,G
            tmp = m_counts[3];
            m_counts[3] = m_counts[2];
            m_counts[2] = tmp;
        }

        // Return the sum of the basecounts for characters lexo. lower than b
        inline size_t getLessThan(char b) const
        {
            size_t out = 0;
            int stop = getBaseRank(b);
            for(int i = 0; i < stop; ++i)
                out += m_counts[i];
            return out;
        }

        // Returns the number of non-zero counts
        uint8_t getNumNonZero() const
        {
            uint8_t count = 0;
            for(int i = 0; i < (isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE); ++i)
            {
                if(m_counts[i] > 0)
                    count += 1;
            }
            return count;
        }

        // Returns true if only one of the DNA characters
        // has a non-zero count
        inline bool hasUniqueDNAChar()
        {
            // index 0 is the '$' character, which we skip
            bool nonzero = false;
            for(int i = 1; i < (isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE); ++i)
            {
                if(m_counts[i] > 0)
                {
                    // if nonzero is set, there is some other nonzero character, return false
                    if(nonzero)
                        return false;
                    else
                        nonzero = true;
                }
            }
            return nonzero;
        }

        // Returns true if any DNA char 
        // has a non-zero count
        inline bool hasDNAChar()
        {
            // index 0 is the '$' character, which we skip
            for(int i = 1; i < (isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE); ++i)
            {
                if(m_counts[i] > 0)
                    return true;
            }
            return false;
        }

        //
        inline char getMaxBase() const
        {
            char base;
            Storage val;
            getMax(base, val);
            return base;
        }

        // 
        inline char getMaxDNABase() const
        {
            Storage max = 0;
            int maxIdx = 0;
            for(int i = 1; i < DNA_ALPHABET_SIZE; ++i)
            {
                if(m_counts[i] > max)
                {
                    max = m_counts[i];
                    maxIdx = i;
                }
            }

            assert(max > 0);
            return DNA_RANK_ALPHABET[maxIdx];
        }
        //
        inline Storage getMaxCount() const
        {
            char base;
            Storage val;
            getMax(base, val);
            return val;
        }

        //
        inline void getMax(char& base, Storage& val) const
        {
            Storage max = 0;
            int maxIdx = 0;
            for(int i = 0; i < isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE; ++i)
            {
                if(m_counts[i] > max)
                {
                    max = m_counts[i];
                    maxIdx = i;
                }
            }
            base = isAminoAcid ? AA_RANK_ALPHABET[maxIdx] : DNA_RANK_ALPHABET[maxIdx];
            val = max;
        }

        //
        inline size_t getSum() const
        {
            size_t sum = m_counts[0];
            sum += m_counts[1];
            sum += m_counts[2];
            sum += m_counts[3];
            sum += m_counts[4];
            return sum;
        }

        // Sort the DNA bases into count order and
        // write them to pOut. This does not include the '$' symbol
        inline void getSorted(char* pOut, size_t len) const
        {
            uint8_t alphabetSize = isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE;
            assert(len >= alphabetSize);
            (void)len;
            for(size_t i = 0; i < alphabetSize; ++i)
                pOut[i] = isAminoAcid ? AA_RANK_ALPHABET[i] : DNA_RANK_ALPHABET[i];
            std::sort(pOut, pOut+alphabetSize, AlphaCountCompareDesc(this));
        }

        // TODO: check whether we need a generalized getUniqueChar function
        // Return the unique DNA character described by the alphacount
        // Returns '$' if no such character exists
        // Asserts if more than one character is described by the alphacount
        inline char getUniqueDNAChar()
        {
            char r = '$';
            for(int i = 1; i < DNA_ALPHABET_SIZE; ++i)
            {
                if(m_counts[i] > 0)
                {
                    assert(r == '$');
                    // lookup the character in the ranked alphabet, since the elements
                    // are stored in lexographic order
                    r = DNA_RANK_ALPHABET[i];
                }
            }
            return r;
        }

        // Operators
        friend std::ostream& operator<<(std::ostream& out, const AlphaCount<Storage, isAminoAcid>& ac)
        {
            uint8_t alphabetSize = isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE;
            std::copy(ac.m_counts, ac.m_counts+alphabetSize, std::ostream_iterator<Storage>(out, " "));
            return out;
        }

        friend std::istream& operator>>(std::istream& in, AlphaCount<Storage, isAminoAcid>& ac)
        {
            for(size_t i = 0; i < isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE; ++i)
                in >> ac.m_counts[i];
            return in;
        }
        
        inline friend AlphaCount operator+(const AlphaCount<Storage, isAminoAcid>& left, const AlphaCount<Storage, isAminoAcid>& right)
        {
            AlphaCount<Storage, isAminoAcid> out;
            for(int i = 0; i < isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE; ++ i)
              out.m_counts[i] = left.m_counts[i] + right.m_counts[i];
            /*
            out.m_counts[0] = left.m_counts[0] + right.m_counts[0];
            out.m_counts[1] = left.m_counts[1] + right.m_counts[1];
            out.m_counts[2] = left.m_counts[2] + right.m_counts[2];
            out.m_counts[3] = left.m_counts[3] + right.m_counts[3];
            out.m_counts[4] = left.m_counts[4] + right.m_counts[4];
            */
            return out;
        }

        inline AlphaCount& operator+=(const AlphaCount<Storage, isAminoAcid>& other)
        {
            for(int i = 0; i < isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE; ++ i)
              m_counts[i] = other.m_counts[i];
            /*
            m_counts[0] += other.m_counts[0];
            m_counts[1] += other.m_counts[1];
            m_counts[2] += other.m_counts[2];
            m_counts[3] += other.m_counts[3];
            m_counts[4] += other.m_counts[4];
            */
            return *this;
        }

        // Specialization for add/subtracing a small AlphaCount to/from a large alphacount
        inline friend void alphacount_add(AlphaCount<uint64_t, 0>& lhs, const AlphaCount<uint8_t, 0>& rhs);
        inline friend void alphacount_subtract(AlphaCount<uint64_t, 0>& lhs, const AlphaCount<uint8_t, 0>& rhs);
        inline friend void alphacount_add16(AlphaCount<uint64_t, 0>& lhs, const AlphaCount<uint16_t, 0>& rhs);
        inline friend void alphacount_subtract16(AlphaCount<uint64_t, 0>& lhs, const AlphaCount<uint16_t, 0>& rhs);
        // the amino acid equivalent
        inline friend void alphacount_add(AlphaCount<uint64_t, 1>& lhs, const AlphaCount<uint8_t, 1>& rhs);
        inline friend void alphacount_subtract(AlphaCount<uint64_t, 1>& lhs, const AlphaCount<uint8_t, 1>& rhs);
        inline friend void alphacount_add16(AlphaCount<uint64_t, 1>& lhs, const AlphaCount<uint16_t, 1>& rhs);
        inline friend void alphacount_subtract16(AlphaCount<uint64_t, 1>& lhs, const AlphaCount<uint16_t, 1>& rhs);
        

        // As the counts are unsigned integers, each value in left
        // must be larger or equal to value in right. The calling function
        // must guarentee this.
        friend AlphaCount operator-(const AlphaCount<Storage, isAminoAcid>& left, const AlphaCount<Storage, isAminoAcid>& right)
        {
            AlphaCount<Storage, isAminoAcid> out;
            for(int i = 0; i < isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE; ++ i)
              out.m_counts[i] = left.m_counts[i] - right.m_counts[i];
            /*
            out.m_counts[0] = left.m_counts[0] - right.m_counts[0];
            out.m_counts[1] = left.m_counts[1] - right.m_counts[1];
            out.m_counts[2] = left.m_counts[2] - right.m_counts[2];
            out.m_counts[3] = left.m_counts[3] - right.m_counts[3];
            out.m_counts[4] = left.m_counts[4] - right.m_counts[4];
            */
            return out;
        }

        inline friend bool operator==(const AlphaCount<Storage, isAminoAcid>& left, const AlphaCount<Storage, isAminoAcid>& right)
        {
            for(int i = 0; i < isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE; ++ i)
            {
              if(left.m_counts[i] != right.m_counts[i]) return false;
            }
            return true;
            /*
            return left.m_counts[0] == right.m_counts[0] && left.m_counts[1] == right.m_counts[1] &&
                   left.m_counts[2] == right.m_counts[2] && left.m_counts[3] == right.m_counts[3] &&
                   left.m_counts[4] == right.m_counts[4];   
            */
        }
        
        inline friend bool operator!=(const AlphaCount<Storage, isAminoAcid>& left, const AlphaCount<Storage, isAminoAcid>& right)
        {
            return !(left == right);
        }
                

    private:
        Storage m_counts[isAminoAcid ? AA_ALPHABET_SIZE : DNA_ALPHABET_SIZE];
        const static size_t maxValue;

        struct AlphaCountCompareDesc
        {
            AlphaCountCompareDesc(const AlphaCount* p) : pAC(p) {}
            bool operator()(char a, char b)
            {
                return pAC->get(a) > pAC->get(b);
            }
            const AlphaCount* pAC;
        };

};

// Typedef commonly used AlphaCounts
typedef AlphaCount<uint64_t, 0> AlphaCount64;
typedef AlphaCount<uint16_t, 0> AlphaCount16;
typedef AlphaCount<uint8_t, 0> AlphaCount8;
typedef AlphaCount<uint64_t, 1> AAAlphaCount64;
typedef AlphaCount<uint16_t, 1> AAAlphaCount16;
typedef AlphaCount<uint8_t, 1> AAAlphaCount8;
//
inline void alphacount_add16(AlphaCount64& lhs, const AlphaCount16& rhs)
{
    for(int i = 0; i < 5; ++ i)
      lhs.m_counts[i] += rhs.m_counts[i];
    /*
    lhs.m_counts[0] += rhs.m_counts[0];
    lhs.m_counts[1] += rhs.m_counts[1];
    lhs.m_counts[2] += rhs.m_counts[2];
    lhs.m_counts[3] += rhs.m_counts[3];
    lhs.m_counts[4] += rhs.m_counts[4];
    */
}

inline void alphacount_subtract16(AlphaCount64& lhs, const AlphaCount16& rhs)
{
    // This function should only be used when lhs is larger the rhs
    // the calling function must guarentee this
    for(int i = 0; i < 5; ++ i)
      lhs.m_counts[i] -= rhs.m_counts[i];
    /*
    lhs.m_counts[0] -= rhs.m_counts[0];
    lhs.m_counts[1] -= rhs.m_counts[1];
    lhs.m_counts[2] -= rhs.m_counts[2];
    lhs.m_counts[3] -= rhs.m_counts[3];
    lhs.m_counts[4] -= rhs.m_counts[4];
    */
}

//
inline void alphacount_add(AlphaCount64& lhs, const AlphaCount8& rhs)
{
    for(int i = 0; i < 5; ++ i)
      lhs.m_counts[i] += rhs.m_counts[i];
    /*    
    lhs.m_counts[0] += rhs.m_counts[0];
    lhs.m_counts[1] += rhs.m_counts[1];
    lhs.m_counts[2] += rhs.m_counts[2];
    lhs.m_counts[3] += rhs.m_counts[3];
    lhs.m_counts[4] += rhs.m_counts[4];
    */
}

inline void alphacount_subtract(AlphaCount64& lhs, const AlphaCount8& rhs)
{
    // This function should only be used when lhs is larger the rhs
    // the calling function must guarentee this
    for(int i = 0; i < 5; ++ i)
      lhs.m_counts[i] -= rhs.m_counts[i];
    /*
    lhs.m_counts[0] -= rhs.m_counts[0];
    lhs.m_counts[1] -= rhs.m_counts[1];
    lhs.m_counts[2] -= rhs.m_counts[2];
    lhs.m_counts[3] -= rhs.m_counts[3];
    lhs.m_counts[4] -= rhs.m_counts[4];
    */
}

// The amino acid equivalent
inline void alphacount_add16(AAAlphaCount64& lhs, const AAAlphaCount16& rhs)
{
    for(int i = 0; i < 21; ++ i)
      lhs.m_counts[i] += rhs.m_counts[i];
    /*
    lhs.m_counts[0] += rhs.m_counts[0];
    lhs.m_counts[1] += rhs.m_counts[1];
    lhs.m_counts[2] += rhs.m_counts[2];
    lhs.m_counts[3] += rhs.m_counts[3];
    lhs.m_counts[4] += rhs.m_counts[4];
    */
}

inline void alphacount_subtract16(AAAlphaCount64& lhs, const AAAlphaCount16& rhs)
{
    // This function should only be used when lhs is larger the rhs
    // the calling function must guarentee this
    for(int i = 0; i < 21; ++ i)
      lhs.m_counts[i] -= rhs.m_counts[i];
    /*
    lhs.m_counts[0] -= rhs.m_counts[0];
    lhs.m_counts[1] -= rhs.m_counts[1];
    lhs.m_counts[2] -= rhs.m_counts[2];
    lhs.m_counts[3] -= rhs.m_counts[3];
    lhs.m_counts[4] -= rhs.m_counts[4];
    */
}

//
inline void alphacount_add(AAAlphaCount64& lhs, const AAAlphaCount8& rhs)
{
    for(int i = 0; i < 21; ++ i)
      lhs.m_counts[i] += rhs.m_counts[i];
    /*
    lhs.m_counts[0] += rhs.m_counts[0];
    lhs.m_counts[1] += rhs.m_counts[1];
    lhs.m_counts[2] += rhs.m_counts[2];
    lhs.m_counts[3] += rhs.m_counts[3];
    lhs.m_counts[4] += rhs.m_counts[4];
    */
}

inline void alphacount_subtract(AAAlphaCount64& lhs, const AAAlphaCount8& rhs)
{
    // This function should only be used when lhs is larger the rhs
    // the calling function must guarentee this
    for(int i = 0; i < 21; ++ i)
      lhs.m_counts[i] -= rhs.m_counts[i];
    /*
    lhs.m_counts[0] -= rhs.m_counts[0];
    lhs.m_counts[1] -= rhs.m_counts[1];
    lhs.m_counts[2] -= rhs.m_counts[2];
    lhs.m_counts[3] -= rhs.m_counts[3];
    lhs.m_counts[4] -= rhs.m_counts[4];
    */
}


#endif
