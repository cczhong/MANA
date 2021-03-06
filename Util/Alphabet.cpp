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
#include "Alphabet.h"
#include <iostream>
#include <iterator>

// IUPAC alphabet
bool IUPAC::isUnambiguous(char c)
{
    switch(c)
    {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return true;
        default:
            return false;
    }
}

// Returns true if c is a valid ambiguity code
bool IUPAC::isAmbiguous(char c)
{
    switch(c)
    {
        case 'M':
        case 'R':
        case 'W':
        case 'S':
        case 'Y':
        case 'K':
        case 'V':
        case 'H':
        case 'D':
        case 'B':
        case 'N':
            return true;
        default:
            return false;
    }
}

// Returns true if c is a valid symbol in this alphabet
bool IUPAC::isValid(char c)
{
    return IUPAC::isUnambiguous(c) || IUPAC::isAmbiguous(c);
}

//
std::string IUPAC::getPossibleSymbols(char c)
{
    switch(c)
    {
        case 'A':
            return "A";
        case 'C':
            return "C";
        case 'G':
            return "G";
        case 'T':
            return "T";
        case 'M':
            return "AC";
        case 'R':
            return "AG";
        case 'W':
            return "AT";
        case 'S':
            return "CG";
        case 'Y':
            return "CT";
        case 'K':
            return "GT";
        case 'V':
            return "ACG";
        case 'H':
            return "ACT";
        case 'D':
            return "AGT";
        case 'B':
            return "CGT";
        case 'N':
            return "ACGT";
        default:
            return "";
    }
}

// IUPAC_AA alphabet
bool IUPAC_AA::isUnambiguous(char c)
{

    switch(c)
    {
        case 'A':
        case 'R':
        case 'N':
        case 'D':
        case 'C':
        case 'Q':
        case 'E':
        case 'G':
        case 'H':
        case 'I':
        case 'L':
        case 'K':
        case 'M':
        case 'F':
        case 'P':
        case 'S':
        case 'T':
        case 'W':
        case 'Y':
        case 'V':
            return true;
        default:
            return false;
    }
}

// Returns true if c is a valid ambiguity code
bool IUPAC_AA::isAmbiguous(char c)
{
    switch(c)
    {
        case 'B':
        case 'Z':
        case 'X':
            return true;
        default:
            return false;
    }
}

// Returns true if c is a valid symbol in this alphabet
bool IUPAC_AA::isValid(char c)
{
    return IUPAC_AA::isUnambiguous(c) || IUPAC_AA::isAmbiguous(c);
}

//
std::string IUPAC_AA::getPossibleSymbols(char c)
{
    switch(c)
    {
        case 'A':
            return "A";
        case 'R':
            return "R";
        case 'N':
            return "N";
        case 'D':
            return "D";
        case 'C':
            return "C";
        case 'Q':
            return "Q";
        case 'E':
            return "E";
        case 'G':
            return "G";
        case 'H':
            return "H";
        case 'I':
            return "I";
        case 'L':
            return "L";
        case 'K':
            return "K";
        case 'M':
            return "M";
        case 'F':
            return "F";
        case 'P':
            return "P";
        case 'S':
            return "S";
        case 'T':
            return "T";
        case 'W':
            return "W";
        case 'Y':
            return "Y";
        case 'V':
            return "V";
        case 'B':
            return "DN";
        case 'Z':
            return "QE";
        case 'X':
            return "ARNDCQEGHILKMFPSTWYV";
        default:
            return "";
    }
}

// AlphaCount
template<> const size_t AlphaCount8::maxValue = std::numeric_limits<uint8_t>::max();
template<> const size_t AlphaCount16::maxValue = std::numeric_limits<uint16_t>::max();
template<> const size_t AlphaCount64::maxValue = std::numeric_limits<uint64_t>::max();
template<> const size_t AAAlphaCount8::maxValue = std::numeric_limits<uint8_t>::max();
template<> const size_t AAAlphaCount16::maxValue = std::numeric_limits<uint16_t>::max();
template<> const size_t AAAlphaCount64::maxValue = std::numeric_limits<uint64_t>::max();
