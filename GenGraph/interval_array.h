/** 
 * @file       IntervalArray.h
 * @brief      Quick interval lookup
 * @details    Internal LCP implementation
 * @date       Modified on Tue 2013-12-17 06:46:57 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __INTERVAL_ARRAY_H__
#define __INTERVAL_ARRAY_H__

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <stack>
#include <boost/unordered_map.hpp>
#include "file.h"

typedef uint16_t VALUETYPE;
typedef boost::unordered_map<size_t, size_t> BigValueMap;

const VALUETYPE MAX_INTERVAL_VALUE = 255;

/**
 * \brief Quick interval lookup for internal LCP values.
 */
class IntervalArray
{
 private:
	size_t size;
	size_t leaf;

	VALUETYPE *leaf_ints;
	VALUETYPE *intervals;

 private:


	long getIndex( size_t l , size_t r );
	VALUETYPE fill( size_t l, size_t r );

 public:
	IntervalArray();
	IntervalArray(VALUETYPE *ints, size_t n);
	~IntervalArray();
	void build(VALUETYPE *ints, size_t n);
	void clear();
	void print( std::ostream &out );
	size_t   getSize() { return size; }
	VALUETYPE* getIntervals() { return intervals; }
	VALUETYPE  getValue(size_t i, size_t j);

	void dump( const char *filename );
	void load( const char *filename, size_t filesize );
};


#endif

