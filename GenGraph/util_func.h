#ifndef _UTIL_FUNC_H_
#define _UTIL_FUNC_H_

#include <iostream>
#include <string>
#include <cstring>
#include <tuple>

class UtilFunc {
public:

    UtilFunc()  {}
    ~UtilFunc() {}

    std::string GetFileStem(const std::string& path)  {
  		// going backward untill the '\/' character
  		int i;
  		for(i = path.length() - 1; i >= 0; -- i) {
    		if(path[i] == '/')  {
      			break;
    		}
  		}
  		return path.substr(i + 1, path.length() - i - 1);
	}

	// The first return value is -1 if a is lexicographically less than b
	//						  is 0 if a is equal to b
	//						  is 1 if a is larger than b
	// The second return value indicates Longest Common Prefix length of the two strings
	// (assuming letters before "begin" are identical)
	std::pair<int, int>	CmpStrWithLCP(char* a, char* b, const int& begin)	{
		int la = strlen(a);
		int lb = strlen(b);
		assert(begin <= la);
		assert(begin <= lb);
		std::pair<int, int> r;
		r.second = begin;
		while(r.second < la && r.second < lb && a[r.second] == b[r.second])	{
			++ r.second;
		}
		if(a[r.second] < b[r.second] || (r.second == la && r.second < lb))	{
			r.first = -1;
		}	else if (a[r.second] > b[r.second] || (r.second == lb && r.second < la))	{
			r.first = 1;
		}	else {
			r.first = 0;
		}
		return r;
	}

};

#endif