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
	std::pair<int, int>	CmpStrWithLCP(const char* a, const char* b, const int& begin)	{
		int la = strlen(a);
		int lb = strlen(b);
		// DEBUG
		//if(begin > la)	{
		//	std::cout << "begin:	" << begin << std::endl;
		//	std::cout << "length a:	" << la << std::endl;
		//	std::cout << "length b:	" << lb << std::endl; 
		//}
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

	// similar purpose as CmpStrWithLCP, but in a reversed direction
	std::pair<int, int> CmpWithLCPRev(const char *a, const char *b, const int& begin)	{
		assert(begin >= 0);
		
		int la = strlen(a);
		int lb = strlen(b);

		std::pair<int, int> r;
		r.second = 0;
		// DEBUG
		//std::cout << a << std::endl;
		//std::cout << b << std::endl;
		//std::cout << la << "	" << lb << std::endl;
		//std::cout << la - begin - 1 - r.second << "	" << lb - begin - 1 - r.second << std::endl;
		//std::cout << a[la - begin - 1 - r.second] << "	" << b[lb - begin - 1 - r.second] << std::endl;
		while((la - begin - 1 - r.second >= 0) && (lb - begin - 1 - r.second >= 0) && a[la - begin - 1 - r.second] == b[lb - begin - 1 - r.second])	{
			// DEBUG
			//std::cout << a[la - begin - 1 - r.second] << "	" << b[lb - begin - 1 - r.second] << std::endl;
			++ r.second;
		}

		if(la - begin - 1 - r.second < 0 && lb - begin - 1 - r.second < 0)	{
			r.first = 0;
		}	else if (la - begin - 1 - r.second < 0)	{
			r.first = -1;
		}	else if (lb - begin - 1 - r.second < 0)	{
			r.first = 1;
		}	else	{
			r.first = a[la - begin - 1 - r.second] < b[lb - begin - 1 - r.second] ? -1 : 1;
		}

		return r;
	}

};

#endif