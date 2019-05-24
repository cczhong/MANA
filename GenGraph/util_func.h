#ifndef _UTIL_FUNC_H_
#define _UTIL_FUNC_H_

#include <iostream>
#include <string>
#include <cstring>

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

};
#endif