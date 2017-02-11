#include "utils.h"
#include <sys/time.h>

namespace pr {
  //! returns the current time in milliseconds
  double getTime(){
    struct timeval tv;
    gettimeofday(&tv, 0);
    return 1e3*tv.tv_sec+tv.tv_usec*1e-3;
  } 

  //! prints a null terminated array od strings
  void printBanner(const char** banner, std::ostream& os){
    while (*banner) {
      os << *banner << std::endl;
      banner++;
    }
  }
}
