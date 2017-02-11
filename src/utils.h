#pragma once
#include <iostream>
namespace pr {
  //! returns the current time in milliseconds
  double getTime(); 

  //! prints a null terminated array od strings
  void printBanner(const char** banner, std::ostream& os=std::cout);
}
