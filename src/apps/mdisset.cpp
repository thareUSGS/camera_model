//MDIS SET STUB
#include <iostream>
#include "Isd.h"

int main(int argc, char *argv[]){
  if(argc<3){
    std::cerr << "Usage: " << argv[0] << " <isd file>" << std::endl;
    return 1;
  }

  return 0;
}
