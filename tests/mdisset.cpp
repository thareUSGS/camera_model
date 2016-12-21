//MDIS SET STUB
#include <iostream>
#include <fstream>
#include <string>
#include "Isd.h"

using namespace std;

int main(int argc, char *argv[]){
  if(argc<2){
    cerr << "Usage: " << argv[0] << " <isd file>" << endl;
    return 1;
  }

  csm::Isd isd;

  //Read the ISD file
  string line;
  string filename(argv[1]);
  ifstream f (argv[1]);

      if (!f.is_open()) {
          perror(("error while opening file " + filename).c_str());
      }
      while (getline(f, line)) {
          cout << line << endl;
          }
      if (f.bad()) {
          perror(("error while reading file " + filename).c_str());
      }
      f.close();


  return 0;
}
