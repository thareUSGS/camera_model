#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <SpiceUsr.h>

#include "SpiceController.h"


//Compiling class
//g++ *.cpp -c -Wall -Wextra -I../../../include -L../../../lib/

using namespace std;


void SpiceController::unload() {

  for (unsigned int i = 0; i < m_kernlist.size();i++) {

    unload_c(m_kernlist[i].c_str() );

  }

}




void SpiceController::load() {

      cout <<"Loading Spice kernels" << endl;


}

void SpiceController::loadKernel(string &kernelFile) {



  fstream inFile(kernelFile.c_str());

  if(!inFile.good() ) {
     cout << "Could not load: " << kernelFile << endl;
     return;
  }

  if(m_furnish) {
    furnsh_c(kernelFile.c_str() );
    m_kernlist.push_back(kernelFile);
    cout << kernelFile << " loaded."  <<endl;
  }






}



