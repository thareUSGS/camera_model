#ifndef SpiceController_h
#define SpiceController_h


#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>


using namespace std;

class SpiceController {


public:

  //Default constructor
  SpiceController():m_kernlist(),m_furnish(true){ }

  /** Returns the number of kernels found and/or loaded */
  int size() const {
    return (m_kernlist.size());
  }

  void load();

  virtual ~SpiceController() {
    unload();
  }

  void loadKernel(string &kernelFile);

private:
  vector<string> m_kernlist;  //!< The list of kernels
  bool m_furnish;             //!< Load the kernels found?
  void unload();

};

#endif

