#include "MdisPlugin.h"
#include "MdisNacSensorModel.h"

#include <Isd.h>
#include <csm.h>

#include <iostream>

using namespace std;

int main() {
  
  csm::Isd isd;
  
  // TODO: create the ISD obj from an ISD file
  
  MdisPlugin plugin;
  MdisNacSensorModel* model;
  
  // Initialize the MdisNacSensorModel from the ISD
  try {
    model = (MdisNacSensorModel*)
            (plugin.constructModelFromISD(isd, "ISIS_MDISNAC_USGSAstro_1_Linux64_csm30.so"));
  }
  catch (csm::Error& e) {
    std::cout << e.what() << endl;
    exit (EXIT_FAILURE);
  }
  
  exit (EXIT_SUCCESS);
}