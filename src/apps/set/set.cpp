//#include "IsdReader.h"
#include "MdisPlugin.h"
#include "MdisNacSensorModel.h"

#include "Isd.h"
#include "csm.h"
#include <json/json.hpp>

#include <string>
#include <fstream>
#include <sstream>



#include <iostream>


using namespace std;
using json = nlohmann::json;


enum DataType {
  INT,
  UINT,
  FLOAT,
  STRING,
  BOOL,
  NULL8,
  UNKNOWN
};



csm::Isd *readISD(string filename);
DataType checkType(json::value_type obj);
void addParam(csm::Isd &isd, json::iterator it, DataType dt, int prec);


int main(int argc,char *argv[]) {


  
  // Create and load in the ISD object from an ISD file
  //csm::Isd isd;

  string filename("json.isd");
  csm::Isd *isd = readISD("json.isd");


#if 0

  
  // Create the plugin
  MdisPlugin plugin;
  
  // Create the Sensor model from the ISD object using the plugin
  //MdisNacSensorModel model =
  //  plugin.constructModelFromISD(&isd, "ISIS_MDISNAC_USGSAstro_1_Linux64_csm30.so");

   MdisNacSensorModel *model =
    (MdisNacSensorModel *)plugin.constructModelFromISD(*isd, "ISIS_MDISNAC_USGSAstro_1_Linux64_csm30.so");




  
  // Test the model's accuracy by visually comparing the results of the function calls
  csm::EcefCoord groundPoint;
  csm::ImageCoord imagePoint;
  
  imagePoint.samp = 300;
  imagePoint.line = 400;
  
  for (int i = 0; i < 10; i++) {
    groundPoint = model->imageToGround(imagePoint, 0);

    std::cout << "Image Point (s, l) : (" << imagePoint.samp << ", " << imagePoint.line << "); "
              << "Ground Point (x, y, z) : (" << groundPoint.x << ", " << groundPoint.y << ", " 
                                              << groundPoint.z << ")" << endl;

    imagePoint = model->groundToImage(groundPoint);
  }
 #endif
  //return 0;
}

csm::Isd *readISD(string filename) {
  //json jsonFile;

  int prec = 12;
  csm::Isd *isd = NULL;


  //Read the ISD file
  stringstream file(filename);
/*
  if (!file.is_open()) {
    perror(("error while opening file " + filename).c_str());
    return NULL;
  }

  else if (file.bad()) {
    perror(("error while reading file " + filename).c_str());
    return NULL;
  }
*/
  //else {
    isd = new csm::Isd();

    json jsonFile(file);
    //file >> jsonFile;

#if 0
    isd->setFilename(filename);

    for (json::iterator i = jsonFile.begin(); i != jsonFile.end(); i++) {
        if (i.value().is_array()){
          DataType arrayType = checkType(i.value()[0]);
          addParam(*isd, i,arrayType,prec);
        }
        else {
          DataType dt = checkType(i.value());
          addParam(*isd,i,dt,prec);
        }
    }//end for
#endif
  //} //end outer-else
  //printISD(*isd);
  //file.close();

  return isd;
}

DataType checkType(json::value_type obj){

  if (obj.is_number()) {
    if (obj.is_number_float())
      return FLOAT;
    else if (obj.is_number_integer())
      return INT;
    else if(obj.is_number_unsigned())
      return UINT;
    else
      return UNKNOWN;
  }
  else if(obj.is_null())
    return NULL8;

  else if(obj.is_string())
    return STRING;

  else if(obj.is_boolean())
    return BOOL;
  else
    return UNKNOWN;

}


/**
 * @brief addParam:
 * @param isd A reference to the ISD object
 * @param it  The iterator to the json file which iterates over the keywords.
 * @param dt  The enum DataType value
 * @param prec The # of decimal places to be written to the ISD (if the value is a float)
 * @author Tyler Wilson
 */
void addParam(csm::Isd &isd, json::iterator it, DataType dt, int prec) {
  ostringstream key;
  //output the key to the ISD
  key << it.key();
  if (it.value().is_array()) {
    if (dt==FLOAT) {
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val << setprecision(prec) << v[j];
        isd.addParam(key.str(),val.str());
      }
    }
    else if(dt==INT){
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val << v[j];
        isd.addParam(key.str(),val.str());
      }
    }

    else if(dt==UINT) {
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val  << v[j];
        isd.addParam(key.str(),val.str());
      }
    }
    else if(dt ==BOOL) {
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val << v[j];
        isd.addParam(key.str(),val.str());
      }
    }
    else if (dt ==STRING) {
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val <<  v[j];
        isd.addParam(key.str(),val.str());
      }
    }
  }
  else {
    if(dt==FLOAT) {
      double v = it.value();
      ostringstream val;
      val << setprecision(prec) << v;
      isd.addParam(key.str(),val.str());
    }
    else if(dt==INT){
      int v = it.value();
      ostringstream val;
      val  << v;
      isd.addParam(key.str(),val.str());
    }
    else if(dt==UINT) {
      unsigned int v = it.value();
      ostringstream val;
      val << v;
      isd.addParam(key.str(),val.str());
    }
    else if(dt ==BOOL) {
      bool v = it.value();
      ostringstream val;
      val  << v;
      isd.addParam(key.str(),val.str());
    }
    else if (dt ==STRING) {
      string v = it.value();
      ostringstream val;
      val  << v;
      isd.addParam(key.str(),val.str());
    }
    else if (dt ==NULL8) {
      ostringstream val;
      val  << "null";
      isd.addParam(key.str(),val.str());
    }
  }//end outer else
}//end addParam
