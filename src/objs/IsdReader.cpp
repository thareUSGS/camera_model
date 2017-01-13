#include <IsdReader.h>

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>

#include <json/json.hpp>
#include <csm/Isd.h>

using namespace std;
using json = nlohmann::json;

csm::Isd *readISD(string filename) {
  json jsonFile;
  int prec = 12;
  csm::Isd *isd = NULL;

  //Read the ISD file
  ifstream file(filename);
  if (!file.is_open()) {
    perror(("error while opening file " + filename).c_str());
    return NULL;
  }

  else if (file.bad()) {
    perror(("error while reading file " + filename).c_str());
    return NULL;
  }

  else {
    isd = new csm::Isd();
    file >> jsonFile;
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
  } //end outer-else
  //printISD(*isd);
  file.close();
  return isd;
}


/**
 * @brief checkType
 * @param obj
 * @return An enum DataType value indicating what the primitive data type of obj is
 * @author Tyler Wilson
 */
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


/**
 * @brief printISD  Display the keyword:value pairs of a CSM::ISD object
 * @param isd
 * @author Tyler Wilson
 */
void printISD(const csm::Isd &isd){
  const multimap<string,string> isdMap = isd.parameters();
  for (auto &i: isdMap)
    cout << i.first << " : " << i.second << endl;
}


