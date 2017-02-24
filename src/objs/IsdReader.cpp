#include <IsdReader.h>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include <csm/Isd.h>
#include <json/json.hpp>

using namespace std;
using json = nlohmann::json;

/**
 * @internal
 *   @todo This should be converted to a C++ class.
 */

/**
 * Reads in a JSON formatted file, parses it, and creates a new csm::Isd in memory.
 * 
 * @param filename JSON file to read.
 * 
 * @return @b csm::Isd* Returns a pointer to the dynamically allocated csm::Isd.
 *                      Returns a null pointer if unsuccessful.
 */
csm::Isd *readISD(string filename) {
  json jsonFile;
  int prec = 15;
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

  // File successfully opened
  else {
    isd = new csm::Isd();
    file >> jsonFile;
    isd->setFilename(filename);
    
    // Parse the JSON and populate the ISD
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
 * Checks the type of the currently parsed JSON token.
 * 
 * @param obj The json::value_type currently being parsed.
 * 
 * @return An enum DataType value indicating what the primitive data type of obj is.
 * 
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
 * Adds a parameter and its value to the ISD object being created.
 * 
 * @param isd A reference to the ISD object being created.
 * @param it  The iterator to the json file which iterates over the keywords.
 * @param dt  The enum DataType value of the value.
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
 * Prints the ISD to standard output.
 * 
 * @brief printISD Display the keyword:value pairs of a CSM::ISD object
 * @param isd Reference to the ISD to output.
 * @author Tyler Wilson
 */
void printISD(const csm::Isd &isd){
  cout.precision(15);
  const multimap<string,string> isdMap = isd.parameters();
  for (auto &i: isdMap)
    cout << i.first << " : " << i.second << endl;
}


