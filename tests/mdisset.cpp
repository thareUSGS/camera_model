//Compile line:
//This needs to incorporated into a Makefile at some future point
//g++ -std=c++11 -I./../include/ -I./../include/json -I./../include/csm/ -L./../lib/ -Wl,-E,-rpath="../lib" -o isdObject isdObject.cpp -lcsmapi
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>

#include <json/json.hpp>
#include <csm/Isd.h>
//#include <json.hpp>
//#include<Isd.h>

using namespace std;
using json = nlohmann::json;
//These are the different data types supported by the json library (with the exception of
//unknown, for handling extraneious input)
enum DataType{INT,UINT,FLOAT,STRING,BOOL,NULL8,UNKNOWN};

DataType checkType(json::value_type obj);

void addParam(csm::Isd * isd, json::iterator, DataType dt, int prec=12);
void printISD(csm::Isd isd);





int main(int argc, char *argv[]) {

  json jsonFile;
  csm::Isd isd;
  int prec = 12;

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <isd file>" << endl;
    return 1;
  }
  //Read the ISD file
  string line;
  string filename(argv[1]);
  ifstream file(argv[1]);
  if (!file.is_open()) {
    perror(("error while opening file " + filename).c_str());
  }

  else if (file.bad()) {
    perror(("error while reading file " + filename).c_str());
  }
  
  else {
    file >> jsonFile;
    isd.setFilename(filename);
    for (json::iterator i = jsonFile.begin(); i != jsonFile.end(); i++) {        
        if (i.value().is_array()){
          DataType arrayType = checkType(i.value()[0]);
          addParam(&isd, i,arrayType,prec);
        }
        else {
          DataType dt = checkType(i.value());
          addParam(&isd,i,dt,prec);
        }
    }//end for
  } //end outer-else
  printISD(isd);
  file.close();

  //return 0;
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
 * @param isd A pointer to the ISD object
 * @param it  The iterator to the json file which iterates over the keywords.
 * @param dt  The enum DataType value
 * @param prec The # of decimal places to be written to the ISD (if the value is a float)
 * @author Tyler Wilson
 */
void addParam(csm::Isd * isd, json::iterator it, DataType dt, int prec) {
  ostringstream key;
  //output the key to the ISD
  key << it.key();
  if (it.value().is_array()) {
    if (dt==FLOAT) {
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val << setprecision(prec) << v[j];
        isd->addParam(key.str(),val.str());
      }
    }
    else if(dt==INT){
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val << v[j];
        isd->addParam(key.str(),val.str());
      }
    }

    else if(dt==UINT) {
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val  << v[j];
        isd->addParam(key.str(),val.str());
      }
    }
    else if(dt ==BOOL) {
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val << v[j];
        isd->addParam(key.str(),val.str());
      }
    }
    else if (dt ==STRING) {
      vector<double> v = it.value();
      for (int j=0;j < v.size(); j++) {
        ostringstream val;
        val <<  v[j];
        isd->addParam(key.str(),val.str());
      }
    }
  }
  else {
    if(dt==FLOAT) {
      double v = it.value();
      ostringstream val;
      val << setprecision(prec) << v;
      isd->addParam(key.str(),val.str());
    }
    else if(dt==INT){
      int v = it.value();
      ostringstream val;
      val  << v;
      isd->addParam(key.str(),val.str());
    }
    else if(dt==UINT) {
      unsigned int v = it.value();
      ostringstream val;
      val << v;
      isd->addParam(key.str(),val.str());
    }
    else if(dt ==BOOL) {
      bool v = it.value();
      ostringstream val;
      val  << v;
      isd->addParam(key.str(),val.str());
    }
    else if (dt ==STRING) {
      string v = it.value();
      ostringstream val;
      val  << v;
      isd->addParam(key.str(),val.str());
    }
    else if (dt ==NULL8) {
      ostringstream val;
      val  << "null";
      isd->addParam(key.str(),val.str());
    }
  }//end outer else
}//end addParam


/**
 * @brief printISD  Display the keyword:value pairs of a CSM::ISD object
 * @param isd
 * @author Tyler Wilson
 */
void printISD(csm::Isd isd){
  const multimap<string,string> isdMap = isd.parameters();
  for (auto &i: isdMap)
    cout << i.first << " : " << i.second << endl;
}


