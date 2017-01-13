#ifndef IsdReader_h
#define IsdReader_h

#include <csm/Isd.h>
#include <json/json.hpp>

using namespace std;
using json = nlohmann::json;

//These are the different data types supported by the json library (with the exception of
//unknown, for handling extraneious input)
enum DataType {
  INT,
  UINT,
  FLOAT,
  STRING,
  BOOL,
  NULL8,
  UNKNOWN
};

void addParam(csm::Isd &isd, json::iterator, DataType dt, int prec=12);
DataType checkType(json::value_type obj);
csm::Isd *readISD(string filename);
void printISD(const csm::Isd &isd);

#endif