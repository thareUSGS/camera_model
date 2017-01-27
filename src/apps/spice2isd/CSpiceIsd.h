#ifndef CSPICEISD_H
#define CSPICEISD_H

#include <utility>
#include <vector>
#include <string>

using namespace std;



class CSpiceIsd
{
  public:
    CSpiceIsd(string cubeFile);
   ~CSpiceIsd();
    void isdJSON(vector<pair<string,double> > * isdData,string sensorModel,
                           string filePath);
    void writeISD();

  private:

    string m_cubeFileString;
    bool m_validCube;
    static const int prec =16;

};


#endif
