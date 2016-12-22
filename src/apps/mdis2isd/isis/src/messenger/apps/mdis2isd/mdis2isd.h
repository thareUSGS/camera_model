#ifndef MDIS2ISD_H
#define MDIS2ISD_H

#include "Cube.h"
#include "CameraPointInfo.h"
#include "FileName.h"

class QString;

namespace Isis{

class mdis2isd
{
  public:

    mdis2isd(QString cubeFile);
    ~mdis2isd();
    void writeISD();


  private:

    QString m_cubeFileString;
    CameraPointInfo m_campt;
    bool m_validCube;
    static const int prec =16;


};

}
#endif // MDIS2ISD_H
