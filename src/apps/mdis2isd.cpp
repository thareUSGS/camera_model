
#include "Isis.h"

#include <cfloat>
#include <cstdio>
#include <iomanip>
#include <QPair>
#include <QList>
#include <QString>
#include <QStringList>

#include "Camera.h"
#include "CameraFocalPlaneMap.h"
#include "Distance.h"
#include "ProcessImportPds.h"
#include "ProcessByLine.h"

#include "UserInterface.h"
#include "FileName.h"
#include "IException.h"
#include "iTime.h"
#include "CameraPointInfo.h"
#include "TextFile.h"
#include "CSVReader.h"
#include "NaifStatus.h"
#include "PvlGroup.h"
#include "PvlKeyword.h"
#include "PvlObject.h"
#include "PvlSequence.h"
#include "LineManager.h"
#include "OriginalLabel.h"
#include "SpecialPixel.h"
#include "Spice.h"



using namespace std;
using namespace Isis;

void writeISD(const UserInterface &ui, PvlGroup * pvl);



void IsisMain() {

  UserInterface &ui = Application::GetUserInterface();
  CameraPointInfo campt;

  campt.SetCube(ui.GetFileName("FROM") + "+" + ui.GetInputAttribute("FROM").toString());
  PvlGroup * caminfo = campt.SetCenter(false,true);
  ProcessImportPds p;
  //FileName inFile = ui.GetFileName("FROM");
  //Cube icube(inFile);



  writeISD(ui,caminfo);


}

void writeISD(const UserInterface &ui,PvlGroup * caminfo){

  //QMap<QString,double> isdMap;
  QList<QPair<QString,double> > isdList;
  int prec =16;
  double spacecraftPosition[3] = {0.0,0.0,0.0};
  double instrumentPosition[3] = {0.0,0.0,0.0};
  double omegaPhiKappa[3] = {0.0,0.0,0.0};
  Distance dRadii[3];
  double isisFocalPlane2SocetPlate[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};


  FileName inFile = ui.GetFileName("FROM");
  Cube icube(inFile);

  if (icube.isProjected()) {
    QString msg = QString("You can only create a CSM from a level 1 image. "
                          "The input image [%1] is a map projected, level "
                          "2, cube.").arg(inFile.expanded());

    throw IException(IException::User, msg, _FILEINFO_);
  }

  // Make sure the image contains the SPICE blobs/tables
  PvlGroup test = icube.label()->findGroup("Kernels", Pvl::Traverse);
  QString instrumentPointing = (QString) test["InstrumentPointing"];
  if (instrumentPointing != "Table") {
    QString msg = QString("Input image [%1] does not contain needed SPICE blobs.  Please run "
                          "spiceinit on the image with attach=yes.").arg(inFile.expanded());

    throw IException(IException::User, msg, _FILEINFO_);
  }

  PvlObject naifKeywords = icube.label()->findObject("NaifKeywords");
  double boresightLine =0.0;
  double boresightSample=0.0;

  Camera *cam = icube.camera();
  CameraFocalPlaneMap *focalMap = cam->FocalPlaneMap();

  boresightLine = focalMap->DetectorSampleOrigin();
  boresightSample = focalMap->DetectorLineOrigin();
  double et = cam->time().Et();

  Spice spice(icube);
  spice.setTime(et);

  //Retrieve instrument position and target body radii in meters
  spice.radii(dRadii);
  double radii[3] = {0.0, 0.0, 0.0};
  radii[0] = dRadii[0].meters();
  radii[1] = dRadii[1].meters();
  radii[2] = dRadii[2].meters();

  //Retrieve Spacecraft position
  spacecraftPosition[0]=((*caminfo)["SpacecraftPosition"][0]).toDouble();
  spacecraftPosition[1]=((*caminfo)["SpacecraftPosition"][0]).toDouble();
  spacecraftPosition[2]=((*caminfo)["SpacecraftPosition"][0]).toDouble();

  spice.instrumentPosition(instrumentPosition);
  for (int i = 0; i < 3; i++) {
    instrumentPosition[i] *= 1000.0;
    spacecraftPosition[i] *= 1000.0;
  }

  // Fetch Bodyfixed -> Camera matrix w cspice
  vector<double>  j2000ToBodyFixedMatrixVector = spice.bodyRotation()->Matrix();
  vector<double>  j2000ToCameraMatrixVector = spice.instrumentRotation()->Matrix();


  // Reformat vector-matrices to 3x3 rotation matricies
  double j2000ToBodyFixedRotationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                                 {0.0, 0.0, 0.0},
                                                 {0.0, 0.0, 0.0}};

  double j2000ToCameraRotationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                              {0.0, 0.0, 0.0},
                                              {0.0, 0.0, 0.0}};

  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      j2000ToBodyFixedRotationMatrix[j][k] = j2000ToBodyFixedMatrixVector[3 * j + k];
      j2000ToCameraRotationMatrix[j][k] = j2000ToCameraMatrixVector[3 * j + k];
    }
  }

  // Compute Camera to Body Fixed rotation matrix
  double cameraToBodyFixedRotationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                                  {0.0, 0.0, 0.0},
                                                  {0.0, 0.0, 0.0}};


  mxmt_c(j2000ToBodyFixedRotationMatrix, j2000ToCameraRotationMatrix,
         cameraToBodyFixedRotationMatrix);



  PvlGroup inst = icube.label()->findGroup("Instrument", Pvl::Traverse);
  QString instrumentId = (QString) inst["InstrumentId"];
  QString spacecraftName = (QString) inst["SpacecraftName"];


  if (spacecraftName == "Messenger") {
    isisFocalPlane2SocetPlate[0][0] = 1.0;
    isisFocalPlane2SocetPlate[1][1] = -1.0;
    isisFocalPlane2SocetPlate[2][2] = -1.0;

  }


    //Calculate ographic coordinates of spacecraft position vector
    double xyzLength = instrumentPosition[0] * instrumentPosition[0] +
                  instrumentPosition[1] * instrumentPosition[1];
    double xyLength = sqrt(xyzLength);
    xyzLength = sqrt (xyzLength + instrumentPosition[2] * instrumentPosition[2]);
    double flattening = (radii[0] - radii[2]) / radii[0];
    double lon = 0.0;
    double lat = 0.0;
    double height = 0.0;
    recgeo_c (instrumentPosition, radii[0], flattening, &lon, &lat, &height);

    // Calculate rotation matrix from Socet Set plate to ocentric ground coordinates
    double socetPlateToOcentricGroundRotationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                                             {0.0, 0.0, 0.0},
                                                             {0.0, 0.0, 0.0}};

    mxmt_c (isisFocalPlane2SocetPlate, cameraToBodyFixedRotationMatrix,
            socetPlateToOcentricGroundRotationMatrix);



    // Populate the ocentric to ographic rotation matrix
    double ocentricToOgraphicRotationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                                     {0.0, 0.0, 0.0},
                                                     {0.0, 0.0, 0.0}};

    double sinLon = instrumentPosition[1] / xyLength;
    double cosLon = instrumentPosition[0] / xyLength;
    double sinLat = instrumentPosition[2] / xyzLength;
    double cosLat = xyLength / xyzLength;
    ocentricToOgraphicRotationMatrix[0][0] = -sinLon;
    ocentricToOgraphicRotationMatrix[1][0] = cosLon;
    ocentricToOgraphicRotationMatrix[2][0] = 0.0;
    ocentricToOgraphicRotationMatrix[0][1] = -sinLat * cosLon;
    ocentricToOgraphicRotationMatrix[1][1] = -sinLat * sinLon;
    ocentricToOgraphicRotationMatrix[2][1] = cosLat;
    ocentricToOgraphicRotationMatrix[0][2] = cosLat * cosLon;
    ocentricToOgraphicRotationMatrix[1][2] = cosLat * sinLon;
    ocentricToOgraphicRotationMatrix[2][2] = sinLat;

    // Compute the Rotation matrix from Socet Set plate to ographic ground coordinates
    // and extract the euler angles to get omega-phi-kappa attidude angles
    double socetPlateToOgrphicGroundRotationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0}};

    mxm_c (socetPlateToOcentricGroundRotationMatrix, ocentricToOgraphicRotationMatrix,
           socetPlateToOgrphicGroundRotationMatrix);


    double omega = 0.0;
    double phi = 0.0;
    double kappa = 0.0;


    m2eul_c (socetPlateToOgrphicGroundRotationMatrix, 3, 2, 1, &kappa, &phi, &omega);

    // Return resulting geographic lat, lon, omega, phi, kappa in decimal degrees
    // height in meters
    //ographicCamPos[0] = lat * RAD2DEG;
    //ographicCamPos[1] = lon * RAD2DEG;
    //ographicCamPos[2] = height;

    omegaPhiKappa[0] = omega * RAD2DEG;
    omegaPhiKappa[1] = phi * RAD2DEG;
    omegaPhiKappa[2] = kappa * RAD2DEG;

  ofstream os;

  QString isdFile = FileName(ui.GetFileName("TO")).expanded();

  os.open(isdFile.toLatin1().data(), ios::out);

  QString modelName("MDIS_SENSOR_MODEL");


  //isdList.append(qMakePair("ISD_LINE_PRINCIPAL_POINT_PIXELS",boresightLine);
  isdList.append(QPair<QString,double>("ISD_LINE_PRINCIPAL_POINT_PIXELS",boresightLine));
  isdList.append(QPair<QString,double>("ISD_SAMPLE_PRINCIPAL_POINT_PIXELS",boresightSample));
  isdList.append(QPair<QString,double>("ISD_FOCAL_LENGTH_PIXELS",double(naifKeywords["TempDependentFocalLength"])));
  isdList.append(QPair<QString,double>("ISD_NUMBER_OF_LINES",(double)(icube.lineCount())));
  isdList.append(QPair<QString,double>("ISD_NUMBER_OF_SAMPLES",(double)(icube.sampleCount())));
  isdList.append(QPair<QString,double>("ISD_SEMI_MAJOR_AXIS_METERS",radii[0]));
  isdList.append(QPair<QString,double>("ISD_SEMI_MINOR_AXIS_METERS",radii[1]));
  isdList.append(QPair<QString,double>("ISD_MIN_ELEVATION_METERS",0.0));
  isdList.append(QPair<QString,double>("ISD_MAX_ELEVATION_METERS",0.0));
  isdList.append(QPair<QString,double>("ISD_X_SENSOR_ORIG_METERS",instrumentPosition[0]));
  isdList.append(QPair<QString,double>("ISD_X_SENSOR_CURR_METERS",instrumentPosition[0]));
  isdList.append(QPair<QString,double>("ISD_Y_SENSOR_ORIG_METERS",instrumentPosition[1]));
  isdList.append(QPair<QString,double>("ISD_Y_SENSOR_CURR_METERS",instrumentPosition[1]));
  isdList.append(QPair<QString,double>("ISD_Z_SENSOR_ORIG_METERS",instrumentPosition[2]));
  isdList.append(QPair<QString,double>("ISD_Z_SENSOR_CURR_METERS",instrumentPosition[2]));
  isdList.append(QPair<QString,double>("ISD_OMEGA_ORIG_RADIANS",omegaPhiKappa[0]));
  isdList.append(QPair<QString,double>("ISD_OMEGA_CURR_RADIANS",omegaPhiKappa[0]));
  isdList.append(QPair<QString,double>("ISD_PHI_ORIG_RADIANS",omegaPhiKappa[1]));
  isdList.append(QPair<QString,double>("ISD_PHI_CURR_RADIANS",omegaPhiKappa[1]));
  isdList.append(QPair<QString,double>("ISD_KAPPA_ORIG_RADIANS",omegaPhiKappa[2]));
  isdList.append(QPair<QString,double>("ISD_KAPPA_CURR_RADIANS",omegaPhiKappa[2]));
  isdList.append(QPair<QString,double>("ISD_ORIGINAL_PARAMETER_COVARIANCE",0.0));
  isdList.append(QPair<QString,double>("ISD_CURRENT_PARAMETER_COVARIANCE",0.0));


  os << "ISD_SENSOR_MODEL_NAME\t";
  os << modelName << endl;
  os << setprecision(prec);

  for (int i =0;i < isdList.size();i++) {

    QPair<QString,double> isdNode=isdList[i];
    os << isdNode.first << "\t\t" << isdNode.second << endl;
  }

  os.close();




}


