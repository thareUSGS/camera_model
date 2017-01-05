#include "mdis2isd.h"

#include <fstream>
#include <iomanip>

#include <QString>


#include "Camera.h"
#include "CameraFocalPlaneMap.h"
#include "CameraPointInfo.h"
#include "Distance.h"


#include "FileName.h"
#include "IException.h"
#include "iTime.h"


#include "NaifStatus.h"
#include "PvlGroup.h"
#include "PvlKeyword.h"
#include "PvlObject.h"
#include "Spice.h"



using namespace std;


namespace Isis {
  mdis2isd::mdis2isd(QString cubeFileName) {
    m_cubeFileString = cubeFileName;
    m_validCube = true;
    try {

    m_campt.SetCube(m_cubeFileString);
    }
    catch (IException &e) {
           //QString msg = "Unable to call CameraPointInfo::SetCube";

           cout << e.toString() << endl;
           m_validCube = false;

           //throw IException(e, IException::Unknown, msg, _FILEINFO_);
         }


    }


  mdis2isd::~mdis2isd(){



  }







  /**
   * @brief mdis2isd::isdJSON  This function outputs the ISD file created by this
   * application in JSON format.  It is very simple, and will need to be modified.
   * The key-value appairs assume <string,double>, which might not be the case as
   * the second value could be anything from a primitive type to an array of primitive
   * types.
   * @param isdList A ptr to a vector of key-value pairs storing the ISD values
   * @param sensorModel The name of the sensor model.
   * @param The output path and name for the JSON file.
   */

  void mdis2isd::isdJSON(std::vector<std::pair<std::string,double> > * isdList,std::string sensorModel,
                         std::string filePath){


    ofstream os;
    os.open(filePath.c_str(), ofstream::out);

    os << "{" << endl;
    os << "\"" <<"ISD_SENSOR_MODEL_NAME"<<"\":" << sensorModel << ",";

    os << setprecision(prec);
    unsigned int nparams = isdList->size();
    for (unsigned int i =0;i < nparams-1;i++) {

      pair<string,double> isdNode=isdList->at(i);
      os << "\"" << isdNode.first << "\":" << isdNode.second << ",";
    }

     pair<string,double> lastNode=isdList->at(nparams-1);
     os << "\"" << lastNode.first << "\":" << lastNode.second << "}";


    os.close();


  }

  /**
   * @brief mdis2isd:writeISD This function grabs the necessary Spice date from ISIS
   * and outputs a simple ISD file.
   */

  void mdis2isd::writeISD(){


    if (m_validCube == false) {
      cout << "Invalid cube" << endl;
      return;

    }
    PvlGroup * caminfo = m_campt.SetCenter(false,true);

    std::vector<std::pair<string,double> > isdList;
    double spacecraftPosition[3] = {0.0,0.0,0.0};
    double instrumentPosition[3] = {0.0,0.0,0.0};
    double omegaPhiKappa[3] = {0.0,0.0,0.0};
    Distance dRadii[3];
    double isisFocalPlane2SocetPlate[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    FileName inFile(m_cubeFileString);
    Cube icube(inFile);

    if (icube.isProjected()) {
      QString msg = QString("You can only create a CSM from a level 1 image. "
                            "The input image [%1] is a map projected, level "
                            "2, cube.").arg(inFile.expanded());

      cout << msg << endl;
      return;
    }

    // Make sure the image contains the SPICE blobs/tables
    PvlGroup test = icube.label()->findGroup("Kernels", Pvl::Traverse);
    QString instrumentPointing = (QString) test["InstrumentPointing"];
    if (instrumentPointing != "Table") {

      QString msg = QString("Input image [%1] does not contain needed SPICE blobs.  Please run "
                            "spiceinit on the image with attach=yes.").arg(inFile.expanded());


      cout << msg << endl;
      return;

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

    QString isdFile = inFile.expanded().split(".",QString::SkipEmptyParts).at(0)+".isd";



    os.open(isdFile.toLatin1().data(), ios::out);

    QString modelName("MDIS_SENSOR_MODEL");


    isdList.push_back(pair<string,double>("ISD_LINE_PRINCIPAL_POINT_PIXELS",boresightLine));
    isdList.push_back(pair<string,double>("ISD_SAMPLE_PRINCIPAL_POINT_PIXELS",boresightSample));
    isdList.push_back(pair<string,double>("ISD_FOCAL_LENGTH_PIXELS",double(naifKeywords["TempDependentFocalLength"])));
    isdList.push_back(pair<string,double>("ISD_NUMBER_OF_LINES",(double)(icube.lineCount())));
    isdList.push_back(pair<string,double>("ISD_NUMBER_OF_SAMPLES",(double)(icube.sampleCount())));
    isdList.push_back(pair<string,double>("ISD_SEMI_MAJOR_AXIS_METERS",radii[0]));
    isdList.push_back(pair<string,double>("ISD_SEMI_MINOR_AXIS_METERS",radii[1]));
    isdList.push_back(pair<string,double>("ISD_MIN_ELEVATION_METERS",0.0));
    isdList.push_back(pair<string,double>("ISD_MAX_ELEVATION_METERS",0.0));
    isdList.push_back(pair<string,double>("ISD_X_SENSOR_ORIG_METERS",instrumentPosition[0]));
    isdList.push_back(pair<string,double>("ISD_X_SENSOR_CURR_METERS",instrumentPosition[0]));
    isdList.push_back(pair<string,double>("ISD_Y_SENSOR_ORIG_METERS",instrumentPosition[1]));
    isdList.push_back(pair<string,double>("ISD_Y_SENSOR_CURR_METERS",instrumentPosition[1]));
    isdList.push_back(pair<string,double>("ISD_Z_SENSOR_ORIG_METERS",instrumentPosition[2]));
    isdList.push_back(pair<string,double>("ISD_Z_SENSOR_CURR_METERS",instrumentPosition[2]));
    isdList.push_back(pair<string,double>("ISD_OMEGA_ORIG_RADIANS",omegaPhiKappa[0]));
    isdList.push_back(pair<string,double>("ISD_OMEGA_CURR_RADIANS",omegaPhiKappa[0]));
    isdList.push_back(pair<string,double>("ISD_PHI_ORIG_RADIANS",omegaPhiKappa[1]));
    isdList.push_back(pair<string,double>("ISD_PHI_CURR_RADIANS",omegaPhiKappa[1]));
    isdList.push_back(pair<string,double>("ISD_KAPPA_ORIG_RADIANS",omegaPhiKappa[2]));
    isdList.push_back(pair<string,double>("ISD_KAPPA_CURR_RADIANS",omegaPhiKappa[2]));
    isdList.push_back(pair<string,double>("ISD_ORIGINAL_PARAMETER_COVARIANCE",0.0));
    isdList.push_back(pair<string,double>("ISD_CURRENT_PARAMETER_COVARIANCE",0.0));


    os << "ISD_SENSOR_MODEL_NAME\t";
    os << modelName << endl;
    os << setprecision(prec);

    for (unsigned int i =0;i < isdList.size();i++) {

      pair<string,double> isdNode=isdList[i];
      os << isdNode.first << "\t\t" << isdNode.second << endl;
    }

    os.close();


    isdJSON(&isdList,modelName.toStdString(),"json.isd");






  }

}
