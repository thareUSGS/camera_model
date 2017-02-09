#include <cmath>
#include <dlfcn.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <csm.h>
#include <Isd.h>
#include <Plugin.h>

#include <gdal/gdal.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>
#include <gdal/cpl_string.h>

#include <IsdReader.h>
#include <MdisPlugin.h>
#include <MdisNacSensorModel.h>


using namespace std;

const double PI = 3.14159265358979323846;

void cubeArray(vector<vector<float> > *cube, GDALRasterBand *poBand);
double computeAngle(const std::vector<double> &v1, const std::vector<double> &v2);
double emissionAngle(const csm::EcefCoord &groundPt,
                     const csm::EcefCoord &sensorPosition);
double incidenceAngle(const csm::EcefCoord &groundPt,
                      const csm::EcefVector &illuminationDirection);
void writeCSV(const string &csvFile,
              const vector<string> &csvHeaders,
              const vector< vector<float> > &cubeData,
              const vector< vector<csm::EcefCoord> > &groundPoints);
void writeCSV(const string &csvFilename,
              const vector<string> &csvHeaders,
              const vector< vector<float> > &cubeData,
              const vector< vector< vector<float> > > &latLonPoints,
              const vector< vector<double> > &emissionAngles,
              const vector< vector<double> > &incidenceAngles);

int main(int argc, char *argv[]) {
  
  // Default cube and isd for now so this will run if user doesn't input these files.
  string isdFile("../../../tests/data/EN1007907102M.json");
  string cubeFile("../../../tests/data/EN1007907102M.cub");
  
  // User can provide ISD and cube if desired.
  if (argc == 3) {
    isdFile = argv[1];
    cubeFile = argv[2];
  }
  else {
    cout << "Usage: set <ISD.json> <cube.cub>\n";
    cout << "Provide an ISD .json file and its associated cube .cub file.\n";  
  }

  csm::Isd *isd = readISD(isdFile);
  if (isd == nullptr) {
    return 1;
  }

  // Find plugins TODO: probably need a dedicated area for plugins to load
  void *pluginFile = dlopen("../../../src/objs/libMdisPlugin.so", RTLD_LAZY);
  if (pluginFile == nullptr) {
    cout << "Could not load plugin." << endl;
    return 1;
  }

  // Choose the correct plugin (e.g. framing v. line-scan?)
  //const csm::PluginList &plugins = csm::Plugin::getList();
  const csm::Plugin *plugin = csm::Plugin::findPlugin("UsgsAstroFrameMdisPluginCSM");
  if (plugin == nullptr) {
    cout << "Could not find plugin: " << "UsgsAstroFrameMdisPluginCSM" << endl;
    return 1;
  }

  MdisNacSensorModel *model;
  
  // Initialize the MdisNacSensorModel from the ISD using the plugin
  try {
    model = dynamic_cast<MdisNacSensorModel*>
            (plugin->constructModelFromISD(*isd, "ISIS_MDISNAC_USGSAstro_1_Linux64_csm30.so"));
  }
  catch (csm::Error &e) {
    cout << e.what() << endl;
    delete isd;
    if (model != nullptr) {
      delete model;
    }
    exit (EXIT_FAILURE);
  }

  if (model == nullptr) {
     cout << "Could not construct the sensor model from the plugin." << endl;
     delete isd;
     return 0;
  }
  
  // Test the model's accuracy by visually comparing the results of the function calls
//   csm::EcefCoord groundPoint;
//   csm::ImageCoord imagePoint;
// 
//   for (int i = 0; i < 10; i++) {
//     groundPoint = model->imageToGround(imagePoint, 0);
// 
//     cout << "Image Point (s, l) : (" << imagePoint.samp << ", " << imagePoint.line << "); "
//               << "Ground Point (x, y, z) : (" << groundPoint.x << ", " << groundPoint.y << ", " 
//                                               << groundPoint.z << ")" << endl;
// 
//     //imagePoint = model->groundToImage(groundPoint);
//   }

  //Read from a cube using the GDAL API, and outputs the DN values to a 2D vector matrix
  //string cubePath("../../../tests/data/CN0108840044M_IF_5_NAC_spiced.cub");
  string cubePath(cubeFile);
  GDALDataset *poDataset;
  GDALRasterBand *poBand;
  int nBlockXSize, nBlockYSize;

  cout << endl;
  cout << "Testing GDAL"  << endl;
  cout << endl;
  GDALAllRegister();

  poDataset = (GDALDataset *)GDALOpen(cubePath.c_str(),GA_ReadOnly);

  if (poDataset == NULL) {
    cout << "Could not open the:  " + cubePath  << endl;
    delete isd;
    delete model;
    return 1;
  }

  else {
    poBand = poDataset->GetRasterBand(1);
    poBand->GetBlockSize(&nBlockXSize, &nBlockYSize);
    
    cout << "Num samples = " << nBlockXSize << endl;
    cout << "Num lines = " <<   nBlockYSize << endl;
    
    //Read a band of data
    vector< vector<float> > cubeMatrix;
    cubeArray(&cubeMatrix,poBand);
    
    // Get ground X,Y,Z for each pixel in image
    vector< vector<csm::EcefCoord> > groundPoints;
    // Get the emission and incidence angles for these XYZ's
    vector< vector<double> > emissionAngles, incidenceAngles;
    
    for (int line = 0; line < cubeMatrix.size(); line++) {
      vector<csm::EcefCoord> groundLine;
      vector<double> emaLine, inaLine;
      
      for (int sample = 0; sample < cubeMatrix[line].size(); sample++) {
        csm::ImageCoord imagePoint(line + 1, sample + 1);
        csm::EcefCoord groundPoint(model->imageToGround(imagePoint, 0.0));
        groundLine.push_back(groundPoint);
        
        // Calculate the emission angle using CSM's getSensorPosition
        csm::EcefCoord sensorPosition(model->getSensorPosition(imagePoint));
        emaLine.push_back(emissionAngle(groundPoint, sensorPosition));
        
        // Calculate the incidence angle using CSM's illumination vector.
        csm::EcefVector illumination(model->getIlluminationDirection(groundPoint));
        inaLine.push_back(incidenceAngle(groundPoint, illumination));
      }
      groundPoints.push_back(groundLine);
      emissionAngles.push_back(emaLine);
      incidenceAngles.push_back(inaLine);
    }
    
    // calculate the radius
    float r = sqrt(groundPoints[0][0].x * groundPoints[0][0].x + 
                    groundPoints[0][0].y * groundPoints[0][0].y + 
                    groundPoints[0][0].z * groundPoints[0][0].z);
    
    // calculate the lat & lon for each ground coordinate
    vector< vector < vector<float> > > latLonPoints;
    for (int i = 0; i < cubeMatrix.size(); i++) {
      vector< vector<float> > latLonLine;
      
      for (int j = 0; j < cubeMatrix[0].size(); j++) {  
        // calculate the latitude & longitude
        vector<float> latLon;
        
        float theta = acos(groundPoints[i][j].z / r); // lat in radians
        float phi = atan2(groundPoints[i][j].y, groundPoints[i][j].x); // lon in radians
        
        latLon.push_back(theta * 180 / M_PI);
        latLon.push_back(phi * 180 / M_PI);
        
        latLonLine.push_back(latLon);
      }
      latLonPoints.push_back(latLonLine);
    }
        
    // Write to csv file
    string csvFilename("ground.csv");
    vector<string> csvHeaders { 
          "Line", "Sample", "DN", "Latitude", "Longitude", "Emission", "Incidence"
    };
    writeCSV(csvFilename, csvHeaders, cubeMatrix, latLonPoints, emissionAngles, incidenceAngles);

  }  //end else
 
  delete isd;
  delete model;
  dlclose(pluginFile);
  exit (EXIT_SUCCESS);
}


void cubeArray(vector <vector<float> > *cube,GDALRasterBand *poBand) {

  vector<float> tempVector;
  float *pafScanline;
  int nsamps = poBand->GetXSize();
  int nlines = poBand->GetYSize();
  for (int j = 0;j<nlines;j++) {

    pafScanline = (float *)CPLMalloc(sizeof(float)*nsamps);
    poBand->RasterIO(GF_Read,0,j,nsamps,1,pafScanline,nsamps,1,GDT_Float32,0,0);

    for (int i = 0;i < nsamps;i++) {
      tempVector.push_back(pafScanline[i]);
    }
    cube->push_back(tempVector);
    tempVector.clear();
    CPLFree(pafScanline);
  }
}


/**
 * Computes an angle between two vectors by using the definition of the dot product.
 * 
 * @param v1 First vector.
 * @param v2 Second vector.
 * 
 * @return @b double Returns the angle between two vectors (in degrees).
 */
double computeAngle(const std::vector<double> &v1, const std::vector<double> &v2) {
  // Recall: dot(v1, v2) = ||v1|| * ||v2|| * cos(theta)
  double dotProduct = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  double magnitude = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]) *
                     sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
  double theta = acos(dotProduct / magnitude);
  
  // Convert to degrees.
  return theta * (180 / PI);
}


/**
 * Calculates the emission angle (degrees) between a ground point's normal and the spacecraft's
 * position (body-fixed).
 * 
 * @param groundPt Body-fixed ground point to use in calculation.
 * @param spacecraftPosition Body-fixed spacecraft position.
 * 
 * @return @b double Returns the emission angle (in degrees) between ground point normal and the
 *                   spacecraft position.
 */
double emissionAngle(const csm::EcefCoord &groundPt,
                     const csm::EcefCoord &sensorPosition) {
  std::vector<double> ground { groundPt.x, groundPt.y, groundPt.z };
  std::vector<double> sensor { sensorPosition.x, sensorPosition.y, sensorPosition.z };
  return computeAngle(ground, sensor);
}


/**
 * Calculates the incidence angle (degrees) between a ground point's normal and the illumination
 * vector from the sun.
 * 
 * @param groundPt Body-fixed ground point to use in calculation.
 * @param illuminationDirection Body-fixed illumination vector of the sun.
 * 
 * @return @b double Returns the incidence angle (in degrees) between ground point normal and the
 *                   illumination vector.
 */
double incidenceAngle(const csm::EcefCoord &groundPt,
                      const csm::EcefVector &illuminationDirection) {
  // First find the body-fixed position of the sun.
  std::vector<double> sun {
    groundPt.x - illuminationDirection.x,
    groundPt.y - illuminationDirection.y,
    groundPt.z - illuminationDirection.z
  };
  
  // Solve for the angle between the sun position and the ground position.
  std::vector<double> ground { groundPt.x, groundPt.y, groundPt.z };
  return computeAngle(ground, sun);
}


/**
 * Writes a CSV file to the given destination.
 * 
 * This CSV will contain a DN, X, Y, and Z for each Line,Sample of the input cube data.
 * 
 * @param csvFilename Name of the output CSV file to write to.
 * @param csvHeaders Vector containing the header elements to write to the CSV file.
 * @param cubeData Matrix of cube DNs.
 * @param groundPoints Matrix of ground points for the image pixels.
 */
void writeCSV(const string &csvFilename,
              const vector<string> &csvHeaders,
              const vector< vector<float> > &cubeData,
              const vector< vector<csm::EcefCoord> > &groundPoints) {
  
  ofstream csvFile(csvFilename);
  if (csvFile.is_open()) {
    
    // Write the csv header
    for (int str = 0; str < csvHeaders.size() - 1; str++) {
      csvFile << csvHeaders[str] << ", ";
    }
    
    // Write the last header element
    csvFile << csvHeaders[csvHeaders.size() - 1] << "\n";
    
    // Write the csv records
    for (int line = 0; line < cubeData.size(); line++) {
      
      for (int sample = 0; sample < cubeData[line].size(); sample++) {
        csm::EcefCoord ground = groundPoints[line][sample];
        csvFile << line + 1 << ", " << sample + 1 << ", " << cubeData[line][sample] << ", "
                << ground.x/1000 << ", " << ground.y/1000 << ", " << ground.z/1000 << "\n";
      }
    }
  }
  
  else {
    cout << "\nUnable to open file \"" << csvFilename << " for writing." << endl;
    return;
  }
}


/**
 * Writes a CSV file to the given destination.
 * 
 * This CSV will contain a DN, Latitude, and Longitude for each Line,Sample of the input cube data.
 * 
 * @param csvFilename Name of the output CSV file to write to.
 * @param csvHeaders Vector containing the header elements to write to the CSV file.
 * @param cubeData Matrix of cube DNs.
 * @param latLonPoints Matrix of latitude/longitude points for the image pixels.
 */
void writeCSV(const string &csvFilename,
              const vector<string> &csvHeaders,
              const vector< vector<float> > &cubeData,
              const vector< vector< vector<float> > > &latLonPoints,
              const vector< vector<double> > &emissionAngles,
              const vector< vector<double> > &incidenceAngles) {
  
  ofstream csvFile(csvFilename);
  if (csvFile.is_open()) {
    
    // Write the csv header
    for (int str = 0; str < csvHeaders.size() - 1; str++) {
      csvFile << csvHeaders[str] << ", ";
    }
    
    // Write the last header element
    csvFile << csvHeaders[csvHeaders.size() - 1] << "\n";
    
    // Write the csv records
    for (int line = 0; line < cubeData.size(); line++) {
      
      for (int sample = 0; sample < cubeData[line].size(); sample++) {
        vector<float> latLon = latLonPoints[line][sample];
        csvFile << line + 1 << ", "
                << sample + 1 << ", "
                << cubeData[line][sample] << ", "
                << latLon[0] << ", " << latLon[1] << ", "
                << emissionAngles[line][sample] << ", "
                << incidenceAngles[line][sample] << "\n";
      }
    }
  }
  
  else {
    cout << "\nUnable to open file \"" << csvFilename << " for writing." << endl;
    return;
  }
}
