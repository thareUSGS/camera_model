

#include "SpiceController.h"
#include "CSpiceIsd.h"
#include <gdal/gdal.h>
#include <SpiceUsr.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>



#include <gdal/gdal.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>
#include <gdal/cpl_string.h>




using namespace std;




void cubeArray(vector<vector<float> > *cube, GDALRasterBand *poBand);

int main(int argc,char *argv[]) {

    int ikid = 236820;

    vector<pair<string,SpiceDouble> >  isdList;


    CSpiceIsd cspice("blah.cub");
    SpiceController sc;

    string kernel1("data/msgr_v231.tf");
    string kernel2("data/msgr_v231.tf");
    string kernel3("data/pck00010_msgr_v23.tpc");
    string kernel4("data/msgr_dyn_v600.tf");
    string kernel5("data/msgr_v231.tf");
    string kernel6("data/msgr_de405_de423s.bsp");
    string kernel7("data/msgr_mdis_sc050727_100302_sub_v1.bc");
    string kernel8("data/msgr_mdis_gm040819_150430v1.bc");
    string kernel9("data/naif0011.tls");
    string kernel10("data/messenger_2548.tsc");


    furnsh_c("data/msgr_mdis_v160.ti");
    furnsh_c("data/msgr_v231.tf");
    furnsh_c("data/pck00010_msgr_v23.tpc");
    furnsh_c("data/msgr_dyn_v600.tf");
    furnsh_c("data/msgr_v231.tf");
    furnsh_c("data/msgr_de405_de423s.bsp");
    furnsh_c("data/msgr_mdis_sc050727_100302_sub_v1.bc");
    furnsh_c("data/msgr_mdis_gm040819_150430v1.bc");
    furnsh_c("data/naif0011.tls");
    furnsh_c("data/messenger_2548.tsc");


    SpiceInt code;
    SpiceInt found;
    bodn2c_c("MESSENGER",&code,&found);


    SpiceDouble focalLength;
    SpiceInt numLines;
    SpiceInt numSamples;
    SpiceInt n;
    SpiceDouble odt_x[9];
    SpiceDouble odt_y[9];
    SpiceDouble flEpsilon;
    SpiceDouble pixelPitch;
    SpiceDouble ccdCenter;
    SpiceDouble ifov;
    SpiceDouble boresight[3];
    SpiceDouble transX[3];
    SpiceDouble transY[3];
    SpiceDouble itranss[3];
    SpiceDouble itransl[3];
    SpiceDouble startSample;
    SpiceDouble startLine;

    int prec = 10;

    gdpool_c ("INS-236820_FOCAL_LENGTH", 0,1, &n, &focalLength, &found );
    gdpool_c ("INS-236820_FL_UNCERTAINTY", 0,1, &n, &flEpsilon, &found );
    gipool_c ("INS-236820_PIXEL_LINES", 0,1, &n, &numLines, &found );
    gipool_c ("INS-236820_PIXEL_SAMPLES", 0,1, &n, &numSamples, &found );
    gdpool_c ("INS-236820_PIXEL_PITCH", 0,1, &n, &pixelPitch, &found );
    gdpool_c ("INS-236820_CCD_CENTER", 0,1, &n, &ccdCenter, &found );
    gdpool_c ("INS-236820_IFOV", 0,1, &n, &ifov, &found );
    gdpool_c ("INS-236820_BORESIGHT", 0,3, &n, boresight, &found );
    gdpool_c ("INS-236820_TRANSX", 0,3, &n, transX, &found );
    gdpool_c ("INS-236820_TRANSY", 0,3, &n, transY, &found );
    gdpool_c ("INS-236820_ITRANSS", 0,3, &n, itranss, &found );
    gdpool_c ("INS-236820_ITRANSL", 0,3, &n, itransl, &found );

    gdpool_c ("INS-236820_OD_T_X", 0,9, &n, odt_x, &found );
    gdpool_c ("INS-236820_OD_T_Y", 0,9, &n, odt_y, &found );

    gdpool_c ("INS-236820_FPUBIN_START_SAMPLE", 0,1, &n, &startSample, &found );
    gdpool_c ("INS-236820_FPUBIN_START_LINE", 0,1, &n, &startLine, &found );


    string cubePath("CN0108840044M_IF_5_NAC_spiced.cub");
    GDALDataset *poDataset;
    GDALRasterBand *poBand;
    int nBlockXSize,nBlockYSize;


    GDALAllRegister();

    poDataset = (GDALDataset *)GDALOpen(cubePath.c_str(),GA_ReadOnly);

    if(poDataset == NULL) {

      cout << "Could not open the:" + cubePath  << endl;

    }

    else {

      poBand = poDataset->GetRasterBand(1);
      poBand->GetBlockSize(&nBlockXSize,&nBlockYSize);

    //Read a band of data

      vector<vector<float> > cubeMatrix;
      cubeArray(&cubeMatrix,poBand);
      for (int i =0; i < cubeMatrix.size(); i++ ) {
        vector<float> v = cubeMatrix[i];

        for (int j = 0; j < v.size(); j++) {
            cout << v[j] << endl;
        }//end inner-for


      }//end outer-for

    }  //end else

}


/**
 * @brief cubeArray:  Translates a GDALRasterBand into a 2D vector matrix
 * @param cube:  The 2D vector matrix which is output by this function.
 * @param poBand:  The GDALRasterBand obtained from the ISIS3 cube.
 */

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
    //free the memory allocated to store the current scanline
    CPLFree(pafScanline);
  }




}



