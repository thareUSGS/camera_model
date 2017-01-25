

#include "SpiceController.h"
#include "cspiceisd.h"
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


    cspiceisd cspice("blah.cub");
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


#if 0
    sc.loadKernel(kernel1);
    sc.loadKernel(kernel2);
    sc.loadKernel(kernel3);
    sc.loadKernel(kernel4);
    sc.loadKernel(kernel5);
    sc.loadKernel(kernel6);
    sc.loadKernel(kernel7);
    sc.loadKernel(kernel8);
    sc.loadKernel(kernel9);
    sc.loadKernel(kernel10);
#endif


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
    int bGotMin,bGotMax;
    double adfMinMax[2];


    GDALAllRegister();

    poDataset = (GDALDataset *)GDALOpen(cubePath.c_str(),GA_ReadOnly);

    if(poDataset == NULL) {

      cout << "Could not open the:" + cubePath  << endl;

    }

    else {

      poBand = poDataset->GetRasterBand(1);
      poBand->GetBlockSize(&nBlockXSize,&nBlockYSize);

      printf( "Block=%dx%d Type=%s, ColorInterp=%s\n",nBlockXSize, nBlockYSize,
              GDALGetDataTypeName(poBand->GetRasterDataType()),
              GDALGetColorInterpretationName(poBand->GetColorInterpretation()) );



    //Read a band of data


      vector<vector<float> > cubeMatrix;
      cubeArray(&cubeMatrix,poBand);

      for (int i =0; i < cubeMatrix.size(); i++ ) {

        vector<float> v = cubeMatrix[i];
        for (int j = 0; j < v.size(); j++)
          cout << v[j] << endl;


      }


#if 0




    float *pafScanline;
    int nXSize = poBand->GetXSize();
    int nYSize = poBand->GetYSize();






    for (int j = 0;j<nYSize;j++) {

    pafScanline = (float *)CPLMalloc(sizeof(float)*nXSize);
    poBand->RasterIO(GF_Read,0,j,nXSize,1,pafScanline,nXSize,1,GDT_Float32,0,0);

    for (int i = 0;i < nXSize;i++)
      cout << setprecision(10) << pafScanline[i] << endl;


    CPLFree(pafScanline);
    }

#endif
    }  //end else

}

void cubeArray(vector <vector<float> > *cube,GDALRasterBand *poBand) {


  vector<float> tempVector;

  float *pafScanline;
  int nsamps = poBand->GetXSize();
  int nlines = poBand->GetYSize();


  for (int j = 0;j<nlines;j++) {

    pafScanline = (float *)CPLMalloc(sizeof(float)*nsamps);
    poBand->RasterIO(GF_Read,0,j,nsamps,1,pafScanline,nsamps,1,GDT_Float32,0,0);

    for (int i = 0;i < nsamps;i++)
      tempVector.push_back(pafScanline[i]);

    cube->push_back(tempVector);
    tempVector.clear();
    CPLFree(pafScanline);
  }




}



