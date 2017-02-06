#include <string>

#include <csm/Isd.h>

#include <gtest/gtest.h>

#include <MdisPlugin.h>
#include <MdisNacSensorModel.h>
#include <IsdReader.h>

#include "MdisNacSensorModelTest.h"

bool MdisNacSensorModelTest::setupFixtureFailed = false;
std::string MdisNacSensorModelTest::setupFixtureError;
csm::Isd *MdisNacSensorModelTest::isd = nullptr;
std::string MdisNacSensorModelTest::dataFile;
MdisPlugin MdisNacSensorModelTest::mdisPlugin;
MdisNacSensorModel *MdisNacSensorModelTest::mdisModel = nullptr;

/*
 * Test imageToGround - truth extracted as follows:
 * setisis isis3
 * qview /work/projects/IAA_camera/data/EN100790102M.cub
 * F (selects "Find Tool")
 * On top toolbar, select "Find Point"
 * Type in 512.5, 512.5 for Sample/Line (ISIS3 pixel center = 1,1)
 * Click "Record Point"
 * Check "XYZ" -> { 1132.18, -1597.75, 1455.66 }
 */
TEST_F(MdisNacSensorModelTest, imageToGround1) {

  // gtest #247 work-around
  if (setupFixtureFailed) {
    FAIL() << setupFixtureError;
  }

  // CSM Line/Sample center = 512, 512
  csm::ImageCoord point(512.0, 512.0);
  double height = 0.0;
  csm::EcefCoord xyz = mdisModel->imageToGround(point, height);
  double truth[] = { 1132.18*1000, -1597.75*1000, 1455.66*1000 };
  EXPECT_EQ(truth[0], xyz.x);
  EXPECT_EQ(truth[1], xyz.y);
  EXPECT_EQ(truth[2], xyz.z);
}


// Test groundToImage
TEST_F(MdisNacSensorModelTest, groundToImage1) {
  // gtest #247 work-around
  if (setupFixtureFailed) {
    FAIL() << setupFixtureError;
  }

  double x = 1132.18*1000;
  double y = -1597.75*1000;
  double z = 1455.66*1000;
  // Override xyz - distortion has some issues, so we will use undistorted xyz for
  // input line,sample 100,100. See python vector ground to image notebook for details.
  x = 1115920.0;
  y = -1603550.0;
  z = 1460830.0;
  csm::EcefCoord xyz(x, y, z);
  csm::ImageCoord pt = mdisModel->groundToImage(xyz);
  // Use 1/2 pixel as tolerance
  EXPECT_NEAR(100.0, pt.line, 0.5);
  EXPECT_NEAR(100.0, pt.samp, 0.5);
}

TEST_F(MdisNacSensorModelTest, distortionModel2) {
  double dx = -6.30;
  double dy = 6.40;
  double udx = 0.0;
  double udy = 0.0;
  double isis3_udx = -6.3036234000160273893698104075156152248383;
  double isis3_udy = 6.3445144408882310216313271666876971721649;
  testMath.setFocalPlane(dx,dy,udx,udy);

  EXPECT_NEAR(udx,isis3_udx,tolerance);
  EXPECT_NEAR(udy,isis3_udy,tolerance);

}
