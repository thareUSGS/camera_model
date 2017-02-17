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
 * Type in 513, 513 for Sample/Line (ISIS3 pixel center = 1,1)
 * Click "Record Point"
 * Check "XYZ" -> { 1132.18, -1597.75, 1455.66 }
 */
TEST_F(MdisNacSensorModelTest, imageToGroundCenter) {

  // gtest #247 work-around
  if (setupFixtureFailed) {
    FAIL() << setupFixtureError;
  }

  csm::ImageCoord point(512.5, 512.5);
  double height = 0.0;
  csm::EcefCoord xyz = mdisModel->imageToGround(point, height);
  double truth[] = { 1129.25*1000, -1599.26*1000, 1455.28*1000 };
  EXPECT_EQ(truth[0], xyz.x);
  EXPECT_EQ(truth[1], xyz.y);
  EXPECT_EQ(truth[2], xyz.z);
}

TEST_F(MdisNacSensorModelTest, imageToGroundOffCenter){
  if (setupFixtureFailed) {
    FAIL() << setupFixtureError;
  }

  csm::ImageCoord point(100, 100);
  double height = 0.0;
  csm::EcefCoord xyz = mdisModel->imageToGround(point, height);
  double truth[] = { 1115.95*1000, -1603.44*1000, 1460.93*1000 };
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

  double x = 1129.25 * 1000;
  double y = -1599.26 * 1000;
  double z = 1455.28 * 1000;
  csm::EcefCoord xyz(x, y, z);
  csm::ImageCoord pt = mdisModel->groundToImage(xyz);
  // Use 1/2 pixel as tolerance
  EXPECT_NEAR(512.5, pt.line, 0.1);
  EXPECT_NEAR(512.5, pt.samp, 0.1);
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


TEST_F(MdisNacSensorModelTest, getImageStart) {
  csm::ImageCoord start = mdisModel->getImageStart();
  EXPECT_EQ(start.line, 1.0);
  EXPECT_EQ(start.samp, 9.0);
}


TEST_F(MdisNacSensorModelTest, getImageSize) {
  csm::ImageVector size = mdisModel->getImageSize();
  EXPECT_EQ(size.line, 1024);
  EXPECT_EQ(size.samp, 1024);
}


TEST_F(MdisNacSensorModelTest, getImageTime) {
  csm::ImageCoord point;
  point.samp = 500;
  point.line = 500;
  double time = mdisModel->getImageTime(point);
  EXPECT_NEAR(time, 418855170.49299997, tolerance);
}


TEST_F(MdisNacSensorModelTest, getSensorPosition) {
  csm::ImageCoord point;
  point.samp = 500;
  point.line = 500;
  csm::EcefCoord position = mdisModel->getSensorPosition(point);
  EXPECT_NEAR(position.x, 1728357.7031238307, tolerance);
  EXPECT_NEAR(position.y, -2088409.0061042644, tolerance);
  EXPECT_NEAR(position.z, 2082873.9280557402, tolerance);
}

