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


// Test imageToProximateImagingLocus
TEST_F(MdisNacSensorModelTest, imageToProximateImagingLocus1) {
  // gtest #247 work-around
  if (setupFixtureFailed) {
    FAIL() << setupFixtureError;
  }
  
  csm::ImageCoord point(512.0, 512.0);
  csm::EcefCoord ground(0,0,0);
  csm::EcefLocus proximateLocus = mdisModel->imageToProximateImagingLocus(point, ground);
  
  double spacecraftX = atof(isd->param("x_sensor_origin").c_str());
  double spacecraftY = atof(isd->param("y_sensor_origin").c_str());
  double spacecraftZ = atof(isd->param("z_sensor_origin").c_str());
  EXPECT_EQ(spacecraftX, proximateLocus.point.x);
  EXPECT_EQ(spacecraftY, proximateLocus.point.y);
  EXPECT_EQ(spacecraftZ, proximateLocus.point.z);
  
  EXPECT_NEAR(-0.601527, proximateLocus.direction.x, tolerance);
  EXPECT_NEAR(0.491036, proximateLocus.direction.y, tolerance);
  EXPECT_NEAR(-0.630118, proximateLocus.direction.z, tolerance);
}


// Test imageToRemoteImagingLocus
TEST_F(MdisNacSensorModelTest, imageToRemoteImagingLocus1) {
  // gtest #247 work-around
  if (setupFixtureFailed) {
    FAIL() << setupFixtureError;
  }
  
  csm::ImageCoord point(512.0, 512.0);
  csm::EcefLocus remoteLocus = mdisModel->imageToRemoteImagingLocus(point);
  
  // Compare the locus origin point with the body-fixed spacecraft position
  double spacecraftX = atof(isd->param("x_sensor_origin").c_str());
  double spacecraftY = atof(isd->param("y_sensor_origin").c_str());
  double spacecraftZ = atof(isd->param("z_sensor_origin").c_str());
  EXPECT_EQ(spacecraftX, remoteLocus.point.x);
  EXPECT_EQ(spacecraftY, remoteLocus.point.y);
  EXPECT_EQ(spacecraftZ, remoteLocus.point.z);
  
  // Truth values pulled from imageToGround's normalized direction vector
  EXPECT_NEAR(-0.601527, remoteLocus.direction.x, tolerance);
  EXPECT_NEAR(0.491036, remoteLocus.direction.y, tolerance);
  EXPECT_NEAR(-0.630118, remoteLocus.direction.z, tolerance);
}


// Tests the getModelState() method with a default constructed MdisNacSensorModel.
TEST_F(MdisNacSensorModelTest, getModelStateDefault) {
  EXPECT_EQ(defaultMdisNac.getModelState(), std::string());
}


// Test getElevation
TEST_F(MdisNacSensorModelTest, computeElevationOnSphere) {
  // (1/4)^2 + (1/2)^2 + z^2 = 1^2; z^2 = 11/16
  double elevation = testMath.computeElevation(0.25, 0.5, sqrt(11)/4.0);
  EXPECT_EQ(0.0, elevation);
}


// Test intersect
TEST_F(MdisNacSensorModelTest, intersectTrivial) {
  std::vector<double> position { 0.0, 0.0, 1.5 };
  std::vector<double> look { 0.0, 0.0, -0.5 };
  csm::EcefCoord intersectGround = testMath.intersect(position, look, 1.0);
  EXPECT_EQ(0.0, intersectGround.x);
  EXPECT_EQ(0.0, intersectGround.y);
  EXPECT_EQ(1.0, intersectGround.z);
}

TEST_F(MdisNacSensorModelTest, intersectLookingAway) {
  std::vector<double> position { 0.0, 0.0, 2.0 };
  std::vector<double> look { 0.0, 0.0, 0.5 };
  csm::EcefCoord ground = testMath.intersect(position, look, 1.0);
  EXPECT_EQ(0.0, ground.x);
  EXPECT_EQ(0.0, ground.y);
  EXPECT_EQ(0.0, ground.z);
}


// Test perpendicular
TEST_F(MdisNacSensorModelTest, perpendicularNonZeros) {
  std::vector<double> v1(3), v2(3);
  v1[0] = 0.0;
  v1[1] = 0.0;
  v1[2] = 1.5;
  v2[0] = 0.0;
  v2[1] = -0.25;
  v2[2] = -0.5;
  std::vector<double> result = testMath.perpendicular(v1, v2);
  EXPECT_NEAR(0.0, result[0], tolerance);
  EXPECT_NEAR(-0.6, result[1], tolerance);
  EXPECT_NEAR(0.3, result[2], tolerance);
}


// Test project
TEST_F(MdisNacSensorModelTest, projectNonZeros) {
  std::vector<double> v1(3), v2(3);
  v1[0] = 0.0;
  v1[1] = 0.0;
  v1[2] = 1.5;
  v2[0] = 0.0;
  v2[1] = -0.25;
  v2[2] = -0.5;
  std::vector<double> result = testMath.project(v1, v2);
  EXPECT_NEAR(0.0, result[0], tolerance);
  EXPECT_NEAR(0.6, result[1], tolerance);
  EXPECT_NEAR(1.2, result[2], tolerance);
}


// Test dot
// TODO: use value-parameterized tests
TEST_F(MdisNacSensorModelTest, dotZeroVectors) {
  std::vector<double> v1(3, 0.0), v2(3, 0.0);
  EXPECT_EQ(0.0, testMath.dot(v1, v2));
}

TEST_F(MdisNacSensorModelTest, dotOneZeroVector) {
  std::vector<double> v1(3, 0.0), v2(3);
  v2[0] = 1.0;
  v2[1] = 2.0;
  v2[2] = 3.0;
  EXPECT_EQ(0.0, testMath.dot(v1, v2));
}

TEST_F(MdisNacSensorModelTest, dotVectors) {
  std::vector<double> v1(3), v2(3);
  v1[0] = 1.0;
  v1[1] = 2.0;
  v1[2] = 3.0;
  v2[0] = -2.0;
  v2[1] = 2.0;
  v2[2] = 3.0;
  EXPECT_EQ(11.0, testMath.dot(v1, v2));
}


// Test magnitude
TEST_F(MdisNacSensorModelTest, magnitudeZero) {
  std::vector<double> v(3, 0.0);
  EXPECT_EQ(0.0, testMath.magnitude(v));
}

TEST_F(MdisNacSensorModelTest, magnitudePositive) {
  std::vector<double>v(3);
  v[0] = 3.0;
  v[1] = 4.0;
  v[2] = 5.0;
  EXPECT_EQ(sqrt(50.0), testMath.magnitude(v));
}

TEST_F(MdisNacSensorModelTest, magnitudeNegative) {
  std::vector<double>v(3);
  v[0] = -3.0;
  v[1] = -4.0;
  v[2] = -5.0;
  EXPECT_EQ(sqrt(50.0), testMath.magnitude(v));
}


// Test normalize
TEST_F(MdisNacSensorModelTest, normalizeVector) {
  std::vector<double>v(3);
  v[0] = 2.0;
  v[1] = 3.0;
  v[2] = 4.0;
  std::vector<double> result = testMath.normalize(v);
  EXPECT_EQ((v[0] / sqrt(29.0)), result[0]);
  EXPECT_EQ((v[1] / sqrt(29.0)), result[1]);
  EXPECT_EQ((v[2] / sqrt(29.0)), result[2]);
}


// Test normalize
TEST_F(MdisNacSensorModelTest, distortionModel1) {
  double dx = 0.0;
  double dy = 0.0;
  double udx = 0.0;
  double udy = 0.0;
  double isis3_udx = 0.0;
  double isis3_udy = 0.0;
  testMath.undistortedFocalCoords(dx,dy,udx,udy);

  EXPECT_NEAR(udx,isis3_udx,tolerance);
  EXPECT_NEAR(udy,isis3_udy,tolerance);

}

TEST_F(MdisNacSensorModelTest, distortionModel2) {
  double dx = -6.30;
  double dy = 6.40;
  double udx = 0.0;
  double udy = 0.0;
  double isis3_udx = -6.3036234000160273893698104075156152248383;
  double isis3_udy = 6.3445144408882310216313271666876971721649;
  testMath.undistortedFocalCoords(dx,dy,udx,udy);

  EXPECT_NEAR(udx,isis3_udx,tolerance);
  EXPECT_NEAR(udy,isis3_udy,tolerance);

}



