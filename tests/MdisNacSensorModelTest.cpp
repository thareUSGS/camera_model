#include <MdisNacSensorModel.h>

#include <gtest/gtest.h>


/**
 * Sub-class MdisNacSensorModel to get test its protected linear algebra methods.
 * 
 * We should be testing the protected methods of MdisNacSensorModel since imageToGround
 * depends on intersect, which depends on project, etc.
 */
class TestableMdisNacSensorModel : public MdisNacSensorModel {
  // Give linear algebra methods public accessing when using instances of this class.
  public:
    using MdisNacSensorModel::intersect;
    using MdisNacSensorModel::perpendicular;
    using MdisNacSensorModel::project;
    using MdisNacSensorModel::dot;
    using MdisNacSensorModel::magnitude;
    using MdisNacSensorModel::normalize;
    using MdisNacSensorModel::setFocalPlane;
};


// Set up a fixture (i.e. objects we can use throughout test)
class MdisNacSensorModelTest : public ::testing::Test {
  protected:
    // Per test-case setup and teardown (e.g. once for this MdisNacSensorModelTest)
    static void SetUpTestCase() {
    }
    
    static void TearDownTestCase() {}
    
    // Per test setup and teardown (e.g. each TEST_F)
    virtual void SetUp() {
      tolerance = 0.00001;
    }

    virtual void TearDown() {}

    MdisNacSensorModel defaultMdisNac;
    double tolerance;
    TestableMdisNacSensorModel testMath;
};


// Tests the getModelState() method with a default constructed MdisNacSensorModel.
TEST_F(MdisNacSensorModelTest, getModelStateDefault) {
  EXPECT_EQ(defaultMdisNac.getModelState(), std::string());
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
  testMath.setFocalPlane(dx,dy,udx,udy);

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
  testMath.setFocalPlane(dx,dy,udx,udy);

  EXPECT_NEAR(udx,isis3_udx,tolerance);
  EXPECT_NEAR(udy,isis3_udy,tolerance);

}



