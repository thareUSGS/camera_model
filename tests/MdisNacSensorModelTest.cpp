#include <MdisPlugin.h>
#include <MdisNacSensorModel.h>
#include <IsdReader.h>

#include <csm/Isd.h>

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
    using MdisNacSensorModel::computeElevation;
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
      dataFile = "./data/EN1007907102M.json";
      isd = readISD(dataFile);
      // printISD(*isd);
      
      // Make sure the isd was read correctly.
      if (isd == nullptr) {
        FAIL() << "Could not create isd from file: " << dataFile;
      }
      
      // Create a model from the ISD so we can test a valid image.
      std::string modelName = MdisNacSensorModel::_SENSOR_MODEL_NAME;
      csm::Model *validModel = mdisPlugin.constructModelFromISD(*isd, modelName);
      
      // We could static_cast, but may be hard to debug if it doesn't correctly cast.
      mdisModel = dynamic_cast<MdisNacSensorModel *>(validModel);
      std::cout << "Construction model: " << mdisModel << "\n";
      
      // Fatal failure if the downcast doesn't work
      if (mdisModel == nullptr) {
        FAIL() << "Could not downcast Model* to MdisNacSensorModel*.";
      }
    }
    
    static void TearDownTestCase() {
      delete isd;
      delete mdisModel;
    }
    
    static csm::Isd *isd;                 // ISD converted from JSON to use for creating model.
    static std::string dataFile;          // JSON data file to be converted to ISD for testing.
    static MdisPlugin mdisPlugin;         // Plugin used to create a model from ISD.
    static MdisNacSensorModel *mdisModel; // MDIS-NAC sensor model created with ISD.
    
    
    // Per test setup and teardown (e.g. each TEST_F)
    virtual void SetUp() {
      tolerance = 0.00001;
    }

    virtual void TearDown() {}
    
    double tolerance;                     // Tolerance to be used for double comparison.
    MdisNacSensorModel defaultMdisNac;    // A default constructed MdisNacSensorModel.
    TestableMdisNacSensorModel testMath;  // Subclassed MdisNacSensorModel for protected methods.
};


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
  double x = 1132.18*1000;
  double y = -1597.75*1000;
  double z = 1455.66*1000;
  csm::EcefCoord xyz(x, y, z);
  csm::ImageCoord pt = mdisModel->groundToImage(xyz);
  EXPECT_EQ(512.0, pt.line);
  EXPECT_EQ(512.0, pt.samp);
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



