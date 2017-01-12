#include <MdisPlugin.h>
#include <MdisNacSensorModel.h>

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
    using MdisNacSensorModel::intersect;
    using MdisNacSensorModel::perpendicular;
    using MdisNacSensorModel::project;
    using MdisNacSensorModel::dot;
    using MdisNacSensorModel::magnitude;
    using MdisNacSensorModel::normalize;
};


// Set up a fixture (i.e. objects we can use throughout test)
class MdisNacSensorModelTest : public ::testing::Test {
  protected:
    // Per test-case setup and teardown (e.g. once for this MdisNacSensorModelTest)
    static void SetUpTestCase() {
      isd = new csm::Isd();
      isd->addParam("boresight", "0.0");
      isd->addParam("boresight", "0.0");
      isd->addParam("boresight", "1.0");
      isd->addParam("ccd_center", "512.5");
      isd->addParam("ephemeris_time", "418855170.493");
      isd->addParam("focal_length", "549.117819537");
      isd->addParam("ifov", "25.44");
      isd->addParam("instrument_id", "MDIS_NAC");
      isd->addParam("itrans_line", "0.0");
      isd->addParam("itrans_line", "0.0");
      isd->addParam("itrans_line", "71.42857143");
      isd->addParam("itrans_sample", "0.0");
      isd->addParam("itrans_sample", "71.42857143");
      isd->addParam("itrans_sample", "0.0");
      isd->addParam("nlines", "1024");
      isd->addParam("nsamples", "1024");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_x", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("odt_y", "0.0");
      isd->addParam("omega", "2.33344166314");
      isd->addParam("phi", "0.606364919159");
      isd->addParam("kappa", "0.704740691838");
      isd->addParam("original_half_samples", "512");
      isd->addParam("original_half_lines", "512");
      isd->addParam("pixel_pitch", "0.014");
      isd->addParam("semi_major_axis", "2439.4");
      isd->addParam("semi_minor_axis", "2439.4");
      isd->addParam("spacecraft_name", "Messenger");
      isd->addParam("starting_detector_line", "1");
      isd->addParam("starting_detector_sample", "9");
      isd->addParam("target_name", "Mercury");
      isd->addParam("transx", "0.0");
      isd->addParam("transx", "0.014");
      isd->addParam("transx", "0.0");
      isd->addParam("transy", "0.0");
      isd->addParam("transy", "0.0");
      isd->addParam("transy", "0.014");
      isd->addParam("x_sensor_origin", "1728357.70312");
      isd->addParam("y_sensor_origin", "-2088409.0061");
      isd->addParam("z_sensor_origin", "2082873.92806");
    }
    
    static void TearDownTestCase() {
      delete isd;
      isd = NULL;
    }
    
    // Per test setup and teardown (e.g. each TEST_F)
    virtual void SetUp() {
      tolerance = 0.00001;
    }

    virtual void TearDown() {}

    static csm::Isd *isd;
    
    double tolerance;
    MdisPlugin mdisPlugin;
    MdisNacSensorModel defaultMdisNac;
    TestableMdisNacSensorModel testMath;
};


csm::Isd *MdisNacSensorModelTest::isd = NULL;


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
  
  // Create a model from the ISD so we can test a valid image.
  std::string modelName = MdisNacSensorModel::_SENSOR_MODEL_NAME;
  csm::Model *validModel = mdisPlugin.constructModelFromISD(*isd, modelName);
  // We could static_cast, but may be hard to debug if it doesn't correctly cast.
  MdisNacSensorModel *mdisModel = dynamic_cast<MdisNacSensorModel *>(validModel);
  
  // Fatal failure if the downcast doesn't work
  if (!mdisModel) {
    FAIL() << "Could not downcast Model* to MdisNacSensorModel*.";
  }
  
  csm::EcefCoord xyz = mdisModel->imageToGround(point, height);
  double truth[] = { 1132.18*1000, -1597.75*1000, 1455.66*1000 };
  EXPECT_EQ(truth[0], xyz.x);
  EXPECT_EQ(truth[1], xyz.y);
  EXPECT_EQ(truth[2], xyz.z);
  
  // Remove the memory we took ownership of.
  delete mdisModel;
}


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
