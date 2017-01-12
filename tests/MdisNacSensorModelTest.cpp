#include <MdisNacSensorModel.h>

#include <gtest/gtest.h>

// Set up a fixture (i.e. objects we can use throughout test)
class MdisNacSensorModelTest : public ::testing::Test {
  // May need this below to test private methods?
  // -- what is best way to test private methods, since they are not part of
  //    the class's public interface?
  // friend class MdisNacSensorModel;
  protected:
    virtual void SetUp() {
    }

    //virtual void TearDown() {}

    MdisNacSensorModel defaultMdisNac;
    //MdisNacSensorModel invalidIntersect;
};

// Tests the getModelState() method with a default constructed MdisNacSensorModel.
TEST_F(MdisNacSensorModelTest, getModelStateDefault) {
  EXPECT_EQ(defaultMdisNac.getModelState(), std::string());
}

// Just a test to see if gtest is set up correctly.
TEST_F(MdisNacSensorModelTest, getModelStateChanged) {
  EXPECT_EQ(defaultMdisNac.getModelState(), std::string("blah"));
}
