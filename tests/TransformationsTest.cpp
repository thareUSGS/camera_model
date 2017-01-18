#include <math.h>
#include <transformations.h>
#include <Eigen/Dense>
#include <iostream>
#include <gtest/gtest.h>

class TransformationsTest : public ::testing::Test {
  protected:
    virtual void SetUp(){
      tolerance = 0.0001;
    }
    double tolerance;
};

TEST_F(TransformationsTest, OpkToRotation){
  // This is right out of Mikhail Modern Photogrammetry p.95
  Eigen::Matrix3d rot;
  rot << 0.9622, 0.2616, -0.0751,
         -0.2578, 0.9645, 0.0562,
         0.0871, -0.0348, 0.9956;

  float o = (2 * M_PI) / 180;
  float p = (5 * M_PI) / 180;
  float k = (15 * M_PI) / 180;
  ASSERT_TRUE(rot.isApprox(opkToRotation(o,p,k), tolerance));
}
