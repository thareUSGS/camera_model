#include <gtest/gtest.h>

TEST(SimpleTest, TrueAssertion) {
  ASSERT_TRUE(0 == 0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
