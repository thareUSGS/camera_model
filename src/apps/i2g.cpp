#include <iostream>

void imageToGround(double sample, double line);

// Simple application to test the CSM imageToGround() method for MDIS-NAC
// i.e. prototyping
int main() {

  double sample = 512.0;
  double line = 512.0;

  std::cout << "Enter a sample: ";
  std::cin >> sample;
  std::cout << "Enter a line: ";
  std::cin >> line;

  imageToGround(sample, line);

  return 0;

}

void imageToGround(double sample, double line) {
// Convert image sample/line to focal plane x/y
//     // Perform affine transformation to get from sample/line to x/y (mm)
  const double transX[3] = {0.0, 0.014, 0.0}; // JSON transx
  const double transY[3] = {0.0, 0.0, 0.014}; // JSON transy

  // center the sample line
  // imgSamps / 2, imgLines / 2: 1024x1024
  sample -= 512.0; // ISD needs a center sample in CSM coord
  line -= 512.0; // ISD needs a center line in CSM coord (.5 .5 pixel centers)

  double focalPlaneX = transX[0] + (transX[1] * sample)
                                 + (transX[2] * line);
  double focalPlaneY = transY[0] + (transY[1] * sample)
                                 + (transY[2] * line);
  std::cout << "\nFOCAL PLANE (X, Y) mm: (" << focalPlaneX << ", " << focalPlaneY << ")\n";
//
}
