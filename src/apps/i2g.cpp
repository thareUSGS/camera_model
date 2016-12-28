#include <cmath>
#include <iostream>
#include <vector>

void imageToGround(double sample, double line);
void groundToImage(double x, double y, double z);

// Simple application to test the CSM imageToGround() method for MDIS-NAC
// i.e. prototyping
int main() {

  double x, y, z = 0.0;

  std::cout << "X: ";
  std::cin >> x;
  std::cout << "Y: ";
  std::cin >> y;
  std::cout << "Z: ";
  std::cin >> z;
  x *= 1000;
  y *= 1000;
  z *= 1000;

  groundToImage(x,y,z);

  double sample = 512.0;
  double line = 512.0;

  std::cout << "Enter a sample: ";
  std::cin >> sample;
  std::cout << "Enter a line: ";
  std::cin >> line;

  imageToGround(sample, line);

  return 0;

}


void groundToImage(double x, double y, double z) {

  double spacecraftX = -6683432.694790688;
  double spacecraftY = -27264735.61516516;
  double spacecraftZ = 8536481.287414147; // SENSOR Z from isd
  const double omega = 2.892112982708664;
  const double phi = 3.263192693345268;
  const double kappa = 149.1628863008221;

  std::vector<double> rotationMatrix(9);
  rotationMatrix[0] = cos(phi) * cos(kappa);
  rotationMatrix[1] = cos(omega) * sin(kappa) + sin(omega) * sin(phi) * cos(kappa);
  rotationMatrix[2] = sin(omega) * sin(kappa) - cos(omega) * sin(phi) * cos(kappa);
  rotationMatrix[3] = -1 * cos(phi) * sin(kappa); 
  rotationMatrix[4] = cos(omega) * cos(kappa) - sin(omega) * sin(phi) * sin(kappa);
  rotationMatrix[5] = sin(omega) * cos(kappa) + cos(omega) * sin(phi) * sin(kappa);
  rotationMatrix[6] = sin(phi);
  rotationMatrix[7] = -1 * sin(omega) * cos(phi);
  rotationMatrix[8] = cos(omega) * cos(phi);

  double f = 549.302734762479; // focal length

  double sample = -f * (
                  ( rotationMatrix[0] * (x - spacecraftX) + rotationMatrix[1] * (y - spacecraftY) + rotationMatrix[2] * (z - spacecraftZ) ) /
                  ( rotationMatrix[6] * (x - spacecraftX) + rotationMatrix[7] * (y - spacecraftY) + rotationMatrix[8] * (z - spacecraftZ) ) )
                + 512.0;

  double line = -f * (
                  ( rotationMatrix[3] * (x - spacecraftX) + rotationMatrix[4] * (y - spacecraftY) + rotationMatrix[5] * (z - spacecraftZ) ) /
                  ( rotationMatrix[6] * (x - spacecraftX) + rotationMatrix[7] * (y - spacecraftY) + rotationMatrix[8] * (z - spacecraftZ) ) )
                + 512.0;

  std::cout << "\nsample, line = " << sample << ", " << line << std::endl;
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

// Rotation stuff
  // angles grabbed from isd (generated from mdis2isd)
  const double omega = 2.892112982708664;
  const double phi = 3.263192693345268;
  const double kappa = 149.1628863008221;

  // Rotation matrix grabbed from Mikhail's "Introduction to Modern Photogrammetry"
  std::vector<double> rotationMatrix(9);
  rotationMatrix[0] = cos(phi) * cos(kappa);
  rotationMatrix[1] = cos(omega) * sin(kappa) + sin(omega) * sin(phi) * cos(kappa);
  rotationMatrix[2] = sin(omega) * sin(kappa) - cos(omega) * sin(phi) * cos(kappa);
  rotationMatrix[3] = -1 * cos(phi) * sin(kappa); 
  rotationMatrix[4] = cos(omega) * cos(kappa) - sin(omega) * sin(phi) * sin(kappa);
  rotationMatrix[5] = sin(omega) * cos(kappa) + cos(omega) * sin(phi) * sin(kappa);
  rotationMatrix[6] = sin(phi);
  rotationMatrix[7] = -1 * sin(omega) * cos(phi);
  rotationMatrix[8] = cos(omega) * cos(phi);

  // Note that elevation will be the height parameter in CSM's imageToGround.
  // For now, assume that ellipsoid has no terrain features
  double elevation = 6051800; // mercury ellipsoid
  double spacecraftZ = 8536481.287414147; // SENSOR Z from isd
  double spacecraftAltitude = elevation - spacecraftZ;
  double z = spacecraftAltitude; // temp

  double spacecraftX = -6683432.694790688;
  double spacecraftY = -27264735.61516516;

  double f = 549.3027347624796 * 1000; // focal length
  double x = z * 
             ( rotationMatrix[0] * (sample - focalPlaneX) + rotationMatrix[3] * (line - focalPlaneY) + rotationMatrix[6] * (-1 * f) ) /
             ( rotationMatrix[2] * (sample - focalPlaneX) + rotationMatrix[5] * (line - focalPlaneY) + rotationMatrix[8] * (-1 * f) )
           + spacecraftX;

  double y = z * 
             ( rotationMatrix[1] * (sample - focalPlaneX) + rotationMatrix[4] * (line - focalPlaneY) + rotationMatrix[7] * (-1 * f) ) /
             ( rotationMatrix[2] * (sample - focalPlaneX) + rotationMatrix[5] * (line - focalPlaneY) + rotationMatrix[8] * (-1 * f) )
           + spacecraftY;


  std::cout << "\n(x, y) = (" << x/1000 << ", " << y/1000 << ", " << z/1000 << ") (km)\n";
}
