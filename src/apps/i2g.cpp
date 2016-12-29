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

  double xl= -6683432.694790688;
  double yl = -27264735.61516516;
  double zl = 8536481.287414147; // SENSOR Z from isd

  std::cout << "\nDEBUG\n";
  std::cout << "x =  " << x;
  std::cout << "\ny =  " << y;
  std::cout << "\nz =  " << z << std::endl;
  std::cout << "xl = " << xl;
  std::cout << "\nyl = " << yl;
  std::cout << "\nzl = " << zl;
  const double omega = 2.892112982708664;
  const double phi = 3.263192693345268;
  const double kappa = 149.1628863008221;

  const double sinw = sin(omega);
  const double sinp = sin(phi);
  const double sink = sin(kappa);
  const double cosw = cos(omega);
  const double cosp = cos(phi);
  const double cosk = cos(kappa);

  std::vector<double> rotationMatrix(9);
  rotationMatrix[0] = cosp * cosk;
  rotationMatrix[1] = cosw * sink + sinw * sinp * cosk;
  rotationMatrix[2] = sinw * sink - cosw * sinp * cosk;
  rotationMatrix[3] = -1 * cosp * sink;
  rotationMatrix[4] = cosw * cosk - sinw * sinp * sink;
  rotationMatrix[5] = sinw * cosk + cosw * sinp * sink;
  rotationMatrix[6] = sinp;
  rotationMatrix[7] = -1 * sinw * cosp;
  rotationMatrix[8] = cosw * cosp;                             

  std::vector<double> &m = rotationMatrix;

  const double f = 549.302734762479; // focal length
  const double x0 = 0.0; // principal point offset x
  const double y0 = 0.0; // principal point offset y
 
  double px = x0 + -f * (
                  ( m[0] * (x - xl) + m[1] * (y - yl) + m[2] * (z - zl) ) /
                  ( m[6] * (x - xl) + m[7] * (y - yl) + m[8] * (z - zl) ) );

  double py = x0 + -f * (
                  ( m[3] * (x - xl) + m[4] * (y - yl) + m[5] * (z - zl) ) /
                  ( m[6] * (x - xl) + m[7] * (y - yl) + m[8] * (z - zl) ) );

  std::cout << "\npx, py = " << px << ", " << py << std::endl;

  // Convert from image plane / focal plane (mm) to sample / line
  const double itranss[3] = { 0.0, 71.42857143, 0.0 };
  const double itransl[3] = { 0.0, 0.0, 71.42857143 };

  px = py = -0.007;

  double sample = itranss[0] + (itranss[1] * px + itranss[2] * py) + 512.0;
  double line = itransl[0] + (itransl[1] * px + itransl[2] * py) + 512.0;

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

// DO WE NEED TO USE THE DETECTOR LINE AND SAMPLE OFFSETS ABOVE?

// Rotation stuff
  // angles grabbed from isd (generated from mdis2isd)
  const double omega = 2.892112982708664;
  const double phi = 3.263192693345268;
  const double kappa = 149.1628863008221;

  const double sinw = sin(omega);
  const double sinp = sin(phi);
  const double sink = sin(kappa);
  const double cosw = cos(omega);
  const double cosp = cos(phi);
  const double cosk = cos(kappa);

  // Rotation matrix grabbed from Mikhail's "Introduction to Modern Photogrammetry"
  std::vector<double> rotationMatrix(9);
  rotationMatrix[0] = cosp * cosk;
  rotationMatrix[1] = cosw * sink + sinw * sinp * cosk;
  rotationMatrix[2] = sinw * sink - cosw * sinp * cosk;
  rotationMatrix[3] = -1 * cosp * sink;
  rotationMatrix[4] = cosw * cosk - sinw * sinp * sink;
  rotationMatrix[5] = sinw * cosk + cosw * sinp * sink;
  rotationMatrix[6] = sinp;
  rotationMatrix[7] = -1 * sinw * cosp;
  rotationMatrix[8] = cosw * cosp;                             

  std::vector<double> &m = rotationMatrix;

  // Note that elevation will be the height parameter in CSM's imageToGround.
  // For now, assume that ellipsoid has no terrain features
  double majorAxis = 6051800; // mercury ellipsoid

  double spacecraftX = -6683432.694790688; // SENSOR X
  double spacecraftY = -27264735.61516516; // SENSOR Y
  double spacecraftZ = 8536481.287414147; // SENSOR Z from isd

  double f = 549.3027347624796; // focal length (mm)
  double x = 1 * 
             ( m[0] * (focalPlaneX) + m[3] * (focalPlaneY) + m[6] * (-1 * f) ) /
             ( m[2] * (focalPlaneX) + m[5] * (focalPlaneY) + m[8] * (-1 * f) )
           + spacecraftX; // what coordinate system is this in?

  double y = 1 * 
             ( m[1] * (focalPlaneX) + m[4] * (focalPlaneY) + m[7] * (-1 * f) ) /
             ( m[2] * (focalPlaneX) + m[5] * (focalPlaneY) + m[8] * (-1 * f) )
           + spacecraftY; // what coordinate system is this in?

  std::cout << "\n(x, y) = (" << x/1000 << ", " << y/1000 << ", " << 1/1000 << ") (km)\n";
}
