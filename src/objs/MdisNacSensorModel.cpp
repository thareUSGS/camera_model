#include "MdisNacSensorModel.h"

#include <iomanip>
#include <iostream>

#include <csm/Error.h>

using namespace std;

const std::string MdisNacSensorModel::_SENSOR_MODEL_NAME 
                                      = "ISIS_MDISNAC_USGSAstro_1_Linux64_csm30.so";


                                      
MdisNacSensorModel::MdisNacSensorModel() {

  m_transX[0] = 0.0;
  m_transX[1] = 0.0;
  m_transX[2] = 0.0;

  m_transY[0] = 0.0;
  m_transY[1] = 0.0;
  m_transY[2] = 0.0;
  
  m_majorAxis = 0.0;
  m_minorAxis = 0.0;
  m_omega = 0.0;
  m_phi = 0.0;
  m_kappa = 0.0;
  m_focalLength = 0.0;
  
  m_spacecraftPosition[0] = 0.0;
  m_spacecraftPosition[1] = 0.0;
  m_spacecraftPosition[2] = 0.0;
  
  m_ccdCenter = 0.0;

  m_startingDetectorSample = 0.0;
  m_startingDetectorLine = 0.0;
  m_targetName = "";
  m_ifov = 0.0;
  m_instrumentID = "";
  m_focalLengthEpsilon = 0.0;
  

#if 0
  //NAC coefficients 
  m_odtX[0]=0.0;
  m_odtX[1]=1.0018542696237999756;
  m_odtX[2]=-0.0;
  m_odtX[3]=-0.0;
  m_odtX[4]=-0.00050944404749411103042;
  m_odtX[5]=0.0;
  m_odtX[6]=1.0040104714688599425e-05;
  m_odtX[7]=0.0;
  m_odtX[8]=1.0040104714688599425e-05;
  m_odtX[9]=0.0;

  m_odtY[0]=0.0;
  m_odtY[1]=0.0;
  m_odtY[2]=1.0;
  m_odtY[3]=0.00090600105949967496381;
  m_odtY[4]=0.0;
  m_odtY[5]=0.00035748426266207598964;
  m_odtY[6]=0.0;
  m_odtY[7]=1.0040104714688599425e-05;
  m_odtY[8]=0.0;
  m_odtY[9]=1.0040104714688599425e-05;
#endif
  m_odtX[0] = 0.0;
  m_odtX[1] = 0.0;
  m_odtX[2] = 0.0;
  m_odtX[3] = 0.0;
  m_odtX[4] = 0.0;
  m_odtX[5] = 0.0;
  m_odtX[6] = 0.0;
  m_odtX[7] = 0.0;
  m_odtX[8] = 0.0;
  m_odtX[9] = 0.0;

  m_odtY[0] = 0.0;
  m_odtY[1] = 0.0;
  m_odtY[2] = 0.0;
  m_odtY[3] = 0.0;
  m_odtY[4] = 0.0;
  m_odtY[5] = 0.0;
  m_odtY[6] = 0.0;
  m_odtY[7] = 0.0;
  m_odtY[8] = 0.0;
  m_odtY[9] = 0.0;

  m_originalHalfLines = 0.0;
  m_spacecraftName = "";
  m_pixelPitch = 0.0; 
    
  m_iTransS[0] = 0.0;
  m_iTransS[1] = 0.0;
  m_iTransS[2] = 0.0;
  
  m_iTransL[0] = 0.0;
  m_iTransL[1] = 0.0;
  m_iTransL[2] = 0.0;

  m_ephemerisTime = 0.0;
  m_originalHalfSamples = 0.0;
  m_boresight[0] = 0.0;
  m_boresight[1] = 0.0;
  m_boresight[2] = 0.0;

  m_nLines = 0;
  m_nSamples = 0;
}


MdisNacSensorModel::~MdisNacSensorModel() {}


/**
 * @brief Compute undistorted focal plane x/y.
 *
 * Computes undistorted focal plane (x,y) coordinates given a distorted focal plane (x,y)
 * coordinate. The undistorted coordinates are solved for using the Newton-Raphson
 * method for root-finding if the distortionFunction method is invoked.
 *
 * @param dx distorted focal plane x in millimeters
 * @param dy distorted focal plane y in millimeters
 * @param undistortedX The undistorted x coordinate, in millimeters.
 * @param undistortedY The undistorted y coordinate, in millimeters.
 *
 * @return if the conversion was successful
 * @todo Review the tolerance and maximum iterations of the root-
 *       finding algorithm.
 * @todo Review the handling of non-convergence of the root-finding
 *       algorithm.
 * @todo Add error handling for near-zero determinant.
*/
bool MdisNacSensorModel::setFocalPlane(double dx,double dy,
                                       double &undistortedX,
                                       double &undistortedY ) const {


  // Solve the distortion equation using the Newton-Raphson method.
  // Set the error tolerance to about one millionth of a NAC pixel.
  const double tol = 1.4E-5;

  // The maximum number of iterations of the Newton-Raphson method.
  const int maxTries = 60;

  double x;
  double y;
  double fx;
  double fy;
  double Jxx;
  double Jxy;
  double Jyx;
  double Jyy;

  // Initial guess at the root
  x = dx;
  y = dy;

  distortionFunction(x, y, fx, fy);

  for (int count = 1; ((fabs(fx) + fabs(fy)) > tol) && (count < maxTries); count++) {

    this->distortionFunction(x, y, fx, fy);

    fx = dx - fx;
    fy = dy - fy;

    distortionJacobian(x, y, Jxx, Jxy, Jyx, Jyy);

    double determinant = Jxx * Jyy - Jxy * Jyx;
    if (determinant < 1E-6) {
      //
      // Near-zero determinant. Add error handling here.
      //
      //-- Just break out and return with no convergence
      break;
    }

    x = x + (Jyy * fx - Jxy * fy) / determinant;
    y = y + (Jxx * fy - Jyx * fx) / determinant;
  }

  if ( (fabs(fx) + fabs(fy)) <= tol) {
    // The method converged to a root.
    undistortedX = x;
    undistortedY = y;

  }
  else {
    // The method did not converge to a root within the maximum
    // number of iterations. Return with no distortion.
    undistortedX = dx;
    undistortedY = dy;
  }

  return true;

}


/**
 * @description Jacobian of the distortion function. The Jacobian was computed
 * algebraically from the function described in the distortionFunction
 * method.
 *
 * @param x
 * @param y
 * @param Jxx  Partial_xx
 * @param Jxy  Partial_xy
 * @param Jyx  Partial_yx
 * @param Jyy  Partial_yy
 */
void MdisNacSensorModel::distortionJacobian(double x, double y, double &Jxx, double &Jxy,
                                            double &Jyx, double &Jyy) const {

  double d_dx[10];
  d_dx[0] = 0;
  d_dx[1] = 1;
  d_dx[2] = 0;
  d_dx[3] = 2 * x;
  d_dx[4] = y;
  d_dx[5] = 0;
  d_dx[6] = 3 * x * x;
  d_dx[7] = 2 * x * y;
  d_dx[8] = y * y;
  d_dx[9] = 0;
  double d_dy[10];
  d_dy[0] = 0;
  d_dy[1] = 0;
  d_dy[2] = 1;
  d_dy[3] = 0;
  d_dy[4] = x;
  d_dy[5] = 2 * y;
  d_dy[6] = 0;
  d_dy[7] = x * x;
  d_dy[8] = 2 * x * y;
  d_dy[9] = 3 * y * y;

  Jxx = 0.0;
  Jxy = 0.0;
  Jyx = 0.0;
  Jyy = 0.0;

  for (int i = 0; i < 10; i++) {
    Jxx = Jxx + d_dx[i] * m_odtX[i];
    Jxy = Jxy + d_dy[i] * m_odtX[i];
    Jyx = Jyx + d_dx[i] * m_odtY[i];
    Jyy = Jyy + d_dy[i] * m_odtY[i];
  }


}


/**
 * @description Compute distorted focal plane (dx,dy) coordinate  given an undistorted focal
 * plane (ux,uy) coordinate. This describes the third order Taylor approximation to the
 * distortion model.
 *
 * @param ux Undistored x
 * @param uy Undistored y
 * @param dx Result distorted x
 * @param dy Result distorted y
 */
void MdisNacSensorModel::distortionFunction(double ux, double uy, double &dx, double &dy) const {

  double f[10];
  f[0] = 1;
  f[1] = ux;
  f[2] = uy;
  f[3] = ux * ux;
  f[4] = ux * uy;
  f[5] = uy * uy;
  f[6] = ux * ux * ux;
  f[7] = ux * ux * uy;
  f[8] = ux * uy * uy;
  f[9] = uy * uy * uy;

  dx = 0.0;
  dy = 0.0;

  for (int i = 0; i < 10; i++) {
    dx = dx + f[i] * m_odtX[i];
    dy = dy + f[i] * m_odtY[i];
  }

}


csm::ImageCoord MdisNacSensorModel::groundToImage(const csm::EcefCoord &groundPt, 
                                                  double desiredPrecision, 
                                                  double *achievedPrecision, 
                                                  csm::WarningList *warnings) const {
  // Get the rotation matrix from ISD supplied omega, phi, kappa
  std::vector<double> rotation = createRotationMatrix(m_omega, m_phi, m_kappa);
  
  // Determine the center look direction (centerX, centerY, f) of the sensor in body-fixed
  // (e.g. the sensor's optical axis rotated into body-fixed frame)
  double centerX = 0.0;
  double centerY = 0.0;
  std::vector<double> sensorLookC { 0.0, 0.0, m_focalLength };
  std::vector<double> sensorLookB = rotate(sensorLookC, rotation);
  
  // Find the look vector from the sensor center look to the ground point in body-fixed
  std::vector<double> lookGroundB {
    groundPt.x - m_spacecraftPosition[0],
    groundPt.y - m_spacecraftPosition[1],
    groundPt.z - m_spacecraftPosition[2],
  };
  
  // Normalize by sacling the sensor-to-ground look to the sensor center look in body-fixed
  double sensorLookBMag = sqrt( dot(sensorLookB, sensorLookB) );
  double lookGroundBMag = sqrt( dot(lookGroundB, lookGroundB) );
  double scale = sensorLookBMag / lookGroundBMag;
  lookGroundB[0] = scale * lookGroundB[0]; 
  lookGroundB[1] = scale * lookGroundB[1]; 
  lookGroundB[2] = scale * lookGroundB[2]; 
    
  // Rotate the sensor-to-ground look vector from body-fixed to sensor frame (inverse rotation)
  std::vector<double> lookGroundC = rotate(lookGroundB, rotation, true);
    
  // Scale the sensor-to-ground sensor frame vector so that it intersects the focal plane
  // (i.e. scale it so its Z component equals the sensor's focal length)
  scale = m_focalLength / lookGroundC[2];
  double focalPlaneX = lookGroundC[0] * scale;
  double focalPlaneY = lookGroundC[1] * scale;
  
  // Distortion
  
  // Convert focal plane mm to pixels
  double pixelX = focalPlaneX * (1.0 / m_transX[1]);
  double pixelY = focalPlaneY * (1.0 / m_transY[2]);
  
  // Convert pixels to line,sample
  double sample = pixelX + m_ccdCenter - 0.5;
  double line = pixelY + m_ccdCenter - 0.5;
  
  // Check if the computed line,sample coordinate is within the image dimensions
  bool outOfBounds = false;
  if (sample > m_nSamples || sample < 0.0 || line > m_nLines || line < 0.0) {
    outOfBounds = true;
  }
  if (outOfBounds) {
    if (warnings != nullptr) {
      std::string msg("The image coordinate is outside the image dimensions.");
      std::string func("MdisNacSensorModel::imageToGround");
      warnings->push_front(csm::Warning(csm::Warning::IMAGE_COORD_OUT_OF_BOUNDS,
                                        msg,
                                        func));
    }
  }
  
  return csm::ImageCoord(line, sample);
}
       
       
csm::ImageCoordCovar MdisNacSensorModel::groundToImage(const csm::EcefCoordCovar &groundPt, 
                                   double desiredPrecision, 
                                   double *achievedPrecision, 
                                   csm::WarningList *warnings) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::groundToImage");
}
                
                
csm::EcefCoord MdisNacSensorModel::imageToGround(const csm::ImageCoord &imagePt, 
                                                 double height, 
                                                 double desiredPrecision, 
                                                 double *achievedPrecision, 
                                                 csm::WarningList *warnings) const {

  double sample = imagePt.samp;
  double line = imagePt.line;

  // center the sample line
  sample -= m_ccdCenter - 0.5; // ISD needs a center sample in CSM coord
  line -= m_ccdCenter - 0.5; // ISD needs a center line in CSM coord (.5 .5 pixel centers)
  
  // Convert from sample/line to focal plane coordinates (in mm)
  double focalPlaneX = m_transX[0] + (m_transX[1] * sample) + (m_transX[2] * line);
  double focalPlaneY = m_transY[0] + (m_transY[1] * sample) + (m_transY[2] * line);

  double undistortedFocalPlaneX = focalPlaneX;
  double undistortedFocalPlaneY = focalPlaneY;

  setFocalPlane(focalPlaneX, focalPlaneY, undistortedFocalPlaneX, undistortedFocalPlaneY);

  // Trigonometric functions for rotation matrix
  const double sinw = std::sin(m_omega);
  const double sinp = std::sin(m_phi);
  const double sink = std::sin(m_kappa);
  const double cosw = std::cos(m_omega);
  const double cosp = std::cos(m_phi);
  const double cosk = std::cos(m_kappa);

  // Rotation matrix taken from Introduction to Mordern Photogrammetry by 
  // Edward M. Mikhail, et al., p. 373
  std::vector<double> m(9);
  m[0] = cosp * cosk;
  m[1] = cosw * sink + sinw * sinp * cosk;
  m[2] = sinw * sink - cosw * sinp * cosk;
  m[3] = -1 * cosp * sink;
  m[4] = cosw * cosk - sinw * sinp * sink;
  m[5] = sinw * cosk + cosw * sinp * sink;
  m[6] = sinp;
  m[7] = -1 * sinw * cosp;
  m[8] = cosw * cosp;
  
  // Multiply the focal vector to the rotation matrix to get the direction of the camera
  std::vector<double> direction(3);
  direction[0] = m[0] * undistortedFocalPlaneX + m[1] * undistortedFocalPlaneY + m[2] * m_focalLength;
  direction[1] = m[3] * undistortedFocalPlaneX + m[4] * undistortedFocalPlaneY + m[5] * m_focalLength;
  direction[2] = m[6] * undistortedFocalPlaneX + m[7] * undistortedFocalPlaneY + m[8] * m_focalLength;
  
  // Save the spacecraft position as a vector
  std::vector<double> spacecraftPosition(3);
  spacecraftPosition[0] = m_spacecraftPosition[0];
  spacecraftPosition[1] = m_spacecraftPosition[1];
  spacecraftPosition[2] = m_spacecraftPosition[2];
  
  // Perform the intersection
  return intersect(spacecraftPosition, direction, m_majorAxis);
}


double MdisNacSensorModel::computeElevation(double x, double y, double z) const {
  // For now, assume a sphere for Mercury. This will change for ellispoids/DEMs.
  // (radius + elevation)^2 = x^2 + y^2 + z^2
  // double localRadius = sqrt(x*x + y*y + z*z);
  // return localRadius - m_majorAxis;
  return 0;
}


csm::EcefCoord MdisNacSensorModel::intersect(const std::vector<double> &sensorPosition, 
                                             const std::vector<double> &direction, 
                                             double radius) const {
  
  std::vector<double> intersectionPoint(3), unitSensorPosition(3), unitDirection(3);
  
  // Transform the look direction vector into unit sphere space.
  unitDirection[0] = direction[0] / radius;
  unitDirection[1] = direction[1] / radius;
  unitDirection[2] = direction[2] / radius;
  
  // Transform the spacecraft position vector into unit sphere space.
  unitSensorPosition[0] = sensorPosition[0] / radius;
  unitSensorPosition[1] = sensorPosition[1] / radius;
  unitSensorPosition[2] = sensorPosition[2] / radius;
  
  // Find this vector:
  // If you extend the look direction infinitely, find the vector that is
  // perpendicular to that look direction and points to the body's origin.
  std::vector<double> perpendicularV = perpendicular(unitSensorPosition, unitDirection);
  
  // The look direction vector extended to the perpendicular vector.
  std::vector<double> positionProj(3);
  positionProj[0] = unitSensorPosition[0] - perpendicularV[0];
  positionProj[1] = unitSensorPosition[1] - perpendicularV[1];
  positionProj[2] = unitSensorPosition[2] - perpendicularV[2];
  
  // Find magnitudes of the "unit" spacecraft position vector and perpendicular vector.
  double positionMag = magnitude(unitSensorPosition);
  double perpendicularMag = magnitude(perpendicularV);
  
  // Use max-norm (infinity-norm) to normalize the "unit" look direction vector.
  std::vector<double> unitDirectionNorm = normalize(unitDirection);
  
  // Positive sign indicates spacecraft is in the target body,
  // negative indicates spacecraft is outside the target body.
  int sign = 1;
  
  // If the spacecraft position is outside the target body
  if (positionMag > 1.0) {
    // If vector perpendicular to look direction has magnitude > 1,
    // we are not looking at the target body (since the body is a "unit sphere")
    if (perpendicularMag > 1.0) {
      csm::EcefCoord coord(intersectionPoint[0], intersectionPoint[1], intersectionPoint[2]);
      return coord; // empty
    }
    
    // If looking away from the target body, there is no intersection.
    if (dot(positionProj, unitDirection) > 0.0) {
      csm::EcefCoord coord(intersectionPoint[0], intersectionPoint[1], intersectionPoint[2]);
      return coord; // empty
    }
    
    // If intersection point is on the limb, then transform back to target body size.
    if (perpendicularMag == 1.0) {
      intersectionPoint[0] = perpendicularV[0] * radius;
      intersectionPoint[1] = perpendicularV[1] * radius;
      intersectionPoint[2] = perpendicularV[2] * radius;
      csm::EcefCoord coord(intersectionPoint[0], intersectionPoint[1], intersectionPoint[2]);
      return coord;
    }
    
    sign = -1;
  }
  // If the spacecraft is on the target body
  else if (positionMag == 1.0) {
    csm::EcefCoord coord(sensorPosition[0], sensorPosition[1], sensorPosition[2]);
    return coord;
  }
  // If the spacecraft is inside the target body??? (target is sky?)
  else {
    sign = 1;
  }
  
  // To visualize:
  // There exists a scalar, scale, such that 
  // ||perpendicularV||^2 + scale^2 = radius^2.
  // Since we are in unit sphere space, the radius = 1.
  double scale = 1 - perpendicularMag * perpendicularMag;
  if (scale < 0.0) {
    scale = 0.0;
  }
  scale = sqrt(scale);
  
  // Find the intersection: 
  /*
   *              perpendicularV 
   * intersect o<-----^
   *          r \     |  scale * normalized look direction * sign
   *           a \    |  
   *            d \   |
   *             i \  |       
   *              u \ ^
   *               s \|  normalized look direction * sign
   *                  o  center of target body
   */
  // We can multiply the scale by the unit look vector and add the perpendicular vector
  // to get the intersection point
  intersectionPoint[0] = perpendicularV[0] + sign * scale * unitDirectionNorm[0];
  intersectionPoint[1] = perpendicularV[1] + sign * scale * unitDirectionNorm[1];
  intersectionPoint[2] = perpendicularV[2] + sign * scale * unitDirectionNorm[2];
  
  // Rescale the intersectionPoint from unit sphere space to ellipsoid space
  intersectionPoint[0] *= radius;
  intersectionPoint[1] *= radius;
  intersectionPoint[2] *= radius;
  
  csm::EcefCoord coord(intersectionPoint[0], intersectionPoint[1], intersectionPoint[2]);

  return coord;
}


std::vector<double> MdisNacSensorModel::perpendicular(const std::vector<double> &v1, 
                                                      const std::vector<double> &v2) const {
  
  double max1 = 0.0;
  double max2 = 0.0;
  
  // Get the infinite norm for each vector.
  for (int i = 0; i < 3; i++) {
    if (std::abs(v1[i]) > max1) {
      max1 = std::abs(v1[i]);
    }
    if (std::abs(v2[i]) > max2) {
      max2 = std::abs(v2[i]);
    }
  }
  
  // Scale the vectors by their max norms (not sure why this is needed, optimization maybe?).
  std::vector<double> newV1(3), newV2(3);
  
  newV1[0] = v1[0] / max1;
  newV1[1] = v1[1] / max1;
  newV1[2] = v1[2] / max1;
  
  newV2[0] = v2[0] / max2;
  newV2[1] = v2[1] / max2;
  newV2[2] = v2[2] / max2;
  
  // Project first vector onto second vector.
  std::vector<double> parallelV = project(newV1, newV2);
  
  std::vector<double> perpendicularV(3);
  
  // Get the perpendicular vector by subtracting the first vector by its projection
  // onto the second.
  perpendicularV[0] = v1[0] - parallelV[0] * max1;
  perpendicularV[1] = v1[1] - parallelV[1] * max1;
  perpendicularV[2] = v1[2] - parallelV[2] * max1;
  
  return perpendicularV;
}


std::vector<double> MdisNacSensorModel::project(const std::vector<double> &v1, 
                                                const std::vector<double> &v2) const {
                                                  
  double v1Dotv2 = dot(v1, v2);
  double v2Dotv2 = dot(v2, v2);
  std::vector<double> projV(3);
  projV[0] = v1Dotv2 / v2Dotv2 * v2[0];
  projV[1] = v1Dotv2 / v2Dotv2 * v2[1];
  projV[2] = v1Dotv2 / v2Dotv2 * v2[2];
  return projV;
}


double MdisNacSensorModel::dot(const std::vector<double> &v1, 
                               const std::vector<double> &v2) const {
                                 
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


double MdisNacSensorModel::magnitude(const std::vector<double> &v) const {
  
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


std::vector<double> MdisNacSensorModel::normalize(const std::vector<double> &v) const {
  
  double mag = magnitude(v);
  std::vector<double> n(3);
  n[0] = v[0] / mag;
  n[1] = v[1] / mag;
  n[2] = v[2] / mag;
  return n;
}


std::vector<double> MdisNacSensorModel::createRotationMatrix(const double omega,
                                                             const double phi,
                                                             const double kappa) const {
  // Trigonometric functions for rotation matrix
  const double sinw = std::sin(omega);
  const double sinp = std::sin(phi);
  const double sink = std::sin(kappa);
  const double cosw = std::cos(omega);
  const double cosp = std::cos(phi);
  const double cosk = std::cos(kappa);
  
  std::vector<double> m(9);
  m[0] = cosp * cosk;
  m[1] = cosw * sink + sinw * sinp * cosk;
  m[2] = sinw * sink - cosw * sinp * cosk;
  m[3] = -1 * cosp * sink;
  m[4] = cosw * cosk - sinw * sinp * sink;
  m[5] = sinw * cosk + cosw * sinp * sink;
  m[6] = sinp;
  m[7] = -1 * sinw * cosp;
  m[8] = cosw * cosp;
  
  return m;
}


std::vector<double> MdisNacSensorModel::rotate(const std::vector<double> &v,
                                               const std::vector<double> &rotationMatrix,
                                               bool invert) const {
  std::vector<double> rotated(3);
  if (!invert) {
    rotated[0] = rotationMatrix[0] * v[0] + rotationMatrix[1] * v[1] + rotationMatrix[2] * v[2];
    rotated[1] = rotationMatrix[3] * v[0] + rotationMatrix[4] * v[1] + rotationMatrix[5] * v[2];
    rotated[2] = rotationMatrix[6] * v[0] + rotationMatrix[7] * v[1] + rotationMatrix[8] * v[2];
  }
  else {
    rotated[0] = rotationMatrix[0] * v[0] + rotationMatrix[3] * v[1] + rotationMatrix[6] * v[2];
    rotated[1] = rotationMatrix[1] * v[0] + rotationMatrix[4] * v[1] + rotationMatrix[7] * v[2];
    rotated[2] = rotationMatrix[2] * v[0] + rotationMatrix[5] * v[1] + rotationMatrix[8] * v[2];
  }
  return rotated;
}


csm::EcefCoordCovar MdisNacSensorModel::imageToGround(const csm::ImageCoordCovar &imagePt, double height, 
                                  double heightVariance, double desiredPrecision, 
                                  double *achievedPrecision, 
                                  csm::WarningList *warnings) const {
    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::imageToGround");
}

csm::EcefLocus MdisNacSensorModel::imageToProximateImagingLocus(const csm::ImageCoord &imagePt, 
                                            const csm::EcefCoord &groundPt, 
                                            double desiredPrecision, 
                                            double *achievedPrecision, 
                                            csm::WarningList *warnings) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::imageToProximateImagingLocus");
}
             
             
csm::EcefLocus MdisNacSensorModel::imageToRemoteImagingLocus(const csm::ImageCoord &imagePt, 
                                         double desiredPrecision, 
                                         double *achievedPrecision, 
                                         csm::WarningList *warnings) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::imageToRemoteImagingLocus");
}


csm::ImageCoord MdisNacSensorModel::getImageStart() const {

    csm::ImageCoord start;
    start.samp = m_startingDetectorSample;
    start.line = m_startingDetectorLine;
    return start;
}

csm::ImageVector MdisNacSensorModel::getImageSize() const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getImageSize");
}

std::pair<csm::ImageCoord, csm::ImageCoord> MdisNacSensorModel::getValidImageRange() const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getValidImageRange");
}

std::pair<double, double> MdisNacSensorModel::getValidHeightRange() const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getValidHeightRange");
}

csm::EcefVector MdisNacSensorModel::getIlluminationDirection(const csm::EcefCoord &groundPt) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getIlluminationDirection");
}

double MdisNacSensorModel::getImageTime(const csm::ImageCoord &imagePt) const {
  
  // check if the image point is in range
  if (imagePt.samp >= m_startingDetectorSample && 
      imagePt.samp <= (m_startingDetectorSample + m_nSamples) &&
      imagePt.line >= m_startingDetectorSample &&
      imagePt.line <= (m_startingDetectorLine + m_nLines)) {
    return m_ephemerisTime;
  }
  else {
    throw csm::Error(csm::Error::BOUNDS,
                     "Image Coordinate out of Bounds",
                     "MdisNacSensorModel::getImageTime");
  }
}

csm::EcefCoord MdisNacSensorModel::getSensorPosition(const csm::ImageCoord &imagePt) const {
  
  // check if the image point is in range
  if (imagePt.samp >= m_startingDetectorSample && 
      imagePt.samp <= (m_startingDetectorSample + m_nSamples) &&
      imagePt.line >= m_startingDetectorSample &&
      imagePt.line <= (m_startingDetectorLine + m_nLines)) {
    csm::EcefCoord sensorPosition;
    sensorPosition.x = m_spacecraftPosition[0];
    sensorPosition.y = m_spacecraftPosition[1];
    sensorPosition.z = m_spacecraftPosition[2];

    return sensorPosition;
  }
  else {
    throw csm::Error(csm::Error::BOUNDS,
                     "Image Coordinate out of Bounds",
                     "MdisNacSensorModel::getSensorPosition");
  }
}

csm::EcefCoord MdisNacSensorModel::getSensorPosition(double time) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getSensorPosition");
}

csm::EcefVector MdisNacSensorModel::getSensorVelocity(const csm::ImageCoord &imagePt) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getSensorVelocity");
}

csm::EcefVector MdisNacSensorModel::getSensorVelocity(double time) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getSensorVelocity");
}

csm::RasterGM::SensorPartials MdisNacSensorModel::computeSensorPartials(int index, const csm::EcefCoord &groundPt, 
                                           double desiredPrecision, 
                                           double *achievedPrecision, 
                                           csm::WarningList *warnings) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::computeSensorPartials");
}

csm::RasterGM::SensorPartials MdisNacSensorModel::computeSensorPartials(int index, const csm::ImageCoord &imagePt, 
                                          const csm::EcefCoord &groundPt, 
                                          double desiredPrecision, 
                                          double *achievedPrecision, 
                                          csm::WarningList *warnings) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::computeSensorPartials");
}
                                              
std::vector<double> MdisNacSensorModel::computeGroundPartials(const csm::EcefCoord &groundPt) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::computeGroundPartials");
}

const csm::CorrelationModel& MdisNacSensorModel::getCorrelationModel() const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getCorrelationModel");
}

std::vector<double> MdisNacSensorModel::getUnmodeledCrossCovariance(const csm::ImageCoord &pt1, 
                                                const csm::ImageCoord &pt2) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getUnmodeledCrossCovariance");
}




csm::Version MdisNacSensorModel::getVersion() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getVersion");
}


std::string MdisNacSensorModel::getModelName() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getModelName");
}


std::string MdisNacSensorModel::getPedigree() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getPedigree");
}


std::string MdisNacSensorModel::getImageIdentifier() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getImageIdentifier");
}


void MdisNacSensorModel::setImageIdentifier(const std::string& imageId,
                                            csm::WarningList* warnings) {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::setImageIdentifier");
}


std::string MdisNacSensorModel::getSensorIdentifier() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getSensorIdentifier");
}


std::string MdisNacSensorModel::getPlatformIdentifier() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getPlatformIdentifier");
}


std::string MdisNacSensorModel::getCollectionIdentifier() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getCollectionIdentifier");
}


std::string MdisNacSensorModel::getTrajectoryIdentifier() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getTrajectoryIdentifier");
}


std::string MdisNacSensorModel::getSensorType() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getSensorType");
}


std::string MdisNacSensorModel::getSensorMode() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getSensorMode");
}


std::string MdisNacSensorModel::getReferenceDateAndTime() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getReferenceDateAndTime");
}


std::string MdisNacSensorModel::getModelState() const {
  // TEMPORARY
  /* commented out for testing the gtest framework
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getModelState");
  */
  return "";
}


void MdisNacSensorModel::replaceModelState(const std::string& argState) {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::replaceModelState");
}




csm::EcefCoord MdisNacSensorModel::getReferencePoint() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getReferencePoint");
}


void MdisNacSensorModel::setReferencePoint(const csm::EcefCoord &groundPt) {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::setReferencePoint");
}


int MdisNacSensorModel::getNumParameters() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getNumParameters");
}


std::string MdisNacSensorModel::getParameterName(int index) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getParameterName");
}


std::string MdisNacSensorModel::getParameterUnits(int index) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getParameterUnits");
}


bool MdisNacSensorModel::hasShareableParameters() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::hasShareableParameters");
}


bool MdisNacSensorModel::isParameterShareable(int index) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::isParameterShareable");
}


csm::SharingCriteria MdisNacSensorModel::getParameterSharingCriteria(int index) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getParameterSharingCriteria");
}


double MdisNacSensorModel::getParameterValue(int index) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getParameterValue");
}


void MdisNacSensorModel::setParameterValue(int index, double value) {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::setParameterValue");
}


csm::param::Type MdisNacSensorModel::getParameterType(int index) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getParameterType");
}


void MdisNacSensorModel::setParameterType(int index, csm::param::Type pType) {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::setParameterType");
}


double MdisNacSensorModel::getParameterCovariance(int index1, int index2) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getParameterCovariance");
}


void MdisNacSensorModel::setParameterCovariance(int index1, int index2, double covariance) {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::setParameterCovariance");
}


int MdisNacSensorModel::getNumGeometricCorrectionSwitches() const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getNumGeometricCorrectionSwitches");
}


std::string MdisNacSensorModel::getGeometricCorrectionName(int index) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getGeometricCorrectionName");
}


void MdisNacSensorModel::setGeometricCorrectionSwitch(int index,
                                                      bool value,
                                                      csm::param::Type pType) {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::setGeometricCorrectionSwitch");
}


bool MdisNacSensorModel::getGeometricCorrectionSwitch(int index) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getGeometricCorrectionSwitch");
}


std::vector<double> MdisNacSensorModel::getCrossCovarianceMatrix(
    const GeometricModel &comparisonModel,
    csm::param::Set pSet,
    const GeometricModelList &otherModels) const {
  throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
                   "Unsupported function",
                   "MdisNacSensorModel::getCrossCovarianceMatrix");
}
