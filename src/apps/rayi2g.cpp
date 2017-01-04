#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

vector<double> imageToGround(double sample, double line, double height);
vector<double> intersect(const vector<double> &point, 
                         const vector<double> &direction, 
                         double radius);
vector<double> perpendicular(const vector<double> &v1, const vector<double> &v2);
vector<double> project(const vector<double> &v1, const vector<double> &v2);
double dot(const vector<double> &v1, const vector<double> &v2);
double magnitude(const vector<double> &v);
vector<double> normalize(const vector<double> &v);

int main() {
  
  vector<double> i2g = imageToGround(500, 500, 1000);
  cout << i2g[0] << endl;
  cout << i2g[1] << endl;
  cout << i2g[2] << endl;
  
  return 0;
}

/**
 * This function determines if a sample, line intersects the target body and if so, where
 * this intersection occurs in body-fixed coordinates.
 * 
 * @param sample Sample of the input image.
 * @param line Line of the input image.
 * @param height ???
 * 
 * @return @b vector<double> Returns the body-fixed X,Y,Z coordinates of the intersection.
 *                           If no intersection, returns a 3-element vector of 0's.
 */
vector<double> imageToGround(double sample, double line, double height) {
  double focalLength = 549.302734762479; // Will be grabbed from ISD.
  
  vector<double> spacecraftPosition(3);
  
  spacecraftPosition[0] = -6683432.694790688; // SENSOR X
  spacecraftPosition[1] = -27264735.61516516; // SENSOR Y
  spacecraftPosition[2] = 8536481.287414147; // SENSOR Z from isd
  
  const double transX[3] = {0.0, 0.014, 0.0}; // Grabbed from ISD.
  const double transY[3] = {0.0, 0.0, 0.014};
  
  // center the sample line
  // imgSamps / 2, imgLines / 2: THIS ASSUMES 1024x1024 image with no offsets
  sample -= 512.0; // ISD needs a center sample in CSM coord
  line -= 512.0; // ISD needs a center line in CSM coord (.5 .5 pixel centers)
  
  // Convert from sample/line to focal plane coordinates
  double focalPlaneX = transX[0] + (transX[1] * sample) + (transX[2] * line);
  double focalPlaneY = transY[0] + (transY[1] * sample) + (transY[2] * line);
  cout << "Focal plane x: " << focalPlaneX << endl;
  cout << "Focal plane y: " << focalPlaneY << endl;

  double majorAxis = 6051800; // mercury ellipsoid (ISD).

  // Determine the look direction of the camera in body-fixed coordinate system
  vector<double> direction(3);
  
  direction[0] = focalPlaneX / 1000 - spacecraftPosition[0];
  direction[1] = focalPlaneY / 1000 - spacecraftPosition[1];
  direction[2] = focalLength / 1000 - spacecraftPosition[2];
  
  cout << "Direction: " << endl;
  cout << direction[0] << endl;
  cout << direction[1] << endl;
  cout << direction[2] << endl;
  cout << endl;
  
  // Try to intersect.
  intersect(spacecraftPosition, direction, majorAxis);
}


/**
 * Given a spacecraft position in body-fixed coordinates, determine if the look direction vector
 * intersects the surface of the target body.
 * 
 * @param point 3-element vector representing the spacecraft's position in body-fixed.
 * @param direction 3-element vector representing the camera's look direction in body-fixed.
 * @param radius Radius for the spherical body.
 * 
 * @return @b vector<double> Returns a 3-element vector representing the intersection in body-fixed.
 *                           If no intersection, the vector will contain 0's.
 * 
 * This code is adapted from the sensorModelRefactor branch's Geometry3D::intersect() method.
 */
vector<double> intersect(const vector<double> &point, 
                         const vector<double> &direction, 
                         double radius) {
  
  vector<double> intersectionPoint(3);
  vector<double> newPoint(3), newDirection(3);
  
  // Transform the look direction vector into unit sphere space.
  newDirection[0] = direction[0] / radius;
  newDirection[1] = direction[1] / radius;
  newDirection[2] = direction[2] / radius;
  
  cout << "New Direction: " << endl;
  cout << newDirection[0] << endl;
  cout << newDirection[1] << endl;
  cout << newDirection[2] << endl;
  cout << endl;
  
  // Transform the spacecraft position vector into unit sphere space.
  newPoint[0] = point[0] / radius;
  newPoint[1] = point[1] / radius;
  newPoint[2] = point[2] / radius;
  
  cout << "New Point: " << endl;
  cout << newPoint[0] << endl;
  cout << newPoint[1] << endl;
  cout << newPoint[2] << endl;
  cout << endl;
  
  // Find this vector:
  // If you extend the look direction infinitely, find the vector that is
  // perpendicular to that look direction and points to the body's origin.
  vector<double> perpendicularV = perpendicular(newPoint, newDirection);
  
  // The look direction vector extended to the perpendicular vector.
  vector<double> newPointProj(3);
  newPointProj[0] = newPoint[0] - perpendicularV[0];
  newPointProj[1] = newPoint[1] - perpendicularV[1];
  newPointProj[2] = newPoint[2] - perpendicularV[2];
  
  // Find magnitudes of the "unit" spacecraft position vector and perpendicular vector.
  double newPointMag = magnitude(newPoint);
  double perpendicularMag = magnitude(perpendicularV);
  
  // Use max-norm (infinity-norm) to normalize the "unit" look direction vector.
  vector<double> newDirectionNorm = normalize(newDirection);
  
  // Positive sign indicates spacecraft is in the target body,
  // negative indicates spacecraft is outside the target body.
  int sign = 1;
  
  // If the spacecraft position is outside the target body
  if (newPointMag > 1.0) {
    // If vector perpendicular to look direction has magnitude > 1,
    // we are not looking at the target body (since the body is a "unit sphere")
    if (perpendicularMag > 1.0) {
      return intersectionPoint; // empty
    }
    
    // If looking away from the target body, there is no intersection.
    if (dot(newPointProj, newDirection) > 0.0) {
      return intersectionPoint; // empty
    }
    
    // If intersection point is on the limb, then transform back to target body size.
    if (perpendicularMag == 1.0) {
      intersectionPoint[0] = perpendicularV[0] * radius;
      intersectionPoint[1] = perpendicularV[1] * radius;
      intersectionPoint[2] = perpendicularV[2] * radius;
      return intersectionPoint;
    }
    
    sign = -1;
  }
  // If the spacecraft is on the target body
  else if (newPointMag == 1.0) {
    return point;
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
  intersectionPoint[0] = perpendicularV[0] + sign * scale * newDirectionNorm[0];
  intersectionPoint[1] = perpendicularV[1] + sign * scale * newDirectionNorm[1];
  intersectionPoint[2] = perpendicularV[2] + sign * scale * newDirectionNorm[2];
  
  return intersectionPoint;
}


/**
 * Projects the first vector onto the second vector, then finds the vector that is
 * perpendicular to the projected vector.
 * 
 * This function is adapted from the sensorModelRefactor branch's
 * LinearAlgebra::perpendicular() method (which is adapted from naif's vperp).
 *
 * @param v1 Vector to project.
 * @param v2 Vector to project on to.
 * 
 * @return @b vector<double> Returns a vector perpendicular to the first vector
 *                           projected on the second vector.
 */
vector<double> perpendicular(const vector<double> &v1, const vector<double> &v2) {
  
  double max1 = 0;
  double max2 = 0;
  
  // Get the infinite norm for each vector.
  for (int i = 0; i < 3; i++) {
    if (abs(v1[i]) > max1) {
      max1 = abs(v1[i]);
    }
    if (abs(v2[i]) > max2) {
      max2 = abs(v2[i]);
    }
  }
  
  // Scale the vectors by their max norms (not sure why this is needed, optimization maybe?).
  vector<double> newV1(3), newV2(3);
  
  newV1[0] = v1[0] / max1;
  newV1[1] = v1[1] / max1;
  newV1[2] = v1[2] / max1;
  
  newV2[0] = v2[0] / max2;
  newV2[1] = v2[1] / max2;
  newV2[2] = v2[2] / max2;
  
  // Project first vector onto second vector.
  vector<double> parallelV = project(newV1, newV2);
  
  vector<double> perpendicularV(3);
  
  // Get the perpendicular vector by subtracting the first vector by its projection
  // onto the second.
  perpendicularV[0] = v1[0] - parallelV[0] * max1;
  perpendicularV[1] = v1[1] - parallelV[1] * max1;
  perpendicularV[2] = v1[2] - parallelV[2] * max1;
  
  return perpendicularV;
}


/**
 * Determines the projection of the first vector onto the second vector.
 * 
 * @param v1 First vector to project.
 * @param v2 Second vector that first vector is being projected onto.
 * 
 * @return @b vector<double> Returns the first vector projected onto the second.
 */
vector<double> project(const vector<double> &v1, const vector<double> &v2) {
  double v1Dotv2 = dot(v1, v2);
  double v2Dotv2 = dot(v2, v2);
  vector<double> projV(3);
  projV[0] = v1Dotv2 / v2Dotv2 * v2[0];
  projV[1] = v1Dotv2 / v2Dotv2 * v2[1];
  projV[2] = v1Dotv2 / v2Dotv2 * v2[2];
  return projV;
}


/**
 * Returns the dot product of two vectors.
 */
double dot(const vector<double> &v1, const vector<double> &v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


/**
 * Returns the magnitude of a vector.
 */
double magnitude(const vector<double> &v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


/**
 * Normalizes the vector (e.g. returns a unit vector).
 */
vector<double> normalize(const vector<double> &v) {
  double magnitude;
  magnitude = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  vector<double> n (3);
  n[0] = v[0] / magnitude;
  n[1] = v[1] / magnitude;
  n[2] = v[2] / magnitude;
  return n;
}

