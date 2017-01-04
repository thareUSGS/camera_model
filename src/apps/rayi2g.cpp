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

vector<double> imageToGround(double sample, double line, double height) {
  double focalLength = 549.302734762479;
  
  vector<double> spacecraftPosition(3);
  
  spacecraftPosition[0] = -6683432.694790688; // SENSOR X
  spacecraftPosition[1] = -27264735.61516516; // SENSOR Y
  spacecraftPosition[2] = 8536481.287414147; // SENSOR Z from isd
  
  const double transX[3] = {0.0, 0.014, 0.0};
  const double transY[3] = {0.0, 0.0, 0.014};
  
  // center the sample line
  // imgSamps / 2, imgLines / 2: 1024x1024
  sample -= 512.0; // ISD needs a center sample in CSM coord
  line -= 512.0; // ISD needs a center line in CSM coord (.5 .5 pixel centers)
  
  double focalPlaneX = transX[0] + (transX[1] * sample) + (transX[2] * line);
  double focalPlaneY = transY[0] + (transY[1] * sample) + (transY[2] * line);
  cout << "Focal plane x: " << focalPlaneX << endl;
  cout << "Focal plane y: " << focalPlaneY << endl;

  double majorAxis = 6051800; // mercury ellipsoid

  vector<double> direction(3);
  
  direction[0] = focalPlaneX / 1000 - spacecraftPosition[0];
  direction[1] = focalPlaneY / 1000 - spacecraftPosition[1];
  direction[2] = focalLength / 1000 - spacecraftPosition[2];
  
  cout << "Direction: " << endl;
  cout << direction[0] << endl;
  cout << direction[1] << endl;
  cout << direction[2] << endl;
  cout << endl;
  
  intersect(spacecraftPosition, direction, majorAxis);
}

vector<double> intersect(const vector<double> &point, 
                         const vector<double> &direction, 
                         double radius) {
  
  vector<double> intersectionPoint(3);
  vector<double> newPoint(3), newDirection(3);
  
  newDirection[0] = direction[0] / radius;
  newDirection[1] = direction[1] / radius;
  newDirection[2] = direction[2] / radius;
  
  cout << "New Direction: " << endl;
  cout << newDirection[0] << endl;
  cout << newDirection[1] << endl;
  cout << newDirection[2] << endl;
  cout << endl;
  
  newPoint[0] = point[0] / radius;
  newPoint[1] = point[1] / radius;
  newPoint[2] = point[2] / radius;
  
  cout << "New Point: " << endl;
  cout << newPoint[0] << endl;
  cout << newPoint[1] << endl;
  cout << newPoint[2] << endl;
  cout << endl;
  
  vector<double> perpendicularV = perpendicular(newPoint, newDirection);
  
  vector<double> newPointProj(3);
  newPointProj[0] = newPoint[0] - perpendicularV[0];
  newPointProj[1] = newPoint[1] - perpendicularV[1];
  newPointProj[2] = newPoint[2] - perpendicularV[2];
  
  double newPointMag = magnitude(newPoint);
  double perpendicularMag = magnitude(perpendicularV);
  
  vector<double> newDirectionNorm = normalize(newDirection);
  
  int sign = 1;
  if (newPointMag > 1.0) {
    if (perpendicularMag > 1.0) {
      return intersectionPoint; // empty
    }
    
    if (dot(newPointProj, newDirection) > 0.0) {
      return intersectionPoint; // empty
    }
    
    if (perpendicularMag == 1.0) {
      intersectionPoint[0] = perpendicularV[0] * radius;
      intersectionPoint[1] = perpendicularV[1] * radius;
      intersectionPoint[2] = perpendicularV[2] * radius;
      return intersectionPoint;
    }
    
    sign = -1;
  }
  else if (newPointMag == 1.0) {
    return point;
  }
  else {
    sign = 1;
  }
  
  double scale = 1 - perpendicularMag * perpendicularMag;
  if (scale < 0.0) {
    scale = 0.0;
  }
  scale = sqrt(scale);
  
  intersectionPoint[0] = perpendicularV[0] + sign * scale * newDirectionNorm[0];
  intersectionPoint[1] = perpendicularV[1] + sign * scale * newDirectionNorm[1];
  intersectionPoint[2] = perpendicularV[2] + sign * scale * newDirectionNorm[2];
  
  return intersectionPoint;
}

vector<double> perpendicular(const vector<double> &v1, const vector<double> &v2) {
  
  double max1 = 0;
  double max2 = 0;
  
  for (int i = 0; i < 3; i++) {
    if (abs(v1[i]) > max1) {
      max1 = abs(v1[i]);
    }
    if (abs(v2[i]) > max2) {
      max2 = abs(v2[i]);
    }
  }
  
  vector<double> newV1(3), newV2(3);
  
  newV1[0] = v1[0]/max1;
  newV1[1] = v1[1]/max1;
  newV1[2] = v1[2]/max1;
  
  newV2[0] = v2[0]/max2;
  newV2[1] = v2[1]/max2;
  newV2[2] = v2[2]/max2;
  
  vector<double> parallelV = project(newV1, newV2);
  
  vector<double> perpendicularV (3);
  
  perpendicularV[0] = v1[0] - parallelV[0] * max1;
  perpendicularV[1] = v1[1] - parallelV[1] * max1;
  perpendicularV[2] = v1[2] - parallelV[2] * max1;
  
  return perpendicularV;
}

vector<double> project(const vector<double> &v1, const vector<double> &v2) {
  double v1Dotv2 = dot(v1, v2);
  double v2Dotv2 = dot(v2, v2);
  vector<double> projV (3);
  projV[0] = v1Dotv2/v2Dotv2 * v2[0];
  projV[1] = v1Dotv2/v2Dotv2 * v2[1];
  projV[2] = v1Dotv2/v2Dotv2 * v2[2];
  return projV;
}

double dot(const vector<double> &v1, const vector<double> &v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double magnitude(const vector<double> &v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

vector<double> normalize(const vector<double> &v) {
  double magnitude;
  magnitude = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  vector<double> n (3);
  n[0] = v[0] / magnitude;
  n[1] = v[1] / magnitude;
  n[2] = v[2] / magnitude;
  return n;
}

