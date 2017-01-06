#ifndef MdisNacSensorModel_h
#define MdisNacSensorModel_h

#include <cmath>
#include <iostream>
#include <vector>

#include "csm/RasterGM.h"


class MdisNacSensorModel : public csm::RasterGM {
  //friend class MdisPlugin;
  
  public:
    MdisNacSensorModel();
    ~MdisNacSensorModel();
    
    virtual csm::ImageCoord groundToImage(const csm::EcefCoord &groundPt, 
                                     double desiredPrecision=0.001, 
                                     double *achievedPrecision=NULL, 
                                     csm::WarningList *warnings=NULL) const;
                                     
    virtual csm::ImageCoordCovar groundToImage(const csm::EcefCoordCovar &groundPt, 
                                          double desiredPrecision=0.001, 
                                          double *achievedPrecision=NULL, 
                                          csm::WarningList *warnings=NULL) const;
                                          
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
    virtual csm::EcefCoord imageToGround(const csm::ImageCoord &imagePt, double height, 
                                    double desiredPrecision=0.001, double *achievedPrecision=NULL, 
                                    csm::WarningList *warnings=NULL) const;
                                    
    virtual csm::EcefCoordCovar imageToGround(const csm::ImageCoordCovar &imagePt, double height, 
                                           double heightVariance, double desiredPrecision=0.001, 
                                           double *achievedPrecision=NULL, 
                                           csm::WarningList *warnings=NULL) const;
                                           
    virtual csm::EcefLocus imageToProximateImagingLocus(const csm::ImageCoord &imagePt, 
                                                      const csm::EcefCoord &groundPt, 
                                                      double desiredPrecision=0.001, 
                                                      double *achievedPrecision=NULL, 
                                                      csm::WarningList *warnings=NULL) const;
                                                      
    virtual csm::EcefLocus imageToRemoteImagingLocus(const csm::ImageCoord &imagePt, 
                                                   double desiredPrecision=0.001, 
                                                   double *achievedPrecision=NULL, 
                                                   csm::WarningList *warnings=NULL) const;
 
    virtual csm::ImageCoord getImageStart() const;
 
    virtual csm::ImageVector getImageSize() const;
 
    virtual std::pair<csm::ImageCoord, csm::ImageCoord> getValidImageRange() const;
  
    virtual std::pair<double, double> getValidHeightRange() const;
 
    virtual csm::EcefVector getIlluminationDirection(const csm::EcefCoord &groundPt) const;
 
    virtual double getImageTime(const csm::ImageCoord &imagePt) const;
 
    virtual csm::EcefCoord getSensorPosition(const csm::ImageCoord &imagePt) const;
 
    virtual csm::EcefCoord getSensorPosition(double time) const;
 
    virtual csm::EcefVector getSensorVelocity(const csm::ImageCoord &imagePt) const;
 
    virtual csm::EcefVector getSensorVelocity(double time) const;
 
    virtual csm::RasterGM::SensorPartials computeSensorPartials(int index, 
                                                                const csm::EcefCoord &groundPt, 
                                                                double desiredPrecision=0.001, 
                                                                double *achievedPrecision=NULL, 
                                                                csm::WarningList *warnings=NULL) const;
 
    virtual csm::RasterGM::SensorPartials computeSensorPartials(int index, 
                                                                const csm::ImageCoord &imagePt, 
                                                                const csm::EcefCoord &groundPt, 
                                                                double desiredPrecision=0.001, 
                                                                double *achievedPrecision=NULL, 
                                                                csm::WarningList *warnings=NULL) const;
                                                 
    virtual std::vector<double> computeGroundPartials(const csm::EcefCoord &groundPt) const;
 
    virtual const csm::CorrelationModel &getCorrelationModel() const;
    
    virtual std::vector<double> getUnmodeledCrossCovariance(const csm::ImageCoord &pt1, 
                                                            const csm::ImageCoord &pt2) const;
                                                            
                                                            
  private:
    
    double m_transX[3];
    double m_transY[3];
    double m_majorAxis;
    double m_omega;
    double m_phi;
    double m_kappa;
    double m_focalLength;
    double m_spacecraftPosition[3];
    double m_ccdCenter;
    
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
    csm::EcefCoord intersect(const std::vector<double> &point, const std::vector<double> &direction, double radius) const;
    
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
    std::vector<double> perpendicular(const std::vector<double> &v1, const std::vector<double> &v2) const;
    
    /**
     * Determines the projection of the first vector onto the second vector.
     * 
     * @param v1 First vector to project.
     * @param v2 Second vector that first vector is being projected onto.
     * 
     * @return @b vector<double> Returns the first vector projected onto the second.
     */
    std::vector<double> project(const std::vector<double> &v1, const std::vector<double> &v2) const;
    
    /**
     * Returns the dot product of two vectors.
     */
    double dot(const std::vector<double> &v1, const std::vector<double> &v2) const;
    
    /**
     * Returns the magnitude of a vector.
     */
    double magnitude(const std::vector<double> &v) const;
    
    /**
     * Normalizes the vector (e.g. returns a unit vector).
     */
    std::vector<double> normalize(const std::vector<double> &v) const;
                                                            
  public:
    static const std::string _SENSOR_MODEL_NAME;
                                                            
    
};

#endif