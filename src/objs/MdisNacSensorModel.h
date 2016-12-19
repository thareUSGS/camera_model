#ifndef MdisNacSensorModel_h
#define MdisNacSensorModel_h

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
                                                            
    
};

#endif