
#include "MdisNacSensorModel.h"

#include "csm/Error.h"

const std::string MdisNacSensorModel::_SENSOR_MODEL_NAME 
                                      = "ISIS_MDISNAC_USGSAstro_1_Linux64_csm30.so";

MdisNacSensorModel::MdisNacSensorModel() {

  
}

MdisNacSensorModel::~MdisNacSensorModel() {
  
}

csm::ImageCoord MdisNacSensorModel::groundToImage(const csm::EcefCoord &groundPt, 
                              double desiredPrecision, 
                              double *achievedPrecision, 
                              csm::WarningList *warnings) const {

  
    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::groundToImage");
}
                                  
csm::ImageCoordCovar MdisNacSensorModel::groundToImage(const csm::EcefCoordCovar &groundPt, 
                                   double desiredPrecision, 
                                   double *achievedPrecision, 
                                   csm::WarningList *warnings) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::groundToImage");
}
                                      
csm::EcefCoord MdisNacSensorModel::imageToGround(const csm::ImageCoord &imagePt, double height, 
                             double desiredPrecision, double *achievedPrecision, 
                             csm::WarningList *warnings) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::imageToGround");
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

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getImageStart");
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

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getImageTime");
}

csm::EcefCoord MdisNacSensorModel::getSensorPosition(const csm::ImageCoord &imagePt) const {

    throw csm::Error(csm::Error::UNSUPPORTED_FUNCTION,
      "Unsupported function",
      "MdisNacSensorModel::getSensorPosition");
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