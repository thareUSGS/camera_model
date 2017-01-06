
#include "MdisNacSensorModel.h"

#include <csm/Error.h>

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