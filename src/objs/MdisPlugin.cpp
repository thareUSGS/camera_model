#include "MdisPlugin.h"

#include <string>

#include <csm/csm.h>
#include <csm/Plugin.h>
#include <csm/Warning.h>

#include "MdisNacSensorModel.h"


MdisPlugin::MdisPlugin() {
}


MdisPlugin::~MdisPlugin() {
}


std::string MdisPlugin::getPluginName() const {
  return "UsgsAstroFrameMdisPluginCSM";
}


std::string MdisPlugin::getManufacturer() const {
  return "UsgsAstrogeology";
}


std::string MdisPlugin::getReleaseDate() const {
  return "TBA";
}


csm::Version MdisPlugin::getCsmVersion() const {
  return csm::Version(3, 1, 0);
}


size_t MdisPlugin::getNumModels() const {
  return 1;
}


std::string MdisPlugin::getModelName(size_t modelIndex) const {
  
  return MdisNacSensorModel::_SENSOR_MODEL_NAME;
}


std::string MdisPlugin::getModelFamily(size_t modelIndex) const {
  return "Raster";
}


csm::Version MdisPlugin::getModelVersion(const std::string &modelName) const {

  return csm::Version(1, 0, 0);
}


bool MdisPlugin::canModelBeConstructedFromState(const std::string &modelName,
                                                const std::string &modelState,
                                                csm::WarningList *warnings) const {
  return false;
}


bool MdisPlugin::canModelBeConstructedFromISD(const csm::Isd &imageSupportData,
                                              const std::string &modelName,
                                              csm::WarningList *warnings) const {
  return true;
}


csm::Model *MdisPlugin::constructModelFromState(const std::string&modelState,
                                                csm::WarningList *warnings) const {
  return NULL;
}


csm::Model *MdisPlugin::constructModelFromISD(const csm::Isd &imageSupportData,
                                              const std::string &modelName,
                                              csm::WarningList *warnings) const {
  return NULL;
}


std::string MdisPlugin::getModelNameFromModelState(const std::string &modelState,
                                                   csm::WarningList *warnings) const {
  return "state";
}


bool MdisPlugin::canISDBeConvertedToModelState(const csm::Isd &imageSupportData,
                                               const std::string &modelName,
                                               csm::WarningList *warnings) const {
  return false;
}


std::string MdisPlugin::convertISDToModelState(const csm::Isd &imageSupportData,
                                               const std::string &modelName,
                                               csm::WarningList *warnings) const {
  return "state";
}
