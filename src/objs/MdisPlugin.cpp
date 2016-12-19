#include "MdisPlugin.h"

#include <string>

#include <csm/csm.h>
#include <csm/Plugin.h>
#include <csm/Warning.h>

#include "MdisNacSensorModel.h"

const std::string m_pluginName = "MDIS_PLUGIN";
const std::string m_manufacturerName = "USGS";
const std::string m_releaseDate = "TBA";
const csm::Version m_csmVersion = csm::Version(3, 1, 0);
const int m_numModels = 1;

MdisPlugin::MdisPlugin() {
}


MdisPlugin::~MdisPlugin() {
}


std::string MdisPlugin::getPluginName() const {
  return m_pluginName;
}


std::string MdisPlugin::getManufacturer() const {
  return m_manufacturerName;
}


std::string MdisPlugin::getReleaseDate() const {
  return m_releaseDate;
}


csm::Version MdisPlugin::getCsmVersion() const {
  return m_csmVersion;
}


size_t MdisPlugin::getNumModels() const {
  return m_numModels;
}


std::string MdisPlugin::getModelName(size_t modelIndex) const {
  // return m_numModels.at(modelIndex);
  return m_pluginName;
}


std::string MdisPlugin::getModelFamily(size_t modelIndex) const {
  return "?";
}


csm::Version MdisPlugin::getModelVersion(const std::string &modelName) const {
  return csm::Version(3, 0, 1);
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
