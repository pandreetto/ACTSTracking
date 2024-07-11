########################################################
# cmake file for building Marlin ACTSTracking Package
# @author Karol Krizka, Samuel Ferraro
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################

# project name
PROJECT( ACTSTracking VERSION 2.0.0 )

### DEPENDENCIES ############################################################
FIND_PACKAGE( Gaudi REQUIRED )
FIND_PACKAGE( EDM4HEP REQUIRED )
FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDCore )
FIND_PACKAGE( ROOT REQUIRED )
FIND_PACKAGE( Acts REQUIRED COMPONENTS Core PluginJson PluginTGeo )
FIND_PACKAGE( k4FWCore REQUIRED )
INCLUDE(GNUInstallDirs)

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### DATA ###########################################################
ADD_SUBDIRECTORY(data)
ADD_SUBDIRECTORY(example)

### LIBRARY #################################################################

# sources
CONFIGURE_FILE(src/config.h.in config.h)
ADD_LIBRARY(ACTSTracking SHARED
	./src/ACTSAlgBase.cxx
	./src/ACTSSeededCKFTrackingAlg.cxx
	./src/ACTSDuplicateRemoval.cxx
	./src/TrackTruthAlg.cxx
	./src/GeometryIdSelector.cxx
	./src/GeometryIdMappingTool.cxx
	./src/Helpers.cxx
  )
TARGET_INCLUDE_DIRECTORIES(ACTSTracking
  PRIVATE
    ./ACTSTracking
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Gaudi_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(ACTSTracking
  ${ROOT_LIBRARIES} 
  DD4hep::DDCore
  Gaudi::GaudiKernel
  Gaudi::GaudiAlg
  EMD4HEP::edm4hep
  k4FWCore::k4FWCore
  ActsCore ActsPluginTGeo ActsPluginJson)

INSTALL( TARGETS ACTSTracking DESTINATION lib )

