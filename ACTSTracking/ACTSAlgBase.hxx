#ifndef ACTSAlgBase_h
#define ACTSAlgBase_h 1

// ACTS
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Plugins/TGeo/TGeoDetectorElement.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

// edm4hep
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>

// Gaudi
#include <Gaudi/Property.h>
#include <GaudiAlg/Transformer.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <k4FWCore/BaseClass.h>

#include <tuple>
#include <string>
#include "GeometryIdMappingTool.hxx"

//! Base processor for ACTS tracking
/**
 * Performs tasks common to all ACTS processors
 *  - loading tracking geometry
 *
 * Assumes that the global TGeoManager with the geometry
 * description is already loaded. For example, via the
 * InitializeDD4hep processor.
 *
 * @author Karol Krizka, Samuel Ferraro
 * @version $Id$
 */
struct ACTSAlgBase : Gaudi::Functional::MultiTransformer<std::tuple<
		    edm4hep::TrackCollection, 
		    edm4hep::TrackCollection>(
		    const edm4hep::TrackerHitPlaneCollection &)> {
	using DetectorElementPtr = std::shared_ptr<const Acts::TGeoDetectorElement>;
	using DetectorStore = std::vector<DetectorElementPtr>;

public:
	ACTSAlgBase(const std::string& name, ISvcLocator* svcLoc);
	StatusCode initialize();
private:
	void buildDetector();
	void buildBfield();

protected:
	//! Path to material file
	Gaudi::Property<std::string> m_matFile{this, "MatFile", std::string(""), "Path to the material description JSON file. Can be empty."};

	//! Path to tracker geometry file
	Gaudi::Property<std::string> m_tgeoFile{this, "TGeoFile", std::string(""), "Path to the tracker geometry file."};

	std::shared_ptr<ACTSTracking::GeometryIdMappingTool> geoIDMappingTool() const;

	const Acts::MagneticFieldContext& magneticFieldContext() const;
	const Acts::GeometryContext& geometryContext() const;
	const Acts::CalibrationContext& calibrationContext() const;

	std::shared_ptr<Acts::MagneticFieldProvider> magneticField() const;
	std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const;

	//! Find surface for hit
	const Acts::Surface* findSurface(const edm4hep::TrackerHit hit) const;

private:
	std::shared_ptr<ACTSTracking::GeometryIdMappingTool> m_geoIDMappingTool;

	Acts::MagneticFieldContext m_magneticFieldContext;
	std::shared_ptr<Acts::MagneticFieldProvider> m_magneticField;

	Acts::GeometryContext m_geometryContext;
	DetectorStore m_detectorStore;
	std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry = nullptr;

	Acts::CalibrationContext m_calibrationContext;
};

#endif
