#ifndef ACTSAlgBase_h
#define ACTSAlgBase_h 1

// ACTS
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Plugins/TGeo/TGeoDetectorElement.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

// emd4hep
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>

// Gaudi
#include <GaudiAlg/GaudiAlgorithm.h>
#include <GaudiAlg/Transformer.h>
#include <k4FWCore/BaseClass.h>

// k4FWCore
#include <k4WFCore/DataHandle.h>

#include <tuple>
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
class ACTSAlgBase : public Gaudi::Functional::MultiTransformer<std::tuple<
		    emd4hep::TrackCollecion, 
		    edm4hep::TrackCollecion>(
		    const edm4hep::TrackerHitPlaneCollecion), BaseClass_t> {
	using DetectorElementPtr = std::shared_ptr<const Acts::TGeoDetectorElement>;
	using DetectorStore = std::vector<DetectorElementPtr>;

public:
	ACTSAlgBase(const std::string& name, ISvcLocator* svcLoc);
private:
	StatusCode initialize();
	void buildDetector();
	void buildBfield();

protected:
	//! Path to material file
	std::string m_matFile{};

	//! Path to tracker geometry file
	std::string m_tgeoFile{};

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
