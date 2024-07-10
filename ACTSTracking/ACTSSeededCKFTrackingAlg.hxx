#ifndef ACTSSeededCKFTrackingAlg_h
#define ACTSSeededCKFTrackingAlg_h 1

#include <edm4hep/TrackerHit.h>
#include <edm4hep/TrackerHitCollection.h>

#include <Acts/Definitions/Units.hpp>

#include "ACTSAlgBase.hxx"
#include "GeometryIdSelector.hxx"

/**
 * This code performs a true pattern recognition by looping over all MC
 * particles and adding all hits associated to them onto a prototrack. This is
 * then fitted and output.
 */
class ACTSSeededCKFTrackingAlg : public ACTSAlgBase {
public:
	ACTSSeededCKFTrackingAlg(const std::string& name, ISvcLocator* svcLoc);
	virtual ~ACTSSeededCKFTrackingAlg();

	virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

protected:
	// Collection names for (in/out)put
	std::vector<std::string> m_inputTrackerHitCollections;
	std::string m_outputSeedCollection;
	std::string m_outputTrackCollection;

	// Run settings
	bool m_runCKF = true;
	bool m_propagateBackward = false;

	// Seed finding configuration
	float m_seedFinding_rMax = 150;
	float m_seedFinding_deltaRMin = 5;
	float m_seedFinding_deltaRMax = 80;
	float m_seedFinding_deltaRMinTop = 0;
	float m_seedFinding_deltaRMaxTop = 0;
	float m_seedFinding_deltaRMinBottom = 0;
	float m_seedFinding_deltaRMaxBottom = 0;
	float m_seedFinding_collisionRegion = 75;
	float m_seedFinding_zMax = 600;
	float m_seedFinding_sigmaScattering = 50;
	float m_seedFinding_radLengthPerSeed = 0.1;
	float m_seedFinding_minPt = 500;
	float m_seedFinding_impactMax = 3 * Acts::UnitConstants::mm;

	StringVec m_seedFinding_zBinEdges;
	int m_zTopBinLen = 1;
	int m_zBottomBinLen = 1;
	int m_phiTopBinLen = 1;
	int m_phiBottomBinLen = 1;

	// Track fit parameters
	double m_initialTrackError_pos;
	double m_initialTrackError_phi;
	double m_initialTrackError_relP;
	double m_initialTrackError_lambda;
	double m_initialTrackError_time =
	100 * Acts::UnitConstants::ns;  // No Marlin default

	double m_CKF_chi2CutOff = 15;
	int32_t m_CKF_numMeasurementsCutOff = 10;

	// Seeding configuration
	std::vector<std::string> m_seedingLayers;
	ACTSTracking::GeometryIdSelector m_seedGeometrySelection;

	uint32_t m_fitFails;
};

#endif
