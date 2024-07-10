#ifndef ACTSTruthTrackingAlg_h
#define ACTSTruthTrackingAlg_h 1

#include <edm4hep/Track.h>
#include <edm4hep/TrackerHit.h>

#include <Acts/Definitions/Units.hpp>
#include "Acts/EventData/ParticleHypothesis.hpp"

#include "ACTSAlgBase.hxx"

/**
 * This code performs a true pattern recognition by looping over all MC
 * particles and adding all hits associated to them onto a prototrack. This is
 * then fitted and output.
 */
class ACTSTruthTrackingAlg : public ACTSAlgBase {
public:

	ACTSTruthTrackingAlg(const std::string& name, ISvcLocator* svcLoc);
        virtual ~ACTSTruthTrackingAlg();

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

private:
	/** Call to get collections
	*/
	template<typename T>
        T* getCollection(const std::string& collectionName);

protected:
	// Encoder
	std::shared_ptr<UTIL::BitField64> m_encoder;

	// Get the subdetector ID from a hit
	int getSubdetector(const edm4hep::TrackerHit* hit);

	// Get the layer ID from a hit
	int getLayer(const edm4hep::TrackerHit* hit);

	// Remove hits in the same layer of the same subdetector
	void removeHitsSameLayer(const std::vector<edmp4hep::TrackerHit*>&,
	std::vector<edm4hep::TrackerHit*>&);

	// Collection names for (in/out)put
	std::vector<std::string> m_inputTrackerHitCollections;
	std::vector<std::string> m_inputTrackerHitRelationCollections;
	std::string m_inputParticleCollection;
	std::string m_outputTrackCollection;

	// Run and event counters
	uint32_t m_eventNumber;
	uint32_t m_runNumber;

	// Track fit parameters
	double m_initialTrackError_d0 = 20 * Acts::UnitConstants::um;  // Marlin: 1.e3
	double m_initialTrackError_phi =
		1 * Acts::UnitConstants::degree;    // Marlin default: 1.e1
	double m_initialTrackError_relP = 0.01;  // Marlin default: 1.e-2
	double m_initialTrackError_lambda =
		1 * Acts::UnitConstants::degree;  // Marlin (tanlambda) default: 1.e1
	double m_initialTrackError_z0 =
		100 * Acts::UnitConstants::um;  // Marlin default: 1.e3
	double m_initialTrackError_time =
		1 * Acts::UnitConstants::ns;  // No Marlin default
	// double m_maxChi2perHit;

	uint32_t m_fitFails;
};

#endif
