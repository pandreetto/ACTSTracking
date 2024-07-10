#ifndef ACTSTruthCKFTrackingAlg_h
#define ACTSTruthCKFTrackingAlg_h 1

#include <edm4hep/TrackerHit.h>

#include <Acts/Definitions/Units.hpp>

#include "ACTSAlgBase.hxx"

/**
 * This code performs a true pattern recognition by looping over all MC
 * particles and adding all hits associated to them onto a prototrack. This is
 * then fitted and output.
 */
class ACTSTruthCKFTrackingAlg : public ACTSAlgBase {
public:
	ACTSTruthCKFTrackingAlg(const std::string& name, ISvcLocator* svcLoc);
	virtual ~ACTSTruthCKFTrackingAlg();

	virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

private:
	/** Call to get collections
	*/
	template<typename T>
        T* getCollection(const std::string& collectionName);

protected:
	// Collection names for (in/out)put
	std::vector<std::string> m_inputTrackerHitCollections;
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

	double m_CKF_chi2CutOff = 15;
	int32_t m_CKF_numMeasurementsCutOff = 10;

	uint32_t m_fitFails;
};

#endif
