#ifndef TrackTruthAlg_h
#define TrackTruthAlg_h 1

//#include <edm4hep.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackerHit.h>

#include <GaudiAlg/GaudiAlgorithm>
#include <GaudiKernel/ServiceHandle.h>
#include <GaudiKernel/ISvcLocator.h>
#include <GaudiKernel/ToolHandle.h>
#include <k4FWCore/DataHandle.h>

#include <string>
#include <vector>

/**
 * Helper processor that creates LCRelation collections for track to hit
 * associations to be used with LCTuple.
 *
 * @param  TrackCollection                Names of Track input collections
 * @param  Track2HitRelationName          Name of output collection for track to
 * hit relations
 */

class TrackTruthAlg : public GaudiAlgorithm {
public:
	TrackTruthAlg(const std::string& name, ISvcLocator* svcLoc);
	virtual ~TrackTruthAlg();

	virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

protected:
	/** Input collection names.
	*/
	std::string m_inColTrack;
	std::string m_inColMCP;
	std::vector<std::string> m_inColH2SH;

	/** Output collection names.
	*/
	std::string m_outColMC2T;
};

#endif
