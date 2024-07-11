#ifndef TrackTruthAlg_h
#define TrackTruthAlg_h 1

// edm4hep
#include <edm4hep/TrackCollecion.h>
#include <edm4hep/MCParticleCollecion.h>
#include <edm4hep/TrackerHitPlaneCollecion.h>
#include <edm4hep/MCRecoParticleAssociationCollecion.h>

// Gaudi
#include <GaudiAlg/GaudiAlgorithm.h>
#include <GaudiAlg/MultiTransformer.h>
#include <k4FWCore/BaseClass.h>

// k4FWCore
#include <k4WFCore/DataHandle.h>

// std
#include <string>
#include <vector>

/**
 * Helper processor that creates LCRelation collections for track to hit
 * associations to be used with LCTuple.
 *
 * @param  TrackCollection                Names of Track input collections
 * @param  Track2HitRelationName          Name of output collection for track to
 * hit relations
 *
 * @author Samuel Ferraro, Unknown
 */

class TrackTruthAlg : public Gaudi::Functional::MultiTransformer<edm4hep::MCRecoParticleAssociationCollecion>(
			const edm4hep::TrackCollecion, 
			const edm4hep::MCParticleCollection, 
			const edm4hep::TrackerHitPlaneCollecion) {
public:
	TrackTruthAlg(const std::string& name, ISvcLocator* svcLoc);

	edm4hep::MCRecoParticleAssociationCollecion operator()(
			const edm4hep::TrackCollecion tracks, 
                        const edm4hep::MCParticleCollection mcParticles,
                        const edm4hep::TrackerHitPlaneCollecion trackerHitRelations);
};

#endif
