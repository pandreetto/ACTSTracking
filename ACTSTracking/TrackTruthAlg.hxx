#ifndef TrackTruthAlg_h
#define TrackTruthAlg_h 1

// edm4hep
#include <edm4hep/TrackCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/MCRecoTrackParticleAssociationCollection.h>
#include <edm4hep/MCRecoTrackerHitPlaneAssociationCollection.h>

// Gaudi
#include <GaudiAlg/GaudiAlgorithm.h>
#include <GaudiAlg/Transformer.h>
#include <k4FWCore/BaseClass.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>

// std
#include <string>
#include <tuple>
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

struct TrackTruthAlg final : Gaudi::Functional::MultiTransformer<std::tuple<edm4hep::MCRecoTrackParticleAssociationCollection>(
			const edm4hep::TrackCollection &, 
			const edm4hep::MCRecoTrackerHitPlaneAssociationCollection &)> {
public:
	TrackTruthAlg(const std::string& name, ISvcLocator* svcLoc);

	std::tuple<edm4hep::MCRecoTrackParticleAssociationCollection> operator()(
			const edm4hep::TrackCollection& tracks,
                        const edm4hep::MCRecoTrackerHitPlaneAssociationCollection& trackerHitRelations) const;
};

#endif
