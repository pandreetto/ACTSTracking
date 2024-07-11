#include "TrackTruthAlg.hxx"
#include "Helpers.hxx"

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackerHit.h>




//------------------------------------------------------------------------------------------------

DECLARE_COMPONENT(TrackTruthAlg)

TrackTruthAlg::TrackTruthAlg(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc, {
		KeyValue("InputTrackCollecionName", "Tracks"),
		KeyValue("InputMCParticleCollecionName", "MCParticle"),
		KeyValue("InputTrackerHit2SimTrackerHitRelationName", "TrackMCRelation") },
		{ KeyValue("OutputParticle2TrackRelationName", "Particle2TrackRelationName") })	{}

edm4hep::MCRecoParticleAssociationCollection TrackTruthAlg::operator()(
			const edm4hep::TrackCollecion tracks,
                        const edm4hep::MCParticleCollection mcParticles,
                        const edm4hep::TrackerHitPlaneCollecion trackerHitRelations) {
	// Map TrackerHits to SimTrackerHits
	std::map<edm4hep::TrackerHit, edm4hep::SimTrackerHit> trackerHit2SimHit;
	for (auto& hitRel : *trackerHitRelations) {
		auto trackerHit = hitRel.getFrom<emd4hep::TrackerHit>();
		auto simTrackerHit = hitRel.getTo<emd4hep::SimTrackerHit>();
		trackerHit2SimHit[trackerHit] = simTrackerHit;
	}

	// Map best matches MCP to Track
	std::map<edm4hep::MCParticle, edm4hep::Track> mcBestMatchTrack;
	std::map<edm4hep::MCParticle, float> mcBestMatchFrac;

	for (auto& track: *tracks) {
		//Get Track
		std::map<edm4hep::MCParticle, uint32_t> trackHit2Mc;
		for (auto& hit : track.getTrackerHits()) {
			//Search for SimHit
			edm4hep::SimTrackerHit* simHit = nullptr;
			for (auto& hit2simhit : *trackerHit2SimHits) {
				const auto simHitIter = hit2simhit.find(hit);
				if (simHitIter != trackerHit2SimHit.end()) { //Found SimHit
					auto simHit = simHitIter->second;
					break;
				}
			}
			if (simHit.getParticle().isAvailable()) {
				trackHit2Mc[simHit.getParticle()]++; //Increment MC Particle counter
			}
		}

		// Update Best Matches
		for (const auto& [mcParticle, hitCount] : trackHit2Mc) {
			float frac = static_cast<float>(hitCount) / track.trackerHits_size();
			bool better = mcBestMatchTrack.count(mcParticle) == 0 || // no best matches exist
				      mcBestMatchFrac[mcParticle] < frac; // this match is better (more hits on track)
			if (better) {
				mcBestMatchTrack[mcParticle] = track;
				mcBestMatchFrac[mcParticle] = frac;
			}
		}
	}

	// Save the best matches
	edm4hep::MCRecoParticleAssociationCollection outColMC2T =  std::make_unique<edm4hep::MCRecoParticleAssociationCollection>();
	for (const auto& [mcParticle, track] : mcBestMatchTrack) {
		edm4hep::MutableMCRecoParticleAssociation association;
		association.setRec(track);
		association.setSim(mcParticle);
		association.setWeight(mcBestMatchFrac[mcParticle]);
		outColMC2T.push_back(association);
	}

	return outColMC2T;
}
