#include "TrackTruthAlg.hxx"
#include "Helpers.hxx"


//------------------------------------------------------------------------------------------------

DECLARE_COMPONENT(TrackTruthAlg)

TrackTruthAlg::TrackTruthAlg(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm {
	declareProperty("TrackCollection", 
			m_inColTrack = std::string("Tracks"), 
			"Name of reconstructed track input collection.");
	declareProperty("MCParticleCollectionName", 
			m_inColMCP = std::string("MCParticle"), 
			"Name of the MCParticle input collection");
	declareProperty("TrackerHit2SimTrackerHitRelationName", 
			m_inColH2SH, 
			"Name of TrackerHit to SimTrackHit relation collection");
	declareProperty("Particle2TrackRelationName", 
			m_outColMC2T = std::string("Particle2TrackRelationName"), 
			"Map from MC particle to reconstructed track.")
}

StatusCode TrackTruthAlg::initialize() {
	info() << "Initializing TrackTruthAlg" << endmsg;
	return StatusCode::SUCCESS;
}

StatusCode TrackTruthAlg::execute() {
	// Load Collections
	edm4hep::TrackCollection* tracks = nullptr;
	edm4hep::MCParticleCollection mcParticles* = nullptr;
	if (ACTSTracking::getCollection(evtSvc(), m_inColTrack, tracks).isFailure() ||
	    ACTSTracking::getCollection(evtSvc(), m_inColMCP, mcParticles).isFailure()) {
		return StatusCode::FAILURE;
	}

	// Map TrackerHits to SimTrackerHits
	std::vector<std::map<edm4hep::TrackerHit, edm4hep::SimTrackerHit>> trackerHit2SimHits;
	for (const std::string& name : m_inColH2SH) {
		//Retrieve Collection
		edm4hep::TrackerHitPlaneCollection* trackerHitRelations = nullptr;
		if (ACTSTracking::getCollection(evtSvc(), name, trackerHitRelations).isFailure()) {
			return StatusCode::FAILURE;
		}
		//Store Relations in a map
		std::map<edm4hep::TrackerHit, edm4hep::SimTrackerHit> trackerHit2SimHit;
		for (auto& hitRel : *trackerHitRelations) {
			auto trackerHit = hitRel.getFrom<emd4hep::TrackerHit>();
			auto simTrackerHit = hitRel.getTo<emd4hep::SimTrackerHit>();
			trackerHit2SimHit[trackerHit] = simTrackerHit;
		}
		//Push map to vector of maps
		trackerHit2SimHits.push_back(trackerHit2SimHit);
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

	put(outColMC2T, m_outColMC2T);
	return StatusCode::SUCCESS;
}
	
	
StatusCode TrackTruthAlg::finalize() {
	info() << "Finalizing TrackTruthAlg" << endmsg;
	return StatusCode::SUCCESS;
}

