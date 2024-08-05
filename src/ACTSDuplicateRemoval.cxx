#include "ACTSDuplicateRemoval.hxx"

// edm4hep
#include <edm4hep/TrackerHit.h>
#include <edm4hep/MutableTrack.h>

// Standard
#include <algorithm>

// ACTSTracking
#include "Helpers.hxx"

namespace ACTSTracking {
/**
 * @brief Return true if `trk1` and `trk2` share at least 50% of hits.
 * @TODO: I do not like the way this had to be implemented. Find doesn't seem to play nice... Maybe with TrackerHits maybe with podio::RelationRange
 */
inline bool tracks_equal(const edm4hep::Track& trk1, const edm4hep::Track& trk2) {
	const auto& hits1 = trk1.getTrackerHits();
	const auto& hits2 = trk2.getTrackerHits();

	// Number of overlapping hist
	uint32_t hitOlap = 0;
	for (const auto& hit1 : hits1) {
		for (const auto& hit2 : hits2) {
			// This is the part I have an issue with
			if (hit1.getCellID() == hit2.getCellID() && hit1.getType() == hit2.getType() &&
                            hit1.getQuality() == hit2.getQuality() && hit1.getTime() == hit2.getTime() &&
                            hit1.getPosition() == hit2.getPosition()) { hitOlap++; }
		}
	}

	// Smaller track count
	uint32_t size = std::min(hits1.size(), hits2.size());

	return 2 * hitOlap > size;  // half of smaller track belong to larger track
}

/**
 * @brief Return true if `trk1` is of higher quality than `trk2`.
 */
bool track_quality_compare(const edm4hep::Track& trk1, const edm4hep::Track& trk2) {
	// If number of hits are different, then the one with more
	// hits should be chosen first
	if (trk1.trackerHits_size() != trk2.trackerHits_size())
		return trk1.trackerHits_size() > trk2.trackerHits_size();

	// Same number of hits means I want smaller chi2
	return trk1.getChi2() < trk2.getChi2();
}

inline bool track_duplicate_compare(const edm4hep::Track& trk1, const edm4hep::Track& trk2) {
	return trk1.getTrackStates(edm4hep::TrackState::AtIP).tanLambda < trk2.getTrackStates(edm4hep::TrackState::AtIP).tanLambda;
}
}  // namespace ACTSTracking

DECLARE_COMPONENT(ACTSDuplicateRemoval)

ACTSDuplicateRemoval::ACTSDuplicateRemoval(const std::string& name, ISvcLocator* svcLoc) : Transformer(name, svcLoc,
				KeyValue("InputTrackCollectionName", "TruthTracks"),
				KeyValue("OutputTrackCollectionName", "DedupedTruthTracks")) {}

edm4hep::TrackCollection ACTSDuplicateRemoval::operator()(const edm4hep::TrackCollection& trackCollection) const{
	// Make output collection
	edm4hep::TrackCollection outputTracks;

	// Insertion sort input tracks
	std::vector<edm4hep::Track> sortedInput;
	for (const auto& track : trackCollection) {
		auto insertion_point = std::upper_bound(sortedInput.begin(), sortedInput.end(), track, ACTSTracking::track_duplicate_compare);
		sortedInput.insert(insertion_point, track);
	}
	
	int total = 0;
	int dupes = 0;
	int added = 0;
	// Loop through all inputs and search for nearby equals
	// Remove if they are too similar
	std::vector<edm4hep::Track> finalTracks;
	for (const auto& track : sortedInput) {
		total++;
		bool foundAnEqual = false;
		int startIdx = (finalTracks.size() >= 10) ? finalTracks.size() - 10 : 0;
		for (int i = startIdx; i < finalTracks.size(); ++i) {
			auto otherTrack = finalTracks[i];
			if (!ACTSTracking::tracks_equal(track, otherTrack)) continue;
			foundAnEqual = true;
			dupes++;
			if (ACTSTracking::track_quality_compare(track, otherTrack)) {
				finalTracks[i] = track;
				break;
			}
		}
		if (!foundAnEqual) {
			added++;
			finalTracks.push_back(track);
		}
	}

	for (const auto& track : finalTracks) {
		auto newTrack = outputTracks.create();
		ACTSTracking::makeMutableTrack(&track, &newTrack);
	}

	return outputTracks;
}
