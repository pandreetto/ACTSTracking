#ifndef ACTSDuplicateRemoval_h
#define ACTSDuplicateRemoval_h 1

#include <edm4hep/Track.h>
#include <edm4hep/TrackCollection.h>

#include <GaudiAlg/GaudiAlgorithm.h>

//! \brief Remove track duplicates
/**
 * If tracks share more than 50% of hits, then
 * remove the best one.
 *
 * @author Karol Krizka
 * @version $Id$
 */
class ACTSDuplicateRemoval : public GaudiAlgorithm {
public:
	ACTSDuplicateRemoval(const std::string& name, ISvcLocator* svcLoc);
	virtual ~ACTSDuplicateRemoval();

	virtual StatusCode initialize();
	virtual StatusCode execute();
	virtual StatusCode finalize();

private:
	Gaudi::Property<std::string> m_inputTrackCollection{this, "InputTrackCollectionName", "TruthTracks", "Name of input track collection"};
	Gaudi::Property<std::string> m_outputTrackCollection{this, "OutputTrackCollection", "DedupedTruthTracks", "Name of output track collection"};
};

namespace ACTSTracking {
//! Compare tracks by quality, best first
/**
 * Decides which track is better using the following algorithm
 *  1. Track with higher number of hits is better
 *  2. Track with smaller chi2 is better
 *
 * 2. is only run with 1. is ambigious.
 */
	bool track_duplicate_compare(const edm4hep::Track& trk1, const edm4hep::Track& trk2);
	bool tracks_quality_compart(const edm4hep::Track& trk1, edm4hep::Track& trk2);
	bol tracks_equal(const edm4hep::Track& trk1, edm4hep::Track& trk2);
}  // namespace ACTSTracking

#endif
