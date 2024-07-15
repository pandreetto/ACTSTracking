#ifndef ACTSDuplicateRemoval_h
#define ACTSDuplicateRemoval_h 1

// edm4hep
#include <edm4hep/TrackCollection.h>
#include <edm4hep/Track.h>

// Gaudi
#include <GaudiAlg/GaudiAlgorithm.h>
#include <GaudiAlg/Transformer.h>
#include <k4FWCore/BaseClass.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>

//! \brief Remove track duplicates
/**
 * If tracks share more than 50% of hits, then
 * remove the best one.
 *
 * @author Karol Krizka, Samuel Ferraro
 * @version $Id$
 */
struct ACTSDuplicateRemoval final : Gaudi::Functional::Transformer <edm4hep::TrackCollection(const edm4hep::TrackCollection&)> {
public:
	ACTSDuplicateRemoval(const std::string& name, ISvcLocator* svcLoc);

	edm4hep::TrackCollection operator()(const edm4hep::TrackCollection& trackCollection) const override;
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
	bool tracks_equal(const edm4hep::Track& trk1, edm4hep::Track& trk2);
}  // namespace ACTSTracking

#endif
