#ifndef ACTSMergeHitCollections_h
#define ACTSMergeHitCollections_h 1

// edm4hep
#include <edm4hep/TrackerHitPlaneCollection.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <GaudiAlg/Transformer.h>

#include <tuple>

//! \brief Combine 6 collections into 1
/**
 * @author Samuel Ferraro
 * @version $Id$
 */
struct ACTSMergeHitCollections final : Gaudi::Functional::MultiTransformer<std::tuple<edm4hep::TrackerHitPlaneCollection>(
		const edm4hep::TrackerHitPlaneCollection &,
		const edm4hep::TrackerHitPlaneCollection &,
		const edm4hep::TrackerHitPlaneCollection &,
		const edm4hep::TrackerHitPlaneCollection &,
		const edm4hep::TrackerHitPlaneCollection &,
		const edm4hep::TrackerHitPlaneCollection &)> {
public:
	ACTSMergeHitCollections(const std::string& name, ISvcLocator* svcLoc);

	std::tuple<edm4hep::TrackerHitPlaneCollection> operator()(
		const edm4hep::TrackerHitPlaneCollection& col1,
                const edm4hep::TrackerHitPlaneCollection& col2,
                const edm4hep::TrackerHitPlaneCollection& col3,
                const edm4hep::TrackerHitPlaneCollection& col4,
                const edm4hep::TrackerHitPlaneCollection& col5,
                const edm4hep::TrackerHitPlaneCollection& col6) const override;
};
#endif // ACTSMergeHitCollections_h
