#ifndef ACTSMergeRelationCollections_h
#define ACTSMergeRelationCollections_h 1

// edm4hep
#include <edm4hep/MCRecoTrackerHitPlaneAssociationCollection.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <GaudiAlg/Transformer.h>

#include <tuple>

//! \brief Combine 6 collections into 1
/**
 * @author Samuel Ferraro
 * @version $Id$
 */
struct ACTSMergeRelationCollections final : Gaudi::Functional::MultiTransformer<std::tuple<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>(
		const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection> &,
		const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection> &,
		const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection> &,
		const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection> &,
		const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection> &,
		const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection> &)> {
public:
	ACTSMergeRelationCollections(const std::string& name, ISvcLocator* svcLoc);

	std::tuple<edm4hep::MCRecoTrackerHitPlaneAssociationCollection> operator()(
		const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col1,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col2,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col3,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col4,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col5,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col6) const override;
};
#endif // ACTSMergeRelationCollections_h
