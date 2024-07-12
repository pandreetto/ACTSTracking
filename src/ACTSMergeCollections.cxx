#include "ACTSMergeCollections.hxx"

#include <edm4hep/TrackerHitPlane.h>
#include <edm4hep/MCRecoTrackParticleAssociation.h>

DECLARE_COMPONENT_WITH_ID(ACTSMergeCollections<edm4hep::TrackerHitPlaneCollection>, "ACTSMergeTrackerHitPlaneCollections")
DECLARE_COMPONENT_WITH_ID(ACTSMergeCollections<edm4hep::MCRecoTrackParticleAssociationCollection>, "ACTSMergeMCRecoTrackParticleAssociationCollections")


ACTSMergeCollections::ACTSMergeCollections(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc, {
		KeyValue("InputCollection1", "Collection1"),
		KeyValue("InputCollection2", "Collection2"),
		KeyValue("InputCollection3", "Collection3"),
		KeyValue("InputCollection4", "Collection4"),
		KeyValue("InputCollection5", "Collection5"),
		KeyValue("InputCollection6", "Collection6") },
	      { KeyValue("OutputCollection", "MergedCollection") }) {}

template<typename CollectionType>
CollectionType ACTSMergeCollections<CollectionType>::operator()(
		const CollectionType& col1,
		const CollectionType& col2,
		const CollectionType& col3,
		const CollectionType& col4,
		const CollectionType& col5,
		const CollectionType& col6) const{
	CollectionType mergedCollection;

	for (const auto& col : {col1, col2, col3, col4, col5, col6}) {
		for (const auto& item : col) { mergedCollection.push_back(item); }
	}

	return mergedCollection;
}
