#include "ACTSMergeHitCollections.hxx"

#include <edm4hep/TrackerHitPlane.h>


DECLARE_COMPONENT(ACTSMergeHitCollections)


ACTSMergeHitCollections::ACTSMergeHitCollections(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc, {
		KeyValue("InputCollection1", "Collection1"),
		KeyValue("InputCollection2", "Collection2"),
		KeyValue("InputCollection3", "Collection3"),
		KeyValue("InputCollection4", "Collection4"),
		KeyValue("InputCollection5", "Collection5"),
		KeyValue("InputCollection6", "Collection6") },
	      { KeyValue("OutputCollection", "MergedCollection") }) {}

std::tuple<edm4hep::TrackerHitPlaneCollection> ACTSMergeHitCollections::operator()(
		const edm4hep::TrackerHitPlaneCollection& col1,
		const edm4hep::TrackerHitPlaneCollection& col2,
		const edm4hep::TrackerHitPlaneCollection& col3,
		const edm4hep::TrackerHitPlaneCollection& col4,
		const edm4hep::TrackerHitPlaneCollection& col5,
		const edm4hep::TrackerHitPlaneCollection& col6) const{
	edm4hep::TrackerHitPlaneCollection mergedCollection;

	for (const auto& col : {&col1, &col2, &col3, &col4, &col5, &col6}) {
		for (const auto& item : *col) { mergedCollection.push_back(item); }
	}

	return std::make_tuple(std::move(mergedCollection));
}
