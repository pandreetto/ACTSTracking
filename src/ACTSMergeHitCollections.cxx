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
		const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col1,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col2,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col3,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col4,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col5,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col6) const{
	edm4hep::TrackerHitPlaneCollection mergedCollection;

	for (const auto& col : {col1.getData(), col2.getData(), col3.getData(), col4.getData(), col5.getData(), col6.getData()}) {
		for (const auto& item : *col) {
			(void)mergedCollection->create(item.getCellID(), item.getType(), item.getQuality(), item.getTime(), item.getEDep(), item.getEDepError(), item.getU(), item.getV(), item.getDu(), item.getDv(), item.getPosition(), item.getCovMatrix());
		}
	}

	return std::make_tuple(std::move(mergedCollection));
}
