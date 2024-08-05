#ifndef ACTSMergeHitCollections_h
#define ACTSMergeHitCollections_h 1

// edm4hep
#include <edm4hep/TrackerHitPlaneCollection.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <GaudiAlg/Transformer.h>

// Standard
#include <tuple>

/**
 * @brief Combine 6 TrackerHitPlaneCollections collections into 1
 * @TODO: I know for a fact that this and ACTSMergeRelationCollections can be the same class. I just cannot get the inhertance structure of edm4hep to work for me...
 * @author Samuel Ferraro
 * @version $Id$
 */
struct ACTSMergeHitCollections final : Gaudi::Functional::MultiTransformer<std::tuple<edm4hep::TrackerHitPlaneCollection>(
		const DataWrapper<edm4hep::TrackerHitPlaneCollection> &,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection> &,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection> &,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection> &,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection> &,
		const DataWrapper<edm4hep::TrackerHitPlaneCollection> &)> {
public:
	/**
         * @brief Constructor for ACTSMergeHitCollections
         * @param name unique string identifier for this instance
         * @param svcLoc a Service Locator passed by the Gaudi AlgManager
         */
	ACTSMergeHitCollections(const std::string& name, ISvcLocator* svcLoc);

	/**
         * @brief ACTSMergeHitCollection operation. The workhorse of this MultiTransformer.
         * @TODO: Note that each param has to be in a DataWrapper. This is because lcio2edm4hep puts things in DataWrappers instead of AnyDataWrappers. If this gets updated, which I'm sure it will, this will have to change.
         * @param col A collection of tracker hits from one section of the detector
         * @return A merged collection with all tracker hits.
         */
	std::tuple<edm4hep::TrackerHitPlaneCollection> operator()(
		const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col1,
                const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col2,
                const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col3,
                const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col4,
                const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col5,
                const DataWrapper<edm4hep::TrackerHitPlaneCollection>& col6) const override;
};
#endif // ACTSMergeHitCollections_h
