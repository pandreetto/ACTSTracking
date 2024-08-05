#ifndef ACTSMergeRelationCollections_h
#define ACTSMergeRelationCollections_h 1

// edm4hep
#include <edm4hep/MCRecoTrackerHitPlaneAssociationCollection.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <GaudiAlg/Transformer.h>

// Standard
#include <tuple>

/**
 * @brief Combine 6 MCRecoTrackerHitPlaneAssociationCollections collections into 1
 * @TODO: I know for a fact that this and ACTSMergeHitCollections can be the same class. I just cannot get the inhertance structure of edm4hep to work for me...
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
	/**
         * @brief Constructor for ACTSMergeRelationCollections
         * @param name unique string identifier for this instance
         * @param svcLoc a Service Locator passed by the Gaudi AlgManager
         */
	ACTSMergeRelationCollections(const std::string& name, ISvcLocator* svcLoc);

	/**
         * @brief ACTSMergeHitCollection operation. The workhorse of this MultiTransformer.
	 * @TODO: Note that each param has to be in a DataWrapper. This is because lcio2edm4hep puts things in DataWrappers instead of AnyDataWrappers. If this gets updated, which I'm sure it will, this will have to change.
         * @param col A collection of tracker hit associations from one section of the detector
         * @return A merged collection with all tracker hits associations.
         */
	std::tuple<edm4hep::MCRecoTrackerHitPlaneAssociationCollection> operator()(
		const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col1,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col2,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col3,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col4,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col5,
                const DataWrapper<edm4hep::MCRecoTrackerHitPlaneAssociationCollection>& col6) const override;
};
#endif // ACTSMergeRelationCollections_h
