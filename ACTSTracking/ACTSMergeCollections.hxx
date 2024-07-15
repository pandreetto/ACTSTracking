#ifndef ACTSMergeCollections_h
#define ACTSMergeCollections_h 1

// edm4hep
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/MCRecoTrackParticleAssociationCollection.h>

// Gaudi
#include <GaudiAlg/GaudiAlgorithm.h>
#include <GaudiAlg/MultiTransformer.h>
#include <k4FWCore/BaseClass.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>

//! \brief Combine 6 collections into 1
/**
 * @author Samuel Ferraro
 * @version $Id$
 */
template<typename CollectionType>
struct ACTSMergeCollection final: Gaudi::Functional::MultiTransformer <CollectionType(
		const CollectionType,
		const CollectionType,
		const CollectionType,
		const CollectionType,
		const CollectionType,
		const CollectionType)> {
public:
	ACTSMergeCollection(const std::string& name, ISvcLocator* svcLoc);

	CollectionType operator()(
		const CollectionType col1,
		const CollectionType col2,
                const CollectionType col3,
                const CollectionType col4,
                const CollectionType col5,
                const CollectionType col6) const override;
};

extern template class ACTSMergeCollections<edm4hep::TrackerHitPlaneCollection>;
extern template class ACTSMergeCollections<edm4hep::MCRecoParticleAssociationCollection>;

#endif // ACTSMergeCollections_h
