#ifndef ACTSMergeRelationCollections_h
#define ACTSMergeRelationCollections_h 1

// edm4hep
#include <edm4hep/MCRecoTrackParticleAssociationCollection.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <GaudiAlg/Transformer.h>

#include <tuple>

//! \brief Combine 6 collections into 1
/**
 * @author Samuel Ferraro
 * @version $Id$
 */
struct ACTSMergeRelationCollections final : Gaudi::Functional::MultiTransformer<std::tuple<edm4hep::MCRecoTrackParticleAssociationCollection>(
		const edm4hep::MCRecoTrackParticleAssociationCollection &,
		const edm4hep::MCRecoTrackParticleAssociationCollection &,
		const edm4hep::MCRecoTrackParticleAssociationCollection &,
		const edm4hep::MCRecoTrackParticleAssociationCollection &,
		const edm4hep::MCRecoTrackParticleAssociationCollection &,
		const edm4hep::MCRecoTrackParticleAssociationCollection &)> {
public:
	ACTSMergeRelationCollections(const std::string& name, ISvcLocator* svcLoc);

	std::tuple<edm4hep::MCRecoTrackParticleAssociationCollection> operator()(
		const edm4hep::MCRecoTrackParticleAssociationCollection& col1,
                const edm4hep::MCRecoTrackParticleAssociationCollection& col2,
                const edm4hep::MCRecoTrackParticleAssociationCollection& col3,
                const edm4hep::MCRecoTrackParticleAssociationCollection& col4,
                const edm4hep::MCRecoTrackParticleAssociationCollection& col5,
                const edm4hep::MCRecoTrackParticleAssociationCollection& col6) const override;
};
#endif // ACTSMergeRelationCollections_h
