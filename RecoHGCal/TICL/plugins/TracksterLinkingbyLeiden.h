#ifndef RecoHGCal_TICL_TracksterLinkingbyLeiden_H__
#define RecoHGCal_TICL_TracksterLinkingbyLeiden_H__

#include <memory>
#include <array>
#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "RecoHGCal/TICL/interface/commons.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "RecoHGCal/TICL/plugins/TICLGraph.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

namespace ticl {
  class TracksterLinkingbyLeiden final : public TracksterLinkingAlgoBase {
  public:
    TracksterLinkingbyLeiden(const edm::ParameterSet &conf, edm::ConsumesCollector iC);
    ~TracksterLinkingbyLeiden() override;

    void initialize(const HGCalDDDConstants *hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;

    void linkTracksters(const Inputs &input,
                        std::vector<Trackster> &resultTracksters,
                        std::vector<std::vector<unsigned int>> &linkedResultTracksters,
                        std::vector<std::vector<unsigned int>> &linkedTracksterIdToInputTracksterId) override;

    static void fillPSetDescription(edm::ParameterSetDescription &desc);

  private:
    void buildLayers();

    std::once_flag initializeGeometry_;

    const HGCalDDDConstants *hgcons_;

    std::unique_ptr<GeomDet> firstDisk_[2];
    std::unique_ptr<GeomDet> interfaceDisk_[2];

    hgcal::RecHitTools rhtools_;

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;

    long long int gamma_{
        1};            //resolution parameter of the algortihm. The higher the gamma, the more communities are yielded
    double theta_{1};  //parameter of the refinement step

    void leidenAlgorithm(TICLGraph &graph, Partition &partition, std::vector<Flat> &flatFinalPartition);

    void TICLGraphProducer(const Inputs &input, edm::Event &evt, const edm::EventSetup &es, TICLGraph &graph);
  };
}  // namespace ticl

bool isAlgorithmDone(ticl::TICLGraph const &graph, ticl::Partition const &partition);

ticl::Partition &removeEmptyCommunities(ticl::Partition &partition);

ticl::Partition &refinePartition(ticl::TICLGraph const &graph,
                                 ticl::Partition &partition,
                                 ticl::Partition &singlePartition,
                                 long long int gamma,
                                 double theta);

ticl::Partition &moveNodesFast(ticl::Partition &partition, long long int gamma);

ticl::Partition &singletonPartition(ticl::TICLGraph const &graph, ticl::Partition &singlePartition);

ticl::Partition &mergeNodesSubset(ticl::Partition &partition, ticl::Community const &subset, long long int gamma);

ticl::TICLGraph &aggregateGraph(ticl::TICLGraph &graph, ticl::Partition const &partition);
#endif
