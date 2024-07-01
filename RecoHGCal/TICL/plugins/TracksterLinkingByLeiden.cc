#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoHGCal/TICL/plugins/TracksterLinkingbyLeiden.h"

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

#include <cmath>
#include <queue>
#include <cassert>
#include <cmath>
#include <random>
#include <string>
#include <utility>

namespace ticl {

  TracksterLinkingbyLeiden::TracksterLinkingbyLeiden(const edm::ParameterSet &conf, edm::ConsumesCollector iC)
      : TracksterLinkingAlgoBase(conf, iC) {}

  TracksterLinkingbyLeiden::~TracksterLinkingbyLeiden() {}

  void TracksterLinkingbyLeiden::buildLayers() {
    // build disks at HGCal front & EM-Had interface for track propagation

    float zVal = hgcons_->waferZ(1, true);
    std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);

    float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
    std::pair<float, float> rMinMax_interface = hgcons_->rangeR(zVal_interface, true);

    for (int iSide = 0; iSide < 2; ++iSide) {
      float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
      firstDisk_[iSide] = std::make_unique<GeomDet>(
          Disk::build(Disk::PositionType(0, 0, zSide),
                      Disk::RotationType(),
                      SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
              .get());

      zSide = (iSide == 0) ? (-1. * zVal_interface) : zVal_interface;
      interfaceDisk_[iSide] = std::make_unique<GeomDet>(
          Disk::build(Disk::PositionType(0, 0, zSide),
                      Disk::RotationType(),
                      SimpleDiskBounds(rMinMax_interface.first, rMinMax_interface.second, zSide - 0.5, zSide + 0.5))
              .get());
    }
  }

  void TracksterLinkingbyLeiden::initialize(const HGCalDDDConstants *hgcons,
                                            const hgcal::RecHitTools rhtools,
                                            const edm::ESHandle<MagneticField> bfieldH,
                                            const edm::ESHandle<Propagator> propH) {
    hgcons_ = hgcons;
    rhtools_ = rhtools;
    buildLayers();

    bfield_ = bfieldH;
    propagator_ = propH;
  }

  void TracksterLinkingbyLeiden::linkTracksters(
      const Inputs &input,
      std::vector<Trackster> &resultTracksters,
      std::vector<std::vector<unsigned int>> &linkedResultTracksters,
      std::vector<std::vector<unsigned int>> &linkedTracksterIdToInputTracksterId) {
    //creating the graph
    TICLGraph graph{};
    TICLGraphProducer(input, graph);
    //std::cout << "GRAPH NODE SIZE " << graph.getNodes().size() << std::endl;
    //std::cout << "TRACKSTER NUMBER " << input.tracksters.size() << std::endl;

    //applying Leiden algo
    Partition partition{std::vector<Community>{}};
    std::vector<Flat> flatFinalPartition;
    singletonPartition(graph, partition);
    //std::cout << "GRAPH NODE SIZE " << partition.getCommunities().size() << std::endl;

    leidenAlgorithm(graph, partition, flatFinalPartition);
    linkedTracksterIdToInputTracksterId.reserve(flatFinalPartition.size());
    //preparing result output
    for (unsigned int i = 0; i < flatFinalPartition.size(); ++i) {
      const auto &flat = flatFinalPartition[i];
      std::vector<unsigned int> linkedTracksters;
      Trackster outTrackster;
      std::vector<unsigned int> InputTracksterIdList{};
      if (!(flat.empty())) {
        for (auto const &elementary : flat) {
          auto tracksterIndex = elementary.getId();
          InputTracksterIdList.push_back(tracksterIndex);
          outTrackster.mergeTracksters(input.tracksters[tracksterIndex]);
        }
        assert(flat.size() == InputTracksterIdList.size());
        linkedTracksterIdToInputTracksterId.push_back(InputTracksterIdList);
        linkedTracksters.push_back(resultTracksters.size());
        resultTracksters.push_back(outTrackster);
        // Store the linked tracksters
        linkedResultTracksters.push_back(linkedTracksters);
      }
    }
  }

  void TracksterLinkingbyLeiden::fillPSetDescription(edm::ParameterSetDescription &desc) {
    desc.add<std::string>("cutTk",
                          "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                          "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
    TracksterLinkingAlgoBase::fillPSetDescription(desc);
  }

  void TracksterLinkingbyLeiden::leidenAlgorithm(TICLGraph &graph,
                                                 Partition &partition,
                                                 std::vector<Flat> &flatFinalPartition) {
    auto const nEdges = totalEdges(graph);
    bool hasNodeBeenMoved{false};
    moveNodesFast(graph, partition, nEdges, hasNodeBeenMoved);
    //moveNodesFast(partition, gamma_);

    if (!(isAlgorithmDone(graph, partition, hasNodeBeenMoved))) {
      Partition refinedPartition = Partition{std::vector<Community>{}};
      assert(refinedPartition.getCommunities().empty());
      refinePartition(graph, partition, refinedPartition, gamma_, nEdges, theta_);
      //create aggregate graph based on refined partition P_ref
      aggregateGraph(graph, partition, refinedPartition);
      //but maintain partition Pf
      auto &communities = partition.getCommunities();
      std::vector<Community> aggregatedCommunities{};
      for (auto const &community : communities) {
        Community aggregatedCommunity{{}, degree(community) + 1};
        for (auto const &aggregateNode : graph.getNodes()) {
          assert(std::holds_alternative<Community>(aggregateNode));
          if (isCommunityContained(std::get<Community>(aggregateNode), community)) {
            aggregatedCommunity.getNodes().push_back(aggregateNode);
          }
        }
        aggregatedCommunities.push_back(aggregatedCommunity);
      }
      communities = aggregatedCommunities;
      leidenAlgorithm(graph, partition, flatFinalPartition);
    }

    else {
      partition.flattenPartition(flatFinalPartition);
    }
  }

  void TracksterLinkingbyLeiden::TICLGraphProducer(const Inputs &input, TICLGraph &graph) {
    auto const &trackstersclue3d = input.tracksters;

    TICLLayerTile tracksterTilePos;
    TICLLayerTile tracksterTileNeg;

    for (size_t id_t = 0; id_t < trackstersclue3d.size(); ++id_t) {
      auto t = trackstersclue3d[id_t];
      if (t.barycenter().eta() > 0.) {
        tracksterTilePos.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
      } else if (t.barycenter().eta() < 0.) {
        tracksterTileNeg.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
      }
    }

    std::vector<Elementary> allElemNodes;

    for (size_t id_t = 0; id_t < trackstersclue3d.size(); ++id_t) {
      auto t = trackstersclue3d[id_t];

      Elementary tNode(id_t);

      auto bary = t.barycenter();
      double del = 0.1;

      double eta_min = std::max(abs(bary.eta()) - del, (double)TileConstants::minEta);
      double eta_max = std::min(abs(bary.eta()) + del, (double)TileConstants::maxEta);

      if (bary.eta() > 0.) {
        std::array<int, 4> search_box =
            tracksterTilePos.searchBoxEtaPhi(eta_min, eta_max, bary.phi() - del, bary.phi() + del);
        if (search_box[2] > search_box[3]) {
          search_box[3] += TileConstants::nPhiBins;
        }

        for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
          for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
            auto &neighbours = tracksterTilePos[tracksterTilePos.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
            for (auto n : neighbours) {
              if (trackstersclue3d[n].barycenter().z() < bary.z()) {
                tNode.addInnerNeighbour(n);
                tNode.addNeighbour(n);
              } else if (trackstersclue3d[n].barycenter().z() > bary.z()) {
                tNode.addOuterNeighbour(n);
                tNode.addNeighbour(n);
              }
            }
          }
        }
      }

      else if (bary.eta() < 0.) {
        std::array<int, 4> search_box =
            tracksterTileNeg.searchBoxEtaPhi(eta_min, eta_max, bary.phi() - del, bary.phi() + del);
        if (search_box[2] > search_box[3]) {
          search_box[3] += TileConstants::nPhiBins;
        }

        for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
          for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
            auto &neighbours = tracksterTileNeg[tracksterTileNeg.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
            for (auto n : neighbours) {
              if (abs(trackstersclue3d[n].barycenter().z()) < abs(bary.z())) {
                tNode.addInnerNeighbour(n);
                tNode.addNeighbour(n);
              } else if (abs(trackstersclue3d[n].barycenter().z()) > abs(bary.z())) {
                tNode.addOuterNeighbour(n);
                tNode.addNeighbour(n);
              }
            }
          }
        }
      }
      allElemNodes.push_back(tNode);
    }
    std::vector<Node> allNodes{};
    allNodes.reserve(allElemNodes.size());
    for (auto const &e : allElemNodes) {
      Node node{e};
      allNodes.push_back(node);
    }
    graph.setNodes(allNodes);
  }

  bool isAlgorithmDone(TICLGraph const &graph, Partition const &partition, bool hasNodeBeenMoved) {
    //empty communities should have been removed in MoveNodesFast
    auto const &communities = partition.getCommunities();
    for (auto const &community : communities) {
      assert(!(community.getNodes().empty()));
    }

    return !(hasNodeBeenMoved) || partition.getCommunities().size() == graph.getNodes().size();
  }

  //quality function, Constant Potts Model
  //tested on godbolt
  auto CPM(Partition const &partition, int gamma) {
    int CPMResult{};
    for (auto const &community : partition.getCommunities()) {
      int n{communitySize(community, 0)};
      CPMResult += (numberOfEdges(community, community) - gamma * n * (n - 1) / 2);
    }
    return CPMResult;
  }

  //interpreting E(C,C) as non null even if C has a single node, if degree(node)>0
  /*auto CPM_contribution_from_new_community(Node const &node, int gamma) {
  Community newCommunity{std::vector<Node>{node}, degree(node) + 1};
  int n = {communitySize(newCommunity, 0)};
  int result{numberOfEdges(newCommunity, newCommunity) - gamma * n * (n - 1) / 2};
  return result;
}*/

  /*auto CPM_after_move(Partition const &partition,
                    int gamma,
                    Community const &communityFrom,
                    Community const &communityTo,
                    Node const &node) {
  int CPMResult{};
  auto const &communities = partition.getCommunities();
  for (auto const &community : communities) {
    if (community == communityFrom && communityFrom.getNodes().size() > 1) {
      std::vector<Node> vectorWithoutNode{};
      std::copy_if(communityFrom.getNodes().begin(),
                   communityFrom.getNodes().end(),
                   std::back_inserter(vectorWithoutNode),
                   [&](Node const &n) { return !(n == node); });
      Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
      CPMResult += (numberOfEdges(communityWithoutNode, communityWithoutNode) -
                    gamma * (communitySize(communityWithoutNode, 0) * communitySize(communityWithoutNode, 0) / 2));
    } else if (community == communityTo) {
      Community communityWithNewNode{community};
      communityWithNewNode.getNodes().push_back(node);
      CPMResult += (numberOfEdges(communityWithNewNode, communityWithNewNode) -
                    gamma * (communitySize(communityWithNewNode, 0) * communitySize(communityWithNewNode, 0) / 2));

    } else {
      CPMResult += (numberOfEdges(community, community) -
                    gamma * (communitySize(community, 0) * communitySize(community, 0) / 2));
      //PROBLEM CAUSED AFTER THIS LINE
    }
  }
  return CPMResult;
}*/

  auto delta_CPM_after_move(int gamma, Community const &communityFrom, Community const &communityTo, Node const &node) {
    std::vector<Node> vectorWithoutNode{};
    std::copy_if(communityFrom.getNodes().begin(),
                 communityFrom.getNodes().end(),
                 std::back_inserter(vectorWithoutNode),
                 [&](Node const &n) { return !(n == node); });
    Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
    int first_term{numberOfEdges(communityWithoutNode, communityWithoutNode) -
                   numberOfEdges(communityFrom, communityFrom) + communitySize(communityFrom, 0) - 1};
    Community communityWithNewNode{communityTo};
    communityWithNewNode.getNodes().push_back(node);
    int second_term{numberOfEdges(communityWithNewNode, communityWithNewNode) -
                    numberOfEdges(communityTo, communityTo) - communitySize(communityTo, 0)};

    return first_term + second_term;
  }

  double delta_modularity_after_move(int totalEdges,
                                     Community const &communityFrom,
                                     Community const &communityTo,
                                     Node const &node) {
    std::vector<Node> vectorWithoutNode{};
    std::copy_if(communityFrom.getNodes().begin(),
                 communityFrom.getNodes().end(),
                 std::back_inserter(vectorWithoutNode),
                 [&](Node const &n) { return !(n == node); });
    Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
    Community communityWithNewNode{communityTo};
    communityWithNewNode.getNodes().push_back(node);
    auto first_term =
        numberOfEdges(communityWithoutNode, communityWithoutNode) - numberOfEdges(communityFrom, communityFrom) +
        numberOfEdges(communityWithNewNode, communityWithNewNode) - numberOfEdges(communityTo, communityTo);
    auto degreeNode = kappa(node);
    double second_term = degreeNode * degreeNode / (2. * totalEdges);

    auto sumDegreeFrom = kappa(communityFrom);
    auto sumDegreeTo = kappa(communityTo);
    double third_term = degreeNode * (sumDegreeFrom - sumDegreeTo) / (2. * totalEdges);

    return (first_term - second_term + third_term) / totalEdges;
  }

  auto delta_CPM_from_empty(Node const &node, int gamma, Community const &communityFrom) {
    std::vector<Node> vectorWithoutNode{};
    std::copy_if(communityFrom.getNodes().begin(),
                 communityFrom.getNodes().end(),
                 std::back_inserter(vectorWithoutNode),
                 [&](Node const &n) { return !(n == node); });
    Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
    int first_term{numberOfEdges(communityWithoutNode, communityWithoutNode) -
                   numberOfEdges(communityFrom, communityFrom) + communitySize(communityFrom, 0) - 1};

    Community newCommunity{std::vector<Node>{node}, degree(node) + 1};
    return first_term + numberOfEdges(newCommunity, newCommunity);
  }

  double delta_modularity_from_empty(int totalEdges, Community const &communityFrom, Node const &node) {
    std::vector<Node> vectorWithoutNode{};
    std::copy_if(communityFrom.getNodes().begin(),
                 communityFrom.getNodes().end(),
                 std::back_inserter(vectorWithoutNode),
                 [&](Node const &n) { return !(n == node); });
    Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
    Community newCommunity{{node}, degree(node) + 1};
    int first_term{numberOfEdges(communityWithoutNode, communityWithoutNode) -
                   numberOfEdges(communityFrom, communityFrom) + numberOfEdges(newCommunity, newCommunity)};
    auto degreeNode = kappa(node);
    double second_term = degreeNode * degreeNode / (2. * totalEdges);
    int sumDegreeFrom{kappa(communityFrom)};
    double third_term = kappa(node) * sumDegreeFrom / (2. * totalEdges);

    return (first_term - second_term + third_term) / totalEdges;
  }

  void moveNode(Community &communityFrom, Community &communityTo, Node const &node) {
    communityTo.getNodes().push_back(node);
    auto it = std::find(communityFrom.getNodes().begin(), communityFrom.getNodes().end(), node);
    assert(it != communityFrom.getNodes().end());
    assert(&(*it) != &node);
    communityFrom.getNodes().erase(it);
  }

  /*void queueCommunity(Community &community, std::queue<Node> &queue) {
    //elements are added to the queue in random order
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(community.getNodes().begin(), community.getNodes().end(), g);

    for (auto const &node : community.getNodes()) {
      queue.push(node);
    }
  }*/

  /*Partition &removeEmptyCommunities(Partition &partition) {
    auto &communities = partition.getCommunities();
    communities.erase(std::remove_if(communities.begin(), communities.end(), [](Community const &community) {
      return community.getNodes().size() == 0;
    }));

    auto const &communitiesAfterRemoval = partition.getCommunities();
    for (auto const &communityAfterRemoval : communitiesAfterRemoval) {
      assert(communityAfterRemoval.getNodes().size() != 0);
    }

    return partition;
  }*/

  /*int bestCommunityIndex(std::vector<Community> const &communities,
                         int &bestDeltaCPM,
                         Community const &currentCommunity,
                         Node const &currentNode,
                         int gamma) {
    //variation of quality function if I move the node to a different community
    assert(bestDeltaCPM == 0);
    int indexBestCommunity{};
    for (unsigned int i = 0; i < communities.size(); ++i) {
      //PROBLEM IS IN THE LINE BELOW!!!!!
      auto deltaCPM = delta_CPM_after_move(gamma, currentCommunity, communities[i], currentNode);
      if (deltaCPM > bestDeltaCPM) {
        bestDeltaCPM = deltaCPM;
        indexBestCommunity = i;
      }
    }
    return indexBestCommunity;
  }*/

  int bestCommunityIndexModularity(std::vector<Community> const &communities,
                                   double &bestDeltaModularity,
                                   Community const &currentCommunity,
                                   Node const &currentNode,
                                   int totalEdges) {
    //variation of quality function if I move the node to a different community
    assert(bestDeltaModularity == 0.);
    int indexBestCommunity{};
    for (unsigned int i = 0; i < communities.size(); ++i) {
      auto deltaModularity = delta_modularity_after_move(totalEdges, currentCommunity, communities[i], currentNode);
      if ((deltaModularity - bestDeltaModularity) > 0.0001) {
        bestDeltaModularity = deltaModularity;
        indexBestCommunity = i;
      }
    }
    return indexBestCommunity;
  }

  /*template <class ADAPTER>
  struct hack : private ADAPTER {
    static auto &get(ADAPTER &a) { return a.*(&hack::c); }
  };

  template <class ADAPTER>
  const auto &get_container(ADAPTER &a) {
    return hack<ADAPTER>::get(a);
  }*/

  /*Partition &moveNodesFast(Partition &partition, int gamma) {
  //all nodes are added to queue in random order
  std::vector<Node> allNodes{};
  std::for_each(
      partition.getCommunities().begin(), partition.getCommunities().end(), [&allNodes](Community const &community) {
        for (auto const &node : community.getNodes())
          allNodes.push_back(node);
      });
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(allNodes.begin(), allNodes.end(), g);
  std::queue<Node> queue{};
  //nodes are added to the queue in random order
  for (auto const &node : allNodes) {
    queue.push(node);
  }

  //std::cout << "PARTITION SIZE: " << partition.getCommunities().size() << std::endl;
  //std::cout << "TOTAL N OF NODES: " << allNodes.size() << std::endl;

  while (!(queue.empty())) {
    Node currentNode{queue.front()};
    queue.pop();
    auto &queue_members = get_container(queue);
    //std::cout << "NODES IN QUEUE: " << queue_members.size() << std::endl;
    auto &currentCommunity = partition.getCommunities()[partition.findCommunityIndex(currentNode)];
    auto &communities = partition.getCommunities();
    //std::cout << "NODE ID " << getId(currentNode) << "COMMUNITY INDEX " << partition.findCommunityIndex(currentNode)
    // << std::endl;
    //std::cout << "N OF COMMUNITIES IN PARTITION: " << communities.size() << std::endl;

    int bestDeltaCPM{0};
    int indexBestCommunity{bestCommunityIndex(communities, bestDeltaCPM, currentCommunity, currentNode, gamma)};

    //variation of quality function if I move the node to an empty community
    int deltaCPMFromEmpty{0};
    if (currentCommunity.getNodes().size() != 1) {
      deltaCPMFromEmpty = delta_CPM_from_empty(currentNode, gamma, currentCommunity);
    }
    //std::cout << "BEST CPM " << bestDeltaCPM << "CPM FROM NEW " << deltaCPMFromEmpty << std::endl;

    if (deltaCPMFromEmpty > bestDeltaCPM) {
      bestDeltaCPM = -1;
      Community newCommunity{std::vector<Node>{}, degree(currentNode) + 1};
      communities.push_back(newCommunity);
      moveNode(currentCommunity, newCommunity, currentNode);
      // making sure all nbrs of currentNode who are not in bestCommunity will be visited
      std::for_each(communities.begin(), communities.end() - 1, [&](auto const &community) {
        std::for_each(community.getNodes().begin(), community.getNodes().end(), [&](auto const &node) {
          if (areNeighbours(currentNode, node) &&
              std::find(queue_members.begin(), queue_members.end(), node) == queue_members.end()) {
            queue.push(node);
          }
        });
      });
    }

    if (bestDeltaCPM > 0) {
      moveNode(currentCommunity, communities[indexBestCommunity], currentNode);
      // making sure all nbrs of currentNode who are not in bestCommunity will be visited
      for (unsigned int i = 0; i < communities.size(); ++i) {
        if (i != static_cast<unsigned int>(indexBestCommunity)) {
          std::for_each(communities[i].getNodes().begin(), communities[i].getNodes().end(), [&](auto const &node) {
            if (areNeighbours(currentNode, node) &&
                std::find(queue_members.begin(), queue_members.end(), node) == queue_members.end()) {
              queue.push(node);
            }
          });
        }
      }
    }

    if (currentCommunity.getNodes().size() == 0) {
      communities.erase(std::remove(communities.begin(), communities.end(), currentCommunity));
    }
  }

  return partition;
}*/

  Partition &moveNodesFast(TICLGraph const &graph, Partition &partition, int nEdges, bool &hasNodeBeenMoved) {
    assert(hasNodeBeenMoved == false);
    auto const nNodes = graph.getNodes().size();
    //all nodes are added to queue in random order
    std::vector<Node> queue{};
    queue.reserve(nNodes * 4);
    auto &communities = partition.getCommunities();
    std::for_each(communities.begin(), communities.end(), [&](Community const &community) {
      auto &nodes = community.getNodes();
      queue.insert(queue.end(), nodes.begin(), nodes.end());
    });
    assert(nNodes == queue.size());

    std::random_device rd;
    std::default_random_engine g(rd());
    std::shuffle(queue.begin(), queue.end(), g);
    //tells me the aggregation degree
    auto deg = degree(queue[0]);
    std::size_t front_index{0};
    while (front_index != queue.size()) {
      auto currentNode = queue[front_index];
      ++front_index;
      auto currentCommunityIndex = partition.findCommunityIndex(currentNode);
      auto &communities = partition.getCommunities();
      auto &currentCommunity = communities[currentCommunityIndex];
      double bestDeltaModularity{0};
      int indexBestCommunity{
          bestCommunityIndexModularity(communities, bestDeltaModularity, currentCommunity, currentNode, nEdges)};

      //variation of quality function if I move the node to an empty community
      int deltaModularityFromEmpty{0};
      if (currentCommunity.getNodes().size() != 1) {
        deltaModularityFromEmpty = delta_modularity_from_empty(nEdges, currentCommunity, currentNode);
      }

      if (deltaModularityFromEmpty > bestDeltaModularity && deltaModularityFromEmpty > 0.001 * (5 * deg + 1)) {
        Community newCommunity{{}, degree(currentNode) + 1};
        moveNode(currentCommunity, newCommunity, currentNode);
        assert(!(newCommunity.getNodes().empty()));
        communities.push_back(newCommunity);
        hasNodeBeenMoved = true;
        // making sure all nbrs of currentNode who are not in bestCommunity will be visited
        std::for_each(communities.begin(), communities.end() - 1, [&](auto const &community) {
          std::for_each(community.getNodes().begin(), community.getNodes().end(), [&](auto const &node) {
            if (areNeighbours(currentNode, node) &&
                std::find(queue.begin() + front_index, queue.end(), node) == queue.end()) {
              queue.push_back(node);
            }
          });
        });
      } else if (bestDeltaModularity > 0.0001 * (5 * deg + 1)) {
        moveNode(currentCommunity, communities[indexBestCommunity], currentNode);
        hasNodeBeenMoved = true;
        // making sure all nbrs of currentNode who are not in bestCommunity will be visited
        for (unsigned int i = 0; i < communities.size(); ++i) {
          if (i != static_cast<unsigned int>(indexBestCommunity)) {
            std::for_each(communities[i].getNodes().begin(), communities[i].getNodes().end(), [&](auto const &node) {
              if (areNeighbours(currentNode, node) &&
                  std::find(queue.begin() + front_index, queue.end(), node) == queue.end()) {
                queue.push_back(node);
              }
            });
          }
        }
      }

      if (currentCommunity.getNodes().empty()) {
        auto it = communities.begin() + currentCommunityIndex;
        assert(it != communities.end());
        communities.erase(it);
      }
    }
    return partition;
  }

  //fills an empty partition with a singleton partition
  Partition &singletonPartition(TICLGraph const &graph, Partition &singlePartition) {
    assert(singlePartition.getCommunities().empty());
    auto const &nodes = graph.getNodes();
    auto &communities = singlePartition.getCommunities();
    for (auto const &node : nodes) {
      Community singletonCommunity{std::vector<Node>{node}, degree(node) + 1};
      communities.push_back(singletonCommunity);
    }
    assert(!(singlePartition.getCommunities().empty()));
    assert(singlePartition.getCommunities().size() == nodes.size());
    for (auto const &community : singlePartition.getCommunities()) {
      assert(community.getNodes().size() == 1);
    }
    return singlePartition;
  }

  bool isNodeWellConnected(Node const &node, Community const &subset, int gamma) {
    auto nodes = subset.getNodes();
    auto it = std::find(nodes.begin(), nodes.end(), node);
    if (it != nodes.end()) {
      nodes.erase(it);
    }
    Community singletonCommunity{{node}, degree(node) + 1};
    Community subsetWithoutNode{nodes, degree(subset)};
    auto edges{numberOfEdges(singletonCommunity, subsetWithoutNode)};
    assert(edges >= 0);
    auto nodeSize{communitySize(singletonCommunity, 0)};
    auto subsetSize{communitySize(subset, 0)};
    return edges >= (gamma * nodeSize * (subsetSize - nodeSize));
  }

  bool isCommunityWellConnected(Community const &community, Community const &subset, int gamma) {
    Community subsetMinuscommunity{{}, degree(subset)};
    for (auto const &node : subset.getNodes()) {
      if (std::find(community.getNodes().begin(), community.getNodes().end(), node) == community.getNodes().end()) {
        subsetMinuscommunity.getNodes().push_back(node);
      }
    }
    auto edges{numberOfEdges(community, subsetMinuscommunity)};
    assert(edges >= 0);
    auto comSize{communitySize(community, 0)};
    auto subsetSize{communitySize(subset, 0)};
    return (edges >= (gamma * comSize * (subsetSize - comSize)));
  }

  int extractRandomCommunityIndex(std::vector<Community> const &communities,
                                  Partition const &partition,
                                  Node const &node,
                                  Community const &nodeCommunity,
                                  Community const &subset,
                                  int nEdges,
                                  double theta) {
    std::vector<double> deltaModularities{};
    //std::random_device rd;
    //std::default_random_engine gen(rd());
    std::default_random_engine gen{};

    //calculating delta_H for all communities
    for (auto const &community : communities) {
      if (isCommunityContained(community, subset) && isCommunityWellConnected(community, subset, 1)) {
        deltaModularities.push_back((delta_modularity_after_move(nEdges, nodeCommunity, community, node)));
      } else {
        // communities not well connected or not within subset are not considered
        deltaModularities.push_back(-1.);
      }
    }

    //creating the discrete probability function
    std::vector<double> distribution{};
    for (auto const &deltaModularity : deltaModularities) {
      if (deltaModularity < 0.) {
        distribution.push_back(0.);
      } else {
        assert(theta > 0.);
        distribution.push_back(std::exp(deltaModularity / theta));
      }
    }

    int resultIndex = 0;
    long int distrubutionSize = distribution.size();
    if (std::count(distribution.begin(), distribution.end(), 0.) != distrubutionSize) {
      //extracting a random community
      std::discrete_distribution<> d(distribution.begin(), distribution.end());
      //extracts a random index
      resultIndex = d(gen);
      assert(isCommunityContained(communities[resultIndex], subset));
    } else {
      resultIndex = -1;
    }
    return resultIndex;
  }

  Partition &mergeNodesSubset(Partition &partition, Community const &subset, int gamma, int nEdges, double theta) {
    auto &communities = partition.getCommunities();
    //consider only nodes that are well connected within subset S
    for (auto const &node : subset.getNodes()) {
      if (isNodeWellConnected(node, subset, gamma)) {
        int index = partition.findCommunityIndex(node);
        auto &nodeCommunity = communities[index];
        assert(communitySize(nodeCommunity, 0) != 0);
        //consider only nodes that have not yet been merged
        if (communitySize(nodeCommunity, 0) == 1) {
          int communityToIndex{
              extractRandomCommunityIndex(communities, partition, node, nodeCommunity, subset, nEdges, theta)};
          if (communityToIndex >= 0 && communityToIndex != index) {
            auto &communityTo = communities[communityToIndex];
            moveNode(nodeCommunity, communityTo, node);
            assert(isCommunityContained(communityTo, subset));
          }
          //removing empty communities
          if (nodeCommunity.getNodes().empty()) {
            auto it = communities.begin() + index;
            assert(it != communities.end());
            communities.erase(it);
          }
        }
      }
    }
    return partition;
  }

  Partition &refinePartition(TICLGraph const &graph,
                             Partition const &partition,
                             Partition &singlePartition,
                             int gamma,
                             int nEdges,
                             double theta) {
    //fills an empty partition with a singleton partition
    auto &refinedPartition = singletonPartition(graph, singlePartition);
    auto const &communities = partition.getCommunities();
    for (auto const &community : communities) {
      mergeNodesSubset(refinedPartition, community, gamma, nEdges, theta);
    }
    assert(refinedPartition.getCommunities().size() >= partition.getCommunities().size());
    return refinedPartition;
  }
  //communities become nodes in aggregate graph
  TICLGraph &aggregateGraph(TICLGraph &graph, Partition const &partition, Partition const &refinedPartition) {
    auto const &refinedCommunities = refinedPartition.getCommunities();
    std::vector<Node> aggregatedNodes{};
    aggregatedNodes.reserve(refinedCommunities.size());
    std::for_each(refinedCommunities.begin(), refinedCommunities.end(), [&](auto const &refinedCommunity) {
      assert(!(refinedCommunity.getNodes().empty()));
     /* for (auto const &node : refinedCommunity.getNodes()) {
        auto const &nodeOriginalCommunity = partition.findCommunity(node);
        assert(isCommunityContained(refinedCommunity, nodeOriginalCommunity));
      }*/
      Node aggregateNode = refinedCommunity;
      //assert(isCommunityContained(std::get<Community>(aggregateNode), nodeOriginalCommunity));
      aggregatedNodes.push_back(aggregateNode);
    });
    assert(aggregatedNodes.size() == refinedCommunities.size());
    graph.setNodes(aggregatedNodes);
    return graph;
  }
}  // namespace ticl