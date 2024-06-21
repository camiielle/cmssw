#include <cmath>
#include <queue>
#include <cassert>
#include <cmath>
#include <random>
#include <string>

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

using namespace ticl;

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
    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
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

  //applying Leiden algo
  Partition partition{std::vector<Community>{}};
  std::vector<ticl::Flat> flatFinalPartition;
  singletonPartition(graph, partition);
  leidenAlgorithm(graph, partition, flatFinalPartition);

  //preparing result output
  for (unsigned int i = 0; i < flatFinalPartition.size(); ++i) {
    const auto& flat = flatFinalPartition[i];
    std::vector<unsigned int> linkedTracksters;
    Trackster outTrackster;
    if (!(flat.empty())) {
      for (auto const &elementary : flat) {
        auto tracksterIndex = elementary.getId();
        linkedTracksterIdToInputTracksterId[i].push_back(tracksterIndex);
        outTrackster.mergeTracksters(input.tracksters[tracksterIndex]);
      }
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

void TracksterLinkingbyLeiden::leidenAlgorithm(ticl::TICLGraph &graph,
                                               ticl::Partition &partition,
                                               std::vector<Flat> &flatFinalPartition) {
  moveNodesFast(partition, gamma_);

  if (!(isAlgorithmDone(graph, partition))) {
    Partition refinedPartition = Partition{std::vector<Community>{}};
    assert(refinedPartition.getCommunities().empty());
    refinePartition(graph, partition, refinedPartition, gamma_, theta_);
    //create aggregate graph based on refined partition P_ref
    aggregateGraph(graph, refinedPartition);
    //but maintain partition P
    auto &communities = partition.getCommunities();
    std::vector<Community> aggregatedCommunities{};
    for (auto const &community : communities) {
      Community aggregatedCommunity{std::vector<Node>{}, degree(community) + 1};
      for (auto const &aggregateNode : graph.getNodes()) {
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
  graph.getNodes() = allNodes;
}

bool isAlgorithmDone(TICLGraph const &graph, Partition const &partition) {
  return partition.getCommunities().size() == graph.getNodes().size();
}

//tested on godbolt
long long int factorial(int n) { return (n == 1 || n == 0) ? 1 : n * factorial(n - 1); }

//tested on godbolt
long long int binomialCoefficient(int n, int k) {
  assert(n >= 0);
  assert(k >= 0);
  if (n < k)
    return 0;
  else if (n == k)
    return 1;
  else
    return factorial(n) / (factorial(k) * factorial(n - k));
}

//quality function, Constant Potts Model
//tested on godbolt
long long int CPM(Partition const &partition, long long int gamma) {
  long long int CPMResult{};
  for (auto const &community : partition.getCommunities()) {
    CPMResult += (static_cast<long long int>(numberOfEdges(community, community)) -
                  (gamma)*binomialCoefficient(communitySize(community), 2));
  }
  return CPMResult;
}

//interpreting E(C,C) as non null even if C has a single node, if degree(node)>0
long long int CPM_contribution_from_new_community(Node const &node, long long int gamma) {
  Community newCommunity{std::vector<Node>{node}, degree(node) + 1};
  long long int result{(static_cast<long long int>(numberOfEdges(newCommunity, newCommunity)) -
                        gamma * binomialCoefficient(communitySize(newCommunity), 2))};
  return result;
}

long long int CPM_after_move(Partition const &partition,
                             long long int gamma,
                             Community const &communityFrom,
                             Community const &communityTo,
                             Node const &node) {
  long long int CPMResult{};
  auto const &communities = partition.getCommunities();
  for (auto const &community : communities) {
    if (community == communityFrom) {
      std::vector<Node> vectorWithoutNode{};
      std::copy_if(communityFrom.getNodes().begin(),
                   communityFrom.getNodes().end(),
                   std::back_inserter(vectorWithoutNode),
                   [&](Node const &n) { return !(n == node); });
      Community communityWithoutNode{vectorWithoutNode, communityFrom.getDegree()};
      CPMResult += (static_cast<long long int>(numberOfEdges(communityWithoutNode, communityWithoutNode)) -
                    gamma * binomialCoefficient(communitySize(communityWithoutNode), 2));
    } else if (community == communityTo) {
      Community communityWithNewNode{community};
      communityWithNewNode.getNodes().push_back(node);
      CPMResult += (static_cast<long long int>(numberOfEdges(communityWithNewNode, communityWithNewNode)) -
                    gamma * binomialCoefficient(communitySize(communityWithNewNode), 2));
    } else {
      CPMResult += (numberOfEdges(community, community) - gamma * binomialCoefficient(communitySize(community), 2));
    }
  }
  return CPMResult;
}

void moveNode(Community &communityFrom, Community &communityTo, Node const &node) {
  communityFrom.getNodes().erase(std::remove(communityFrom.getNodes().begin(), communityFrom.getNodes().end(), node));
  communityTo.getNodes().push_back(node);
}

auto queueCommunity(Community &community, std::queue<Node> &queue) {
  //elements are added to the queue in random order
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(community.getNodes().begin(), community.getNodes().end(), g);

  for (auto const &node : community.getNodes()) {
    queue.push(node);
  }
  return queue;
}

Partition &removeEmptyCommunities(Partition &partition) {
  auto &communities = partition.getCommunities();
  communities.erase(std::remove_if(communities.begin(), communities.end(), [](Community const &community) {
    return community.getNodes().size() == 0;
  }));

  auto const &communitiesAfterRemoval = partition.getCommunities();
  for (auto const &communityAfterRemoval : communitiesAfterRemoval) {
    assert(communityAfterRemoval.getNodes().size() != 0);
  }

  return partition;
}

int bestCommunityIndex(Partition const &partition,
                       std::vector<Community> const &communities,
                       long long int &bestDeltaCPM,
                       Community const &currentCommunity,
                       Node const &currentNode,
                       long long int currentCPM,
                       long long int gamma) {
  //variation of quality function if I move the node to a different community
  assert(bestDeltaCPM == 0);
  int indexBestCommunity{};
  int iterationIndex{-1};
  for (auto const &community : communities) {
    ++iterationIndex;
    auto AfterMoveCPM = CPM_after_move(partition, gamma, currentCommunity, community, currentNode);
    auto deltaCPM = AfterMoveCPM - currentCPM;
    if (deltaCPM > bestDeltaCPM) {
      bestDeltaCPM = deltaCPM;
      indexBestCommunity = iterationIndex;
    }
  }
  return indexBestCommunity;
}

Partition &moveNodesFast(Partition &partition, long long int gamma) {
  //all nodes are added to queue in random order
  auto shuffledCommunities = partition.getCommunities();
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(shuffledCommunities.begin(), shuffledCommunities.end(), g);
  std::queue<Node> queue{};
  for (auto &community : shuffledCommunities) {
    queueCommunity(community, queue);
  }

  while (!queue.empty()) {
    Node const &currentNode{queue.front()};
    auto currentCPM = CPM(partition, gamma);
    auto &currentCommunity = partition.getCommunities()[partition.findCommunityIndex(currentNode)];
    auto &communities = partition.getCommunities();

    long long int bestDeltaCPM{0};
    int indexBestCommunity{
        bestCommunityIndex(partition, communities, bestDeltaCPM, currentCommunity, currentNode, currentCPM, gamma)};

    //variation of quality function if I move the node to an empty community
    auto deltaCPMFromEmpty = CPM_contribution_from_new_community(currentNode, gamma) - currentCPM;

    // contains nbrs of currentNode that are not in bestCommunity
    std::vector<Node> currentNeighbours{};

    if (deltaCPMFromEmpty > bestDeltaCPM) {
      bestDeltaCPM = -1;
      Community newCommunity{std::vector<Node>{}, degree(currentNode) + 1};
      communities.push_back(newCommunity);
      moveNode(currentCommunity, newCommunity, currentNode);
      std::for_each(communities.begin(), communities.end() - 1, [&](auto const &community) {
        std::for_each(community.getNodes().begin(), community.getNodes().end(), [&](auto const &node) {
          if (areNeighbours(currentNode, node)) {
            currentNeighbours.push_back(node);
          }
        });
      });
    }

    if (bestDeltaCPM > 0) {
      moveNode(currentCommunity, communities[indexBestCommunity], currentNode);
      std::for_each(communities.begin(), communities.end(), [&](auto const &community) {
        if (!(community == communities[indexBestCommunity])) {
          std::for_each(community.getNodes().begin(), community.getNodes().end(), [&](auto const &node) {
            if (areNeighbours(currentNode, node)) {
              currentNeighbours.push_back(node);
            }
          });
        }
      });
    }

    // making sure all nbrs of currentNode who are not in bestCommunity will be visited
    for (auto const &neighbour : currentNeighbours) {
      queue.push(neighbour);
    }
  }

  //remove communities that, after node moving, are empty and then returns the result
  return removeEmptyCommunities(partition);
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

  return singlePartition;
}

bool isNodeWellConnected(Node const &node, Community const &subset, long long int gamma) {
  auto nodes = subset.getNodes();
  if (std::find(nodes.begin(), nodes.end(), node) != nodes.end()) {
    nodes.erase(std::remove(nodes.begin(), nodes.end(), node));
  }
  Community singletonCommunity{std::vector<Node>{node}, degree(node) + 1};
  Community subsetWithoutNode{nodes, degree(subset)};
  long long int edges{numberOfEdges(singletonCommunity, subsetWithoutNode)};
  assert(edges >= 0);
  long long int nodeSize{communitySize(singletonCommunity)};
  long long int subsetSize{communitySize(subset)};
  return (edges >= (gamma * nodeSize * (subsetSize - nodeSize)));
}

bool isCommunityWellConnected(Community const &community, Community const &subset, long long int gamma) {
  Community subsetMinuscommunity{std::vector<Node>{}, degree(subset)};
  for (auto const &node : subset.getNodes()) {
    if (std::find(community.getNodes().begin(), community.getNodes().end(), node) == community.getNodes().end()) {
      subsetMinuscommunity.getNodes().push_back(node);
    }
  }
  long long int edges{numberOfEdges(community, subsetMinuscommunity)};
  assert(edges >= 0);
  long long int comSize{communitySize(community)};
  long long int subsetSize{communitySize(subset)};
  return (edges >= (gamma * comSize * (subsetSize - comSize)));
}

int extractRandomCommunityIndex(std::vector<Community> const &communities,
                                Partition const &partition,
                                Node const &node,
                                Community const &nodeCommunity,
                                Community const &subset,
                                long long int gamma,
                                double theta) {
  auto currentCPM = CPM(partition, gamma);
  std::vector<long long int> deltaCPMs{};

  //calculating delta_H for all communities
  for (auto const &community : communities) {
    if (isCommunityWellConnected(community, subset, gamma)) {
      deltaCPMs.push_back((CPM_after_move(partition, gamma, nodeCommunity, community, node) - currentCPM));
    } else {
      // communities not well connected are not considered
      deltaCPMs.push_back(-1);
    }
  }

  //creating the discrete probability function
  std::vector<double> distribution{};
  for (auto const &deltaCPM : deltaCPMs) {
    if (deltaCPM < 0) {
      distribution.push_back(0.);
    } else {
      assert(theta > 0);
      distribution.push_back(std::exp(deltaCPM / theta));
    }
  }

  //extracting a random community
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(distribution.begin(), distribution.end());
  //extracts a random index
  int resultIndex = d(gen);

  return resultIndex;
}

Partition &mergeNodesSubset(Partition &partition, Community const &subset, long long int gamma, double theta) {
  auto &communities = partition.getCommunities();

  //consider only nodes that are well connected within subset S
  for (auto const &node : subset.getNodes()) {
    if (isNodeWellConnected(node, subset, gamma)) {
      int index{static_cast<int>(partition.findCommunityIndex(node))};
      auto &nodeCommunity = communities[index];
      assert((communitySize(nodeCommunity)) != 0);
      //consider only nodes that have not yet been merged
      if (communitySize(nodeCommunity) == 1) {
        int communityToIndex{
            extractRandomCommunityIndex(communities, partition, node, nodeCommunity, subset, gamma, theta)};
        auto &communityTo = communities[communityToIndex];
        moveNode(nodeCommunity, communityTo, node);
      }
    }
  }
  return partition;
}

Partition &refinePartition(
    TICLGraph const &graph, Partition &partition, Partition &singlePartition, long long int gamma, double theta) {
  //fills an empty partition with a singleton partition
  auto &refinedPartition = singletonPartition(graph, singlePartition);
  auto const &communities = partition.getCommunities();
  for (auto const &community : communities) {
    mergeNodesSubset(refinedPartition, community, gamma, theta);
  }
  return refinedPartition;
}

//communities become nodes in aggregate graph
TICLGraph &aggregateGraph(TICLGraph &graph, Partition const &partition) {
  std::vector<Community> const &communities{partition.getCommunities()};
  std::vector<Node> aggregatedNodes{};
  aggregatedNodes.reserve(communities.size());
  std::for_each(communities.begin(), communities.end(), [&aggregatedNodes](auto const &community) {
    aggregatedNodes.push_back(Node{community});
  });

  assert(aggregatedNodes.size() == communities.size());
  auto &oldNodes = graph.getNodes();
  oldNodes = aggregatedNodes;

  return graph;
}