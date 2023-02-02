#ifndef fitterCUDA_h
#define fitterCUDA_h
//#include "CUDADataFormats/Track/interface/TrackForPVHeterogeneous.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
//#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerCUDA.h"
#include "RecoVertex/PrimaryVertexProducer/interface/WeightedMeanFitter.h"
//#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"
//#include <cstddef>
//#include <cstdint>
//#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"
//#include <thrust/device_vector.h>
//#include <thrust/host_vector.h>
//#include <thrust/sort.h>
//#include "HeterogeneousCore/CUDAUtilities/interface/radixSort.h"


namespace fitterCUDA {

  struct algo { //JS_EDIT: moved from PrimaryVertexProducer
    VertexFitter<5>* fitter;
    VertexCompatibleWithBeam* vertexSelector;
    std::string label;
    bool useBeamConstraint;
    double minNdof;
  };

  //void wrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, cudaStream_t stream);
  std::vector<TransientVertex> wrapper(
    algo algorithm,
    std::vector<std::vector<reco::TransientTrack> >&& clusters,
    reco::BeamSpot beamSpot,
    VertexState beamVertexState,
    bool f4D,
    bool validBS,
    bool weightFit,
    bool fVerbose );
}

#endif
