// CUDA include files
#include <cuda_runtime.h>

// CMSSW include files
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include <stdio.h>

#include <vector>
#include <math.h> // JS_EDIT
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

namespace WeightedMeanFitterCUDA {

  __global__ inline std::pair<GlobalPoint, double> nearestPoint(
      const GlobalPoint& vertex, reco::Track iclus
  );
  __global__ inline TransientVertex weightedMeanOutlierRejection(
      const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points,
      std::vector<reco::TransientTrack> iclus
  );
  __global__ inline TransientVertex weightedMeanOutlierRejectionBeamSpot(
      const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points,
      std::vector<reco::TransientTrack> iclus,
      const reco::BeamSpot& beamSpot
  );
  __global__ inline TransientVertex weightedMeanOutlierRejectionVarianceAsError(
      const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points,
      std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus
  );



}
