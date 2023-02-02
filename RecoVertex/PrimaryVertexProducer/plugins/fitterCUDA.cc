// CUDA include files
#include <cuda_runtime.h>

// CMSSW include files
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "RecoVertex/PrimaryVertexProducer/interface/fitterCUDA.h"
#include "CUDADataFormats/Track/interface/TrackForPVHeterogeneous.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
//#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerCUDA.h" //JS
//#include "RecoVertex/PrimaryVertexProducer/interface/WeightedMeanFitter.h" //JS
//#include "DataFormats/Math/interface/Error.h"
//#include <cstddef>
//#include <cstdint>
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"
//#include <thrust/device_vector.h>
//#include <thrust/host_vector.h>
//#include <thrust/sort.h>
//#include "HeterogeneousCore/CUDAUtilities/interface/radixSort.h"

namespace fitterCUDA {

/*
__global__ void fitterKernel(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks){

}
*/


//#ifdef __CUDACC__
std::vector<TransientVertex> wrapper(
    algo algorithm,
    std::vector<std::vector<reco::TransientTrack> >&& clusters,
    reco::BeamSpot beamSpot,
    VertexState beamVertexState,
    bool f4D,
    bool validBS,
    bool weightFit,
    bool fVerbose ){
    /*
    unsigned int blockSize = 512;
    unsigned int gridSize  = 1;
    fitterKernel<<<gridSize, blockSize, 0, stream>>>(ntracks, tracks);
    cudaCheck(cudaGetLastError());
    */


    std::vector<TransientVertex> pvs;

    for (std::vector<std::vector<reco::TransientTrack> >::const_iterator iclus = clusters.begin();
         iclus != clusters.end();
         iclus++) {
      double sumwt = 0.;
      double sumwt2 = 0.;
      double sumw = 0.;
      double meantime = 0.;
      double vartime = 0.;
      if (f4D) {
        for (const auto& tk : *iclus) {
          const double time = tk.timeExt();
          const double err = tk.dtErrorExt();
          const double inverr = err > 0. ? 1.0 / err : 0.;
          const double w = inverr * inverr;
          sumwt += w * time;
          sumwt2 += w * time * time;
          sumw += w;
        }
        meantime = sumwt / sumw;
        double sumsq = sumwt2 - sumwt * sumwt / sumw;
        double chisq = iclus->size() > 1 ? sumsq / double(iclus->size() - 1) : sumsq / double(iclus->size());
        vartime = chisq / sumw;
      }

      TransientVertex v;
      if (!weightFit) {
        if (algorithm.useBeamConstraint && validBS && (iclus->size() > 1)) {
          v = algorithm.fitter->vertex(*iclus, beamSpot);
        } else if (!(algorithm.useBeamConstraint) && (iclus->size() > 1)) {
          v = algorithm.fitter->vertex(*iclus);
        }
      } else if (weightFit) {
        std::vector<std::pair<GlobalPoint, GlobalPoint>> points;
        if (algorithm.useBeamConstraint && validBS && (iclus->size() > 1)) {
            for (const auto& itrack : *iclus){
                   GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
                   GlobalPoint err(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.track().dzError());
                   std::pair<GlobalPoint, GlobalPoint> p2(p, err);
                   points.push_back(p2);
            }
            v = WeightedMeanFitter::weightedMeanOutlierRejectionBeamSpot(points, *iclus, beamSpot);
            if ((v.positionError().matrix())(2,2) != (WeightedMeanFitter::startError*WeightedMeanFitter::startError)) pvs.push_back(v);

        }
        else if (!(algorithm.useBeamConstraint) && (iclus->size() > 1)) {
           for (const auto& itrack : *iclus){
                   GlobalPoint p = itrack.impactPointState().globalPosition();
                   GlobalPoint err(itrack.track().dxyError(), itrack.track().dxyError(), itrack.track().dzError());
                   std::pair<GlobalPoint, GlobalPoint> p2(p, err);
                   points.push_back(p2);
           }
           v = WeightedMeanFitter::weightedMeanOutlierRejection(points, *iclus);
           if ((v.positionError().matrix())(2,2) != (WeightedMeanFitter::startError*WeightedMeanFitter::startError)) pvs.push_back(v); //FIX with constants

        }
      }

      // 4D vertices: add timing information
      if (f4D and v.isValid()) {
        auto err = v.positionError().matrix4D();
        err(3, 3) = vartime;
        auto trkWeightMap3d = v.weightMap();  // copy the 3d-fit weights
        v = TransientVertex(v.position(), meantime, err, v.originalTracks(), v.totalChiSquared(), v.degreesOfFreedom());
        v.weightMap(trkWeightMap3d);
      }

      if (fVerbose) {
        if (v.isValid()) {
          std::cout << "x,y,z";
          if (f4D)
            std::cout << ",t";
          std::cout << "=" << v.position().x() << " " << v.position().y() << " " << v.position().z();
          if (f4D)
            std::cout << " " << v.time();
          std::cout << " cluster size = " << (*iclus).size() << std::endl;
        } else {
          std::cout << "Invalid fitted vertex,  cluster size=" << (*iclus).size() << std::endl;
        }
      }

      if (v.isValid() && (v.degreesOfFreedom() >= algorithm.minNdof) &&
          (!validBS || (*(algorithm.vertexSelector))(v, beamVertexState)))
        pvs.push_back(v);
    }

    return pvs;
}
//#endif
}
