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
#include <math.h>

namespace fitterCUDA {



__global__ void fitterKernel(
    unsigned int ntracks,
    TrackForPV::TrackForPVSoA* tracks,
    TrackForPV::VertexForPVSoA* vertices,
    algo algorithm
){
//      return;
// RECOMMENT IF NEEDED FROM HERE
  

  size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x; // This is going to be the vertex index
  size_t gridSize = blockDim.x * gridDim.x;

  
  //1 block for each vertex (for block in BLOCKS)
  for (unsigned int k = firstElement; k < vertices->nTrueVertex(0); k += gridSize) {
    if (!vertices->isGood(k)) continue; //skip if not good
    unsigned int ivertex = vertices->order(k);
    //z weighted averaging
    //currently approxs Gaussian with linear (triangular) PDF
    //only applies PDF to z-axis (ignores x, y, t)
    unsigned int iavg_z = 0;
    //1 thread / track (for each thread in block)
    for (unsigned int kk = 0; kk < vertices->ntracks(ivertex); kk++){
      unsigned int itrack = tracks->order(kk);
      //-----crit load, no thread dependence-----
      double influence = 0;
      double dist = pow(tracks->z(itrack) - vertices->z(ivertex), 2) / (tracks->dz(itrack) + tracks->dzError(itrack));
      if (dist <= 9.0) {
        influence = 1.0 - dist / 9.0;
      }
      iavg_z += tracks->z(itrack) * influence;
    }
    vertices->z(ivertex) = iavg_z / vertices->ntracks(ivertex);
  }
  __syncthreads();

  




//    for (std::vector<std::vector<reco::reco::TransientTrack> >::const_iterator iclus = clusters.begin();
//         iclus != clusters.end();
//         iclus++) {

/*
      for (int i = 0; i < clusters_size; i++) {
      double sumwt = 0.;
      double sumwt2 = 0.;
      double sumw = 0.;
      double meantime = 0.;
      double vartime = 0.;
      if (f4D) {
        for (int j = 0; j < max_size; j++) {
          const double time = cuda_clusters[i][j].timeExt();
          const double err = cuda_clusters[i][j].dtErrorExt();
          const double inverr = err > 0. ? 1.0 / err : 0.;
          const double w = inverr * inverr;
          sumwt += w * time;
          sumwt2 += w * time * time;
          sumw += w;
        }
        meantime = sumwt / sumw;
        double sumsq = sumwt2 - sumwt * sumwt / sumw;
        double chisq = max_size > 1 ? sumsq / double(max_size - 1) : sumsq / double(max_size);
        vartime = chisq / sumw;
      }

      TransientVertex v;

      if (weightFit) {
        std::vector<std::pair<GlobalPoint, GlobalPoint>> points;
        if (algorithm.useBeamConstraint && validBS && (max_size > 1)) {
            for (int j = 0; j < max_size; j++){
                   GlobalPoint p =  cuda_clusters[i][j].stateAtBeamLine().trackStateAtPCA().position();
                   GlobalPoint err(cuda_clusters[i][j].stateAtBeamLine().transverseImpactParameter().error(), cuda_clusters[i][j].stateAtBeamLine().transverseImpactParameter().error(), cuda_clusters[i][j].track().dzError());
                   std::pair<GlobalPoint, GlobalPoint> p2(p, err);
                   points.push_back(p2);
            }
            v = WeightedMeanFitter::weightedMeanOutlierRejectionBeamSpot(points, *cuda_clusters[i], beamSpot);
            if ((v.positionError().matrix())(2,2) != (WeightedMeanFitter::startError*WeightedMeanFitter::startError)) cuda_pvs[idx] = v; //pvs.push_back(v);

        }
        else if (!(algorithm.useBeamConstraint) && (cuda_clusters[i]->size() > 1)) {
           for (int j = 0; j < max_size; j++){
                   GlobalPoint p = cuda_clusters[i][j].impactPointState().globalPosition();
                   GlobalPoint err(cuda_clusters[i][j].track().dxyError(), cuda_clusters[i][j].track().dxyError(), cuda_clusters[i][j].track().dzError());
                   std::pair<GlobalPoint, GlobalPoint> p2(p, err);
                   points.push_back(p2);
           }
           v = WeightedMeanFitter::weightedMeanOutlierRejection(points, *cuda_clusters[i]);
           if ((v.positionError().matrix())(2,2) != (WeightedMeanFitter::startError*WeightedMeanFitter::startError)) cuda_pvs[idx] = v; //pvs.push_back(v); //FIX with constants

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

      if (v.isValid() && (v.degreesOfFreedom() >= algorithm.minNdof) &&
          (!validBS || (*(algorithm.vertexSelector))(v, beamVertexState)))
        cuda_pvs[idx] = v; //pvs.push_back(v);
    }
*/
}


#ifdef __CUDACC__
void wrapper(
    unsigned int ntracks,
    TrackForPV::TrackForPVSoA* GPUtracksObject,
    TrackForPV::VertexForPVSoA* GPUverticesObject,
    algo algorithm
){

    //defines grid
    unsigned int blockSize = 1; //optimal size depends, probably 1 block, multiple threads per vertex
    unsigned int gridSize  = 1; //might need experimental determination

    //action!
    fitterKernel<<<gridSize, blockSize>>>(
    	ntracks,
        GPUtracksObject,
    	GPUverticesObject,
    	algorithm
    );

    cudaCheck(cudaGetLastError());
}
#endif
}




