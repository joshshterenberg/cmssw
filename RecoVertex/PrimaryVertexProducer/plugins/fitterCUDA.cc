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
__global__ void fitterKernel(
    TransientVertex cuda_pvs[],
    algo algorithm,
    reco::TransientTrack cuda_clusters[][max_size],
    reco::BeamSpot beamSpot,
    VertexState beamVertexState,
    int clusters_size,
    int max_size,
    bool f4D,
    bool validBS,
    bool weightFit,
    bool fVerbose
){

      int idx = blockIdx.x * blockDim.x + threadIdx.x;


//    for (std::vector<std::vector<reco::reco::TransientTrack> >::const_iterator iclus = clusters.begin();
//         iclus != clusters.end();
//         iclus++) {


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

      if (fVerbose) {
        if (v.isValid()) {
          std::cout << "x,y,z";
          if (f4D)
            std::cout << ",t";
          std::cout << "=" << v.position().x() << " " << v.position().y() << " " << v.position().z();
          if (f4D)
            std::cout << " " << v.time();
          std::cout << " cluster size = " << (*cuda_clusters[i]).size() << std::endl;
        } else {
          std::cout << "Invalid fitted vertex,  cluster size=" << (*cuda_clusters[i]).size() << std::endl;
        }
      }

      if (v.isValid() && (v.degreesOfFreedom() >= algorithm.minNdof) &&
          (!validBS || (*(algorithm.vertexSelector))(v, beamVertexState)))
        cuda_pvs[idx] = v; //pvs.push_back(v);
    }
}
*/

#ifdef __CUDACC__
std::vector<TransientVertex> wrapper(
    algo algorithm,
    std::vector<std::vector<reco::TransientTrack> >&& clusters,
    reco::BeamSpot beamSpot,
    VertexState beamVertexState,
    bool f4D,
    bool validBS,
    bool weightFit,
    bool fVerbose ){

    std::cout << "\n\n\n\n\ngot to the wrapper\n";

    //defines grid
    unsigned int blockSize = 1;
    unsigned int gridSize  = 1;
    std::cout << "defined grid size\n";

    //create and allocate all host memory (only pvs needed)
    TransientVertex *cpu_pvs;
    cpu_pvs = (TransientVertex *) malloc(clusters.size() * sizeof(TransientVertex));
    std::cout << "created and allocated all host memory\n";


    //create and allocate all device memory (pvs and clusters)
    reco::TransientTrack **cuda_clusters;
    TransientVertex *cuda_pvs;
    std::cout << "initial defs good\n";
    cudaCheck(cudaMalloc(&cuda_pvs, clusters.size() * sizeof(TransientVertex)));
    cudaCheck(cudaMalloc(&cuda_clusters, clusters.size() * sizeof(reco::TransientTrack *)));
    std::cout << "initial cuda mallocs good\n";

    reco::TransientTrack sample_track;
    long unsigned int max_size = 0;
    for (long unsigned int i = 0; i < clusters.size(); i++) if (max_size < clusters[i].size()) max_size = clusters[i].size();
    for (long unsigned int i = 0; i < clusters.size(); i++) {
        cudaCheck(cudaMalloc(&(cuda_clusters[i]), (max_size+1) * sizeof(reco::TransientTrack))); //this is segfaulting rn
        std::cout << "cudamalloc of clusters 2d is good\n";
        for (long unsigned int j = 0; j < clusters[i].size(); j++) {
            std::cout << "everything until population good\n";
            cuda_clusters[i][j] = clusters[i][j]; //populating rectangular 2D cluster array on valid vals
        }
        cuda_clusters[i][clusters[i].size()] = sample_track; //SET FOR FINAL VALUES, COMPARE IN KERNEL
    }
    std::cout << "created and allocated all device memory\n";

    //host to device memory copy (NONE NEEDED)
    //cudaCheck(cudaMemcpy(cuda_clusters, cpu_clusters, memSize1, cudaMemcpyHostToDevice));
    //cudaCheck(cudaMemcpy(cuda_pvs, cpu_pvs, memSize2, cudaMemcpyHostToDevice));
    //std::cout << "host memory copied to device memory\n";


    //action!
    /*
    fitterKernel<<<gridSize, blockSize>>>(
      cuda_pvs,
      algorithm,
      cuda_clusters,
      beamSpot,
      beamVertexState,
      clusters.size(),
      max_size,
      f4D,
      validBS,
      weightFit,
      fVerbose
    );
    std::cout << "main action complete\n";
    */

    //wait for device to complete / error check
    cudaDeviceSynchronize();
    cudaCheck(cudaGetLastError());
    std::cout << "sync / error check complete\n";

    //device to host memory copy
    cudaCheck(cudaMemcpy(cpu_pvs, cuda_pvs, clusters.size() * sizeof(TransientVertex), cudaMemcpyDeviceToHost));
    std::cout << "device memory copied to host memory\n";

    //clear device memory
    for (long unsigned int i = 0; i < clusters.size(); i++) {
        cudaCheck(cudaFree(cuda_clusters[i]));
    }
    cudaCheck(cudaFree(cuda_clusters));
    cudaCheck(cudaFree(cuda_pvs));
    std::cout << "cleared device memory\n\n";

    //move output array to vector
    std::vector<TransientVertex> pvs;
    for (long unsigned int i = 0; i < clusters.size(); i++) pvs.push_back(cpu_pvs[i]); //if statement needed
    std::cout << "pvs vector created\n";

    //done
    return pvs;
}
#endif
}




