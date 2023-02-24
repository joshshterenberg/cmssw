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



__global__ void fitterKernel(
    unsigned int ntracks,
    TrackForPV::TrackForPVSoA* GPUtracksObject,
    TrackForPV::VertexForPVSoA* GPUverticesObject,
    algo algorithm
){
//      return;
// RECOMMENT IF NEEDED FROM HERE
      int idx = blockIdx.x * blockDim.x + threadIdx.x;

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
    unsigned int blockSize = 1;
    unsigned int gridSize  = 1;
    std::cout << "defined grid size\n";


    /*
    //create and allocate all host memory (only pvs needed)
    TransientVertex *cpu_pvs;
    TrackForPV::TrackForPVSoA* cpu_clusters;

    cpu_pvs = (TransientVertex *) malloc(ntracks * sizeof(TransientVertex));
    cpu_clusters = (reco::TransientTrack *) malloc(clusters.size() * max_size * sizeof(reco::TransientTrack)); //1D flattened


    /////////////////////////////////////////////////////////
    ////// MODIFY TYPE / CONVERT TO WHAT'S NEEDED HERE //////
    /////////////////////////////////////////////////////////

    long unsigned int n_iter[clusters.size()];
    std::cout << "size of clusters array: " << clusters.size() << "X" << max_size << "\n";
    //std::cout << "size of cpu_clusters: " << end(cpu_clusters) - begin(cpu_clusters) << "\n";
    for (long unsigned int i = 0; i < clusters.size(); i++) {
        std::cout << "i is " << i << ". current size is " << clusters[i].size() << "\n";
        n_iter[i] = clusters[i].size(); //use in kernel to keep track of bounds on 2D array
        for (long unsigned int j = 0; j < clusters[i].size(); j++) {
            std::cout << "j is " << j << "\t";
            reco::TransientTrack test = clusters[i][j];
            std::cout << "read ok\t";
            cpu_clusters[i * max_size + j] = test; //populating rectangular 2D cluster array on valid vals
            std::cout << "write ok\n";
        }
    }

    std::cout << "created and allocated all host memory\n";


    //create and allocate all device memory (pvs and clusters)
    reco::TransientTrack *cuda_clusters;
    TransientVertex *cuda_pvs;
    cudaCheck(cudaMalloc(&cuda_pvs, clusters.size() * sizeof(TransientVertex)));
    cudaCheck(cudaMalloc(&cuda_clusters, clusters.size() * max_size * sizeof(reco::TransientTrack)));

    std::cout << "created and allocated all device memory\n";

    //host to device memory copy
    cudaCheck(cudaMemcpy(cuda_clusters, cpu_clusters, clusters.size() * max_size * sizeof(reco::TransientTrack), cudaMemcpyHostToDevice));
    std::cout << "host memory copied to device memory\n";
    */

    //action!
    fitterKernel<<<gridSize, blockSize>>>(
    	ntracks,
        GPUtracksObject,
    	GPUverticesObject,
    	algorithm
    );
    std::cout << "main action complete\n";

    //wait for device to complete / error check
    cudaDeviceSynchronize();
    cudaCheck(cudaGetLastError());
    std::cout << "sync / error check complete\n";


    /*
    //device to host memory copy
    cudaCheck(cudaMemcpy(cpu_pvs, cuda_pvs, clusters.size() * sizeof(TransientVertex), cudaMemcpyDeviceToHost));
    std::cout << "device memory copied to host memory\n";

    //clear device memory
    cudaCheck(cudaFree(cuda_clusters));
    cudaCheck(cudaFree(cuda_pvs));
    std::cout << "cleared device memory\n\n";

    //move output array to vector
    std::vector<TransientVertex> pvs;
    for (long unsigned int i = 0; i < clusters.size(); i++) pvs.push_back(cpu_pvs[i]); //if statement needed
    std::cout << "pvs vector created\n";

    //done
    return pvs;

    */
}
#endif
}




