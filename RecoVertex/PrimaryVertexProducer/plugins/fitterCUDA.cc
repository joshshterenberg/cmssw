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

using Vector512d = Eigen::Matrix<double, 1024, 1>;

namespace fitterCUDA {

__global__ void fitterKernel(
    unsigned int ntracks,
    TrackForPV::TrackForPVSoA* tracks,
    TrackForPV::VertexForPVSoA* vertices,
    algo algorithm
){

  /*
	OUTPUTS:
		DONE nTrueVertex (filled in clusterizer)
		DONE isGood flag (filled in clusterizer, can be modified (for now no mod))
		DONE x, y, z, t
		DONE chi2 (dep. on errs)
		DONE ndof (+1 per track if track has any influence)
		DONE errx, erry, errz
		DONE ntracks (same as ndof)
		DONE track_id (type Vector512d, save position here. ignoring...)
		DONE track_weight (type Vector512d, save weight in order here. ignoring...)
  */


  size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x;
  size_t gridSize = blockDim.x * gridDim.x;
  float precision = 1e-24;

  //1 block for each vertex (for block in BLOCKS)
  for (unsigned int k = firstElement; k < vertices->nTrueVertex(0); k += gridSize) {
    if (!vertices->isGood(k)) continue; //skip if not good
    unsigned int ivertex = vertices->order(k);

    //position/vector loop
    //1 thread / track (for each thread in block)
    unsigned int iavg_x = 0, iavg_y = 0, iavg_z = 0, iavgt = 0;
    Vector512d track_ids, track_weights;
    for (unsigned int kk = 0; kk < vertices->ntracks(ivertex); kk++){
      unsigned int itrack = tracks->order(kk);
      track_ids[itrack] = itrack;

      //GAUSSEAN WEIGHTING (just in z for now)
      double influence = 0;
      float wz = (tracks->dzError(itrack) <= precision) ? pow(precision,2) : pow(tracks->dzError(itrack),2);
      double dist = pow(tracks->dz(itrack) - vertices->z(ivertex), 2) / //change to x/dx?
                    (wz); //add err_z?
      if (dist <= 9.0) {
        influence = 1.0 - dist / 9.0;
        vertices->ndof(ivertex) += 1;
        vertices->ntracks(ivertex) += 1;
      }

      track_weights[itrack] = influence;

      iavg_x += tracks->x(itrack) * influence;
      iavg_y += tracks->y(itrack) * influence;
      iavg_z += tracks->z(itrack) * influence;
      //iavgt  += tracks->t(itrack) * influence; //not using
    }
    vertices->x(ivertex) = iavg_x / vertices->ntracks(ivertex);
    vertices->y(ivertex) = iavg_y / vertices->ntracks(ivertex);
    vertices->z(ivertex) = iavg_z / vertices->ntracks(ivertex);
    //vertices->t(ivertex) = iavgt  / vertices->ntracks(ivertex);

    vertices->track_id(ivertex) = track_ids;
    vertices->track_weight(ivertex) = track_weights;

    __syncthreads(); //not sure if needed, just in case RAW

    //errs loop, similar deal
    double s1x = 0, s1y = 0, s1z = 0, s2x = 0, s2y = 0, s2z = 0;
    for (unsigned int kk = 0; kk < vertices->ntracks(ivertex); kk++){
      unsigned int itrack = tracks->order(kk);
      s1x += tracks->dxError(itrack); //maybe change to dx?
      s1y += tracks->dyError(itrack);
      s1z += tracks->dzError(itrack);
      s2x += (s1x + tracks->weight(itrack));
      s2y += (s1y + tracks->weight(itrack));
      s2z += (s1z + tracks->weight(itrack));
    }
    vertices->errx(ivertex) = s1x / s2x;
    vertices->erry(ivertex) = s1y / s2y;
    vertices->errz(ivertex) = s1z / s2z;

    __syncthreads();


    //chi2 loop, similar deal
    double dist = 0;
    for (unsigned int kk = 0; kk < vertices->ntracks(ivertex); kk++){
      unsigned int itrack = tracks->order(kk);
      float wx = (tracks->dxError(itrack) <= precision) ? precision : tracks->dxError(itrack); //maybe change to dx?
      float wy = (tracks->dyError(itrack) <= precision) ? precision : tracks->dyError(itrack);
      float wz = (tracks->dzError(itrack) <= precision) ? precision : tracks->dzError(itrack);
      dist +=       pow(tracks->dx(itrack) - vertices->x(ivertex), 2) / //maybe change to x?
                    (pow(wx, 2) + vertices->errx(ivertex));
      dist +=       pow(tracks->dy(itrack) - vertices->y(ivertex), 2) /
                    (pow(wy, 2) + vertices->erry(ivertex));
      dist +=       pow(tracks->dz(itrack) - vertices->z(ivertex), 2) /
                    (pow(wz, 2) + vertices->errz(ivertex));
    }
    vertices->chi2(ivertex) = dist; //likely eliminated RAW dependence here

    __syncthreads();

  }





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
    unsigned int blockSize = 512; //optimal size depends, probably 1 block, multiple threads per vertex
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




