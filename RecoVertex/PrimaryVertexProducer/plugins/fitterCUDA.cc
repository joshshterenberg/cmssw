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
  double precision = 1.0; //FLOAT CHANGE

  //1 block for each vertex (for block in BLOCKS)
  for (unsigned int k = firstElement; k < vertices->nTrueVertex(0); k += gridSize) {
    unsigned int ivertex = vertices->order(k);
    if (!vertices->isGood(ivertex)) continue; //skip if not good

    /*
		WEIGHTEDMEANFITTER ORDER

	fill err(x,y,z) with 4,4,400, corr(x,y,z) with 1.2,1.2,1.4
	x = sum(xi * (1/dxi^2)) / sum((1/dxi^2)), continue with y,z (check for precision)
	define err_x = 1 / sum((1/dxi^2)), continue with y,z
	for every point calc weight (xi - x) / (dxi^2 + err_x) < 9 ?
		set (ntracks, ndof) if > 0
	x = sum(xi * weight / (dxi^2)) / sum(dxi^2)
	err(x,y,z) = sum(dxi^2 * weight) * corr(x,y,z)^2 / (dxi^2)^2
	chi2 = sum((xi-x)^2/(dxi^2+err(x)) + y,z)

    */


    //ASSUMING THROUGHOUT THAT tracks->x == p.first.x() and tracks->dx == p.second.x()


    //-----------------------------position/vector loop-----------------------------
    //1 thread / track (for each thread in block)
    double err_x = 4.0, err_y = 4.0, err_z = 400.0, corr_x = 1.2, corr_y = 1.2, corr_z = 1.4;
    Vector512d track_ids, track_weights;
    int track_id_counter = 0, track_weight_counter = 0;

    //loop to calculate z weighted average ("old_z")
    double old_z = 0, old_x = 0, old_y = 0;
    double s_wz = 0;
    for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){
      unsigned int itrack = tracks->order(kk);
      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == k) {
        //this is a valid track that is associated with this vertex
        track_ids[track_id_counter] = itrack;
        track_id_counter++;
        //double wx = (tracks->dx2(itrack) < precision) ? 1./pow(tracks->dx2(itrack),2) : 1./pow(precision,2);
        //double wy = (tracks->dy2(itrack) < precision) ? 1./pow(tracks->dy2(itrack),2) : 1./pow(precision,2);
        double wz = (tracks->dz2(itrack) < precision) ? 1./tracks->dz2(itrack) : 1./precision; //dz2 is already squared, precision is too
        old_z += tracks->z(itrack) * wz;
        //old_x += tracks->x(itrack) * wx;
        //old_y += tracks->y(itrack) * wy;
        s_wz += wz;
        //s_wx += wx;
        //s_wy += wy;
      }
    }
    //old_x /= s_wx;
    //old_y /= s_wy;
    old_z /= s_wz;

    //err_x = 1/s_wx;
    //err_y = 1/s_wy;
    err_z = 1/s_wz;

    vertices->track_id(ivertex) = track_ids;

    //loop to calculate vertex positions

    double iavg_x = 0, iavg_y = 0, iavg_z = 0;
    s_wz = 0;
    int s1x = 0, s1y = 0, s1z = 0;
    for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){
      unsigned int itrack = tracks->order(kk);
      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == k) {
        ////GAUSSEAN WEIGHTING (just triangular in z for now)
        double influence = 0;
        //double wx = (tracks->dx(itrack) <= precision) ? pow(precision,2) : pow(tracks->dx(itrack),2);
        //double wy = (tracks->dy(itrack) <= precision) ? pow(precision,2) : pow(tracks->dy(itrack),2);
        double wz = (tracks->dz2(itrack) < precision) ? tracks->dz2(itrack) : precision; //???
        double distz = pow(tracks->z(itrack) - old_z, 2) / (wz + err_z);
        //double distx = pow(tracks->x(itrack) - old_x, 2) / (wx + err_x);
        //double disty = pow(tracks->y(itrack) - old_y, 2) / (wy + err_y);
        if (distz < 9.0) {
          influence = 1.0; //straight in z for now (- abs(distz)/3.0)
          vertices->ndof(ivertex) += 1;
          vertices->ntracks(ivertex) += 1;
          track_weights[track_weight_counter] = influence;
        } else {
          track_weights[track_weight_counter] = 0;
          //continue;
        }
        track_weight_counter++;
        ////

        //iavg_x += tracks->x(itrack) * (influence / wx);
        //iavg_y += tracks->y(itrack) * (influence / wy);
        iavg_z += tracks->z(itrack) * (influence / wz);

        //s_wx += wx;
        //s_wy += wy;
        s_wz += influence / wz; //CHECK

        //s1x += wx * influence;
        //s1y += wy * influence;
        s1z += influence * influence / wz; //CHECK
      }
    }
    //vertices->x(ivertex) = iavg_x / s_wx;
    //vertices->y(ivertex) = iavg_y / s_wy;
    vertices->z(ivertex) = iavg_z / s_wz;

    vertices->track_weight(ivertex) = track_weights;

    //vertices->errx(ivertex) = s1x * pow(corr_x,2) / pow(s_wx,2);
    //vertices->erry(ivertex) = s1y * pow(corr_y,2) / pow(s_wy,2);
    vertices->errz(ivertex) = s1z * pow(corr_z,2) / pow(s_wz,2);

    printf("current vertex: %d\t old_z: %f\t vertices->z: %f\n", ivertex, old_z, vertices->z(ivertex));

    //----------------------------------chi2 loop-------------------------------------
    //creates error: GammaContinuedFraction::a too large, ITMAX too small
    /*
    double dist = 0, chi2 = 0;
    for (unsigned int kk = 0; kk < vertices->ntracks(ivertex); kk++){
      unsigned int itrack = tracks->order(kk);
      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == k) {
        double wx = (tracks->dx(itrack) <= precision) ? precision : tracks->dx(itrack);
        double wy = (tracks->dy(itrack) <= precision) ? precision : tracks->dy(itrack);
        double wz = (tracks->dz(itrack) <= precision) ? precision : tracks->dz(itrack);
        dist =        pow(tracks->x(itrack) - old_x, 2) /
                      (pow(wx, 2) + vertices->errx(ivertex));
        dist +=       pow(tracks->y(itrack) - old_y, 2) /
                      (pow(wy, 2) + vertices->erry(ivertex));
        dist +=       pow(tracks->z(itrack) - old_z, 2) /
                      (pow(wz, 2) + vertices->errz(ivertex));
        chi2 += dist;
      }
    }
    vertices->chi2(ivertex) = chi2;
    */

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




