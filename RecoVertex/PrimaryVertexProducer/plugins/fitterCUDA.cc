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
  double precision = 1e-24; //FLOAT CHANGE

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
    double s_wx = 0, s_wy = 0, s_wz = 0;

    for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){
      unsigned int itrack = tracks->order(kk);
      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == k) {
        //this is a valid track that is associated with this vertex
        track_ids[track_id_counter] = itrack;
        track_id_counter++;
        //double wx = (tracks->dxy2(itrack) > precision) ? 1./tracks->dxy2(itrack) : 1./precision;
        //double wy = (tracks->dxy2(itrack) > precision) ? 1./tracks->dxy2(itrack) : 1./precision;
        double wz = (tracks->dz2(itrack) > precision) ? 1./tracks->dz2(itrack) : 1./precision;
        old_z += tracks->z(itrack) * wz;
        old_x += tracks->x(itrack) * 1; //mod
        old_y += tracks->y(itrack) * 1;
        s_wz += wz;
        s_wx += 1;
        s_wy += 1;
      }
    }


    old_x /= s_wx;
    old_y /= s_wy;
    old_z /= s_wz;

    err_x = 1/s_wx;
    err_y = 1/s_wy;
    err_z = 1/s_wz;

    vertices->track_id(ivertex) = track_ids;

    //loop to calculate vertex positions

    double iavg_x = 0, iavg_y = 0, iavg_z = 0;
    s_wz = 0; s_wy = 0; s_wx = 0;
    double s2x = 0, s2y = 0, s2z = 0;

    for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){
      unsigned int itrack = tracks->order(kk);
      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == k) {
        ////WEIGHTING
        double influence = 0;
        double wx = (tracks->dxy2(itrack) <= precision) ? precision : tracks->dxy2(itrack);
        double wy = (tracks->dxy2(itrack) <= precision) ? precision : tracks->dxy2(itrack);
        double wz = (tracks->dz2(itrack) <= precision) ? precision : tracks->dz2(itrack);
        double distz = pow(tracks->z(itrack) - old_z, 2) / (wz + err_z);
        double distx = pow(tracks->x(itrack) - old_x, 2) / (wx + err_x);
        double disty = pow(tracks->y(itrack) - old_y, 2) / (wy + err_y);

        if (distz + distx + disty < 9.0) {
          influence = 1;
          vertices->ndof(ivertex) += 1;
          vertices->ntracks(ivertex) += 1;
          track_weights[track_weight_counter] = influence;
        } else {
          track_weights[track_weight_counter] = 0;
        }
        track_weight_counter++;
        ////

        iavg_x += tracks->x(itrack) * (influence / wx);
        iavg_y += tracks->y(itrack) * (influence / wy);
        iavg_z += tracks->z(itrack) * (influence / wz);

        s_wx += influence / wx;
        s_wy += influence / wy;
        s_wz += influence / wz;

        s2x += influence * influence / wx;
        s2y += influence * influence / wy;
        s2z += influence * influence / wz;
      }
    }

    vertices->x(ivertex) = old_x; //mod_ just the average pos
    vertices->y(ivertex) = old_y;
    vertices->z(ivertex) = iavg_z / s_wz;

    vertices->track_weight(ivertex) = track_weights;

    vertices->errx(ivertex) = s2x * pow(corr_x,2) / pow(s_wx,2);
    vertices->erry(ivertex) = s2y * pow(corr_y,2) / pow(s_wy,2);
    vertices->errz(ivertex) = s2z * pow(corr_z,2) / pow(s_wz,2);

    //----------------------------------chi2 loop-------------------------------------
    //creates error: GammaContinuedFraction::a too large, ITMAX too small
    //maybe just average chi2s from the tracks? those are being calculated anyway.

    double dist = 0, chi2 = 0;
    for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){
      unsigned int itrack = tracks->order(kk);
      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == k) {

        double wx = (tracks->dxy2(itrack) <= precision) ? precision : tracks->dxy2(itrack);
        double wy = (tracks->dxy2(itrack) <= precision) ? precision : tracks->dxy2(itrack);
        double wz = (tracks->dz2(itrack) <= precision) ? precision : tracks->dz2(itrack);

        printf("x,y,z: %f, %f, %f\n", tracks->x(itrack), tracks->y(itrack), tracks->z(itrack));
        printf("vx,vy,vz: %f, %f, %f\n",  vertices->x(ivertex),  vertices->y(ivertex),  vertices->z(ivertex));
        printf("wx,wy,wz: %f, %f, %f\n", wx, wy, wz);
        printf("errx,erry,errz: %f, %f, %f\n", vertices->errx(ivertex), vertices->erry(ivertex), vertices->errz(ivertex));

        dist =        pow(tracks->x(itrack) - vertices->x(ivertex), 2) /
                      (wx + vertices->errx(ivertex));
        printf("dist post x: %f\n", dist);
        dist +=       pow(tracks->y(itrack) - vertices->y(ivertex), 2) /
                      (wy + vertices->erry(ivertex));
        printf("dist post y: %f\n", dist);
        dist +=       pow(tracks->z(itrack) - vertices->z(ivertex), 2) /
                      (wz + vertices->errz(ivertex));
        printf("dist post z: %f\n", dist);
        chi2 += dist;
        //chi2 += tracks->chi2(itrack); //TODO: TAKE OUT OF DATAFORMAT
      }
    }


    vertices->chi2(ivertex) = chi2; //requirement
    //weight tracks by chi2?


    printf("%d, %f, %f, %f, %d, %f, %f\n",
      ivertex,
      vertices->x(ivertex),
      vertices->y(ivertex),
      vertices->z(ivertex),
      vertices->ntracks(ivertex),
      vertices->errz(ivertex),
      vertices->chi2(ivertex)
    );

  }
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




