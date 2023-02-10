#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerDumbFitter.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"

#include "RecoVertex/LinearizationPointFinders/interface/CrossingPtBasedLinearizationPointFinder.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include <cmath>
#include <algorithm>

namespace {
  inline GlobalPoint operator-(const GlobalPoint& a, const GlobalPoint& b) {
    return GlobalPoint(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
  }

  inline GlobalPoint operator+(const GlobalPoint& a, const GlobalPoint& b) {
    return GlobalPoint(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
  }

  inline GlobalPoint operator/(const GlobalPoint& a, const double b) {
    return GlobalPoint(a.x() / b, a.y() / b, a.z() / b);
  }

#ifndef __clang__
  inline GlobalPoint operator*(const GlobalPoint& a, const double b) {
    return GlobalPoint(a.x() * b, a.y() * b, a.z() * b);
  }

  inline GlobalPoint operator*(const double b, const GlobalPoint& a) {
    return GlobalPoint(a.x() * b, a.y() * b, a.z() * b);
  }
#endif
}


struct CompareTwoTracks {
    int operator()(const reco::TransientTrack &a, const reco::TransientTrack &b) {
      return a.initialFreeState().momentum().mag() > b.initialFreeState().momentum().mag();
      //       return a.p() > b.p();
    };
};

PrimaryVertexProducerDumbFitter::PrimaryVertexProducerDumbFitter(const edm::ParameterSet& conf)
    : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))), theConfig(conf) {
  fVerbose = conf.getUntrackedParameter<bool>("verbose", false);

  trkToken = consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("TrackLabel"));
  bsToken = consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpotLabel"));
  f4D = false;

  // select and configure the track selection
  std::string trackSelectionAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
  if (trackSelectionAlgorithm == "filter") {
    theTrackFilter = new TrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else if (trackSelectionAlgorithm == "filterWithThreshold") {
    theTrackFilter = new HITrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else {
    throw VertexException("PrimaryVertexProducer: unknown track selection algorithm: " + trackSelectionAlgorithm);
  }

  // select and configure the track clusterizer
  std::string clusteringAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<std::string>("algorithm");
  if (clusteringAlgorithm == "gap") {
    theTrackClusterizer = new GapClusterizerInZ(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkGapClusParameters"));
  } else if (clusteringAlgorithm == "DA") {
    theTrackClusterizer = new DAClusterizerInZ(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  }
  // provide the vectorized version of the clusterizer, if supported by the build
  else if (clusteringAlgorithm == "DA_vect") {
    theTrackClusterizer = new DAClusterizerInZ_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  } else if (clusteringAlgorithm == "DASub_vect") {
    theTrackClusterizer = new DAClusterizerInZSubCluster_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDASubClusParameters"));
  } else if (clusteringAlgorithm == "DA2D_vect") {
    theTrackClusterizer = new DAClusterizerInZT_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
    f4D = true;
  }

  else {
    throw VertexException("PrimaryVertexProducer: unknown clustering algorithm: " + clusteringAlgorithm);
  }

  if (f4D) {
    trkTimesToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("TrackTimeResosLabel"));
  }

  // select and configure the vertex fitters
  if (conf.exists("vertexCollections")) {
    std::vector<edm::ParameterSet> vertexCollections =
        conf.getParameter<std::vector<edm::ParameterSet> >("vertexCollections");

    for (std::vector<edm::ParameterSet>::const_iterator algoconf = vertexCollections.begin();
         algoconf != vertexCollections.end();
         algoconf++) {
      algo algorithm;
      std::string fitterAlgorithm = algoconf->getParameter<std::string>("algorithm");
      if (fitterAlgorithm == "KalmanVertexFitter") {
        algorithm.fitter = new KalmanVertexFitter();
      } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
        algorithm.fitter = new AdaptiveVertexFitter(GeometricAnnealing(algoconf->getParameter<double>("chi2cutoff")));
      } else if (fitterAlgorithm == "weightedMean") {
        std::cout << "Using weighted mean as fitter" << std::endl;
      } else {
        throw VertexException("PrimaryVertexProducer: unknown algorithm: " + fitterAlgorithm);
      }
      algorithm.label = algoconf->getParameter<std::string>("label");
      algorithm.minNdof = algoconf->getParameter<double>("minNdof");
      algorithm.useBeamConstraint = algoconf->getParameter<bool>("useBeamConstraint");
      algorithm.vertexSelector =
          new VertexCompatibleWithBeam(VertexDistanceXY(), algoconf->getParameter<double>("maxDistanceToBeam"));
      algorithms.push_back(algorithm);

      produces<reco::VertexCollection>(algorithm.label);
    }
  } else {
    edm::LogWarning("MisConfiguration")
        << "this module's configuration has changed, please update to have a vertexCollections=cms.VPSet parameter.";

    algo algorithm;
    std::string fitterAlgorithm = conf.getParameter<std::string>("algorithm");
    if (fitterAlgorithm == "KalmanVertexFitter") {
      algorithm.fitter = new KalmanVertexFitter();
    } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
      algorithm.fitter = new AdaptiveVertexFitter();
    } else if (fitterAlgorithm == "weightedMean") {
        std::cout << "Using weighted mean as fitter" << std::endl;
    } else {
      throw VertexException("PrimaryVertexProducerAlgorithm: unknown algorithm: " + fitterAlgorithm);
    }
    algorithm.label = "";
    algorithm.minNdof = conf.getParameter<double>("minNdof");
    algorithm.useBeamConstraint = conf.getParameter<bool>("useBeamConstraint");

    algorithm.vertexSelector = new VertexCompatibleWithBeam(
        VertexDistanceXY(),
        conf.getParameter<edm::ParameterSet>("PVSelParameters").getParameter<double>("maxDistanceToBeam"));

    algorithms.push_back(algorithm);
    produces<reco::VertexCollection>(algorithm.label);
  }

  //check if this is a recovery iteration
  fRecoveryIteration = conf.getParameter<bool>("isRecoveryIteration");
  if (fRecoveryIteration) {
    if (algorithms.empty()) {
      throw VertexException("PrimaryVertexProducer: No algorithm specified. ");
    } else if (algorithms.size() > 1) {
      throw VertexException(
          "PrimaryVertexProducer: Running in Recovery mode and more than one algorithm specified.  Please "
          "only one algorithm.");
    }
    recoveryVtxToken = consumes<reco::VertexCollection>(conf.getParameter<edm::InputTag>("recoveryVtxCollection"));
  }
}

TransientVertex weightedMean(const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points, std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus){
     float x=0, y=0, z=0, s_wxy=0, s_wz=0, s2_wxy=0, s2_wz=0, wxy=0, wz=0, chi2=0;
     float ndof_z = 0, ndof_xy = 0;
     AlgebraicSymMatrix33 err;
     err(0,0) = 2 * 2;
     err(1,1) = 2 * 2;
     err(2,2) = 20 * 20; // error is 20 cm, so cov -> is 20 ^ 2
     //bool isValid = false;
     for (const auto& p : points){ 
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
            //GlobalError e(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(),  itrack.track().dzError());

            wxy = p.second.x();
            wxy = wxy <=  1e-4 ? 1e+8 : 1. / std::pow(wxy,2);

            wz = p.second.z();
            wz = wz <=  1e-4 ? 1e+8 : 1. / std::pow(wz,2);

            x += p.first.x() * wxy;
            y += p.first.y() * wxy;
            z += p.first.z() * wz;

            s_wxy += wxy;
            s_wz  += wz;
            s2_wxy += wxy * wxy;
            s2_wz  += wz  * wz;
     }

     if ( s_wxy == 0. || s_wz == 0. ){
        std::cout << "Vertex fitting failed at beginning" << std::endl;
        return TransientVertex(GlobalPoint(0,0,0), err, *iclus, 0, 0);
     }

     ndof_xy = (s_wxy * s_wxy)/s2_wxy;
     ndof_z  = (s_wz  * s_wz )/s2_wz;

     x /= s_wxy;     
     y /= s_wxy;     
     z /= s_wz;     
     float err_x=0, err_y=0, err_z=0;


     for (const auto& p : points){ 
        wxy = p.second.x();
        wxy = wxy <=  1e-4 ? 1e+8 : 1. / std::pow(wxy,2);

        wz = p.second.z();
        wz = wz <=  1e-4 ? 1e+8 : 1. / std::pow(wz,2);

        err_x += wxy * pow(p.first.x() - x, 2);
        err_y += wxy * pow(p.first.y() - y, 2);
        err_z += wz  * pow(p.first.z() - z, 2);
     }
//     err(0,0) = ( err_x / (s_wxy * (points.size() -1) ) );
//     err(1,1) = ( err_y / (s_wxy * (points.size() -1) ) );
//     err(2,2) = ( err_z / (s_wz  * (points.size() -1) ) );

     err(0,0) = ( err_x / (s_wxy * (ndof_xy-1) ) );
     err(1,1) = ( err_y / (s_wxy * (ndof_xy-1) ) );
     err(2,2) = ( err_z / (s_wz  * (ndof_z -1) ) );

     //isValid = true;
     float dist = 0; 
     //float varikk = 0;
     //float var_x = 0, var_y = 0, var_z = 0;
     //std::vector<float> weights;
     for (const auto& p : points){ 
          //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
        //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
        wxy = p.second.x();
        wz = p.second.z();
        dist =  std::pow(p.first.x() - x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0), 2) );
        dist += std::pow(p.first.y() - y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1), 2) );
        dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz , 2) + std::pow( err(2,2), 2) ); 
//        dist =  std::pow(p.first.x() - x, 2) / (std::pow( err(0,0), 2) );
//        dist += std::pow(p.first.y() - y, 2) / (std::pow( err(1,1), 2) );
//        dist += std::pow(p.first.z() - z, 2) / (std::pow( err(2,2), 2) ); 
        //weights.push_back(dist)

        /*
        dist = std::pow(p.first.x() - x, 2) / p.first.x();
        dist += std::pow(p.first.y() - y, 2) / p.first.y();
        dist += std::pow(p.first.z() - z, 2) / p.first.z(); 
        */
        chi2 += dist;
        //ndof ++;
        
       //var_x += wxy * pow(p.first.x() - x, 2)
     }
     TransientVertex v(GlobalPoint(x,y,z), err, *iclus, chi2, int((2 * ndof_xy + ndof_z)));
     return v;
}



TransientVertex weightedMeanManyIter_originalIP(const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points, std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus){
     float x=0, y=0, z=0, s_wxy=0, s_wz=0, s2_wxy=0, s2_wz=0, wxy=0, wz=0, chi2=0;
     float ndof_z = 0, ndof_xy = 0;
     AlgebraicSymMatrix33 err;
     err(0,0) = 2 * 2;
     err(1,1) = 2 * 2;
     err(2,2) = 20 * 20; // error is 20 cm, so cov -> is 20 ^ 2

     for (const auto& p : points){ 
            //std::vector<reco::TransientTrack>::const_iterator itrack = iclus.begin(); itrack!= iclus.end(); itrack++) {
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
            //GlobalError e(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(),  itrack.track().dzError());

            wxy = p.second.x();
            wxy = wxy <=  1e-4 ? 1e+8 : 1. / std::pow(wxy,2);

            wz = p.second.z();
            wz = wz <=  1e-4 ? 1e+8 : 1.  / std::pow(wz,2);

            x += p.first.x() * wxy;
            y += p.first.y() * wxy;
            z += p.first.z() * wz;

            s_wxy += wxy;
            s_wz  += wz;
     }

     
     if ( s_wxy == 0. || s_wz == 0. ){
        std::cout << "Vertex fitting failed at beginning" << std::endl;
        return TransientVertex(GlobalPoint(0,0,0), err, *iclus, 0, 0);
     }


      x /= s_wxy;     
      y /= s_wxy;     
      z /= s_wz;
      float err_x, err_y, err_z;
      float old_x, old_y, old_z;
      float dist; 
      for (unsigned int iter=0; iter < 10 ; iter++){

          s_wxy = 0; s_wz = 0;
          old_x = x, old_y = y, old_z = z;
          x=0; y=0; z=0;

          for (const auto& p : points){ 
          //for (unsigned int i=0; i < (unsigned int) iclus->size(); i++){ 
            //dist = IPTools::absoluteImpactParameter3D((*iclus)[i],(reco::Vertex) TransientVertex(GlobalPoint(old_x, old_y, old_z), err, *iclus, 0, 0)).second();
//            TrajectoryStateClosestToPoint traj = (*iclus)[i].trajectoryStateClosestToPoint(GlobalPoint(old_x, old_y, old_z));
//            if (! traj.isValid() ) continue;
//            std::pair<GlobalPoint, GlobalPoint> p(traj.position(), GlobalPoint(0,0,0));
            dist =  std::pow(p.first.x() - old_x, 2); 
            dist += std::pow(p.first.y() - old_y, 2); 
            dist += std::pow(p.first.z() - old_z, 2);       
            dist = dist <= 1e-4 ? 1e+8 : 1. / std::pow(dist,2);

            x += p.first.x() * dist;
            y += p.first.y() * dist;
            z += p.first.z() * dist;

            s_wxy += dist ;
            s_wz  += dist ;

            s2_wxy += dist * dist;
            s2_wz  += dist * dist;
          }
          if ( s_wxy == 0. || s_wz == 0. ){
            std::cout << "Vertex fitting failed at iteration " << iter << std::endl;
            return TransientVertex(GlobalPoint(0,0,0), err, *iclus, 0, 0);
          }
         ndof_xy = (s_wxy * s_wxy)/s2_wxy;
         ndof_z  = (s_wz  * s_wz )/s2_wz;
          x /= s_wxy;     
          y /= s_wxy;     
          z /= s_wz;     

          err_x=0.; err_y=0.; err_z=0.;

          for (const auto& p : points){ 
          //for (unsigned int i=0; i < (unsigned int) iclus->size(); i++){ 
//            TrajectoryStateClosestToPoint traj = (*iclus)[i].trajectoryStateClosestToPoint(GlobalPoint(old_x, old_y, old_z));
//            if (! traj.isValid() ) continue;
//            std::pair<GlobalPoint, GlobalPoint> p(traj.position(), GlobalPoint(0,0,0));
            dist =  std::pow(p.first.x() - old_x, 2); 
            dist += std::pow(p.first.y() - old_y, 2);
            dist += std::pow(p.first.z() - old_z, 2);
            dist = dist <= 1e-4 ? 1e+8 : 1. / std::pow(dist,2);
            //dist = dist == 0. ? 10000. : 1. / std::pow(dist,2);
            
            err_x += dist * std::pow(p.first.x() - x, 2);
            err_y += dist * std::pow(p.first.y() - y, 2);
            err_z += dist * std::pow(p.first.z() - z, 2);
          }
          /*
          if (err_x <= 0 || err_y <= 0 || err_z <= 0) std::cout << "\none err negative\n" << std::endl;
          if (s_wxy <= 0 || s_wz <= 0) std::cout << "\none sum negative\n" << std::endl;
          if ((points.size() - 1) <= 0.5) std::cout << "\npoints size <= 0.5\n" << std::endl;
          if ( (err_x / (s_wxy * (points.size() -1))) <= 0.){
            std::cout << "\nsqrt content negative!\n" << std::endl;  
          }
          else{
            std::cout << "pt_size,errx,sqrt_arg " << (points.size()-1) << " , " << err_x << " , " << (err_x / (s_wxy * (points.size() -1))) << std::endl;
            std::cout << "x,y,z" << x << " , " << y << " , " << z << std::endl;
          }
          */
          //std::cout << err(0,0) << " , " << err(2,1) << " , " << err(2,2) << std::endl;
      }   
      
//      err(0,0) = ( err_x / (s_wxy * (points.size() -1) )) / 1.;
//      err(1,1) = ( err_y / (s_wxy * (points.size() -1) )) / 1.;
//      err(2,2) = ( err_z / (s_wz  * (points.size() -1) )) / 1.;
//      err(0,0) = ( err_x / (s_wxy)) / 1.;
//      err(1,1) = ( err_y / (s_wxy)) / 1.;
//      err(2,2) = ( err_z / (s_wz )) / 1.;
     err(0,0) = ( err_x / (s_wxy * (ndof_xy-1) ) );
     err(1,1) = ( err_y / (s_wxy * (ndof_xy-1) ) );
     err(2,2) = ( err_z / (s_wz  * (ndof_z -1) ) );

      if (err(0,0) < 1e-8){
        err(0,0) = 1e-8;
      }
      if (err(1,1) < 1e-8){
        err(1,1) = 1e-8;
      }
      if (err(2,2) < 1e-8){
        err(2,2) = 1e-8;
      }

      //std::cout << err(0,0) << " , " << err(1,1) << " , " << err(2,2) << std::endl;
      
      for (const auto& p : points){
          //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
        //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
//            TrajectoryStateClosestToPoint traj = (*iclus)[i].trajectoryStateClosestToPoint(GlobalPoint(old_x, old_y, old_z));
//            if (! traj.isValid() ) continue;
//            std::pair<GlobalPoint, GlobalPoint> p(traj.position(), GlobalPoint(0,0,0));
//        wxy = p.second.x();
//        wz = p.second.z();
//        dist = std::pow(p.first.x() - x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0),2) );
//        dist += std::pow(p.first.y() - y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1),2) );
//        dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2),2) ); 
//        dist =  std::pow(p.first.x() - x, 2) / ( std::pow(wxy, 2)  );
//        dist += std::pow(p.first.y() - y, 2) / ( std::pow(wxy, 2) );
//        dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) ); 
        dist =  std::pow(p.first.x() - x, 2) / ( (err(0,0), 2)  );
        dist += std::pow(p.first.y() - y, 2) / ( (err(1,1), 2) );
        dist += std::pow(p.first.z() - z, 2) / ( (err(2,2), 2) ); 
        chi2 += dist;
        //ndof ++;
        
      }
    TransientVertex v(GlobalPoint(x,y,z), err, *iclus, chi2, ndof_z);
    return v;
}




TransientVertex weightedMeanManyIter_IP(const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points, std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus){
     float x=0, y=0, z=0, s_wxy=0, s_wz=0, wxy=0, wz=0, chi2=0;
     unsigned int ndof = 0;
     AlgebraicSymMatrix33 err;
     err(0,0) = 2 * 2;
     err(1,1) = 2 * 2;
     err(2,2) = 20 * 20; // error is 20 cm, so cov -> is 20 ^ 2

     for (const auto& p : points){ 
            //std::vector<reco::TransientTrack>::const_iterator itrack = iclus.begin(); itrack!= iclus.end(); itrack++) {
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
            //GlobalError e(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(),  itrack.track().dzError());

            wxy = p.second.x();
            wxy = wxy <=  1e-4 ? 1e+8 : 1. / std::pow(wxy,2);

            x += p.first.x() * wxy;
            y += p.first.y() * wxy;
            s_wxy += wxy;

            wz = p.second.z();
            wz = wz <=  1e-4 ? 1e+8 : 1. / std::pow(wz,2);
            z += p.first.z() * wz;
            s_wz += wz;
     }

     //auto err = v.positionError().matrix();
     
     if ( s_wxy == 0. || s_wz == 0. ){
        std::cout << "Vertex fitting failed at beginning" << std::endl;
        return TransientVertex(GlobalPoint(0,0,0), err, *iclus, 0, 0);
     }
      x /= s_wxy;     
      y /= s_wxy;     
      z /= s_wz;
      float err_x=0, err_y=0, err_z=0;
      float old_x = x, old_y = y, old_z = z;
      float dist = 0; 
      //std::cout << "begin for loop \n\n";
      for (unsigned int iter=0; iter < 10; iter++){

          s_wxy = 0; s_wz = 0;
          old_x = x, old_y = y, old_z = z;
          x=0; y=0; z=0;

          //for (const auto& p : points){ 
          for (unsigned int i=0; i < (unsigned int) iclus->size(); i++){ 
            //dist = IPTools::absoluteImpactParameter3D((*iclus)[i],(reco::Vertex) TransientVertex(GlobalPoint(old_x, old_y, old_z), err, *iclus, 0, 0)).second();
            TrajectoryStateClosestToPoint traj = (*iclus)[i].trajectoryStateClosestToPoint(GlobalPoint(old_x, old_y, old_z));
            if (! traj.isValid() ) continue;
            std::pair<GlobalPoint, GlobalPoint> p(traj.position(), GlobalPoint(0,0,0));
            dist =  std::pow(p.first.x() - old_x, 2); 
            dist += std::pow(p.first.y() - old_y, 2); 
            dist += std::pow(p.first.z() - old_z, 2);       
            dist = dist <= 1e-4 ? 1e+8 : 1. / std::pow(dist,2);

            x += p.first.x() * dist;
            y += p.first.y() * dist;
            z += p.first.z() * dist;

            s_wxy += dist ;
            s_wz  += dist ;
          }
          if ( s_wxy == 0. || s_wz == 0. ){
            std::cout << "Vertex fitting failed at iteration " << iter << std::endl;
            return TransientVertex(GlobalPoint(0,0,0), err, *iclus, 0, 0);
          }
          x /= s_wxy;     
          y /= s_wxy;     
          z /= s_wz;     

          err_x=0.; err_y=0.; err_z=0.;

          //for (const auto& p : points){ 
          for (unsigned int i=0; i < (unsigned int) iclus->size(); i++){ 
            TrajectoryStateClosestToPoint traj = (*iclus)[i].trajectoryStateClosestToPoint(GlobalPoint(old_x, old_y, old_z));
            if (! traj.isValid() ) continue;
            std::pair<GlobalPoint, GlobalPoint> p(traj.position(), GlobalPoint(0,0,0));
            dist =  std::pow(p.first.x() - old_x, 2); 
            dist += std::pow(p.first.y() - old_y, 2);
            dist += std::pow(p.first.z() - old_z, 2);
            dist = dist <= 1e-4 ? 1e+8 : 1. / std::pow(dist,2);
            //dist = dist == 0. ? 10000. : 1. / std::pow(dist,2);
            
            err_x += dist * std::pow(p.first.x() - x, 2);
            err_y += dist * std::pow(p.first.y() - y, 2);
            err_z += dist * std::pow(p.first.z() - z, 2);
          }
          /*
          if (err_x <= 0 || err_y <= 0 || err_z <= 0) std::cout << "\none err negative\n" << std::endl;
          if (s_wxy <= 0 || s_wz <= 0) std::cout << "\none sum negative\n" << std::endl;
          if ((points.size() - 1) <= 0.5) std::cout << "\npoints size <= 0.5\n" << std::endl;
          if ( (err_x / (s_wxy * (points.size() -1))) <= 0.){
            std::cout << "\nsqrt content negative!\n" << std::endl;  
          }
          else{
            std::cout << "pt_size,errx,sqrt_arg " << (points.size()-1) << " , " << err_x << " , " << (err_x / (s_wxy * (points.size() -1))) << std::endl;
            std::cout << "x,y,z" << x << " , " << y << " , " << z << std::endl;
          }
          */
          //std::cout << err(0,0) << " , " << err(2,1) << " , " << err(2,2) << std::endl;
      }   
      
//      err(0,0) = ( err_x / (s_wxy * (points.size() -1) )) / 1.;
//      err(1,1) = ( err_y / (s_wxy * (points.size() -1) )) / 1.;
//      err(2,2) = ( err_z / (s_wz  * (points.size() -1) )) / 1.;
      err(0,0) = ( err_x / (s_wxy)) / 1.;
      err(1,1) = ( err_y / (s_wxy)) / 1.;
      err(2,2) = ( err_z / (s_wz )) / 1.;
      if (err(0,0) < 1e-8){
        err(0,0) = 1e-8;
      }
      if (err(1,1) < 1e-8){
        err(1,1) = 1e-8;
      }
      if (err(2,2) < 1e-8){
        err(2,2) = 1e-8;
      }

      //std::cout << err(0,0) << " , " << err(1,1) << " , " << err(2,2) << std::endl;
      
      //for (const auto& p : points){
          for (unsigned int i=0; i < (unsigned int) iclus->size(); i++){ 
          //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
        //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
            TrajectoryStateClosestToPoint traj = (*iclus)[i].trajectoryStateClosestToPoint(GlobalPoint(old_x, old_y, old_z));
            if (! traj.isValid() ) continue;
            std::pair<GlobalPoint, GlobalPoint> p(traj.position(), GlobalPoint(0,0,0));
//        wxy = p.second.x();
//        wz = p.second.z();
//        dist = std::pow(p.first.x() - x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0),2) );
//        dist += std::pow(p.first.y() - y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1),2) );
//        dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2),2) ); 
//        dist =  std::pow(p.first.x() - x, 2) / ( std::pow(wxy, 2)  );
//        dist += std::pow(p.first.y() - y, 2) / ( std::pow(wxy, 2) );
//        dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) ); 
        dist =  std::pow(p.first.x() - x, 2) / ( (err(0,0), 2)  );
        dist += std::pow(p.first.y() - y, 2) / ( (err(1,1), 2) );
        dist += std::pow(p.first.z() - z, 2) / ( (err(2,2), 2) ); 
        chi2 += dist;
        ndof ++;
        
      }
    TransientVertex v(GlobalPoint(x,y,z), err, *iclus, chi2, ndof);
    return v;
}


PrimaryVertexProducerDumbFitter::~PrimaryVertexProducerDumbFitter() {
  if (theTrackFilter)
    delete theTrackFilter;
  if (theTrackClusterizer)
    delete theTrackClusterizer;
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    if (algorithm->fitter)
      delete algorithm->fitter;
    if (algorithm->vertexSelector)
      delete algorithm->vertexSelector;
  }
}

void PrimaryVertexProducerDumbFitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get the BeamSpot, it will always be needed, even when not used as a constraint
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from EventSetup";
  }

  bool validBS = true;
  VertexState beamVertexState(beamSpot);
  if ((beamVertexState.error().cxx() <= 0.) || (beamVertexState.error().cyy() <= 0.) ||
      (beamVertexState.error().czz() <= 0.)) {
    validBS = false;
    edm::LogError("UnusableBeamSpot") << "Beamspot with invalid errors " << beamVertexState.error().matrix();
  }

  //if this is a recovery iteration, check if we already have a valid PV
  if (fRecoveryIteration) {
    auto const& oldVertices = iEvent.get(recoveryVtxToken);
    //look for the first valid (not-BeamSpot) vertex
    for (auto const& old : oldVertices) {
      if (!(old.isFake())) {
        //found a valid vertex, write the first one to the collection and return
        //otherwise continue with regular vertexing procedure
        auto result = std::make_unique<reco::VertexCollection>();
        result->push_back(old);
        iEvent.put(std::move(result), algorithms.begin()->label);
        return;
      }
    }
  }

  // get RECO tracks from the event
  // `tks` can be used as a ptr to a reco::TrackCollection
  edm::Handle<reco::TrackCollection> tks;
  iEvent.getByToken(trkToken, tks);

  // interface RECO tracks to vertex reconstruction
  const auto& theB = &iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> t_tks;

  if (f4D) {
    edm::Handle<edm::ValueMap<float> > trackTimesH;
    edm::Handle<edm::ValueMap<float> > trackTimeResosH;
    iEvent.getByToken(trkTimesToken, trackTimesH);
    iEvent.getByToken(trkTimeResosToken, trackTimeResosH);
    t_tks = (*theB).build(tks, beamSpot, *(trackTimesH.product()), *(trackTimeResosH.product()));
  } else {
    t_tks = (*theB).build(tks, beamSpot);
  }
  if (fVerbose) {
    std::cout << "RecoVertex/PrimaryVertexProducer"
              << "Found: " << t_tks.size() << " reconstructed tracks"
              << "\n";
  }

  // select tracks
  std::vector<reco::TransientTrack>&& seltks = theTrackFilter->select(t_tks);
//  if (seltks.size() > 4096) {
//    std::cout << "RecoVertex/PrimaryVertexProducer"
//              << "Found: " << seltks.size() << " reconstructed tracks"
//              << "\n";
//    }

  // clusterize tracks in Z
  /*
  std::sort(seltks.begin(), seltks.end(), 
    [](const reco::TransientTrack & a, const reco::TransientTrack & b) -> bool
    { 
        return (a.stateAtBeamLine().trackStateAtPCA()).position().z() > (b.stateAtBeamLine().trackStateAtPCA()).position().z(); 
    });
  std::vector<std::vector<reco::TransientTrack> > clusters;
  for (unsigned int block = 0; block < (unsigned int) std::floor(seltks.size()/512) ; block ++){
    std::vector<reco::TransientTrack> subtracks;
    //std::copy(seltks.begin() + (size_t)( block * 512 ), seltks.begin() + std::min(seltks.size()-1, (size_t)((block+1)*512)), subtracks.begin()); 
    for (unsigned int i = (block * 512); i < (unsigned int) std::min(seltks.size(),(size_t) (block+1)*512); i++) subtracks.push_back(seltks[i]);

    std::vector<std::vector<reco::TransientTrack> >&& subclusters = theTrackClusterizer->clusterize(subtracks);
    for (auto const& cluster: subclusters)  clusters.push_back(cluster); 
  }
  */
  std::vector<std::vector<reco::TransientTrack> >&& clusters = theTrackClusterizer->clusterize(seltks);

  if (fVerbose) {
    std::cout << " clustering returned  " << clusters.size() << " clusters  from " << seltks.size()
              << " selected tracks" << std::endl;
  }

  // vertex fits
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    auto result = std::make_unique<reco::VertexCollection>();
    reco::VertexCollection& vColl = (*result);

    std::vector<TransientVertex> pvs;
    for (std::vector<std::vector<reco::TransientTrack> >::const_iterator iclus = clusters.begin();
         iclus != clusters.end();
         iclus++) {
         /*
         if (iclus->size() <= 1){
             std::cout << "Cluster size <= 2, not using it" << std::endl;     
             continue;
         }
         //if (validBS && !validBS) std::cout << "ciao" << std::endl;
         std::vector<std::pair<GlobalPoint, GlobalPoint>> points;
         //std::vector<GlobalPoint> errors;
         if (algorithm->useBeamConstraint && validBS && (iclus->size() > 1)) {
             std::cout << "Using BS constraint" << std::endl;
             for (const auto& itrack : *iclus){ 
                    GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
                    GlobalPoint err(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.track().dzError());
                    std::pair<GlobalPoint, GlobalPoint> p2(p, err);
                    points.push_back(p2);
             }

            TransientVertex v = weightedMean(points, iclus);
            if ((v.positionError().matrix())(2,2) != (20*20)) pvs.push_back(v);
         } 
         else if (!(algorithm->useBeamConstraint) && (iclus->size() > 1)) {
             std::cout << "Not using BS constraint" << std::endl;
            for (const auto& itrack : *iclus){ 
                    GlobalPoint p = itrack.impactPointState().globalPosition();
                    GlobalPoint err(itrack.track().dxyError(), itrack.track().dxyError(), itrack.track().dzError());
                    std::pair<GlobalPoint, GlobalPoint> p2(p, err);
                    points.push_back(p2);
            }

            TransientVertex v = weightedMean(points, iclus);
            if ((v.positionError().matrix())(2,2) != (20*20)) pvs.push_back(v);

         }

    */    


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
      if (algorithm->useBeamConstraint && validBS && (iclus->size() > 1)) {
        v = algorithm->fitter->vertex(*iclus, beamSpot);
      } else if (!(algorithm->useBeamConstraint) && (iclus->size() > 1)) {
        v = algorithm->fitter->vertex(*iclus);
      }  // else: no fit ==> v.isValid()=False

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

      if (v.isValid() && (v.degreesOfFreedom() >= algorithm->minNdof) &&
          (!validBS || (*(algorithm->vertexSelector))(v, beamVertexState)))
        pvs.push_back(v);


    }  // end of cluster loop

    if (fVerbose) {
      std::cout << "PrimaryVertexProducerAlgorithm::vertices  candidates =" << pvs.size() << std::endl;
    }

    if (clusters.size() > 2 && clusters.size() > 2 * pvs.size())
      edm::LogWarning("PrimaryVertexProducer")
          << "more than half of candidate vertices lost " << pvs.size() << ' ' << clusters.size();

    if (pvs.empty() && seltks.size() > 5)
      edm::LogWarning("PrimaryVertexProducer")
          << "no vertex found with " << seltks.size() << " tracks and " << clusters.size() << " vertex-candidates";

    // sort vertices by pt**2  vertex (aka signal vertex tagging)
    if (pvs.size() > 1) {
      sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
    }

    // convert transient vertices returned by the theAlgo to (reco) vertices
    for (std::vector<TransientVertex>::const_iterator iv = pvs.begin(); iv != pvs.end(); iv++) {
      reco::Vertex v = *iv;
      vColl.push_back(v);
    }

    if (vColl.empty()) {
      GlobalError bse(beamSpot.rotatedCovariance3D());
      if ((bse.cxx() <= 0.) || (bse.cyy() <= 0.) || (bse.czz() <= 0.)) {
        AlgebraicSymMatrix33 we;
        we(0, 0) = 10000;
        we(1, 1) = 10000;
        we(2, 2) = 10000;
        vColl.push_back(reco::Vertex(beamSpot.position(), we, 0., 0., 0));
        if (fVerbose) {
          std::cout << "RecoVertex/PrimaryVertexProducer: "
                    << "Beamspot with invalid errors " << bse.matrix() << std::endl;
          std::cout << "Will put Vertex derived from dummy-fake BeamSpot into Event.\n";
        }
      } else {
        vColl.push_back(reco::Vertex(beamSpot.position(), beamSpot.rotatedCovariance3D(), 0., 0., 0));
        if (fVerbose) {
          std::cout << "RecoVertex/PrimaryVertexProducer: "
                    << " will put Vertex derived from BeamSpot into Event.\n";
        }
      }
    }

    if (fVerbose) {
      int ivtx = 0;
      for (reco::VertexCollection::const_iterator v = vColl.begin(); v != vColl.end(); ++v) {
        std::cout << "recvtx " << ivtx++ << "#trk " << std::setw(3) << v->tracksSize() << " chi2 " << std::setw(4)
                  << v->chi2() << " ndof " << std::setw(3) << v->ndof() << " x " << std::setw(6) << v->position().x()
                  << " dx " << std::setw(6) << v->xError() << " y " << std::setw(6) << v->position().y() << " dy "
                  << std::setw(6) << v->yError() << " z " << std::setw(6) << v->position().z() << " dz " << std::setw(6)
                  << v->zError();
        if (f4D) {
          std::cout << " t " << std::setw(6) << v->t() << " dt " << std::setw(6) << v->tError();
        }
        std::cout << std::endl;
      }
    }
    /*
      int ivtx = 0;
        std::cout << "recvtx,#trk,chi2,ndof,x,dx,y,dy,z,dz" << std::endl;
      for (reco::VertexCollection::const_iterator v = vColl.begin(); v != vColl.end(); ++v) {
        std::cout << ivtx++ << "," << v->tracksSize() << "," << v->chi2() << ","  << v->ndof() << ","  << v->position().x()
                  << ","  << v->xError() << ","  << v->position().y() << ","
                   << v->yError() << ","  << v->position().z() << "," 
                  << v->zError();
        std::cout << std::endl;

      }
    */
    iEvent.put(std::move(result), algorithm->label);
  }
}


void PrimaryVertexProducerDumbFitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // offlinePrimaryVertices
  edm::ParameterSetDescription desc;
  {
    edm::ParameterSetDescription vpsd1;
    vpsd1.add<double>("maxDistanceToBeam", 1.0);
    vpsd1.add<std::string>("algorithm", "AdaptiveVertexFitter");
    vpsd1.add<bool>("useBeamConstraint", false);
    vpsd1.add<std::string>("label", "");
    vpsd1.add<double>("chi2cutoff", 2.5);
    vpsd1.add<double>("minNdof", 0.0);
    std::vector<edm::ParameterSet> temp1;
    temp1.reserve(2);
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", false);
      temp2.addParameter<std::string>("label", "");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("minNdof", 0.0);
      temp1.push_back(temp2);
    }
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", true);
      temp2.addParameter<std::string>("label", "WithBS");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("minNdof", 2.0);
      temp1.push_back(temp2);
    }
    desc.addVPSet("vertexCollections", vpsd1, temp1);
  }
  desc.addUntracked<bool>("verbose", false);
  {
    edm::ParameterSetDescription psd0;
    TrackFilterForPVFinding::fillPSetDescription(psd0);
    psd0.add<int>("numTracksThreshold", 0);  // HI only
    desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
  }
  desc.add<edm::InputTag>("beamSpotLabel", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("TrackLabel", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("TrackTimeResosLabel", edm::InputTag("dummy_default"));  // 4D only
  desc.add<edm::InputTag>("TrackTimesLabel", edm::InputTag("dummy_default"));      // 4D only

  {
    edm::ParameterSetDescription psd0;
    {
      edm::ParameterSetDescription psd1;
      DAClusterizerInZT_vect::fillPSetDescription(psd1);
      psd0.add<edm::ParameterSetDescription>("TkDAClusParameters", psd1);

      edm::ParameterSetDescription psd2;
      GapClusterizerInZ::fillPSetDescription(psd2);
      psd0.add<edm::ParameterSetDescription>("TkGapClusParameters", psd2);

      edm::ParameterSetDescription psd3;
      DAClusterizerInZSubCluster_vect::fillPSetDescription(psd3);
      psd0.add<edm::ParameterSetDescription>("TkDASubClusParameters", psd3);
    }
    psd0.add<std::string>("algorithm", "DA_vect");
    desc.add<edm::ParameterSetDescription>("TkClusParameters", psd0);
  }

  desc.add<bool>("isRecoveryIteration", false);
  desc.add<edm::InputTag>("recoveryVtxCollection", {""});

  descriptions.add("primaryVertexProducerDumbFitter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrimaryVertexProducerDumbFitter);
