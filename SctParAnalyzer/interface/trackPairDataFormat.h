#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

//#include "DataFormats/Math/interface/Error.h"
//#include "DataFormats/Math/interface/Point3D.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"



#include "TLorentzVector.h"
#include "TVector3.h"

#include <TROOT.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>


struct trackPair
{

   int type = 99;  // 0:electrons 1:muons
   bool canFitVertex = false;
   bool hasValidVertex = false;
   double Lxy_PV = 0;
   double Ixy_PV = 0;
   double Lxy_BS = 0;
   double Ixy_BS = 0;
   double Lxy_0 = 0;
   double Ixy_0 = 0;
   double normalizedChi2 = 0;
   double trackDxy = 0;
   double trackIxy = 0; 
   double trackDxy_PV = 0;
   double trackIxy_PV = 0; 
   double trackDxy_0 = 0; 
   double trackIxy_0 = 0; 
   double trackDxy_BS = 0;
   double trackIxy_BS = 0;
   double etaA = 0;
   double etaB = 0;
   double leadingPt = 0;
   double subleadingPt = 0;
   double leadingEt = 0;
   double subleadingEt = 0;
   double mass = 0;   
   double ptll = 0;
   double cosAlpha = 0;
   double dPhi = 0;
   double lldPhi = 0;
   double dR = 0;
   double relisoA = 0;
   double relisoB = 0;
   double vx = 0;       
   double vy = 0;       

   int fromPVA = 0;
   int fromPVB = 0;
   int PVAssociation = 0;

   trackPair(const reco::Vertex &pv, const reco::BeamSpot &bs, const TransientTrackBuilder &theTransientTrackBuilder, const reco::TrackRef &tr_A, const reco::TrackRef &tr_B, bool isEE)
   { 
      Init(theTransientTrackBuilder, pv, bs, tr_A, tr_B, isEE); 
   };
   ~trackPair(){};

   void Init(const TransientTrackBuilder &theTransientTrackBuilder, const reco::Vertex &pv, const reco::BeamSpot &bs,const reco::TrackRef &tr_A, const reco::TrackRef &tr_B, bool isEE)
   {
      if (isEE) { type = 0; }
      else { type = 1; }

      // Get tracks:
      std::vector<reco::TransientTrack> vec_refitTracks;
      reco::TransientTrack isotransienttrackA = theTransientTrackBuilder.build(tr_A);
      reco::TransientTrack isotransienttrackB = theTransientTrackBuilder.build(tr_B);
      vec_refitTracks.push_back(isotransienttrackA); vec_refitTracks.push_back(isotransienttrackB);

      // Fit tracks:
      //AdaptiveVertexFitter  thefitterll(GeometricAnnealing(2.5));
      KalmanVertexFitter thefitterll;
      TransientVertex myVertex = thefitterll.vertex(vec_refitTracks);
      const reco::Vertex secV = myVertex;


      // If the vertex is valid get the geometric information:
      if (secV.isValid()) {

          hasValidVertex = true;

          // Define the axis along the direction of the distance is defined:
          GlobalVector axis(0,0,0);
          axis = GlobalVector(secV.x(),secV.y(),secV.z());

          // Vertex position and fit details:
          normalizedChi2 = myVertex.normalisedChiSquared();         
          vx = secV.x();
          vy = secV.y();

          // Kinematics: 
          leadingPt = (tr_A->pt()>tr_B->pt())? tr_A->pt(): tr_B->pt();
          subleadingPt = (tr_A->pt()<tr_B->pt())? tr_A->pt(): tr_B->pt();

          etaA = tr_A->eta();
          etaB = tr_B->eta();   

          // Vector angles: 
	  /* 
          TVector3 vec3A(tr_A.px(), tr_A.py(), tr_A.pz());
          TVector3 vec3B(tr_B.px(), tr_B.py(), tr_B.pz());
          TVector3 divec3 = vec3A + vec3B;
          TVector3 vtxvec3(secV.x() - pv.x(), secV.y() - pv.y(), secV.z() - pv.z());
          cosAlpha = TMath::Cos(vec3A.Angle(vec3B));
          dPhi = divec3.DeltaPhi(vtxvec3);
          dR = vec3A.DeltaR(vec3B);
	  */

          // Physical magnitudes:
          TLorentzVector la;
          TLorentzVector lb;

          if (isEE) {
             la.SetPtEtaPhiM(tr_A->pt(), tr_A->eta(), tr_A->phi(), 0.510/1000.0);
             lb.SetPtEtaPhiM(tr_B->pt(), tr_B->eta(), tr_B->phi(), 0.510/1000.0);
          } else {
             la.SetPtEtaPhiM(tr_A->pt(), tr_A->eta(), tr_A->phi(), 105.658/1000.0);
             lb.SetPtEtaPhiM(tr_B->pt(), tr_B->eta(), tr_B->phi(), 105.658/1000.0);
          }

          mass = (la + lb).M();
          ptll = (la + lb).Pt();

      } else {
         hasValidVertex = false;
      }
   } 
};
