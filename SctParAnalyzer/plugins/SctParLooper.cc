// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1F.h"
#include <cmath>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "ScoutingParking/SctParAnalyzer/interface/trackPairDataFormat.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"



//Header
class SctParLooper : public edm::one::EDAnalyzer<> {
  public:
    explicit SctParLooper(const edm::ParameterSet&);
    ~SctParLooper() override;
    void beginJob() override;
    void endJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:
    triggerExpression::Data triggerCache_;
    const std::vector<std::string> vtriggerSelection_;
    
    edm::EDGetToken algToken_;
    std::vector<triggerExpression::Evaluator*> vtriggerSelector_;
    std::shared_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
    std::vector<std::string> l1Seeds_;
    TString l1Names[100] = {""};
    Bool_t l1Result[100] = {false};
    //Scouting
    edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> pvTokenScouting_;
    edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> svTokenScouting_;
    edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muTokenScouting_;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTransientTrackBuilderToken_;
    //PAT
    edm::EDGetTokenT<std::vector<pat::Muon>> PATMuonCollection_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> PATVertexCollection_;
    edm::EDGetTokenT<std::vector<reco::Track>> muon_TrackCollection_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotCollection_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertexCollection_;
    //GEN
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;

    Bool_t hltResult[100] = {false}; 
    TFile* fout;
    TTree* tout;
    unsigned int run, lumi, evtn;
    bool passL1, passHLT;
    float sct_PV_x, sct_PV_y, sct_PV_z;
    Float_t BeamSpot_x0;
    Float_t BeamSpot_y0;
    Float_t BeamSpot_z0;
    Float_t BeamSpot_BeamWidthX;
    Float_t BeamSpot_BeamWidthY;
    //Scouting SV variables
    std::vector<unsigned int> sct_index, sct_ndof;
    std::vector<float> sct_x, sct_y, sct_z;
    std::vector<float> sct_xe, sct_ye, sct_ze;
    std::vector<float> sct_chi2, sct_prob, sct_SV_chi2Ndof;
    std::vector<float> sct_lxy, sct_l3d;
    std::vector<float> sct_mindx, sct_mindy, sct_mindz, sct_mindxy, sct_mind3d;
    std::vector<float> sct_maxdx, sct_maxdy, sct_maxdz, sct_maxdxy, sct_maxd3d;
    std::vector<int> sct_SV_selected;
    std::vector<bool> sct_onModule, sct_onModuleWithinUnc;
    std::vector<float> sct_closestDet_x, sct_closestDet_y, sct_closestDet_z;
    std::vector<float> sct_minDistanceFromDet, sct_minDistanceFromDet_x, sct_minDistanceFromDet_y, sct_minDistanceFromDet_z;
    //Scouting muon variables
    std::vector<std::vector<int>> sct_vtxIdxs;
    std::vector<unsigned int> sct_saHits, sct_saMatchedStats;
    std::vector<unsigned int> sct_muHits, sct_muChambs, sct_muCSCDT, sct_muMatch, sct_muMatchedStats, sct_muExpMatchedStats, sct_muMatchedRPC;
    std::vector<unsigned int> sct_pixHits, sct_stripHits;
    std::vector<unsigned int> sct_pixLayers, sct_trkLayers;
    std::vector<int> sct_bestAssocSVIdx, sct_bestAssocSVOverlapIdx;
    std::vector<int> sct_ch;
    std::vector<bool> sct_isGlobal, sct_isTracker;
    std::vector<float> sct_pt, sct_eta, sct_phi;
    std::vector<float> sct_mu_chi2Ndof;
    std::vector<float> sct_ecalIso, sct_hcalIso, sct_trackIso;
    std::vector<float> sct_ecalRelIso, sct_hcalRelIso, sct_trackRelIso;
    std::vector<float> sct_dxy, sct_dxye, sct_dz, sct_dze;
    std::vector<float> sct_dxysig, sct_dzsig;
    std::vector<float> sct_phiCorr, sct_dxyCorr;
    std::vector<float> sct_PFIsoChg0p3, sct_PFIsoChg0p4, sct_PFIsoAll0p3, sct_PFIsoAll0p4;
    std::vector<float> sct_PFRelIsoChg0p3, sct_PFRelIsoChg0p4, sct_PFRelIsoAll0p3, sct_PFRelIsoAll0p4;
    std::vector<float> sct_mindrPF0p3, sct_mindrPF0p4;
    std::vector<float> sct_mindr, sct_maxdr;
    std::vector<float> sct_invMass;
    std::vector<float> sct_mindrJet, sct_mindphiJet, sct_mindetaJet;
    std::vector<TLorentzVector> sct_vec;
    std::vector<bool> sct_mu_selected;

    std::vector<int> sct_nhitsbeforesv, sct_ncompatible, sct_ncompatibletotal, sct_nexpectedhits, sct_nexpectedhitsmultiple, sct_nexpectedhitsmultipletotal, sct_nexpectedhitstotal, nsct_muons;
    //PAT muon variables
    std::vector<float> PAT_Muon_pt, PAT_Muon_phi, PAT_Muon_eta, PAT_vx, PAT_vy, PAT_lxy, PAT_invMass; 
    std::vector<int> nPAT_muons, PAT_pair_index, PAT_selected_index;
    std::vector<bool> PAT_muons_selected;
    //std::vector<std::vector<bool>> PAT_muons_selected;

    //GEN particle variables
    std::vector<float> gen_pt, gen_eta, gen_phi, gen_m;
    std::vector<float> gen_vx, gen_vy, gen_vz, gen_lxy;
    std::vector<int> gen_status, gen_index, gen_pdgId, gen_motherIndex, gen_motherPdgId, nGenMuons, gen_matchIndex;
    std::vector<float> gen_ct;
    std::vector<float> gen_selected;
  };

//Constructor
SctParLooper::SctParLooper(const edm::ParameterSet& iConfig) :
  triggerCache_{triggerExpression::Data(iConfig.getParameterSet("triggerConfiguration"), consumesCollector())},
  vtriggerSelection_{iConfig.getParameter<vector<string>>("triggerSelection")},
  algToken_{consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("AlgInputTag"))},

  pvTokenScouting_{consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("hltScoutingPrimaryVertexPacker_primaryVtx"))},
  svTokenScouting_{consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("hltScoutingMuonPacker_displacedVtx"))},
  muTokenScouting_{consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("Scoutingmuons"))},
  theTransientTrackBuilderToken_{esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))},
  PATMuonCollection_{consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("PATMuonCollection"))},
  PATVertexCollection_{consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("PATVertexCollection"))},
  muon_TrackCollection_{consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("muon_TrackCollection"))},
  beamSpotCollection_{consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"))},
  primaryVertexCollection_{consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertexCollection"))},
  genToken_{consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticlesInputTag"))}

  {
    vtriggerSelector_.reserve(vtriggerSelection_.size());
    for (auto const& vt : vtriggerSelection_)
      vtriggerSelector_.push_back(triggerExpression::parse(vt));
    l1GtUtils_ = std::make_shared<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), l1t::UseEventSetupIn::RunAndEvent);
    l1Seeds_ = iConfig.getParameter<std::vector<std::string>>("l1Seeds");
    for (unsigned int i = 0; i < l1Seeds_.size(); i++) {
      const auto& l1seed(l1Seeds_.at(i));
      l1Names[i] = TString(l1seed);
    }
  }

//Destructor
SctParLooper::~SctParLooper() = default;

void SctParLooper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  evtn = iEvent.id().event();
  passL1 = false;
  passHLT = false;
  const auto& theTransientTrackBuilder = iSetup.getData(theTransientTrackBuilderToken_);
  edm::Handle<std::vector<Run3ScoutingVertex>> ScoutingprimaryVertices;
  edm::Handle<std::vector<Run3ScoutingVertex>> ScoutingdisplacedVertices;
  edm::Handle<std::vector<Run3ScoutingMuon>> Scoutingmuons;
  edm::Handle<std::vector<pat::Muon>> PATMuon;
  edm::Handle<std::vector<reco::Vertex>> PATVertex;
  edm::Handle<reco::BeamSpot> beamSpot;
  edm::Handle<std::vector<reco::Vertex>> primaryvertices;
  edm::Handle<std::vector<reco::Track>> muon_Track;
  edm::Handle<std::vector<reco::GenParticle>> GenParts;

  iEvent.getByToken(pvTokenScouting_, ScoutingprimaryVertices);
  iEvent.getByToken(svTokenScouting_, ScoutingdisplacedVertices);
  iEvent.getByToken(muTokenScouting_, Scoutingmuons);
  iEvent.getByToken(PATMuonCollection_, PATMuon);
  iEvent.getByToken(PATVertexCollection_, PATVertex);
  iEvent.getByToken(muon_TrackCollection_, muon_Track);
  iEvent.getByToken(beamSpotCollection_, beamSpot);
  iEvent.getByToken(primaryVertexCollection_, primaryvertices);
  iEvent.getByToken(genToken_, GenParts);
  

  reco::BeamSpot beamSpotObject = *beamSpot;
  BeamSpot_x0 = beamSpotObject.x0();
  BeamSpot_y0 = beamSpotObject.y0();
  BeamSpot_z0 = beamSpotObject.z0();
  BeamSpot_BeamWidthX = beamSpotObject.BeamWidthX();
  BeamSpot_BeamWidthY = beamSpotObject.BeamWidthY();
  GlobalPoint _BSpoint(beamSpotObject.x0(), beamSpotObject.y0(), beamSpotObject.z0());
  GlobalPoint _0point(0.0, 0.0, 0.0);


  //////////////////////////
  //////GEN particles///////
  /////////////////////////
  gen_pt.clear(); gen_eta.clear(); gen_phi.clear(); gen_m.clear();
  gen_vx.clear(); gen_vy.clear(); gen_vz.clear(); gen_lxy.clear();
  gen_status.clear(); gen_index.clear(); gen_pdgId.clear(); gen_motherIndex.clear(); gen_motherPdgId.clear(); gen_ct.clear(), nGenMuons.clear();

  auto GenPart = *GenParts;
  unsigned int nGenParts = GenPart.size();
  int nMuonDaughters = 0;
  bool isSelected = false;
  for (unsigned int iGen = 0; iGen < nGenParts; ++iGen) {
    auto genpart = GenPart[iGen];

    int motherIdx = -1, motherPdgId = 0;
    reco::GenParticle lastCopy=genpart; // Default value
    if (genpart.numberOfMothers() > 0) {
      motherIdx = genpart.motherRef().index();
      while (GenPart[motherIdx].pdgId() == genpart.pdgId()) {
        if (GenPart[motherIdx].numberOfMothers() == 0) break;
        lastCopy = GenPart[motherIdx];
        motherIdx = GenPart[motherIdx].motherRef().index();
      }
      motherPdgId = GenPart[motherIdx].pdgId();
    }
    
    if (genpart.pdgId()==13) {
      if (motherPdgId==9900015 ) {   // Heavy neutral lepton
        ++nMuonDaughters;
        isSelected = true;
        float vx1 = lastCopy.vx() * 10.; // mm
        float vy1 = lastCopy.vy() * 10.; // mm
        float vx0 = GenPart[motherIdx].vx() * 10.; // mm
        float vy0 = GenPart[motherIdx].vy() * 10.; // mm
        gen_matchIndex.push_back(motherIdx);
        gen_ct.push_back ( ((vx1 - vx0)*GenPart[motherIdx].px() + (vy1 - vy0)*GenPart[motherIdx].py())*GenPart[motherIdx].mass()/(GenPart[motherIdx].pt()*GenPart[motherIdx].pt()) ); // in mm
      }
    }
    gen_selected.push_back(isSelected);
    gen_pt.push_back(genpart.pt());
    gen_eta.push_back(genpart.eta());
    gen_phi.push_back(genpart.phi());
    gen_m.push_back(genpart.mass());
    gen_vx.push_back(genpart.vx());
    gen_vy.push_back(genpart.vy());
    gen_vz.push_back(genpart.vz());
    gen_lxy.push_back(TMath::Hypot(genpart.vx(), genpart.vy()));
    gen_status.push_back(genpart.status());
    gen_index.push_back(iGen);
    gen_pdgId.push_back(genpart.pdgId());
    //gen_motherIndex.push_back(motherIdx);
    //gen_motherPdgId.push_back(motherPdgId);
  }
  nGenMuons.push_back(nMuonDaughters);

  if (nMuonDaughters < 2) return;

  /////////////////////////
  //Triggers:
  /////////////////////////
  const reco::Vertex &thePrimaryVertex = (*primaryvertices)[0];
  if (!ScoutingprimaryVertices.isValid()) {
    std::cout << "Primary vertex not valid" << std::endl;	  
    return;
  }
  if (!ScoutingdisplacedVertices.isValid()) {
    std::cout << "Scouting displaced vertex not valid" << std::endl;
    return;
  }
  if (!Scoutingmuons.isValid()) {
    std::cout << "Scouting muons not valid" << std::endl;
    return;
  }
  l1GtUtils_->retrieveL1(iEvent, iSetup, algToken_);

  passL1 = false;

  for (unsigned int i = 0; i < l1Seeds_.size(); i++) {
    const auto& l1seed(l1Seeds_.at(i));
    bool l1htbit = false;
    double prescale = -1;
    l1GtUtils_->getFinalDecisionByName(l1seed, l1htbit);
    l1GtUtils_->getPrescaleByName(l1seed, prescale);
    l1Result[i] = l1htbit;
    if (l1Result[i] == 1){
      passL1 = true;
    }
  }
  //if (!passL1) {
  //  return;
  //}
  

  if (triggerCache_.setEvent(iEvent, iSetup)) {
    for (unsigned int i = 0; i < vtriggerSelector_.size(); i++) {
      auto& vts(vtriggerSelector_.at(i));
      bool result = false;
      if (vts) {
        if (triggerCache_.configurationUpdated())
          vts->init(triggerCache_);
        result = (*vts)(triggerCache_);
      }
      hltResult[i] = result;
      if (result)
        passHLT = true;
    }
  }
  //if (!passHLT) {
  //  return;
  //}

  

  /////////////////////////
  //Scouting muons:
  /////////////////////////
  //bool relaxedSVSel = false;
  //float sfSVsel = (relaxedSVSel) ? 5.0 : 1.0;
  //float maxXerr=0.05*sfSVsel, maxYerr=0.05*sfSVsel, maxZerr=0.10*sfSVsel, maxChi2=3.0*sfSVsel;
  
  sct_PV_x = 0;
  sct_PV_y = 0;
  sct_PV_z = 0;

  if (ScoutingprimaryVertices.isValid() && !ScoutingprimaryVertices->empty() && (*ScoutingprimaryVertices)[0].isValidVtx()) {
    sct_PV_x = (*ScoutingprimaryVertices)[0].x();
    sct_PV_y = (*ScoutingprimaryVertices)[0].y();
    sct_PV_z = (*ScoutingprimaryVertices)[0].z();
  }



  //SV selection
  sct_index.clear(); sct_ndof.clear();
  sct_x.clear(); sct_y.clear(); sct_z.clear();
  sct_xe.clear(); sct_ye.clear(); sct_ze.clear();
  sct_chi2.clear(); sct_prob.clear(); sct_SV_chi2Ndof.clear();
  sct_lxy.clear(); sct_l3d.clear();
  sct_mindx.clear(); sct_mindy.clear(); sct_mindz.clear(); sct_mindxy.clear(); sct_mind3d.clear();
  sct_maxdx.clear(); sct_maxdy.clear(); sct_maxdz.clear(); sct_maxdxy.clear(); sct_maxd3d.clear();
  sct_SV_selected.clear();
  sct_onModule.clear(); sct_onModuleWithinUnc.clear();
  sct_closestDet_x.clear(); sct_closestDet_y.clear(); sct_closestDet_z.clear();
  sct_minDistanceFromDet.clear(), sct_minDistanceFromDet_x.clear(), sct_minDistanceFromDet_y.clear(), sct_minDistanceFromDet_z.clear();
  const auto& sct_svCollection = *ScoutingdisplacedVertices;
  unsigned int nSVs = sct_svCollection.size();

  //if (nSVs < 1){
  //  return;
  //}
  //for (unsigned int iSV = 0; iSV < nSVs; ++iSV) {
  //  const auto& sv = sct_svCollection[iSV];
  //  if (!sv.isValidVtx()) continue;	
  //  sct_index.push_back(iSV);
  //  sct_ndof.push_back(sv.ndof());
  //  sct_x.push_back(sv.x());
  //  sct_y.push_back(sv.y());
  //  sct_z.push_back(sv.z());
  //  sct_xe.push_back(sv.xError());
  //  sct_ye.push_back(sv.yError());
  //  sct_ze.push_back(sv.zError());
  //  sct_chi2.push_back(sv.chi2());
  //  sct_prob.push_back(TMath::Prob(sv.chi2(), sv.ndof()));
  //  sct_SV_chi2Ndof.push_back(sv.chi2()/sv.ndof());
  //  float dx = sv.x() - sct_PV_x;
  //  float dy = sv.y() - sct_PV_y;
  //  float dz = sv.z() - sct_PV_z;
  //  float sct_lxy_int = std::sqrt(dx * dx + dy * dy);
  //  float sct_l3d_int = std::sqrt(sct_lxy_int * sct_lxy_int + dz * dz);
//
  //  sct_lxy.push_back(sct_lxy_int);
  //  sct_l3d.push_back(sct_l3d_int);
//
  //  //bool passSel = (sv.xError() < maxXerr && sv.yError() < maxYerr && sv.zError() < maxZerr && sv.chi2() / sv.ndof() < maxChi2);
  //  //sct_SV_selected.push_back(passSel);
  //}
  
  //Muon selection  
  sct_vtxIdxs.clear();
  sct_saHits.clear(); sct_saMatchedStats.clear();
  sct_muHits.clear(); sct_muChambs.clear(); sct_muCSCDT.clear(); sct_muMatch.clear(); sct_muMatchedStats.clear(); sct_muExpMatchedStats.clear(); sct_muMatchedRPC.clear();
  sct_pixHits.clear(); sct_stripHits.clear();
  sct_pixLayers.clear(); sct_trkLayers.clear();
  sct_bestAssocSVIdx.clear(); sct_bestAssocSVOverlapIdx.clear();
  sct_ch.clear();
  sct_isGlobal.clear(); sct_isTracker.clear();
  sct_pt.clear(); sct_eta.clear(); sct_phi.clear();
  sct_mu_chi2Ndof.clear();
  sct_ecalIso.clear(); sct_hcalIso.clear(); sct_trackIso.clear();
  sct_ecalRelIso.clear(); sct_hcalRelIso.clear(); sct_trackRelIso.clear();
  sct_dxy.clear(); sct_dxye.clear(); sct_dz.clear(); sct_dze.clear();
  sct_dxysig.clear(); sct_dzsig.clear();
  sct_phiCorr.clear(); sct_dxyCorr.clear();
  sct_PFIsoChg0p3.clear(); sct_PFIsoAll0p3.clear();
  sct_PFRelIsoChg0p3.clear(); sct_PFRelIsoAll0p3.clear();
  sct_PFIsoChg0p4.clear(); sct_PFIsoAll0p4.clear();
  sct_PFRelIsoChg0p4.clear(); sct_PFRelIsoAll0p4.clear();
  sct_mindrPF0p3.clear(); sct_mindrPF0p4.clear();
  sct_mindr.clear(); sct_maxdr.clear();
  sct_mindrJet.clear(); sct_mindphiJet.clear(); sct_mindetaJet.clear();
  sct_vec.clear();
  sct_mu_selected.clear();
  sct_nhitsbeforesv.clear(); sct_ncompatible.clear(); sct_ncompatibletotal.clear(); sct_nexpectedhits.clear(); sct_nexpectedhitsmultiple.clear(); sct_nexpectedhitsmultipletotal.clear(); sct_nexpectedhitstotal.clear();
  sct_invMass.clear(); nsct_muons.clear();

  if (!Scoutingmuons.isValid()) {
    return;
  }
  const auto& muonCollection = *Scoutingmuons;
  unsigned int nMus = muonCollection.size();
  nsct_muons.push_back(nMus);
  std::vector<bool> sct_mu_select(nMus, false);
  for (unsigned int iMu = 0; iMu < nMus; ++iMu) {
    if (nMus < 1){
      break;
    }
    if (nSVs < 1){
      break;
    }
    const auto& mu = muonCollection[iMu];
    std::vector<int> sct_matchedAndSelVtxIdxs;
    for (unsigned int matchedVtxIdx : mu.vtxIndx()) {
      for (size_t iSV = 0; iSV < sct_index.size(); ++iSV) {
        if (matchedVtxIdx == sct_index[iSV]) {
          sct_matchedAndSelVtxIdxs.push_back(matchedVtxIdx);
        }
      }
    }
    //Pair muons if they have the same vertex index and find their invariant mass
    for (unsigned int iMu2 = iMu + 1; iMu2 < nMus; ++iMu2) {
      const auto& mu2 = muonCollection[iMu2];
      const std::vector<int> vtxIndx_1 = mu.vtxIndx();
      const std::vector<int> vtxIndx_2 = mu2.vtxIndx();
      math::PtEtaPhiMLorentzVector sct_mu1(mu.pt(), mu.eta(), mu.phi(), mu.m());
      math::PtEtaPhiMLorentzVector sct_mu2(mu2.pt(), mu2.eta(), mu2.phi(), mu2.m());
      for (const auto& commonIdx : vtxIndx_2) {
        if (std::find(vtxIndx_1.begin(), vtxIndx_1.end(), commonIdx) != vtxIndx_1.end()) {
          sct_invMass.push_back((sct_mu1 + sct_mu2).mass());
          if ((sct_mu1 + sct_mu2).mass() > 1.7 && (sct_mu1 + sct_mu2).mass() < 2.3 && mu.charge() != mu2.charge()) {  //If the dimuon system comes from Zd
            sct_mu_select[iMu] = true;
            sct_mu_select[iMu2] = true;
            for (unsigned int iSV = 0; iSV < nSVs; ++iSV) {
              const auto& sv = sct_svCollection[iSV];
              if (!sv.isValidVtx()) continue;	
              sct_index.push_back(iSV);
              sct_ndof.push_back(sv.ndof());
              sct_x.push_back(sv.x());
              sct_y.push_back(sv.y());
              sct_z.push_back(sv.z());
              sct_xe.push_back(sv.xError());
              sct_ye.push_back(sv.yError());
              sct_ze.push_back(sv.zError());
              sct_chi2.push_back(sv.chi2());
              sct_prob.push_back(TMath::Prob(sv.chi2(), sv.ndof()));
              sct_SV_chi2Ndof.push_back(sv.chi2()/sv.ndof());
              float dx = sv.x() - sct_PV_x;
              float dy = sv.y() - sct_PV_y;
              float dz = sv.z() - sct_PV_z;
              float sct_lxy_int = std::sqrt(dx * dx + dy * dy);
              float sct_l3d_int = std::sqrt(sct_lxy_int * sct_lxy_int + dz * dz);
              sct_lxy.push_back(sct_lxy_int);
              sct_l3d.push_back(sct_l3d_int);
            }      
          }
          break;
        }
      }
    }
    sct_pt.push_back(mu.pt());
    sct_eta.push_back(mu.eta());
    sct_phi.push_back(mu.phi());
    sct_ch.push_back(mu.charge());
    sct_isGlobal.push_back(mu.isGlobalMuon());
    sct_isTracker.push_back(mu.isTrackerMuon());
    sct_mu_chi2Ndof.push_back(mu.normalizedChi2());
    sct_mu_selected.push_back(sct_mu_select[iMu]);
  }
  /////////////////////////
  //PAT muons:
  /////////////////////////
  int nPatMuons      = PATMuon->size();
  PAT_Muon_pt.clear(); PAT_Muon_phi.clear(); PAT_Muon_eta.clear(); PAT_lxy.clear(); PAT_vx.clear(); PAT_vy.clear(); 
  PAT_invMass.clear(); nPAT_muons.clear(); PAT_muons_selected.clear(); PAT_pair_index.clear(); PAT_selected_index.clear();

  float minChi2 = 10000;
  nPAT_muons.push_back(nPatMuons);
  int selected_index = 0;
  int selected_pair_index = 0;
  std::vector<bool> PAT_mu_select(nPatMuons, false);
  
  if (nPatMuons > 1) {
    for (int i = 0; i < nPatMuons; i++) {
      int min_i = 99;
      int min_j = 99;
      for (int j = i+1; j < nPatMuons; j++) {
        if (i == j) continue;
        const pat::Muon& mu_i = (*PATMuon)[i];
        const pat::Muon& mu_j = (*PATMuon)[j];
        reco::TrackRef tr_i;
        reco::TrackRef tr_j;
        if (mu_i.isGlobalMuon() && mu_j.isGlobalMuon()){
          tr_i = mu_i.globalTrack();
          tr_j = mu_j.globalTrack();
        }
        else if (mu_i.isTrackerMuon() && mu_j.isTrackerMuon()){
          tr_i = mu_i.innerTrack();
          tr_j = mu_j.innerTrack();
        }
        else continue;

        trackPair testcandidate(thePrimaryVertex, beamSpotObject, theTransientTrackBuilder, tr_i, tr_j, false);
        if (!testcandidate.hasValidVertex) continue ;
        if (testcandidate.normalizedChi2 < minChi2) {
          minChi2 = testcandidate.normalizedChi2;
          min_i = i;
          min_j = j;
        }
      }

      if (min_i != 99 && min_j != 99) {
        const pat::Muon& mu_1 = (*PATMuon)[min_i];
        const pat::Muon& mu_2 = (*PATMuon)[min_j];
        reco::TrackRef tr_1;
        reco::TrackRef tr_2;
        if (mu_1.isGlobalMuon() && mu_2.isGlobalMuon()){
          tr_1 = mu_1.globalTrack();
          tr_2 = mu_2.globalTrack();
        }
        else if (mu_1.isTrackerMuon() && mu_2.isTrackerMuon()){
          tr_1 = mu_1.innerTrack();
          tr_2 = mu_2.innerTrack();
        }
        trackPair muonPair(thePrimaryVertex, beamSpotObject, theTransientTrackBuilder, tr_1, tr_2, false);
        std::vector<pat::Muon> sorted_mu;
        sorted_mu.push_back(mu_1);
        sorted_mu.push_back(mu_2);
        std::sort(std::begin(sorted_mu), std::end(sorted_mu), [&](const pat::Muon mu1, const pat::Muon mu2) {
          return mu_1.pt() > mu_2.pt();
        });
        pat::Muon leading_mu = sorted_mu.at(0);
        pat::Muon subleading_mu = sorted_mu.at(1);
        PAT_Muon_pt.push_back(leading_mu.pt());
        PAT_Muon_pt.push_back(subleading_mu.pt());
        PAT_Muon_phi.push_back(leading_mu.phi());
        PAT_Muon_phi.push_back(subleading_mu.phi());
        PAT_Muon_eta.push_back(leading_mu.eta());
        PAT_Muon_eta.push_back(subleading_mu.eta());
        float vx = muonPair.vx;
        float vy = muonPair.vy;
        float lxy = std::sqrt(vx * vx + vy * vy);
        float invMass_pat = muonPair.mass;
        if (invMass_pat > 1.7 && invMass_pat < 2.3) {  //If the dimuon system comes from Zd
          PAT_mu_select[min_i] = true;
          PAT_mu_select[min_j] = true;
          PAT_selected_index.push_back(selected_index);
          PAT_selected_index.push_back(selected_index + 1);
          PAT_pair_index.push_back(selected_pair_index);
        }
        PAT_vx.push_back(vx);
        PAT_vy.push_back(vy);
        PAT_lxy.push_back(lxy);
        PAT_invMass.push_back(invMass_pat);
        selected_index = selected_index + 2;
        ++selected_pair_index;
      }  
    PAT_muons_selected.push_back(PAT_mu_select[i]);
    }
  }


  tout->Fill();
      
}

void SctParLooper::beginJob() {
  fout = new TFile("output_ctau-100-mA-2p00-mpi-10.root", "RECREATE");
  tout = new TTree("tout","Run3ScoutingTree");

  for (size_t iL1 = 0; iL1 < l1Seeds_.size(); ++iL1) {
    tout->Branch(l1Seeds_[iL1].c_str(), &l1Result[iL1], (l1Seeds_[iL1] + "/B").c_str());
  }

  for (size_t i = 0; i < vtriggerSelector_.size(); ++i) {
    tout->Branch(vtriggerSelection_[i].c_str(), &hltResult[i], (vtriggerSelection_[i] + "/B").c_str());
  }

  tout->Branch("run", &run);
  tout->Branch("lumi", &lumi);
  tout->Branch("evtn", &evtn);
  tout->Branch("passL1", &passL1);
  tout->Branch("passHLT", &passHLT);
  tout->Branch("sct_PV_x", &sct_PV_x);
  tout->Branch("sct_PV_y", &sct_PV_y);
  tout->Branch("sct_PV_z", &sct_PV_z);

  //Scouting SV branches
  tout->Branch("sct_SV_xe", &sct_xe);
  tout->Branch("sct_SV_ye", &sct_ye);
  tout->Branch("sct_SV_ze", &sct_ze);
  tout->Branch("sct_SV_chi2", &sct_chi2);
  tout->Branch("sct_SV_prob", &sct_prob);
  tout->Branch("sct_SV_chi2Ndof", &sct_SV_chi2Ndof);
  tout->Branch("sct_SV_lxy", &sct_lxy);
  tout->Branch("sct_SV_l3d", &sct_l3d);
  tout->Branch("sct_SV_selected", &sct_SV_selected);
  
  //Scouting Muon branches

  tout->Branch("sct_Muon_pt", &sct_pt);
  tout->Branch("sct_Muon_eta", &sct_eta);
  tout->Branch("sct_Muon_phi", &sct_phi);
  tout->Branch("sct_Muon_isGlobal", &sct_isGlobal);
  tout->Branch("sct_Muon_isTracker", &sct_isTracker);
  tout->Branch("sct_invMass", &sct_invMass);
  tout->Branch("nsct_muons", &nsct_muons);
  tout->Branch("sct_Muon_selected", &sct_mu_selected);

  //PAT Muon branches
  tout->Branch("PAT_Muon_pt", &PAT_Muon_pt);
  tout->Branch("PAT_Muon_eta", &PAT_Muon_eta);
  tout->Branch("PAT_Muon_phi", &PAT_Muon_phi);
  tout->Branch("PAT_lxy", &PAT_lxy);
  tout->Branch("PAT_vx", &PAT_vx);
  tout->Branch("PAT_vy", &PAT_vy);
  tout->Branch("PAT_invMass", &PAT_invMass);
  tout->Branch("nPAT_muons", &nPAT_muons);
  tout->Branch("PAT_Muon_selected", &PAT_muons_selected);
  tout->Branch("PAT_pair_index", &PAT_pair_index);
  tout->Branch("PAT_selected_index", &PAT_selected_index);

  //GEN branches
  tout->Branch("Gen_pt", &gen_pt);
  tout->Branch("Gen_eta", &gen_eta);
  tout->Branch("Gen_phi", &gen_phi);
  tout->Branch("Gen_m", &gen_m);
  tout->Branch("Gen_vx", &gen_vx);
  tout->Branch("Gen_vy", &gen_vy);
  tout->Branch("Gen_vz", &gen_vz);
  tout->Branch("Gen_lxy", &gen_lxy);
  tout->Branch("Gen_status", &gen_status);
  tout->Branch("Gen_index", &gen_index);
  tout->Branch("Gen_pdgId", &gen_pdgId);
  tout->Branch("Gen_motherIndex", &gen_motherIndex);
  tout->Branch("Gen_motherPdgId", &gen_motherPdgId);
  tout->Branch("Gen_ct", &gen_ct);
  tout->Branch("Number_GenMuons", &nGenMuons);

}


void SctParLooper::endJob() {
  fout->Write();
  if (tout) {
    std::cout << "SctParLooper: total entries in tout = " << tout->GetEntries() << std::endl;
  } else {
    std::cout << "SctParLooper: tout is null!" << std::endl;
  }
  fout->Close();
}
DEFINE_FWK_MODULE(SctParLooper);



