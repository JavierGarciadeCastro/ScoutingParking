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
    edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muTokenScoutingNoVtx_;
    edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muTokenScoutingVtx_;
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
    bool selected_sv;
    float sct_PV_x, sct_PV_y, sct_PV_z;
    Float_t BeamSpot_x0;
    Float_t BeamSpot_y0;
    Float_t BeamSpot_z0;
    Float_t BeamSpot_BeamWidthX;
    Float_t BeamSpot_BeamWidthY;
    //Scouting SV variables
    std::vector<unsigned int> sct_ndof, nsct_SV;
    std::vector<float> sct_x, sct_y, sct_z;
    std::vector<float> sct_xe, sct_ye, sct_ze;
    std::vector<float> sct_chi2, sct_prob, sct_SV_chi2Ndof;
    std::vector<float> sct_lxy, sct_l3d;
    
    //Scouting muon variables
    std::vector<int> sct_ch_NoVtx, nsct_muons_NoVtx;
    std::vector<bool> sct_isGlobal_NoVtx, sct_isTracker_NoVtx;
    std::vector<float> sct_pt_NoVtx, sct_eta_NoVtx, sct_phi_NoVtx;
    std::vector<float> sct_mu_chi2Ndof_NoVtx;
    std::vector<std::vector<int>> sct_mu_vtxIdx_NoVtx;

    std::vector<int> sct_ch_Vtx, nsct_muons_Vtx;
    std::vector<bool> sct_isGlobal_Vtx, sct_isTracker_Vtx;
    std::vector<float> sct_pt_Vtx, sct_eta_Vtx, sct_phi_Vtx;
    std::vector<float> sct_mu_chi2Ndof_Vtx;
    std::vector<std::vector<int>> sct_mu_vtxIdx_Vtx;

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
  muTokenScoutingNoVtx_{consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("ScoutingmuonsNoVtx"))},
  muTokenScoutingVtx_{consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("ScoutingmuonsVtx"))},
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
  edm::Handle<std::vector<Run3ScoutingMuon>> ScoutingmuonsNoVtx;
  edm::Handle<std::vector<Run3ScoutingMuon>> ScoutingmuonsVtx;
  edm::Handle<std::vector<pat::Muon>> PATMuon;
  edm::Handle<std::vector<reco::Vertex>> PATVertex;
  edm::Handle<reco::BeamSpot> beamSpot;
  edm::Handle<std::vector<reco::Vertex>> primaryvertices;
  edm::Handle<std::vector<reco::Track>> muon_Track;
  edm::Handle<std::vector<reco::GenParticle>> GenParts;

  iEvent.getByToken(pvTokenScouting_, ScoutingprimaryVertices);
  iEvent.getByToken(svTokenScouting_, ScoutingdisplacedVertices);
  iEvent.getByToken(muTokenScoutingNoVtx_, ScoutingmuonsNoVtx);
  iEvent.getByToken(muTokenScoutingVtx_, ScoutingmuonsVtx);
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
    
    if (abs(genpart.pdgId())==13 && genpart.status() == 1) {
      if (abs(motherPdgId)==9900015 ) {   
        ++nMuonDaughters;
        isSelected = true;
        float vx1 = lastCopy.vx(); // mm
        float vy1 = lastCopy.vy(); // mm
        float vx0 = GenPart[motherIdx].vx(); // mm
        float vy0 = GenPart[motherIdx].vy(); // mm
        gen_matchIndex.push_back(motherIdx);
        gen_ct.push_back ( ((vx1 - vx0)*GenPart[motherIdx].px() + (vy1 - vy0)*GenPart[motherIdx].py())*GenPart[motherIdx].mass()/(GenPart[motherIdx].pt()*GenPart[motherIdx].pt()) ); // in mm
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
      }
    }

  }
  nGenMuons.push_back(nMuonDaughters);

  //if (nMuonDaughters < 2) return;

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
  if (!ScoutingmuonsNoVtx.isValid()) {
    std::cout << "Scouting muons No Vtx not valid" << std::endl;
    return;
  }

  if (!ScoutingmuonsVtx.isValid()) {
    std::cout << "Scouting muons Vtx not valid" << std::endl;
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


  ///////////////////////////////////////////////////
  ////////////// Scouting Muons ///////////////////////
  ///////////////////////////////////////////////////
  
  sct_pt_NoVtx.clear(); sct_eta_NoVtx.clear(); sct_phi_NoVtx.clear(); sct_ch_NoVtx.clear(); sct_isGlobal_NoVtx.clear(); sct_isTracker_NoVtx.clear(); sct_mu_chi2Ndof_NoVtx.clear(); sct_mu_vtxIdx_NoVtx.clear();
  sct_pt_Vtx.clear(); sct_eta_Vtx.clear(); sct_phi_Vtx.clear(); sct_ch_Vtx.clear(); sct_isGlobal_Vtx.clear(); sct_isTracker_Vtx.clear(); sct_mu_chi2Ndof_Vtx.clear(); sct_mu_vtxIdx_Vtx.clear();

  const auto& muonCollectionNoVtx = *ScoutingmuonsNoVtx;
  unsigned int nMusNoVtx = muonCollectionNoVtx.size();
  nsct_muons_NoVtx.push_back(nMusNoVtx);

  for (unsigned int iMu = 0; iMu < nMusNoVtx; ++iMu) {
    const auto& mu = muonCollectionNoVtx[iMu];
    sct_pt_NoVtx.push_back(mu.pt());
    sct_eta_NoVtx.push_back(mu.eta());
    sct_phi_NoVtx.push_back(mu.phi());
    sct_ch_NoVtx.push_back(mu.charge());
    sct_isGlobal_NoVtx.push_back(mu.isGlobalMuon());
    sct_isTracker_NoVtx.push_back(mu.isTrackerMuon());
    sct_mu_chi2Ndof_NoVtx.push_back(mu.normalizedChi2());
    sct_mu_vtxIdx_NoVtx.push_back(mu.vtxIndx());
  }

  const auto& muonCollectionVtx = *ScoutingmuonsVtx;
  unsigned int nMus_Vtx = muonCollectionVtx.size();
  nsct_muons_Vtx.push_back(nMus_Vtx);

  for (unsigned int iMu_Vtx = 0; iMu_Vtx < nMus_Vtx; ++iMu_Vtx) {
    const auto& mu_Vtx = muonCollectionVtx[iMu_Vtx];
    sct_pt_Vtx.push_back(mu_Vtx.pt());
    sct_eta_Vtx.push_back(mu_Vtx.eta());
    sct_phi_Vtx.push_back(mu_Vtx.phi());
    sct_ch_Vtx.push_back(mu_Vtx.charge());
    sct_isGlobal_Vtx.push_back(mu_Vtx.isGlobalMuon());
    sct_isTracker_Vtx.push_back(mu_Vtx.isTrackerMuon());
    sct_mu_chi2Ndof_Vtx.push_back(mu_Vtx.normalizedChi2());
    sct_mu_vtxIdx_Vtx.push_back(mu_Vtx.vtxIndx());
  }

  ///////////////////////////////////////////////////
  ////////////// Scouting SVs ///////////////////////
  ///////////////////////////////////////////////////

  sct_ndof.clear(); sct_x.clear(); sct_y.clear(); sct_z.clear(); sct_xe.clear(); sct_ye.clear(); sct_ze.clear(); 
  sct_chi2.clear(); sct_prob.clear(); sct_SV_chi2Ndof.clear(); sct_lxy.clear(); sct_l3d.clear();

  const auto& sct_svCollection = *ScoutingdisplacedVertices;
  unsigned int nSVs = sct_svCollection.size();
  nsct_SV.push_back(nSVs);

  for (unsigned int iSV = 0; iSV < nSVs; ++iSV) {
    const auto& sv = (*ScoutingdisplacedVertices)[iSV];
    //sct_invMass.push_back((sct_mu1 + sct_mu2).mass());
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
  


  /////////////////////////
  ////// PAT muons ////////
  /////////////////////////
  PAT_Muon_pt.clear(); PAT_Muon_phi.clear(); PAT_Muon_eta.clear(); PAT_lxy.clear(); PAT_vx.clear(); PAT_vy.clear(); 
  PAT_invMass.clear(); nPAT_muons.clear(); PAT_muons_selected.clear(); PAT_pair_index.clear(); PAT_selected_index.clear();

  unsigned int nPatMuons = PATMuon->size();
  nPAT_muons.push_back(nPatMuons);

  for (unsigned int iPAT = 0; iPAT < nPatMuons; ++iPAT) {
    const pat::Muon& Pat = (*PATMuon)[iPAT];
    PAT_Muon_pt.push_back(Pat.pt());
    PAT_Muon_phi.push_back(Pat.phi());
    PAT_Muon_eta.push_back(Pat.eta());
  }


  tout->Fill();
      
}

void SctParLooper::beginJob() {
  fout = new TFile("output_test.root", "RECREATE");
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
  tout->Branch("nsct_SV", &nsct_SV);
  tout->Branch("sct_SV_x", &sct_x);
  tout->Branch("sct_SV_y", &sct_y);
  tout->Branch("sct_SV_z", &sct_z);
  tout->Branch("sct_SV_xe", &sct_xe);
  tout->Branch("sct_SV_ye", &sct_ye);
  tout->Branch("sct_SV_ze", &sct_ze);
  tout->Branch("sct_SV_chi2", &sct_chi2);
  tout->Branch("sct_SV_prob", &sct_prob);
  tout->Branch("sct_SV_chi2Ndof", &sct_SV_chi2Ndof);
  tout->Branch("sct_SV_lxy", &sct_lxy);
  tout->Branch("sct_SV_l3d", &sct_l3d);
  
  //Scouting Muon branches
  tout->Branch("nsct_muons_NoVtx", &nsct_muons_NoVtx);
  tout->Branch("sct_Muon_pt_NoVtx", &sct_pt_NoVtx);
  tout->Branch("sct_Muon_eta_NoVtx", &sct_eta_NoVtx);
  tout->Branch("sct_Muon_phi_NoVtx", &sct_phi_NoVtx);
  tout->Branch("sct_Muon_isGlobal_NoVtx", &sct_isGlobal_NoVtx);
  tout->Branch("sct_Muon_isTracker_NoVtx", &sct_isTracker_NoVtx);
  tout->Branch("sct_Muon_vtxIdx_NoVtx", &sct_mu_vtxIdx_NoVtx);

  tout->Branch("nsct_muons_Vtx", &nsct_muons_Vtx);
  tout->Branch("sct_Muon_pt_Vtx", &sct_pt_Vtx);
  tout->Branch("sct_Muon_eta_Vtx", &sct_eta_Vtx);
  tout->Branch("sct_Muon_phi_Vtx", &sct_phi_Vtx);
  tout->Branch("sct_Muon_isGlobal_Vtx", &sct_isGlobal_Vtx);
  tout->Branch("sct_Muon_isTracker_Vtx", &sct_isTracker_Vtx);
  tout->Branch("sct_Muon_vtxIdx_Vtx", &sct_mu_vtxIdx_Vtx);

  //PAT Muon branches
  tout->Branch("PAT_Muon_pt", &PAT_Muon_pt);
  tout->Branch("PAT_Muon_eta", &PAT_Muon_eta);
  tout->Branch("PAT_Muon_phi", &PAT_Muon_phi);
  tout->Branch("nPAT_muons", &nPAT_muons);


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



