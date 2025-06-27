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
    //edm::EDGetTokenT<std::vector<bool>> l1ResultToken_;
    //edm::EDGetTokenT<std::vector<bool>> hltResultToken_;
    //edm::EDGetTokenT<std::vector<std::string>> l1NameToken_;
    //edm::EDGetTokenT<std::vector<std::string>> hltNameToken_;
    
    edm::EDGetToken algToken_;
    std::vector<triggerExpression::Evaluator*> vtriggerSelector_;
    std::shared_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
    std::vector<std::string> l1Seeds_;
    TString l1Names[100] = {""};
    Bool_t l1Result[100] = {false};
    edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> pvToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> svToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muToken_;
    Bool_t hltResult[100] = {false}; 
    TFile* fout;
    TTree* tout;
    unsigned int run, lumi, evtn;
    bool passL1, passHLT;
    float PV_x, PV_y, PV_z;
    //SV variables
    std::vector<unsigned int> index, ndof;
    std::vector<float> x, y, z;
    std::vector<float> xe, ye, ze;
    std::vector<float> chi2, prob, SV_chi2Ndof;
    std::vector<float> lxy, l3d;
    std::vector<float> mindx, mindy, mindz, mindxy, mind3d;
    std::vector<float> maxdx, maxdy, maxdz, maxdxy, maxd3d;
    std::vector<bool> SV_selected;
    std::vector<bool> onModule, onModuleWithinUnc;
    std::vector<float> closestDet_x, closestDet_y, closestDet_z;
    std::vector<float> minDistanceFromDet, minDistanceFromDet_x, minDistanceFromDet_y, minDistanceFromDet_z;
    //Muon variables
    std::vector<std::vector<int>> vtxIdxs;
    std::vector<unsigned int> saHits, saMatchedStats;
    std::vector<unsigned int> muHits, muChambs, muCSCDT, muMatch, muMatchedStats, muExpMatchedStats, muMatchedRPC;
    std::vector<unsigned int> pixHits, stripHits;
    std::vector<unsigned int> pixLayers, trkLayers;
    std::vector<int> bestAssocSVIdx, bestAssocSVOverlapIdx;
    std::vector<int> ch;
    std::vector<int> isGlobal, isTracker, isStandAlone;
    std::vector<float> pt, eta, phi;
    std::vector<float> mu_chi2Ndof;
    std::vector<float> ecalIso, hcalIso, trackIso;
    std::vector<float> ecalRelIso, hcalRelIso, trackRelIso;
    std::vector<float> dxy, dxye, dz, dze;
    std::vector<float> dxysig, dzsig;
    std::vector<float> phiCorr, dxyCorr;
    std::vector<float> PFIsoChg0p3, PFIsoChg0p4, PFIsoAll0p3, PFIsoAll0p4;
    std::vector<float> PFRelIsoChg0p3, PFRelIsoChg0p4, PFRelIsoAll0p3, PFRelIsoAll0p4;
    std::vector<float> mindrPF0p3, mindrPF0p4;
    std::vector<float> mindr, maxdr;
    std::vector<float> mindrJet, mindphiJet, mindetaJet;
    std::vector<TLorentzVector> vec;
    std::vector<bool> mu_selected;
    std::vector<int> nhitsbeforesv, ncompatible, ncompatibletotal, nexpectedhits, nexpectedhitsmultiple, nexpectedhitsmultipletotal, nexpectedhitstotal;
  };

//Constructor
SctParLooper::SctParLooper(const edm::ParameterSet& iConfig) :
  triggerCache_{triggerExpression::Data(iConfig.getParameterSet("triggerConfiguration"), consumesCollector())},
  vtriggerSelection_{iConfig.getParameter<vector<string>>("triggerSelection")},
  algToken_{consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("AlgInputTag"))},

  //l1ResultToken_{consumes<std::vector<bool>>(iConfig.getParameter<edm::InputTag>("l1Result"))},
  //hltResultToken_{consumes<std::vector<bool>>(iConfig.getParameter<edm::InputTag>("hltResult"))},
  //l1NameToken_{consumes<std::vector<std::string>>(iConfig.getParameter<edm::InputTag>("l1Names"))},
  //hltNameToken_{consumes<std::vector<std::string>>(iConfig.getParameter<edm::InputTag>("hltNames"))},
  pvToken_{consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("hltScoutingPrimaryVertexPacker_primaryVtx"))},
  svToken_{consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("hltScoutingMuonPacker_displacedVtx"))},
  muToken_{consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muons"))}
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
  edm::Handle<std::vector<Run3ScoutingVertex>> primaryVertices;
  edm::Handle<std::vector<Run3ScoutingVertex>> displacedVertices;
  edm::Handle<std::vector<Run3ScoutingMuon>> muons;

  iEvent.getByToken(pvToken_, primaryVertices);
  iEvent.getByToken(svToken_, displacedVertices);
  iEvent.getByToken(muToken_, muons);
  if (!primaryVertices.isValid()) {
    std::cout << "Primary vertex not valid" << std::endl;	  
    return;
  }
  if (!displacedVertices.isValid()) {
    std::cout << "Displaced vertex not valid" << std::endl;
    return;
  }
  if (!muons.isValid()) {
    std::cout << "Muons not valid" << std::endl;
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
  if (!passL1) {
    return;
  }
  

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
  if (!passHLT) {
    return;
  }
  //int nMuonAssoc;
  bool relaxedSVSel = false;
  float sfSVsel = (relaxedSVSel) ? 5.0 : 1.0;
  float maxXerr=0.05*sfSVsel, maxYerr=0.05*sfSVsel, maxZerr=0.10*sfSVsel, maxChi2=3.0*sfSVsel;
  
  PV_x = 0;
  PV_y = 0;
  PV_z = 0;

  if (primaryVertices.isValid() && !primaryVertices->empty() && (*primaryVertices)[0].isValidVtx()) {
    PV_x = (*primaryVertices)[0].x();
    PV_y = (*primaryVertices)[0].y();
    PV_z = (*primaryVertices)[0].z();
  }

  //SV selection
  index.clear(); ndof.clear();
  x.clear(); y.clear(); z.clear();
  xe.clear(); ye.clear(); ze.clear();
  chi2.clear(); prob.clear(); SV_chi2Ndof.clear();
  lxy.clear(); l3d.clear();
  mindx.clear(); mindy.clear(); mindz.clear(); mindxy.clear(); mind3d.clear();
  maxdx.clear(); maxdy.clear(); maxdz.clear(); maxdxy.clear(); maxd3d.clear();
  SV_selected.clear();
  onModule.clear(); onModuleWithinUnc.clear();
  closestDet_x.clear(); closestDet_y.clear(); closestDet_z.clear();
  minDistanceFromDet.clear(), minDistanceFromDet_x.clear(), minDistanceFromDet_y.clear(), minDistanceFromDet_z.clear();
  const auto& svCollection = *displacedVertices;
  unsigned int nSVs = svCollection.size();
  if (nSVs < 1){
    return;
  }
  for (unsigned int iSV = 0; iSV < nSVs; ++iSV) {
    const auto& sv = svCollection[iSV];
    if (!sv.isValidVtx()) continue;	
    index.push_back(iSV);
    ndof.push_back(sv.ndof());
    x.push_back(sv.x());
    y.push_back(sv.y());
    z.push_back(sv.z());
    xe.push_back(sv.xError());
    ye.push_back(sv.yError());
    ze.push_back(sv.zError());
    chi2.push_back(sv.chi2());
    prob.push_back(TMath::Prob(sv.chi2(), sv.ndof()));
    SV_chi2Ndof.push_back(sv.chi2()/sv.ndof());
    float dx = sv.x() - PV_x;
    float dy = sv.y() - PV_y;
    float dz = sv.z() - PV_z;
    float lxy_int = std::sqrt(dx * dx + dy * dy);
    float l3d_int = std::sqrt(lxy_int * lxy_int + dz * dz);

    lxy.push_back(lxy_int);
    l3d.push_back(l3d_int);

    bool passSel = (sv.xError() < maxXerr && sv.yError() < maxYerr && sv.zError() < maxZerr && sv.chi2() / sv.ndof() < maxChi2);
    SV_selected.push_back(passSel);
  }
  
  //Muon selection  
  vtxIdxs.clear();
  saHits.clear(); saMatchedStats.clear();
  muHits.clear(); muChambs.clear(); muCSCDT.clear(); muMatch.clear(); muMatchedStats.clear(); muExpMatchedStats.clear(); muMatchedRPC.clear();
  pixHits.clear(); stripHits.clear();
  pixLayers.clear(); trkLayers.clear();
  bestAssocSVIdx.clear(); bestAssocSVOverlapIdx.clear();
  ch.clear();
  isGlobal.clear(); isTracker.clear(); isStandAlone.clear();
  pt.clear(); eta.clear(); phi.clear();
  mu_chi2Ndof.clear();
  ecalIso.clear(); hcalIso.clear(); trackIso.clear();
  ecalRelIso.clear(); hcalRelIso.clear(); trackRelIso.clear();
  dxy.clear(); dxye.clear(); dz.clear(); dze.clear();
  dxysig.clear(); dzsig.clear();
  phiCorr.clear(); dxyCorr.clear();
  PFIsoChg0p3.clear(); PFIsoAll0p3.clear();
  PFRelIsoChg0p3.clear(); PFRelIsoAll0p3.clear();
  PFIsoChg0p4.clear(); PFIsoAll0p4.clear();
  PFRelIsoChg0p4.clear(); PFRelIsoAll0p4.clear();
  mindrPF0p3.clear(); mindrPF0p4.clear();
  mindr.clear(); maxdr.clear();
  mindrJet.clear(); mindphiJet.clear(); mindetaJet.clear();
  vec.clear();
  mu_selected.clear();
  nhitsbeforesv.clear(); ncompatible.clear(); ncompatibletotal.clear(); nexpectedhits.clear(); nexpectedhitsmultiple.clear(); nexpectedhitsmultipletotal.clear(); nexpectedhitstotal.clear();

  if (!muons.isValid()) {
    return;
  }
  const auto& muonCollection = *muons;
  unsigned int nMus = muonCollection.size();
  if (nMus < 1){
    return;
  }
  for (unsigned int iMu = 0; iMu < nMus; ++iMu) {
    const auto& mu = muonCollection[iMu];
    std::vector<int> matchedAndSelVtxIdxs;
    for (unsigned int matchedVtxIdx : mu.vtxIndx()) {
      for (size_t iSV = 0; iSV < index.size(); ++iSV) {
        if (matchedVtxIdx == index[iSV]) {
          matchedAndSelVtxIdxs.push_back(matchedVtxIdx);
        }
      }
    }

    vtxIdxs.push_back(matchedAndSelVtxIdxs);
    saHits.push_back(mu.nValidStandAloneMuonHits());
    saMatchedStats.push_back(mu.nStandAloneMuonMatchedStations());
    muHits.push_back(mu.nValidRecoMuonHits());
    muChambs.push_back(mu.nRecoMuonChambers());
    muCSCDT.push_back(mu.nRecoMuonChambersCSCorDT());
    muMatch.push_back(mu.nRecoMuonMatches());
    muMatchedStats.push_back(mu.nRecoMuonMatchedStations());
    muExpMatchedStats.push_back(mu.nRecoMuonExpectedMatchedStations());
    muMatchedRPC.push_back(mu.nRecoMuonMatchedRPCLayers());
    pixHits.push_back(mu.nValidPixelHits());
    stripHits.push_back(mu.nValidStripHits());
    pixLayers.push_back(mu.nPixelLayersWithMeasurement());
    trkLayers.push_back(mu.nTrackerLayersWithMeasurement());
    pt.push_back(mu.pt());
    eta.push_back(mu.eta());
    phi.push_back(mu.phi());
    ch.push_back(mu.charge());
    //isGlobal.push_back(isGlobalMuon(mu.type()));
    //isTracker.push_back(isTrackerMuon(mu.type()));
    //isStandAlone.push_back(isStandAloneMuon(mu.type()));

    mu_chi2Ndof.push_back(mu.normalizedChi2());

    ecalIso.push_back(mu.ecalIso());
    hcalIso.push_back(mu.hcalIso());
    trackIso.push_back(mu.trackIso());


    //Muons.ecalRelIso.push_back(mu.ecalIso() / pt);
    //Muons.hcalRelIso.push_back(mu.hcalIso() / pt);
    //Muons.trackRelIso.push_back(mu.trackIso() / pt);

    //Muons.dxy.push_back(mu.trk_dxy());
    //Muons.dxye.push_back(mu.trk_dxyError());
    //Muons.dz.push_back(mu.trk_dz());
    //Muons.dze.push_back(mu.trk_dzError());

    //Muons.dxysig.push_back(mu.trk_dxy() / mu.trk_dxyError());
    //Muons.dzsig.push_back(mu.trk_dz() / mu.trk_dzError());
    bool mu_selection = (mu.pt() > 3.0 && std::abs(mu.eta()) < 2.4 && mu.normalizedChi2() < 3.0);
    mu_selected.push_back(mu_selection);
  }

  tout->Fill();
      
}

void SctParLooper::beginJob() {
  fout = new TFile("output.root", "RECREATE");
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
  tout->Branch("PV_x", &PV_x);
  tout->Branch("PV_y", &PV_y);
  tout->Branch("PV_z", &PV_z);

  //SV branches
  tout->Branch("SV_xe", &xe);
  tout->Branch("SV_ye", &ye);
  tout->Branch("SV_ze", &ze);
  tout->Branch("SV_chi2", &chi2);
  tout->Branch("SV_prob", &prob);
  tout->Branch("SV_chi2Ndof", &SV_chi2Ndof);
  tout->Branch("SV_lxy", &lxy);
  tout->Branch("SV_l3d", &l3d);
  tout->Branch("SV_selected", &SV_selected);
  
  //Muon branches
  tout->Branch("Muon_vtxIdxs", &vtxIdxs);
  tout->Branch("Muon_saHits", &saHits);
  tout->Branch("Muon_saMatchedStats", &saMatchedStats);
  tout->Branch("Muon_muHits", &muHits);
  tout->Branch("Muon_muChambs", &muChambs);
  tout->Branch("Muon_muCSCDT", &muCSCDT);
  tout->Branch("Muon_muMatch", &muMatch);
  tout->Branch("Muon_muMatchedStats", &muMatchedStats);
  tout->Branch("Muon_muExpMatchedStats", &muExpMatchedStats);
  tout->Branch("Muon_muMatchedRPC", &muMatchedRPC);
  tout->Branch("Muon_pixHits", &pixHits);
  tout->Branch("Muon_stripHits", &stripHits);
  tout->Branch("Muon_pixLayers", &pixLayers);
  tout->Branch("Muon_trkLayers", &trkLayers);
  tout->Branch("Muon_pt", &pt);
  tout->Branch("Muon_eta", &eta);
  tout->Branch("Muon_phi", &phi);
  tout->Branch("Muon_selected", &mu_selected);
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



