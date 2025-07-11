import FWCore.ParameterSet.Config as cms

DoubleMuL1 = [
    "L1_DoubleMu0_Upt8_SQ_er2p0",
    "L1_DoubleMu0_Upt7_SQ_er2p0",
    "L1_DoubleMu_15_7",
    "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7",
    "L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18",
    "L1_DoubleMu8_SQ",
    "L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
    "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
    "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
    "L1_DoubleMu0_Upt15_Upt7",
    "L1_DoubleMu0_Upt6_IP_Min1_Upt4",
    "L1_DoubleMu0_Upt6_SQ_er2p0"
]
SingleMuL1 = ["L1_SingleMu11_SQ14_BMTF","L1_SingleMu10_SQ14_BMTF"]

DoubleMuHLT = ["DST_PFScouting_DoubleMuon_v*", "HLT_DoubleMu4_3_LowMass_v*"]
SingleMuHLT = ["DST_PFScouting_SingleMuon_v*", "HLT_Mu9_Barrel_L1HP10_IP6_v*", "HLT_Mu10_Barrel_L1HP11_IP6_v*", "HLT_Mu7_Barrel_L1HP8_IP6_v*", "HLT_Mu8_Barrel_L1HP9_IP6_v*", "HLT_Mu0_Barrel_L1HP6_IP6_v*", "HLT_Mu6_Barrel_L1HP7_IP6_v*"]
SctParLooper = cms.EDAnalyzer("SctParLooper",
    hltScoutingPrimaryVertexPacker_primaryVtx = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    hltScoutingMuonPacker_displacedVtx = cms.InputTag("hltScoutingMuonPackerNoVtx","displacedVtx","HLT"),
    Scoutingmuons = cms.InputTag("hltScoutingMuonPackerNoVtx"),
    PATMuonCollection = cms.InputTag('slimmedMuons'),
    PATVertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
    muon_TrackCollection = cms.InputTag('displacedGlobalMuons'),
    beamSpotCollection = cms.InputTag('offlineBeamSpot'),
    primaryVertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
    genParticlesInputTag = cms.InputTag('prunedGenParticles', '', 'PAT'),
    triggerSelection = cms.vstring(SingleMuHLT + DoubleMuHLT),
    triggerConfiguration = cms.PSet(
        hltResults            = cms.InputTag('TriggerResults','','HLT'),
        l1tResults            = cms.InputTag('','',''),
        l1tIgnoreMaskAndPrescale = cms.bool(False),
        throw                 = cms.bool(True),
        usePathStatus = cms.bool(False),
    ),
    AlgInputTag = cms.InputTag("gtStage2Digis"),
    l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
    l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
    ReadPrescalesFromFile = cms.bool(False),
    l1Seeds = cms.vstring(DoubleMuL1 + SingleMuL1) 
)
nTuplizer = cms.Sequence(SctParLooper)
