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
#SingleMuL1 = ["L1_SingleMu11_SQ14_BMTF","L1_SingleMu10_SQ14_BMTF"]
SctParLooper = cms.EDAnalyzer("SctParLooper",
    #l1Result = cms.InputTag("triggerMaker", "l1result"),
    #hltResult = cms.InputTag("triggerMaker", "hltresult"),
    #l1Names  = cms.InputTag("triggerMaker", "l1name"),
    #hltNames = cms.InputTag("triggerMaker", "hltname"),    
    hltScoutingPrimaryVertexPacker_primaryVtx = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
    hltScoutingMuonPacker_displacedVtx = cms.InputTag("hltScoutingMuonPackerNoVtx","displacedVtx","HLT"),
    muons = cms.InputTag("hltScoutingMuonPackerNoVtx"),
    triggerSelection = cms.vstring(["DST_PFScouting_ZeroBias_v*", "DST_PFScouting_DoubleEG_v*", "DST_PFScouting_JetHT_v*"]),
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
    l1Seeds = cms.vstring(DoubleMuL1), #Full list of double muon L1 seeds
)
nTuplizer = cms.Sequence(SctParLooper)
